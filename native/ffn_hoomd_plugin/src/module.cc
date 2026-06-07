#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "baoab_kernel.cuh"
#include "constrained_kernel.cuh"
#include "fixman_kernel.cuh"
#include "position_kick_kernel.cuh"
#include "shake_kernel.cuh"

#include <hoomd/BoxDim.h>
#include <hoomd/ExecutionConfiguration.h>
#include <hoomd/GPUArray.h>
#include <hoomd/ParticleData.h>
#include <hoomd/SystemDefinition.h>
#include <hoomd/Trigger.h>
#include <hoomd/Updater.h>

#include <cmath>
#include <cstdint>
#include <cstring>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace py = pybind11;

class FFNNoOpUpdater : public hoomd::Updater
    {
    public:
    FFNNoOpUpdater(std::shared_ptr<hoomd::SystemDefinition> sysdef,
                   std::shared_ptr<hoomd::Trigger> trigger)
        : hoomd::Updater(sysdef, trigger), m_update_count(0), m_last_timestep(0)
        {
        }

    void update(uint64_t timestep) override
        {
        ++m_update_count;
        m_last_timestep = timestep;
        }

    uint64_t getUpdateCount() const
        {
        return m_update_count;
        }

    uint64_t getLastTimestep() const
        {
        return m_last_timestep;
        }

    private:
    uint64_t m_update_count;
    uint64_t m_last_timestep;
    };

class FFNPositionKickUpdater : public hoomd::Updater
    {
    public:
    FFNPositionKickUpdater(std::shared_ptr<hoomd::SystemDefinition> sysdef,
                           std::shared_ptr<hoomd::Trigger> trigger,
                           hoomd::Scalar dx,
                           hoomd::Scalar dy,
                           hoomd::Scalar dz)
        : hoomd::Updater(sysdef, trigger), m_dx(dx), m_dy(dy), m_dz(dz), m_update_count(0),
          m_last_timestep(0), m_block_size(256)
        {
        }

    void update(uint64_t timestep) override
        {
        ++m_update_count;
        m_last_timestep = timestep;

        if (m_exec_conf->isCUDAEnabled())
            {
            m_exec_conf->setDevice();
            hoomd::ArrayHandle<hoomd::Scalar4> d_pos(m_pdata->getPositions(),
                                                     hoomd::access_location::device,
                                                     hoomd::access_mode::readwrite);
#ifdef ENABLE_HIP
            const hipError_t kernel_status = ffn_native::gpu_position_kick(
                d_pos.data, m_pdata->getN(), m_dx, m_dy, m_dz, m_block_size);
            if (kernel_status != hipSuccess)
                {
                throw std::runtime_error(std::string("FFNPositionKickUpdater GPU kernel failed: ")
                                         + hipGetErrorString(kernel_status));
                }
            if (m_exec_conf->isCUDAErrorCheckingEnabled())
                {
                CHECK_CUDA_ERROR();
                }
#else
            throw std::runtime_error("FFNPositionKickUpdater GPU path requires ENABLE_HIP");
#endif
            return;
            }

        hoomd::ArrayHandle<hoomd::Scalar4> h_pos(m_pdata->getPositions(),
                                                 hoomd::access_location::host,
                                                 hoomd::access_mode::readwrite);
        for (unsigned int i = 0; i < m_pdata->getN(); ++i)
            {
            h_pos.data[i].x += m_dx;
            h_pos.data[i].y += m_dy;
            h_pos.data[i].z += m_dz;
            }
        }

    uint64_t getUpdateCount() const
        {
        return m_update_count;
        }

    uint64_t getLastTimestep() const
        {
        return m_last_timestep;
        }

    private:
    hoomd::Scalar m_dx;
    hoomd::Scalar m_dy;
    hoomd::Scalar m_dz;
    uint64_t m_update_count;
    uint64_t m_last_timestep;
    unsigned int m_block_size;
    };

// Native overdamped Leimkuhler-Matthews BAOAB-limit integrator (Stage 1a).
// Mirrors ffn_sim.integrator.baoab_device.OverdampedBAOABDevice: attach to an
// md.Integrator(forces=[...], methods=[]); HOOMD computes net force each step,
// this Updater reads it and advances positions device-resident. gamma_by_tag is
// the per-tag Stokes drag (dense [0,N) tag space); bd_prefactor = sqrt(kT/(2 g dt)).
class FFNBaoabUpdater : public hoomd::Updater
    {
    public:
    FFNBaoabUpdater(std::shared_ptr<hoomd::SystemDefinition> sysdef,
                    std::shared_ptr<hoomd::Trigger> trigger,
                    hoomd::Scalar dt,
                    uint64_t seed,
                    hoomd::Scalar kT,
                    std::vector<hoomd::Scalar> gamma_by_tag)
        : hoomd::Updater(sysdef, trigger), m_dt(dt), m_seed(seed), m_block_size(256),
          m_n_tag(gamma_by_tag.size()), m_gamma(gamma_by_tag.size(), m_exec_conf),
          m_bd_pref(gamma_by_tag.size(), m_exec_conf), m_prv(3 * gamma_by_tag.size(), m_exec_conf)
        {
        if (gamma_by_tag.empty())
            throw std::runtime_error("FFNBaoabUpdater: gamma_by_tag is empty");
        if (!(dt > hoomd::Scalar(0)))
            throw std::runtime_error("FFNBaoabUpdater: dt must be > 0");
            {
            hoomd::ArrayHandle<hoomd::Scalar> h_g(m_gamma, hoomd::access_location::host,
                                                  hoomd::access_mode::overwrite);
            hoomd::ArrayHandle<hoomd::Scalar> h_p(m_bd_pref, hoomd::access_location::host,
                                                  hoomd::access_mode::overwrite);
            for (size_t t = 0; t < gamma_by_tag.size(); ++t)
                {
                if (!(gamma_by_tag[t] > hoomd::Scalar(0)))
                    throw std::runtime_error("FFNBaoabUpdater: gamma must be > 0");
                h_g.data[t] = gamma_by_tag[t];
                h_p.data[t]
                    = std::sqrt(kT / (hoomd::Scalar(2) * gamma_by_tag[t] * dt));
                }
            }
            {
            hoomd::ArrayHandle<hoomd::Scalar> h_prv(m_prv, hoomd::access_location::host,
                                                    hoomd::access_mode::overwrite);
            for (size_t k = 0; k < 3 * gamma_by_tag.size(); ++k)
                h_prv.data[k] = hoomd::Scalar(0);
            }
        }

    void update(uint64_t timestep) override
        {
        const unsigned int N = m_pdata->getN();
        if (N == 0)
            return;
        if (N > m_n_tag)
            throw std::runtime_error("FFNBaoabUpdater: N exceeds the gamma/prv buffer "
                                     "(tag space grew; this Action is fixed-N).");
        const hoomd::BoxDim& box = m_pdata->getGlobalBox();
        const hoomd::Scalar3 L = box.getL();
        const hoomd::Scalar xy = box.getTiltFactorXY();
        const hoomd::Scalar xz = box.getTiltFactorXZ();
        const hoomd::Scalar yz = box.getTiltFactorYZ();

        if (!m_exec_conf->isCUDAEnabled())
            throw std::runtime_error(
                "FFNBaoabUpdater currently requires a GPU device (ENABLE_HIP).");

#ifdef ENABLE_HIP
        m_exec_conf->setDevice();
        hoomd::ArrayHandle<hoomd::Scalar4> d_pos(
            m_pdata->getPositions(), hoomd::access_location::device, hoomd::access_mode::readwrite);
        hoomd::ArrayHandle<hoomd::Scalar4> d_force(
            m_pdata->getNetForce(), hoomd::access_location::device, hoomd::access_mode::read);
        hoomd::ArrayHandle<int3> d_image(
            m_pdata->getImages(), hoomd::access_location::device, hoomd::access_mode::readwrite);
        hoomd::ArrayHandle<unsigned int> d_tag(
            m_pdata->getTags(), hoomd::access_location::device, hoomd::access_mode::read);
        hoomd::ArrayHandle<hoomd::Scalar> d_g(
            m_gamma, hoomd::access_location::device, hoomd::access_mode::read);
        hoomd::ArrayHandle<hoomd::Scalar> d_p(
            m_bd_pref, hoomd::access_location::device, hoomd::access_mode::read);
        hoomd::ArrayHandle<hoomd::Scalar> d_prv(
            m_prv, hoomd::access_location::device, hoomd::access_mode::readwrite);
        const hipError_t st = ffn_native::gpu_baoab_step(
            d_pos.data, d_force.data, d_image.data, d_tag.data, d_g.data, d_p.data, d_prv.data, N,
            m_seed, timestep, m_dt, L.x, L.y, L.z, xy, xz, yz, m_block_size);
        if (st != hipSuccess)
            throw std::runtime_error(std::string("FFNBaoabUpdater kernel failed: ")
                                     + hipGetErrorString(st));
        if (m_exec_conf->isCUDAErrorCheckingEnabled())
            {
            CHECK_CUDA_ERROR();
            }
#else
        throw std::runtime_error("FFNBaoabUpdater GPU path requires ENABLE_HIP");
#endif
        }

    private:
    hoomd::Scalar m_dt;
    uint64_t m_seed;
    unsigned int m_block_size;
    size_t m_n_tag;
    hoomd::GPUArray<hoomd::Scalar> m_gamma;
    hoomd::GPUArray<hoomd::Scalar> m_bd_pref;
    hoomd::GPUArray<hoomd::Scalar> m_prv;
    };

#ifdef ENABLE_HIP
// Standalone native M-SHAKE projection (Stage 1b validation handle). Takes numpy
// arrays, returns (projected_pos, lambda, nonconverged) — for bit-parity testing
// against ffn_sim.integrator.constrained_baoab.shake_project_chains.
py::tuple shake_project(
    py::array_t<double, py::array::c_style | py::array::forcecast> pos,
    py::array_t<double, py::array::c_style | py::array::forcecast> ref,
    py::array_t<double, py::array::c_style | py::array::forcecast> inv_mass,
    py::array_t<int, py::array::c_style | py::array::forcecast> chains,
    double rest_length,
    double Lx,
    double Ly,
    double Lz,
    double tol,
    unsigned int max_iter)
    {
    const auto posb = pos.request();
    const auto chb = chains.request();
    if (posb.ndim != 2 || posb.shape[1] != 3)
        throw std::runtime_error("pos must be (N, 3)");
    if (chb.ndim != 2)
        throw std::runtime_error("chains must be (F, m+1)");
    const unsigned int N = static_cast<unsigned int>(posb.shape[0]);
    const unsigned int F = static_cast<unsigned int>(chb.shape[0]);
    const unsigned int m = static_cast<unsigned int>(chb.shape[1]) - 1u;

    py::array_t<double> out_pos(std::vector<py::ssize_t>{(py::ssize_t)N, 3});
    std::memcpy(out_pos.mutable_data(), pos.data(), (size_t)3 * N * sizeof(double));
    py::array_t<double> out_lam(std::vector<py::ssize_t>{(py::ssize_t)F, (py::ssize_t)m});
    py::array_t<int> out_nc(std::vector<py::ssize_t>{(py::ssize_t)F});

    const hipError_t err = ffn_native::shake_project_host(
        out_pos.mutable_data(), ref.data(), inv_mass.data(), chains.data(),
        out_lam.mutable_data(), N, F, m, rest_length, Lx, Ly, Lz, tol, max_iter,
        out_nc.mutable_data());
    if (err != hipSuccess)
        throw std::runtime_error(std::string("shake_project kernel: ") + hipGetErrorString(err));
    return py::make_tuple(out_pos, out_lam, out_nc);
    }

// Standalone native Fixman pseudo-force (Stage 1b validation handle). Returns
// (force, U_F, bad_sign) for parity vs fixman_logdet_and_force.
py::tuple fixman_force(
    py::array_t<double, py::array::c_style | py::array::forcecast> pos,
    py::array_t<double, py::array::c_style | py::array::forcecast> inv_gamma,
    py::array_t<int, py::array::c_style | py::array::forcecast> chains,
    double kT,
    double Lx,
    double Ly,
    double Lz)
    {
    const auto posb = pos.request();
    const auto chb = chains.request();
    if (posb.ndim != 2 || posb.shape[1] != 3)
        throw std::runtime_error("pos must be (N, 3)");
    if (chb.ndim != 2)
        throw std::runtime_error("chains must be (F, m+1)");
    const unsigned int N = static_cast<unsigned int>(posb.shape[0]);
    const unsigned int F = static_cast<unsigned int>(chb.shape[0]);
    const unsigned int m = static_cast<unsigned int>(chb.shape[1]) - 1u;

    py::array_t<double> force(std::vector<py::ssize_t>{(py::ssize_t)N, 3});
    py::array_t<double> logdet(std::vector<py::ssize_t>{(py::ssize_t)F});
    py::array_t<int> bad(std::vector<py::ssize_t>{(py::ssize_t)F});
    const hipError_t err = ffn_native::fixman_force_host(
        pos.data(), inv_gamma.data(), chains.data(), force.mutable_data(),
        logdet.mutable_data(), bad.mutable_data(), N, F, m, Lx, Ly, Lz);
    if (err != hipSuccess)
        throw std::runtime_error(std::string("fixman_force kernel: ") + hipGetErrorString(err));

    const double half_kT = 0.5 * kT; // FIXMAN_SIGN = +1
    double* fd = force.mutable_data();
    for (size_t i = 0; i < (size_t)3 * N; ++i)
        fd[i] *= half_kT;
    const double* ld = logdet.data();
    double U_F = 0.0;
    for (unsigned int f = 0; f < F; ++f)
        U_F += ld[f];
    U_F *= half_kT;
    return py::make_tuple(force, U_F, bad);
    }
#endif

// Native device-resident constrained L-M BAOAB integrator (Stage 1b assembly).
// Mirrors ffn_sim.integrator.constrained_baoab.ConstrainedLeimkuhlerMatthewsBAOAB:
// methods=[] integrator; HOOMD computes net force (angle/LJ/ERM/…, the stretch
// bond is REPLACED by the constraint); this updater adds Fixman + predictor +
// M-SHAKE each step on the GPU. inv_gamma_by_tag = per-tag mobility (1/γ);
// chains_tag = flat F*(m+1) bead TAGS (uniform length m+1).
class FFNConstrainedBaoabUpdater : public hoomd::Updater
    {
    public:
    FFNConstrainedBaoabUpdater(std::shared_ptr<hoomd::SystemDefinition> sysdef,
                               std::shared_ptr<hoomd::Trigger> trigger,
                               hoomd::Scalar dt,
                               uint64_t seed,
                               hoomd::Scalar kT,
                               std::vector<hoomd::Scalar> inv_gamma_by_tag,
                               std::vector<int> chains_tag,
                               unsigned int F,
                               unsigned int m,
                               double rest_length,
                               double tol,
                               unsigned int max_iter)
        : hoomd::Updater(sysdef, trigger), m_dt(dt), m_half_kT(0.5 * (double)kT), m_seed(seed),
          m_rest_length(rest_length), m_tol(tol), m_max_iter(max_iter), m_block(128), m_F(F), m_m(m),
          m_N(m_pdata->getN()),
          m_inv_gamma_by_tag(inv_gamma_by_tag.size(), m_exec_conf),
          m_bd_pref_by_tag(inv_gamma_by_tag.size(), m_exec_conf), m_prv(3 * m_N, m_exec_conf),
          m_chains_tag(F * (m + 1) > 0 ? F * (m + 1) : 1, m_exec_conf),
          m_pos_d(3 * m_N, m_exec_conf), m_ref_d(3 * m_N, m_exec_conf),
          m_ffix_d(3 * m_N, m_exec_conf), m_pred_d(3 * m_N, m_exec_conf),
          m_inv_gamma_row(m_N, m_exec_conf), m_row_of_tag(m_N, m_exec_conf),
          m_chains_row(F * (m + 1) > 0 ? F * (m + 1) : 1, m_exec_conf),
          m_lambda(F * m > 0 ? F * m : 1, m_exec_conf), m_logdet(F > 0 ? F : 1, m_exec_conf),
          m_nonconv(F > 0 ? F : 1, m_exec_conf), m_bad_sign(F > 0 ? F : 1, m_exec_conf)
        {
        if (inv_gamma_by_tag.size() < m_N)
            throw std::runtime_error("FFNConstrainedBaoabUpdater: inv_gamma_by_tag shorter than N");
        if (!(dt > 0))
            throw std::runtime_error("FFNConstrainedBaoabUpdater: dt must be > 0");
        const double half_dt = 2.0 * (double)dt;
            {
            hoomd::ArrayHandle<hoomd::Scalar> h_ig(m_inv_gamma_by_tag,
                hoomd::access_location::host, hoomd::access_mode::overwrite);
            hoomd::ArrayHandle<hoomd::Scalar> h_bp(m_bd_pref_by_tag,
                hoomd::access_location::host, hoomd::access_mode::overwrite);
            for (size_t t = 0; t < inv_gamma_by_tag.size(); ++t)
                {
                h_ig.data[t] = inv_gamma_by_tag[t];
                h_bp.data[t] = std::sqrt((double)kT * (double)inv_gamma_by_tag[t] / half_dt);
                }
            }
            {
            hoomd::ArrayHandle<hoomd::Scalar> h_prv(m_prv, hoomd::access_location::host,
                hoomd::access_mode::overwrite);
            for (size_t k = 0; k < 3 * (size_t)m_N; ++k)
                h_prv.data[k] = 0.0;
            }
        if (F * (m + 1) > 0)
            {
            hoomd::ArrayHandle<int> h_ct(m_chains_tag, hoomd::access_location::host,
                hoomd::access_mode::overwrite);
            for (size_t k = 0; k < (size_t)F * (m + 1); ++k)
                h_ct.data[k] = chains_tag[k];
            }
        }

    void update(uint64_t timestep) override
        {
        const unsigned int N = m_pdata->getN();
        if (N == 0)
            return;
        if (N != m_N)
            throw std::runtime_error("FFNConstrainedBaoabUpdater: N changed (fixed-N Action).");
        const hoomd::BoxDim& box = m_pdata->getGlobalBox();
        const hoomd::Scalar3 L = box.getL();
        const hoomd::Scalar xy = box.getTiltFactorXY();
        const hoomd::Scalar xz = box.getTiltFactorXZ();
        const hoomd::Scalar yz = box.getTiltFactorYZ();
        if (!m_exec_conf->isCUDAEnabled())
            throw std::runtime_error("FFNConstrainedBaoabUpdater requires a GPU device.");
#ifdef ENABLE_HIP
        m_exec_conf->setDevice();
        hoomd::ArrayHandle<hoomd::Scalar4> d_pos(m_pdata->getPositions(),
            hoomd::access_location::device, hoomd::access_mode::readwrite);
        hoomd::ArrayHandle<hoomd::Scalar4> d_force(m_pdata->getNetForce(),
            hoomd::access_location::device, hoomd::access_mode::read);
        hoomd::ArrayHandle<int3> d_image(m_pdata->getImages(),
            hoomd::access_location::device, hoomd::access_mode::readwrite);
        hoomd::ArrayHandle<unsigned int> d_tag(m_pdata->getTags(),
            hoomd::access_location::device, hoomd::access_mode::read);
        hoomd::ArrayHandle<hoomd::Scalar> d_ig(m_inv_gamma_by_tag,
            hoomd::access_location::device, hoomd::access_mode::read);
        hoomd::ArrayHandle<hoomd::Scalar> d_bp(m_bd_pref_by_tag,
            hoomd::access_location::device, hoomd::access_mode::read);
        hoomd::ArrayHandle<hoomd::Scalar> d_prv(m_prv,
            hoomd::access_location::device, hoomd::access_mode::readwrite);
        hoomd::ArrayHandle<int> d_ct(m_chains_tag,
            hoomd::access_location::device, hoomd::access_mode::read);
        hoomd::ArrayHandle<hoomd::Scalar> d_pd(m_pos_d, hoomd::access_location::device,
            hoomd::access_mode::overwrite);
        hoomd::ArrayHandle<hoomd::Scalar> d_rd(m_ref_d, hoomd::access_location::device,
            hoomd::access_mode::overwrite);
        hoomd::ArrayHandle<hoomd::Scalar> d_fd(m_ffix_d, hoomd::access_location::device,
            hoomd::access_mode::overwrite);
        hoomd::ArrayHandle<hoomd::Scalar> d_prd(m_pred_d, hoomd::access_location::device,
            hoomd::access_mode::overwrite);
        hoomd::ArrayHandle<hoomd::Scalar> d_igr(m_inv_gamma_row,
            hoomd::access_location::device, hoomd::access_mode::overwrite);
        hoomd::ArrayHandle<int> d_rot(m_row_of_tag, hoomd::access_location::device,
            hoomd::access_mode::overwrite);
        hoomd::ArrayHandle<int> d_cr(m_chains_row, hoomd::access_location::device,
            hoomd::access_mode::overwrite);
        hoomd::ArrayHandle<hoomd::Scalar> d_lam(m_lambda, hoomd::access_location::device,
            hoomd::access_mode::overwrite);
        hoomd::ArrayHandle<hoomd::Scalar> d_ld(m_logdet, hoomd::access_location::device,
            hoomd::access_mode::overwrite);
        hoomd::ArrayHandle<int> d_nc(m_nonconv, hoomd::access_location::device,
            hoomd::access_mode::overwrite);
        hoomd::ArrayHandle<int> d_bs(m_bad_sign, hoomd::access_location::device,
            hoomd::access_mode::overwrite);
        const hipError_t st = ffn_native::gpu_constrained_step(
            d_pos.data, d_force.data, d_image.data, d_tag.data, d_ig.data, d_bp.data, d_prv.data,
            d_ct.data, d_pd.data, d_rd.data, d_fd.data, d_prd.data, d_igr.data, d_rot.data,
            d_cr.data, d_lam.data, d_nc.data, d_bs.data, d_ld.data, N, m_F, m_m, m_seed, timestep,
            m_dt, m_half_kT, m_rest_length, L.x, L.y, L.z, xy, xz, yz, m_tol, m_max_iter, m_block);
        if (st != hipSuccess)
            throw std::runtime_error(std::string("FFNConstrainedBaoabUpdater kernel: ")
                                     + hipGetErrorString(st));
        if (m_exec_conf->isCUDAErrorCheckingEnabled())
            {
            CHECK_CUDA_ERROR();
            }
#endif
        }

    // Accumulated per-bond Lagrange multipliers from the LAST step's M-SHAKE
    // (shape (F, m)). Per RIGID_LAGRANGE_TENSION_DESIGN.md the consumer converts
    // each λ to scalar bond tension T_bond = λ·r0/dt for the method-of-planes γ.
    pybind11::array_t<double> get_lambda() const
        {
        pybind11::array_t<double> out(
            std::vector<pybind11::ssize_t>{(pybind11::ssize_t)m_F, (pybind11::ssize_t)m_m});
        if (m_F * m_m == 0)
            return out;
        hoomd::ArrayHandle<hoomd::Scalar> h(m_lambda, hoomd::access_location::host,
                                            hoomd::access_mode::read);
        std::memcpy(out.mutable_data(), h.data, (size_t)m_F * m_m * sizeof(double));
        return out;
        }

    // Number of chains whose M-SHAKE did not converge last step (0 = all OK).
    int nonconverged_count() const
        {
        if (m_F == 0)
            return 0;
        hoomd::ArrayHandle<int> h(m_nonconv, hoomd::access_location::host,
                                  hoomd::access_mode::read);
        int s = 0;
        for (unsigned int f = 0; f < m_F; ++f)
            s += h.data[f];
        return s;
        }

    private:
    hoomd::Scalar m_dt;
    double m_half_kT;
    uint64_t m_seed;
    double m_rest_length;
    double m_tol;
    unsigned int m_max_iter;
    unsigned int m_block;
    unsigned int m_F;
    unsigned int m_m;
    unsigned int m_N;
    hoomd::GPUArray<hoomd::Scalar> m_inv_gamma_by_tag, m_bd_pref_by_tag, m_prv;
    hoomd::GPUArray<int> m_chains_tag;
    hoomd::GPUArray<hoomd::Scalar> m_pos_d, m_ref_d, m_ffix_d, m_pred_d, m_inv_gamma_row;
    hoomd::GPUArray<int> m_row_of_tag, m_chains_row;
    hoomd::GPUArray<hoomd::Scalar> m_lambda, m_logdet;
    hoomd::GPUArray<int> m_nonconv, m_bad_sign;
    };

std::string build_info()
    {
    return "ffn_hoomd_plugin compile/load probe: native hot-loop scaffold";
    }

PYBIND11_MODULE(_ffn_native, m)
    {
    m.doc() = "FFN native HOOMD hot-loop extension scaffold.";
    m.def("build_info", &build_info);
#ifdef ENABLE_HIP
    m.def("shake_project", &shake_project, py::arg("pos"), py::arg("ref"),
          py::arg("inv_mass"), py::arg("chains"), py::arg("rest_length"), py::arg("Lx"),
          py::arg("Ly"), py::arg("Lz"), py::arg("tol") = 1e-10, py::arg("max_iter") = 100);
    m.def("fixman_force", &fixman_force, py::arg("pos"), py::arg("inv_gamma"),
          py::arg("chains"), py::arg("kT"), py::arg("Lx"), py::arg("Ly"), py::arg("Lz"));
#endif
    py::class_<FFNNoOpUpdater, hoomd::Updater, std::shared_ptr<FFNNoOpUpdater>>(
        m, "FFNNoOpUpdater")
        .def(py::init<std::shared_ptr<hoomd::SystemDefinition>, std::shared_ptr<hoomd::Trigger>>())
        .def("get_update_count", &FFNNoOpUpdater::getUpdateCount)
        .def("get_last_timestep", &FFNNoOpUpdater::getLastTimestep);
    py::class_<FFNPositionKickUpdater, hoomd::Updater, std::shared_ptr<FFNPositionKickUpdater>>(
        m, "FFNPositionKickUpdater")
        .def(py::init<std::shared_ptr<hoomd::SystemDefinition>,
                      std::shared_ptr<hoomd::Trigger>,
                      hoomd::Scalar,
                      hoomd::Scalar,
                      hoomd::Scalar>())
        .def("get_update_count", &FFNPositionKickUpdater::getUpdateCount)
        .def("get_last_timestep", &FFNPositionKickUpdater::getLastTimestep);
    py::class_<FFNBaoabUpdater, hoomd::Updater, std::shared_ptr<FFNBaoabUpdater>>(
        m, "FFNBaoabUpdater")
        .def(py::init<std::shared_ptr<hoomd::SystemDefinition>,
                      std::shared_ptr<hoomd::Trigger>,
                      hoomd::Scalar,
                      uint64_t,
                      hoomd::Scalar,
                      std::vector<hoomd::Scalar>>());
    py::class_<FFNConstrainedBaoabUpdater, hoomd::Updater,
               std::shared_ptr<FFNConstrainedBaoabUpdater>>(m, "FFNConstrainedBaoabUpdater")
        .def(py::init<std::shared_ptr<hoomd::SystemDefinition>,
                      std::shared_ptr<hoomd::Trigger>,
                      hoomd::Scalar,
                      uint64_t,
                      hoomd::Scalar,
                      std::vector<hoomd::Scalar>,
                      std::vector<int>,
                      unsigned int,
                      unsigned int,
                      double,
                      double,
                      unsigned int>())
        .def("get_lambda", &FFNConstrainedBaoabUpdater::get_lambda)
        .def("nonconverged_count", &FFNConstrainedBaoabUpdater::nonconverged_count);
    }
