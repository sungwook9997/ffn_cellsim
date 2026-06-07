# ffn_hoomd_plugin

Native HOOMD-blue extension sandbox for the FFN VRAM-resident hot loop.

This directory is intentionally outside the normal Python package build. Normal
`pytest` and editable installs do not build it. Use it to develop compiled
C++/CUDA operations before wiring them into `ffn_sim`.

## First Target

1. Compile/load probe (`_ffn_native.build_info`).
2. Native HOOMD updater attach probe.
3. Native GPU particle-state mutation probe.
4. Native Leimkuhler-Matthews BAOAB operation.
5. Fixed-pool myosin/xlink attachment force kernels.
6. Native custom forces for ERM, turgor, membrane, and manifold contact.

## Build Probe

On gbook:

```bash
cd ~/ffn_cellsim/native/ffn_hoomd_plugin
source ~/miniconda3/etc/profile.d/conda.sh
conda activate ffn_sim

# One-time toolchain bootstrap for the conda env:
conda install -y --override-channels -c conda-forge \
  cmake ninja cxx-compiler pybind11 eigen \
  cuda-nvcc=12.9 cuda-cudart-dev=12.9 cuda-libraries-dev=12.9

cmake -S . -B build -DCMAKE_PREFIX_PATH="$CONDA_PREFIX"
cmake --build build -j
PYTHONPATH="$PWD/build" python - <<'PY'
import _ffn_native
print(_ffn_native.build_info())
PY
```

The probe does not replace any runtime path. It only proves the external HOOMD
plugin toolchain can compile and load against the active conda HOOMD package.

## Smoke Tests

Use both the compiled extension and the experimental Python wrapper package:

```bash
export PYTHONPATH="$PWD/build:$PWD/python"
python smoke_noop_updater.py --device gpu --steps 16
python smoke_position_kick.py --device gpu --steps 8 --dx 0.125
python smoke_position_kick.py --device cpu --steps 8 --dx 0.125
```

Expected output:

```text
native-noop-updater smoke PASS device=gpu steps=16 count=16 last=15
native-position-kick smoke PASS device=gpu steps=8 x_final=1 count=8
native-position-kick smoke PASS device=cpu steps=8 x_final=1 count=8
```

The position-kick smoke proves that the native updater can mutate
`ParticleData::getPositions()` directly. On GPU it uses a device `ArrayHandle`
and a CUDA/HIP kernel.

## Overhead Probe

```bash
export PYTHONPATH="$PWD/build:$PWD/python"
for mode in none python-noop native-noop native-position-kick; do
  python benchmark_updater_overhead.py \
    --device gpu --mode "$mode" --steps 200000 --n-particles 1
done
```

Gbook snapshot on 2026-06-07:

| Mode | Steps/s |
|---|---:|
| no updater | 859,979 |
| Python no-op `CustomUpdater` | 763,130 |
| native no-op updater | 829,522 |
| native position-kick CUDA kernel | 330,985 |

## Native BAOAB Contract

The native operation must preserve the Python `OverdampedBAOABDevice` update:

```text
r_i <- wrap_box(
    r_i
    + F_i / gamma[type_i] * dt
    + sqrt(kT / gamma[type_i]) * (W_i,n + W_i,n-1) * dt
)
```

Required state:

- previous Gaussian vector per tag
- gamma by particle type/tag
- deterministic seed/timestep random stream
- GPU and CPU implementations, or a clear GPU-only experimental guard

Validation before production:

- existing `test_baoab_device.py`
- GPU-vs-CPU diffusion/drift/harmonic statistics
- H.7 `n_fil=150` and `n_fil=1000` speed/correctness smoke
- no gate loosening
