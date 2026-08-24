"""Independent cross-checks for the cortex F-actin filament count (70,686).

Reproducer for docs/v2_audit/CORTEX_COUNT_CROSSCHECK_2026-07-21.md. All inputs are KB-cited; no fitting.
Shows that the current count (areal-density route) is ~9-15x above what the cortex mesh/concentration imply
at the model's 3 um filament length -- the (N, L, xi) consistency triangle can satisfy at most 2 of 3.
"""

from __future__ import annotations

import numpy as np

# ---- KB-sourced inputs (cited) ------------------------------------------------------------------------
R_UM = 7.5          # cell radius [um]  (Wagner2011 Coulter dia 14.8 um -> R~7.4; setpoint 7.5)
H_UM = 0.20         # cortex thickness [um] (KB-3.1/3.5: ~200 nm)
L_UM = 3.0          # MODEL filament contour length [um] ((beads-1)*l0 = 6*0.5)
XI_ACTIN_UM = 0.100  # actin meshwork pore [um] (KB-3.1/3.18: ~100 nm)
RISE_UM = 2.7e-3    # axial rise per strand [um/monomer] (KB-DRAFT-3-26)
STRANDS = 2         # F-actin = 2 protofilaments
R_FIL_UM = 0.0035   # F-actin radius [um] (7 nm diameter, KB-3.18)
NA = 6.022e23
MW_ACTIN = 42000.0  # g/mol

A_UM2 = 4.0 * np.pi * R_UM**2           # cortex surface [um^2]
V_UM3 = A_UM2 * H_UM                     # cortex shell volume [um^3]
LIN_MON_PER_UM = STRANDS / RISE_UM       # monomers per um of filament (~741)
VOL_PER_UM = np.pi * R_FIL_UM**2         # F-actin volume per um length [um^3]


def conc_uM(total_len_um: float) -> float:
    """Volumetric F-actin concentration [uM] implied by a total contour length in the cortex shell."""
    mol = (total_len_um * LIN_MON_PER_UM) / NA
    return mol / (V_UM3 * 1e-15) * 1e6


def vol_frac_pct(total_len_um: float) -> float:
    """F-actin volume fraction [%] implied by a total contour length."""
    return (total_len_um * VOL_PER_UM) / V_UM3 * 100.0


def n_from_len(total_len_um: float) -> float:
    """Filament count from total contour length at the model filament length L_UM."""
    return total_len_um / L_UM


def main() -> None:
    print(f"cortex A = {A_UM2:.1f} um^2 | V = {V_UM3:.1f} um^3 | linear density {LIN_MON_PER_UM:.0f} mon/um\n")

    # Route 1 -- areal density x surface (current baseline)
    n1 = 100.0 * A_UM2
    lt1 = n1 * L_UM
    la1 = lt1 / A_UM2
    print("ROUTE 1  areal density x surface (CURRENT)")
    print(f"  N = 100/um^2 x {A_UM2:.0f} = {n1:.0f}")
    print(f"  total len {lt1:.0f} um | conc {conc_uM(lt1):.0f} uM | vol-frac {vol_frac_pct(lt1):.1f}%"
          f" | implied mesh (square 2/L_A) {2/la1*1000:.0f} nm\n")

    # Route 2 -- mesh size (square-grid geometric)
    la2 = 2.0 / XI_ACTIN_UM
    lt2 = la2 * A_UM2
    print("ROUTE 2  actin mesh xi=100nm -> L_A=2/xi (KB-3.1/3.18)")
    print(f"  N = {n_from_len(lt2):.0f} | conc {conc_uM(lt2):.0f} uM | vol-frac {vol_frac_pct(lt2):.2f}%\n")

    # Route 3 -- mesh via Schmidt/Kaes xi[um] = 0.3/sqrt(c[mg/mL])
    c_mgml = (0.3 / XI_ACTIN_UM) ** 2
    c_molL = c_mgml / MW_ACTIN            # mg/mL = g/L; /(g/mol) = mol/L
    mon3 = c_molL * (V_UM3 * 1e-15) * NA
    lt3 = mon3 / LIN_MON_PER_UM
    print("ROUTE 3  mesh -> concentration via xi-c (Schmidt/Kaes)")
    print(f"  xi=100nm -> c = {c_mgml:.1f} mg/mL = {c_molL*1e6:.0f} uM")
    print(f"  N = {n_from_len(lt3):.0f} | vol-frac {vol_frac_pct(lt3):.2f}%\n")

    # Consistency triangle for the baseline
    print("CONSISTENCY TRIANGLE (N, L, xi) -- at most 2 of 3")
    lt_b = 70686 * L_UM
    la_b = lt_b / A_UM2
    print(f"  N=70686 & L=3um  -> mesh {2/la_b*1000:.0f} nm (square), conc {conc_uM(lt_b):.0f} uM,"
          f" vol-frac {vol_frac_pct(lt_b):.1f}%")
    print(f"  N=70686 & xi=100nm -> filament length must be {lt2/70686*1000:.0f} nm (not 3um)")


if __name__ == "__main__":
    main()
