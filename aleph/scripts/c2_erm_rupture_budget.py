"""Does the ERM preload ask each tether for more tension than its rupture threshold allows?

`erm_tether_force_kernel` returns WITHOUT accumulating when tension > f_rupt, so a tether the preload
has stretched past threshold contributes exactly zero force — which is indistinguishable, in f_in,
from never having been preloaded at all. That is the bit-identical f_in observed across both sweeps.
"""
from aleph.components.incumbent.compartments import MEMBRANE_DEFAULTS as P
from aleph.laws.membrane_surface import GAMMA_MCA_PN_UM, erm_rupture_force

kappa_m = P["kappa_m"]
gamma_mem = P["gamma_mem"]
f_rupt = erm_rupture_force(kappa_m, gamma_mem)

print(f"kappa_m    = {kappa_m} pN*um")
print(f"gamma_mem  = {gamma_mem} pN/um")
print(f"gamma_MCA  = {GAMMA_MCA_PN_UM} pN/um")
print(f"f_rupt     = {f_rupt:.4f} pN     (docstring band 5-40 pN)")
print()
print("PRELOAD DEMAND, measured job 69 (full native):")
for label, need in (("p50", 37.28847443893447), ("absmax", 47.39767963656912)):
    verdict = "EXCEEDS f_rupt -> tether contributes ZERO force" if need > f_rupt else "within f_rupt"
    print(f"  |f_in| {label:>6} = {need:8.4f} pN   {verdict}")
print()
print("sweep 2 doubled the demand to ~74.6 pN, which is further past it.")
