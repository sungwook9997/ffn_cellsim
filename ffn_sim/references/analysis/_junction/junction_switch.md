# Junction switch (cadherin -> integrin clutch) + bulk pressure — calibrated values

## Bulk (compressive) pressure in spheroids

Cells in the interior of a growing spheroid experience **compressive solid stress**
that builds toward the core as the tissue is confined. Reported magnitudes:

- Externally applied / internally generated compressive stress in the
  **0.5 – 10 kPa** range mechanically suppresses proliferation and triggers
  motility/jamming transitions (Helmlinger et al. 1997 Nat Biotech; Montel et al.
  2011 PRL ~5–10 kPa; Cheng et al. 2009 PLoS ONE).
- **Onset ≈ 0.5 kPa** for the growth-suppression / mechanotransduction response;
  **saturation ≈ 5 kPa** (response plateaus). We map our per-cell crowding proxy
  onto this **[0.5, 5] kPa** band.

### Pressure proxy (mechanistic, not a paper closed-form)

Per cell we count **neighbour crowding** = number of OTHER cells whose centroid is
within a contact radius r_contact = 2.2·R_patch, and convert to a compressive
stress via a linear map calibrated so an interior cell with a full coordination
shell (~12 FCC neighbours) sits near the saturation stress and a free-surface cell
(~4–6 neighbours) sits near/below onset:

```
crowd       = #neighbour cells within r_contact
P_kPa       = P_min + (P_max - P_min) * clip((crowd - c_lo)/(c_hi - c_lo), 0, 1)
             with c_lo = 4, c_hi = 12, P_min = 0.0, P_max = 6.0 kPa
```

So core (crowd ~12) -> ~6 kPa (above the 5-kPa saturation), rim (crowd ~5) -> ~0.5 kPa.

## Junction switch: cadherin -> integrin clutch

Under elevated compressive stress / at the unjamming free surface, cells
**down-regulate cadherin (cell–cell) adhesion and up-regulate integrin (cell–ECM /
substrate) adhesion** — the cadherin→integrin clutch switch that lets rim cells
unjam and migrate outward (EMT-like; Friedl & Alexander 2011 Cell; Trepat &
Sahai 2018; mechanically-gated E-cad/integrin reciprocity, Mui et al. 2016 JCS).

### Switch rule used by JunctionSwitchUpdater

```
P_switch = 0.5 kPa            # onset threshold
if P_cell > P_switch:
    cadherin (cell-cell W_cc)      *= cadherin_weak_factor   (0.3, weaken)
    integrin (cell-substrate W_cs) *= integrin_strong_factor (3.0, strengthen)
    flag cell as 'switched'
```

Switched RIM cells (high pressure that is RELEASED at the free surface — they have
open space + are now strongly substrate-adherent) spread outward, the physiological
unjamming-and-migrate behaviour. Switched CORE cells are also necrotic (depth) so
their switch is moot — the visible junction-switch population is the high-pressure
rim ring, which is exactly the spreading front.
