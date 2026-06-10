# Necrosis / 3-zone spheroid — literature-calibrated values

Tumour spheroids develop a **layered radial architecture** governed by the
diffusion-limited penetration of oxygen and nutrients from the free surface
inward. From the surface inward the canonical three zones are:

| Zone | Depth band below surface | State |
|------|--------------------------|-------|
| Proliferating rim | 0 – ~40 µm (outer 1–3 cell layers) | actively dividing (O₂ + nutrient replete) |
| Quiescent shell | ~40 – ~150 µm | viable but G0 arrested (hypoxic, no division) |
| Necrotic core | > ~150 µm | dead (anoxic; below the O₂ penetration depth) |

## Calibration constants (trusted literature bands)

- **Oxygen penetration depth ≈ 150–200 µm.** The viable rim thickness of tumour
  spheroids is repeatedly reported at ~100–200 µm, set by the O₂ diffusion length
  L = sqrt(D·C0 / k) (D ≈ 2e-9 m²/s, cellular O₂ consumption). We use
  **d_necrotic = 150 µm** as the depth at which cells cross into anoxia → necrosis.
  (Grimes et al. 2014 J R Soc Interface; Mueller-Klieser 1987 reviews; Sutherland 1988 Science.)
- **Proliferating rim ≈ outer 40 µm.** Ki-67⁺ proliferation is confined to the
  outermost ~2–3 cell diameters (~30–50 µm). We use **d_prolif = 40 µm**.
- **Necrotic-core onset diameter ≈ 400–500 µm.** A spheroid develops a visible
  necrotic core only once its diameter exceeds ~2× the viable-rim thickness,
  i.e. ~400 µm; below ~300 µm spheroids are fully viable. A 5-day MCF7/HCT116
  spheroid is typically 300–600 µm diameter.
- **Timing.** The necrotic core forms gradually over the **days-long aggregation**
  (cells compact, the centre crosses the penetration depth, then dies) — it is
  already present before any substrate-spreading assay begins. We therefore model
  an AGGREGATION phase (ball compacts + ages, core crosses into necrosis by depth)
  BEFORE the spreading phase.

## Depth thresholds used by NecrosisUpdater

```
d = R_cluster - r_cell          # depth below the spheroid surface (centroid radius)
d <  40 µm   -> PROLIFERATING   (rim; may divide)
40 <= d < 150 µm -> QUIESCENT   (viable, no division)
d >= 150 µm  -> NECROTIC        (dead; division off, marked grey/black)
```

## Coarse-graining note (REQUIRED for visibility, documented)

A real >400 µm spheroid resolved with 7.5-µm fine DCM cells needs hundreds of
cells (each cell = ~42 nodes) → thousands–tens-of-thousands of particles, which is
a GPU job. To make the 3-zone split VISIBLE in a ~5–10 min CPU demo we COARSE-GRAIN:
each DCM shell represents a **tissue patch**, the cell radius is scaled up so a
~40–60-cell ball spans the necrotic-onset scale (cluster radius ~150–250 µm,
diameter ~300–500 µm). The depth thresholds (40 / 150 µm) are kept at their
literature values; with R_patch ≈ 30 µm a ~50-cell ball reaches R_cluster ≈ 180 µm
so the inner cells cross the 150-µm necrotic threshold and an outer rim stays
proliferating — the 3 zones appear. This scaling is stated in every figure caption.
