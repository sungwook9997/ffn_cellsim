#!/bin/bash
# Stage figures for Claude Design: curated deck set + movie posters/gifs + thumbnails for browse gallery.
set -u
cd /Users/sw1/ffn_cellsim/ffn_sim/outputs || exit 1

DA=deck_assets
mkdir -p "$DA/curated" "$DA/_thumbs"

echo "=== 1) Copy curated still figures (slide-numbered) ==="
# format: SRC|DSTNAME
printf '%s\n' \
"warp_decohesion/figs/ipc_sphere_clean_n1000_mesh.png|slide01_hero_spheroid_n1138.png" \
"layer2/figs/fig_layer2_pi_overlay.png|slide02_pi_experiment_overlay.png" \
"h3/figs/fig_h3_evidence_montage.png|slide07-08_cortex_validation_montage.png" \
"warp_decohesion/figs/young_dupre_doublet_gate.png|slide08_young_dupre_doublet.png" \
"layer2/figs/fig_layer2_prod_rlaw.png|slide09_rlaw_r2_0998.png" \
"layer2/figs/fig_layer2_ligand_production.png|slide09_ligand_ordering.png" \
"layer2/figs/fig_layer2_magnitude_gap.png|slide10_magnitude_gap.png" \
"layer2/figs/fig_layer2_decisive_traction.png|slide10_decisive_traction_ejection.png" \
"h_dcm_two_stage/figs/two_stage_n400.png|slide12_two_stage_storyboard_n400.png" \
"warp_decohesion/figs/cohesion_fix_peeling_vs_cohesive.png|slide13_retraction_peeling_vs_cohesive.png" \
"warp_decohesion/figs/A1_actual_cells_3d.png|slide13_actual_cells_7mobile_93frozen.png" \
"warp_decohesion/figs/faceting_comparison.png|slide14_faceting_comparison.png" \
"warp_decohesion/figs/RECONSIDER_xsec_surfcore.png|slide14_xsection_surf_core.png" \
"warp_decohesion/figs/ipc_cleave_opt_n1000_mesh.png|slide16_engine_scale_n1128.png" \
"warp_decohesion/figs/CLEAVE_mechanism.png|backup_cleave_mechanism.png" \
"h_dcm_two_stage/figs/a2_decohesion_n100_montage.png|backup_decohesion_montage_n100.png" \
"h_dcm_gpu_lod/figs/necrosis_3zone_spatial.png|backup_necrosis_3zone.png" \
| while IFS='|' read -r src dst; do
    if [ -f "$src" ]; then cp -f "$src" "$DA/curated/$dst"; echo "  + $dst"; else echo "  ! MISSING $src"; fi
  done

echo "=== 2) Movie -> poster PNG + compressed GIF ==="
# format: SRC|BASENAME
printf '%s\n' \
"h_dcm_two_stage/figs/two_stage_n400_agg.mp4|movie_aggregation_n400" \
"h_dcm_two_stage/figs/two_stage_n800_spread.mp4|movie_spread_compaction_n800" \
"h_dcm_gpu_lod/figs/junction_switch_spatial.mp4|movie_decohesion_junction_switch" \
"h_dcm_gpu_lod/figs/necrosis_3zone_spatial.mp4|movie_necrosis_3zone" \
| while IFS='|' read -r src base; do
    if [ ! -f "$src" ]; then echo "  ! MISSING $src"; continue; fi
    # poster: representative frame (last ~1s = end-state)
    ffmpeg -y -sseof -1.0 -i "$src" -frames:v 1 -update 1 "$DA/curated/${base}_poster.png" >/dev/null 2>&1 \
      || ffmpeg -y -i "$src" -vf "thumbnail" -frames:v 1 -update 1 "$DA/curated/${base}_poster.png" >/dev/null 2>&1
    # gif: 12fps, 720px wide, palette for quality
    ffmpeg -y -i "$src" -vf "fps=12,scale=720:-1:flags=lanczos,palettegen" "$DA/_thumbs/_pal_${base}.png" >/dev/null 2>&1
    ffmpeg -y -i "$src" -i "$DA/_thumbs/_pal_${base}.png" -lavfi "fps=12,scale=720:-1:flags=lanczos [x]; [x][1:v] paletteuse" "$DA/curated/${base}.gif" >/dev/null 2>&1
    rm -f "$DA/_thumbs/_pal_${base}.png"
    echo "  + ${base}_poster.png + ${base}.gif"
  done

echo "=== 3) Thumbnails for ALL pngs (parallel sips, 360px) ==="
find . -path "./$DA" -prune -o -name '*.png' -print | sed 's|^\./||' \
  | xargs -P 8 -I {} sh -c '
      src="$1"; flat=$(printf "%s" "$src" | tr "/ " "__");
      sips -Z 360 "$src" --out "deck_assets/_thumbs/$flat" >/dev/null 2>&1
    ' _ {}
echo "  png thumbs: $(find $DA/_thumbs -name '*.png' | wc -l | tr -d ' ')"

echo "=== 4) Posters/thumbnails for ALL movies (parallel ffmpeg, 360px) ==="
find . -path "./$DA" -prune -o -name '*.mp4' -print | sed 's|^\./||' \
  | xargs -P 4 -I {} sh -c '
      src="$1"; flat=$(printf "%s" "$src" | tr "/ " "__"); out="deck_assets/_thumbs/${flat%.mp4}__MOV.png";
      ffmpeg -y -i "$src" -vf "thumbnail,scale=360:-1" -frames:v 1 -update 1 "$out" >/dev/null 2>&1
    ' _ {}
echo "  mp4 posters: $(find $DA/_thumbs -name '*__MOV.png' | wc -l | tr -d ' ')"

echo "=== 5) Build browse gallery HTML ==="
python3 /Users/sw1/ffn_cellsim/ffn_sim/outputs/deck_assets/build_gallery.py

echo "=== DONE ==="
du -sh "$DA" 2>/dev/null
echo "curated files: $(find $DA/curated -type f | wc -l | tr -d ' ')"
