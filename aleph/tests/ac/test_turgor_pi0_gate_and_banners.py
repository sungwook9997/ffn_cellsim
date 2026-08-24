"""Π₀ gate (PI 2026-07-25) + retired-cortex-build banner presence — the executable half of the
2026-07-25 claim-hygiene pass.

Two independent things are asserted, and each one is written so a REGRESSION MAKES IT FAIL (the point
of the whole pass — a gate a broken build cannot fail is worse than no gate):

A. **Π₀ is REQUIRED with no silent default.** ``aleph/laws/params_turgor.yaml`` holds
   ``value: null`` / ``evidence_status: "GAP — PI"``; ``resolve_pi0_pa()`` with no value raises
   :class:`TurgorPi0Unspecified` naming all three competing in-tree values and pointing at
   ``PARAM_PROVENANCE_AUDIT_2026-07-24.md``. Passing 40.0 still runs, and stamps
   ``provenance="CONVENIENCE"`` on the artifact record. The physics entry points that consume Π₀
   (``ff.gamma_floor.turgor_pressure`` / ``measure_gamma``) are exercised through the gate too, so
   re-introducing a default in EITHER the ledger or a function signature fails here.

B. **The retired-cortex-build banner is present** in the docs/artifacts that carry the γ-floor
   headline, and the three still-live 2026-07-23 structural defects on
   ``gamma_floor.build_crosslinked_cortex`` are asserted to be exactly as reported (so the banner
   cannot be quietly removed while the defects remain, and cannot be silently *kept* after they are
   fixed — the "defects are still live" assertions fail once the builder is fixed, which is the
   signal to re-run and drop the banner).

Runs on the dev Mac (no CUDA): pure host-side parameter/parse checks, no Warp kernels, no HOOMD.
"""

from __future__ import annotations

import ast
import inspect
import re
import textwrap
from pathlib import Path

import numpy as np
import pytest

from aleph.laws.turgor_pi0 import (
    AUDIT_DOC,
    PI0_CLAIM_BAND_IMPLIED_133PA,
    PI0_CLAIM_HELA_PROXY_40PA,
    PI0_CLAIM_MCF7_GEOMETRY_72PA,
    TurgorPi0Invalid,
    TurgorPi0Unspecified,
    convenience_label,
    ledger_gated_value,
    resolve_pi0_pa,
)

REPO = Path(__file__).resolve().parents[3]


# ─────────────────────────────── A. the Π₀ gate ───────────────────────────────

def test_ledger_holds_no_value_so_there_is_nothing_to_default_to():
    """The gate's foundation: the ledger's Π₀ is `null`. If someone fills it in, this fails LOUDLY —
    a value in the ledger would silently become a default again for every caller."""
    assert ledger_gated_value() is None, (
        "params_turgor.yaml::parameters.PI_0.value must stay null (GAP — PI). A value here restores "
        "the silent default the 2026-07-25 gate removed; Π₀ selection is a PI decision."
    )


def test_missing_pi0_raises_with_all_three_competing_values_and_the_audit():
    """NEGATIVE CONTROL: no Π₀ ⇒ refuse, and say why in a way a reader can act on."""
    with pytest.raises(TurgorPi0Unspecified) as ei:
        resolve_pi0_pa(caller="test")
    msg = str(ei.value)
    for pa in (PI0_CLAIM_HELA_PROXY_40PA, PI0_CLAIM_BAND_IMPLIED_133PA, PI0_CLAIM_MCF7_GEOMETRY_72PA):
        assert f"{pa:g} Pa" in msg, f"refusal must name the {pa:g} Pa claim"
    assert "params_turgor.yaml" in msg
    assert "PARAM_PROVENANCE_AUDIT_2026-07-24" in msg
    assert AUDIT_DOC in msg


def test_explicit_40_pa_runs_and_is_recorded_as_convenience():
    """POSITIVE CONTROL: an explicit 40.0 still runs — and the artifact says CONVENIENCE, not sourced."""
    r = resolve_pi0_pa(40.0, caller="test")
    assert r.value_pa == 40.0
    assert r.provenance == "CONVENIENCE"
    assert r.claim_id == "claim_hela_proxy_40pa"
    assert r.is_sourced_for_mcf7 is False
    rec = r.as_artifact_record()
    assert rec["turgor_PI_0_Pa"] == 40.0
    assert rec["turgor_PI_0_provenance"] == "CONVENIENCE"
    assert rec["turgor_PI_0_sourced_for_MCF7"] is False
    assert "HeLa" in rec["turgor_PI_0_source"]


def test_133_pa_is_recorded_as_circular_band_implied_not_sourced():
    r = resolve_pi0_pa(133.0, caller="test")
    assert r.provenance == "PI_RATIFIED_BAND_IMPLIED"
    assert r.is_sourced_for_mcf7 is False


def test_pi0_sign_and_finiteness_sanity_gate():
    """Sanity Gate: Π₀ ≥ 0 (positive = net OUTWARD); Π₀=0 is accepted but flagged unphysical."""
    with pytest.raises(TurgorPi0Invalid):
        resolve_pi0_pa(-1.0, caller="test")
    with pytest.raises(TurgorPi0Invalid):
        resolve_pi0_pa(float("nan"), caller="test")
    zero = resolve_pi0_pa(0.0, caller="test")
    assert zero.provenance == "UNPHYSICAL_ZERO_BASELINE"


def test_unregistered_value_is_accepted_but_flagged():
    """A sweep point / new selection is not blocked, but it is not laundered as sourced either."""
    r = resolve_pi0_pa(57.0, caller="test")
    assert r.provenance == "UNREGISTERED"
    assert r.is_sourced_for_mcf7 is False


def test_gamma_floor_turgor_entry_points_have_no_silent_pi0_default():
    """The physics entry points must not carry a numeric Π₀ default in their signature."""
    from aleph.laws import gamma_floor as gf

    for fn in (gf.turgor_pressure, gf.external_force, gf.measure_gamma,
               gf.gamma_floor_run, gf.gamma_floor_production, gf.equilibrate):
        default = inspect.signature(fn).parameters["dP0"].default
        assert default is None, (
            f"{fn.__name__} re-introduced a Π₀ default ({default!r}); Π₀ is gated (PI 2026-07-25)."
        )


def test_turgor_pressure_refuses_without_pi0_and_runs_with_it():
    """End-to-end through the real physics path: negative then positive control."""
    from aleph.laws.gamma_floor import (TURGOR_DP0, CortexParams, build_crosslinked_cortex,
                                        turgor_pressure)

    cx = build_crosslinked_cortex(CortexParams(), n_filaments=30, n_xl=30, n_myo=3,
                                  rng=np.random.default_rng(0))
    with pytest.raises(TurgorPi0Unspecified):
        turgor_pressure(cx, cx.net.pos)                    # no Π₀ → refuse
    dP, R_mean = turgor_pressure(cx, cx.net.pos, dP0=TURGOR_DP0)   # explicit 40 Pa → run
    assert np.isfinite(dP) and dP > 0.0 and R_mean > 0.0


def test_measure_gamma_requires_pi0_only_when_the_passive_channel_is_measured():
    """turgor=False (the γ_active-only channel) must stay usable without Π₀ — the gate is scoped."""
    from aleph.laws.gamma_floor import (TURGOR_DP0, build_crosslinked_cortex, equilibrate,
                                        measure_gamma)

    cx = build_crosslinked_cortex(n_filaments=30, n_xl=90, n_myo=20, rng=np.random.default_rng(3))
    equilibrate(cx, 0.0, n_steps=40, turgor=False)
    off = measure_gamma(cx, 5.0, turgor=False)              # no Π₀ needed
    assert off["gamma_passive"] == 0.0
    assert off["turgor_PI_0_provenance"] == "TURGOR_OFF"
    with pytest.raises(TurgorPi0Unspecified):
        measure_gamma(cx, 5.0, turgor=True)                # passive channel without Π₀ → refuse
    on = measure_gamma(cx, 5.0, turgor=True, dP0=TURGOR_DP0)
    assert on["gamma_passive"] == pytest.approx(0.5 * TURGOR_DP0 * cx.R_um)
    assert on["turgor_PI_0_provenance"] == "CONVENIENCE"    # provenance travels with the number


def test_pi0_units_identity_pa_equals_pn_per_um2():
    """Sanity Gate (dimensional): 1 Pa == 1 pN/µm², so Π₀ needs no conversion between conventions."""
    from aleph.components.incumbent import assemble
    from aleph.laws.gamma_floor import TURGOR_DP0

    assert assemble.PI_0_PA == TURGOR_DP0
    assert assemble.PI_0_PROVENANCE["turgor_PI_0_provenance"] == "CONVENIENCE"


def test_ac_cell_convenience_labels_are_registered_and_match_the_code():
    """The three labeled-not-gated params: ledger value must equal the live code default."""
    from aleph.components.incumbent.assemble import CellConfig
    from aleph.laws.architecture_spec import CORTEX

    cfg = CellConfig()
    assert convenience_label("cortex_seg_um")["value"] == cfg.cortex_seg_um == CORTEX.filament.seg_um
    assert convenience_label("cortex_length_um")["value"] == cfg.cortex_length_um == CORTEX.filament.length_um
    assert (convenience_label("cortex_density_per_fil")["value"]
            == cfg.cortex_density_per_fil == CORTEX.crosslinker.density_per_fil)
    for name in ("cortex_seg_um", "cortex_length_um", "cortex_density_per_fil"):
        row = convenience_label(name)
        assert row["provenance_class"] == "CONVENIENCE"
        assert row["gated"] is False                       # PI 2026-07-25: labeled, NOT gated
        assert row["deviation_factor"]                     # a non-empty statement of the deviation
        assert row["owner"] == "PI" and row["expiry_trigger"]


# ─────────────── B. the retired-cortex-build banner + the still-live defects ───────────────

BANNER_KEY = "pre-2026-07-23 cortex build"
BANNERED_DOCS = [
    "aleph/laws/ENGINE.md",
    "aleph/dcm/ENGINE.md",
    "aleph/docs/v2_audit/FF_STAGE6O_NATIVE_GAMMA_2026-07-01.md",
    "aleph/outputs/_archive/mech_hier/REPORT.md",
    "aleph/outputs/_archive/engine_reval/REPORT.md",
]


@pytest.mark.parametrize("rel", BANNERED_DOCS)
def test_retired_cortex_build_banner_present(rel: str):
    """Every doc carrying a full-native ff/ headline must say the cortex it was measured on is retired."""
    text = (REPO / rel).read_text(encoding="utf-8", errors="ignore")
    assert BANNER_KEY in text, f"{rel} is missing the retired-cortex-build banner"
    assert "CORTEX_STRUCTURE_AUDIT_2026-07-23" in text, f"{rel} banner must cite the audit"
    assert re.search(r"not (?:been )?re-run", text), f"{rel} banner must state it has not been re-run"


def test_three_2026_07_23_defects_are_still_live_in_the_gamma_floor_builder():
    """The banner's factual basis, asserted against the CODE (not a doc claim).

    When these assertions start FAILING, the defects have been fixed — that is the signal to re-run
    the γ-floor on the fixed cortex and then remove the banner. Until then the banner stands.
    """
    from aleph.laws import cortex_assembly, gamma_floor

    fn = gamma_floor.build_crosslinked_cortex
    sig = inspect.signature(fn)
    # Inspect the executable BODY only — the docstring itself names these tokens (it carries the
    # banner), so a naive whole-source grep would self-trip and the gate would be vacuous. AST-strip
    # the docstring rather than string-replacing it (getdoc() re-indents and would not match).
    tree = ast.parse(textwrap.dedent(inspect.getsource(fn)))
    fdef = tree.body[0]
    assert isinstance(fdef, ast.FunctionDef)
    stmts = fdef.body[1:] if ast.get_docstring(fdef) is not None else fdef.body
    body = "\n".join(ast.unparse(s) for s in stmts)

    # (1) n_xl is a caller argument and the spec is never consulted.
    assert "n_xl" in sig.parameters, "n_xl is still a caller argument, not read from architecture_spec"
    assert "architecture_spec" not in body and "CORTEX" not in body, (
        "build_crosslinked_cortex now reads the architecture spec — re-run the γ-floor and drop the banner"
    )
    # ... and the production preset still passes density 1.0 (n_xl == n_filaments).
    assert gamma_floor.PROD_N_XL == gamma_floor.PROD_N_FIL, (
        "PROD_N_XL != PROD_N_FIL — crosslink density changed; re-run the γ-floor and drop the banner"
    )

    # (2) resolve_overlaps is never passed, so build_cortex_network runs at its default False.
    assert "resolve_overlaps" not in body, "resolve_overlaps is now wired — re-run and drop the banner"
    assert (inspect.signature(cortex_assembly.build_cortex_network)
            .parameters["resolve_overlaps"].default is False)

    # (3) length_dist defaults to "mono".
    assert sig.parameters["length_dist"].default == "mono"

    # and the banner text the module ships states all three.
    banner = gamma_floor._CORTEX_BUILD_BANNER
    assert BANNER_KEY in banner
    for token in ("n_xl", "resolve_overlaps", "length_dist", "CORTEX_STRUCTURE_AUDIT_2026-07-23"):
        assert token in banner


def test_engine_md_going_forward_claims_are_qualified():
    """ff/ENGINE.md and dcm/ENGINE.md both called themselves "the going-forward engine"; the canonical
    composition layer is `engine/` (PI 2026-07-22; `ac/engine/` until the 2026-08-09 rename).
    The banner must qualify the claim."""
    for rel in ("aleph/laws/ENGINE.md", "aleph/dcm/ENGINE.md"):
        text = (REPO / rel).read_text(encoding="utf-8", errors="ignore")
        assert "engine/" in text, f"{rel} must name engine/ as the canonical composition layer"
