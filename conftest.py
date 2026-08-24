"""Repo-root pytest configuration — machine-enforced import invariants.

The 2026-07-25 adversarial audit found that two of this project's HARD rules were
enforced only by convention, and that one of them was already silently violated on
the GATE-A path. Both are now enforced by the interpreter, at conftest import time,
before a single test module is collected.

--------------------------------------------------------------------------------
1. ONE module identity per source file  (audit finding R4)
--------------------------------------------------------------------------------
``aleph/`` is importable under two names whenever both the repo root *and*
``aleph/`` itself are on ``sys.path`` — which is exactly the documented gbook
runtime contract (``PYTHONPATH=~/ffn_ac_native:~/ffn_ac_native/ffn_sim`` as it was
written on the retired gbook host — the path is quoted as history, not as a target; kept that
way to defeat the stale-editable-install shadow). Under that contract:

    from ac.cell.assemble import CellConfig        # 22 scripts, incl. the GATE-A
                                                   # ac_resting_converge.py
    from aleph.components.incumbent.compartments import ...   # what assemble.py itself does

loads the SAME FILES TWICE as two unrelated module objects. Measured consequences,
not hypothetical ones: duplicated module-level state (caches, registries, RNG),
duplicated ``@wp.kernel`` registration, and ``isinstance`` / dataclass identity
failing across the split (``isinstance(ac...CORTEX_REGION, aleph.ac...RegionSpec)``
is ``False``). A bug of that shape presents as "the parameter I set did not take
effect", which is indistinguishable from a physics result.

``_CanonicalIdentityFinder`` makes ``aleph.<pkg>`` the single canonical identity
and resolves every bare ``<pkg>...`` import to the *same object*, so the two spellings
are aliases rather than copies. ``aleph.*`` is canonical (not the reverse) because
that is what the library's own internal absolute imports use.

Scope of the alias set is deliberately narrow and evidence-based: only the four
sibling packages actually imported bare somewhere in ``aleph/{scripts,tests}``
(``ac``, ``ff``) plus their two peers that could be next (``dcm``, ``common``).
``archive`` (HOOMD port spec) and ``validation`` (oracles, runtime-import-forbidden)
are deliberately NOT aliased — making them easier to import top-level would work
against two other HARD rules.

--------------------------------------------------------------------------------
2. No HOOMD inside a pytest process  (I0-A contract)
--------------------------------------------------------------------------------
"Warp-CUDA is the only simulation runtime; HOOMD is never imported or executed" was
enforced by (a) a static AST scan of ``ac/`` only and (b) HOOMD not happening to be
installed. Neither is enforcement: the CI env used to *create* ``hoomd=7.0.1`` and
pass ``aleph/tests`` explicitly, so ~35 HOOMD-importing legacy test modules ran on
every push. ``_HoomdImportBlocker`` refuses the import in-process regardless of what
is installed in the environment.

It raises an ``ImportError`` subclass on purpose: ``except ImportError`` sites and
``pytest.importorskip("hoomd")`` then degrade to an honest SKIP instead of an
opaque crash, while an unguarded ``import hoomd`` still fails loudly with the reason.

--------------------------------------------------------------------------------
3. No stale-editable shadow
--------------------------------------------------------------------------------
Also checked: the ``aleph`` that actually imports is the one in THIS working tree.
A stale editable/site-packages install silently testing other code is the
"editable-shadow" trap that already cost this project a debugging session.

Failure mode of this file is loud and total: any violated invariant raises during
conftest import, which aborts the whole pytest session. There is no env-var escape
hatch — an escape hatch on an import-identity invariant is gate-loosening.
"""

from __future__ import annotations

import importlib
import importlib.abc
import importlib.machinery
import importlib.util
import sys
from pathlib import Path
from types import ModuleType

_REPO_ROOT = Path(__file__).resolve().parent
_CANONICAL_ROOT = "aleph"

# Sibling packages of aleph/ that a bare `import <name>` must resolve to the SAME
# object as `import aleph.<name>`. The bare-import habit came from 22 scripts that
# double-loaded the package; the aliases keep a stray one from producing a second copy.
#
# Updated by the 2026-08-09 restructure: `ac` split into `engine` + `components`, `ff`
# became `laws`, and `dcm` moved under `archive`. NOT aliased on purpose: `archive`
# (nothing live may reach it) and `validation` (oracles must never sit on a runtime
# import path — a runtime that can import its own oracle can be made to agree with it).
_ALIASED_PACKAGES = ("engine", "components", "laws")


# --------------------------------------------------------------------------- #
# 1. one canonical module identity
# --------------------------------------------------------------------------- #
class _AliasLoader(importlib.abc.Loader):
    """Loader that re-publishes an already-imported canonical module under an alias.

    ``create_module`` hands back the *existing* module object, so the alias name and
    the canonical name map to one identity in ``sys.modules``. ``exec_module`` then
    restores the canonical ``__name__`` / ``__spec__`` / ``__package__`` that the
    import machinery re-stamps with the alias identity on the way through, so the
    single shared module keeps telling the truth about where it came from.
    """

    def __init__(self, canonical_name: str, canonical_spec: importlib.machinery.ModuleSpec | None,
                 canonical_package: str | None) -> None:
        self._canonical_name = canonical_name
        self._canonical_spec = canonical_spec
        self._canonical_package = canonical_package

    def create_module(self, spec: importlib.machinery.ModuleSpec) -> ModuleType:
        return importlib.import_module(self._canonical_name)

    def exec_module(self, module: ModuleType) -> None:
        # Already executed under its canonical name — re-executing would duplicate
        # module-level side effects (the very thing this finder exists to prevent).
        module.__name__ = self._canonical_name
        if self._canonical_spec is not None:
            module.__spec__ = self._canonical_spec
            module.__loader__ = self._canonical_spec.loader
        module.__package__ = self._canonical_package


class _CanonicalIdentityFinder(importlib.abc.MetaPathFinder):
    """Resolve ``<alias>.a.b`` to the identical object as ``aleph.<alias>.a.b``."""

    def __init__(self, aliases: tuple[str, ...], canonical_root: str) -> None:
        self._aliases = frozenset(aliases)
        self._canonical_root = canonical_root

    def find_spec(
        self,
        fullname: str,
        path: object = None,
        target: ModuleType | None = None,
    ) -> importlib.machinery.ModuleSpec | None:
        if fullname.partition(".")[0] not in self._aliases:
            return None
        canonical = f"{self._canonical_root}.{fullname}"

        # Existence probe first: a genuinely missing module must fall through to the
        # normal machinery (so the error names the module the user asked for), while a
        # real import error INSIDE an existing canonical module must propagate rather
        # than be swallowed into "no such module".
        try:
            if importlib.util.find_spec(canonical) is None:
                return None
        except (ImportError, AttributeError, ValueError):
            return None

        module = importlib.import_module(canonical)
        loader = _AliasLoader(canonical, getattr(module, "__spec__", None),
                             getattr(module, "__package__", None))
        spec = importlib.machinery.ModuleSpec(
            fullname, loader, origin=getattr(module, "__file__", None)
        )
        search = getattr(module, "__path__", None)
        if search is not None:
            spec.submodule_search_locations = list(search)
        return spec


def _install_canonical_identity_finder() -> None:
    """Install the alias finder at the FRONT of ``sys.meta_path`` (before PathFinder)."""
    for existing in sys.meta_path:
        if isinstance(existing, _CanonicalIdentityFinder):
            return
    # A duplicate already resident in sys.modules cannot be repaired by aliasing:
    # references to the wrong copy are already held. Refuse rather than pretend.
    for alias in _ALIASED_PACKAGES:
        stale = sys.modules.get(alias)
        if stale is None:
            continue
        canonical = sys.modules.get(f"{_CANONICAL_ROOT}.{alias}")
        if canonical is not None and stale is not canonical:
            raise RuntimeError(
                f"module-identity split already established before conftest.py ran: "
                f"'{alias}' and '{_CANONICAL_ROOT}.{alias}' are two distinct module "
                f"objects for {getattr(stale, '__file__', '?')}. Something imported "
                f"'{alias}' before pytest loaded the root conftest (a plugin, a "
                f"sitecustomize, or -p). Import via '{_CANONICAL_ROOT}.{alias}'."
            )
    sys.meta_path.insert(0, _CanonicalIdentityFinder(_ALIASED_PACKAGES, _CANONICAL_ROOT))


def _assert_single_module_identity() -> None:
    """Self-check with teeth: prove the alias actually yields ONE object, per package.

    If this ever regresses (a Python import-machinery change, a competing meta-path
    finder, a stale sys.modules entry) the whole pytest session aborts here. A guard
    that cannot fail is worse than no guard, so this runs unconditionally.
    """
    for alias in _ALIASED_PACKAGES:
        canonical_name = f"{_CANONICAL_ROOT}.{alias}"
        canonical = importlib.import_module(canonical_name)
        aliased = importlib.import_module(alias)
        if aliased is not canonical:
            raise RuntimeError(
                f"module identity NOT canonical: '{alias}' is {id(aliased):#x} but "
                f"'{canonical_name}' is {id(canonical):#x} — the same source file is "
                f"loaded twice in one process (duplicated state + duplicated Warp "
                f"kernel registration + cross-identity isinstance failures)."
            )


def _assert_no_duplicate_loaded_files() -> None:
    """No file under this repo may be resident in ``sys.modules`` as 2+ objects."""
    by_path: dict[str, dict[int, list[str]]] = {}
    for name, module in list(sys.modules.items()):
        origin = getattr(module, "__file__", None)
        if not origin:
            continue
        try:
            resolved = Path(origin).resolve()
        except OSError:
            continue
        if _REPO_ROOT not in resolved.parents:
            continue
        by_path.setdefault(str(resolved), {}).setdefault(id(module), []).append(name)
    duplicated = {p: ids for p, ids in by_path.items() if len(ids) > 1}
    if duplicated:
        detail = "; ".join(
            f"{p} loaded as {sorted(n for names in ids.values() for n in names)}"
            for p, ids in sorted(duplicated.items())
        )
        raise RuntimeError(f"dual module identity detected: {detail}")


# --------------------------------------------------------------------------- #
# 2. HOOMD is not importable inside a pytest process (I0-A)
# --------------------------------------------------------------------------- #
class HoomdImportForbidden(ImportError):
    """Raised when anything in a pytest process tries to import HOOMD."""


class _HoomdImportBlocker(importlib.abc.MetaPathFinder):
    """Refuse ``import hoomd`` regardless of what is installed in the environment."""

    def find_spec(
        self,
        fullname: str,
        path: object = None,
        target: ModuleType | None = None,
    ) -> None:
        if fullname == "hoomd" or fullname.startswith("hoomd."):
            raise HoomdImportForbidden(
                f"import of '{fullname}' is forbidden: Warp-CUDA is the only simulation "
                "runtime (I0-A, PI 2026-07-16). The retired HOOMD tree was DELETED on "
                "2026-07-29 (PI: 'HOOMD 관련 내용은 전부 삭제'); it is recoverable from git "
                "history at 1a9ded66^ if a port needs re-deriving. Nothing in this repo "
                "may import it.",
                name=fullname,
            )
        return None


def _install_hoomd_blocker() -> None:
    for existing in sys.meta_path:
        if isinstance(existing, _HoomdImportBlocker):
            return
    sys.meta_path.insert(0, _HoomdImportBlocker())


# --------------------------------------------------------------------------- #
# 3. the aleph under test is THIS working tree
# --------------------------------------------------------------------------- #
def _assert_no_editable_shadow() -> None:
    if str(_REPO_ROOT) not in sys.path:
        sys.path.insert(0, str(_REPO_ROOT))
    package = importlib.import_module(_CANONICAL_ROOT)
    origin = getattr(package, "__file__", None)
    expected = _REPO_ROOT / _CANONICAL_ROOT
    if origin is None or Path(origin).resolve().parent != expected:
        raise RuntimeError(
            f"'{_CANONICAL_ROOT}' resolves to {origin} but this working tree is "
            f"{expected} — a stale editable/site-packages install is shadowing the "
            "repo, so the tests would validate code that is not the code you edited. "
            "Uninstall the shadow (pip uninstall ffn_cellsim) or fix PYTHONPATH."
        )


# --------------------------------------------------------------------------- #
# import-time enforcement (runs before any test module is collected)
# --------------------------------------------------------------------------- #
_install_hoomd_blocker()
_install_canonical_identity_finder()
_assert_no_editable_shadow()
_assert_single_module_identity()
_assert_no_duplicate_loaded_files()
