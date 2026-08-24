"""Derived render buffers, and the registry that is the only way to get one.

Manuscript §11: *"Viewer receives only derived render buffers; no authoritative array or writable IPC
handle is exposed."* That sentence can be honoured by convention — everybody agrees not to hand the
viewer the real array — and a convention is exactly what fails once the code is large enough that
nobody reads all of it. So it is honoured here by construction, in three layers:

1. :class:`DerivedBuffer` refuses to hold a writeable array and refuses to hold a *view*. An array
   with a non-``None`` ``base`` is a window onto somebody else's memory, and freezing the window does
   not freeze what it looks at: the owner can still write through its own handle, and the viewer's
   "immutable" buffer changes under it. So a view is refused even when it is read-only.
2. :class:`ObservationRegistry` keeps the authoritative arrays in a name-mangled private attribute and
   exposes only metadata about them — owners, array names, shapes, dtypes. The derivation functions
   registered against it *do* receive the authoritative arrays, because that is where the reduction
   has to run; what the registry guarantees is that nothing the viewer can reach came back out.
3. :func:`assert_no_authoritative_reach` walks the registry's public surface and every buffer it has
   produced, and fails if any reachable object shares memory with an authoritative array or is a
   writeable array. It is the guard that closes the gap the first two layers cannot see: a derivation
   that returns a slice, an accessor added later that forwards the wrong object.

The distinction between "read-only" and "derived" is worth stating once. A read-only view of
authoritative memory satisfies "the viewer cannot write", and it still breaks the contract, because
the viewer is then reading physics that is being mutated underneath it mid-step and will report torn
state as an observation. Derived means *detached*: its own memory, its own lifetime, and a stated
derivation naming what reduction produced it.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Callable, Iterable, Mapping, Sequence

import numpy as np

__all__ = [
    "AuthoritativeReachError",
    "DerivationSpec",
    "DerivedBuffer",
    "ObservationRegistry",
    "assert_no_authoritative_reach",
    "freeze_derived",
]


class AuthoritativeReachError(RuntimeError):
    """Something reachable from the viewer side touches authoritative memory.

    Raised at the moment the reach is created rather than when it is used, so the traceback names the
    derivation or the accessor that leaked rather than the innocent consumer that read it.
    """


def freeze_derived(values: Any, *, dtype: Any = np.float64) -> np.ndarray:
    """Return a detached, read-only copy of ``values``.

    ``np.array(..., copy=True)`` rather than ``asarray``: the copy is the point. The returned array
    owns its memory (``base is None``), so nothing that happens to the source afterwards can be
    observed through it.

    Args:
        values: Anything array-like.
        dtype: Element type of the copy.

    Returns:
        A C-contiguous, non-writeable array with ``base is None``.
    """
    frozen = np.array(values, dtype=dtype, copy=True, order="C")
    frozen.setflags(write=False)
    return frozen


@dataclass(frozen=True, slots=True)
class DerivedBuffer:
    """One render buffer: detached values plus everything needed to read them honestly.

    Attributes:
        name: Registration name of the derivation that produced it.
        owner: The authoritative owner the derivation read from. Recorded so a scene node can be
            traced back to whose physics it depicts, without holding a handle on that physics.
        units: Physical units of the values, as text. Required — an unlabelled render buffer is a
            picture, and a picture is not a measurement.
        derivation: What reduction was applied, in words. ``"per-node force magnitude"`` is a
            different observable from ``"per-node force z-component"`` and the buffer must say which.
        values: The detached, read-only values.
        active_count: How many entries are live. A buffer sized for the maximum population and
            partly filled is normal; a viewer that colours the unfilled tail is reporting zeros as
            physics, so the count travels with the values.
        accepted_step: The accepted step this buffer was derived from.
        topology_epoch: The topology generation at that step.
    """

    name: str
    owner: str
    units: str
    derivation: str
    values: np.ndarray
    active_count: int
    accepted_step: int
    topology_epoch: int

    def __post_init__(self) -> None:
        for field_name in ("name", "owner", "units", "derivation"):
            if not str(getattr(self, field_name)).strip():
                raise ValueError(f"DerivedBuffer.{field_name} is required and must be non-empty")
        if not isinstance(self.values, np.ndarray):
            raise TypeError(
                f"DerivedBuffer.values must be an ndarray; got {type(self.values).__name__}"
            )
        if self.values.flags.writeable:
            raise AuthoritativeReachError(
                f"derived buffer {self.name!r} was handed a writeable array. A writeable buffer is a "
                "write path into whatever it aliases, and §11 gives the viewer no write access of "
                "any kind. Build it with freeze_derived()."
            )
        if self.values.base is not None:
            raise AuthoritativeReachError(
                f"derived buffer {self.name!r} was handed a view (base is not None), not a derived "
                "array. Freezing a view freezes the window and not the memory: the owner still "
                "writes through its own handle and the viewer reads torn mid-step state. Copy it."
            )
        if isinstance(self.active_count, bool) or not isinstance(self.active_count, int):
            raise TypeError(f"active_count must be an integer; got {self.active_count!r}")
        if not (0 <= self.active_count <= int(self.values.shape[0] if self.values.ndim else 1)):
            raise ValueError(
                f"derived buffer {self.name!r} declares active_count={self.active_count} over "
                f"{self.values.shape} values; an active count larger than the buffer is a count of "
                "entries that were never written"
            )
        for field_name in ("accepted_step", "topology_epoch"):
            value = getattr(self, field_name)
            if isinstance(value, bool) or not isinstance(value, int):
                raise TypeError(f"DerivedBuffer.{field_name} must be an integer; got {value!r}")
            if value < 0:
                raise ValueError(f"DerivedBuffer.{field_name} cannot be negative; got {value}")

    @property
    def active_values(self) -> np.ndarray:
        """The live prefix of the buffer, still read-only.

        This is a slice and therefore a view — of the *derived* buffer, which nothing authoritative
        aliases, so the hazard :class:`DerivedBuffer` refuses at construction does not apply here.
        """
        sliced = self.values[: self.active_count]
        sliced.setflags(write=False)
        return sliced

    @property
    def extrema(self) -> tuple[float, float]:
        """``(min, max)`` over the active prefix, as actually measured.

        Raises:
            ValueError: When nothing is active. An extremum over an empty population is not zero.
        """
        if self.active_count == 0:
            raise ValueError(
                f"derived buffer {self.name!r} has no active entries, so it has no extrema; "
                "reporting 0.0 here is how an empty population becomes a measured value"
            )
        active = self.active_values
        return (float(np.min(active)), float(np.max(active)))

    def as_json_obj(self) -> dict[str, Any]:
        """Return the buffer's metadata, without the values.

        The values are bulk data and belong in a frame payload; this is what goes in a manifest.
        """
        return {
            "name": self.name,
            "owner": self.owner,
            "units": self.units,
            "derivation": self.derivation,
            "shape": list(self.values.shape),
            "dtype": str(self.values.dtype),
            "active_count": self.active_count,
            "accepted_step": self.accepted_step,
            "topology_epoch": self.topology_epoch,
        }


@dataclass(frozen=True, slots=True)
class DerivationSpec:
    """A registered reduction from one owner's authoritative arrays to one render buffer.

    Attributes:
        name: Unique registration name.
        owner: Which authority the reduction reads.
        units: Units of the result.
        derivation: What the reduction is, in words.
        reduce: ``reduce(arrays) -> array_like``. Runs on the producer side, with the authoritative
            arrays, and must return something new. Returning an input is caught by the registry.
        active_count: ``active_count(arrays) -> int``, or ``None`` to mean "all of it".
    """

    name: str
    owner: str
    units: str
    derivation: str
    reduce: Callable[[Mapping[str, np.ndarray]], Any]
    active_count: Callable[[Mapping[str, np.ndarray]], int] | None = None

    def __post_init__(self) -> None:
        for field_name in ("name", "owner", "units", "derivation"):
            if not str(getattr(self, field_name)).strip():
                raise ValueError(f"DerivationSpec.{field_name} is required")
        if not callable(self.reduce):
            raise TypeError(f"DerivationSpec.reduce for {self.name!r} must be callable")


class ObservationRegistry:
    """The producer's side of the viewer boundary: the only source of derived buffers.

    The authoritative arrays go in through :meth:`declare_authority` and never come back out. What
    comes out of :meth:`derive` is a :class:`DerivedBuffer`, which by its own constructor cannot be a
    handle on the arrays that produced it.

    The private store is name-mangled (``self.__authority``) rather than single-underscore. That is
    not security — Python has none — it is a statement about intent that makes an accidental
    ``registry._authority["membrane"]["positions"]`` in viewer code fail rather than work.
    """

    def __init__(self) -> None:
        self.__authority: dict[str, dict[str, np.ndarray]] = {}
        self._specs: dict[str, DerivationSpec] = {}
        self._derive_count: int = 0

    # -- authority side ---------------------------------------------------------------------
    def declare_authority(self, owner: str, arrays: Mapping[str, np.ndarray]) -> None:
        """Register one owner's authoritative arrays.

        Args:
            owner: The owning participant's name.
            arrays: Its live authoritative arrays, by name. Held by reference on purpose — the whole
                point is that a derivation reads the *current* state, not a stale copy.

        Raises:
            ValueError: If the owner is already declared, or an array name is empty.
            TypeError: If any value is not an ndarray.
        """
        if not owner.strip():
            raise ValueError("an authority needs an owner name")
        if owner in self.__authority:
            raise ValueError(
                f"authority {owner!r} is already declared; re-declaring it would leave two array "
                "sets claiming to be the same owner's state"
            )
        held: dict[str, np.ndarray] = {}
        for name, array in arrays.items():
            if not str(name).strip():
                raise ValueError(f"authority {owner!r} declares an unnamed array")
            if not isinstance(array, np.ndarray):
                raise TypeError(
                    f"authority {owner!r} array {name!r} must be an ndarray; got "
                    f"{type(array).__name__}"
                )
            held[str(name)] = array
        self.__authority[owner] = held

    def register(self, spec: DerivationSpec) -> None:
        """Register a reduction.

        Args:
            spec: The derivation.

        Raises:
            ValueError: If the name is taken, or the owner is not declared. Registering against an
                undeclared owner is refused rather than deferred: a derivation whose authority
                appears later is a derivation nobody checked the array names of.
        """
        if not isinstance(spec, DerivationSpec):
            raise TypeError("register() takes a DerivationSpec")
        if spec.name in self._specs:
            raise ValueError(f"derivation {spec.name!r} is already registered")
        if spec.owner not in self.__authority:
            raise ValueError(
                f"derivation {spec.name!r} reads authority {spec.owner!r}, which is not declared; "
                f"declared authorities are {sorted(self.__authority)}"
            )
        self._specs[spec.name] = spec

    # -- metadata only ----------------------------------------------------------------------
    @property
    def authority_owners(self) -> tuple[str, ...]:
        """Names of the declared authorities. Names, not arrays."""
        return tuple(sorted(self.__authority))

    @property
    def derivation_names(self) -> tuple[str, ...]:
        """Names of the registered derivations."""
        return tuple(sorted(self._specs))

    @property
    def derive_count(self) -> int:
        """How many buffers this registry has produced. A cheap cost witness."""
        return self._derive_count

    def authority_shapes(self) -> dict[str, dict[str, tuple[Any, ...]]]:
        """Per-owner, per-array ``(shape..., dtype)`` metadata.

        Deliberately returns plain tuples of ints and a dtype *string*, never the arrays and never
        anything that holds a reference to them.
        """
        return {
            owner: {
                name: (*(int(n) for n in array.shape), str(array.dtype))
                for name, array in sorted(arrays.items())
            }
            for owner, arrays in sorted(self.__authority.items())
        }

    def spec(self, name: str) -> DerivationSpec:
        """The registered spec for ``name``.

        Returns the spec, which carries a callable that closes over nothing authoritative unless the
        caller made it so. The registry cannot police what a caller's closure captured, which is why
        :func:`assert_no_authoritative_reach` checks results rather than trusting registration.
        """
        try:
            return self._specs[name]
        except KeyError:
            raise KeyError(
                f"no derivation named {name!r}; registered: {list(self.derivation_names)}"
            ) from None

    # -- the boundary crossing --------------------------------------------------------------
    def derive(self, name: str, *, accepted_step: int, topology_epoch: int) -> DerivedBuffer:
        """Run one registered reduction and return its detached buffer.

        Args:
            name: Registered derivation name.
            accepted_step: The accepted step being observed.
            topology_epoch: The topology generation at that step.

        Returns:
            The derived buffer.

        Raises:
            AuthoritativeReachError: If the reduction returned an object that shares memory with any
                authoritative array. This is the guard that catches ``lambda arrays:
                arrays["positions"]`` and every slice of it — a reduction that forgot to reduce.
            KeyError: If ``name`` is not registered.
        """
        spec = self.spec(name)
        arrays = self.__authority[spec.owner]
        raw = spec.reduce(arrays)

        candidate = raw if isinstance(raw, np.ndarray) else None
        if candidate is not None:
            for owner, owned in self.__authority.items():
                for array_name, array in owned.items():
                    if np.shares_memory(candidate, array):
                        raise AuthoritativeReachError(
                            f"derivation {name!r} returned memory shared with authoritative array "
                            f"{owner}/{array_name}. A viewer holding that array reads physics as it "
                            "is being written and can be handed a write path by any later refactor. "
                            "Return a reduction, not the state."
                        )

        values = freeze_derived(raw)
        if spec.active_count is None:
            active = int(values.shape[0]) if values.ndim else 1
        else:
            active = int(spec.active_count(arrays))
        self._derive_count += 1
        return DerivedBuffer(
            name=spec.name,
            owner=spec.owner,
            units=spec.units,
            derivation=spec.derivation,
            values=values,
            active_count=active,
            accepted_step=int(accepted_step),
            topology_epoch=int(topology_epoch),
        )

    def derive_all(
        self, *, accepted_step: int, topology_epoch: int
    ) -> dict[str, DerivedBuffer]:
        """Run every registered reduction. Order is by name, so a manifest is reproducible."""
        return {
            name: self.derive(name, accepted_step=accepted_step, topology_epoch=topology_epoch)
            for name in self.derivation_names
        }

    # -- the guard's own hook ----------------------------------------------------------------
    def _authoritative_arrays_for_audit(self) -> tuple[tuple[str, np.ndarray], ...]:
        """Every authoritative array, labelled, for :func:`assert_no_authoritative_reach`.

        Named with a leading underscore and a suffix nobody would type by accident. The audit needs
        the arrays in order to test whether anything *else* aliases them, so there has to be one door
        — and the audit is the only caller, which is why it is not part of the public surface above.
        """
        return tuple(
            (f"{owner}/{name}", array)
            for owner, arrays in sorted(self.__authority.items())
            for name, array in sorted(arrays.items())
        )


_AUDIT_HOOK = "_authoritative_arrays_for_audit"

#: The private store, spelled the way ``dir()`` reports it. Exempt because it *is* the authority
#: rather than a copy of it, and reaching it requires typing the mangled name — which is what the
#: mangling is for. Everything else that holds the same memory is a second holder, and a second holder
#: is the leak.
_AUTHORITY_STORE = "_ObservationRegistry__authority"

#: Attribute names the audit is allowed to skip, with the reason.
_AUDIT_EXEMPT: Mapping[str, str] = {
    _AUDIT_HOOK: "the audit's own door onto the authority; not reachable by name from viewer code",
    _AUTHORITY_STORE: "the private authority store itself, which is the thing being protected",
    "spec": "returns the producer-side reduction callable, which runs before the boundary",
    "declare_authority": "the authority ingress; takes arrays, returns None",
    "register": "registration ingress",
}


def _reachable_arrays(obj: object, *, depth: int = 3) -> Iterable[tuple[str, np.ndarray]]:
    """Yield ``(path, array)`` for every ndarray reachable from ``obj`` within ``depth`` hops."""
    seen: set[int] = set()

    def walk(node: object, path: str, level: int) -> Iterable[tuple[str, np.ndarray]]:
        if level < 0 or id(node) in seen:
            return
        seen.add(id(node))
        if isinstance(node, np.ndarray):
            yield (path, node)
            return
        if isinstance(node, Mapping):
            for key, value in node.items():
                yield from walk(value, f"{path}[{key!r}]", level - 1)
            return
        if isinstance(node, (list, tuple, set, frozenset)):
            for index, value in enumerate(node):
                yield from walk(value, f"{path}[{index}]", level - 1)
            return
        slots = getattr(type(node), "__slots__", None)
        names: Sequence[str]
        if slots is not None and not isinstance(node, type):
            names = tuple(slots)
        elif hasattr(node, "__dict__") and not isinstance(node, type):
            names = tuple(vars(node))
        else:
            return
        for attr in names:
            if attr.startswith("__"):
                continue
            try:
                value = getattr(node, attr)
            except Exception:  # noqa: BLE001 - an accessor that raises exposes nothing
                continue
            yield from walk(value, f"{path}.{attr}", level - 1)

    return list(walk(obj, "", depth))


def assert_no_authoritative_reach(
    registry: ObservationRegistry,
    *,
    extra_objects: Sequence[object] = (),
) -> None:
    """Fail if anything on the viewer side of ``registry`` reaches authoritative memory.

    What it checks, and in each case what the failure would look like in production:

    * Every public attribute of the registry, and every zero-argument public method's return value,
      is walked for ndarrays. Any that shares memory with an authoritative array is a leak — the
      shape of a convenience accessor added later that returns the wrong object.
    * Any *writeable* ndarray reachable from that surface is a leak even when it aliases nothing,
      because a writeable buffer handed across the boundary is a write path waiting for an aliasing
      change upstream.
    * ``extra_objects`` gets the same treatment, which is how derived buffers and frames that have
      already left the registry are included.

    Args:
        registry: The registry to audit.
        extra_objects: Objects that came out of it, or out of a transport fed by it.

    Raises:
        AuthoritativeReachError: On the first leak, naming the path and the array.
        TypeError: If ``registry`` is not an :class:`ObservationRegistry`.
    """
    if not isinstance(registry, ObservationRegistry):
        raise TypeError("assert_no_authoritative_reach audits an ObservationRegistry")

    authority = dict(getattr(registry, _AUDIT_HOOK)())
    if not authority:
        raise ValueError(
            "the registry declares no authority, so this audit would pass by vacuity; declare the "
            "authoritative arrays before auditing, or the test proves nothing"
        )

    surface: list[tuple[str, object]] = []
    for attr in dir(registry):
        if attr.startswith("__") or attr in _AUDIT_EXEMPT:
            continue
        try:
            value = getattr(registry, attr)
        except Exception:  # noqa: BLE001
            continue
        if callable(value):
            try:
                surface.append((f"registry.{attr}()", value()))
            except TypeError:
                continue  # needs arguments; covered through derive() results in extra_objects
            except Exception:  # noqa: BLE001
                continue
        else:
            surface.append((f"registry.{attr}", value))
    for index, obj in enumerate(extra_objects):
        surface.append((f"extra[{index}]", obj))

    for label, obj in surface:
        for path, array in _reachable_arrays(obj):
            for name, authoritative in authority.items():
                if np.shares_memory(array, authoritative):
                    raise AuthoritativeReachError(
                        f"{label}{path} shares memory with authoritative array {name}. The viewer "
                        "side must hold derived buffers only; this is a window onto live physics."
                    )
            if array.flags.writeable:
                raise AuthoritativeReachError(
                    f"{label}{path} is a writeable array on the viewer side of the boundary. §11 "
                    "gives the viewer no write access; freeze it or do not export it."
                )
