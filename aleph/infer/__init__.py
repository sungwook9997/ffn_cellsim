"""Inference layer: likelihoods, Fisher information, identifiability, posteriors, typed refusal.

Everything here is written against a generic forward model, ``Callable[[NDArray], NDArray]`` mapping
parameters to a predicted observable.  The layer never imports an observation operator or an
analytic oracle: the oracles judge the runtime and must not be reachable from the code they judge,
and the observation operator is a separate lane whose interface is exactly the callable above.  A
consequence worth stating plainly --- every test in ``tests/infer`` supplies its own closed-form toy
model, so this layer is verifiable on its own, today, without any physics attached.

Modules:
    :mod:`~aleph.infer.refusal`: the typed refusal vocabulary, shared by everything below.
    :mod:`~aleph.infer.numdiff`: finite differences, the arbiter for every analytic gradient here.
    :mod:`~aleph.infer.likelihood`: Gaussian (known and a-priori per-bin variance) and
        Gamma/Wishart likelihoods for spectra, with analytic gradients.
    :mod:`~aleph.infer.fisher`: Fisher information, named stiff and sloppy directions, Cramer-Rao.
    :mod:`~aleph.infer.identifiability`: structural null space, degeneracy geometry, refusal.
    :mod:`~aleph.infer.posterior`: Laplace approximation, Metropolis check, coverage campaign.

Only the refusal vocabulary is re-exported at package level.  It is the part other lanes need in
order to *handle* an answer this layer declines to give, and it is the only module here that pulls
in no numerical dependencies, so importing it costs nothing.
"""

from .refusal import (
    InsufficientData,
    InsufficientDataError,
    ModelInadequate,
    ModelInadequateError,
    NativeQueryRequest,
    OutOfDistribution,
    OutOfDistributionError,
    ProvenanceIncomplete,
    ProvenanceIncompleteError,
    Refusal,
    RefusalError,
    Unidentifiable,
    UnidentifiableDirectionError,
    is_refusal,
)

__all__ = [
    "InsufficientData",
    "InsufficientDataError",
    "ModelInadequate",
    "ModelInadequateError",
    "NativeQueryRequest",
    "OutOfDistribution",
    "OutOfDistributionError",
    "ProvenanceIncomplete",
    "ProvenanceIncompleteError",
    "Refusal",
    "RefusalError",
    "Unidentifiable",
    "UnidentifiableDirectionError",
    "is_refusal",
]
