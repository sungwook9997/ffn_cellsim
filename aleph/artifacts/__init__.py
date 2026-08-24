"""Schemas, digests, and stamps — the layer that makes an Aleph number auditable.

* :mod:`aleph.artifacts.stamp` captures which build a run was measured on, refusing to resolve a
  disagreement between a declared commit and a git commit, and hashes configurations under a
  documented canonicalization.
* :mod:`aleph.artifacts.digest` content-addresses accepted state, with coverage declared rather
  than incidental.
* :mod:`aleph.artifacts.record` is the ``aleph-run-record@1`` schema, including the census guard
  that refuses a completeness claim without a per-owner breakdown.
"""

from aleph.artifacts.digest import (
    DIGEST_COVERAGE,
    OwnerStateBlock,
    RngState,
    SimulationClock,
    accepted_state_digest,
)
from aleph.artifacts.record import (
    ARTIFACT_SCHEMA,
    BANNED_CENSUS_KEYS,
    ArtifactRef,
    DeviceDescription,
    OverclaimingCensusError,
    RecordIdentity,
    RunRecord,
    TimingBlock,
    Verdict,
    reject_overclaiming_census,
)
from aleph.artifacts.stamp import (
    BuildSource,
    BuildStamp,
    CommitDisagreementError,
    RepositoryState,
    build_stamp,
    canonical_json,
    config_hash,
)

__all__ = [
    "ARTIFACT_SCHEMA",
    "BANNED_CENSUS_KEYS",
    "DIGEST_COVERAGE",
    "ArtifactRef",
    "BuildSource",
    "BuildStamp",
    "CommitDisagreementError",
    "DeviceDescription",
    "OverclaimingCensusError",
    "OwnerStateBlock",
    "RecordIdentity",
    "RepositoryState",
    "RngState",
    "RunRecord",
    "SimulationClock",
    "TimingBlock",
    "Verdict",
    "accepted_state_digest",
    "build_stamp",
    "canonical_json",
    "config_hash",
    "reject_overclaiming_census",
]
