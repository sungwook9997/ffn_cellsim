"""Evidence capture: the only path by which a viewer's output becomes citable.

Manuscript §11, evidence-capture row: *"Predeclared accepted cadence, build/config/census/frame/topology
IDs, units, exact cost, and required views."* And, from the dependency-boundary row: *"Evidence capture
may publish only from the unified accepted-step predicate."*

Each clause is a refusal here rather than a field to fill in.

* **Predeclared cadence.** :class:`CaptureCadence` records the accepted-step index the producer had
  reached when the cadence was fixed. :meth:`EvidenceCapture.mint` refuses a frame from a step at or
  before that index. Choosing a cadence after seeing which steps looked good is the visualisation
  equivalent of choosing a threshold from the runs that motivated it, and PLAN §11 already contains one
  measurement showing how much that can move an answer — twenty stopping points on one mesh spanning
  603x in the error, with five sign changes.
* **Publish only from an accepted step.** :meth:`mint` takes ``accepted`` explicitly and raises
  :class:`UnacceptedStepError` when it is false. There is no default, so a caller who does not know
  whether the step was accepted cannot get a frame by omission.
* **Required views.** The set of derivations that must be present in every frame is declared up front.
  A capture missing a declared view is refused, because "the frame we happened to have" and "the frame
  we said we would take" are different objects and only the second one is a protocol.
* **Exact cost.** Accumulated from an injectable clock, and reported as a
  :class:`~aleph.artifacts.record.TimingBlock` so it is the same cost object the run record uses. When
  the accumulated cost is exactly zero the manifest says the timing block is absent *and why*, rather
  than fabricating a positive number to satisfy the constructor.

The manifest's evidence label is produced by :func:`aleph.evidence.ladder.classify_claim` and it is
deliberately unflattering: a viewer's derived buffer is a projection of accepted state, not an
independent measurement of it, so the quantitative axis starts ``BLOCKED`` and stays there. What capture
buys is *retrievability* — a later reader can find the accepted state a picture came from — and that is
worth having on its own.
"""

from __future__ import annotations

import json
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable, Mapping, Sequence

from aleph.artifacts.record import TimingBlock, reject_overclaiming_census
from aleph.artifacts.stamp import BuildStamp, config_hash
from aleph.evidence.ladder import EvidenceLabel, EvidenceRung, QuantitativeStatus, classify_claim
from aleph.viz.buffers import DerivedBuffer
from aleph.viz.channels import ColourBar
from aleph.viz.frames import AcceptedFrame, AcceptedFrameStream, FrameIdentity
from aleph.viz.modes import ViewerMode
from aleph.viz.transport import transport_report

__all__ = [
    "CaptureCadence",
    "CaptureManifest",
    "EvidenceCapture",
    "MissingRequiredViewError",
    "OffCadenceError",
    "PrematureCadenceError",
    "UnacceptedStepError",
]


class UnacceptedStepError(RuntimeError):
    """A frame was requested for a step the acceptance predicate did not accept."""


class OffCadenceError(RuntimeError):
    """A frame was requested at a step the predeclared cadence does not sample."""


class PrematureCadenceError(RuntimeError):
    """A frame was requested from a step at or before the cadence's own declaration."""


class MissingRequiredViewError(RuntimeError):
    """A declared required view is absent from the offered buffers."""


@dataclass(frozen=True, slots=True)
class CaptureCadence:
    """A sampling rule fixed before the steps it samples.

    Attributes:
        every_n_accepted_steps: Sample every N accepted steps.
        declared_at_accepted_step: The producer's accepted-step count at the moment the cadence was
            written down. Frames from at or before this index are refused: they were already in the
            past when the rule was chosen.
        required_views: Derivation names every sampled frame must carry.
        note: Why this cadence, in words. Required — a cadence with no stated reason is a number
            somebody liked.
    """

    every_n_accepted_steps: int
    declared_at_accepted_step: int
    required_views: tuple[str, ...]
    note: str

    def __post_init__(self) -> None:
        if isinstance(self.every_n_accepted_steps, bool) or not isinstance(
            self.every_n_accepted_steps, int
        ):
            raise TypeError("every_n_accepted_steps must be an integer")
        if self.every_n_accepted_steps < 1:
            raise ValueError(
                f"cadence must be at least 1 accepted step; got {self.every_n_accepted_steps}"
            )
        if self.declared_at_accepted_step < 0:
            raise ValueError("declared_at_accepted_step cannot be negative")
        if not self.required_views:
            raise ValueError(
                "a capture with no required views has declared nothing it must show, so no frame "
                "can fail to satisfy it; §11 requires the views to be predeclared"
            )
        if len(set(self.required_views)) != len(self.required_views):
            raise ValueError(f"required_views repeats a name: {list(self.required_views)}")
        if not self.note.strip():
            raise ValueError(
                "a cadence needs a stated reason; without one there is no way to tell a designed "
                "sampling rate from the one that produced a convenient picture"
            )

    def samples(self, accepted_step: int) -> bool:
        """Whether this cadence samples ``accepted_step``."""
        return int(accepted_step) % self.every_n_accepted_steps == 0

    def check(self, accepted_step: int) -> None:
        """Raise unless a frame at ``accepted_step`` is inside this cadence's remit.

        Raises:
            PrematureCadenceError: If the step is at or before the declaration.
            OffCadenceError: If the cadence does not sample this step.
        """
        step = int(accepted_step)
        if step <= self.declared_at_accepted_step:
            raise PrematureCadenceError(
                f"accepted step {step} is at or before this cadence's declaration point "
                f"({self.declared_at_accepted_step}), so the rule was chosen after the step it "
                "would sample. A cadence fixed with the data already in hand is a selection, not a "
                "protocol."
            )
        if not self.samples(step):
            raise OffCadenceError(
                f"accepted step {step} is not sampled by the declared cadence of every "
                f"{self.every_n_accepted_steps} accepted steps. Publishing here would make the "
                "series' effective rate different from the one the manifest states."
            )

    def as_json_obj(self) -> dict[str, Any]:
        """Return the JSON-able cadence block."""
        return {
            "every_n_accepted_steps": self.every_n_accepted_steps,
            "declared_at_accepted_step": self.declared_at_accepted_step,
            "required_views": list(self.required_views),
            "note": self.note,
        }


@dataclass(frozen=True, slots=True)
class CaptureManifest:
    """Everything a later reader needs in order to trust — or distrust — a capture.

    Attributes:
        run_label: Free-text identity of the capture.
        mode: Always ``EVIDENCE_CAPTURE``.
        cadence: The predeclared sampling rule.
        identity: Build/config/census/schema identifiers.
        census: The population census, validated against the overclaiming guard.
        frames: One header per published frame.
        evidence: The two-axis evidence label.
        timing: Exact cost, when it was non-zero.
        timing_absent_reason: Why there is no timing block, when there is none.
        transport: What the frames went through, including its drop count.
        attachment: The viewer attachment's own report.
        refusals: Every frame request that was refused, with the reason. Non-empty is normal and
            healthy: it is the record of the guards firing.
    """

    run_label: str
    mode: ViewerMode
    cadence: CaptureCadence
    identity: FrameIdentity
    census: Mapping[str, Any]
    frames: tuple[Mapping[str, Any], ...]
    evidence: EvidenceLabel
    timing: TimingBlock | None
    timing_absent_reason: str
    transport: Mapping[str, Any]
    attachment: Mapping[str, Any]
    refusals: tuple[str, ...]

    def as_json_obj(self) -> dict[str, Any]:
        """Return the whole manifest as a JSON-able object."""
        block: dict[str, Any] = {
            "schema": "aleph-viz-capture@1",
            "run_label": self.run_label,
            "mode": self.mode.value,
            "cadence": self.cadence.as_json_obj(),
            "identity": self.identity.as_json_obj(),
            "census": dict(self.census),
            "n_frames": len(self.frames),
            "frames": [dict(frame) for frame in self.frames],
            "transport": dict(self.transport),
            "attachment": dict(self.attachment),
            "refusals": list(self.refusals),
        }
        if self.timing is not None:
            block["timing"] = self.timing.as_json_obj()
            # Mirrors the rule in `aleph.artifacts.record.RunRecord`: a cost is a property of a run at
            # a stated convergence. A capture evaluates no gate, so its cost is never comparable
            # against another capture's, and saying so here stops someone lifting the number.
            block["timing"]["comparable"] = False
            block["timing"]["not_comparable_reason"] = (
                "a capture evaluates no acceptance gate of its own, so this cost is the cost of "
                "observing whatever the producer happened to do; comparing it against another "
                "capture would compare two different trajectories"
            )
        else:
            block["timing"] = None
            block["timing_absent_reason"] = self.timing_absent_reason
        block.update(self.evidence.as_artifact_fields())
        return block

    def to_json(self, *, indent: int = 2) -> str:
        """Return the manifest as JSON text, newline-terminated."""
        return json.dumps(self.as_json_obj(), indent=indent, sort_keys=False) + "\n"

    def write_json(self, path: Path | str) -> Path:
        """Write the manifest to ``path``."""
        destination = Path(path)
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_text(self.to_json(), encoding="utf-8")
        return destination


class EvidenceCapture:
    """Mints accepted frames, and refuses every way of getting one that §11 forbids.

    Constructed *before* the steps it will sample, because the cadence it holds is predeclared by
    definition. The census is validated at construction by
    :func:`aleph.artifacts.record.reject_overclaiming_census`, so a capture whose population claims a
    completeness it cannot show does not exist rather than existing and being caught at write time.
    """

    def __init__(
        self,
        *,
        run_label: str,
        cadence: CaptureCadence,
        build: BuildStamp,
        config: Mapping[str, Any],
        census: Mapping[str, Any],
        schema_hash: str,
        transport: Any,
        accepted_state_digest_for: Callable[[int], str | None] | None = None,
        clock: Callable[[], float] = time.perf_counter,
    ) -> None:
        """Open a capture.

        Args:
            run_label: Identity of this capture.
            cadence: The predeclared sampling rule.
            build: The build stamp, from :func:`aleph.artifacts.stamp.build_stamp`.
            config: Run configuration, hashed into the frame identity.
            census: Population census.
            schema_hash: Digest of the state schema.
            transport: Where frames go. A lossless transport is expected; a lossy one is accepted and
                its drop count travels into every frame, which is the honest outcome rather than a
                refusal — a capture over a lossy channel is a capture with a stated hole in it.
            accepted_state_digest_for: ``f(accepted_step) -> digest | None``, supplied by the producer.
                Without it, frames carry no accepted-state digest and
                :attr:`AcceptedFrame.supports_quantitative_claim` is ``False`` for every frame — which
                is correct: a frame that cannot name the state it depicts supports no number.
            clock: Monotonic seconds source.

        Raises:
            OverclaimingCensusError: If the census overclaims completeness.
            TypeError: On a wrong argument type.
        """
        if not run_label.strip():
            raise ValueError("a capture needs a run_label")
        if not isinstance(cadence, CaptureCadence):
            raise TypeError("EvidenceCapture.cadence must be a CaptureCadence")
        if not isinstance(build, BuildStamp):
            raise TypeError("EvidenceCapture.build must be a BuildStamp from build_stamp()")
        reject_overclaiming_census(census)

        self.run_label = run_label
        self.cadence = cadence
        self.mode = ViewerMode.EVIDENCE_CAPTURE
        self._build = build
        self._config = dict(config)
        self._census = dict(census)
        self._schema_hash = schema_hash
        self._transport = transport
        self._digest_for = accepted_state_digest_for
        self._clock = clock

        self._stream = AcceptedFrameStream(
            cadence_every_n_accepted_steps=cadence.every_n_accepted_steps
        )
        self._refusals: list[str] = []
        self._cost_seconds = 0.0
        self._opened_at = clock()

    # -- introspection ------------------------------------------------------------------------
    @property
    def stream(self) -> AcceptedFrameStream:
        """The accepted stream this capture is filling."""
        return self._stream

    @property
    def frames(self) -> tuple[AcceptedFrame, ...]:
        """Published frames, in order."""
        return self._stream.frames

    @property
    def refusals(self) -> tuple[str, ...]:
        """Every refused request, with its reason."""
        return tuple(self._refusals)

    @property
    def cost_seconds(self) -> float:
        """Accumulated wall time inside :meth:`mint`."""
        return self._cost_seconds

    @property
    def census_digest(self) -> str:
        """Canonical hash of the census."""
        return config_hash(self._census)

    @property
    def config_digest(self) -> str:
        """Canonical hash of the configuration."""
        return config_hash(self._config)

    def identity_for(self, accepted_step: int) -> FrameIdentity:
        """Build the identity block for a frame at ``accepted_step``."""
        commit = self._build.git_commit or self._build.declared_commit
        if not commit:
            commit = f"unavailable:{self._build.source.value}"
        return FrameIdentity(
            build_commit=commit,
            config_digest=self.config_digest,
            census_digest=self.census_digest,
            schema_hash=self._schema_hash,
            accepted_state_digest=(
                self._digest_for(int(accepted_step)) if self._digest_for is not None else None
            ),
        )

    # -- the one way to get a frame ------------------------------------------------------------
    def mint(
        self,
        *,
        accepted: bool,
        accepted_step: int,
        attempted_step: int,
        physical_time_s: float,
        topology_epoch: int,
        active_count: int,
        colour_bar: ColourBar,
        buffers: Mapping[str, DerivedBuffer],
        dropped_frame_count: int,
        cost_seconds: float = 0.0,
        events: Sequence[str] = (),
        notes: Mapping[str, Any] | None = None,
    ) -> AcceptedFrame:
        """Mint and file one accepted frame, or refuse and say why.

        Args:
            accepted: Whether the acceptance predicate accepted this step. **No default.**
            accepted_step: Accepted-step index.
            attempted_step: Attempted-step index.
            physical_time_s: Simulated time.
            topology_epoch: Topology generation.
            active_count: Live population.
            colour_bar: Display mapping, with units and clamps.
            buffers: The derived buffers.
            dropped_frame_count: Publishable frames lost before this one.
            cost_seconds: Cost of deriving this frame.
            events: Typed events at this step.
            notes: Anything else.

        Returns:
            The filed frame.

        Raises:
            UnacceptedStepError: If ``accepted`` is false. This is the §11 dependency boundary: a
                viewer may publish only from the accepted-step predicate, and a rejected candidate's
                buffers depict state that no authoritative array ever held.
            PrematureCadenceError: If the step predates the cadence's declaration.
            OffCadenceError: If the cadence does not sample this step.
            MissingRequiredViewError: If a declared required view is absent.
        """
        started = self._clock()
        try:
            if not accepted:
                why = (
                    f"refusing to publish from attempted step {attempted_step}: the acceptance "
                    "predicate did not accept it, and §11 permits evidence capture to publish only "
                    "from the unified accepted-step predicate. A rejected candidate's buffers depict "
                    "state that was rolled back; they belong in the diagnostic layer."
                )
                self._refusals.append(why)
                raise UnacceptedStepError(why)
            try:
                self.cadence.check(accepted_step)
            except (PrematureCadenceError, OffCadenceError) as exc:
                self._refusals.append(str(exc))
                raise
            missing = [view for view in self.cadence.required_views if view not in buffers]
            if missing:
                why = (
                    f"accepted step {accepted_step} is missing required view(s) {missing}; the "
                    f"capture declared {list(self.cadence.required_views)} and a frame that omits "
                    "one is not the frame the protocol promised"
                )
                self._refusals.append(why)
                raise MissingRequiredViewError(why)

            frame = AcceptedFrame(
                mode=self.mode,
                identity=self.identity_for(accepted_step),
                accepted_step=int(accepted_step),
                attempted_step=int(attempted_step),
                physical_time_s=float(physical_time_s),
                topology_epoch=int(topology_epoch),
                active_count=int(active_count),
                colour_bar=colour_bar,
                buffers=dict(buffers),
                dropped_frame_count=int(dropped_frame_count),
                cost_seconds=float(cost_seconds),
                cadence_every_n_accepted_steps=self.cadence.every_n_accepted_steps,
                events=tuple(events),
                notes=dict(notes or {}),
            )
            self._stream.append(frame)
            return frame
        finally:
            self._cost_seconds += max(0.0, self._clock() - started)

    # -- the manifest -------------------------------------------------------------------------
    def blocking_reasons(self) -> tuple[str, ...]:
        """Why this capture's magnitudes may not be quoted.

        Always non-empty. The first reason is structural and does not go away: a render buffer is a
        reduction of accepted state chosen for legibility, and a reduction chosen for legibility is not
        a measurement protocol. The rest appear only when they apply.
        """
        reasons = [
            "a derived render buffer is a projection of accepted state selected for legibility, not "
            "an independent observation of it; quoting a number from it would quote the projection"
        ]
        dropped = int(getattr(self._transport, "dropped_total", 0) or 0)
        if dropped:
            reasons.append(
                f"the transport dropped {dropped} frame(s), so the series is not the cadence the "
                "manifest declares"
            )
        failures = tuple(getattr(self._transport, "capture_failures", ()) or ())
        if failures:
            reasons.append(
                f"{len(failures)} frame(s) could not be accepted by the transport: {failures[0]}"
            )
        clamped = [
            frame.scalar_channel for frame in self.frames if frame.colour_bar.is_clamped
        ]
        if clamped:
            reasons.append(
                f"{len(clamped)} frame(s) render a clamped channel, so their saturated bands cannot "
                "be inverted to a value"
            )
        untraced = [
            frame.accepted_step
            for frame in self.frames
            if not frame.identity.traces_to_accepted_state
        ]
        if untraced:
            reasons.append(
                f"{len(untraced)} frame(s) name no accepted-state digest, so the state they depict "
                "cannot be retrieved and re-measured"
            )
        if not self.frames:
            reasons.append("the capture published no frames, so there is nothing to quote")
        return tuple(reasons)

    def evidence_label(self, *, rung: EvidenceRung = EvidenceRung.CPU_UNIT) -> EvidenceLabel:
        """Grade the capture on both axes.

        ``requested_status`` is ``PROVISIONAL`` and never higher, and the blocking reasons will demote
        it to ``BLOCKED`` — which is the point. Asking for ``QUOTABLE`` here and letting the policy
        refuse would be asking the ladder to say no on the caller's behalf.
        """
        return classify_claim(
            rung=rung,
            basis=(
                f"{len(self.frames)} accepted frame(s) published from the accepted-step predicate at "
                f"a cadence of every {self.cadence.every_n_accepted_steps} accepted steps, on the "
                "CPU path with no device transport"
            ),
            requested_status=QuantitativeStatus.PROVISIONAL,
            blocking_reasons=self.blocking_reasons(),
        )

    def manifest(
        self,
        *,
        attachment_report: Mapping[str, Any] | None = None,
        rung: EvidenceRung = EvidenceRung.CPU_UNIT,
    ) -> CaptureManifest:
        """Assemble the manifest.

        Args:
            attachment_report: :meth:`aleph.viz.publisher.ViewerAttachment.report`, when there was one.
            rung: The structural rung to grade at. ``CPU_UNIT`` is correct for this lane and a caller
                cannot buy a stronger quantitative status with a stronger rung — that is rule 1 of
                :func:`~aleph.evidence.ladder.classify_claim`.

        Returns:
            The manifest.
        """
        total = float(self._clock() - self._opened_at)
        timing: TimingBlock | None = None
        absent = ""
        if total > 0.0:
            timing = TimingBlock(
                wall_seconds=total,
                n_attempted_steps=None,
                n_accepted_steps=len(self.frames) or None,
            )
        else:
            absent = (
                "the injected clock reported zero elapsed seconds, so there is no measured cost; a "
                "TimingBlock requires a positive wall time and fabricating one to satisfy the "
                "constructor would put an invented number where a measurement belongs"
            )
        return CaptureManifest(
            run_label=self.run_label,
            mode=self.mode,
            cadence=self.cadence,
            identity=self.identity_for(
                self.frames[-1].accepted_step if self.frames else 0
            ),
            census=dict(self._census),
            frames=tuple(frame.as_json_obj() for frame in self.frames),
            evidence=self.evidence_label(rung=rung),
            timing=timing,
            timing_absent_reason=absent,
            transport=transport_report(self._transport),
            attachment=dict(attachment_report or {}),
            refusals=self.refusals,
        )

    def write(
        self,
        directory: Path | str,
        *,
        attachment_report: Mapping[str, Any] | None = None,
        rung: EvidenceRung = EvidenceRung.CPU_UNIT,
    ) -> Path:
        """Write the manifest and every frame's values under ``directory``.

        Layout::

            <directory>/capture_manifest.json
            <directory>/frames/accepted_<step:08d>_<channel>.npy   (one per buffer per frame)

        Values are written as ``.npy`` rather than into the JSON, because a render buffer is bulk data
        and inlining it makes the manifest unreadable by the person it exists for.

        Returns:
            The manifest path.
        """
        root = Path(directory)
        (root / "frames").mkdir(parents=True, exist_ok=True)
        import numpy as np  # local: keeps the module importable where only manifests are read

        for frame in self.frames:
            for name in sorted(frame.buffers):
                buffer = frame.buffers[name]
                target = root / "frames" / f"accepted_{frame.accepted_step:08d}_{name}.npy"
                np.save(target, buffer.active_values)
        manifest = self.manifest(attachment_report=attachment_report, rung=rung)
        return manifest.write_json(root / "capture_manifest.json")
