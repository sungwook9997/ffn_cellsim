"""hoomd_legacy — the retired HOOMD-blue fine-grained runtime (archived 2026-06-29).

Demoted from production engine to frozen reference / parity oracle when the
going-forward engines moved to NVIDIA Warp (``ff/`` + ``dcm/``). These packages
— cell, cortex, bridge, ecm, junction, integrator, spheroid, gpu_opt — still run
under HOOMD-blue and are imported only by parity-fixture generators and legacy
scripts/tests, never by the active ``ff/`` or ``dcm/`` runtime.
"""
