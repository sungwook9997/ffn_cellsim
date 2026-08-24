# AFM protocol selection audit - Chugh et al. 2017

## Decision

Use Chugh et al. 2017 (DOI `10.1038/ncb3525`) as the HeLa cortical-tension magnitude and rounded-cell
state oracle. Do **not** describe the experiment-only spherical indenter as a reproduction of Chugh. The
paper's tension apparatus is a tipless, approximately flat cantilever compressing a rounded cell against a
glass-bottom dish; it is not a spherical probe.

The spherical indenter therefore remains a separate in-silico perturbation experiment. Its force-depth curve
can test the engine and rate dependence after C-2 passes, but Chugh cannot supply its radius, contact distance,
penalty stiffness, or spherical-contact inversion.

## What the source directly fixes

The source Methods and Supplementary Fig. 1 specify:

- rounded HeLa states: trypsinized adherent HeLa and suspension HeLa, in interphase and mitosis;
- a Nanowizard1 AFM with a tipless silicon-nitride cantilever (`HQ:CSC38/tipless/NoAl`);
- cantilever spring constant calibrated per experiment, reported as `0.07-0.12 N/m`;
- lowering speed `0.5 um/s`;
- force set point `10 nN`, producing `1-4 um` average whole-cell compression;
- constant-height hold for `300 s`, with force read after the initial relaxation;
- cortical-tension inversion from cantilever force, maximum cross-section radius, contact radius and cell
  height (Methods equations 1-2), assuming negligible adhesion to dish and cantilever.

These values are direct evidence for a future **flat-cantilever Chugh reproduction**. They are not silently
transferred to the spherical-indenter schedule.

## What the source does not fix

Chugh does not report the engine's resting osmotic setpoint `Pi_0`, exterior-medium viscosity, or a numerical
contact penalty. A rigid experimental cantilever calibration is not a calibration of the connector penalty.
Consequently:

- `Pi_0 = 40 Pa` remains an explicitly cross-protocol HeLa proxy and licenses only a mechanism demonstration;
- spherical contact distance and stiffness remain PI gates;
- exterior viscosity remains a PI gate;
- the `0.5 um/s` loading speed is direct for the flat-cantilever protocol only;
- no quantitative spherical-indenter cortical-tension result is licensed by this paper.

## Implementation split

1. First close C-2 on the full-native rounded cell without accepting a rejected candidate.
2. Keep the spherical sixth CONTACT edge as the engine perturbation harness and require its own sourced or
   PI-approved calibration.
3. Implement a separate flat-cantilever apparatus before claiming a Chugh reproduction; it must include the
   lower dish contact and the paper's geometry-based tension inversion rather than relabeling the spherical
   force-depth curve.
4. For either apparatus, use at least three distinct seeds and report seed scatter before resolving a tension
   difference.

## Source files inspected

- `/Users/sw1/ffn_cellsim/aleph/references/emss-72183.pdf`, Methods page 12.
- `/Users/sw1/ffn_cellsim/aleph/references/ncb3525 (1).pdf`, Supplementary Fig. 1 (PDF page 12).

