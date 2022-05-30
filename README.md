# Beam-Prop
**Beam Propagation Method (BPM)** for planar **photonic integrated circuits (PIC)**, implemented in MATLAB with finite-differences.

## Brief description of the BPM-FD-2D

The BPM is a **paraxial spectral-domain stepping algorithm**, propagating a transversal spatial excitation at the input facet of structure (e.g. at one port of a PIC waveguide component) down to its output facet, along a well-defined (straight) optical axis; only 'slow' longitudinal perturbations of the structure are allowed, such as tapers, S-type bends or MMI-sections when dealing with PICs. The BPM is implemented with finite differences in 2D, where one axis is the waveguide cross-section (e.g. x), the other axis is the propagation direction (e.g. z). Perfectly matched layers (PML) truncate the cross-section in the lateral x-direction. Both scalar and semi-vector variants are included, corresponding to TE and TM polarized modes in the PIC waveguides. A variable wide-angle (Pade) correction can be used in conjunction with the Crank-Nicolson scheme, for the z-propagation algorithm.

Note that a **multi-layer slab waveguide (MLSWG) characteristic equation solver** is included, that can be used to calculate the modal excitation fed at the input port of the PIC. It can solve for arbitrary number of slab layers, with complex-valued indices, for both TE and TM polarizations. It can be used to find *all* modes inside a specified effective-index range, for a given wavelength.

## Examples included

Four scripts/examples are included, to ease the understanding of these routines and how to define the various input variables.
1. Solver used for Coupler Supermodes -- Extracts the (super)modes supported by a waveguide coupler, assumed as a "super-waveguide".
2. BPM along a simple coupler -- Excites one port of a waveguide coupler and propagates the field down to the end of the coupler.
3. BPM tap coupler -- Models a realistic 10 dB tap coupler, including S-bends to uncouple the I/O waveguide ports.
4. BPM MZI modulator -- Models an amplitude modulator, including a pair of MZI arms with input/output Y-splitter/combiners, respectively.

## Files

Here is a list of the m-files included in this repository
```
\BPM\BPMFD2D_DoProp.m
    \BPMFD2D_DrawLayout.m
    \BPMFD2D_PreProcLayout.m
\Solver\interpinv.m
       \MLSWG.m
       \MLSWG_CharEq.m
       \myNewtonRaphson.m
\Misc\flwcs.m
     \fmfp.m
     \LVCMv2.m
myScriptExample1_Solver_Supermodes.m
myScriptExample2_BPM_Simple_coupler.m
myScriptExample3_BPM_Tap_coupler.m
myScriptExample4_BPM_MZI_modulator.m
```

## References

If this code is used for research papers, be a pal and please cite some of my relevant works, e.g., [this one (thermo-optically tunable plasmonic components)](https://ieeexplore.ieee.org/abstract/document/5955059) or [that one (nonlinear graphene-comprising isolator)](https://ieeexplore.ieee.org/abstract/document/9395480). Both relate to finite-element method implementations, which I plan to soon<sup>TM<\sup> add, too.
		
