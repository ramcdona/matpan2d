# MATPAN2D -- Matlab Panel Code 2D

MATPAN2D is an attempt to create a fast, simple, and hackable design tool
for interesting aerodynamic phenemona.  These include power-effects,
inteference between bodies, etc.

MATPAN2D is released under the MIT License.

Rob McDonald

# What is MATPAN2D?

* Suite of tools for aerodynamic analysis and design of concepts with
  propulsion and/or interaction
* Mix of custom-developed and off-the-shelf code
* Simplifications aggressively accepted/pursued
  * Capture essential physics
  * Speed development & execution
  * Reduce design space
  * Clarify concept benefits or drawbacks
* Informed design preferred to optimization
* 2D Axisymmetric panel code (planar eventually)
  * Lumped ring-vortex with Dirichlet BC (V=0 on surface)
  * Based on: R. I. Lewis, _Vortex Element Methods for fluid Dynamic
    Analysis of Engineering Systems_ , 1991
* Actuator disk
  * Uniform loading, no swirl, contracted streamtube
  * Ducted -- zero or large tip-gap (> Rd/10) can be modeled
* Flexible & Fast
  * Custom Matlab implementation
  * Arbitrary combinations of bodies, ducts & disks
  * NACA 4-Digit duct generated on-the-fly
  * ~0.1 to 3 seconds for solution, ~10-15 seconds for visualization

# Planned Work

* Compressibility Correction -- Gothert's Rule
* BL2D Integration
  * NASA-TM-83207 Harris & Blanchard 1982
  * Loose coupling
    * Prescribed inviscid pressure distributions
  * Strong coupling (may be possible)
    * Momentum thickness feedback to inviscid solver
    * Converge on combined viscous/inviscid flow
* Design
  * Inviscid inverse design
  * Inviscid general pressure distribution matching
  * Design pressure distribution to tailor boundary layer

# Examples

[Unpowered Component Buildup](./CompBuildup.md)

[Powered Body](./PoweredBody.md)

[Powered Duct](./PoweredDuct.md)

[Powered Configuration](./PoweredBuildup.md)
