# amr1d
A simple 1D AMR code, it is part of a numerical Hydrodynamics course at UNAM.

It solves the 1D Euler equations on an adaptive block based cartesian mesh, 
based in the blockIDs idea used in the walicxe3d code.  
At this point the solution is obtained with a 1st order Godunov scheme with 
the HLLC fluxes.

*The code is work in progress as the course evolves.*

The strategy is not focussed on performance but rather on being easy to
understand.
