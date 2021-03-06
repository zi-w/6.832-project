6.832-project
=============

This repository contains the code for trajectory optimization for underactuated systems via TrajOpt [1], an approach based on sequential convex optimization (SCO). SCO solves non-convex problems by solving a local approximated convex sub-problem iteratively. We implemented the algorithm using a matlab-based convex optimizer CVX [2,3] for the sub-problem within the framework of Drake [4]. Different dynamics constraints were considered including forward Euler (trajopt_forward.m), backward Euler (trajopt_backward.m), mid-point Euler (trajopt_midpoint.m), and collocation constraints (trajopt_collocation.m).

The simulation results are partially available here http://people.csail.mit.edu/ziw/6.832/.

To run the matlab files, Drake [4] and CVX [3] are required to be installed. An example showing how to run trajectory optimization with our code is shown in run_nlink.m. run_nlink.m should be put into the directory in drake, drake/examples/PlanarNLink/, together with the functions we provided.
 

References

[1] J. Schulman, Y. Duan, J. Ho, A. Lee, I. Awwal, H. Bradlow, J. Pan,
S. Patil, K. Goldberg, and P. Abbeel. Motion planning with sequential
convex optimization and convex collision checking. The International
Journal of Robotics Research, 33(9):1251–1270, 2014.

[2] M. Grant and S. Boyd. Graph implementations for nonsmooth convex
programs. In V. Blondel, S. Boyd, and H. Kimura, editors, Recent
Advances in Learning and Control, Lecture Notes in Control and
Information Sciences, pages 95–110. Springer-Verlag Limited, 2008.
http://stanford.edu/ ̃boyd/graph_dcp.html.

[3] M. Grant and S. Boyd. CVX: Matlab software for disciplined convex
programming, version 2.1. http://cvxr.com/cvx, Mar. 2014.

[4] R. Tedrake. Drake: A planning, control, and analysis toolbox for
nonlinear dynamical systems, 2014.
