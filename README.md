# LQ_Kernel_StateConstraints_SICON
Matlab code to play with the kernel solution to the state-constrained Linear Quadratic Regulator (LQR) and reproduce results from P.-C. Aubin-Frankowski "Linearly-constrained Linear Quadratic Regulator from the viewpoint of kernel methods", SICON 2021 [article]https://epubs.siam.org/doi/abs/10.1137/20M1348765 [arxiv]https://arxiv.org/abs/2011.02196# One can find in https://pcaubin.github.io/ several presentations related to the topic (or freely contact the author).

The key message is that state-constrained continuous-time LQR can be solved easily with matrix-valued kernels. Owing to Van Loan's technique, the latter are computed at lightspeed in the time-invariant case and for a quadratic Lagrangian only depending on the control. Note that using the eta parameter and SOCP is not compulsory and that the discretized constraints already give a very satisfactory sparse solution for fine grids.

The source files in this repository can of course be used for implementing problems. They could otherwise be used as a fruitful inspiration since several files need to be manually adapted to fit the user's specific problem.

## Description of the code /AUBIN_SICON2021_code used for the linear pendulum experiment

StateConstrained_pendulum_SICON_generation.m : __where to manually define the inputs of the problem__, gives a formatted output with unique time points. This ouput is an input of StateConstrained_pendulum_SICON_solver.m (uses ConvertingConsProblem_toUniquePoints.m as a subroutine)

StateConstrained_pendulum_SICON_solver.m : __where to manually add extra inputs of the optimization problem__, solves the optimization problem with CVX+MOSEK. This part could be replaced by the user with his favourite QPQC or SOCP solver (uses GramianComputingVanLoan.m to compute the Gram matrices as suggested by https://ieeexplore.ieee.org/document/1101743; and EtaComputingVanLoan_wVaryingC.m as a subroutine to approximate the values of the eta parameter of the SOCP technique)

StateConstrained_pendulum_SICON_minimal.m: __where to manually define some global variables of the problem__, provides a figure and an example of how to run and visualize the matrix-valued kernel solution to the Linear Quadratic Regulator.

StateConstrained_pendulum_SICON_figures.m: __where to manually define some global variables of the problem__, provides all the figures of the article which can be found in the directory /AUBIN_SICON2021_figures for comparison.

GramianComputingVanLoan.m: *the most important file*, implementation that you can use to quickly compute kernel Gram matrices (for instance Gramians of controllability) in the case Q=0
