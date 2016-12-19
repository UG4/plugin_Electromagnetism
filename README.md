# Electromagnetism #

**UG4-Plugin** implementing the FE discretization of the Maxwell equations using the Nedelec-type-1 (Whitney-1) elements.

The basic tools for the Nedelec elements are located in the root directory of the plugin.

Subdirectory EddyCurrent_E_Nedelec contains the implementation of the time-harmonic E-based formulation of the eddy current model, i.e. the discretization, the hybrid smoother by R. Hiptmair and auxiliary utilities.

For the mathematical details, we refere to

* A. Bossavit. Computational Electromagnetism. Academic Press (Boston), 1998
* R. Hiptmair. Multigrid method for Maxwell’s equations. SIAM J. Numer. Anal., 36(1): pp. 204–225 (1998), DOI 10.1137/S0036142997326203
* O. Sterz. Modellierung und Numerik zeitharmonischer Wirbelstromprobleme in 3D. PhD thesis, Heidelberg University, 2003
* O. Sterz, A. Hauser, G. Wittum. Adaptive local multigrid methods for solving time-harmonic eddy-current problems, IEEE Transactions on Magnetics, Volume 42, Issue 2, Feb. 2006, DOI: 10.1109/TMAG.2005.859064

Copyright 2009-2016 Goethe Center for Scientific Computing, University Frankfurt

Please install/clone this repository through UG4's package manager
[ughub](https://github.com/UG4/ughub).
