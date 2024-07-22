# Finite-size analysis of prepare-and-measure and decoy-state QKD via entropy accumulation

This is a public version of the code used in *Finite-size analysis of prepare-and-measure and decoy-state QKD via entropy accumulation* [arxiv link](https://arxiv.org/abs/2406.10198). This code was designed during the early stages of the [Open QKD Security package](https://github.com/Optical-Quantum-Communication-Theory/openQKDsecurity) version 2 rewrite. As such, we provide it here as a standalone vesion, separate from the main repository. A small number of required functions from the main repository are included.

\<Write description of code here. For example, what figures and data it generates, and the methods it uses.\>


## Install instructions
> [!CAUTION]
> This repository is for archival and transparency purposes; we do not guarantee compatibility with updated versions of the Open QKD Security package.


These installation instructions were adapted from the main repository.
Before installing the software, ensure you have the latest version of MATLAB installed on your machine.
Our software requires *at least version 2020b* for full functionality, but installing the latest version is preferable.

Our software has the following dependencies for its default settings:

- [CVX](https://cvxr.com/cvx/download/) v2.2, a library for solving convex optimization problems in MATLAB.
- [QETLAB](https://github.com/nathanieljohnston/QETLAB) *above* v0.9, a MATLAB toolbox for operations involving quantum channels and entanglement. Note that you cannot use the version from their website as it has a bugs associated with implementing Choi matrices. *You must use their latest copy on GitHub*. At the time of writing, this means downloading their code with the green "code" button and *not* the v0.9 release.
- [MOSEK](https://www.mosek.com/) (optional), a more advanced semidefinite programming (SDP) solver than the default (SDPT3) used in CVX. Note that the MOSEK solver can be downloaded together with CVX, but requires a separate license to use. See [this page](https://cvxr.com/cvx/doc/mosek.html) for more information.

Please refer to the documentation of each of these software packages for installation instructions.

To install the openQKDSecurity package, download the repository on GitHub as a zip file and unzip to a preferred directory. This can be done on the main page of the repository by pressing the green “Code” button at the top of the page and then selecting "Download ZIP".

Next, ensure that our software and the above dependencies are on your MATLAB path.To place a folder on your path, navigate to it in MATLAB, then right-click and select "Add to Path\textgreater Selected folder and Subfolders". Make sure you do this for OpenQKDSecurity, QETLAB and CVX. We also recommend you run "cvx_setup" to check if CVX is properly configured and which solvers are available.

Before you close MATLAB, go to "Home\>Environment\>Set Path" and click "save" at the bottom. If you don't, then you have to add everything to your path each time you restart MATLAB.

\<Install directions for this repository. For example, add this folder to the Matlab path and save. Run this test function. Etc.\>
