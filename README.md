# ROM_for_geomNL
A repo containing the FE codes and Manlab scripts for reduced order modelling of geometrically nonlinear structures with Normal Form method or Quadratic Manifold method.


The folders are divided in the conservative and the non conservative case.

Inside each folder, three methods are available:
- DNF SO: second order direct normal form
- QM SMD: quadratic manifold with static modal derivatives
- QM MD: quadratic manifold with modal derivatives

Inside each of them:
- the folder "astk" contains the code aster files for producing the ROM from FE
- the folder "manlab" contains the ManLab files for continuation with HBM


A short guide "HowToRunAstkJob" is included.
The algorithm of DNF is detailed in the pdf "DirectNormalFormAlgorithm"
A mesh of a clamped clamped beam (.med file) and some matlab figures with the results are included for verification.


Useful papers:
- Direct normal form implementation
https://doi.org/10.1016/j.cma.2021.113957
- Quadratic manifold implementation
https://doi.org/10.1007/s11071-020-05813-1

