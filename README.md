# euler-2D-kfvs
Euler equations unstructured solver using kinetic flux vector splitting (KFVS) scheme.

Some features:
- Finite volume method.
- Unstructured quadrilateral meshes.
- Kinetic flux vector splitting (KFVS) scheme. 
- Written in Modern Fortran (Fortran 95 + some Fortran 2003). No real OOP has been implemented yet.
- Use core code of gmsh-to-vtk-and-tecplot as a mesh reader.

euler-2D-kfvs is inspired by the structured finite volume solver of [Prof. Luc Mieussens](https://www.math.u-bordeaux.fr/~lmieusse/PAGE_WEB/ENSEIGNEMENT/MMK3/SIMULATION_NUMERIQUE_ECOULEMENTS_FLUIDES/simulations.html)

 
[More detailed description](https://github.com/truongd8593/euler-2D-kfvs/wiki)

Remark:
 - This Fortran 95 code is no longer maintained. I have upgraded it to Fortran 2003, following an object-oriented model. The code is [here](https://github.com/truongd8593/euler2D-kfvs-Fortran2003).
