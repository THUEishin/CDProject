# CDProject
The project for the course 'Computational Dynamics'. This README.md is used to record the modification to the origin [FEM code](https://www.github.com/xzhang66/stappp) 

### Modification by Ruichen Ni
#### Date: 2019/4/25
1. Add element type: plain strain Quadrilateral(4 nodes). The number of this _ElementType for input is 2.
2. Replace element stress with nodal stress in output phase. Apply Average SPR Method to calculate nodal stress.
#### Date: 2019/4/26
3. Postproccess phase for [Tecplot](https://www.tecplot.com/) Display.
#### Date: 2019/4/28
4. Add [Sparse Storage Method](https://software.intel.com/en-us/mkl-developer-reference-c-sparse-matrix-storage-formats) and [MKL PARDISO](https://software.intel.com/en-us/mkl-developer-reference-c-pardiso) Solver to the code. To use the Solver, set MODEX to be 2.
5. There are some bugs in Tecplot output which have been revised.