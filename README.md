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
#### Date: 2019/5/16
6. Add nodal mass as EFEP90 code.
#### Date: 2019/5/21
7. Add a Solver for eigenvalue calculation with MKL [dfeast_scsrgv](https://software.intel.com/en-us/mkl-developer-reference-c-feast-scsrgv-feast-hcsrgv) function.
#### Date: 2019/5/22
8. Add visualization module for eigenvalue and eigenvector in Tecplot
#### Date: 2019/5/24
9. Add SUBSPACE Eigen Solver

### Modification by Cong Chen
#### Date:2019/4/28
1. Add element type: HexT （20 nodes）.The number of this _ElementType for input is 5;
#### Date:2019/5/5
2. Preprocess for truss .inp file from Abaqus.
3. Revise the bug in the HexT file

### Modification by Jia Sun
#### date:2019/05/07
1. Add element type: H8 (8 nodes). The number of this _ElementType for input is 4;
#### date:2019/05/22
2. H8.cpp has been modified.

### Modification by Weiyin Huang