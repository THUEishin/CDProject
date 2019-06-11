#pragma once
/*****************************************************************************/
/*  SparseMatrix: used to store the sparse matrix for MKL PARDISO Solver     */
/*     Added by Ruichen Ni, 2018311066                                       */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Date: April 27, 2019                                                  */
/*****************************************************************************/


class CSparseMatrix
{
private:
	//! Gloabl Stiffness Matrix
	double* _K;	

	//! Row index
	int* _iK;	

	//! Coloumn index
	int* _jK;	

	//! Dimension of the stiffness matrix
	unsigned int NEQ_;

	//! Total number of non-zero elements in the stiffness matrix
	unsigned int NNZ_;

public:
	//! constructors
	CSparseMatrix();

	//! destructor
	~CSparseMatrix();

	//! Return Dimension of the stiffness matrix
	inline unsigned int GetNEQ() { return NEQ_; }

	//! Return Total number of non-zero elements in the stiffness matrix
	inline unsigned int GetNNZ() { return NNZ_; }

	//! Assemble the element stiffness matrix to the global stiffness matrix
	void Assembly(double* Matrix, unsigned int* LocationMatrix, size_t ND);

	friend class CPARDISOSolver;
	friend class CFEASTGVSolver;
	friend class CSUBSPACESolver;
	friend class CSUBSPACESolver_CG;
};