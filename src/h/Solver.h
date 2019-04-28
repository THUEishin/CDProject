/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.11, November 22, 2017                                       */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#pragma once

#include "SkylineMatrix.h"
#include "SparseMatrix.h"
#include "mkl_types.h"

//!	Base class for a solver
/*	New solver should be derived from this base class, and match the storage scheme
	of the global stiffness matrix employed in Domain class. */
class CSolver
{
public:

	CSolver() {}

//!	Perform L*D*L(T) factorization of the stiffness matrix
	virtual void LDLT() = 0;
    
//!	Reduce right-hand-side load vector and back substitute
	virtual void BackSubstitution(double* Force) = 0;

//! Release Additional memory
	virtual void ReleasePhase() {}
};

//!	LDLT solver: A in core solver using skyline storage  and column reduction scheme
class CLDLTSolver : public CSolver
{
private:

	CSkylineMatrix<double>* K;

public:

//!	Constructor
	CLDLTSolver(CSkylineMatrix<double>* K) : CSolver() { this->K = K; };

//!	Perform L*D*L(T) factorization of the stiffness matrix
	virtual void LDLT();

//!	Reduce right-hand-side load vector and back substitute
	virtual void BackSubstitution(double* Force); 
};

//! PARDISO solver: A in core solver using sparse storage and MKL PARDISO
class CPARDISOSolver: public CSolver
{
private:

	CSparseMatrix* SparseM;

	MKL_INT NEQ;

//! Type of stiffness matrix
	MKL_INT mtype;

//! Number of right hand sides.
	MKL_INT nrhs;

	void *pt[64];

	MKL_INT iparm[64];
	MKL_INT maxfct, mnum, phase, error, msglvl;

	// Auxiliary variables.
	double ddum;			//Double dummy 
	MKL_INT idum;			//Integer dummy

public:

	CPARDISOSolver(CSparseMatrix* M);

//!	Perform L*D*L(T) factorization of the stiffness matrix
	virtual void LDLT();

//!	Reduce right-hand-side load vector and back substitute
	virtual void BackSubstitution(double* Force);

//! Release Additional memory
	virtual void ReleasePhase();
};