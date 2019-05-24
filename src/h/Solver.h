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

//! Calculate the first m0 Eigenvalues and their Eigenvectors
	virtual void Calculate_GV() {};

	virtual double* GetEigen() { return nullptr; }

	virtual double* GetEigenV() { return nullptr; }
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

//! Set value of nrhs
	void Set_NRHS(MKL_INT num) { nrhs = num; }
};

//! Generalized Eigenvalue solver: A in core solver using sparse storage and MKL FEAST
class CFEASTGVSolver
{
private:

	//! Stiffness Matrix
	CSparseMatrix* SparseK;

	//! Mass Matrix
	double* _M;

	MKL_INT* _iM;

	MKL_INT* _jM;

	MKL_INT NEQ;

public:

	CFEASTGVSolver(CSparseMatrix* K);

	//! Calculate the eigenvalue in interval [emin, emax] with a guess of total number of eigenvalue: m0
	//! m is actual number of eigenvalue in interval [emin, emax] and m should be smaller than m0
	void Calculate_GV_FEAST(double emin, double emax, MKL_INT& m0, MKL_INT& m, double* lambda, double* res, double* Q);
};

//! Generalized Eigenvalue solver: A in core solver using sparse storage and subspace method
class CSUBSPACESolver: public CSolver
{
public:
	CPARDISOSolver* Pardiso;

	MKL_INT NEQ;
	MKL_INT m0;

	double* _K;

	double* _M;

	double* Q;
	double* E;
	double* X;
	double* Y1;
	double* Y2;

	double* e;

public:
	CSUBSPACESolver(CSparseMatrix* M, MKL_INT Num_eig);

	//!	Perform L*D*L(T) factorization of the stiffness matrix
	virtual void LDLT() { Pardiso->LDLT(); };

	//!	Reduce right-hand-side load vector and back substitute
	virtual void BackSubstitution(double* Force) { Pardiso->BackSubstitution(Force); };

	//! Release Additional memory
	virtual void ReleasePhase() { Pardiso->ReleasePhase(); }

	//! Calculate the first m0 Eigenvalues and their Eigenvectors
	virtual void Calculate_GV();

	virtual double* GetEigen() { return E; }

	virtual double* GetEigenV() { return Q; }
};