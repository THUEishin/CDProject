/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.11, November 22, 2017                                       */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#include "Solver.h"

#include <cmath>
#include <cfloat>
#include <iostream>
#include <algorithm>
#include <mkl.h>
#include "../h/Domain.h"
#include "mkl_pardiso.h"
#include "stdlib.h"
#include <fstream>
#include <iomanip>

using namespace std;

// LDLT facterization
void CLDLTSolver::LDLT()
{
	unsigned int N = K->dim();
    unsigned int* ColumnHeights = K->GetColumnHeights();   // Column Hights

	for (unsigned int j = 2; j <= N; j++)      // Loop for column 2:n (Numbering starting from 1)
	{
        // Row number of the first non-zero element in column j (Numbering starting from 1)
		unsigned int mj = j - ColumnHeights[j-1];
        
		for (unsigned int i = mj+1; i <= j-1; i++)	// Loop for mj+1:j-1 (Numbering starting from 1)
		{
            // Row number of the first nonzero element in column i (Numbering starting from 1)
			unsigned int mi = i - ColumnHeights[i-1];

			double C = 0.0;
			for (unsigned int r = max(mi, mj); r <= i-1; r++)
				C += (*K)(r,i) * (*K)(r,j);		// C += L_ri * U_rj

			(*K)(i,j) -= C;	// U_ij = K_ij - C
		}

		for (unsigned int r = mj; r <= j-1; r++)	// Loop for mj:j-1 (column j)
		{
			double Lrj = (*K)(r,j) / (*K)(r,r);	// L_rj = U_rj / D_rr
			(*K)(j,j) -= Lrj * (*K)(r,j);	// D_jj = K_jj - sum(L_rj*U_rj, r=mj:j-1)
			(*K)(r,j) = Lrj;
		}

        if (fabs((*K)(j,j)) <= FLT_MIN)
        {
            cerr << "*** Error *** Stiffness matrix is not positive definite !" << endl
            	 << "    Euqation no = " << j << endl
            	 << "    Pivot = " << (*K)(j,j) << endl;
            
            exit(4);
        }
    }
};

// Solve displacement by back substitution
void CLDLTSolver::BackSubstitution(double* Force)
{
	unsigned int N = K->dim();
    unsigned int* ColumnHeights = K->GetColumnHeights();   // Column Hights

//	Reduce right-hand-side load vector (LV = R)
	for (unsigned int i = 2; i <= N; i++)	// Loop for i=2:N (Numering starting from 1)
	{
        unsigned int mi = i - ColumnHeights[i-1];

		for (unsigned int j = mi; j <= i-1; j++)	// Loop for j=mi:i-1
			Force[i-1] -= (*K)(j,i) * Force[j-1];	// V_i = R_i - sum_j (L_ji V_j)
	}

//	Back substitute (Vbar = D^(-1) V, L^T a = Vbar)
	for (unsigned int i = 1; i <= N; i++)	// Loop for i=1:N
		Force[i-1] /= (*K)(i,i);	// Vbar = D^(-1) V

	for (unsigned int j = N; j >= 2; j--)	// Loop for j=N:2
	{
        unsigned int mj = j - ColumnHeights[j-1];

		for (unsigned int i = mj; i <= j-1; i++)	// Loop for i=mj:j-1
			Force[i-1] -= (*K)(i,j) * Force[j-1];	// a_i = Vbar_i - sum_j(L_ij Vbar_j)
	}
};


//! Constructer
CPARDISOSolver::CPARDISOSolver(CSparseMatrix* M) : CSolver()
{
	this->SparseM = M;
	NEQ = SparseM->NEQ_;
	mtype = -2;  //Real symmetric matrix
	nrhs = 1;
}

// LDLT facterization
void CPARDISOSolver::LDLT()
{
	MKL_INT i;

//-------------------------------------------------------------------- 
// Setup Pardiso control parameters. 
//--------------------------------------------------------------------
	for (i = 0; i < 64; i++)
	{
		iparm[i] = 0;
	}
	iparm[0] = 1;			/* No solver default */
	iparm[1] = 2;			/* Fill-in reordering from METIS */
							/* Numbers of processors, value of OMP_NUM_THREADS */
	iparm[2] = 1;
	iparm[3] = 0;			/* No iterative-direct algorithm */
	iparm[4] = 0;			/* No user fill-in reducing permutation */
	iparm[5] = 0;			/* Write solution into x */
	iparm[6] = 0;			/* Not in use */
	iparm[7] = 2;			/* Max numbers of iterative refinement steps */
	iparm[8] = 0;			/* Not in use */
	iparm[9] = 13;		/* Perturb the pivot elements with 1E-13 */
	iparm[10] = 1;		/* Use nonsymmetric permutation and scaling MPS */
	iparm[11] = 0;		/* Not in use */
	iparm[12] = 0;		/* Maximum weighted matching algorithm is switched-off (default for symmetric). Try iparm[12] = 1 in case of inappropriate accuracy */
	iparm[13] = 0;		/* Output: Number of perturbed pivots */
	iparm[14] = 0;		/* Not in use */
	iparm[15] = 0;		/* Not in use */
	iparm[16] = 0;		/* Not in use */
	iparm[17] = -1;		/* Output: Number of nonzeros in the factor LU */
	iparm[18] = -1;		/* Output: Mflops for LU factorization */
	iparm[19] = 0;		/* Output: Numbers of CG Iterations */
	maxfct = 1;			/* Maximum number of numerical factorizations. */
	mnum = 1;			/* Which factorization to use. */
	msglvl = 1;			/* Print statistical information in file */
	error = 0;			/* Initialize error flag */

//-------------------------------------------------------------------- 
// Initialize the internal solver memory pointer. This is only 
// necessary for the FIRST call of the PARDISO solver.
//--------------------------------------------------------------------
	for (i = 0; i < 64; i++)
	{
		pt[i] = 0;
	}

//--------------------------------------------------------------------
// Reordering and Symbolic Factorization. This step also allocates
// all memory that is necessary for the factorization
//--------------------------------------------------------------------
	phase = 11;
	PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &NEQ, SparseM->_K, SparseM->_iK, SparseM->_jK, &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
	if (error != 0)
	{
		printf("\nERROR during symbolic factorization: %d", error);
		system("Pause");
		exit(1);
	}
	printf("\nReordering completed ... ");
	printf("\nNumber of nonzeros in factors = %d", iparm[17]);
	printf("\nNumber of factorization MFLOPS = %d", iparm[18]);

//--------------------------------------------------------------------
// Numerical factorization
//--------------------------------------------------------------------
	phase = 22;
	PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &NEQ, SparseM->_K, SparseM->_iK, SparseM->_jK, &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
	if (error != 0)
	{
		printf("\nERROR during numerical factorization: %d", error);
		system("Pause");
		exit(2);
	}
	printf("\nFactorization completed ... ");
}

// Solve displacement by back substitution
void CPARDISOSolver::BackSubstitution(double* Force)
{
	double* dis = new double[NEQ*nrhs];
//--------------------------------------------------------------------
// Back substitution and iterative refinement
//--------------------------------------------------------------------
	phase = 33;
	iparm[7] = 2;			// Max numbers of iterative refinement steps.
	PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &NEQ, SparseM->_K, SparseM->_iK, SparseM->_jK, &idum, &nrhs, iparm, &msglvl, Force, dis, &error);
	if (error != 0)
	{
		printf("\nERROR during solution: %d", error);
		system("Pause");
		exit(3);
	}
	printf("\nSolve completed ... ");
	printf("\nThe solution of the system is: ");

	for (int i = 0; i < NEQ*nrhs; i++)
	{
		Force[i] = dis[i];
	}
	delete[] dis;
}

// Release Additional memory
void CPARDISOSolver::ReleasePhase()
{
//--------------------------------------------------------------------
// Termination and release of memory.
// -------------------------------------------------------------------
	phase = -1;
	PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &NEQ, SparseM->_K, SparseM->_iK, SparseM->_jK, &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
	return;
}

CFEASTGVSolver::CFEASTGVSolver(CSparseMatrix* K)
{
	SparseK = K;
	NEQ = SparseK->NEQ_;

	CDomain* FEMData = CDomain::Instance();

	int NUMNP = FEMData->GetNUMNP();
	CNode* nodelist = FEMData->GetNodeList();

	_M = new double[NEQ];
	_jM = new MKL_INT[NEQ];
	_iM = new MKL_INT[NEQ + 1];

	for (int i = 0; i < NEQ; i++)
	{
		_iM[i] = i + 1;
		_jM[i] = i + 1;
	}

	_iM[NEQ] = NEQ + 1;

	int NDF = CNode::NDF;
	for (int i = 0; i < NUMNP; i++)
	{
		CNode* node = &nodelist[i];
		for (int j = 0; j < NDF; j++)
		{
			if (node->bcode[j] > 0)
			{
				_M[node->bcode[j] - 1] = node->mass;
			}
		}
	}
}

void CFEASTGVSolver::Calculate_GV_FEAST(double emin, double emax, MKL_INT& m0, MKL_INT& m, double* lambda, double* res, double* Q)
{
	char uplo = 'U';
	MKL_INT fpm[128];
	double epsout;
	MKL_INT loop;
	MKL_INT info;

	feastinit(fpm);

	dfeast_scsrgv(&uplo, &NEQ, SparseK->_K, SparseK->_iK, SparseK->_jK, _M, _iM, _jM, fpm, &epsout, &loop, &emin, &emax, &m0, lambda, Q, &m, res, &info);

	if (info != 0)
	{
		cout << "There is error in FEAST Solver. Error Code is " << info << endl;
		exit(1);
	}
}

CSUBSPACESolver::CSUBSPACESolver(CSparseMatrix* M, MKL_INT Num_eig)
{
	Pardiso = new CPARDISOSolver(M);
	m0 = Num_eig;
	NEQ = M->NEQ_;

	_K = new double[m0*m0];

	_M = new double[m0*m0];

	Q = new double[NEQ*m0];
	E = new double[m0];
	X = new double[NEQ*m0];
	Y1 = new double[NEQ*m0];
	Y2 = new double[NEQ*m0];
	e = new double[m0];

	for (int i = 0; i < m0; i++)
	{
		for (int j = 0; j < NEQ; j++)
		{
			if (i == 0)
			{
				Q[j+i*NEQ] = 1.0;
			}
			else
			{
				if ((j + 1) == i)
					Q[j+i*NEQ] = 1.0;
				else
					Q[j+i*NEQ] = 0.0;
			}
		}
	}

	for (int i = 0; i < m0; i++)
	{
		E[i] = 0.0;
	}
}

void CSUBSPACESolver::Calculate_GV()
{
	CDomain* FEMData = CDomain::Instance();

	int NUMNP = FEMData->GetNUMNP();
	CNode* nodelist = FEMData->GetNodeList();

	int NDF = CNode::NDF;
	int bcode;
	for (int i = 0; i < NUMNP; i++)
	{
		CNode* node = &nodelist[i];
		for (int j = 0; j < NDF; j++)
		{
			bcode = node->bcode[j];
			if (bcode > 0)
			{
				for (int k = 0; k < m0; k++)
				{
					Y1[k*NEQ + bcode - 1] = node->mass*Q[k*NEQ + bcode - 1];
				}
			}
		}
	}

	double tol = 1.0e-5;
	double error = 1.0;
	MKL_INT info;
	Pardiso->Set_NRHS(m0);
	while (error > tol)
	{
		for (int i = 0; i < NEQ*m0; i++)
			X[i] = Y1[i];

		//! Calculate X by solve K*X = Y1
		BackSubstitution(X);

		//! Calculate Y2 by multiplying Y2 = M*X
		for (int i = 0; i < NUMNP; i++)
		{
			CNode* node = &nodelist[i];
			for (int j = 0; j < NDF; j++)
			{
				bcode = node->bcode[j];
				if (bcode > 0)
				{
					for (int k = 0; k < m0; k++)
					{
						Y2[k*NEQ + bcode - 1] = node->mass*X[k*NEQ + bcode - 1];
					}
				}
			}
		}

		//! Calculate _K by multiplying _K = X^T*Y1
		//! Calculate _M by multiplying _M = X^T*Y2
		clear(_K, m0*m0);
		clear(_M, m0*m0);
		for (int i = 0; i < m0; i++)
		{
			for (int j = 0; j < m0; j++)
			{
				for (int k = 0; k < NEQ; k++)
				{
					_K[i*m0+j] += X[i*NEQ + k] * Y1[j*NEQ + k];
					_M[i*m0+j] += X[i*NEQ + k] * Y2[j*NEQ + k];
				}
			}
		}

		info = LAPACKE_dsygv(LAPACK_ROW_MAJOR, 1, 'V', 'U', m0, _K, m0, _M, m0, e);

		error = 0;
		
		for (int i = 0; i < m0; i++)
		{
			error = max(error, abs(e[i] - E[i]) / e[i]);
			E[i] = e[i];
		}

		clear(Y1, NEQ*m0);
		for (int i = 0; i < m0; i++)
		{
			for (int j = 0; j < NEQ; j++)
			{
				for (int k = 0; k < m0; k++)
				{
					Y1[i*NEQ + j] += Y2[k*NEQ + j] * _K[k*m0 + i];
				}
			}
		}
	}

	for (int i = 0; i < NEQ*m0; i++)
		X[i] = Y1[i];

	//! Calculate X by solve K*X = Y1
	BackSubstitution(X);

	for (int i = 0; i < m0; i++)
	{
		for (int j = 0; j < NEQ; j++)
		{
			for (int k = 0; k < m0; k++)
			{
				Q[i*NEQ + j] += X[k*NEQ + j] * _K[k*m0 + i];
			}
		}
	}

	Pardiso->Set_NRHS(1);
	return;
}