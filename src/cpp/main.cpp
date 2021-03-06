/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.11, November 22, 2017                                       */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#include "Domain.h"
#include "Bar.h"
#include "Outputter.h"
#include "TecOutputter.h"
#include "Clock.h"
#include <iomanip>

using namespace std;

int main(int argc, char *argv[])
{
	if (argc != 2) //  Print help message
	{
	    cout << "Usage: stap++ InputFileName\n";
		exit(1);
	}

	string filename(argv[1]);
    size_t found = filename.find_last_of('.');

    // If the input file name is provided with an extension
    if (found != std::string::npos) {
        if (filename.substr(found) == ".dat")
            filename = filename.substr(0, found);
        else {
            // The input file name must has an extension of 'dat'
            cout << "*** Error *** Invalid file extension: "
                 << filename.substr(found+1) << endl;
            exit(1);
        }
    }

    string InFile = filename + ".dat";
	string OutFile = filename + ".out";
	string TecOutFile = filename + "result.dat";

	CDomain* FEMData = CDomain::Instance();
	CTecOutputter* TecOutput = CTecOutputter::Instance(TecOutFile);

    Clock timer;
    timer.Start();

//  Read data and define the problem domain
	if (!FEMData->ReadData(InFile, OutFile))
	{
		cerr << "*** Error *** Data input failed!" << endl;
		exit(1);
	}
    
    double time_input = timer.ElapsedTime();

//  Allocate global vectors and matrices, such as the Force, ColumnHeights,
//  DiagonalAddress and StiffnessMatrix, and calculate the column heights
//  and address of diagonal elements
	FEMData->AllocateMatrices();
    
//  Assemble the gloabl stiffness matrix
	FEMData->AssembleStiffnessMatrix();
    
    double time_assemble = timer.ElapsedTime();

//  Solve the linear equilibrium equations for displacements
	//! Calculate the Eigenvalue by MKL FEAST
	if (FEMData->GetMODEX() == 4)
	{
		CFEASTGVSolver* Solver = new CFEASTGVSolver(FEMData->GetSparseStiffnessMatrix());
		int NEQ = FEMData->GetNEQ();
		double emin = 0.0;
		double emax = 100000000.0;
		int m0 = 5, m;
		double* lambda = new double[m0];
		double* res = new double[m0];
		double* Q = new double[NEQ*m0];

		Solver->Calculate_GV_FEAST(emin, emax, m0, m, lambda, res, Q);
		ofstream eig;
		eig.open(filename + ".eig");
		eig << setiosflags(ios::scientific) << setprecision(10);
		eig << "The number of eigenvalue between the interval [" << emin << ", " << emax << "] is " << m << endl;
		for (int i = 0; i < m; i++)
		{
			eig << "Eigen Value " << i + 1 << " : " << lambda[i] << endl;
			eig << "The Eigen Vector is: " << endl;
			for (int j = 0; j < NEQ; j++)
			{
				eig << Q[i*NEQ + j] << "    ";
			}
			eig << endl;
		}
		TecOutput->OutputEIGModule(Q, lambda, m);
		return 0;
	}

	if (FEMData->GetMODEX() == 5)
	{
		int Num_eig=8;
		CSUBSPACESolver_CG* Solver = new CSUBSPACESolver_CG(FEMData->GetSparseStiffnessMatrix(),Num_eig);
		int NEQ = FEMData->GetNEQ();
		Solver->CGSUBSPACEIteration();
		double* EigenV = Solver->GetEigenV();
		double* Eigen = Solver->GetEigen();
		ofstream eig;
		eig.open(filename + ".eig");
		eig << setiosflags(ios::scientific) << setprecision(10);
		eig << "The number of eigenvalue is " << Num_eig << endl;
		for (int i = 0; i < Num_eig; i++)
		{
			eig << "Eigen Value " << i + 1 << " : " << Eigen[i] << endl;
			eig << "The Eigen Vector is: " << endl;
			for (int j = 0; j < NEQ; j++)
			{
				eig << EigenV[i*NEQ + j] << "    ";
			}
			eig << endl;
		}
		TecOutput->OutputEIGModule(EigenV, Eigen, Num_eig);
		return 0;
	}


	CSolver* Solver;
	if (FEMData->GetMODEX() == 2)
	{
		Solver = new CPARDISOSolver(FEMData->GetSparseStiffnessMatrix());
	}
	else if(FEMData->GetMODEX() == 1)
	{
		Solver = new CLDLTSolver(FEMData->GetStiffnessMatrix());
	}
	else if (FEMData->GetMODEX() == 3)
	{
		int Num_eig = 8;
		Solver = new CSUBSPACESolver(FEMData->GetSparseStiffnessMatrix(), Num_eig);
	}
    
//  Perform L*D*L(T) factorization of stiffness matrix
    Solver->LDLT();

    COutputter* Output = COutputter::Instance();

#ifdef _DEBUG_
    Output->PrintStiffnessMatrix();
#endif
//! Calculate Eigenvalues and Eigenvectors
	if (FEMData->GetMODEX() == 3)
	{
		int Num_eig = 8;
		int NEQ = FEMData->GetNEQ();
		Solver->Calculate_GV();
		double* EigenV = Solver->GetEigenV();
		double* Eigen = Solver->GetEigen();
		ofstream eig;
		eig.open(filename + ".eig");
		eig << setiosflags(ios::scientific) << setprecision(10);
		eig << "The number of eigenvalue is " << Num_eig << endl;
		for (int i = 0; i < Num_eig; i++)
		{
			eig << "Eigen Value " << i + 1 << " : " << Eigen[i] << endl;
			eig << "The Eigen Vector is: " << endl;
			for (int j = 0; j < NEQ; j++)
			{
				eig << EigenV[i*NEQ + j] << "    ";
			}
			eig << endl;
		}
		TecOutput->OutputEIGModule(EigenV, Eigen, Num_eig);
		return 0;
	}

//  Loop over for all load cases
    for (unsigned int lcase = 0; lcase < FEMData->GetNLCASE(); lcase++)
    {
//      Assemble righ-hand-side vector (force vector)
        FEMData->AssembleForce(lcase + 1);
            
//      Reduce right-hand-side force vector and back substitute
        Solver->BackSubstitution(FEMData->GetForce());
            
#ifdef _DEBUG_
        Output->PrintDisplacement(lcase);
#endif
            
        Output->OutputNodalDisplacement(lcase);

		Output->OutputNodalStress();

		//! Write result to Tecplot file in initial phase
		TecOutput->OutputResult(1, lcase);

		//! Write result to Tecplot file in deformed phase
		TecOutput->OutputResult(2, lcase);
    }

	Solver->ReleasePhase();

    double time_solution = timer.ElapsedTime();

//  Calculate and output stresses of all elements
//	Output->OutputElementStress();
    
    double time_stress = timer.ElapsedTime();
    
    timer.Stop();
    
    *Output << "\n S O L U T I O N   T I M E   L O G   I N   S E C \n\n"
            << "     TIME FOR INPUT PHASE = " << time_input << endl
            << "     TIME FOR CALCULATION OF STIFFNESS MATRIX = " << time_assemble - time_input << endl
            << "     TIME FOR FACTORIZATION AND LOAD CASE SOLUTIONS = " << time_solution - time_assemble << endl << endl
            << "     T O T A L   S O L U T I O N   T I M E = " << time_stress << endl;

	return 0;
}
