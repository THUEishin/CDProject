/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.11, November 22, 2017                                       */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/


#include "../h/T4.h"

#include <iostream>
#include <iomanip>
#include <cmath>


//	Constructor
CT4::CT4()
{
	NEN_ = 4;	// Each element has 20 nodes
	nodes_ = new CNode * [NEN_];


	ND_ = 12;	// Each node has 3 dimension
	LocationMatrix_ = new unsigned int[ND_];


	ElementMaterial_ = nullptr;
}


//	Deconstructor
CT4::~CT4()
{
	if (nodes_)
		delete[] nodes_;

	if (LocationMatrix_)
		delete[] LocationMatrix_;
}


//	Read element data from stream Input
bool CT4::Read(ifstream& Input, unsigned int Ele, CMaterial* MaterialSets, CNode* NodeList)
{
	unsigned int N;


	Input >> N;	// element number

	if (N != Ele + 1)
	{
		cerr << "*** Error *** Elements must be inputted in order !" << endl
			<< "    Expected element : " << Ele + 1 << endl
			<< "    Provided element : " << N << endl;


		return false;
	}


	unsigned int MSet;	// Material property set number
	unsigned int N1, N2, N3, N4;	// Four nodes number

	Input >> N1 >> N2 >> N3 >> N4  >> MSet;
	ElementMaterial_ = dynamic_cast<CT4Material*>(MaterialSets) + MSet - 1;
	nodes_[0] = &NodeList[N1 - 1];
	nodes_[1] = &NodeList[N2 - 1];
	nodes_[2] = &NodeList[N3 - 1];
	nodes_[3] = &NodeList[N4 - 1];

	for (int i = 0; i < 4; i++)
	{
		nodes_[i]->Tec_flag = true;
	}

	// MASS ERROR
	//!< Calculate nodal mass
	double rho = ElementMaterial_->density_0;
	double Jacobi, SHP[4];
	double B[6][12]; // No use here

	StrainMatrix(B, Jacobi);
	double Gmass = rho *1 * abs(Jacobi); // we can get volumn by integrating the function: f(x)=1
	
	SHPFunction(SHP);
	for (int i = 0; i < 4; i++)
	{
		nodes_[i]->mass += Gmass * SHP[i];
	}



	return true;
}


//	Write element data to stream
void CT4::Write(COutputter& output, unsigned int Ele)
{
	output << setw(5) << Ele + 1 << setw(11) << nodes_[0]->NodeNumber
		<< setw(9) << nodes_[1]->NodeNumber << setw(9) << nodes_[2]->NodeNumber << setw(9) << nodes_[3]->NodeNumber
		<< setw(12) << ElementMaterial_->nset << endl;
}


//	Write element data to stream
void CT4::Write(CTecOutputter& output)
{
	output << nodes_[0]->NodeTecNumber << " " << nodes_[1]->NodeTecNumber << " " << nodes_[2]->NodeTecNumber << " " << nodes_[2]->NodeTecNumber << " "
		<< nodes_[3]->NodeTecNumber << " " << nodes_[3]->NodeTecNumber << " " << nodes_[3]->NodeTecNumber << " " << nodes_[3]->NodeTecNumber << endl;
}


//  Generate location matrix: the global equation number that corresponding to each DOF of the element
//	Caution:  Equation number is numbered from 1 !
void CT4::GenerateLocationMatrix()
{
	unsigned int i = 0;
	for (unsigned int N = 0; N < NEN_; N++)
		for (unsigned int D = 0; D < 3; D++)
			LocationMatrix_[i++] = nodes_[N]->bcode[D];
}


//	Return the size of the element stiffness matrix (stored as an array column by column)
//	For 4 node T4 element, element stiffness is a 12x12 matrix, whose upper triangular part
//	has 150 elements
unsigned int CT4::SizeOfStiffnessMatrix() { return 150; }


//	Calculate element stiffness matrix 
//	Upper triangular matrix, stored as an array column by colum starting from the diagonal element
void CT4::ElementStiffness(double* Matrix)
{
	clear(Matrix, SizeOfStiffnessMatrix());


	


	//! Stiffness Matrix
	//! loop through the guass points to calculate strain matrix
	double B[6][12], BTMP[6][12];	//!< BTMP = D*B, K = B'*D*B = B'*BTMP
	double D[6][6];	//!< Constitutive relation
	double Jacobi;
	Constitutive(D);
	
			
	StrainMatrix(B, Jacobi);
	for (int L = 0; L < 12; L++)
		{
			BTMP[0][L] = D[0][0] * B[0][L] + D[0][1] * B[1][L] + D[0][2] * B[2][L] + D[0][3] * B[3][L] + D[0][4] * B[4][L] + D[0][5] * B[5][L];
			BTMP[1][L] = D[1][0] * B[0][L] + D[1][1] * B[1][L] + D[1][2] * B[2][L] + D[1][3] * B[3][L] + D[1][4] * B[4][L] + D[1][5] * B[5][L];
			BTMP[2][L] = D[2][0] * B[0][L] + D[2][1] * B[1][L] + D[2][2] * B[2][L] + D[2][3] * B[3][L] + D[2][4] * B[4][L] + D[2][5] * B[5][L];
			BTMP[3][L] = D[3][0] * B[0][L] + D[3][1] * B[1][L] + D[3][2] * B[2][L] + D[3][3] * B[3][L] + D[3][4] * B[4][L] + D[3][5] * B[5][L];
			BTMP[4][L] = D[4][0] * B[0][L] + D[4][1] * B[1][L] + D[4][2] * B[2][L] + D[4][3] * B[3][L] + D[4][4] * B[4][L] + D[4][5] * B[5][L];
			BTMP[5][L] = D[5][0] * B[0][L] + D[5][1] * B[1][L] + D[5][2] * B[2][L] + D[5][3] * B[3][L] + D[5][4] * B[4][L] + D[5][5] * B[5][L];
		}

		//! Stiffness Matrix
		for (int L = 0; L < 12; L++)
		{
			unsigned int Diag = L * (L + 1) / 2 + 1;
			for (int M = 0; M <= L; M++)
			{//error
			Matrix[Diag + L - M - 1] += 1 * abs(Jacobi) * (B[0][M] * BTMP[0][L] + B[1][M] * BTMP[1][L] + B[2][M] * BTMP[2][L] + B[3][M] * BTMP[3][L] + B[4][M] * BTMP[4][L] + B[5][M] * BTMP[5][L]);
			}
		}
}


//	Calculate element stress and nodal stress
//ERROR
void CT4::ElementStress(double* stress, double* Displacement)
{
	clear(stress, 6);
	
	double B[6][12], C[6][12];	//!< C= D*B
	double D[6][6];
	double Jacob;
	// SINGLE gAUSS ERROR
	
	Constitutive(D);
	StrainMatrix(B,  Jacob);
	for (int k = 0; k < 12; k++)
	{
		C[0][k] = D[0][0] * B[0][k] + D[0][1] * B[1][k] + D[0][2] * B[2][k] + D[0][3] * B[3][k] + D[0][4] * B[4][k] + D[0][5] * B[5][k];
		C[1][k] = D[1][0] * B[0][k] + D[1][1] * B[1][k] + D[1][2] * B[2][k] + D[1][3] * B[3][k] + D[1][4] * B[4][k] + D[1][5] * B[5][k];
		C[2][k] = D[2][0] * B[0][k] + D[2][1] * B[1][k] + D[2][2] * B[2][k] + D[2][3] * B[3][k] + D[2][4] * B[4][k] + D[2][5] * B[5][k];
		C[3][k] = D[3][0] * B[0][k] + D[3][1] * B[1][k] + D[3][2] * B[2][k] + D[3][3] * B[3][k] + D[3][4] * B[4][k] + D[3][5] * B[5][k];
		C[4][k] = D[4][0] * B[0][k] + D[4][1] * B[1][k] + D[4][2] * B[2][k] + D[4][3] * B[3][k] + D[4][4] * B[4][k] + D[4][5] * B[5][k];
		C[5][k] = D[5][0] * B[0][k] + D[5][1] * B[1][k] + D[5][2] * B[2][k] + D[5][3] * B[3][k] + D[5][4] * B[4][k] + D[5][5] * B[5][k];
	}

	for (int k = 0; k < 12; k++)
	{
		if (LocationMatrix_[k])
		{
			stress[0] += C[0][k] * Displacement[LocationMatrix_[k] - 1];
			stress[1] += C[1][k] * Displacement[LocationMatrix_[k] - 1];
			stress[2] += C[2][k] * Displacement[LocationMatrix_[k] - 1];
			stress[3] += C[3][k] * Displacement[LocationMatrix_[k] - 1];
			stress[4] += C[4][k] * Displacement[LocationMatrix_[k] - 1];
			stress[5] += C[5][k] * Displacement[LocationMatrix_[k] - 1];
		}
	}
			

	double str[6];
	str[0] = stress[0];
	str[1] = stress[1];
	str[2] = stress[2];
	str[3] = stress[3];
	str[4] = stress[4];
	str[5] = stress[5];

	for (unsigned int N = 0; N < NEN_; N++)
	{
		//error
		nodes_[N]->count += 1;
		nodes_[N]->Stress[0] += str[0];		// stress xx
		nodes_[N]->Stress[1] += str[1];		// stress yy
		nodes_[N]->Stress[2] += str[2];		// stress zz
		nodes_[N]->Stress[3] += str[3];		// stress xy
		nodes_[N]->Stress[4] += str[4];		// stress xz
		nodes_[N]->Stress[5] += str[5];		// stress yz
	}
}

//! Return the constitutive relation matrix of plain strain
void CT4::Constitutive(double D[6][6])
{
	CT4Material* mat = dynamic_cast<CT4Material*>(ElementMaterial_);
	double F = mat->E / (1 + mat->Poisson);
	double G = F * (mat->Poisson) / (1 - 2 * mat->Poisson);
	double H = F + G;

	D[0][0] = H;
	D[0][1] = G;
	D[0][2] = G;
	D[0][3] = 0;
	D[0][4] = 0;
    D[0][5] = 0;
	D[1][0] = G;
	D[1][1] = H;
	D[1][2] = G;
	D[1][3] = 0;
	D[1][4] = 0;
	D[1][5] = 0;
	D[2][0] = G;
	D[2][1] = G;
	D[2][2] = H;
	D[2][3] = 0;
	D[2][4] = 0;
	D[2][5] = 0;
	D[3][0] = 0;
	D[3][1] = 0;
	D[3][2] = 0;
	D[3][3] = F / 2.0;
	D[3][4] = 0;
	D[3][5] = 0;
	D[4][0] = 0;
	D[4][1] = 0;
	D[4][2] = 0;
	D[4][3] = 0;
	D[4][4] = F / 2.0;
	D[4][5] = 0;
	D[5][0] = 0;
	D[5][1] = 0;
	D[5][2] = 0;
	D[5][3] = 0;
	D[5][4] = 0;
	D[5][5] = F / 2.0;
}

//! Return the shape function value of point with parent coordinate (R, S, T)
void CT4::SHPFunction(double SHP[4])
{
	SHP[0] = 0.25;
	SHP[1] = 0.25;
	SHP[2] = 0.25;
	SHP[3] = 0.25;
}


//! Return the strain matrix value of point with parent coordinate (R, S, T)
void CT4::StrainMatrix(double B[6][12], double& Jacob)
{
	double P[3][4];
	//With respect to R
	P[0][0] = 1;
	P[0][1] = 0;
	P[0][2] = 0;
	P[0][3] = -1;
	//With respect to S
	P[1][0] = 0;
	P[1][1] = 1;
	P[1][2] = 0;
	P[1][3] = -1;
	//With respect to T
	P[2][0] = 0;
	P[2][1] = 0;
	P[2][2] = 1;
	P[2][3] = -1;
	

	//Evaluate the jacobi matrix at point (R,S)
	double XJ[3][3];
	double XJI[3][3];
	double DUM;
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			DUM = 0;
			for (int k = 0; k < 4; k++)
			{
				DUM = DUM + P[i][k] * (nodes_[k]->XYZ[j]);
				XJ[i][j] = DUM;
			}
		}
	}


	//Caculate the determinant of the jacobi matrix at point (R,S,T)
	Jacob = XJ[0][0] * (XJ[1][1] * XJ[2][2] - XJ[1][2] * XJ[2][1]) - XJ[0][1] * (XJ[1][0] * XJ[2][2] - XJ[1][2] * XJ[2][0]) + XJ[0][2] * (XJ[1][0] * XJ[2][1] - XJ[1][1] * XJ[2][0]);


	if (Jacob == 0)
	{
		cout << "The Jacobi of an element equals zeros!" << endl;
		exit(1);
	}


	//Cacaulate inverse of the jacobi matrix
	DUM = 1 / Jacob;
	XJI[0][0] = (XJ[1][1] * XJ[2][2] - XJ[1][2] * XJ[2][1]) * DUM;
	XJI[0][1] = -(XJ[0][1] * XJ[2][2] - XJ[2][1] * XJ[0][2]) * DUM;
	XJI[0][2] = (XJ[0][1] * XJ[1][2] - XJ[1][1] * XJ[0][2]) * DUM;
	XJI[1][0] = -(XJ[1][0] * XJ[2][2] - XJ[1][2] * XJ[2][0]) * DUM;
	XJI[1][1] = (XJ[0][0] * XJ[2][2] - XJ[2][0] * XJ[0][2]) * DUM;
	XJI[1][2] = -(XJ[0][0] * XJ[1][2] - XJ[1][0] * XJ[0][2]) * DUM;
	XJI[2][0] = (XJ[1][0] * XJ[2][1] - XJ[1][1] * XJ[2][0]) * DUM;
	XJI[2][1] = -(XJ[0][0] * XJ[2][1] - XJ[0][1] * XJ[2][0]) * DUM;
	XJI[2][2] = (XJ[0][0] * XJ[1][1] - XJ[0][1] * XJ[1][0]) * DUM;


	//Caculate global derivative operator B
	int K2 = 2;
	for (int k = 0; k <4; k++)
	{
		B[0][K2 - 2] = 0;
		B[0][K2 - 1] = 0;
		B[0][K2] = 0;
		B[1][K2 - 2] = 0;
		B[1][K2 - 1] = 0;
     	B[1][K2] = 0;
		B[2][K2 - 2] = 0;
		B[2][K2 - 1] = 0;
		B[2][K2] = 0;
		B[3][K2 - 2] = 0;
		B[3][K2 - 1] = 0;
		B[3][K2] = 0;
		B[4][K2 - 2] = 0;
		B[4][K2 - 1] = 0;
		B[4][K2] = 0;
		B[5][K2 - 2] = 0;
		B[5][K2 - 1] = 0;
		B[5][K2] = 0;
		for (int i = 0; i < 3; i++)
		{
			B[0][K2 - 2] = B[0][K2 - 2] + XJI[0][i] * P[i][k];
			B[1][K2 - 1] = B[1][K2 - 1] + XJI[1][i] * P[i][k];
			B[2][K2] = B[2][K2] + XJI[2][i] * P[i][k];
		}
		B[3][K2 - 2] = B[1][K2 - 1];
		B[3][K2 - 1] = B[0][K2 - 2];
		B[4][K2 - 2] = B[2][K2];
		B[4][K2] = B[0][K2 - 2];
		B[5][K2 - 1] = B[2][K2];
		B[5][K2] = B[1][K2 - 1];
		K2 = K2 + 3;
	}
}