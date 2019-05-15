/*****************************************************************************/
/*  Quad: Element Quadratic for Plain Problem                                */
/*     Added by Ruichen Ni, 2018311066                                       */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Date: March 16, 2019                                                  */
/*****************************************************************************/

#include "../h/Quad.h"

#include <iostream>
#include <iomanip>
#include <cmath>

//	Constructor
CQuad::CQuad()
{
	NEN_ = 4;	// Each element has 4 nodes
	nodes_ = new CNode*[NEN_];

	ND_ = 8;	// Each node has 2 dimension for plain problem
	LocationMatrix_ = new unsigned int[ND_];

	ElementMaterial_ = nullptr;
}

//	Deconstructor
CQuad::~CQuad()
{
	if (nodes_)
		delete[] nodes_;

	if (LocationMatrix_)
		delete[] LocationMatrix_;
}

//	Read element data from stream Input
bool CQuad::Read(ifstream& Input, unsigned int Ele, CMaterial* MaterialSets, CNode* NodeList)
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
	unsigned int N1, N2, N3, N4;	// Left node number and right node number

	Input >> N1 >> N2 >> N3 >> N4 >> MSet;
	ElementMaterial_ = dynamic_cast<CQuadMaterial*>(MaterialSets) + MSet - 1;
	nodes_[0] = &NodeList[N1 - 1];
	nodes_[1] = &NodeList[N2 - 1];
	nodes_[2] = &NodeList[N3 - 1];
	nodes_[3] = &NodeList[N4 - 1];

	for (int i = 0; i < 4; i++)
	{
		nodes_[i]->Tec_flag = true;
	}

	//!< Calculate nodal mass
	double rho = ElementMaterial_->density_0;
	//! two point Guass Quadrature
	double GP[2], weight[2], Jacobi, SHP[4];
	Guassian(2, GP, weight);
	double B[3][8]; // No use here

	for (int I = 0; I < 2; I++)
	{
		for (int J = 0; J < 2; J++)
		{
			StrainMatrix(B, GP[I], GP[J], Jacobi);
			double Gmass = rho * weight[I] * weight[J] * abs(Jacobi); // we can get volumn by integrating the function: f(x)=1
			// Expolate the mass from Gauss point to nodes
			SHPFunction(SHP, GP[I], GP[J]);
			for (int i = 0; i < 4; i++)
			{
				nodes_[i]->mass += Gmass * SHP[i];
			}
		}
	}

	return true;
}

//	Write element data to stream
void CQuad::Write(COutputter& output, unsigned int Ele)
{
	output << setw(5) << Ele + 1 << setw(11) << nodes_[0]->NodeNumber << setw(9) << nodes_[1]->NodeNumber
		<< setw(11) << nodes_[2]->NodeNumber << setw(9) << nodes_[3]->NodeNumber << setw(12) << ElementMaterial_->nset << endl;
}

//	Write element data to stream
void CQuad::Write(CTecOutputter& output)
{
	output << nodes_[0]->NodeTecNumber << " " << nodes_[1]->NodeTecNumber << " " << nodes_[2]->NodeTecNumber << " " << nodes_[3]->NodeTecNumber << endl;
}

//  Generate location matrix: the global equation number that corresponding to each DOF of the element
//	Caution:  Equation number is numbered from 1 !
void CQuad::GenerateLocationMatrix()
{
	unsigned int i = 0;
	for (unsigned int N = 0; N < NEN_; N++)
		for (unsigned int D = 0; D < 2; D++)	// 2-dimension and z-direction must be fixed
			LocationMatrix_[i++] = nodes_[N]->bcode[D];
}

//	Return the size of the element stiffness matrix (stored as an array column by column)
//	For 4 node quadratic element, element stiffness is a 8x8 matrix, whose upper triangular part
//	has 36 elements
unsigned int CQuad::SizeOfStiffnessMatrix() { return 36; }

//	Calculate element stiffness matrix 
//	Upper triangular matrix, stored as an array column by colum starting from the diagonal element
void CQuad::ElementStiffness(double* Matrix)
{
	clear(Matrix, SizeOfStiffnessMatrix());

	//! two point Guass Quadrature
	double GP[2], weight[2];
	Guassian(2, GP, weight);

	//! Stiffness Matrix
	//! loop through the guass points to calculate strain matrix
	double B[3][8], BTMP[3][8];	//!< BTMP = D*B, K = B'*D*B = B'*BTMP
	double D[3][3];	//!< Constitutive relation
	double Jacobi;
	Constitutive(D);
	for (int I = 0; I < 2; I++)
	{
		for (int J = 0; J < 2; J++)
		{
			StrainMatrix(B, GP[I], GP[J], Jacobi);
			for (int K = 0; K < 8; K++)
			{
				BTMP[0][K] = D[0][0] * B[0][K] + D[0][1] * B[1][K] + D[0][2] * B[2][K];
				BTMP[1][K] = D[1][0] * B[0][K] + D[1][1] * B[1][K] + D[1][2] * B[2][K];
				BTMP[2][K] = D[2][0] * B[0][K] + D[2][1] * B[1][K] + D[2][2] * B[2][K];
			}

			//! Stiffness Matrix
			for (int K = 0; K < 8; K++)
			{
				unsigned int Diag = K * (K + 1) / 2 + 1;
				for (int L = 0; L <= K; L++)
				{
					Matrix[Diag + K - L - 1] += weight[I] * weight[J] * abs(Jacobi)*(B[0][L] * BTMP[0][K] + B[1][L] * BTMP[1][K] + B[2][L] * BTMP[2][K]);
				}
			}
		}
	}
}

//	Calculate element stress and nodal stress
void CQuad::ElementStress(double* stress, double* Displacement)
{
	clear(stress, 12);
	double GP[2], weight[2];
	Guassian(2, GP, weight);
	double B[3][8], C[3][8];	//!< C= D*B
	double D[3][3];
	double Jacob;
	for (int i = 0; i < 2; i++)
	{
		for (int j = 0; j < 2; j++)
		{
			Constitutive(D);
			StrainMatrix(B, GP[i], GP[j], Jacob);
			for (int k = 0; k < 8; k++)
			{
				C[0][k] = D[0][0] * B[0][k] + D[0][1] * B[1][k] + D[0][2] * B[2][k];
				C[1][k] = D[1][0] * B[0][k] + D[1][1] * B[1][k] + D[1][2] * B[2][k];
				C[2][k] = D[2][0] * B[0][k] + D[2][1] * B[1][k] + D[2][2] * B[2][k];
			}

			for (int k = 0; k < 8; k++)
			{
				if (LocationMatrix_[k])
				{
					stress[(i * 2 + j) * 3] += C[0][k] * Displacement[LocationMatrix_[k] - 1];
					stress[(i * 2 + j) * 3 + 1] += C[1][k] * Displacement[LocationMatrix_[k] - 1];
					stress[(i * 2 + j) * 3 + 2] += C[2][k] * Displacement[LocationMatrix_[k] - 1];
				}
			}
		}
	}

	double str[3];
	str[0] = stress[0] + stress[3] + stress[6] + stress[9];		// stress xx
	str[1] = stress[1] + stress[4] + stress[7] + stress[10];	// stress yy
	str[2] = stress[2] + stress[5] + stress[8] + stress[11];	// stress xy

	for (unsigned int N = 0; N < NEN_; N++)
	{
		nodes_[N]->count += 4;				// 4 is the number of Gauss Point in an element
		nodes_[N]->Stress[0] += str[0];		// stress xx
		nodes_[N]->Stress[1] += str[1];		// stress yy
		nodes_[N]->Stress[3] += str[3];		// stress xy
	}
}

//! Return the constitutive relation matrix of plain strain
void CQuad::Constitutive(double D[3][3])
{
	CQuadMaterial* mat = dynamic_cast<CQuadMaterial*>(ElementMaterial_);
	double temp;
	temp = mat->E / ((1 + mat->Poisson)*(1 - 2 * mat->Poisson));
	D[0][0] = temp * (1 - mat->Poisson);
	D[0][1] = temp * mat->Poisson;
	D[0][2] = 0;
	D[1][0] = temp * mat->Poisson;
	D[1][1] = temp * (1 - mat->Poisson);
	D[1][2] = 0;
	D[2][0] = 0;
	D[2][1] = 0;
	D[2][2] = temp * (1 - 2 * mat->Poisson) / 2.0;
}

//! Return the shape function value of point with parent coordinate (xi, eta)
void CQuad::SHPFunction(double SHP[4], double xi, double eta)
{
	SHP[0] = (1 - xi)*(1 - eta) / 4.0;
	SHP[1] = (1 + xi)*(1 - eta) / 4.0;
	SHP[2] = (1 + xi)*(1 + eta) / 4.0;
	SHP[3] = (1 - xi)*(1 + eta) / 4.0;
}

//! Return the strain matrix value of point with parent coordinate (xi, eta)
void CQuad::StrainMatrix(double B[3][8], double xi, double eta, double& Jacob)
{
	double GN[2][4];
	GN[0][0] = (eta - 1) / 4.0;
	GN[1][0] = (xi - 1) / 4.0;
	GN[0][1] = (1 - eta) / 4.0;
	GN[1][1] = (-xi - 1) / 4.0;
	GN[0][2] = (1 + eta) / 4.0;
	GN[1][2] = (1 + xi) / 4.0;
	GN[0][3] = (-eta - 1) / 4.0;
	GN[1][3] = (1 - xi) / 4.0;

	double J[2][2];
	J[0][0] = 0; J[0][1] = 0; J[1][0] = 0; J[1][1] = 0;

	for (int I = 0; I < 4; I++)
	{
		J[0][0] += GN[0][I] * nodes_[I]->XYZ[0];
		J[0][1] += GN[0][I] * nodes_[I]->XYZ[1];
		J[1][0] += GN[1][I] * nodes_[I]->XYZ[0];
		J[1][1] += GN[1][I] * nodes_[I]->XYZ[1];
	}

	Jacob = J[0][0] * J[1][1] - J[0][1] * J[1][0];
	if (Jacob == 0)
	{
		cout << "The Jacobi of an element equals zeros!" << endl;
		exit(1);
	}
	double temp = J[0][0] / Jacob;
	J[0][0] = J[1][1] / Jacob;
	J[1][1] = temp;
	J[0][1] = -J[0][1] / Jacob;
	J[1][0] = -J[1][0] / Jacob;

	double DPHI[2][4];
	for (int I = 0; I < 4; I++)
	{
		DPHI[0][I] = J[0][0] * GN[0][I] + J[0][1] * GN[1][I];
		DPHI[1][I] = J[1][0] * GN[0][I] + J[1][1] * GN[1][I];
	}

	for (int I = 0; I < 4; I++)
	{
		B[0][2 * I] = DPHI[0][I];
		B[0][2 * I + 1] = 0;
		B[1][2 * I] = 0;
		B[1][2 * I + 1] = DPHI[1][I];
		B[2][2 * I] = DPHI[1][I];
		B[2][2 * I + 1] = DPHI[0][I];
	}
}