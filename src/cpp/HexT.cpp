/*****************************************************************************/
/*  Quad: Element Quadratic for Plain Problem                                */
/*     Added by Ruichen Ni, 2018311066                                       */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Date: March 16, 2019                                                  */
/*****************************************************************************/

#include "HexT.h"

#include <iostream>
#include <iomanip>
#include <cmath>

//	Constructor
CHexT::CHexT()
{
	NEN_ = 20;	// Each element has 20 nodes
	nodes_ = new CNode*[NEN_];

	ND_ = 60;	// Each node has 3 dimension
	LocationMatrix_ = new unsigned int[ND_];

	ElementMaterial_ = nullptr;
}

//	Deconstructor
CHexT::~CHexT()
{
	if (nodes_)
		delete[] nodes_;

	if (LocationMatrix_)
		delete[] LocationMatrix_;
}

//	Read element data from stream Input
bool CHexT::Read(ifstream& Input, unsigned int Ele, CMaterial* MaterialSets, CNode* NodeList)
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
	unsigned int N1, N2, N3, N4, N5, N6, N7, N8, N9, N10, N11, N12, N13, N14, N15, N16, N17, N18, N19, N20;	// Twenty nodes number

	Input >> N1 >> N2 >> N3 >> N4 >> N5 >> N6 >> N7 >> N8 >> N9 >> N10 >> N11 >> N12 >> N13 >> N14 >> N15 >> N16 >> N17 >> N18 >> N19 >> N20 >> MSet;
	ElementMaterial_ = dynamic_cast<CHexTMaterial*>(MaterialSets) + MSet - 1;
	nodes_[0] = &NodeList[N1 - 1];
	nodes_[1] = &NodeList[N2 - 1];
	nodes_[2] = &NodeList[N3 - 1];
	nodes_[3] = &NodeList[N4 - 1];
	nodes_[4] = &NodeList[N5 - 1];
	nodes_[5] = &NodeList[N6 - 1];
	nodes_[6] = &NodeList[N7 - 1];
	nodes_[7] = &NodeList[N8 - 1];
	nodes_[8] = &NodeList[N9 - 1];
	nodes_[9] = &NodeList[N10 - 1];
	nodes_[10] = &NodeList[N11 - 1];
	nodes_[11] = &NodeList[N12 - 1];
	nodes_[12] = &NodeList[N13 - 1];
	nodes_[13] = &NodeList[N14 - 1];
	nodes_[14] = &NodeList[N15 - 1];
	nodes_[15] = &NodeList[N16 - 1];
	nodes_[16] = &NodeList[N17 - 1];
	nodes_[17] = &NodeList[N18 - 1];
	nodes_[18] = &NodeList[N19 - 1];
	nodes_[19] = &NodeList[N20 - 1];

	for (int i = 0; i < 8; i++)
	{
		nodes_[i]->Tec_flag = true;
	}

	for (int i = 8; i < NEN_; i++)
	{
		nodes_[i]->Tec_flag = false;
	}

	return true;
}

//	Write element data to stream
void CHexT::Write(COutputter& output, unsigned int Ele)
{
	output << setw(5) << Ele + 1 << setw(11) << nodes_[0]->NodeNumber
		<< setw(9) << nodes_[1]->NodeNumber << setw(9) << nodes_[2]->NodeNumber << setw(9) << nodes_[3]->NodeNumber
		<< setw(9) << nodes_[4]->NodeNumber << setw(9) << nodes_[5]->NodeNumber << setw(9) << nodes_[6]->NodeNumber
		<< setw(9) << nodes_[7]->NodeNumber << setw(9) << nodes_[8]->NodeNumber << setw(9) << nodes_[9]->NodeNumber
		<< setw(9) << nodes_[10]->NodeNumber << setw(9) << nodes_[11]->NodeNumber << setw(9) << nodes_[12]->NodeNumber
		<< setw(9) << nodes_[13]->NodeNumber << setw(9) << nodes_[14]->NodeNumber << setw(9) << nodes_[15]->NodeNumber
		<< setw(9) << nodes_[16]->NodeNumber << setw(9) << nodes_[17]->NodeNumber << setw(9) << nodes_[18]->NodeNumber
		<< setw(9) << nodes_[19]->NodeNumber
		<< setw(12) << ElementMaterial_->nset << endl;
}

//	Write element data to stream
void CHexT::Write(CTecOutputter& output)
{
	output << nodes_[0]->NodeTecNumber << " " << nodes_[1]->NodeTecNumber << " " << nodes_[2]->NodeTecNumber << " " << nodes_[3]->NodeTecNumber << " "<< nodes_[4]->NodeTecNumber << " " 
		<< nodes_[5]->NodeTecNumber << " " << nodes_[6]->NodeTecNumber << " " << nodes_[7]->NodeTecNumber << endl;
}

//  Generate location matrix: the global equation number that corresponding to each DOF of the element
//	Caution:  Equation number is numbered from 1 !
void CHexT::GenerateLocationMatrix()
{
	unsigned int i = 0;
	for (unsigned int N = 0; N < NEN_; N++)
		for (unsigned int D = 0; D < 3; D++)
			LocationMatrix_[i++] = nodes_[N]->bcode[D];
}

//	Return the size of the element stiffness matrix (stored as an array column by column)
//	For 20 node HexT20 element, element stiffness is a 60x60 matrix, whose upper triangular part
//	has 1830 elements
unsigned int CHexT::SizeOfStiffnessMatrix() { return 1830; }

//	Calculate element stiffness matrix 
//	Upper triangular matrix, stored as an array column by colum starting from the diagonal element
void CHexT::ElementStiffness(double* Matrix)
{
	clear(Matrix, SizeOfStiffnessMatrix());

	//! two point Guass Quadrature
	double GP[2], weight[2];
	Guassian(2, GP, weight);

	//! Stiffness Matrix
	//! loop through the guass points to calculate strain matrix
	double B[6][60], BTMP[6][60];	//!< BTMP = D*B, K = B'*D*B = B'*BTMP
	double D[6][6];	//!< Constitutive relation
	double Jacobi;
	Constitutive(D);
	for (int I = 0; I < 2; I++)
	{
		for (int J = 0; J < 2; J++)
			for (int K = 0; K < 2; K++)
			{
				StrainMatrix(B, GP[I], GP[J], GP[K], Jacobi);
				for (int L = 0; L < 60; L++)
				{
					BTMP[0][L] = D[0][0] * B[0][L] + D[0][1] * B[1][L] + D[0][2] * B[2][L] + D[0][3] * B[3][L] + D[0][4] * B[4][L] + D[0][5] * B[5][L];
					BTMP[1][L] = D[1][0] * B[0][L] + D[1][1] * B[1][L] + D[1][2] * B[2][L] + D[1][3] * B[3][L] + D[1][4] * B[4][L] + D[1][5] * B[5][L];
					BTMP[2][L] = D[2][0] * B[0][L] + D[2][1] * B[1][L] + D[2][2] * B[2][L] + D[2][3] * B[3][L] + D[2][4] * B[4][L] + D[2][5] * B[5][L];
					BTMP[3][L] = D[3][0] * B[0][L] + D[3][1] * B[1][L] + D[3][2] * B[2][L] + D[3][3] * B[3][L] + D[3][4] * B[4][L] + D[3][5] * B[5][L];
					BTMP[4][L] = D[4][0] * B[0][L] + D[4][1] * B[1][L] + D[4][2] * B[2][L] + D[4][3] * B[3][L] + D[4][4] * B[4][L] + D[4][5] * B[5][L];
					BTMP[5][L] = D[5][0] * B[0][L] + D[5][1] * B[1][L] + D[5][2] * B[2][L] + D[5][3] * B[3][L] + D[5][4] * B[4][L] + D[5][5] * B[5][L];
				}

				//! Stiffness Matrix
				for (int L = 0; L < 60; L++)
				{
					unsigned int Diag = L * (L + 1) / 2 + 1;
					for (int M = 0; M <= L; M++)
					{
						Matrix[Diag + L - M - 1] += weight[I] * weight[J] * weight[K] * abs(Jacobi)*(B[0][M] * BTMP[0][L] + B[1][M] * BTMP[1][L] + B[2][M] * BTMP[2][L]+ B[3][M] * BTMP[3][L] + B[4][M] * BTMP[4][L] + B[5][M] * BTMP[5][L]);
					}
				}
			}

	}
}

//	Calculate element stress and nodal stress
void CHexT::ElementStress(double* stress, double* Displacement)
{
	clear(stress, 48);
	double GP[2], weight[2];
	Guassian(2, GP, weight);
	double B[6][60], C[6][60];	//!< C= D*B
	double D[6][6];
	double Jacob;
	for (int i = 0; i < 2; i++)
	{
		for (int j = 0; j < 2; j++)
		{
			for (int h = 0; h < 2; h++)
			{
				Constitutive(D);
				StrainMatrix(B, GP[i], GP[j], GP[h], Jacob);
				for (int k = 0; k < 60; k++)
				{
					C[0][k] = D[0][0] * B[0][k] + D[0][1] * B[1][k] + D[0][2] * B[2][k] + D[0][3] * B[3][k] + D[0][4] * B[4][k] + D[0][5] * B[5][k];
					C[1][k] = D[1][0] * B[0][k] + D[1][1] * B[1][k] + D[1][2] * B[2][k] + D[1][3] * B[3][k] + D[1][4] * B[4][k] + D[1][5] * B[5][k];
					C[2][k] = D[2][0] * B[0][k] + D[2][1] * B[1][k] + D[2][2] * B[2][k] + D[2][3] * B[3][k] + D[2][4] * B[4][k] + D[2][5] * B[5][k];
					C[3][k] = D[3][0] * B[0][k] + D[3][1] * B[1][k] + D[3][2] * B[2][k] + D[3][3] * B[3][k] + D[3][4] * B[4][k] + D[3][5] * B[5][k];
					C[4][k] = D[4][0] * B[0][k] + D[4][1] * B[1][k] + D[4][2] * B[2][k] + D[4][3] * B[3][k] + D[4][4] * B[4][k] + D[4][5] * B[5][k];
					C[5][k] = D[5][0] * B[0][k] + D[5][1] * B[1][k] + D[5][2] * B[2][k] + D[5][3] * B[3][k] + D[5][4] * B[4][k] + D[5][5] * B[5][k];
				}

				for (int k = 0; k < 60; k++)
				{
					if (LocationMatrix_[k])
					{
						stress[((i * 2 + j) * 2 + h) * 6] += C[0][k] * Displacement[LocationMatrix_[k] - 1];
						stress[((i * 2 + j) * 2 + h) * 6 + 1] += C[1][k] * Displacement[LocationMatrix_[k] - 1];
						stress[((i * 2 + j) * 2 + h) * 6 + 2] += C[2][k] * Displacement[LocationMatrix_[k] - 1];
						stress[((i * 2 + j) * 2 + h) * 6 + 3] += C[3][k] * Displacement[LocationMatrix_[k] - 1];
						stress[((i * 2 + j) * 2 + h) * 6 + 4] += C[4][k] * Displacement[LocationMatrix_[k] - 1];
						stress[((i * 2 + j) * 2 + h) * 6 + 5] += C[5][k] * Displacement[LocationMatrix_[k] - 1];
					}
				}
			}
		}

	}

	double str[6];
	str[0] = stress[0] + stress[6] + stress[12] + stress[18] + stress[24] + stress[30] + stress[36] + stress[42];
	str[1] = stress[1] + stress[7] + stress[13] + stress[19] + stress[25] + stress[31] + stress[37] + stress[43];
	str[2] = stress[2] + stress[8] + stress[14] + stress[20] + stress[26] + stress[32] + stress[38] + stress[44];
	str[3] = stress[3] + stress[9] + stress[15] + stress[21] + stress[27] + stress[33] + stress[39] + stress[45];
	str[4] = stress[4] + stress[10] + stress[16] + stress[22] + stress[28] + stress[34] + stress[40] + stress[46];
	str[5] = stress[5] + stress[11] + stress[17] + stress[23] + stress[29] + stress[35] + stress[41] + stress[47];

	for (unsigned int N = 0; N < NEN_; N++)
	{
		nodes_[N]->count += 8;
		nodes_[N]->Stress[0] += str[0];		// stress xx
		nodes_[N]->Stress[1] += str[1];		// stress yy
		nodes_[N]->Stress[2] += str[2];		// stress zz
		nodes_[N]->Stress[3] += str[3];		// stress xy
		nodes_[N]->Stress[4] += str[4];		// stress xz
		nodes_[N]->Stress[5] += str[5];		// stress yz
	}
}

//! Return the constitutive relation matrix of plain strain
void CHexT::Constitutive(double D[6][6])
{
	CHexTMaterial* mat = dynamic_cast<CHexTMaterial*>(ElementMaterial_);
	double F = mat->E / (1 + mat->Poisson);
	double G = F*(mat->Poisson) / (1 - 2 * mat->Poisson);
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

//! Return the shape function value of point with parent coordinate (xi, eta)
void CHexT::SHPFunction(double SHP[20], double R, double S , double T)
{
	double R1, R2, S1, S2, T1, T2, N9, N10, N11, N12, N13, N14, N15, N16, N17, N18, N19, N20;
	double N1P, N2P, N3P, N4P, N5P, N6P, N7P, N8P;
	double N1, N2, N3, N4, N5, N6, N7, N8;
	R1 = 1 + R;
	R2 = 1 - R;
	S1 = 1 + S;
	S2 = 1 - S;
	T1 = 1 + T;
	T2 = 1 - T;
	N1P = 0.125*R2*S2*T2;
	N2P = 0.125*R1*S2*T2;
	N3P = 0.125*R1*S1*T2;
	N4P = 0.125*R2*S1*T2;
	N5P = 0.125*R2*S2*T1;
	N6P = 0.125*R1*S2*T1;
	N7P = 0.125*R1*S1*T1;
	N8P = 0.125*R2*S1*T1;
	SHP[8] = 0.25*T2*R1*R2*S2;
	SHP[9] = 0.25*T2*R1*S1*S2;
	SHP[10] = 0.25*T2*R1*R2*S1;
	SHP[11] = 0.25*T2*R2*S1*S2;
	SHP[12] = 0.25*T2*T1*R2*S2;
	SHP[13] = 0.25*T2*T1*R1*S2;
	SHP[14] = 0.25*T2*T1*R1*S1;
	SHP[15] = 0.25*T2*T1*R2*S1;
	SHP[16] = 0.25*T1*R1*R2*S2;
	SHP[17] = 0.25*T1*R1*S1*S2;
	SHP[18] = 0.25*T1*R1*R2*S1;
	SHP[19] = 0.25*T1*R2*S1*S2;
	SHP[0] = N1P - 0.5*N9 - 0.5*N12 - 0.5*N13;
	SHP[1] = N2P - 0.5*N9 - 0.5*N10 - 0.5*N14;
	SHP[2] = N3P - 0.5*N10 - 0.5*N11 - 0.5*N15;
	SHP[3] = N4P - 0.5*N11 - 0.5*N12 - 0.5*N16;
	SHP[4] = N5P - 0.5*N13 - 0.5*N17 - 0.5*N20;
	SHP[5] = N6P - 0.5*N17 - 0.5*N18 - 0.5*N14;
	SHP[6] = N7P - 0.5*N18 - 0.5*N19 - 0.5*N15;
	SHP[7] = N8P - 0.5*N19 - 0.5*N20 - 0.5*N16;
}

//! Return the strain matrix value of point with parent coordinate (R, S, T)
void CHexT::StrainMatrix(double B[6][60], double R, double S, double T, double& Jacob)
{
	double P[3][20];
	//With respect to R
	P[0][0] = -0.125*(1 - S)*(1 - T) - 0.125*(1 - R)*(1 - S)*(1 - T) +
		0.125 *(1 + R)*(1 - S)*(1 - T) + 0.125*(1 - S)*(1 + S)*(1 - T) +
		0.125 *(1 - S)*(1 - T)*(1 + T);
	P[0][1] = 0.125 *(1 - S) *(1 - T) - 0.125 *(1 - R)* (1 - S) *(1 - T) +
		0.125 *(1 + R) *(1 - S) *(1 - T) - 0.125* (1 - S)* (1 + S)* (1 - T) -
		0.125* (1 - S) *(1 - T) *(1 + T);
	P[0][2] = 0.125 *(1 + S) *(1 - T) - 0.125 *(1 - R) *(1 + S) *(1 - T) +
		0.125 *(1 + R) *(1 + S) *(1 - T) - 0.125 *(1 - S) *(1 + S) *(1 - T) -
		0.125 *(1 + S) *(1 - T) *(1 + T);
	P[0][3] = -0.125 *(1 + S) *(1 - T) - 0.125 *(1 - R) *(1 + S) *(1 - T) +
		0.125 *(1 + R) *(1 + S) *(1 - T) + 0.125 *(1 - S) *(1 + S) *(1 - T) +
		0.125 *(1 + S) *(1 - T) *(1 + T);
	P[0][4] = -0.125 *(1 - S) *(1 + T) - 0.125 *(1 - R) *(1 - S) *(1 + T) +
		0.125 *(1 + R) *(1 - S) *(1 + T) + 0.125 *(1 - S) *(1 + S) *(1 + T) +
		0.125 *(1 - S) *(1 - T) *(1 + T);
	P[0][5] = 0.125 *(1 - S) *(1 + T) - 0.125 *(1 - R) *(1 - S) *(1 + T) +
		0.125 *(1 + R) *(1 - S) *(1 + T) - 0.125 *(1 - S) *(1 + S) *(1 + T) -
		0.125 *(1 - S) *(1 - T) *(1 + T);
	P[0][6] = 0.125 *(1 + S) *(1 + T) - 0.125 *(1 - R) *(1 + S) *(1 + T) +
		0.125 *(1 + R) *(1 + S) *(1 + T) - 0.125 *(1 - S) *(1 + S) *(1 + T) -
		0.125 *(1 + S) *(1 - T) *(1 + T);
	P[0][7] = -0.125 *(1 + S) *(1 + T) - 0.125 *(1 - R) *(1 + S) *(1 + T) +
		0.125 *(1 + R) *(1 + S) *(1 + T) + 0.125 *(1 - S) *(1 + S) *(1 + T) +
		0.125 *(1 + S) *(1 - T) *(1 + T);
	P[0][8] = 0.25 *(1 - R) *(1 - S) *(1 - T) - 0.25 *(1 + R) *(1 - S) *(1 - T);
	P[0][9] = 0.25 *(1 - S) *(1 + S) *(1 - T);
	P[0][10] = 0.25 *(1 - R) *(1 + S) *(1 - T) - 0.25 *(1 + R)* (1 + S) *(1 - T);
	P[0][11] = -0.25 *(1 - S) *(1 + S)* (1 - T);
	P[0][12] = -0.25 *(1 - S) *(1 - T) *(1 + T);
	P[0][13] = 0.25 *(1 - S) *(1 - T) *(1 + T);
	P[0][14] = 0.25 *(1 + S) *(1 - T) *(1 + T);
	P[0][15] = -0.25 *(1 + S) *(1 - T) *(1 + T);
	P[0][16] = 0.25 *(1 - R) *(1 - S)* (1 + T) - 0.25 *(1 + R) *(1 - S) *(1 + T);
	P[0][17] = 0.25 *(1 - S) *(1 + S) *(1 + T);
	P[0][18] = 0.25 *(1 - R) *(1 + S) *(1 + T) - 0.25 *(1 + R) *(1 + S) *(1 + T);
	P[0][19] = -0.25 *(1 - S) *(1 + S) *(1 + T);

	//With respect to S
	P[1][0] = -0.125 *(1 - R) *(1 - T) + 0.125 *(1 - R) *(1 + R) *(1 - T) -
		0.125 *(1 - R) *(1 - S) *(1 - T) + 0.125 *(1 - R) *(1 + S) *(1 - T) +
		0.125 *(1 - R) *(1 - T) *(1 + T);
	P[1][1] = -0.125 *(1 + R) *(1 - T) + 0.125 *(1 - R) *(1 + R) *(1 - T) -
		0.125 *(1 + R) *(1 - S) *(1 - T) + 0.125 *(1 + R) *(1 + S) *(1 - T) +
		0.125 *(1 + R) *(1 - T) *(1 + T);
	P[1][2] = 0.125 *(1 + R) *(1 - T) - 0.125 *(1 - R) *(1 + R) *(1 - T) -
		0.125 *(1 + R) *(1 - S) *(1 - T) + 0.125 *(1 + R) *(1 + S) *(1 - T) -
		0.125 *(1 + R) *(1 - T) *(1 + T);
	P[1][3] = 0.125 *(1 - R) *(1 - T) - 0.125 *(1 - R) *(1 + R) *(1 - T) -
		0.125 *(1 - R) *(1 - S) *(1 - T) + 0.125 *(1 - R) *(1 + S) *(1 - T) -
		0.125 *(1 - R) *(1 - T) *(1 + T);
	P[1][4] = -0.125 *(1 - R) *(1 + T) + 0.125 *(1 - R) *(1 + R) *(1 + T) -
		0.125 *(1 - R) *(1 - S) *(1 + T) + 0.125 *(1 - R) *(1 + S) *(1 + T) +
		0.125 *(1 - R) *(1 - T) *(1 + T);
	P[1][5] = -0.125 *(1 + R) *(1 + T) + 0.125 *(1 - R) *(1 + R) *(1 + T) -
		0.125 *(1 + R) *(1 - S) *(1 + T) + 0.125 *(1 + R) *(1 + S) *(1 + T) +
		0.125 *(1 + R) *(1 - T) *(1 + T);
	P[1][6] = 0.125 *(1 + R) *(1 + T) - 0.125 *(1 - R) *(1 + R) *(1 + T) -
		0.125 *(1 + R) *(1 - S) *(1 + T) + 0.125 *(1 + R) *(1 + S) *(1 + T) -
		0.125 *(1 + R) *(1 - T) *(1 + T);
	P[1][7] = 0.125 *(1 - R) *(1 + T) - 0.125 *(1 - R) *(1 + R) *(1 + T) -
		0.125 *(1 - R) *(1 - S) *(1 + T) + 0.125 *(1 - R) *(1 + S) *(1 + T) -
		0.125 *(1 - R) *(1 - T) *(1 + T);
	P[1][8] = -0.25 *(1 - R) *(1 + R) *(1 - T);
	P[1][9] = 0.25 *(1 + R) *(1 - S) *(1 - T) - 0.25 *(1 + R) *(1 + S) *(1 - T);
	P[1][10] = 0.25 *(1 - R) *(1 + R)* (1 - T);
	P[1][11] = 0.25 *(1 - R) *(1 - S)* (1 - T) - 0.25 *(1 - R) *(1 + S) *(1 - T);
	P[1][12] = -0.25 *(1 - R) *(1 - T) *(1 + T);
	P[1][13] = -0.25 *(1 + R) *(1 - T) *(1 + T);
	P[1][14] = 0.25 *(1 + R) *(1 - T) *(1 + T);
	P[1][15] = 0.25 *(1 - R) *(1 - T) *(1 + T);
	P[1][16] = -0.25 *(1 - R) *(1 + R)* (1 + T);
	P[1][17] = 0.25 *(1 + R) *(1 - S) *(1 + T) - 0.25 *(1 + R) *(1 + S) *(1 + T);
	P[1][18] = 0.25 *(1 - R) *(1 + R) *(1 + T);
	P[1][19] = 0.25 *(1 - R) *(1 - S) *(1 + T) - 0.25 *(1 - R) *(1 + S) *(1 + T);

	//With respect to T
	P[2][0] = -0.125 *(1 - R) *(1 - S) + 0.125 *(1 - R) *(1 + R) *(1 - S) +
		0.125 *(1 - R) *(1 - S) *(1 + S) - 0.125 *(1 - R) *(1 - S) *(1 - T) +
		0.125 *(1 - R) *(1 - S) *(1 + T);
	P[2][1] = -0.125 *(1 + R) *(1 - S) + 0.125 *(1 - R) *(1 + R) *(1 - S) +
		0.125 *(1 + R) *(1 - S) *(1 + S) - 0.125 *(1 + R) *(1 - S) *(1 - T) +
		0.125 *(1 + R) *(1 - S) *(1 + T);
	P[2][2] = -0.125 *(1 + R) *(1 + S) + 0.125 *(1 - R) *(1 + R) *(1 + S) +
		0.125 *(1 + R) *(1 - S) *(1 + S) - 0.125 *(1 + R) *(1 + S) *(1 - T) +
		0.125 *(1 + R) *(1 + S) *(1 + T);
	P[2][3] = -0.125 *(1 - R) *(1 + S) + 0.125 *(1 - R) *(1 + R) *(1 + S) +
		0.125 *(1 - R) *(1 - S) *(1 + S) - 0.125 *(1 - R) *(1 + S) *(1 - T) +
		0.125 *(1 - R) *(1 + S) *(1 + T);
	P[2][4] = 0.125 *(1 - R) *(1 - S) - 0.125 *(1 - R) *(1 + R) *(1 - S) -
		0.125 *(1 - R) *(1 - S) *(1 + S) - 0.125 *(1 - R) *(1 - S) *(1 - T) +
		0.125 *(1 - R) *(1 - S) *(1 + T);
	P[2][5] = 0.125 *(1 + R)* (1 - S) - 0.125 *(1 - R) *(1 + R) *(1 - S) -
		0.125 *(1 + R) *(1 - S) *(1 + S) - 0.125 *(1 + R) *(1 - S) *(1 - T) +
		0.125 *(1 + R) *(1 - S) *(1 + T);
	P[2][6] = 0.125 *(1 + R) *(1 + S) - 0.125 *(1 - R) *(1 + R) *(1 + S) -
		0.125 *(1 + R) *(1 - S) *(1 + S) - 0.125 *(1 + R) *(1 + S) *(1 - T) +
		0.125 *(1 + R) *(1 + S) *(1 + T);
	P[2][7] = 0.125 *(1 - R) *(1 + S) - 0.125 *(1 - R) *(1 + R) *(1 + S) -
		0.125 *(1 - R) *(1 - S) *(1 + S) - 0.125 *(1 - R) *(1 + S) *(1 - T) +
		0.125 *(1 - R) *(1 + S) *(1 + T);
	P[2][8] = -0.25 *(1 - R) *(1 + R) *(1 - S);
	P[2][9] = -0.25 *(1 + R) *(1 - S) *(1 + S);
	P[2][10] = -0.25 *(1 - R) *(1 + R) *(1 + S);
	P[2][11] = -0.25 *(1 - R) *(1 - S) *(1 + S);
	P[2][12] = 0.25 *(1 - R) *(1 - S) *(1 - T) - 0.25 *(1 - R) *(1 - S) *(1 + T);
	P[2][13] = 0.25 *(1 + R) *(1 - S) *(1 - T) - 0.25 *(1 + R) *(1 - S) *(1 + T);
	P[2][14] = 0.25 *(1 + R) *(1 + S) *(1 - T) - 0.25 *(1 + R) *(1 + S) *(1 + T);
	P[2][15] = 0.25 *(1 - R) *(1 + S) *(1 - T) - 0.25 *(1 - R) *(1 + S) *(1 + T);
	P[2][16] = 0.25 *(1 - R) *(1 + R) *(1 - S);
	P[2][17] = 0.25 *(1 + R) *(1 - S) *(1 + S);
	P[2][18] = 0.25 *(1 - R) *(1 + R) *(1 + S);
	P[2][19] = 0.25 *(1 - R) *(1 - S) *(1 + S);

	//Evaluate the jacobi matrix at point (R,S)
	double XJ[3][3];
	double XJI[3][3];
	double DUM;
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			DUM = 0;
			for (int k = 0; k < 20; k++)
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
	for (int k = 0; k < 20; k++)
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