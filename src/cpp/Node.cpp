/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.11, November 22, 2017                                       */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#include <iostream>
#include <iomanip>

#include "Node.h"

CNode::CNode(double X, double Y, double Z)
{
    XYZ[0] = X;		// Coordinates of the node
    XYZ[1] = Y;
    XYZ[2] = Z;

	for (int i = 0; i < 6; i++)
		Stress[i] = 0.0;

	count = 0;
	Tec_flag = false;
	NodeTecNumber = 0;
    
    bcode[0] = 0;	// Boundary codes
    bcode[1] = 0;
    bcode[2] = 0;
};

//	Read element data from stream Input
bool CNode::Read(ifstream& Input, unsigned int np)
{
	unsigned int N;

	Input >> N;	// node number
	if (N != np + 1) 
	{
		cerr << "*** Error *** Nodes must be inputted in order !" << endl 
			 << "   Expected node number : " << np + 1 << endl
			 << "   Provided node number : " << N << endl;

		return false;
	}

	NodeNumber = N;

	Input >> bcode[0] >> bcode[1] >> bcode[2]
		  >> XYZ[0] >> XYZ[1] >> XYZ[2];

	return true;
}

//	Output nodal point data to stream
void CNode::Write(COutputter& output, unsigned int np)
{
	output << setw(9) << np + 1 << setw(5) << bcode[0] << setw(5) << bcode[1] << setw(5) << bcode[2]
		   << setw(18) << XYZ[0] << setw(15) << XYZ[1] << setw(15) << XYZ[2] << endl;
}

//	Output nodal point data to Tecplot stream
void CNode::Write(CTecOutputter& output, unsigned int PTYPE, unsigned int flag,double* Displacement)
{
	//! flag: 0-nodal position of Initial phase
	//!       1-nodal position of initial phase after calculation
	//!       2-nodal position of deformed phase after calculation
	//! Displacement = nullptr before calculation
	if (!Tec_flag) return;

	double UXYZ[NDF];
	for (unsigned int i = 0; i < NDF; i++) UXYZ[i] = 0.0;

	if (Displacement)
	{
		for (unsigned int j = 0; j < NDF; j++)
		{
			if (bcode[j] == 0)
			{
				UXYZ[j] = 0.0;
			}
			else
			{
				UXYZ[j] = Displacement[bcode[j] - 1];
			}
		}
	}

	if (PTYPE)
	{
		if(flag == 0)
			output << XYZ[0] << " " << XYZ[1] << " " << XYZ[2] <<" ";
		else if(flag == 2)
		{
			if (!Displacement)
			{
				cout << "***Error*** Can not write nodal position of deformed phase before calculation" << endl;
				exit(1);
			}
			
			output << XYZ[0] + UXYZ[0] << " " << XYZ[1] + UXYZ[1] << " " << XYZ[2] + UXYZ[2] << " ";
		}

		for (unsigned int i = 0; i < NDF; i++) output<< UXYZ[i] << " ";

		for (unsigned int i = 0; i < 6; i++) output << Stress[i] << " ";

		output << endl;
	}
	else
	{
		if (flag == 0)
			output << XYZ[0] << " " << XYZ[1] << " ";
		else if(flag == 2)
		{
			if (!Displacement)
			{
				cout << "***Error*** Can not write nodal position of deformed phase before calculation" << endl;
				exit(1);
			}

			output << XYZ[0] + UXYZ[0] << " " << XYZ[1] + UXYZ[1] << " ";
		}

		for (unsigned int i = 0; i < 2; i++) output << UXYZ[i] << " ";

		output << Stress[0] << " " << Stress[1] << " " << Stress[3] << " ";

		output << endl;
	}
}

//	Output equation numbers of nodal point to stream
void CNode::WriteEquationNo(COutputter& output, unsigned int np)
{
	output << setw(9) << np+1 << "       ";

	for (unsigned int dof = 0; dof < CNode::NDF; dof++)	// Loop over for DOFs of node np
	{
		output << setw(5) << bcode[dof];
	}

	output << endl;
}

//	Write nodal displacement
void CNode::WriteNodalDisplacement(COutputter& output, unsigned int np, double* Displacement)
{
	output << setw(5) << np + 1 << "        ";

	for (unsigned int j = 0; j < NDF; j++)
	{
		if (bcode[j] == 0)
		{
			output << setw(18) << 0.0;
		}
		else
		{
			output << setw(18) << Displacement[bcode[j] - 1];
		}
	}

	output << endl;
}

//  Set nodal stress to 0 for every LOAD CASE
void CNode::ResetNodalStress()
{
	for (int i = 0; i < 6; i++)
		Stress[i] = 0.0;

	count = 0;
}
