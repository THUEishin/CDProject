#pragma once
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

#include "Element.h"
using namespace std;

class CT4 : public CElement
{
public:

	//!	Constructor
	CT4();

	//!	Desconstructor
	~CT4();

	//!	Read element data from stream Input
	virtual bool Read(ifstream& Input, unsigned int Ele, CMaterial* MaterialSets, CNode* NodeList);

	//!	Write element data to stream
	virtual void Write(COutputter& output, unsigned int Ele);

	//!	Write element data to Tecplot stream
	virtual void Write(CTecOutputter& output);

	//! Generate location matrix: the global equation number that corresponding to each DOF of the element
	//	Caution:  Equation number is numbered from 1 !
	virtual void GenerateLocationMatrix();

	//!	Calculate element stiffness matrix
	virtual void ElementStiffness(double* Matrix);

	//!	Calculate element stress
	virtual void ElementStress(double* stress, double* Displacement);

	//!	Return the size of the element stiffness matrix (stored as an array column by column)
	virtual unsigned int SizeOfStiffnessMatrix();

	//! Return the constitutive relation matrix of plain strain or stress
	void Constitutive(double D[6][6]);

	//! Return the shape function value of point with parent coordinate (R,S,T)
	void SHPFunction(double SHP[4] );

	//! Return the strain matrix value of point with parent coordinate (R,S,T)
	void StrainMatrix(double B[6][12],  double& Jacob);
};