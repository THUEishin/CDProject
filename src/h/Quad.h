#pragma once
/*****************************************************************************/
/*  Quad: Element Quadratic for Plain Problem                                */
/*     Added by Ruichen Ni, 2018311066                                       */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Date: March 16, 2019                                                  */
/*****************************************************************************/

#pragma once

#include "Element.h"
using namespace std;

class CQuad : public CElement
{
public:

	//!	Constructor
	CQuad();

	//!	Desconstructor
	~CQuad();

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
	void Constitutive(double D[3][3]);

	//! Return the shape function value of point with parent coordinate (xi, eta)
	void SHPFunction(double SHP[4], double xi, double eta);

	//! Return the strain matrix value of point with parent coordinate (xi, eta)
	void StrainMatrix(double B[3][8], double xi, double eta, double& Jacob);
};