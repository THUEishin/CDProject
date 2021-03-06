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

#include "Outputter.h"

using namespace std;

//!	Material base class which only define one data member
/*!	All type of material classes should be derived from this base class */
class CMaterial
{
public:

	unsigned int nset;	//!< Number of set
	
	double E;  //!< Young's modulus

	double density_0;  //!< initial density 

public:

//! Virtual deconstructor
    virtual ~CMaterial() {density_0 = 0.0;};

//!	Read material data from stream Input
	virtual bool Read(ifstream& Input, unsigned int mset) = 0;

//!	Write material data to Stream
    virtual void Write(COutputter& output, unsigned int mset) = 0;

};

//!	Material class for bar element
class CBarMaterial : public CMaterial
{
public:

	double Area;	//!< Sectional area of a bar element

public:
	
//!	Read material data from stream Input
	virtual bool Read(ifstream& Input, unsigned int mset);

//!	Write material data to Stream
	virtual void Write(COutputter& output, unsigned int mset);
};

//! Material class for quadratic element
class CQuadMaterial : public CMaterial
{
public:

	double Poisson;		//!< Poisson rate of a quadratic element

public:

	//!	Read material data from stream Input
	virtual bool Read(ifstream& Input, unsigned int mset);

	//!	Write material data to Stream
	virtual void Write(COutputter& output, unsigned int mset);
};

//! Material class for HexT20 element
class CHexTMaterial : public CMaterial
{
public:

	double Poisson;		//!< Poisson rate of a quadratic element

public:

	//!	Read material data from stream Input
	virtual bool Read(ifstream& Input, unsigned int mset);

	//!	Write material data to Stream
	virtual void Write(COutputter& output, unsigned int mset);
};

class CH8Material : public CMaterial
{
public:

	double Poisson;		//!< Poisson rate of a quadratic element

public:

	//!	Read material data from stream Input
	virtual bool Read(ifstream& Input, unsigned int mset);

	//!	Write material data to Stream
	virtual void Write(COutputter& output, unsigned int mset);
};

class CT4Material : public CMaterial
{
public:

	double Poisson;		//!< Poisson rate of a quadratic element

public:

	//!	Read material data from stream Input
	virtual bool Read(ifstream& Input, unsigned int mset);

	//!	Write material data to Stream
	virtual void Write(COutputter& output, unsigned int mset);
};