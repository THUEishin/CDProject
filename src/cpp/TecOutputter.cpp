#include "Domain.h"
#include "TecOutputter.h"
#include <iomanip>

using namespace std;

CTecOutputter* CTecOutputter::_instance = nullptr;

//	Constructor
CTecOutputter::CTecOutputter(string FileName)
{
	OutputFile.open(FileName);

	NUMTP = 0;

	if (!OutputFile)
	{
		cerr << "*** Error *** File " << FileName << " does not exist !" << endl;
		exit(3);
	}
}

//	Return the single instance of the class
CTecOutputter* CTecOutputter::Instance(string FileName)
{
	if (!_instance)
		_instance = new CTecOutputter(FileName);
	return _instance;
}

// 	Output title of Tecplot file
void CTecOutputter::OutputHeading()
{
	CDomain* FEMData = CDomain::Instance();

	*this << "TITLE = " << FEMData->GetTitle() << endl;
}

//	Print Initial data
void CTecOutputter::OutputInitInfo()
{
	CDomain* FEMData = CDomain::Instance();

	CNode* NodeList = FEMData->GetNodeList();

	CElementGroup* EleGrpList = FEMData->GetEleGrpList();

	if (FEMData->GetPTYPE()) //! 3D
	{
		*this << "VARIABLES = \"X\", \"Y\", \"Z\", \"UX\", \"UY\", \"UZ\", \"SXX\", \"SYY\", \"SZZ\", \"SXY\", \"SXZ\", \"SYZ\"" << endl;
	}
	else  //£¡2D
	{
		*this << "VARIABLES = \"X\", \"Y\", \"UX\", \"UY\", \"SXX\", \"SYY\", \"SXY\"" << endl;
	}

	unsigned int NUMNP = FEMData->GetNUMNP();
	unsigned int PTYPE = FEMData->GetPTYPE();
	unsigned int NUMEG = FEMData->GetNUMEG();
	unsigned int NUMET = 0;

	for (unsigned int N = 0; N < NUMEG; N++)
	{
		NUMET += EleGrpList[N].GetNUME();
	}
	
	for (unsigned int i = 0; i < NUMNP; i++)
	{
		if (NodeList[i].Tec_flag)
			NUMTP++;
	}

	*this << "ZONE T=\"P_Initial\", F=FEPOINT, N=" << NUMTP << ", E=" << NUMET;
	
	if (PTYPE)
	{
		*this << ", ET=BRICK" << endl;
	}
	else
	{
		*this << ", ET=QUADRILATERAL" << endl;
	}

	for (unsigned int np = 0; np < NUMNP; np++)
		NodeList[np].Write(*this, PTYPE);

	for (unsigned int neg = 0; neg < NUMEG; neg++)
	{
		CElementGroup& EleGrp = EleGrpList[neg];
		unsigned int NET = EleGrp.GetNUME();
		for (unsigned int ne = 0; ne < NET; ne++)
			EleGrp[ne].Write(*this);
	}
}

//	Print Result data
void CTecOutputter::OutputResult(unsigned int flag, unsigned int lcase)
{
	//! flag: 1-result of initial phase
	//!       2-result of deformed phase

	CDomain* FEMData = CDomain::Instance();

	CNode* NodeList = FEMData->GetNodeList();

	CElementGroup* EleGrpList = FEMData->GetEleGrpList();

	unsigned int NUMNP = FEMData->GetNUMNP();
	unsigned int PTYPE = FEMData->GetPTYPE();
	unsigned int NUMEG = FEMData->GetNUMEG();
	double* Displacement = FEMData->GetDisplacement();
	unsigned int NUMET = 0;

	for (unsigned int N = 0; N < NUMEG; N++)
	{
		NUMET += EleGrpList[N].GetNUME();
	}

	if (flag == 1)
	{
		*this << "ZONE T=\"P_lcase_" << to_string(lcase) << "_initphase\", ";
	}
	else
	{
		*this << "ZONE T=\"P_lcase_" << to_string(lcase) << "_deformphase\", ";
	}
	*this << "F=FEPOINT, N=" << NUMTP << ", E=" << NUMET;

	if (PTYPE)
	{
		*this << ", ET=BRICK, ";
	}
	else
	{
		*this << ", ET=QUADRILATERAL, ";
	}

	if (flag == 1)
	{
		if(PTYPE)
			*this << "D=(1,2,3,FECONNECT)" << endl;
		else
			*this << "D=(1,2,FECONNECT)" << endl;
	}
	else
	{
		*this << "D=(FECONNECT)" << endl;
	}

	for (unsigned int np = 0; np < NUMNP; np++)
		NodeList[np].Write(*this, PTYPE, flag, Displacement);
}