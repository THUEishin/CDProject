#pragma once
/*****************************************************************************/
/*  TecOutputter: postprocess phase of writting results to Tecplot file      */
/*     Added by Ruichen Ni, 2018311066                                       */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Date: April 25, 2019                                                  */
/*****************************************************************************/

#include <fstream>
#include <iostream>

using namespace std;

//! TecOutputer class is used to output results to Tecplot file for visualization
class CTecOutputter
{
private:

	//!	File stream for output
	ofstream OutputFile;

protected:

	//!	Constructor
	CTecOutputter(string FileName);

	//!	Designed as a single instance class
	static CTecOutputter* _instance;

public:

	//!	Return pointer to the output file stream
	inline ofstream* GetOutputFile() { return &OutputFile; }

	//!	Return the single instance of the class
	static CTecOutputter* Instance(string FileName = " ");

	//!	Output title of Tecplot file and Variables' name
	void OutputHeading();

	//!	Output nodal point initial position data
	void OutputInitInfo();

	//!	Output nodal point initial position data
	void OutputResult(unsigned int flag, unsigned int lcase);

	//! Overload the operator <<
	template <typename T>
	CTecOutputter& operator<<(const T& item)
	{
		OutputFile << item;
		return *this;
	}

	typedef std::basic_ostream<char, std::char_traits<char> > CharOstream;
	CTecOutputter& operator<<(CharOstream& (*op)(CharOstream&))
	{
		op(OutputFile);
		return *this;
	}
};