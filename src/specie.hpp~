#ifndef specie_HPP_
#define specie_HPP_

#include "specie.hpp"
#include <string>
#include "gmm/gmm.h"
#include <iostream>
#include <string>
#include <iostream>
#include <fstream>
#include <vector>

using namespace gmm;
using namespace std;

//definisco quelle standard, e se sono fluide, e se sono gas

const string allnames[]={"KER1","KER2","KER3","KER4","KER5","KER6","KER7","KER8", "RES","OLIO", "CHX", "CH4"};
const int allisfluid[]={0,0,0,0,0,0,0,0,0,1,1,1};
const int allisgas[]={0,0,0,0,0,0,0,0,0,0,1,1};


class Species
{
public:


	//! @name Constructor & Destructor
	//@{
		
	//! Empty constructor
	Species();

	explicit Species(string);
	
	//! Destructor
	virtual ~Species();
	
	inline string nome() const {return M_nome;}

	inline bool isfluid() const {return M_isfluid;}

	inline void setCiniSp(double val) {M_Cini=val;}

	inline double getCiniSp() const {return M_Cini;}

	inline void setposition(size_type n) {M_position=n;}

	inline size_type getposition() {return M_position;}

	inline void setNcanali(size_type n) {M_ncanali=n;}

	inline void setCanali(double frazione) {
		if (frazione>0) 
			{M_canali.push_back(frazione);}
		else
			{M_canali.pop_back();}
	}

	inline size_type getNcanali() {return M_ncanali;}

	inline vector<double> getCanali() {return M_canali;}

	inline void setGlobalIndex(size_type n) {M_GlobalIndex=n;}
	
	inline size_type getGlobalIndex() const {return M_GlobalIndex;}

	void update_timeHis(std::vector<double> &);

	bool export_timeHis(string, vector<double>);

private:	
	
	string M_nome;
	bool M_isfluid;
	bool M_isgas;
	size_type M_ncanali;
	size_type M_GlobalIndex;
	vector<double> M_canali;
	size_type M_position;
	double M_Cini;
	std::vector<double> M_timeHis;
};


#endif /* GEOMTRIANGLE_HPP_ */

