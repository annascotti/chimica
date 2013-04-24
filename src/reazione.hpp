#ifndef reazione_HPP_
#define reazione_HPP_

#include "specie.hpp"
#include <string>
#include "gmm/gmm.h"
#include <iostream>
#include <string>
#include <iostream>
#include <fstream>

using namespace gmm;
using namespace std;

const double Rgas(1.9872);  //in calorie 
class React
{
public:
	//! @name Constructor & Destructor
	//@{
		
	//! Empty constructor
	React();
	
	//! Destructor
	virtual ~React();

	inline void setReactant(string nome) {M_reagente=nome;}

	inline string getReactant() const {return M_reagente;}

	inline string getProduct(size_type n) const {return M_prodotti[n];}
	
	inline void setNprod(size_type num) {M_nprod=num;}

	inline size_type getNprod() {return M_nprod;}
	
	inline void addProduct(string nome) {M_prodotti.push_back(nome);}

	inline void setNcanali(size_type n) {M_ncanali=n;}

	inline size_type getNcanali() const {return M_ncanali;}

	inline void setEact(double energ) {M_Eact.push_back(energ);}

	inline void setAf(double energ) {M_Af.push_back(energ);}

	inline void setFr(double energ) {M_Fr.push_back(energ);}

	inline void setStcoefficients(double coeff) {if (coeff<0) M_stReagente=coeff; else M_stprodotti.push_back(coeff);}

	inline double getStRcoefficient() {return M_stReagente;}

	inline vector<double> getStProdcoefficients() {return M_stprodotti;}

	inline void setGlobalIndex(size_type n) {M_GlobalIndex=n;}
	
	inline size_type getGlobalIndex() const {return M_GlobalIndex;}

	vector<double> getRvel(double, vector<double> &,  map<string,Species*>);
private:	
	string M_reagente;    //nome del reagente
	vector<string> M_prodotti;  //nomi dei prodotti
	size_type M_nprod;
	size_type M_ncanali;
	size_type M_GlobalIndex;
	vector<double> M_Eact;  //energie di attivazione  
	vector<double> M_Af;    //fattori preesponenziali
	vector<double> M_Fr;    //frazione di reazione che va con quella velocità
	double M_stReagente;    //coefficiente stechio del reagente
	vector<double> M_stprodotti; //coefficienti stechio dei prodotti
};


#endif 

