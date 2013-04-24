
#include "chimica.hpp"

using namespace std;
using namespace gmm;


//----------------costruttore vuoto-------------

React::React(): M_reagente(), M_prodotti(0) {};

//distruttore
React::~React() {};

//---------------calcola e restituisce la velocità di reazione------------

vector<double> React::getRvel(double Temp, vector<double> & conc, map<string,Species*> mappa){
	vector<double> vel;
	size_type who(mappa[M_reagente]->getGlobalIndex());

	for (size_type i=0; i<M_ncanali; ++i){  //ciclo sui vari canali

		double provv;
		provv=M_Af[i]*exp(-M_Eact[i]/(Rgas*Temp))*conc[who+i]; //formula di Arrhenius
		vel.push_back(provv);
	}
	
	return vel;
}



