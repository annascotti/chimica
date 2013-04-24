
#include "specie.hpp"


using namespace std;
using namespace gmm;


//----------------costruttore vuoto-------------

Species::Species(): M_nome(), M_isfluid(0), M_isgas(0) {};

//----------------costruttore--------------------

Species::Species(string nome): M_nome(nome), M_ncanali(1) 
{
	for (size_type i=0; i<sizeof(allnames)/sizeof(string);++i){  //capisco qual è di quelle predefinite 
		if (strcmp(allnames[i].c_str(), nome.c_str())==0){  //di conseguenza setto se è gas, se è solida o fluida
			M_isfluid=allisfluid[i];
			M_isgas=allisgas[i];
		}
	}
	this->setCanali(1);
}

//distruttore
Species::~Species() {};

//------------dato che una specie può essere divisa in vari canali, per poi esportarne la concentrazione li sommo tutti

void Species::update_timeHis(std::vector<double> & conc){

	double totale(0);
	for (size_type i=0; i<M_ncanali; ++i){
		totale=totale+conc[M_GlobalIndex + i];
	}
	M_timeHis.push_back(totale);
}

//------------------scrivo sul file tempo e concentrazione----------------

bool Species::export_timeHis(string path, vector<double> tempi)
{
	string nomefile(path);
	nomefile.append(M_nome.c_str());
	ofstream myfile;
	myfile.open(nomefile.c_str());
	for (gmm::size_type i=0; i<M_timeHis.size();++i){
		myfile<<tempi[i]<<"  " <<M_timeHis[i] <<endl;
	}

	myfile.close();
	//scrittura su file

	return true;
}
