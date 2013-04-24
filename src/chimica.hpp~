#include "reazione.hpp"
#include "specie.hpp"
#include <string>
#include "gmm/gmm.h"
#include <iostream>
#include <string>
#include <iostream>
#include <fstream>
#include <map>

#ifndef chimica_HPP_
#define chimica_HPP_


using namespace gmm;
using namespace std;

//classe che contiene le info relaive alla chimica, quante reazioni, quante specie...

class Kem
{
public:
	//! @name Constructor & Destructor
	//@{
		
	//! Empty constructor
	Kem();
	
	explicit Kem(string);

	//! Destructor
	virtual ~Kem();
	
	inline gmm::size_type nReazioni() const {return M_nreazioni;}

	void setSpecies();

	void setCini();

	void setCiniCanali();

	void setGlobalIndeces();

	

	size_type contaCanaliSpecie();

	size_type contaCanaliReazioni();

	void assemblyS();

//------------getters-------------------------
	gmm::dense_matrix<double> getS() const {return M_S;}
	gmm::dense_matrix<double> getSp() const {return M_Spiu;}
	gmm::dense_matrix<double> getSmt() const {return M_Smenotilde;}
	vector<double> getCiniCanali() const {return M_CiniCanali;}
	inline vector <Species> getAllSpecies() const {return M_species;}
	inline vector <React> getAllReact() const {return M_reazioni;}
	map<string,Species*> getIndexMap() const {return indexMap;}

	void update_timehis(std::vector<double> &);

private:	
	string M_nomefile;   //file .kem
	size_type M_nspecie;  //numero di specie coinvolte
	size_type M_nreazioni; //numero reazioni
	size_type M_ncanaliR;  //canali di reazione
	size_type M_ncanaliS;  //frazioni di specie
	gmm::dense_matrix<double> M_S;   //matrice stechiometrica
	gmm::dense_matrix<double> M_Spiu;  //parte positiva
	gmm::dense_matrix<double> M_Smenotilde;  //parte negativa normalizzata
	vector<React> M_reazioni;    //vettore contenente le reazioni 
	vector <Species> M_species;  //vettore delle specie
	vector <Species*> M_species_pointer;  //vettore di puntatori alle specie
        map<string,Species*> indexMap;   //mappa stra il nome della specie e il puntatore 
	vector<double> M_CiniCanali; //condizioni iniziali per le specie

//un paio di metodi privati
	size_type contaReazioni(std::string);  
	React buildReact(string, string, size_type);
};


#endif /* GEOMTRIANGLE_HPP_ */

