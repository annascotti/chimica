// i/o example

/**
 * Interior Penalty for Darcy's equation
 *  
*/

// #include <utility>
#include "getfem/getfem_regular_meshes.h"
#include "getfem/getfem_interpolation.h"
#include "getfem/getfem_config.h"
#include "getfem/getfem_assembling.h" // import assembly methods (and comp. of 
#include "getfem/getfem_import.h"
#include "getfem/getfem_export.h"     // export functions (save the solution in a file)
#include "gmm/gmm.h"
#include "gmm/gmm_inoutput.h"
#include "gmm/gmm_MUMPS_interface.h"
#include "gmm/gmm_superlu_interface.h"

// Level Set and Xfem stuff:
#include "getfem/getfem_mesh_im_level_set.h"
#include "getfem/getfem_mesh_fem_level_set.h"
#include "getfem/getfem_mesh_fem_product.h"
#include "getfem/getfem_mesh_fem_global_function.h"


/* try to enable the SIGFPE if something evaluates to a Not-a-number
 * of infinity during computations
 */      // System Matrix
#ifdef GETFEM_HAVE_FEENABLEEXCEPT
#  include <fenv.h>
#endif

#include <iostream>
#include <string>
#include <cstring>
#include <sstream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <vector>

//#include "src/ricerche.hpp"

#include "src/chimica.hpp"
#include "src/integrate.hpp"


using namespace gmm;
using namespace std;

int main ()
{

	ifstream myfile;
	myfile.open ("KEMFILE1.kem");  //apro il file di dati che contiene le info sulla chimica

	Kem bla("KEMFILE1.kem");    //costruisco l'oggetto della classe kem 

	size_type n_canaliR(bla.contaCanaliReazioni());  //conto il numero totale di canali di reazione--> il numero di velocità di reazione
	size_type n_canaliS(bla.contaCanaliSpecie());  //conto il numero di specie o frazioni di specie che voglio tracciare->la dimensione del sistema ODE

	dense_matrix<double> S(n_canaliS, n_canaliR); //matrice stechiometrica
	
	bla.assemblyS();	//assemblo la matrice stechiometrica
	bla.setCiniCanali();	//setto le condizioni iniziali per le specie (o frazioni di ogni specie: in pratica la c_0)

	vector<double> conc(bla.getCiniCanali());  //inizializzo il vettore soluzione

	TimeInt timeloop;  //costruisco l'oggetto timeloop che si occupa dell'integrazione in tempo 
	timeloop.setTini(0);
	timeloop.setTfin(350);   //tutti i tempi sono in Ma
	timeloop.setDeltaT(0.01);

	cout << scaled(conc,1.)<<endl;  //scrivo a schermo le c_0
	
	double tt;   //variabile tempo
	timeloop.addstep(tt);  //inizialmente zero

	while (tt<timeloop.getTfin()) //ciclo
	{   
		//cout << "tempo  "<<tt<<endl;
		tt=tt+timeloop.getDeltaT();   //aggiorno il tempo
		timeloop.addstep(tt);
		std::vector<double> conc_old(conc);	//la concentrazione vecchia
		timeloop.onestepPATEE(&bla,conc_old, conc );    //faccio un passo in tempo con il metodo che voglio
		bla.update_timehis(conc);   //aggiorno la soluzione
	}

	cout << scaled(conc,1.)<<endl;   //scrivo a schermo le concentrazioni finali

	//timeloop.esportaTempi("out/");   //esporto se voglio il vettore di tempi

	for (size_type i=0; i<bla.getAllSpecies().size();++i)  //esporto le verie specie su file
	{
		bla.getAllSpecies()[i].export_timeHis("out/", timeloop.getTempi());
	}


  return 0;
}
