#include "integrate.hpp"

using namespace std;
using namespace gmm;


//----------------costruttore vuoto-------------

TimeInt::TimeInt() {};

//--------------------esportazione dei tempi su file------------------------

void TimeInt::esportaTempi(string path)
{
	string nomefile(path);
	nomefile.append("time");
	ofstream myfile;
	myfile.open(nomefile.c_str());
	for (gmm::size_type i=0; i<M_tempi.size();++i)
	{
		myfile<< M_tempi[i] <<endl;
	}
}

//--------------------integrazione con eulero esplicito------------------------

void TimeInt::onestepEE(Kem *chimica, std::vector<double> & c_old,std::vector<double> & c_new )
{
	vector<React>  allreact(chimica->getAllReact());
	double Temp(273+0.5*M_tempi[M_tempi.size()-2]); //calcolo la temperatura con una legge, in realtà la dovrei caricare da file
	// dato che è esplito la valuto al tempo n
	
	size_type n_canaliR(chimica->contaCanaliReazioni()); 
	size_type n_canaliS(chimica->contaCanaliSpecie());

	vector<double> r(n_canaliR,0.0);    //vettore velocità di reazione
	vector<double> rhs(n_canaliS,0.0);  //forzante del sistema ottenuto come S per r

	for (size_type i=0; i<chimica->nReazioni();++i)  //per ogni reazione
	{
		for(size_type ii=0; ii<allreact[i].getNcanali();++ii)  //per ogni canale di questa reazione
		{ 
			r[allreact[i].getGlobalIndex()+ii]=allreact[i].getRvel(Temp,c_old, chimica->getIndexMap())[ii];  //calcolo la rate
		}
	}

	mult(chimica->getS(),r,rhs);  //moltiplico la matrice stechiometrica per le reaction rate

	gmm::scale(rhs,conv_t*M_DeltaT);   //moltiplico per dt
	
	for (gmm::size_type i=0;i<c_new.size();++i)
	{ 
		c_new[i]=c_old[i]+rhs[i];    //aggiorno le concentrazioni

	}
}

//--------------------integrazione con eulero implicito------------------------

void TimeInt::onestepEI(Kem *chimica, std::vector<double> & c_old,std::vector<double> & c_new )
{
	vector<React>  allreact(chimica->getAllReact());
	double Temp(273+0.5*M_tempi[M_tempi.size()-1]);
	size_type n_canaliR(chimica->contaCanaliReazioni());
	size_type n_canaliS(chimica->contaCanaliSpecie());
	for (gmm::size_type i=0;i<c_new.size();++i)
	{ 
		c_new[i]=c_old[i];  //guess iniziale per la soluzione (è un problema implicito potenzialmente non lineare)
	}
	
	size_type cont(0);  //contatore n. iterazioni
	double err(1);      //indicatore dell'errore, inizializzato a 1 (>toll)
	double normc;
	
	while (err>1.0e-15)
	{
		vector<double> r(n_canaliR,0.0);  //vettore velocità di reazione
		vector<double> rhs(n_canaliS,0.0);  //forzante del sistema ottenuto come S per r

		for (size_type i=0; i<chimica->nReazioni();++i) //per ogni reazione
		{
			for(size_type ii=0; ii<allreact[i].getNcanali();++ii)//per ogni canale di questa reazione
			{ 
				r[allreact[i].getGlobalIndex()+ii]=allreact[i].getRvel(Temp,c_new, chimica->getIndexMap())[ii]; //calcolo la rate
			}
		}
		mult(chimica->getS(),r,rhs);
		gmm::scale(rhs,conv_t*M_DeltaT);
		err=0;
		normc=0;
		for (gmm::size_type i=0;i<c_new.size();++i)
		{ 
			err=err+(c_new[i]-(c_old[i]+rhs[i]))*(c_new[i]-(c_old[i]+rhs[i]));   //calcolo la norma del residuo, sommando i quadrati di ogni comp.
			c_new[i]=c_old[i]+rhs[i];    //aggiorno la soluzione
			normc=normc+c_new[i]*c_new[i];   //  calcolo allo stesso modo la norma della soluzione
		}
		cont=cont+1;   //aggiorno il contatore 
		err=pow(err,0.5)/pow(normc,0.5);  // calcolo residuo normalizzato
	}//fine while

}

//--------------------integrazione con Runge Kutta 2-----------------------

void TimeInt::onestepRK2(Kem *chimica, std::vector<double> & c_old,std::vector<double> & c_new )
{
	vector<React>  allreact(chimica->getAllReact());
	double Temp(273+0.5*M_tempi[M_tempi.size()-2]);   //temperatura al tempo n e n+1
	double Temp2(273+0.5*M_tempi[M_tempi.size()-1]);

	size_type n_canaliR(chimica->contaCanaliReazioni());  
	size_type n_canaliS(chimica->contaCanaliSpecie());
	
	vector<double> r(n_canaliR,0.0);
	vector<double> rhs(n_canaliS,0.0),rhs2(n_canaliS,0.0);
	std::vector<double> cprovv(c_old);  //vettore per la soluzione intermedia

	for (size_type i=0; i<chimica->nReazioni();++i)
	{
		for(size_type ii=0; ii<allreact[i].getNcanali();++ii)
		{ 
			r[allreact[i].getGlobalIndex()+ii]=allreact[i].getRvel(Temp,c_old, chimica->getIndexMap())[ii];
		}
	}
	mult(chimica->getS(),r,rhs);  //.... come sopra

	gmm::scale(rhs,conv_t*M_DeltaT);
	for (gmm::size_type i=0;i<c_new.size();++i)
	{ 
		cprovv[i]=c_old[i]+rhs[i];  //questo è il primo stadio del RK
	}

	//adesso calcolo le velocità di reazione corrispondenti a cprovv, e con la temperatura al tempo n+1
	for (size_type i=0; i<chimica->nReazioni();++i){
		for(size_type ii=0; ii<allreact[i].getNcanali();++ii){ 
			r[allreact[i].getGlobalIndex()+ii]=allreact[i].getRvel(Temp2,cprovv, chimica->getIndexMap())[ii];
		}
	}
	mult(chimica->getS(),r,rhs2);

	gmm::scale(rhs2,conv_t*M_DeltaT);
	for (gmm::size_type i=0;i<c_new.size();++i)
	{ 
		c_new[i]=c_old[i]+0.5*rhs[i]+0.5*rhs2[i];   //con i due rhs mediati ottengo la soluzione al tempo n+1
	}
}

//--------------------integrazione con lo schema di Patankar (1 ordine)-----------------------

void TimeInt::onestepPATEE(Kem *chimica, std::vector<double> & c_old,std::vector<double> & c_new )
{
	vector<React>  allreact(chimica->getAllReact());
	double Temp(273+0.5*M_tempi[M_tempi.size()-2]);
	size_type n_canaliR(chimica->contaCanaliReazioni());
	size_type n_canaliS(chimica->contaCanaliSpecie());

	vector<double> r(n_canaliR,0.0);
	vector<double> rhs(n_canaliS,0.0),rhs2(n_canaliS,0.0);

	// devo costruire la matrice P che mi dice "quanto" la produzione della specie i è dovuta alla specie j...
	// formula nella tesi di Anna e nel paper Scotti Formaggia Positivity...
	gmm::dense_matrix<double> RR(n_canaliR,n_canaliR); 
	gmm::dense_matrix<double> P(n_canaliS,n_canaliS);
	gmm::dense_matrix<double> provv(n_canaliR,n_canaliS);
	gmm::dense_matrix<double> M(n_canaliS,n_canaliS);

	for (size_type i=0; i<chimica->nReazioni();++i)
	{
		for(size_type ii=0; ii<allreact[i].getNcanali();++ii)
		{ 
			r[allreact[i].getGlobalIndex()+ii]=allreact[i].getRvel(Temp,c_old, chimica->getIndexMap())[ii];   //rate di reazione
			//le metto in una matrice, sulla diagonale. inefficiente, lo so!
			RR(allreact[i].getGlobalIndex()+ii,allreact[i].getGlobalIndex()+ii)=r[allreact[i].getGlobalIndex()+ii];
		}
	}

	mult(RR, transposed(chimica->getSmt()), provv);  //le moltiplico per i coefficienti stechiometrici DEI REAGENTI (quelli =-1)
	mult(chimica->getSp(),provv, P);   // moltiplico per i coefficienti stechiometrici dei prodotti (quelli >0)
	//ottengo con questa formula la matrice P

	// nel metodo di patankar devo risolvere un sistema lineare Mc_new=c_old
	// costruisco M
	for (gmm::size_type i=0;i<c_new.size();++i)
	{ 
		M(i,i)=1.;
		double s(0);
		for (gmm::size_type j=0;j<c_new.size();++j)
		{ 	
			if (c_old[j]>0)
			{
				M(i,j)=M(i,j)-M_DeltaT*conv_t*P(i,j)/c_old[j];
			}
			s=s+P(j,i);
		}
		if (c_old[i]>0)
		{
			M(i,i)=M(i,i)+M_DeltaT*conv_t*s/c_old[i];
		}
	}

double rcond;
lu_solve(M, c_new, c_old);  //risolvo e ottengo la nuova c

}

//--------------------integrazione con lo schema di Patankar (2 ordine) INCOMPLETO-----------------------

void TimeInt::onestepPATRK2(Kem *chimica, std::vector<double> & c_old,std::vector<double> & c_new ){
vector<React>  allreact(chimica->getAllReact());
double Temp(273+0.5*M_tempi[M_tempi.size()-2]);
double Temp2(273+0.5*M_tempi[M_tempi.size()-1]);
size_type n_canaliR(chimica->contaCanaliReazioni());
size_type n_canaliS(chimica->contaCanaliSpecie());

vector<double> r(n_canaliR,0.0);
vector<double> rhs(n_canaliS,0.0),rhs2(n_canaliS,0.0);
std::vector<double> cprovv(c_old);
gmm::dense_matrix<double> RR(n_canaliR,n_canaliR);
gmm::dense_matrix<double> P(n_canaliS,n_canaliS),P2(n_canaliS,n_canaliS);
gmm::dense_matrix<double> provv(n_canaliR,n_canaliS);
gmm::dense_matrix<double> M(n_canaliS,n_canaliS);


for (size_type i=0; i<chimica->nReazioni();++i){
	for(size_type ii=0; ii<allreact[i].getNcanali();++ii){ 
		r[allreact[i].getGlobalIndex()+ii]=allreact[i].getRvel(Temp,c_old, chimica->getIndexMap())[ii];
		RR(allreact[i].getGlobalIndex()+ii,allreact[i].getGlobalIndex()+ii)=r[allreact[i].getGlobalIndex()+ii];
	}
}

mult(RR, transposed(chimica->getSmt()), provv);
mult(chimica->getSp(),provv, P);

//primo step pat
for (gmm::size_type i=0;i<c_new.size();++i)
{ 
	M(i,i)=1.;
	double s(0);

	for (gmm::size_type j=0;j<c_new.size();++j)
	{ 	
		if (c_old[j]>0){
		M(i,j)=M(i,j)-M_DeltaT*conv_t*P(i,j)/c_old[j];}
		s=s+P(j,i);
	}
	if (c_old[i]>0){
	M(i,i)=M(i,i)+M_DeltaT*conv_t*s/c_old[i];
}
}
double rcond;
lu_solve(M, c_new, c_old);
 //secondo step pat

/*
for (size_type i=0; i<chimica->nReazioni();++i){
	for(size_type ii=0; ii<allreact[i].getNcanali();++ii){ 
		r[allreact[i].getGlobalIndex()+ii]=allreact[i].getRvel(Temp2,cprovv, chimica->getIndexMap())[ii];
		RR(allreact[i].getGlobalIndex()+ii,allreact[i].getGlobalIndex()+ii)=r[allreact[i].getGlobalIndex()+ii];
	}
}
clear(provv);
mult(RR, transposed(chimica->getSmt()), provv);
mult(chimica->getSp(),provv, P2);

//primo step pat
for (gmm::size_type i=0;i<c_new.size();++i)
{ 
	M(i,i)=1.;
	double s(0);

	for (gmm::size_type j=0;j<c_new.size();++j)
	{ 	
		if (c_old[j]>0){
		M(i,j)=M(i,j)-M_DeltaT*conv_t*P(i,j)/c_old[j];}
		s=s+P(j,i);
	}
	if (c_old[i]>0){
	M(i,i)=M(i,i)+M_DeltaT*conv_t*s/c_old[i];
}
}
double rcond;
lu_solve(M, c_new, c_old);

*/
}

