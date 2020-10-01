#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "Random64.h"
#include "bases.h"
#include "trabajadores.h"
#include "dynamics.h"

typedef void(*reactions) (grupo &Sa, grupo &Sb, grupo &Ea, grupo &Eb, grupo &Pa, grupo &Pb, grupo &PTa, grupo &PTb, grupo &PTAa, grupo &PTAb, grupo &La, grupo &Lb, grupo &LTa, grupo &LTb, grupo &LTAa, grupo &LTAb, grupo &LAa, grupo &LAb, grupo &IAa, grupo &IAb, grupo &Ra, grupo &Rb, Crandom &ran, trabajadores *altos, trabajadores *bajos);

int main(void)
{
  //Defino las proporciones de tipo bajo y tipo alto
  double pa = 0.18, pb = 0.82;
  const unsigned int Na = pa*N, Nb = pb*N;

  //Creo los arreglos de cada tipo
  trabajadores altos[Na], bajos[Nb];

  //Creo los vectores en donde están las personas de cada estadío
  grupo susal, susba;
  grupo expal, expba;
  grupo preal, preba;
  grupo preTal, preTba;
  grupo preTAal, preTAba;
  grupo leval, levba;
  grupo levTal, levTba;
  grupo levTAal, levTAba;
  grupo levAal, levAba;
  grupo infAal, infAba;
  grupo recal, recba;

  //Defino la prevalencia externa
  double prev = 0.22;

  //Creo el generador de semillas
  Crandom gseed(83554);

  //Defino la cantidad de tiempo de la corrida
  int T = 140;
  double t;

  //Defino las variables del acordeón
  double nu = 2, delta = 5;
  unsigned int loops = T/(nu+delta);

  //Defino las variables para el testeo masivo
  unsigned int tests = 500;
  double dt = nu/((double)tests);

  //Defino la variable de tiempo propio de cada reacción(tj)
  double tj[14];

  //Defino el número de corridas
  unsigned int ensemble = 1000;

  //Creo el arreglo de las funciones de reacción
  reactions react[14] = {reaction0, reaction1, reaction2, reaction3, reaction4, reaction5, reaction6, reaction7, reaction8, reaction9,
			 reaction10, reaction11, reaction12, reaction13};

  //Variables auxiliares
  std::vector<double> ti_in;
  unsigned int n1, n2, index;
  double aux;
  
  std::ofstream fout;
  std::string name;

  //Genero las corridas
  for(unsigned int i=0; i<ensemble; i++){
    std::cout << i << std::endl;

    //Le asigno a cada entrada el índice como valor
    susal.resize(Na);    susba.resize(Nb);
    for(unsigned int j=0; j<Na; j++){susal[j] = j;    altos[j].init();}
    for(unsigned int j=0; j<Nb; j++){susba[j] = j;    bajos[j].init();}

    //name = "Data/datos_" + std::to_string(i) + ".csv";
    name = "prueba.csv";
    fout.open(name);

    //Inicio el tiempo
    t = 0.0;
    fout << t << '\t' << susal.size() << '\t' << susba.size() << '\t';
    fout << expal.size() << '\t' << expba.size() << '\t';
    fout << preal.size() << '\t' << preba.size() << '\t';
    fout << preTal.size() << '\t' << preTba.size() << '\t';
    fout << preTAal.size() << '\t' << preTAba.size() << '\t';
    fout << leval.size() << '\t' << levba.size() << '\t';
    fout << levTal.size() << '\t' << levTba.size() << '\t';
    fout << levTAal.size() << '\t' << levTAba.size() << '\t';
    fout << levAal.size() << '\t' << levAba.size() << '\t';
    fout << infAal.size() << '\t' << infAba.size() << '\t';
    fout << recal.size() << '\t' << recba.size() << std::endl;
    
    //Inicio los tiempos propios de cada reacción
    for(unsigned int j=0; j<14; j++){tj[j] = 0.0;}   

    for(unsigned int j=0; j<loops; j++){
      
      //Región de testeo masivo
      aux = 0.0;
      n1 = 0;
      while(aux < nu){
	//Obtengo el tiempo e índice de la reacción
	ti_in = contagio(susal.size(), susba.size(), expal.size(), expba.size(), preal.size(), preba.size(), preTal.size(), preTba.size(), preTAal.size(), preTAba.size(), leval.size(), levba.size(), levTal.size(), levTba.size(), levTAal.size(), levTAba.size(), levAal.size(), levAba.size(), infAal.size(), infAba.size(), Na, Nb, prev, gseed, t, tj);
      
	//Si se tiene el tiempo máximo como tiempo mínimo, entonces termino la simulación
	if(ti_in[0] == 1e6){break;}	
	
	//Actualizo los tiempos de los testeados, y si ya les dieron resultado los aislo
	tested_isolated(preTal, preTAal, altos, ti_in[0], 3, 4);
	tested_isolated(preTba, preTAba, bajos, ti_in[0], 3, 4);
	tested_isolated(levTal, levTAal, altos, ti_in[0], 5, 6);
	tested_isolated(levTba, levTAba, bajos, ti_in[0], 5, 6);
      
	//Genero la reacción según el índice que acabo de obtener
	react[(int)ti_in[1]](susal, susba, expal, expba, preal, preba, preTal, preTba, preTAal, preTAba, leval, levba, levTal, levTba, levTAal, levTAba, levAal, levAba, infAal, infAba, recal, recba, gseed, altos, bajos);

	//Sumo el tiempo de la reacción
	t += ti_in[0];
	aux += ti_in[0];	

	//Genero lo tests continuos para los leves
	if((int)ti_in[1] == 4){continue_reaction(leval, levTal, altos, gseed);}
	else if((int)ti_in[1] == 5){continue_reaction(levba, levTba, bajos, gseed);}

	//Genero los tests masivos
	n2 = (int)(aux/dt);
	for(unsigned int k=n1; k<n2 && k<tests; k++){
	  if(gseed.r()*N < Na){massive_reaction(susal, expal, preal, preTal, leval, levTal, recal, gseed, altos);}
	  else{massive_reaction(susba, expba, preba, preTba, levba, levTba, recba, gseed, bajos);}
	}
	n1 = n2;
	
	fout << t << '\t' << susal.size() << '\t' << susba.size() << '\t';
	fout << expal.size() << '\t' << expba.size() << '\t';
	fout << preal.size() << '\t' << preba.size() << '\t';
	fout << preTal.size() << '\t' << preTba.size() << '\t';
	fout << preTAal.size() << '\t' << preTAba.size() << '\t';
	fout << leval.size() << '\t' << levba.size() << '\t';
	fout << levTal.size() << '\t' << levTba.size() << '\t';
	fout << levTAal.size() << '\t' << levTAba.size() << '\t';
	fout << levAal.size() << '\t' << levAba.size() << '\t';
	fout << infAal.size() << '\t' << infAba.size() << '\t';
	fout << recal.size() << '\t' << recba.size() << std::endl;
	
	//Borro el vector de tiempo e índice
	ti_in.clear();
      }

      //Si el vector de tiempo e índice no se borró, es porque se rompió el ciclo
      if(ti_in.size() != 0){break;}

      //Región sin testeo masivo
      aux = 0.0;
      while(aux < delta){
	//Obtengo el tiempo e índice de la reacción
	ti_in = contagio(susal.size(), susba.size(), expal.size(), expba.size(), preal.size(), preba.size(), preTal.size(), preTba.size(), preTAal.size(), preTAba.size(), leval.size(), levba.size(), levTal.size(), levTba.size(), levTAal.size(), levTAba.size(), levAal.size(), levAba.size(), infAal.size(), infAba.size(), Na, Nb, prev, gseed, t, tj);
      
	//Si se tiene el tiempo máximo como tiempo mínimo, entonces termino la simulación
	if(ti_in[0] == 1e6){break;}	

	//Actualizo los tiempos de los testeados, y si ya les dieron resultado los aislo
	tested_isolated(preTal, preTAal, altos, ti_in[0], 3, 4);
	tested_isolated(preTba, preTAba, bajos, ti_in[0], 3, 4);
	tested_isolated(levTal, levTAal, altos, ti_in[0], 5, 6);
	tested_isolated(levTba, levTAba, bajos, ti_in[0], 5, 6);
      
	//Genero la reacción según el índice que acabo de obtener
	react[(int)ti_in[1]](susal, susba, expal, expba, preal, preba, preTal, preTba, preTAal, preTAba, leval, levba, levTal, levTba, levTAal, levTAba, levAal, levAba, infAal, infAba, recal, recba, gseed, altos, bajos);	

	//Genero lo tests continuos para los leves
	if((int)ti_in[1] == 4){continue_reaction(leval, levTal, altos, gseed);}
	else if((int)ti_in[1] == 5){continue_reaction(levba, levTba, bajos, gseed);}
      
	//Sumo el tiempo de la reacción
	t += ti_in[0];
	aux += ti_in[0];
	
	fout << t << '\t' << susal.size() << '\t' << susba.size() << '\t';
	fout << expal.size() << '\t' << expba.size() << '\t';
	fout << preal.size() << '\t' << preba.size() << '\t';
	fout << preTal.size() << '\t' << preTba.size() << '\t';
	fout << preTAal.size() << '\t' << preTAba.size() << '\t';
	fout << leval.size() << '\t' << levba.size() << '\t';
	fout << levTal.size() << '\t' << levTba.size() << '\t';
	fout << levTAal.size() << '\t' << levTAba.size() << '\t';
	fout << levAal.size() << '\t' << levAba.size() << '\t';
	fout << infAal.size() << '\t' << infAba.size() << '\t';
	fout << recal.size() << '\t' << recba.size() << std::endl;
	
	//Borro el vector de tiempo e índice
	ti_in.clear();
      }
      
      //Si el vector de tiempo e índice no se borró, es porque se rompió el ciclo
      if(ti_in.size() != 0){break;}
    }
    fout.close();

    //Borro los vectores
    susal.clear();    susba.clear();
    expal.clear();    expba.clear();
    preal.clear();    preba.clear();
    preTal.clear();    preTba.clear();
    preTAal.clear();    preTAba.clear();
    leval.clear();    levba.clear();
    levTal.clear();    levTba.clear();
    levTAal.clear();    levTAba.clear();
    levAal.clear();    levAba.clear();
    infAal.clear();    infAba.clear();
    recal.clear();    recba.clear();
  }

  return 0;
}
