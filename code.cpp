#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "Random64.h"
#include "bases.h"
#include "trabajadores.h"
#include "dynamics.h"

typedef void(*reactions) (grupo &Sa, grupo &Sb, grupo &Ea, grupo &Eb, grupo &Pa, grupo &Pb, grupo &PTAa, grupo &PTAb, grupo &La, grupo &Lb, grupo &LTAa, grupo &LTAb, grupo &LAa, grupo &LAb, grupo &IAa, grupo &IAb, grupo &Ra, grupo &Rb, Crandom &ran);

int main(void)
{
  //Defino las proporciones de tipo bajo y tipo alto
  double pa = 0.18, pb = 0.82;
  const int Na = pa*N, Nb = pb*N;

  //Creo los arreglos de cada tipo
  trabajadores altos[Na], bajos[Nb];

  //Creo los vectores en donde están las personas de cada estadío
  grupo susal, susba;
  grupo expal, expba;
  grupo preal, preba;
  grupo preTAal, preTAba;
  grupo leval, levba;
  grupo levTAal, levTAba;
  grupo levAal, levAba;
  grupo infAal, infAba;
  grupo recal, recba;

  //Defino la prevalencia externa
  double prev = 0.03;

  //Creo el generador de semillas
  Crandom gseed(1717171);

  //Defino la cantidad de tiempo de la corrida
  int T = 140;
  double t;

  //Defino las variables del acordeón
  double r;
  double nu = 2, delta = 5;
  int loops = T/(nu+delta);

  //Defino el número de corridas
  int ensemble = 1;

  //Creo el arreglo de las funciones de reacción
  reactions react[26] = {reaction0, reaction1, reaction2, reaction3, reaction4, reaction5, reaction6, reaction7, reaction8, reaction9,
			 reaction10, reaction11, reaction12, reaction13, reaction14, reaction15, reaction16, reaction17, reaction18,
			 reaction19, reaction20, reaction21, reaction22, reaction23, reaction24, reaction25};

  //Variables auxiliares
  std::vector<double> ti_in;
  int n;
  double aux;
  
  std::ofstream fout;
  std::string name;

  //Genero las corridas
  for(unsigned int i=0; i<ensemble; i++){

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
    fout << preTAal.size() << '\t' << preTAba.size() << '\t';
    fout << leval.size() << '\t' << levba.size() << '\t';
    fout << levTAal.size() << '\t' << levTAba.size() << '\t';
    fout << levAal.size() << '\t' << levAba.size() << '\t';
    fout << infAal.size() << '\t' << infAba.size() << '\t';
    fout << recal.size() << '\t' << recba.size() << std::endl;

    for(unsigned j=0; j<loops; j++){
      aux = 0.0;
      r = 1.0;
      while(aux < nu){
	//Obtengo el tiempo e índice de la reacción
	ti_in = contagio(susal.size(), susba.size(), expal.size(), expba.size(), preal.size(), preba.size(), preTAal.size(), preTAba.size(), leval.size(), levba.size(), levTAal.size(), levTAba.size(), levAal.size(), levAba.size(), infAal.size(), infAba.size(), Na, Nb, prev, gseed, r);

	//Si se tiene el tiempo máximo como tiempo mínimo, entonces termino la simulación
	if(ti_in[0] == 1e6){break;}

	//Genero la reacción según el índice que acabo de obtener
	react[(int)ti_in[1]](susal, susba, expal, expba, preal, preba, preTAal, preTAba, leval, levba, levTAal, levTAba, levAal, levAba, infAal, infAba, recal, recba, gseed);

	//Sumo el tiempo de la reacción
	t += ti_in[0];
	aux += ti_in[0];

	fout << t << '\t' << susal.size() << '\t' << susba.size() << '\t';
	fout << expal.size() << '\t' << expba.size() << '\t';
	fout << preal.size() << '\t' << preba.size() << '\t';
	fout << preTAal.size() << '\t' << preTAba.size() << '\t';
	fout << leval.size() << '\t' << levba.size() << '\t';
	fout << levTAal.size() << '\t' << levTAba.size() << '\t';
	fout << levAal.size() << '\t' << levAba.size() << '\t';
	fout << infAal.size() << '\t' << infAba.size() << '\t';
	fout << recal.size() << '\t' << recba.size() << std::endl;
	
	//Borro el vector de tiempo e índice
	ti_in.clear();
      }

      aux = 0.0;
      r = 0.0;
      while(aux < delta){
	//Obtengo el tiempo e índice de la reacción
	ti_in = contagio(susal.size(), susba.size(), expal.size(), expba.size(), preal.size(), preba.size(), preTAal.size(), preTAba.size(), leval.size(), levba.size(), levTAal.size(), levTAba.size(), levAal.size(), levAba.size(), infAal.size(), infAba.size(), Na, Nb, prev, gseed, r);

	//Si se tiene el tiempo máximo como tiempo mínimo, entonces termino la simulación
	if(ti_in[0] == 1e6){break;}

	//Genero la reacción según el índice que acabo de obtener
	react[(int)ti_in[1]](susal, susba, expal, expba, preal, preba, preTAal, preTAba, leval, levba, levTAal, levTAba, levAal, levAba, infAal, infAba, recal, recba, gseed);

	//Sumo el tiempo de la reacción
	t += ti_in[0];
	aux += ti_in[0];

	fout << t << '\t' << susal.size() << '\t' << susba.size() << '\t';
	fout << expal.size() << '\t' << expba.size() << '\t';
	fout << preal.size() << '\t' << preba.size() << '\t';
	fout << preTAal.size() << '\t' << preTAba.size() << '\t';
	fout << leval.size() << '\t' << levba.size() << '\t';
	fout << levTAal.size() << '\t' << levTAba.size() << '\t';
	fout << levAal.size() << '\t' << levAba.size() << '\t';
	fout << infAal.size() << '\t' << infAba.size() << '\t';
	fout << recal.size() << '\t' << recba.size() << std::endl;
	
	//Borro el vector de tiempo e índice
	ti_in.clear();
      }
    }
    fout.close();

    //Borro los vectores
    susal.clear();    susba.clear();
    expal.clear();    expba.clear();
    preal.clear();    preba.clear();
    preTAal.clear();    preTAba.clear();
    leval.clear();    levba.clear();
    levTAal.clear();    levTAba.clear();
    levAal.clear();    levAba.clear();
    infAal.clear();    infAba.clear();
    recal.clear();    recba.clear();
  }


  return 0;
}
