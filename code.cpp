#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "Random64.h"
#include "bases.h"
#include "trabajadores.h"
#include "dynamics.h"

typedef void(*reactions) (grupo &Sa, grupo &Sb, grupo &STa, grupo &Stb, grupo &Ea, grupo &Eb, grupo &ETa, grupo &Etb, grupo &Pa, grupo &Pb, grupo &PTa, grupo &PTb, grupo &PTAa, grupo &PTAb, grupo &La, grupo &Lb, grupo &LTa, grupo &LTb, grupo &LTAa, grupo &LTAb, grupo &IAa, grupo &IAb, grupo &RTa, grupo &RTb, grupo &RIa, grupo &RIb, grupo &RAa, grupo &RAb, Crandom &ran, trabajadores *altos, trabajadores *bajos, std::vector<weib_d> &dist);

int main(void)
{
  //Defino las proporciones de tipo bajo y tipo alto
  double pa = 0.18, pb = 0.82;
  const unsigned int Na = pa*N, Nb = pb*N;

  //Creo los arreglos de cada tipo
  trabajadores altos[Na], bajos[Nb];

  //Creo los vectores en donde están las personas de cada estadío
  grupo susal, susba;  grupo susTal, susTba;
  grupo expal, expba;  grupo expTal, expTba;
  grupo preal, preba;  grupo preTal, preTba;  grupo preTAal, preTAba;
  grupo leval, levba;  grupo levTal, levTba;  grupo levTAal, levTAba;
  grupo infAal, infAba;
  grupo recTal, recTba;  grupo recIal, recIba;  grupo recAal, recAba;

  //Defino la prevalencia externa
  double prev = 0.22;

  //Creo el generador de semillas
  Crandom gseed(83554);

  //Defino la cantidad de tiempo de la corrida
  int T = 420;
  double t;

  //Defino las variables del acordeón
  double nu = 2, delta = 0;
  unsigned int loops = T/(nu+delta);

  //Defino las variables para el testeo masivo
  unsigned int tests = (int)(N*theta);
  double dt = nu/((double)tests);
  unsigned int value_tests;

  //Defino la variable de tiempo propio de cada reacción(tj)
  double tj[14];

  //Defino el número de corridas
  unsigned int ensemble = 1;

  //Creo el arreglo de las funciones de reacción
  reactions react[14] = {reaction0, reaction1, reaction2, reaction3, reaction4, reaction5, reaction6, reaction7, reaction8, reaction9,
			 reaction10, reaction11, reaction12, reaction13};

  //Creo los arreglos para los argumentos de forma y orden de las distribuciones
  double d_order[12] = {De, De, Dpl, Dpl, Dpg, Dpg, Dpl+Dil, Dpl+Dil, Dil, Dil, Dig, Dig};
  double d_shape = 2.0;
  std::vector<weib_d> dist;
  for(unsigned int i=0; i<12; i++){
    weib_d my_dist(d_shape, d_order[i]);
    dist.push_back(my_dist);
  }
  
  //Variables auxiliares
  std::vector<double> ti_in;
  unsigned int n1, n2;
  //unsigned int index;
  double aux;
  unsigned int aux2;
  
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
    fout << susTal.size() << '\t' << susTba.size() << '\t';
    fout << expal.size() << '\t' << expba.size() << '\t';
    fout << expTal.size() << '\t' << expTba.size() << '\t';
    fout << preal.size() << '\t' << preba.size() << '\t';
    fout << preTal.size() << '\t' << preTba.size() << '\t';
    fout << preTAal.size() << '\t' << preTAba.size() << '\t';
    fout << leval.size() << '\t' << levba.size() << '\t';
    fout << levTal.size() << '\t' << levTba.size() << '\t';
    fout << levTAal.size() << '\t' << levTAba.size() << '\t';
    fout << infAal.size() << '\t' << infAba.size() << '\t';
    fout << recTal.size() << '\t' << recTba.size() << '\t';
    fout << recIal.size() << '\t' << recIba.size() << '\t';
    fout << recAal.size() << '\t' << recAba.size() << std::endl;
    
    //Inicio los tiempos propios de cada reacción
    for(unsigned int j=0; j<14; j++){tj[j] = 0.0;}

    for(unsigned int j=0; j<loops; j++){
      
      //Región de testeo masivo
      aux = 0.0;
      n1 = 0;
      while(aux < nu){
	//Obtengo el tiempo e índice de la reacción
	ti_in = contagio(susal.size(), susba.size(), susTal.size(), susTba.size(), expal.size(), expba.size(), expTal.size(), expTba.size(), preal.size(), preba.size(), preTal.size(), preTba.size(), preTAal.size(), preTAba.size(), leval.size(), levba.size(), levTal.size(), levTba.size(), levTAal.size(), levTAba.size(), infAal.size(), infAba.size(), Na, Nb, prev, gseed, t, tj);
      
	//Si se tiene el tiempo máximo como tiempo mínimo, entonces termino la simulación
	if(ti_in[0] == 1e6){break;}

	//Actualizo los tiempos de los estados que pueden transitar
	update_times(expal, altos, ti_in[0]);	update_times(expba, bajos, ti_in[0]);
	update_times(expTal, altos, ti_in[0]);	update_times(expTba, bajos, ti_in[0]);
	update_times(preal, altos, ti_in[0]);	update_times(preba, bajos, ti_in[0]);
	update_times(preTal, altos, ti_in[0]);	update_times(preTba, bajos, ti_in[0]);
	update_times(preTAal, altos, ti_in[0]);	update_times(preTAba, bajos, ti_in[0]);
	update_times(leval, altos, ti_in[0]);	update_times(levba, bajos, ti_in[0]);
	update_times(levTal, altos, ti_in[0]);	update_times(levTba, bajos, ti_in[0]);
	update_times(levTAal, altos, ti_in[0]);	update_times(levTAba, bajos, ti_in[0]);
	update_times(infAal, altos, ti_in[0]);	update_times(infAba, bajos, ti_in[0]);
	
	//Actualizo los tiempos de los testeados, y si ya les dieron resultado los aislo
	tested_isolated_inf(preTal, preTAal, preal, altos, ti_in[0], 3, 4, 4, gseed);
	tested_isolated_inf(preTba, preTAba, preba, bajos, ti_in[0], 3, 4, 4, gseed);
	tested_isolated_inf(levTal, levTAal, leval, altos, ti_in[0], 5, 6, 6, gseed);
	tested_isolated_inf(levTba, levTAba, levba, bajos, ti_in[0], 5, 6, 6, gseed);

	//Actualizo el tiempo de los testeados que van a dar negativo, y si ya cumplieron el tiempo los devuelvo al estado sin testear	
	tested_isolated(susTal, susal, altos, ti_in[0], 1, 1, gseed);
	tested_isolated(susTba, susba, bajos, ti_in[0], 1, 1, gseed);
	tested_isolated(expTal, expal, altos, ti_in[0], 1, 1, gseed);
	tested_isolated(expTba, expba, bajos, ti_in[0], 1, 1, gseed);
	tested_isolated(recTal, recIal, altos, ti_in[0], 1, 1, gseed);
	tested_isolated(recTba, recIba, bajos, ti_in[0], 1, 1, gseed);
      
	//Genero la reacción según el índice que acabo de obtener
	react[(int)ti_in[1]](susal, susba, susTal, susTba, expal, expba, expTal, expTba, preal, preba, preTal, preTba, preTAal, preTAba, leval, levba, levTal, levTba, levTAal, levTAba, infAal, infAba, recTal, recTba, recIal, recIba, recAal, recAba, gseed, altos, bajos, dist);

	//Sumo el tiempo de la reacción
	t += ti_in[0];
	aux += ti_in[0];

	//Genero los tests masivos
	aux2 = susal.size() + susba.size() + expal.size() + expba.size() + preal.size() + preba.size();
	aux2 += leval.size() + levba.size() + recIal.size() + recIba.size();
	n2 = (int)(aux/dt);
	for(unsigned int k=n1; k<n2 && k<tests && 0<aux2; k++){
	  massive_reaction(susal, susba, susTal, susTba, expal, expba, expTal, expTba, preal, preba, preTal, preTba, leval, levba, levTal, levTba, recIal, recIba, recTal, recTba, gseed, altos, bajos);
	  aux2 = susal.size() + susba.size() + expal.size() + expba.size() + preal.size() + preba.size();
	  aux2 += leval.size() + levba.size() + recIal.size() + recIba.size();
	}
	n1 = n2;
	
	fout << t << '\t' << susal.size() << '\t' << susba.size() << '\t';
	fout << susTal.size() << '\t' << susTba.size() << '\t';
	fout << expal.size() << '\t' << expba.size() << '\t';
	fout << expTal.size() << '\t' << expTba.size() << '\t';
	fout << preal.size() << '\t' << preba.size() << '\t';
	fout << preTal.size() << '\t' << preTba.size() << '\t';
	fout << preTAal.size() << '\t' << preTAba.size() << '\t';
	fout << leval.size() << '\t' << levba.size() << '\t';
	fout << levTal.size() << '\t' << levTba.size() << '\t';
	fout << levTAal.size() << '\t' << levTAba.size() << '\t';
	fout << infAal.size() << '\t' << infAba.size() << '\t';
	fout << recTal.size() << '\t' << recTba.size() << '\t';
	fout << recIal.size() << '\t' << recIba.size() << '\t';
	fout << recAal.size() << '\t' << recAba.size() << std::endl;
	
	//Borro el vector de tiempo e índice
	ti_in.clear();
      }      

      //Si el vector de tiempo e índice no se borró, es porque se rompió el ciclo
      if(ti_in.size() != 0){break;}

      //Región sin testeo masivo
      aux = 0.0;
      while(aux < delta){
	//Obtengo el tiempo e índice de la reacción
	ti_in = contagio(susal.size(), susba.size(), susTal.size(), susTba.size(), expal.size(), expba.size(), expTal.size(), expTba.size(), preal.size(), preba.size(), preTal.size(), preTba.size(), preTAal.size(), preTAba.size(), leval.size(), levba.size(), levTal.size(), levTba.size(), levTAal.size(), levTAba.size(), infAal.size(), infAba.size(), Na, Nb, prev, gseed, t, tj);
      
	//Si se tiene el tiempo máximo como tiempo mínimo, entonces termino la simulación
	if(ti_in[0] == 1e6){break;}

	//Actualizo los tiempos de los estados que pueden transitar
	update_times(expal, altos, ti_in[0]);	update_times(expba, bajos, ti_in[0]);
	update_times(expTal, altos, ti_in[0]);	update_times(expTba, bajos, ti_in[0]);
	update_times(preal, altos, ti_in[0]);	update_times(preba, bajos, ti_in[0]);
	update_times(preTal, altos, ti_in[0]);	update_times(preTba, bajos, ti_in[0]);
	update_times(preTAal, altos, ti_in[0]);	update_times(preTAba, bajos, ti_in[0]);
	update_times(leval, altos, ti_in[0]);	update_times(levba, bajos, ti_in[0]);
	update_times(levTal, altos, ti_in[0]);	update_times(levTba, bajos, ti_in[0]);
	update_times(levTAal, altos, ti_in[0]);	update_times(levTAba, bajos, ti_in[0]);
	update_times(infAal, altos, ti_in[0]);	update_times(infAba, bajos, ti_in[0]);

	//Actualizo los tiempos de los testeados, y si ya les dieron resultado los aislo
	tested_isolated_inf(preTal, preTAal, preal, altos, ti_in[0], 3, 4, 4, gseed);
	tested_isolated_inf(preTba, preTAba, preba, bajos, ti_in[0], 3, 4, 4, gseed);
	tested_isolated_inf(levTal, levTAal, leval, altos, ti_in[0], 5, 6, 6, gseed);
	tested_isolated_inf(levTba, levTAba, levba, bajos, ti_in[0], 5, 6, 6, gseed);

	//Actualizo el tiempo de los testeados que van a dar negativo, y si ya cumplieron el tiempo los devuelvo al estado sin testear	
	tested_isolated(susTal, susal, altos, ti_in[0], 1, 1, gseed);
	tested_isolated(susTba, susba, bajos, ti_in[0], 1, 1, gseed);
	tested_isolated(expTal, expal, altos, ti_in[0], 1, 1, gseed);
	tested_isolated(expTba, expba, bajos, ti_in[0], 1, 1, gseed);
	tested_isolated(recTal, recIal, altos, ti_in[0], 1, 1, gseed);
	tested_isolated(recTba, recIba, bajos, ti_in[0], 1, 1, gseed);
      
	//Genero la reacción según el índice que acabo de obtener
	react[(int)ti_in[1]](susal, susba, susTal, susTba, expal, expba, expTal, expTba, preal, preba, preTal, preTba, preTAal, preTAba, leval, levba, levTal, levTba, levTAal, levTAba, infAal, infAba, recTal, recTba, recIal, recIba, recAal, recAba, gseed, altos, bajos, dist);
      
	//Sumo el tiempo de la reacción
	t += ti_in[0];
	aux += ti_in[0];
	
	fout << t << '\t' << susal.size() << '\t' << susba.size() << '\t';
	fout << susTal.size() << '\t' << susTba.size() << '\t';
	fout << expal.size() << '\t' << expba.size() << '\t';
	fout << expTal.size() << '\t' << expTba.size() << '\t';
	fout << preal.size() << '\t' << preba.size() << '\t';
	fout << preTal.size() << '\t' << preTba.size() << '\t';
	fout << preTAal.size() << '\t' << preTAba.size() << '\t';
	fout << leval.size() << '\t' << levba.size() << '\t';
	fout << levTal.size() << '\t' << levTba.size() << '\t';
	fout << levTAal.size() << '\t' << levTAba.size() << '\t';
	fout << infAal.size() << '\t' << infAba.size() << '\t';
	fout << recTal.size() << '\t' << recTba.size() << '\t';
	fout << recIal.size() << '\t' << recIba.size() << '\t';
	fout << recAal.size() << '\t' << recAba.size() << std::endl;
	
	//Borro el vector de tiempo e índice
	ti_in.clear();
      }
      
      //Si el vector de tiempo e índice no se borró, es porque se rompió el ciclo
      if(ti_in.size() != 0){break;}
    }
    fout.close();    

    //Borro los vectores
    susal.clear();    susba.clear();    susTal.clear();    susTba.clear();
    expal.clear();    expba.clear();    expTal.clear();    expTba.clear();
    preal.clear();    preba.clear();    preTal.clear();    preTba.clear();    preTAal.clear();    preTAba.clear();
    leval.clear();    levba.clear();    levTal.clear();    levTba.clear();    levTAal.clear();    levTAba.clear();
    infAal.clear();    infAba.clear();
    recTal.clear();    recTba.clear();    recIal.clear();    recIba.clear();    recAal.clear();    recAba.clear();
  }

  return 0;
}
