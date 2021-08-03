/**
 * @file code.cpp
 * @author Daniel Felipe
 * @date 2020
 * @brief File containing the main function of the Hospital_AM model
 */


#include <iostream>
#include <fstream>
#include <time.h>

#include <dynamics.h>
#include <test.h>
#include <trace.h>
#include <reaction.h>
#include <other_functions.h>

typedef void(*reactions) (std::map<std::string, group> &Val, std::map<std::string, group> &Vba, Crandom &ran, Workers *altos, Workers *bajos, std::vector<lognormal_d> &dist, int agentI);

int main(void)
{
  //Creo los arreglos de cada tipo
  Workers altos[Na], bajos[Nb];

  //Creo los vectores en donde están las personas de cada estadío
  std::map<std::string, group> vecal{
				     { "SUS", {} }, { "SUST", {} }, { "SUSA", {} },
				     { "EXP", {} }, { "EXPT", {} }, { "EXPA", {} },
				     { "PRE", {} }, { "PRET", {} }, { "PREA", {} },
				     { "MSYM", {} }, { "MSYMT", {} }, { "MSYMA", {} },
				     { "SSYMA", {} },
				     { "RECI", {} }, { "RECT", {} }, { "RECA", {} }
  };

  std::map<std::string, group> vecba{
				     { "SUS", {} }, { "SUST", {} }, { "SUSA", {} },
				     { "EXP", {} }, { "EXPT", {} }, { "EXPA", {} },
				     { "PRE", {} }, { "PRET", {} }, { "PREA", {} },
				     { "MSYM", {} }, { "MSYMT", {} }, { "MSYMA", {} },
				     { "SSYMA", {} },
				     { "RECI", {} }, { "RECT", {} }, { "RECA", {} }
  }; //En el hospital

  //Creo el generador de semillas
  Crandom gseed(difftime(time(0),0));

  //Defino la cantidad de tiempo de la corrida
  int T = 720;
  double t;

  //Defino las variables del acordeón
  //double nu = myconstants["Nu"], delta = myconstants["Delta"];
  //unsigned int loops = T/(nu+delta);

  //Defino las variables para el testeo masivo
  unsigned int tests = (int)(N*MyCons.theta);
  double dt = MyCons.nu/((double)tests);

  //Defino la variable de tiempo propio de cada reacción(tj)
  double tj[14];

  //Defino el número de corridas
  unsigned int ensemble = 1e3;

  //Creo el arreglo de las funciones de reacción
  reactions react[14] = {reaction0, reaction1, reaction2, reaction3, reaction4, reaction5, reaction6, reaction7, reaction8, reaction9,
			 reaction10, reaction11, reaction12, reaction13};

  //Creo los arreglos para los argumentos de forma y orden de las distribuciones
  double proms[12] = {De, De, Dpl, Dpl, Dpg, Dpg, Dpl+Dil, Dpl+Dil, Dil, Dil, Dig, Dig};
  double d_shape, d_order, CoefVar = 3.9/5.2;
  d_shape = std::sqrt(std::log((CoefVar*CoefVar) + 1.0));
  std::vector<lognormal_d> dist;
  for(unsigned int i=0; i<12; i++){
    d_order = std::log(proms[i]) - 0.5*(d_shape*d_shape);
    lognormal_d my_dist(d_order, d_shape);
    dist.push_back(my_dist);
  }

  //Variables auxiliares
  std::vector<double> ti_in;
  unsigned int n1, n2;
  double aux, contador, promrec;
  unsigned int aux2, num;
  double ARi, ARh, ARc;

  std::ofstream fout;
  std::string name;

  name = "res/Data/datos_" + MyCons.name + ".csv";
  fout.open(name);

  contador = 0;
  promrec = 0;
  num = 0;

  //Genero las corridas
  while(contador < ensemble){
    std::cout << '\r';
    std::cout << "Vamos en la simulacion numero " << contador << " de " << ensemble << ".";

    //Le asigno a cada entrada el índice como valor
    vecal["SUS"].resize(Na);    vecba["SUS"].resize(Nb);
    for(unsigned int j=0; j<Na; j++){vecal["SUS"][j] = j;    altos[j].init();}
    for(unsigned int j=0; j<Nb; j++){vecba["SUS"][j] = j;    bajos[j].init();}

    name = "res/Data_Prevalencia/datos_" + MyCons.name + "_" + std::to_string(num) + ".csv";
    //name = "prueba.csv";
    //fout.open(name);
    //fout.close();

    //Inicio el tiempo
    t = 0.0;
    ARi = 0.0;
    ARh = 0.0;
    ARc = 0.0;

    //Imprimo los datos
    print_all(vecal, vecba, t, name);

    //Inicio los tiempos propios de cada reacción
    for(unsigned int j=0; j<14; j++){tj[j] = 0.0;}

    while(t < T){
      //Región de testeo masivo
      aux = 0.0;
      n1 = 0;
      while(aux < MyCons.nu){
	//Obtengo el tiempo e índice de la reacción
	ti_in = contagio(vecal, vecba, gseed, t, tj);

	//Si se tiene el tiempo máximo como tiempo mínimo, entonces termino la simulación
	if(ti_in[0] == 1e6){break;}

	//Voy contado de que tipo fue el contagio
	if((int)ti_in[1] == 0 || (int)ti_in[1] == 1){
	  if((int)ti_in[2] == -2){ARc++;}
	  else if((int)ti_in[2] == -1){ARh++;}
	  else{ARi++;}
	}

	//Actualizo los tiempos de los estados que pueden transitar
	update_times_all(vecal, vecba, altos, bajos, ti_in[0]);

	//Actualizo los tiempos de los testeados masivamente que dan negativo
	update_massive_all(vecal, vecba, altos, bajos, ti_in[0]);

	//Actualizo los tiempos de los leves aislado
	if(MyCons.AisLev){result_lev_ais(vecal, vecba, altos, bajos, ti_in[0], gseed);}

	//Actualizo los tiempos de los rastreados
	trace_massive_all(vecal, vecba, altos, bajos, ti_in[0]);

	//Actualizo los tiempos de los testeados y hago el rastreo de los nuevos aislados
	main_trace(vecal, vecba, altos, bajos, ti_in[0], gseed);

	//Genero los tests masivos
	/* Cuento los sus, exp, pre, lev y recI */
	aux2 = vecal["SUS"].size() + vecba["SUS"].size() + vecal["EXP"].size() + vecba["EXP"].size() + vecal["PRE"].size() + vecba["PRE"].size();
	aux2 += vecal["MSYM"].size() + vecba["MSYM"].size() + vecal["RECI"].size() + vecba["RECI"].size();
	n2 = (int)(aux/dt);
	for(unsigned int k=n1; k<n2 && k<tests && 0<aux2; k++){
	  massive_reaction(vecal, vecba, gseed, altos, bajos);
	  aux2 = vecal["SUS"].size() + vecba["SUS"].size() + vecal["EXP"].size() + vecba["EXP"].size() + vecal["PRE"].size() + vecba["PRE"].size();
	  aux2 += vecal["MSYM"].size() + vecba["MSYM"].size() + vecal["RECI"].size() + vecba["RECI"].size();
	}
	n1 = n2;

	//Genero la reacción según el índice que acabo de obtener
	react[(int)ti_in[1]](vecal, vecba, gseed, altos, bajos, dist, ti_in[2]);

	//Sumo el tiempo de la reacción
	t += ti_in[0];
	aux += ti_in[0];

	//Imprimo los datos
	print_all(vecal, vecba, t, name);

	//Borro el vector de tiempo e índice
	ti_in.clear();
      }

      //Muevo los testeados masivos a su respectivo lugar
      move_massive_all(vecal, vecba, altos, bajos);

      //Si el vector de tiempo e índice no se borró, es porque se rompió el ciclo
      if(ti_in.size() != 0){break;}

      //Región sin testeo masivo
      aux = 0.0;
      while(aux < MyCons.delta){
	//Obtengo el tiempo e índice de la reacción
	ti_in = contagio(vecal, vecba, gseed, t, tj);

	//Si se tiene el tiempo máximo como tiempo mínimo, entonces termino la simulación
	if(ti_in[0] == 1e6){break;}

	//Voy contado de que tipo fue el contagio
	if((int)ti_in[1] == 0 || (int)ti_in[1] == 1){
	  if((int)ti_in[2] == -2){ARc++;}
	  else if((int)ti_in[2] == -1){ARh++;}
	  else{ARi++;}
	}

	//Actualizo los tiempos de los estados que pueden transitar
	update_times_all(vecal, vecba, altos, bajos, ti_in[0]);

	//Actualizo los tiempos de los leves aislado
	if(MyCons.AisLev){result_lev_ais(vecal, vecba, altos, bajos, ti_in[0], gseed);}

	//Actualizo los tiempos de los rastreados
	trace_massive_all(vecal, vecba, altos, bajos, ti_in[0]);

	//Actualizo los tiempos de los testeados y hago el rastreo de los nuevos aislados
	main_trace(vecal, vecba, altos, bajos, ti_in[0], gseed);

	//Actualizo los tiempos de los testeados masivamente, y si ya cumplieron tiempo, los devuelvo
	tested_massive_all(vecal, vecba, altos, bajos, ti_in[0]);

	//Genero la reacción según el índice que acabo de obtener
	react[(int)ti_in[1]](vecal, vecba, gseed, altos, bajos, dist, ti_in[2]);

	//Sumo el tiempo de la reacción
	t += ti_in[0];
	aux += ti_in[0];

	//Imprimo los datos
        print_all(vecal, vecba, t, name);

	//Borro el vector de tiempo e índice
	ti_in.clear();
      }

      //Si el vector de tiempo e índice no se borró, es porque se rompió el ciclo
      if(ti_in.size() != 0){break;}
    }
    ti_in.clear();

    /* Cuento los recI, recT y recA */
    promrec += vecal["RECI"].size() + vecba["RECI"].size() + vecal["RECT"].size() + vecba["RECT"].size() + vecal["RECA"].size() + vecba["RECA"].size();
    fout << vecal["RECI"].size() + vecba["RECI"].size() + vecal["RECT"].size() + vecba["RECT"].size() + vecal["RECA"].size() + vecba["RECA"].size() << '\t';
    fout << ARi << '\t' << ARh << '\t' << ARc << std::endl;

    /* Imprimo la Red */
    name = "res/Data_Grafos/datos_" + MyCons.name + "_" + std::to_string(num) + ".csv";
    print_net(vecal, vecba, altos, bajos, name);

    //Borro los vectores
    for(std::map<std::string, group>::iterator it=vecal.begin(); it != vecal.end(); it++){
      vecal[it->first].clear();
      vecba[it->first].clear();
    }

    contador++;
    num++;
  }
  fout.close();

  promrec /= (contador*N);

  std::cout << "\nAR " << MyCons.name << ": " << promrec << std::endl;

  return 0;
}
