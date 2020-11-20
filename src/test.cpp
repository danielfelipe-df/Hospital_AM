#include <test.h>
#include <trace.h>


void tested_reaction(grupo &Out, grupo &In, int index, trabajadores *family, int typeout, int typein, double delta){
  int agent = Out[index];
  Out.erase(Out.begin() + index);
  In.push_back(agent);
  family[agent].change(typein, typeout);
  family[agent].time = 0.0;
  family[agent].tmax = delta;
}


void massive_reaction(std::vector<grupo> &Val, std::vector<grupo> &Vba, Crandom &ran, trabajadores *altos, trabajadores *bajos){
  unsigned int M = Val[0].size() + Vba[0].size() + Val[2].size() + Vba[2].size() + Val[5].size() + Vba[5].size() + Val[8].size() + Vba[8].size() + Val[12].size() + Vba[12].size();
  unsigned int num = (int)(ran.r()*M);

  if(num < Val[0].size()){
    tested_reaction(Val[0], Val[1], (int)(ran.r()*Val[0].size()), altos, 0, 1, Tt); //Susceptible alto
  }
  else if(num < Val[0].size() + Vba[0].size()){
    tested_reaction(Vba[0], Vba[1], (int)(ran.r()*Vba[0].size()), bajos, 0, 1, Tt); //Suscpetible bajo
  }
  else if(num < Val[0].size() + Vba[0].size() + Val[2].size()){
    tested_reaction(Val[2], Val[3], (int)(ran.r()*Val[2].size()), altos, 2, 3, Tt); //Expuesto alto
  }
  else if(num < Val[0].size() + Vba[0].size() + Val[2].size() + Vba[2].size()){
    tested_reaction(Vba[2], Vba[3], (int)(ran.r()*Vba[2].size()), bajos, 2, 3, Tt); //Expuesto bajo
  }
  else if(num < Val[0].size() + Vba[0].size() + Val[2].size() + Vba[2].size() + Val[5].size()){
    tested_reaction(Val[5], Val[6], (int)(ran.r()*Val[5].size()), altos, 5, 6, Tt); //Presintomático alto
  }
  else if(num < Val[0].size() + Vba[0].size() + Val[2].size() + Vba[2].size() + Val[5].size() + Vba[5].size()){
    tested_reaction(Vba[5], Vba[6], (int)(ran.r()*Vba[5].size()), bajos, 5, 6, Tt); //Presintomático bajo
  }
  else if(num < Val[0].size() + Vba[0].size() + Val[2].size() + Vba[2].size() + Val[5].size() + Vba[5].size() + Val[8].size()){
    tested_reaction(Val[8], Val[9], (int)(ran.r()*Val[8].size()), altos, 8, 9, Tt); //Leve alto
  }
  else if(num < Val[0].size() + Vba[0].size() + Val[2].size() + Vba[2].size() + Val[5].size() + Vba[5].size() + Val[8].size() + Vba[8].size()){
    tested_reaction(Vba[8], Vba[9], (int)(ran.r()*Vba[8].size()), bajos, 8, 9, Tt); //Leve bajo
  }
  else if(num < Val[0].size() + Vba[0].size() + Val[2].size() + Vba[2].size() + Val[5].size() + Vba[5].size() + Val[8].size() + Vba[8].size() + Val[12].size()){
    tested_reaction(Val[12], Val[13], (int)(ran.r()*Val[12].size()), altos, 12, 13, Tt); //Recuperado no-detectado alto
  }
  else{
    tested_reaction(Vba[12], Vba[13], (int)(ran.r()*Vba[12].size()), bajos, 12, 13, Tt); //Recuperado no-detectado bajo
  }
}


int tested_isolated_inf(grupo &T, grupo &TA, grupo &G, trabajadores *family, double time, int typeout, int typein1, int typein2, Crandom &ran){
  int index, contador = 0;
  for(unsigned int i=0; i<T.size(); i++){
    index = T[i];      family[index].time += time;
    if(family[index].time > family[index].tmax){
      if(ran.r() < xi){
	family[index].change(typein1, typeout);	  TA.push_back(index);	contador++;
	if(typein1 == 10){tested_lev_ais(index, family, 1e6);} //Si es Leve aislado a donde pasa pues se le pone el tlevmax en 1e6
      }
      else{
	family[index].change(typein2, typeout);	  G.push_back(index);
      }
      T.erase( T.begin() + i);
      family[index].time = 0.0;
      family[index].tmax = 0.0;
      i--;
    }
  }
  return contador;
}


void tested_lev_ais(int agent, trabajadores *family, double Tmax){
  family[agent].tlev = 0.0;
  family[agent].tlevmax = Tmax;
}


void tested_massive(grupo &T, grupo &G, trabajadores *family, double time, int typeout, int typein){
  for(unsigned int i=0; i<T.size(); i++){
    family[T[i]].time += time;
    if(family[T[i]].time > family[T[i]].tmax){tested_reaction(T, G, i, family, typeout, typein, 0.0);      i--;}
  }
}


void move_massive(grupo &T, grupo &G, trabajadores *family, unsigned int typeout, unsigned int typein){
  for(unsigned int i=0; i<T.size(); i++){
    if(family[T[i]].time > family[T[i]].tmax){tested_reaction(T, G, i, family, typeout, typein, 0.0);      i--;}
  }
}


void result_lev_ais(std::vector<grupo> &Val, std::vector<grupo> &Vba, trabajadores *altos, trabajadores *bajos, double time, Crandom &ran){
  int index;
  unsigned int contador = 0;
  for(unsigned int i=0; i<Val[10].size()-contador; i++){
    index = Val[10][i];
    altos[index].tlev += time;
    if(altos[index].tlev > altos[index].tlevmax){
      if(ran.r() < xi){
	tested_lev_ais(index, altos, 1e6);
	Val[10].erase(Val[10].begin() + i);	Val[10].push_back(index);
	aux_main(1, Val, Vba, altos, bajos, ran, phi1, mu, true, 10);
	contador++;
      }
      else{
	tested_reaction(Val[10], Val[8], i, altos, 10, 8, 0.0);
      }
      i--;
    }
  }

  contador = 0;
  for(unsigned int i=0; i<Vba[10].size()-contador; i++){
    index = Vba[10][i];
    bajos[index].tlev += time;
    if(bajos[index].tlev > bajos[index].tlevmax){
      if(ran.r() < xi){
	tested_lev_ais(index, bajos, 1e6);
	Vba[10].erase(Vba[10].begin() + i);	Vba[10].push_back(index);
	aux_main(1, Val, Vba, altos, bajos, ran, mu, chi, false, 10);
	contador++;
      }
      else{
	tested_reaction(Vba[10], Vba[8], i, bajos, 10, 8, 0.0);
      }
      i--;
    }
  }
}


void move_massive_all(std::vector<grupo> &Val, std::vector<grupo> &Vba, trabajadores *altos, trabajadores *bajos){
  move_massive(Val[1], Val[0], altos, 1, 0);      move_massive(Vba[1], Vba[0], bajos, 1, 0); //Susceptible testeado a susceptible
  move_massive(Val[3], Val[2], altos, 3, 2);      move_massive(Vba[3], Vba[2], bajos, 3, 2); //Expuesto testeado a expuesto
  move_massive(Val[13], Val[12], altos, 13, 12);      move_massive(Vba[13], Vba[12], bajos, 13, 12); //Recuperado testeado a recuperado
}


void tested_massive_all(std::vector<grupo> &Val, std::vector<grupo> &Vba, trabajadores *altos, trabajadores *bajos, double time){
  tested_massive(Val[1], Val[0], altos, time, 1, 0);  tested_massive(Vba[1], Vba[0], bajos, time, 1, 0); //Susceptible testeado a susceptible
  tested_massive(Val[3], Val[2], altos, time, 3, 2);  tested_massive(Vba[3], Vba[2], bajos, time, 3, 2); //Expuesto testeado a expuesto
  tested_massive(Val[13], Val[12], altos, time, 13, 12);  tested_massive(Vba[13], Vba[12], bajos, time, 13, 12); //Recuperado testeado a recuperado
}
