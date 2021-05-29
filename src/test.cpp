#include <test.h>
#include <trace.h>


void tested_reaction(grupo &Out, grupo &In, int index, Workers *family, int typeout, int typein, double value, bool dist){
  int agent = Out[index];
  Out.erase(Out.begin() + index);
  In.push_back(agent);
  family[agent].change(typein, typeout);
  family[agent].time = 0.0;
  if(dist){family[agent].tmax = quantile(dist_Tt, value);}
  else{family[agent].tmax = value;}
}


void tested_lev_ais(int agent, Workers *family, double value, bool dist){
  family[agent].tlev = 0.0;
  if(dist){family[agent].tlevmax = quantile(dist_Tt, value);}
  else{family[agent].tlevmax = value;}
}


void massive_reaction(std::vector<grupo> &Val, std::vector<grupo> &Vba, Crandom &ran, Workers *altos, Workers *bajos){
  unsigned int M = Val[0].size() + Vba[0].size() + Val[3].size() + Vba[3].size() + Val[6].size() + Vba[6].size() + Val[9].size() + Vba[9].size() + Val[13].size() + Vba[13].size();
  unsigned int num = (int)(ran.r()*M);

  if(num < Val[0].size()){
    tested_reaction(Val[0], Val[1], (int)(ran.r()*Val[0].size()), altos, 0, 1, ran.r(), true); //Susceptible alto
  }
  else if(num < Val[0].size() + Vba[0].size()){
    tested_reaction(Vba[0], Vba[1], (int)(ran.r()*Vba[0].size()), bajos, 0, 1, ran.r(), true); //Suscpetible bajo
  }
  else if(num < Val[0].size() + Vba[0].size() + Val[3].size()){
    tested_reaction(Val[3], Val[4], (int)(ran.r()*Val[3].size()), altos, 3, 4, ran.r(), true); //Expuesto alto
  }
  else if(num < Val[0].size() + Vba[0].size() + Val[3].size() + Vba[3].size()){
    tested_reaction(Vba[3], Vba[4], (int)(ran.r()*Vba[3].size()), bajos, 3, 4, ran.r(), true); //Expuesto bajo
  }
  else if(num < Val[0].size() + Vba[0].size() + Val[3].size() + Vba[3].size() + Val[6].size()){
    tested_reaction(Val[6], Val[7], (int)(ran.r()*Val[6].size()), altos, 6, 7, ran.r(), true); //Presintomático alto
  }
  else if(num < Val[0].size() + Vba[0].size() + Val[3].size() + Vba[3].size() + Val[6].size() + Vba[6].size()){
    tested_reaction(Vba[6], Vba[7], (int)(ran.r()*Vba[6].size()), bajos, 6, 7, ran.r(), true); //Presintomático bajo
  }
  else if(num < Val[0].size() + Vba[0].size() + Val[3].size() + Vba[3].size() + Val[6].size() + Vba[6].size() + Val[9].size()){
    tested_reaction(Val[9], Val[10], (int)(ran.r()*Val[9].size()), altos, 9, 10, ran.r(), true); //Leve alto
  }
  else if(num < Val[0].size() + Vba[0].size() + Val[3].size() + Vba[3].size() + Val[6].size() + Vba[6].size() + Val[9].size() + Vba[9].size()){
    tested_reaction(Vba[9], Vba[10], (int)(ran.r()*Vba[9].size()), bajos, 9, 10, ran.r(), true); //Leve bajo
  }
  else if(num < Val[0].size() + Vba[0].size() + Val[3].size() + Vba[3].size() + Val[6].size() + Vba[6].size() + Val[9].size() + Vba[9].size() + Val[13].size()){
    tested_reaction(Val[13], Val[14], (int)(ran.r()*Val[13].size()), altos, 13, 14, ran.r(), true); //Recuperado no-detectado alto
  }
  else{
    tested_reaction(Vba[13], Vba[14], (int)(ran.r()*Vba[13].size()), bajos, 13, 14, ran.r(), true); //Recuperado no-detectado bajo
  }
}


int tested_isolated_inf(grupo &T, grupo &TA, grupo &G, Workers *family, double time, int typeout, int typein1, int typein2, Crandom &ran){
  int index, contador = 0;
  for(unsigned int i=0; i<T.size(); i++){
    index = T[i];      family[index].time += time;
    if(family[index].time > family[index].tmax){
      if(ran.r() < MyCons.xi){
	family[index].change(typein1, typeout);	  TA.push_back(index);	contador++;
	if(typein1 == 11){tested_lev_ais(index, family, 1e6, false);} //Si es Leve aislado a donde pasa pues se le pone el tlevmax en 1e6
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


void tested_massive(grupo &T, grupo &G, Workers *family, double time, int typeout, int typein){
  for(unsigned int i=0; i<T.size(); i++){
    family[T[i]].time += time;
    if(family[T[i]].time > family[T[i]].tmax){tested_reaction(T, G, i, family, typeout, typein, 0.0, false);      i--;}
  }
}


void move_massive(grupo &T, grupo &G, Workers *family, unsigned int typeout, unsigned int typein){
  for(unsigned int i=0; i<T.size(); i++){
    if(family[T[i]].time > family[T[i]].tmax){tested_reaction(T, G, i, family, typeout, typein, 0.0, false);      i--;}
  }
}


void result_lev_ais(std::vector<grupo> &Val, std::vector<grupo> &Vba, Workers *altos, Workers *bajos, double time, Crandom &ran){
  int index;
  unsigned int contador = 0;
  for(unsigned int i=0; i<Val[11].size()-contador; i++){
    index = Val[11][i];
    altos[index].tlev += time;
    if(altos[index].tlev > altos[index].tlevmax){
      if(ran.r() < MyCons.xi){
	tested_lev_ais(index, altos, 1e6, false);
	Val[11].erase(Val[11].begin() + i);	Val[11].push_back(index);
	aux_main(1, Val, Vba, altos, bajos, ran, MyCons.phi1, MyCons.mu, true, 11, false);
	contador++;
      }
      else{
	tested_reaction(Val[11], Val[9], i, altos, 11, 9, 0.0, false);
      }
      i--;
    }
  }

  contador = 0;
  for(unsigned int i=0; i<Vba[11].size()-contador; i++){
    index = Vba[11][i];
    bajos[index].tlev += time;
    if(bajos[index].tlev > bajos[index].tlevmax){
      if(ran.r() < MyCons.xi){
	tested_lev_ais(index, bajos, 1e6, false);
	Vba[11].erase(Vba[11].begin() + i);	Vba[11].push_back(index);
	aux_main(1, Val, Vba, altos, bajos, ran, MyCons.mu, MyCons.chi, false, 11, false);
	contador++;
      }
      else{
	tested_reaction(Vba[11], Vba[9], i, bajos, 11, 9, 0.0, false);
      }
      i--;
    }
  }
}


void move_massive_all(std::vector<grupo> &Val, std::vector<grupo> &Vba, Workers *altos, Workers *bajos){
  move_massive(Val[1], Val[0], altos, 1, 0);      move_massive(Vba[1], Vba[0], bajos, 1, 0); //Susceptible testeado a susceptible
  move_massive(Val[4], Val[3], altos, 4, 3);      move_massive(Vba[4], Vba[3], bajos, 4, 3); //Expuesto testeado a expuesto
  move_massive(Val[14], Val[13], altos, 14, 13);      move_massive(Vba[14], Vba[13], bajos, 14, 13); //Recuperado testeado a recuperado
}


void tested_massive_all(std::vector<grupo> &Val, std::vector<grupo> &Vba, Workers *altos, Workers *bajos, double time){
  tested_massive(Val[1], Val[0], altos, time, 1, 0);  tested_massive(Vba[1], Vba[0], bajos, time, 1, 0); //Susceptible testeado a susceptible
  tested_massive(Val[4], Val[3], altos, time, 4, 3);  tested_massive(Vba[4], Vba[3], bajos, time, 4, 3); //Expuesto testeado a expuesto
  tested_massive(Val[14], Val[13], altos, time, 14, 13);  tested_massive(Vba[14], Vba[13], bajos, time, 14, 13); //Recuperado testeado a recuperado
}
