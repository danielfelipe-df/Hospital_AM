#include <test.h>
#include <trace.h>


void tested_reaction(grupo &Out, grupo &In, int index, Workers *family, Stages typeout, Stages typein, double value, bool dist){
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


void massive_reaction(std::map<std::string, grupo> &Val, std::map<std::string, grupo> &Vba, Crandom &ran, Workers *altos, Workers *bajos){
  unsigned int M = Val["SUS"].size() + Vba["SUS"].size() + Val["EXP"].size() + Vba["EXP"].size() + Val["PRE"].size() + Vba["PRE"].size() + Val["MSYM"].size() + Vba["MSYM"].size() + Val["RECI"].size() + Vba["RECI"].size();
  unsigned int num = (int)(ran.r()*M);

  if(num < Val["SUS"].size()){
    tested_reaction(Val["SUS"], Val["SUST"], (int)(ran.r()*Val["SUS"].size()), altos, SUS, SUST, ran.r(), true); //Susceptible alto
  }
  else if(num < Val["SUS"].size() + Vba["SUS"].size()){
    tested_reaction(Vba["SUS"], Vba["SUST"], (int)(ran.r()*Vba["SUS"].size()), bajos, SUS, SUST, ran.r(), true); //Suscpetible bajo
  }
  else if(num < Val["SUS"].size() + Vba["SUS"].size() + Val["EXP"].size()){
    tested_reaction(Val["EXP"], Val["EXPT"], (int)(ran.r()*Val["EXP"].size()), altos, EXP, EXPT, ran.r(), true); //Expuesto alto
  }
  else if(num < Val["SUS"].size() + Vba["SUS"].size() + Val["EXP"].size() + Vba["EXP"].size()){
    tested_reaction(Vba["EXP"], Vba["EXPT"], (int)(ran.r()*Vba["EXP"].size()), bajos, EXP, EXPT, ran.r(), true); //Expuesto bajo
  }
  else if(num < Val["SUS"].size() + Vba["SUS"].size() + Val["EXP"].size() + Vba["EXP"].size() + Val["PRE"].size()){
    tested_reaction(Val["PRE"], Val["PRET"], (int)(ran.r()*Val["PRE"].size()), altos, PRE, PRET, ran.r(), true); //Presintomático alto
  }
  else if(num < Val["SUS"].size() + Vba["SUS"].size() + Val["EXP"].size() + Vba["EXP"].size() + Val["PRE"].size() + Vba["PRE"].size()){
    tested_reaction(Vba["PRE"], Vba["PRET"], (int)(ran.r()*Vba["PRE"].size()), bajos, PRE, PRET, ran.r(), true); //Presintomático bajo
  }
  else if(num < Val["SUS"].size() + Vba["SUS"].size() + Val["EXP"].size() + Vba["EXP"].size() + Val["PRE"].size() + Vba["PRE"].size() + Val["MSYM"].size()){
    tested_reaction(Val["MSYM"], Val["MSYMT"], (int)(ran.r()*Val["MSYM"].size()), altos, MSYM, MSYMT, ran.r(), true); //Leve alto
  }
  else if(num < Val["SUS"].size() + Vba["SUS"].size() + Val["EXP"].size() + Vba["EXP"].size() + Val["PRE"].size() + Vba["PRE"].size() + Val["MSYM"].size() + Vba["MSYM"].size()){
    tested_reaction(Vba["MSYM"], Vba["MSYMT"], (int)(ran.r()*Vba["MSYM"].size()), bajos, MSYM, MSYMT, ran.r(), true); //Leve bajo
  }
  else if(num < Val["SUS"].size() + Vba["SUS"].size() + Val["EXP"].size() + Vba["EXP"].size() + Val["PRE"].size() + Vba["PRE"].size() + Val["MSYM"].size() + Vba["MSYM"].size() + Val["RECI"].size()){
    tested_reaction(Val["RECI"], Val["RECT"], (int)(ran.r()*Val["RECI"].size()), altos, RECI, RECT, ran.r(), true); //Recuperado no-detectado alto
  }
  else{
    tested_reaction(Vba["RECI"], Vba["RECT"], (int)(ran.r()*Vba["RECI"].size()), bajos, RECI, RECT, ran.r(), true); //Recuperado no-detectado bajo
  }
}


int tested_isolated_inf(grupo &T, grupo &TA, grupo &G, Workers *family, double time, Stages typeout, Stages typein1, Stages typein2, Crandom &ran){
  int index, contador = 0;
  for(size_t i=0; i<T.size(); i++){
    index = T[i];      family[index].time += time;
    if(family[index].time > family[index].tmax){
      if(ran.r() < MyCons.xi){
	family[index].change(typein1, typeout);	  TA.push_back(index);	contador++;
	if(typein1 == MSYMA){tested_lev_ais(index, family, 1e6, false);} //Si es Leve aislado a donde pasa pues se le pone el tlevmax en 1e6
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


void tested_massive(grupo &T, grupo &G, Workers *family, double time, Stages typeout, Stages typein){
  for(size_t i=0; i<T.size(); i++){
    family[T[i]].time += time;
    if(family[T[i]].time > family[T[i]].tmax){tested_reaction(T, G, i, family, typeout, typein, 0.0, false);      i--;}
  }
}


void move_massive(grupo &T, grupo &G, Workers *family, Stages typeout, Stages typein){
  for(size_t i=0; i<T.size(); i++){
    if(family[T[i]].time > family[T[i]].tmax){tested_reaction(T, G, i, family, typeout, typein, 0.0, false);      i--;}
  }
}


void result_lev_ais(std::map<std::string, grupo> &Val, std::map<std::string, grupo> &Vba, Workers *altos, Workers *bajos, double time, Crandom &ran){
  int index;
  unsigned int contador = 0;
  for(size_t i=0; i<Val["MSYMA"].size()-contador; i++){
    index = Val["MSYMA"][i];
    altos[index].tlev += time;
    if(altos[index].tlev > altos[index].tlevmax){
      if(ran.r() < MyCons.xi){
	tested_lev_ais(index, altos, 1e6, false);
	Val["MSYMA"].erase(Val["MSYMA"].begin() + i);
	Val["MSYMA"].push_back(index);
	aux_main(1, Val, Vba, altos, bajos, ran, MyCons.phi1, MyCons.mu, true, "MSYMA", false);
	contador++;
      }
      else{
	tested_reaction(Val["MSYMA"], Val["MSYM"], i, altos, MSYMA, MSYM, 0.0, false);
      }
      i--;
    }
  }

  contador = 0;
  for(size_t i=0; i<Vba["MSYMA"].size()-contador; i++){
    index = Vba["MSYMA"][i];
    bajos[index].tlev += time;
    if(bajos[index].tlev > bajos[index].tlevmax){
      if(ran.r() < MyCons.xi){
	tested_lev_ais(index, bajos, 1e6, false);
	Vba["MSYMA"].erase(Vba["MSYMA"].begin() + i);
	Vba["MSYMA"].push_back(index);
	aux_main(1, Val, Vba, altos, bajos, ran, MyCons.mu, MyCons.chi, false, "MSYMA", false);
	contador++;
      }
      else{
	tested_reaction(Vba["MSYMA"], Vba["MSYM"], i, bajos, MSYMA, MSYM, 0.0, false);
      }
      i--;
    }
  }
}


void move_massive_all(std::map<std::string, grupo> &Val, std::map<std::string, grupo> &Vba, Workers *altos, Workers *bajos){
  move_massive(Val["SUST"], Val["SUS"], altos, SUST, SUS);      move_massive(Vba["SUST"], Vba["SUS"], bajos, SUST, SUS); //Susceptible testeado a susceptible
  move_massive(Val["EXPT"], Val["EXP"], altos, EXPT, EXP);      move_massive(Vba["EXPT"], Vba["EXP"], bajos, EXPT, EXP); //Expuesto testeado a expuesto
  move_massive(Val["RECT"], Val["RECI"], altos, RECT, RECI);      move_massive(Vba["RECT"], Vba["RECI"], bajos, RECT, RECI); //Recuperado testeado a recuperado
}


void tested_massive_all(std::map<std::string, grupo> &Val, std::map<std::string, grupo> &Vba, Workers *altos, Workers *bajos, double time){
  tested_massive(Val["SUST"], Val["SUS"], altos, time, SUST, SUS);  tested_massive(Vba["SUST"], Vba["SUS"], bajos, time, SUST, SUS); //Susceptible testeado a susceptible
  tested_massive(Val["EXPT"], Val["EXP"], altos, time, EXPT, EXP);  tested_massive(Vba["EXPT"], Vba["EXP"], bajos, time, EXPT, EXP); //Expuesto testeado a expuesto
  tested_massive(Val["RECT"], Val["RECI"], altos, time, RECT, RECI);  tested_massive(Vba["RECT"], Vba["RECI"], bajos, time, RECT, RECI); //Recuperado testeado a recuperado
}
