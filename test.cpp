#include "test.h"


void tested_reaction(grupo &Out, grupo &In, int index, trabajadores *family, int typeout, int typein, double delta){
  int agent = Out[index];
  Out.erase(Out.begin() + index);
  In.push_back(agent);
  family[agent].change(typein, typeout);
  family[agent].time = 0.0;
  family[agent].tmax = delta;
}


void massive_reaction(std::vector<grupo> &Val, std::vector<grupo> &Vba, Crandom &ran, trabajadores *altos, trabajadores *bajos){
  unsigned int M = Val[0].size() + Vba[0].size() + Val[2].size() + Vba[2].size() + Val[5].size() + Vba[5].size() + Val[9].size() + Vba[9].size() + Val[14].size() + Vba[14].size();
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
  else if(num < Val[0].size() + Vba[0].size() + Val[2].size() + Vba[2].size() + Val[5].size() + Vba[5].size() + Val[9].size()){
    tested_reaction(Val[9], Val[10], (int)(ran.r()*Val[9].size()), altos, 9, 10, Tt); //Leve alto
  }
  else if(num < Val[0].size() + Vba[0].size() + Val[2].size() + Vba[2].size() + Val[5].size() + Vba[5].size() + Val[9].size() + Vba[9].size()){
    tested_reaction(Vba[9], Vba[10], (int)(ran.r()*Vba[9].size()), bajos, 9, 10, Tt); //Leve bajo
  }
  else if(num < Val[0].size() + Vba[0].size() + Val[2].size() + Vba[2].size() + Val[5].size() + Vba[5].size() + Val[9].size() + Vba[9].size() + Val[14].size()){
    tested_reaction(Val[14], Val[15], (int)(ran.r()*Val[14].size()), altos, 14, 15, Tt); //Recuperado no-detectado alto
  }
  else{
    tested_reaction(Vba[14], Vba[15], (int)(ran.r()*Vba[14].size()), bajos, 14, 15, Tt); //Recuperado no-detectado bajo
  }
}


int tested_isolated_inf(grupo &T, grupo &TA, grupo &G, trabajadores *family, double time, int typeout, int typein1, int typein2, Crandom &ran){
  int index, contador = 0;
  for(unsigned int i=0; i<T.size(); i++){
    index = T[i];      family[index].time += time;
    if(family[index].time > family[index].tmax){
      if(ran.r() < xi){family[index].change(typein1, typeout);	  TA.push_back(index);	contador++;}
      else{family[index].change(typein2, typeout);	  G.push_back(index);}
      T.erase( T.begin() + i);
      family[index].time = 0.0;
      family[index].tmax = 0.0;
      i--;
    }
  }
  return contador;
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


void move_massive_all(std::vector<grupo> &Val, std::vector<grupo> &Vba, trabajadores *altos, trabajadores *bajos){
  move_massive(Val[1], Val[0], altos, 1, 0);      move_massive(Vba[1], Vba[0], bajos, 1, 0); //Susceptible testeado a susceptible
  move_massive(Val[3], Val[2], altos, 3, 2);      move_massive(Vba[3], Vba[2], bajos, 3, 2); //Expuesto testeado a expuesto
  move_massive(Val[8], Val[5], altos, 8, 5);      move_massive(Vba[8], Vba[5], bajos, 8, 5); //Presintomático masivo a presintomático
  move_massive(Val[12], Val[9], altos, 12, 9);      move_massive(Vba[12], Vba[9], bajos, 12, 9); //Leve masivo a leve
  move_massive(Val[15], Val[14], altos, 15, 14);      move_massive(Vba[15], Vba[14], bajos, 15, 14); //Recuperado testeado a recuperado
}


void tested_massive_all(std::vector<grupo> &Val, std::vector<grupo> &Vba, trabajadores *altos, trabajadores *bajos, double time){
  tested_massive(Val[1], Val[0], altos, time, 1, 0);  tested_massive(Vba[1], Vba[0], bajos, time, 1, 0); //Susceptible testeado a susceptible
  tested_massive(Val[3], Val[2], altos, time, 3, 2);  tested_massive(Vba[3], Vba[2], bajos, time, 3, 2); //Expuesto testeado a expuesto
  tested_massive(Val[8], Val[5], altos, time, 8, 5);  tested_massive(Vba[8], Vba[5], bajos, time, 8, 5); //Presintomático masivo a presintomátivo
  tested_massive(Val[12], Val[9], altos, time, 12, 9);  tested_massive(Vba[12], Vba[9], bajos, time, 12, 9); //Leve masivo a masivo
  tested_massive(Val[14], Val[15], altos, time, 14, 15);  tested_massive(Vba[14], Vba[15], bajos, time, 14, 15); //Recuperado testeado a recuperado
}
