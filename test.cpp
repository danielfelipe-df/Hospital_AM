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
  unsigned int M = Val[0].size() + Vba[0].size() + Val[2].size() + Vba[2].size() + Val[4].size() + Vba[4].size() + Val[7].size() + Vba[7].size() + Val[11].size() + Vba[11].size();
  unsigned int num = (int)(ran.r()*M);
  
  if(num < Val[0].size()){
    tested_reaction(Val[0], Val[1], (int)(ran.r()*Val[0].size()), altos, 0, 1, Tt);
  }
  else if(num < Val[0].size() + Vba[0].size()){
    tested_reaction(Vba[0], Vba[1], (int)(ran.r()*Vba[0].size()), bajos, 0, 1, Tt);
  }
  else if(num < Val[0].size() + Vba[0].size() + Val[2].size()){
    tested_reaction(Val[2], Val[3], (int)(ran.r()*Val[2].size()), altos, 2, 3, Tt);
  }
  else if(num < Val[0].size() + Vba[0].size() + Val[2].size() + Vba[2].size()){
    tested_reaction(Vba[2], Vba[3], (int)(ran.r()*Vba[2].size()), bajos, 2, 3, Tt);
  }
  else if(num < Val[0].size() + Vba[0].size() + Val[2].size() + Vba[2].size() + Val[4].size()){
    tested_reaction(Val[4], Val[5], (int)(ran.r()*Val[4].size()), altos, 4, 5, Tt);
  }
  else if(num < Val[0].size() + Vba[0].size() + Val[2].size() + Vba[2].size() + Val[4].size() + Vba[4].size()){
    tested_reaction(Vba[4], Vba[5], (int)(ran.r()*Vba[4].size()), bajos, 4, 5, Tt);
  }
  else if(num < Val[0].size() + Vba[0].size() + Val[2].size() + Vba[2].size() + Val[4].size() + Vba[4].size() + Val[7].size()){
    tested_reaction(Val[7], Val[8], (int)(ran.r()*Val[7].size()), altos, 7, 8, Tt);
  }
  else if(num < Val[0].size() + Vba[0].size() + Val[2].size() + Vba[2].size() + Val[4].size() + Vba[4].size() + Val[7].size() + Vba[7].size()){
    tested_reaction(Vba[7], Vba[8], (int)(ran.r()*Vba[7].size()), bajos, 7, 8, Tt);
  }
  else if(num < Val[0].size() + Vba[0].size() + Val[2].size() + Vba[2].size() + Val[4].size() + Vba[4].size() + Val[7].size() + Vba[7].size() + Val[11].size()){
    tested_reaction(Val[11], Val[12], (int)(ran.r()*Val[11].size()), altos, 11, 12, Tt);
  }
  else{
    tested_reaction(Vba[11], Vba[12], (int)(ran.r()*Vba[11].size()), bajos, 11, 12, Tt);
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
  move_massive(Val[1], Val[0], altos, 1, 0);      move_massive(Vba[1], Vba[0], bajos, 1, 0);
  move_massive(Val[3], Val[2], altos, 3, 2);      move_massive(Vba[3], Vba[2], bajos, 3, 2);
  move_massive(Val[12], Val[11], altos, 12, 11);      move_massive(Vba[12], Vba[11], bajos, 12, 11);
}


void tested_massive_all(std::vector<grupo> &Val, std::vector<grupo> &Vba, trabajadores *altos, trabajadores *bajos, double time){
  tested_massive(Val[1], Val[0], altos, time, 1, 0);  tested_massive(Vba[1], Vba[0], bajos, time, 1, 0);
  tested_massive(Val[3], Val[2], altos, time, 3, 2);  tested_massive(Vba[3], Vba[2], bajos, time, 3, 2);
  tested_massive(Val[12], Val[11], altos, time, 12, 11);  tested_massive(Vba[12], Vba[11], bajos, time, 12, 11);
}
