#include "test.h"

void tested_reaction(grupo &Out, grupo &In, int index, trabajadores *family, int typeout, int typein, double delta){
  int agent = Out[index];
  Out.erase(Out.begin() + index);
  In.push_back(agent);
  family[agent].change(typein, typeout);
  family[agent].time = 0.0;
  family[agent].tmax = delta;
}

void massive_reaction(grupo &Sa, grupo &Sb, grupo &SMa, grupo &SMb, grupo &Ea, grupo &Eb, grupo &EMa, grupo &EMb, grupo &Pa, grupo &Pb, grupo &PTa, grupo &PTb, grupo &La, grupo &Lb, grupo &LTa, grupo &LTb, grupo &RIa, grupo &RIb, grupo &RMa, grupo &RMb, Crandom &ran, trabajadores *altos, trabajadores *bajos){
  unsigned int M = Sa.size() + Sb.size() + Ea.size() + Eb.size() + Pa.size() + Pb.size() + La.size() + Lb.size() + RIa.size() + RIb.size();
  unsigned int num = (int)(ran.r()*M);
  
  if(num < Sa.size()){
    tested_reaction(Sa, SMa, (int)(ran.r()*Sa.size()), altos, 0, 1, Tt);
  }
  else if(num < Sa.size() + Sb.size()){
    tested_reaction(Sb, SMb, (int)(ran.r()*Sb.size()), bajos, 0, 1, Tt);
  }
  else if(num < Sa.size() + Sb.size() + Ea.size()){
    tested_reaction(Ea, EMa, (int)(ran.r()*Ea.size()), altos, 2, 3, Tt);
  }
  else if(num < Sa.size() + Sb.size() + Ea.size() + Eb.size()){
    tested_reaction(Eb, EMb, (int)(ran.r()*Eb.size()), bajos, 2, 3, Tt);
  }
  else if(num < Sa.size() + Sb.size() + Ea.size() + Eb.size() + Pa.size()){
    tested_reaction(Pa, PTa, (int)(ran.r()*Pa.size()), altos, 4, 5, Tt);
  }
  else if(num < Sa.size() + Sb.size() + Ea.size() + Eb.size() + Pa.size() + Pb.size()){
    tested_reaction(Pb, PTb, (int)(ran.r()*Pb.size()), bajos, 4, 5, Tt);
  }
  else if(num < Sa.size() + Sb.size() + Ea.size() + Eb.size() + Pa.size() + Pb.size() + La.size()){
    tested_reaction(La, LTa, (int)(ran.r()*La.size()), altos, 7, 8, Tt);
  }
  else if(num < Sa.size() + Sb.size() + Ea.size() + Eb.size() + Pa.size() + Pb.size() + La.size() + Lb.size()){
    tested_reaction(Lb, LTb, (int)(ran.r()*Lb.size()), bajos, 7, 8, Tt);
  }
  else if(num < Sa.size() + Sb.size() + Ea.size() + Eb.size() + Pa.size() + Pb.size() + La.size() + Lb.size() + RIa.size()){
    tested_reaction(RIa, RMa, (int)(ran.r()*RIa.size()), altos, 11, 12, Tt);
  }
  else{
    tested_reaction(RIb, RMb, (int)(ran.r()*RIb.size()), bajos, 11, 12, Tt);
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
