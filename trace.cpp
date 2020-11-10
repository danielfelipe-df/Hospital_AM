#include <iostream>
#include <vector>
#include <algorithm>
#include <iterator>
#include "trace.h"
#include "test.h"


void main_trace(std::vector<grupo> &Val, std::vector<grupo> &Vba, trabajadores *altos, trabajadores *bajos, double time, Crandom &ran){
  int num;
  num = tested_isolated_inf(Val[6], Val[7], Val[5], altos, time, 6, 7, 5, ran); //Presintomáticos altos
  aux_main(num, Val, Vba, altos, bajos, ran, phi1, mu, true, 7);

  num = tested_isolated_inf(Vba[6], Vba[7], Vba[5], bajos, time, 6, 7, 5, ran); //Presintomáticos bajos
  aux_main(num, Val, Vba, altos, bajos, ran, mu, chi, false, 7);

  num = tested_isolated_inf(Val[9], Val[10], Val[8], altos, time, 9, 10, 8, ran); //Leves altos
  aux_main(num, Val, Vba, altos, bajos, ran, phi1, mu, true, 10);

  num = tested_isolated_inf(Vba[9], Vba[10], Vba[8], bajos, time, 9, 10, 8, ran); //Leves bajos
  aux_main(num, Val, Vba, altos, bajos, ran, mu, chi, false, 10);
}


void aux_main(int num, std::vector<grupo> &Val, std::vector<grupo> &Vba, trabajadores *altos, trabajadores *bajos, Crandom &ran, double cons1, double cons2, bool type, unsigned int index){
  unsigned int contador, aux, ind, my_size, Gsize;
  double value1, value2, value3;
  if(type){Gsize = Val[index].size();}
  else{Gsize = Vba[index].size();}

  for(unsigned int i=Gsize-num; i<Gsize; i++){
    contador = 0;
    if(type){my_size = altos[Val[index][i]].my_inf.size();}
    else{my_size = bajos[Vba[index][i]].my_inf.size();}

    for(unsigned int j=0; j<my_size && contador<trace; j++){
      if(type){ind = altos[Val[index][i]].my_inf[j];}
      else{ind = bajos[Vba[index][i]].my_inf[j];}

      if(ind < Na){aux = reaction_trace(Val, altos, ran, ind);}
      else{aux = reaction_trace(Vba, bajos, ran, ind-Na);}
      contador += aux;
    }

    while(contador<trace){
      value1 = Val[4].size() + Val[7].size() + Val[10].size() + Val[11].size() + Val[14].size(); //Todos los aislados altos
      value2 = Vba[4].size() + Vba[7].size() + Vba[10].size() + Vba[11].size() + Vba[14].size(); //Todos los aislados bajos
      value3 = ran.r()*(cons1*(Na-value1) + cons2*(Nb-value2));
      if(value3 < cons1*(Na-value1)){aux = reaction_trace(Val, altos, ran, -1);}
      else{aux = reaction_trace(Vba, bajos, ran, -1);}
      contador += aux;
    }
  }
}


int reaction_trace(std::vector<grupo> &V, trabajadores *family, Crandom &ran, int index){
  int agent;
  if(index == -1){
    grupo aux;
    int ind;
    for(unsigned int i=0; i<=3; i++){std::copy(V[i].begin(), V[i].end(), std::back_inserter(aux));} //Susceptibles y expuestos
    std::copy(V[5].begin(), V[5].end(), std::back_inserter(aux)); //Presintomáticos
    std::copy(V[6].begin(), V[6].end(), std::back_inserter(aux)); //Presintomáticos testeados
    std::copy(V[8].begin(), V[8].end(), std::back_inserter(aux)); //Leves
    std::copy(V[9].begin(), V[9].end(), std::back_inserter(aux)); //Leves testeados
    std::copy(V[12].begin(), V[12].end(), std::back_inserter(aux)); //Recuperados no-detectados
    std::copy(V[13].begin(), V[13].end(), std::back_inserter(aux)); //Recuperados testeados

    ind = (int)(ran.r()*aux.size());
    agent = aux[ind];
    aux.clear();
  }
  else{
    agent = index;
  }

  if(family[agent].kind == 0){aux_trace(V[0], V[1], family, 0, 1, agent, true);    return 1;} //Susceptible
  else if(family[agent].kind == 1){aux_trace(V[0], V[1], family, 0, 1, agent, false);    return 1;} //Susceptible testeado
  else if(family[agent].kind == 2){aux_trace(V[2], V[4], family, 2, 4, agent, true);    return 1;} //Expuesto
  else if(family[agent].kind == 3){aux_trace(V[3], V[4], family, 3, 4, agent, true);    return 1;} //Expuesto testeado
  else if(family[agent].kind == 5){aux_trace(V[5], V[7], family, 5, 7, agent, true);    return 1;} //Presintomático
  else if(family[agent].kind == 6){aux_trace(V[6], V[7], family, 6, 7, agent, true);    return 1;} //Presintomático testeado
  else if(family[agent].kind == 8){aux_trace(V[8], V[10], family, 8, 10, agent, true);    return 1;} //Leve
  else if(family[agent].kind == 9){aux_trace(V[9], V[10], family, 9, 10, agent, true);    return 1;} //Leve testeado
  else if(family[agent].kind == 12){aux_trace(V[12], V[13], family, 12, 13, agent, true);    return 1;} //Recuperado no-detectado
  else if(family[agent].kind == 13){aux_trace(V[12], V[13], family, 12, 13, agent, false);    return 1;} //Recuperado testeado
  else{return 0;}
}


void aux_trace(grupo &G, grupo &T, trabajadores *family, int typeout, int typein, int index, bool normal){
  std::vector<int>::iterator it;
  unsigned int ind;

  if(normal){
    it = std::find(G.begin(), G.end(), index);  ind = std::distance(G.begin(), it);
    tested_reaction(G, T, ind, family, typeout, typein, 0.0);
  }
  else{
    it = std::find(T.begin(), T.end(), index);  ind = std::distance(T.begin(), it);
    family[T[ind]].time = 0.0;
  }
}
