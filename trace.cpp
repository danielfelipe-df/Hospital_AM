#include <iostream>
#include <vector>
#include <algorithm>
#include <iterator>
#include "trace.h"
#include "test.h"


void main_trace(std::vector<grupo> &Val, std::vector<grupo> &Vba, trabajadores *altos, trabajadores *bajos, double time, Crandom &ran){
  int num;
  num = tested_isolated_inf(Val[5], Val[6], Val[4], altos, time, 5, 6, 4, ran);
  aux_main(num, Val, Vba, altos, bajos, ran, phi1, mu, true, 6);

  num = tested_isolated_inf(Vba[5], Vba[6], Vba[4], bajos, time, 5, 6, 4, ran);
  aux_main(num, Val, Vba, altos, bajos, ran, mu, chi, false, 6);

  num = tested_isolated_inf(Val[8], Val[9], Val[7], altos, time, 8, 9, 7, ran);
  aux_main(num, Val, Vba, altos, bajos, ran, phi1, mu, true, 9);

  num = tested_isolated_inf(Vba[8], Vba[9], Vba[7], bajos, time, 8, 9, 7, ran);
  aux_main(num, Val, Vba, altos, bajos, ran, mu, chi, false, 9);
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
      value1 = Val[6].size() + Val[9].size() + Val[10].size() + Val[13].size();
      value2 = Vba[6].size() + Vba[9].size() + Vba[10].size() + Vba[13].size();
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
    for(unsigned int i=0; i<6; i++){std::copy(V[i].begin(), V[i].end(), std::back_inserter(aux));}
    std::copy(V[7].begin(), V[7].end(), std::back_inserter(aux));
    std::copy(V[8].begin(), V[8].end(), std::back_inserter(aux));
    std::copy(V[11].begin(), V[11].end(), std::back_inserter(aux));
    std::copy(V[12].begin(), V[12].end(), std::back_inserter(aux));
    
    ind = (int)(ran.r()*aux.size());
    agent = aux[ind];
    aux.clear();
  }
  else{
    agent = index;
  }

  if(family[agent].kind == 0){std::cout << 1 << std::endl;    aux_trace(V[0], V[1], family, 0, 1, agent, true);    return 1;}
  else if(family[agent].kind == 1){std::cout << 1 << std::endl;    aux_trace(V[0], V[1], family, 0, 1, agent, false);    return 1;}
  else if(family[agent].kind == 2){std::cout << 2 << std::endl;    aux_trace(V[2], V[3], family, 2, 3, agent, true);    return 1;}
  else if(family[agent].kind == 3){std::cout << 2 << std::endl;    aux_trace(V[2], V[3], family, 2, 3, agent, false);    return 1;}
  else if(family[agent].kind == 4){std::cout << 3 << std::endl;    aux_trace(V[4], V[6], family, 4, 6, agent, true);    return 1;}
  else if(family[agent].kind == 5){std::cout << 3 << std::endl;    aux_trace(V[5], V[6], family, 5, 6, agent, true);    return 1;}
  else if(family[agent].kind == 7){std::cout << 4 << std::endl;    aux_trace(V[7], V[9], family, 7, 9, agent, true);    return 1;}
  else if(family[agent].kind == 8){std::cout << 4 << std::endl;    aux_trace(V[8], V[9], family, 8, 9, agent, true);    return 1;}
  else if(family[agent].kind == 11){std::cout << 5 << std::endl;    aux_trace(V[11], V[12], family, 11, 12, agent, true);    return 1;}
  else if(family[agent].kind == 12){std::cout << 5 << std::endl;    aux_trace(V[11], V[12], family, 11, 12, agent, false);    return 1;}
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
