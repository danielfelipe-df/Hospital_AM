#include <vector>
#include <algorithm>
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
  if(type){Gsize = Val[index].size();}
  else{Gsize = Vba[index].size();}
  
  for(unsigned int i=Gsize; i<Gsize; i++){
    contador = 0;
    if(type){my_size = altos[Val[index][i]].my_inf.size();}
    else{my_size = bajos[Vba[index][i]].my_inf.size();}
    
    for(unsigned int j=0; j<my_size && contador<trace; j++){
      if(type){ind = altos[Val[index][i]].my_inf[j];}
      else{ind = bajos[Vba[index][i]].my_inf[j];}
      
      if(ind < Na){aux = reaction_trace(ind, Val, altos);}
      else{aux = reaction_trace(ind-Na, Vba, bajos);}
      contador += aux;
    }
    
    while(contador<trace){
      ind = (int)(ran.r()*(cons1*Na + cons2*Nb));
      if(ind < cons1*Na){aux = reaction_trace(ind, Val, altos);}
      else{aux = reaction_trace(ind-Na, Vba, bajos);}
      contador += aux;
    }
  }
}


int reaction_trace(int index, std::vector<grupo> &V, trabajadores *family){
  int num;

  num = aux_trace(V[0], V[1], family, 0, 1, index);
  if(num == 1){return num;}

  num = aux_trace(V[2], V[3], family, 2, 3, index);
  if(num == 1){return num;}

  num = aux_trace(V[4], V[5], family, 4, 5, index);
  if(num == 1){return num;}

  num = aux_trace(V[7], V[8], family, 7, 8, index);
  if(num == 1){return num;}

  num = aux_trace(V[11], V[12], family, 11, 12, index);
  if(num == 1){return num;}

  return 0;
}


int aux_trace(grupo &G, grupo &T, trabajadores *family, int typeout, int typein, int index){
  std::vector<int>::iterator it;
  unsigned int ind;
  
  if(G.size() != 0){
    it = std::find(G.begin(), G.end(), index);  ind = std::distance(G.begin(), it);
    if(ind < G.size()){tested_reaction(G, T, ind, family, typeout, typein, Tt);      return 1;}
  }
  if(T.size() != 0){
    it = std::find(T.begin(), T.end(), index);  ind = std::distance(T.begin(), it);
    if(ind < T.size()){family[T[ind]].time = 0.0;      return 1;}
  }

  return 0;
}
