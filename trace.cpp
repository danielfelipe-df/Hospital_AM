#include <vector>
#include <algorithm>
#include "trace.h"
#include "test.h"


void main_trace(std::vector<grupo> &Val, std::vector<grupo> &Vba, trabajadores *altos, trabajadores *bajos, double time, Crandom &ran){
  int num;
  num = tested_isolated_inf(Val[5], Val[6], Val[4], altos, time, 5, 6, 4, ran);
  aux_main(num, Val[6], Val[0], Vba[0], Val[1], Vba[1], Val[2], Vba[2], Val[3], Vba[3], Val[4], Vba[4], Val[5], Vba[5], Val[7], Vba[7], Val[8], Vba[8], Val[11], Vba[11], Val[12], Vba[12], altos, bajos, ran, phi1, mu, true);

  num = tested_isolated_inf(Vba[5], Vba[6], Vba[4], bajos, time, 5, 6, 4, ran);
  aux_main(num, Vba[6], Val[0], Vba[0], Val[1], Vba[1], Val[2], Vba[2], Val[3], Vba[3], Val[4], Vba[4], Val[5], Vba[5], Val[7], Vba[7], Val[8], Vba[8], Val[11], Vba[11], Val[12], Vba[12], altos, bajos, ran, mu, chi, false);

  num = tested_isolated_inf(Val[8], Val[9], Val[7], altos, time, 8, 9, 7, ran);
  aux_main(num, Val[9], Val[0], Vba[0], Val[1], Vba[1], Val[2], Vba[2], Val[3], Vba[3], Val[4], Vba[4], Val[5], Vba[5], Val[7], Vba[7], Val[8], Vba[8], Val[11], Vba[11], Val[12], Vba[12], altos, bajos, ran, phi1, mu, true);

  num = tested_isolated_inf(Vba[8], Vba[9], Vba[7], bajos, time, 8, 9, 7, ran);
  aux_main(num, Vba[9], Val[0], Vba[0], Val[1], Vba[1], Val[2], Vba[2], Val[3], Vba[3], Val[4], Vba[4], Val[5], Vba[5], Val[7], Vba[7], Val[8], Vba[8], Val[11], Vba[11], Val[12], Vba[12], altos, bajos, ran, mu, chi, false);
}


void aux_main(int num, grupo &G, grupo &Sa, grupo &Sb, grupo &STa, grupo &STb, grupo &Ea, grupo &Eb, grupo &ETa, grupo &ETb, grupo &Pa, grupo &Pb, grupo &PTa, grupo &PTb, grupo &La, grupo &Lb, grupo &LTa, grupo &LTb, grupo &RIa, grupo &RIb, grupo &RTa, grupo &RTb, trabajadores *altos, trabajadores *bajos, Crandom &ran, double cons1, double cons2, bool type){
  unsigned int contador, aux, ind, my_size;
  for(unsigned int i=G.size()-num; i<G.size(); i++){
    contador = 0;
    if(type){my_size = altos[G[i]].my_inf.size();}
    else{my_size = bajos[G[i]].my_inf.size();}
    for(unsigned int j=0; j<my_size && contador<trace; j++){
      if(type){ind = altos[G[i]].my_inf[j];}
      else{ind = bajos[G[i]].my_inf[j];}
      if(ind < Na){aux = reaction_trace(ind, Sa, STa, Ea, ETa, Pa, PTa, La, LTa, RIa, RTa, altos);}
      else{aux = reaction_trace(ind-Na, Sb, STb, Eb, ETb, Pb, PTb, Lb, LTb, RIb, RTb, bajos);}
      contador += aux;
    }
    
    while(contador<trace){
      ind = (int)(ran.r()*(cons1*Na + cons2*Nb));
      if(ind < cons1*Na){aux = reaction_trace(ind, Sa, STa, Ea, ETa, Pa, PTa, La, LTa, RIa, RTa, altos);}
      else{aux = reaction_trace(ind-Na, Sb, STb, Eb, ETb, Pb, PTb, Lb, LTb, RIb, RTb, bajos);}
      contador += aux;
    }
  }
}


int reaction_trace(int index, grupo &S, grupo &ST, grupo &E, grupo &ET, grupo &P, grupo &PT, grupo &L, grupo &LT, grupo &RI, grupo &RT, trabajadores *family){
  int num;

  num = aux_trace(S, ST, family, 0, 1, index);
  if(num == 1){return num;}

  num = aux_trace(E, ET, family, 2, 3, index);
  if(num == 1){return num;}

  num = aux_trace(P, PT, family, 4, 5, index);
  if(num == 1){return num;}

  num = aux_trace(L, LT, family, 7, 8, index);
  if(num == 1){return num;}

  num = aux_trace(RI, RT, family, 11, 12, index);
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
