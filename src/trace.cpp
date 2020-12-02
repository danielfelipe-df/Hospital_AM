#include <iostream>
#include <vector>
#include <algorithm>
#include <iterator>
#include <trace.h>
#include <test.h>


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
  unsigned int contador, aux, aux2, ind, my_size, Gsize, j;
  double value;
  grupo vaux;
  if(type){Gsize = Val[index].size();}
  else{Gsize = Vba[index].size();}

  for(unsigned int i=Gsize-num; i<Gsize; i++){
    if(type){my_size = altos[Val[index][i]].my_inf.size();}
    else{my_size = bajos[Vba[index][i]].my_inf.size();}

    for(j=0; j<my_size && j<trace_net; j++){
      if(type){ind = altos[Val[index][i]].my_inf[j];}
      else{ind = bajos[Vba[index][i]].my_inf[j];}
      vaux.push_back(ind);
    }

    while(j<trace_net){
      value = ran.r()*(cons1*Na + cons2*Nb);
      if(value < cons1*Na){ind = (int)(ran.r()*Na);}
      else{ind = (int)(ran.r()*Nb) + Na;}
      vaux.push_back(ind);
      j++;
    }

    contador = 0;
    while(vaux.size() > 0 && contador < trace){
      aux2 = (int)(ran.r()*vaux.size());
      ind = vaux[aux2];
      if(ind < Na){aux = reaction_trace(Val, altos, ran, ind);}
      else{aux = reaction_trace(Vba, bajos, ran, ind-Na);}
      contador += aux;
      vaux.erase(vaux.begin() + aux2);
    }
    vaux.clear();
  }
}


int reaction_trace(std::vector<grupo> &V, trabajadores *family, Crandom &ran, int index){
  if(family[index].kind == 0){aux_trace(V[0], V[1], family, 0, 1, index, true);    return 1;} //Susceptible
  else if(family[index].kind == 1){aux_trace(V[0], V[1], family, 0, 1, index, false);    return 1;} //Susceptible testeado
  else if(family[index].kind == 2){aux_trace(V[2], V[4], family, 2, 4, index, true);    return 1;} //Expuesto
  else if(family[index].kind == 3){aux_trace(V[3], V[4], family, 3, 4, index, true);    return 1;} //Expuesto testeado
  else if(family[index].kind == 5){aux_trace(V[5], V[7], family, 5, 7, index, true);    return 1;} //Presintomático
  else if(family[index].kind == 6){aux_trace(V[6], V[7], family, 6, 7, index, true);    return 1;} //Presintomático testeado
  else if(family[index].kind == 8){aux_trace(V[8], V[10], family, 8, 10, index, true);    tested_lev_ais(V[10].back(), family, 1e6, false);    return 1;} //Leve
  else if(family[index].kind == 9){aux_trace(V[9], V[10], family, 9, 10, index, true);    tested_lev_ais(V[10].back(), family, 1e6, false);    return 1;} //Leve testeado
  else if(family[index].kind == 12){aux_trace(V[12], V[13], family, 12, 13, index, true);    return 1;} //Recuperado no-detectado
  else if(family[index].kind == 13){aux_trace(V[12], V[13], family, 12, 13, index, false);    return 1;} //Recuperado testeado
  else{return 0;}
}


void aux_trace(grupo &G, grupo &T, trabajadores *family, int typeout, int typein, int index, bool normal){
  std::vector<int>::iterator it;
  unsigned int ind;

  if(normal){
    it = std::find(G.begin(), G.end(), index);  ind = std::distance(G.begin(), it);
    tested_reaction(G, T, ind, family, typeout, typein, 0.0, false);
  }
  else{
    it = std::find(T.begin(), T.end(), index);  ind = std::distance(T.begin(), it);
    family[T[ind]].time = 0.0;
  }
}
