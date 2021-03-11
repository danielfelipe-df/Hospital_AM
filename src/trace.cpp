#include <iostream>
#include <vector>
#include <algorithm>
#include <iterator>
#include <trace.h>
#include <test.h>

void trace_reaction(grupo &Out, grupo &In, int index, trabajadores *family, int typeout, int typein, double time, bool is_traced){
  int agent = Out[index];
  Out.erase(Out.begin() + index);
  In.push_back(agent);
  family[agent].change(typein, typeout);
  family[agent].ttrace = time;
  family[agent].btrace = is_traced;
}


void main_trace(std::vector<grupo> &Val, std::vector<grupo> &Vba, trabajadores *altos, trabajadores *bajos, double time, Crandom &ran){
  int num;
  num = tested_isolated_inf(Val[7], Val[8], Val[6], altos, time, 7, 8, 6, ran); //Presintomáticos altos
  aux_main(num, Val, Vba, altos, bajos, ran, MyCons.phi1, MyCons.mu, true, 8, true);

  num = tested_isolated_inf(Vba[7], Vba[8], Vba[6], bajos, time, 7, 8, 6, ran); //Presintomáticos bajos
  aux_main(num, Val, Vba, altos, bajos, ran, MyCons.mu, MyCons.chi, false, 8, true);

  num = tested_isolated_inf(Val[10], Val[11], Val[9], altos, time, 10, 11, 9, ran); //Leves altos
  aux_main(num, Val, Vba, altos, bajos, ran, MyCons.phi1, MyCons.mu, true, 11, false);

  num = tested_isolated_inf(Vba[10], Vba[11], Vba[9], bajos, time, 10, 11, 9, ran); //Leves bajos
  aux_main(num, Val, Vba, altos, bajos, ran, MyCons.mu, MyCons.chi, false, 11, false);
}


void aux_main(int num, std::vector<grupo> &Val, std::vector<grupo> &Vba, trabajadores *altos, trabajadores *bajos, Crandom &ran, double cons1, double cons2, bool type, unsigned int index, bool pre){
  unsigned int contador, aux, aux2, ind, my_size, Gsize, j;
  double value, aux3;
  grupo vaux;
  if(type){Gsize = Val[index].size();}
  else{Gsize = Vba[index].size();}

  for(unsigned int i=Gsize-num; i<Gsize; i++){
    if(type){my_size = altos[Val[index][i]].my_inf.size();}
    else{my_size = bajos[Vba[index][i]].my_inf.size();}

    for(j=0; j<my_size; j++){
      if(type){ind = altos[Val[index][i]].my_inf[j];}
      else{ind = bajos[Vba[index][i]].my_inf[j];}
      vaux.push_back(ind);
    }

    if(pre){
      if(type){aux3 = (int)(3*(MyCons.trace_net - my_size));}
      else{aux3 = (int)(3*(MyCons.trace_net - my_size));}
    }
    else{
      if(type){aux3 = (int)((altos[Val[index][i]].tstate + 2)*(MyCons.trace_net - my_size));}
      else{aux3 = (int)((bajos[Vba[index][i]].tstate + 2)*(MyCons.trace_net - my_size));}
    }

    while((0 < j) && (j < aux3)){
      value = ran.r()*(cons1*Na + cons2*Nb);
      if(value < cons1*Na){ind = (int)(ran.r()*Na);}
      else{ind = (int)(ran.r()*Nb) + Na;}
      vaux.push_back(ind);
      j++;
    }

    contador = 0;
    eliminar_repetidos(vaux);
    while(vaux.size() > 0 && contador < MyCons.trace){
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
  if(family[index].kind == 0){aux_trace(V[0], V[2], family, 0, 2, index, true);    return 1;} //Susceptible
  else if(family[index].kind == 1){aux_trace(V[1], V[2], family, 1, 2, index, true);    return 1;} //Susceptible testeado
  else if(family[index].kind == 3){aux_trace(V[3], V[5], family, 3, 5, index, true);    return 1;} //Expuesto
  else if(family[index].kind == 4){aux_trace(V[4], V[5], family, 4, 5, index, true);    return 1;} //Expuesto testeado
  else if(family[index].kind == 6){aux_trace(V[6], V[8], family, 6, 8, index, true);    return 1;} //Presintomático
  else if(family[index].kind == 7){aux_trace(V[7], V[8], family, 7, 8, index, true);    return 1;} //Presintomático testeado
  else if(family[index].kind == 9){aux_trace(V[9], V[11], family, 9, 11, index, true);    if(MyCons.AisLev){tested_lev_ais(V[11].back(), family, 1e6, false);}    return 1;} //Leve
  else if(family[index].kind == 10){aux_trace(V[10], V[11], family, 10, 11, index, true);    if(MyCons.AisLev){tested_lev_ais(V[11].back(), family, 1e6, false);}    return 1;} //Leve testeado
  else if(family[index].kind == 13){aux_trace(V[13], V[14], family, 13, 14, index, true);    return 1;} //Recuperado no-detectado
  else if(family[index].kind == 14){aux_trace(V[14], V[14], family, 14, 14, index, false);    return 1;} //Recuperado testeado
  else{return 0;}
}


void aux_trace(grupo &G, grupo &T, trabajadores *family, int typeout, int typein, int index, bool normal){
  std::vector<int>::iterator it;
  unsigned int ind;

  if(normal){
    it = std::find(G.begin(), G.end(), index);  ind = std::distance(G.begin(), it);
    trace_reaction(G, T, ind, family, typeout, typein, 0.0, true);
  }
  else{
    it = std::find(T.begin(), T.end(), index);  ind = std::distance(T.begin(), it);
    family[T[ind]].time = 0.0;
  }
}


void eliminar_repetidos(std::vector<int> &y)
{
  auto end = y.end();
  for(auto it = y.begin(); it != end; it++){
    end = std::remove(it + 1, end, *it);
  }
  y.erase(end, y.end());
}


void trace_massive(grupo &R, grupo &G, trabajadores *family, double time, int typeout, int typein){
  for(size_t i=0; i<R.size(); i++){
    if(family[R[i]].btrace){
      family[R[i]].ttrace += time;
      if(family[R[i]].ttrace > TtraceMax){trace_reaction(R, G, i, family, typeout, typein, 0.0, false);	i--;}
    }
  }
}


void trace_massive_all(std::vector<grupo> &Val, std::vector<grupo> &Vba, trabajadores *altos, trabajadores *bajos, double time){
  trace_massive(Val[2], Val[0], altos, time, 2, 0);  trace_massive(Vba[2], Vba[0], bajos, time, 2, 0); // Susceptible aislado a susceptible
  trace_massive(Val[5], Val[3], altos, time, 5, 3);  trace_massive(Vba[5], Vba[3], bajos, time, 5, 3); // Expuesto aislado a expuesto
  trace_massive(Val[8], Val[6], altos, time, 8, 6);  trace_massive(Vba[8], Vba[6], bajos, time, 8, 6); // Presintomático aislado a presintomático
  trace_massive(Val[11], Val[9], altos, time, 11, 9);  trace_massive(Vba[11], Vba[9], bajos, time, 11, 9); // Leve aislado a leve
  trace_massive(Val[14], Val[13], altos, time, 14, 13);  trace_massive(Vba[14], Vba[13], bajos, time, 14, 13); // Recuperado aislado a recuperado
}
