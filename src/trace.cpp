#include <algorithm>
#include <iterator>
#include <trace.h>
#include <test.h>

void trace_reaction(grupo &Out, grupo &In, int index, Workers *family, Stages typeout, Stages typein, double time, bool is_traced){
  int agent = Out[index];
  Out.erase(Out.begin() + index);
  In.push_back(agent);
  family[agent].change(typein, typeout);
  family[agent].ttrace = time;
  family[agent].btrace = is_traced;
}


void main_trace(std::map<std::string, grupo> &Val, std::map<std::string, grupo> &Vba, Workers *altos, Workers *bajos, double time, Crandom &ran){
  int num;
  num = tested_isolated_inf(Val["PRET"], Val["PREA"], Val["PRE"], altos, time, PRET, PREA, PRE, ran); //Presintomáticos altos
  aux_main(num, Val, Vba, altos, bajos, ran, MyCons.phi1, MyCons.mu, true, "PREA", true);

  num = tested_isolated_inf(Vba["PRET"], Vba["PREA"], Vba["PRE"], bajos, time, PRET, PREA, PRE, ran); //Presintomáticos bajos
  aux_main(num, Val, Vba, altos, bajos, ran, MyCons.mu, MyCons.chi, false, "PREA", true);

  num = tested_isolated_inf(Val["MSYMT"], Val["MSYMA"], Val["MSYM"], altos, time, MSYMT, MSYMA, MSYM, ran); //Leves altos
  aux_main(num, Val, Vba, altos, bajos, ran, MyCons.phi1, MyCons.mu, true, "MSYMA", false);

  num = tested_isolated_inf(Vba["MSYMT"], Vba["MSYMA"], Vba["MSYM"], bajos, time, MSYMT, MSYMA, MSYM, ran); //Leves bajos
  aux_main(num, Val, Vba, altos, bajos, ran, MyCons.mu, MyCons.chi, false, "MSYMA", false);
}


void aux_main(int num, std::map<std::string, grupo> &Val, std::map<std::string, grupo> &Vba, Workers *altos, Workers *bajos, Crandom &ran, double cons1, double cons2, bool type, std::string s, bool pre){
  unsigned int contador, aux, aux2, ind, my_size, Gsize, j;
  double value, aux3;
  grupo vaux;
  if(type){Gsize = Val[s].size();}
  else{Gsize = Vba[s].size();}

  for(unsigned int i=Gsize-num; i<Gsize; i++){
    if(type){my_size = altos[Val[s][i]].my_inf.size();}
    else{my_size = bajos[Vba[s][i]].my_inf.size();}

    for(j=0; j<my_size; j++){
      if(type){ind = altos[Val[s][i]].my_inf[j];}
      else{ind = bajos[Vba[s][i]].my_inf[j];}
      vaux.push_back(ind);
    }

    if(pre){
      if(type){aux3 = (int)(3*(MyCons.trace_net - my_size));}
      else{aux3 = (int)(3*(MyCons.trace_net - my_size));}
    }
    else{
      if(type){aux3 = (int)((altos[Val[s][i]].tstate + 2)*(MyCons.trace_net - my_size));}
      else{aux3 = (int)((bajos[Vba[s][i]].tstate + 2)*(MyCons.trace_net - my_size));}
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


int reaction_trace(std::map<std::string, grupo> &V, Workers *family, Crandom &ran, int index){
  if(family[index].kind == SUS){
    aux_trace(V["SUS"], V["SUSA"], family, SUS, SUSA, index, true);
    return 1;
  } //Susceptible
  else if(family[index].kind == SUST){
    aux_trace(V["SUST"], V["SUSA"], family, SUST, SUSA, index, true);
    return 1;
  } //Susceptible testeado
  else if(family[index].kind == EXP){
    aux_trace(V["EXP"], V["EXPA"], family, EXP, EXPA, index, true);
    return 1;
  } //Expuesto
  else if(family[index].kind == EXPT){
    aux_trace(V["EXPT"], V["EXPA"], family, EXPT, EXPA, index, true);
    return 1;
  } //Expuesto testeado
  else if(family[index].kind == PRE){
    aux_trace(V["PRE"], V["PREA"], family, PRE, PREA, index, true);
    return 1;
  } //Presintomático
  else if(family[index].kind == PRET){
    aux_trace(V["PRET"], V["PREA"], family, PRET, PREA, index, true);
    return 1;
  } //Presintomático testeado
  else if(family[index].kind == MSYM){
    aux_trace(V["MSYM"], V["MSYMA"], family, MSYM, MSYMA, index, true);
    if(MyCons.AisLev){tested_lev_ais(V["MSYMA"].back(), family, 1e6, false);}
    return 1;
  } //Leve
  else if(family[index].kind == MSYMT){
    aux_trace(V["MSYMT"], V["MSYMA"], family, MSYMT, MSYMA, index, true);
    if(MyCons.AisLev){tested_lev_ais(V["MSYMA"].back(), family, 1e6, false);}
    return 1;
  } //Leve testeado
  else if(family[index].kind == RECI){
    aux_trace(V["RECI"], V["RECT"], family, RECI, RECT, index, true);
    return 1;
  } //Recuperado no-detectado
  else if(family[index].kind == RECT){
    aux_trace(V["RECT"], V["RECT"], family, RECT, RECT, index, false);
    return 1;
  } //Recuperado testeado
  else{
    return 0;
  }
}


void aux_trace(grupo &G, grupo &T, Workers *family, Stages typeout, Stages typein, int index, bool normal){
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


void trace_massive(grupo &R, grupo &G, Workers *family, double time, Stages typeout, Stages typein){
  for(size_t i=0; i<R.size(); i++){
    if(family[R[i]].btrace){
      family[R[i]].ttrace += time;
      if(family[R[i]].ttrace > TtraceMax){
	trace_reaction(R, G, i, family, typeout, typein, 0.0, false);
	i--;
      }
    }
  }
}


void trace_massive_all(std::map<std::string, grupo> &Val, std::map<std::string, grupo> &Vba, Workers *altos, Workers *bajos, double time){
  trace_massive(Val["SUSA"], Val["SUS"], altos, time, SUSA, SUS);  trace_massive(Vba["SUSA"], Vba["SUS"], bajos, time, SUSA, SUS); // Susceptible aislado a susceptible
  trace_massive(Val["EXPA"], Val["EXP"], altos, time, EXPA, EXP);  trace_massive(Vba["EXPA"], Vba["EXP"], bajos, time, EXPA, EXP); // Expuesto aislado a expuesto
  trace_massive(Val["PREA"], Val["PRE"], altos, time, PREA, PRE);  trace_massive(Vba["PREA"], Vba["PRE"], bajos, time, PREA, PRE); // Presintomático aislado a presintomático
  trace_massive(Val["MSYMA"], Val["MSYM"], altos, time, MSYMA, MSYM);  trace_massive(Vba["MSYMA"], Vba["MSYM"], bajos, time, MSYMA, MSYM); // Leve aislado a leve
  trace_massive(Val["RECT"], Val["RECI"], altos, time, RECT, RECI);  trace_massive(Vba["RECT"], Vba["RECI"], bajos, time, RECT, RECI); // Recuperado aislado a recuperado
}
