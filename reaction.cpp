#include "reaction.h"
#include "trace.h"
#include "test.h"


void mother_reaction(grupo &Out, grupo &In, int index, trabajadores *family, int typeout, int typein){
  int agent = Out[index];
  Out.erase(Out.begin() + index);
  In.push_back(agent);
  family[agent].change(typein, typeout);
  family[agent].tstate = 0.0;
}


int index_time(grupo &Out, trabajadores *family, lognormal_d &dist, double value){
  //Jose
  unsigned int n = Out.size(), agent;
  double times[n];
  double param = value;

  //Hallo la diferencia del tiempo con el tiempo que lleve en el estado
  for(unsigned int i=0; i<n; i++){agent = Out[i];    times[i] = std::abs(cdf(dist, family[agent].tstate) - param);}

  //Retorno el índice donde está el tiempo mínimo
  return std::distance(times, std::min_element(times, times+n));
}


int who_infected(grupo &Pa, grupo &Pb, grupo &PTa, grupo &PTb, grupo &PTAa, grupo &PTAb, grupo &La, grupo &Lb, grupo &LTa, grupo &LTb, grupo &LTAa, grupo &LTAb, grupo &IAa, grupo &IAb, double cons1, double cons2, Crandom &ran, int index, trabajadores *altos, trabajadores *bajos){
  double num[4];
  num[0] = cons1*(Pa.size() + PTa.size() + La.size() + LTa.size())/(double)Na;
  num[1] = cons2*(Pb.size() + PTb.size() + Lb.size() + LTb.size())/(double)Nb;
  num[2] = (1-alpha)*cons1*(IAa.size() + PTAa.size() + LTAa.size())/(double)Na;
  num[3] = (1-alpha)*cons2*(IAb.size() + PTAb.size() + LTAb.size())/(double)Nb;
  grupo aux;
  if(num[0]+num[1]+num[2]+num[3] > 0.0){
    double num2 = ran.r()*(num[0] + num[1] + num[2] + num[3]);
    if(num2 < num[0]){return selection_infectious(Pa, PTa, La, LTa, ran, index, altos);}
    else if(num2 < num[0] + num[1]){return selection_infectious(Pb, PTb, Lb, LTb, ran, index, bajos) + Na;}
    else if(num2 < num[0] + num[1] + num[2]){return selection_infectious(IAa, PTAa, LTAa, aux, ran, index, altos);}
    else{return selection_infectious(IAb, PTAb, LTAb, aux, ran, index, bajos) + Na;}
  }
  else{
    return -1;
  }
}


int selection_infectious(grupo &Ga, grupo &Gb, grupo &Gc, grupo &Gd, Crandom &ran, int index, trabajadores *family){
  double num = ran.r()*(Ga.size() + Gb.size() + Gc.size() + Gd.size());
  int ind, agent;
  if(num < Ga.size()){ind = (int)(ran.r()*Ga.size());    agent = Ga[ind];    family[agent].my_inf.push_back(index);}
  else if(num < Ga.size() + Gb.size()){ind = (int)(ran.r()*Gb.size());    agent = Gb[ind];    family[agent].my_inf.push_back(index);}
  else if(num < Ga.size() + Gb.size() + Gc.size()){ind = (int)(ran.r()*Gc.size());    agent = Gc[ind];    family[agent].my_inf.push_back(index);}
  else{ind = (int)(ran.r()*Gd.size());    agent = Gd[ind];    family[agent].my_inf.push_back(index);}
  return agent;
}


void reaction0(std::vector<grupo> &Val, std::vector<grupo> &Vba, Crandom &ran, trabajadores *altos, trabajadores *bajos, std::vector<lognormal_d> &dist){
  unsigned int index = (int)(ran.r()*(Val[0].size() + Val[1].size()));
  int agentS, value1 = Val[2].size(), value2;
  if(index < Val[0].size()){agentS = Val[0][index];    mother_reaction(Val[0], Val[2], index, altos, 0, 2);}
  else{agentS = Val[1][index-Val[0].size()];    mother_reaction(Val[1], Val[3], index-Val[0].size(), altos, 1, 3);}
  value2 = Val[2].size();

  int agentI = who_infected(Val[4], Vba[4], Val[5], Vba[5], Val[6], Vba[6], Val[7], Vba[7], Val[8], Vba[8], Val[9], Vba[9], Val[10], Vba[10], phi1, mu, ran, agentS, altos, bajos);
  if(value2 > value1){altos[Val[2].back()].DF = agentI;}
  else{altos[Val[3].back()].DF = agentI;}
}


void reaction1(std::vector<grupo> &Val, std::vector<grupo> &Vba, Crandom &ran, trabajadores *altos, trabajadores *bajos, std::vector<lognormal_d> &dist){
  unsigned int index = (int)(ran.r()*(Vba[0].size() + Vba[1].size()));
  int agentS, value1 = Vba[2].size(), value2;
  if(index < Vba[0].size()){agentS = Vba[0][index];    mother_reaction(Vba[0], Vba[2], index, bajos, 0, 2);}
  else{agentS = Vba[1][index-Vba[0].size()];    mother_reaction(Vba[1], Vba[3], index-Vba[0].size(), bajos, 1, 3);}
  value2 = Vba[2].size();

  int agentI = who_infected(Val[4], Vba[4], Val[5], Vba[5], Val[6], Vba[6], Val[7], Vba[7], Val[8], Vba[8], Val[9], Vba[9], Val[10], Vba[10], mu, chi, ran, agentS + Na, altos, bajos);
  if(value2 > value1){bajos[Vba[2].back()].DF = agentI;}
  else{bajos[Vba[3].back()].DF = agentI;}
}


void reaction2(std::vector<grupo> &Val, std::vector<grupo> &Vba, Crandom &ran, trabajadores *altos, trabajadores *bajos, std::vector<lognormal_d> &dist){
  std::vector<int> aux;

  for(unsigned int i=0; i<Val[2].size(); i++){if(altos[Val[2][i]].tstate > TM){aux.push_back(Val[2][i]);}}
  for(unsigned int i=0; i<Val[3].size(); i++){if(altos[Val[3][i]].tstate > TM){aux.push_back(Val[3][i]);}}

  if(aux.size() != 0){
    double value = (double)ran.r();
    unsigned int index = index_time(aux, altos, dist[0], value);
    unsigned int agent = aux[index];
    std::vector<int>::iterator it;
    unsigned int ind;
    if(Val[2].size() != 0){
      it = std::find(Val[2].begin(), Val[2].end(), agent);      ind = std::distance(Val[2].begin(), it);
      if(ind < Val[2].size()){mother_reaction(Val[2], Val[4], ind, altos, 2, 4);}
    }
    if(Val[3].size() != 0){
      it = std::find(Val[3].begin(), Val[3].end(), agent);      ind = std::distance(Val[3].begin(), it);
      if(ind < Val[3].size()){mother_reaction(Val[3], Val[4], ind, altos, 3, 4);}
    }
  }
}


void reaction3(std::vector<grupo> &Val, std::vector<grupo> &Vba, Crandom &ran, trabajadores *altos, trabajadores *bajos, std::vector<lognormal_d> &dist){
  std::vector<int> aux;

  for(unsigned int i=0; i<Vba[2].size(); i++){if(bajos[Vba[2][i]].tstate > TM){aux.push_back(Vba[2][i]);}}
  for(unsigned int i=0; i<Vba[3].size(); i++){if(bajos[Vba[3][i]].tstate > TM){aux.push_back(Vba[3][i]);}}  

  if(aux.size() != 0){
    double value = (double)ran.r();
    unsigned int index = index_time(aux, bajos, dist[1], value);
    unsigned int agent = aux[index];
    std::vector<int>::iterator it;
    unsigned int ind;
    if(Vba[2].size() != 0){
      it = std::find(Vba[2].begin(), Vba[2].end(), agent);      ind = std::distance(Vba[2].begin(), it);
      if(ind < Vba[2].size()){mother_reaction(Vba[2], Vba[4], ind, bajos, 2, 4);}
    }
    if(Vba[3].size() != 0){
      it = std::find(Vba[3].begin(), Vba[3].end(), agent);      ind = std::distance(Vba[3].begin(), it);
      if(ind < Vba[3].size()){mother_reaction(Vba[3], Vba[4], ind, bajos, 3, 4);}
    }    
  }
}


void reaction4(std::vector<grupo> &Val, std::vector<grupo> &Vba, Crandom &ran, trabajadores *altos, trabajadores *bajos, std::vector<lognormal_d> &dist){
  std::vector<int> aux;

  for(unsigned int i=0; i<Val[4].size(); i++){if(altos[Val[4][i]].tstate > TM){aux.push_back(Val[4][i]);}}
  for(unsigned int i=0; i<Val[5].size(); i++){if(altos[Val[5][i]].tstate > TM){aux.push_back(Val[5][i]);}}
  for(unsigned int i=0; i<Val[6].size(); i++){if(altos[Val[6][i]].tstate > TM){aux.push_back(Val[6][i]);}}

  if(aux.size() != 0){
    double value = (double)ran.r();
    unsigned int index = index_time(aux, altos, dist[2], value);
    unsigned int agent = aux[index];
    std::vector<int>::iterator it;
    unsigned int ind;
    if(Val[4].size() != 0){
      it = std::find(Val[4].begin(), Val[4].end(), agent);      ind = std::distance(Val[4].begin(), it);
      if(ind < Val[4].size()){
	mother_reaction(Val[4], Val[7], ind, altos, 4, 7);
	if(ran.r() < iota){
	  tested_reaction(Val[7], Val[9], Val[7].size()-1, altos, 7, 9, 0.0);
	  aux_main(1, Val[9], Val[0], Vba[0], Val[1], Vba[1], Val[2], Vba[2], Val[3], Vba[3], Val[4], Vba[4], Val[5], Vba[5], Val[7], Vba[7], Val[8], Vba[8], Val[11], Vba[11], Val[12], Vba[12], altos, bajos, ran, phi1, mu, true);
	}
      }
    }
    if(Val[5].size() != 0){
      it = std::find(Val[5].begin(), Val[5].end(), agent);      ind = std::distance(Val[5].begin(), it);
      if(ind < Val[5].size()){mother_reaction(Val[5], Val[8], ind, altos, 5, 8);}
    }
    if(Val[6].size() != 0){
      it = std::find(Val[6].begin(), Val[6].end(), agent);      ind = std::distance(Val[6].begin(), it);
      if(ind < Val[6].size()){mother_reaction(Val[6], Val[9], ind, altos, 6, 9);}
    }
  }
}


void reaction5(std::vector<grupo> &Val, std::vector<grupo> &Vba, Crandom &ran, trabajadores *altos, trabajadores *bajos, std::vector<lognormal_d> &dist){
  std::vector<int> aux;

  for(unsigned int i=0; i<Vba[4].size(); i++){if(bajos[Vba[4][i]].tstate > TM){aux.push_back(Vba[4][i]);}}
  for(unsigned int i=0; i<Vba[5].size(); i++){if(bajos[Vba[5][i]].tstate > TM){aux.push_back(Vba[5][i]);}}
  for(unsigned int i=0; i<Vba[6].size(); i++){if(bajos[Vba[6][i]].tstate > TM){aux.push_back(Vba[6][i]);}}

  if(aux.size() != 0){
    double value = (double)ran.r();
    unsigned int index = index_time(aux, bajos, dist[3], value);
    unsigned int agent = aux[index];
    std::vector<int>::iterator it;
    unsigned int ind;
    if(Vba[4].size() != 0){
      it = std::find(Vba[4].begin(), Vba[4].end(), agent);      ind = std::distance(Vba[4].begin(), it);
      if(ind < Vba[4].size()){
	mother_reaction(Vba[4], Vba[7], ind, bajos, 4, 7);
	if(ran.r() < iota){
	  tested_reaction(Vba[7], Vba[9], Vba[7].size()-1, bajos, 7, 9, 0.0);
	  aux_main(1, Vba[9], Val[0], Vba[0], Val[1], Vba[1], Val[2], Vba[2], Val[3], Vba[3], Val[4], Vba[4], Val[5], Vba[5], Val[7], Vba[7], Val[8], Vba[8], Val[11], Vba[11], Val[12], Vba[12], altos, bajos, ran, mu, chi, false);
	}
      }
    }
    if(Vba[5].size() != 0){
      it = std::find(Vba[5].begin(), Vba[5].end(), agent);      ind = std::distance(Vba[5].begin(), it);
      if(ind < Vba[5].size()){mother_reaction(Vba[5], Vba[8], ind, bajos, 5, 8);}
    }
    if(Vba[6].size() != 0){
      it = std::find(Vba[6].begin(), Vba[6].end(), agent);      ind = std::distance(Vba[6].begin(), it);
      if(ind < Vba[6].size()){mother_reaction(Vba[6], Vba[9], ind, bajos, 6, 9);}
    }
  }
}


void reaction6(std::vector<grupo> &Val, std::vector<grupo> &Vba, Crandom &ran, trabajadores *altos, trabajadores *bajos, std::vector<lognormal_d> &dist){
  std::vector<int> aux;

  for(unsigned int i=0; i<Val[4].size(); i++){if(altos[Val[4][i]].tstate > TM){aux.push_back(Val[4][i]);}}
  for(unsigned int i=0; i<Val[5].size(); i++){if(altos[Val[5][i]].tstate > TM){aux.push_back(Val[5][i]);}}
  for(unsigned int i=0; i<Val[6].size(); i++){if(altos[Val[6][i]].tstate > TM){aux.push_back(Val[6][i]);}}

  if(aux.size() != 0){
    double value = (double)ran.r();
    unsigned int index = index_time(aux, altos, dist[4], value);
    unsigned int agent = aux[index];
    std::vector<int>::iterator it;
    unsigned int ind;
    if(Val[4].size() != 0){
      it = std::find(Val[4].begin(), Val[4].end(), agent);      ind = std::distance(Val[4].begin(), it);
      if(ind < Val[4].size()){mother_reaction(Val[4], Val[10], ind, altos, 4, 10);}
    }
    if(Val[5].size() != 0){
      it = std::find(Val[5].begin(), Val[5].end(), agent);      ind = std::distance(Val[5].begin(), it);
      if(ind < Val[5].size()){mother_reaction(Val[5], Val[10], ind, altos, 5, 10);}
    }
    if(Val[6].size() != 0){
      it = std::find(Val[6].begin(), Val[6].end(), agent);      ind = std::distance(Val[6].begin(), it);
      if(ind < Val[6].size()){mother_reaction(Val[6], Val[10], ind, altos, 6, 10);}
    }
    aux_main(1, Val[10], Val[0], Vba[0], Val[1], Vba[1], Val[2], Vba[2], Val[3], Vba[3], Val[4], Vba[4], Val[5], Vba[5], Val[7], Vba[7], Val[8], Vba[8], Val[11], Vba[11], Val[12], Vba[12], altos, bajos, ran, phi1, mu, true);
  }
}


void reaction7(std::vector<grupo> &Val, std::vector<grupo> &Vba, Crandom &ran, trabajadores *altos, trabajadores *bajos, std::vector<lognormal_d> &dist){
  std::vector<int> aux;

  for(unsigned int i=0; i<Vba[4].size(); i++){if(bajos[Vba[4][i]].tstate > TM){aux.push_back(Vba[4][i]);}}
  for(unsigned int i=0; i<Vba[5].size(); i++){if(bajos[Vba[5][i]].tstate > TM){aux.push_back(Vba[5][i]);}}
  for(unsigned int i=0; i<Vba[6].size(); i++){if(bajos[Vba[6][i]].tstate > TM){aux.push_back(Vba[6][i]);}}

  if(aux.size() != 0){
    double value = (double)ran.r();
    unsigned int index = index_time(aux, bajos, dist[5], value);
    unsigned int agent = aux[index];
    std::vector<int>::iterator it;
    unsigned int ind;
    if(Vba[4].size() != 0){
      it = std::find(Vba[4].begin(), Vba[4].end(), agent);      ind = std::distance(Vba[4].begin(), it);
      if(ind < Vba[4].size()){mother_reaction(Vba[4], Vba[10], ind, bajos, 4, 10);}
    }
    if(Vba[5].size() != 0){
      it = std::find(Vba[5].begin(), Vba[5].end(), agent);      ind = std::distance(Vba[5].begin(), it);
      if(ind < Vba[5].size()){mother_reaction(Vba[5], Vba[10], ind, bajos, 5, 10);}
    }
    if(Vba[6].size() != 0){
      it = std::find(Vba[6].begin(), Vba[6].end(), agent);      ind = std::distance(Vba[6].begin(), it);
      if(ind < Vba[6].size()){mother_reaction(Vba[6], Vba[10], ind, bajos, 6, 10);}
    }
    aux_main(1, Vba[10], Val[0], Vba[0], Val[1], Vba[1], Val[2], Vba[2], Val[3], Vba[3], Val[4], Vba[4], Val[5], Vba[5], Val[7], Vba[7], Val[8], Vba[8], Val[11], Vba[11], Val[12], Vba[12], altos, bajos, ran, mu, chi, false);
  }
}


void reaction8(std::vector<grupo> &Val, std::vector<grupo> &Vba, Crandom &ran, trabajadores *altos, trabajadores *bajos, std::vector<lognormal_d> &dist){
  std::vector<int> aux;

  for(unsigned int i=0; i<Val[4].size(); i++){if(altos[Val[4][i]].tstate > TM){aux.push_back(Val[4][i]);}}
  for(unsigned int i=0; i<Val[5].size(); i++){if(altos[Val[5][i]].tstate > TM){aux.push_back(Val[5][i]);}}
  for(unsigned int i=0; i<Val[6].size(); i++){if(altos[Val[6][i]].tstate > TM){aux.push_back(Val[6][i]);}}

  if(aux.size() != 0){
    double value = (double)ran.r();
    unsigned int index = index_time(aux, altos, dist[6], value);
    unsigned int agent = aux[index];
    std::vector<int>::iterator it;
    unsigned int ind;
    if(Val[4].size() != 0){
      it = std::find(Val[4].begin(), Val[4].end(), agent);      ind = std::distance(Val[4].begin(), it);
      if(ind < Val[4].size()){mother_reaction(Val[4], Val[11], ind, altos, 4, 11);}
    }
    if(Val[5].size() != 0){
      it = std::find(Val[5].begin(), Val[5].end(), agent);      ind = std::distance(Val[5].begin(), it);
      if(ind < Val[5].size()){mother_reaction(Val[5], Val[13], ind, altos, 5, 13);}
    }
    if(Val[6].size() != 0){
      it = std::find(Val[6].begin(), Val[6].end(), agent);      ind = std::distance(Val[6].begin(), it);
      if(ind < Val[6].size()){mother_reaction(Val[6], Val[13], ind, altos, 6, 13);}
    }
  }
}


void reaction9(std::vector<grupo> &Val, std::vector<grupo> &Vba, Crandom &ran, trabajadores *altos, trabajadores *bajos, std::vector<lognormal_d> &dist){
  std::vector<int> aux;

  for(unsigned int i=0; i<Vba[4].size(); i++){if(bajos[Vba[4][i]].tstate > TM){aux.push_back(Vba[4][i]);}}
  for(unsigned int i=0; i<Vba[5].size(); i++){if(bajos[Vba[5][i]].tstate > TM){aux.push_back(Vba[5][i]);}}
  for(unsigned int i=0; i<Vba[6].size(); i++){if(bajos[Vba[6][i]].tstate > TM){aux.push_back(Vba[6][i]);}}

  if(aux.size() != 0){
    double value = (double)ran.r();
    unsigned int index = index_time(aux, bajos, dist[7], value);
    unsigned int agent = aux[index];
    std::vector<int>::iterator it;
    unsigned int ind;
    if(Vba[4].size() != 0){
      it = std::find(Vba[4].begin(), Vba[4].end(), agent);      ind = std::distance(Vba[4].begin(), it);
      if(ind < Vba[4].size()){mother_reaction(Vba[4], Vba[11], ind, bajos, 4, 11);}
    }
    if(Vba[5].size() != 0){
      it = std::find(Vba[5].begin(), Vba[5].end(), agent);      ind = std::distance(Vba[5].begin(), it);
      if(ind < Vba[5].size()){mother_reaction(Vba[5], Vba[13], ind, bajos, 5, 13);}
    }
    if(Vba[6].size() != 0){
      it = std::find(Vba[6].begin(), Vba[6].end(), agent);      ind = std::distance(Vba[6].begin(), it);
      if(ind < Vba[6].size()){mother_reaction(Vba[6], Vba[13], ind, bajos, 6, 13);}
    }
  }
}


void reaction10(std::vector<grupo> &Val, std::vector<grupo> &Vba, Crandom &ran, trabajadores *altos, trabajadores *bajos, std::vector<lognormal_d> &dist){
  std::vector<int> aux;

  for(unsigned int i=0; i<Val[7].size(); i++){if(altos[Val[7][i]].tstate > TM){aux.push_back(Val[7][i]);}}
  for(unsigned int i=0; i<Val[8].size(); i++){if(altos[Val[8][i]].tstate > TM){aux.push_back(Val[8][i]);}}
  for(unsigned int i=0; i<Val[9].size(); i++){if(altos[Val[9][i]].tstate > TM){aux.push_back(Val[9][i]);}}

  if(aux.size() != 0){
    double value = (double)ran.r();
    unsigned int index = index_time(aux, altos, dist[8], value);
    unsigned int agent = aux[index];
    std::vector<int>::iterator it;
    unsigned int ind;
    if(Val[7].size() != 0){
      it = std::find(Val[7].begin(), Val[7].end(), agent);      ind = std::distance(Val[7].begin(), it);
      if(ind < Val[7].size()){mother_reaction(Val[7], Val[11], ind, altos, 7, 11);}
    }
    if(Val[8].size() != 0){
      it = std::find(Val[8].begin(), Val[8].end(), agent);      ind = std::distance(Val[8].begin(), it);
      if(ind < Val[8].size()){mother_reaction(Val[8], Val[13], ind, altos, 8, 13);}
    }
    if(Val[9].size() != 0){
      it = std::find(Val[9].begin(), Val[9].end(), agent);      ind = std::distance(Val[9].begin(), it);
      if(ind < Val[9].size()){mother_reaction(Val[9], Val[13], ind, altos, 9, 13);}
    }
  }
}


void reaction11(std::vector<grupo> &Val, std::vector<grupo> &Vba, Crandom &ran, trabajadores *altos, trabajadores *bajos, std::vector<lognormal_d> &dist){
  std::vector<int> aux;

  for(unsigned int i=0; i<Vba[7].size(); i++){if(bajos[Vba[7][i]].tstate > TM){aux.push_back(Vba[7][i]);}}
  for(unsigned int i=0; i<Vba[8].size(); i++){if(bajos[Vba[8][i]].tstate > TM){aux.push_back(Vba[8][i]);}}
  for(unsigned int i=0; i<Vba[9].size(); i++){if(bajos[Vba[9][i]].tstate > TM){aux.push_back(Vba[9][i]);}}

  if(aux.size() != 0){
    double value = (double)ran.r();
    unsigned int index = index_time(aux, bajos, dist[9], value);
    unsigned int agent = aux[index];
    std::vector<int>::iterator it;
    unsigned int ind;
    if(Vba[7].size() != 0){
      it = std::find(Vba[7].begin(), Vba[7].end(), agent);      ind = std::distance(Vba[7].begin(), it);
      if(ind < Vba[7].size()){mother_reaction(Vba[7], Vba[11], ind, bajos, 7, 11);}
    }
    if(Vba[8].size() != 0){
      it = std::find(Vba[8].begin(), Vba[8].end(), agent);      ind = std::distance(Vba[8].begin(), it);
      if(ind < Vba[8].size()){mother_reaction(Vba[8], Vba[13], ind, bajos, 8, 13);}
    }
    if(Vba[9].size() != 0){
      it = std::find(Vba[9].begin(), Vba[9].end(), agent);      ind = std::distance(Vba[9].begin(), it);
      if(ind < Vba[9].size()){mother_reaction(Vba[9], Vba[13], ind, bajos, 9, 13);}
    }
  }
}


void reaction12(std::vector<grupo> &Val, std::vector<grupo> &Vba, Crandom &ran, trabajadores *altos, trabajadores *bajos, std::vector<lognormal_d> &dist){
  std::vector<int> aux;
  
  for(unsigned int i=0; i<Val[10].size(); i++){if(altos[Val[10][i]].tstate > TM){aux.push_back(Val[10][i]);}}

  if(aux.size() != 0){    
    unsigned int index = index_time(aux, altos, dist[10], (double)ran.r());
    unsigned int agent = aux[index];
    std::vector<int>::iterator it = std::find(Val[10].begin(), Val[10].end(), agent);
    unsigned int ind = std::distance(Val[10].begin(), it);
    if(ind < Val[10].size()){mother_reaction(Val[10], Val[13], ind, altos, 10, 13);}
  }
}


void reaction13(std::vector<grupo> &Val, std::vector<grupo> &Vba, Crandom &ran, trabajadores *altos, trabajadores *bajos, std::vector<lognormal_d> &dist){
  std::vector<int> aux;
  
  for(unsigned int i=0; i<Vba[10].size(); i++){if(bajos[Vba[10][i]].tstate > TM){aux.push_back(Vba[10][i]);}}

  if(aux.size() != 0){
    unsigned int index = index_time(aux, bajos, dist[11], (double)ran.r());
    unsigned int agent = aux[index];
    std::vector<int>::iterator it = std::find(Vba[10].begin(), Vba[10].end(), agent);
    unsigned int ind = std::distance(Vba[10].begin(), it);
    if(ind < Vba[10].size()){mother_reaction(Vba[10], Vba[13], ind, bajos, 10, 13);}
  }  
}
