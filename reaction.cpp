#include <iostream>
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

int who_infected(grupo &Pa, grupo &Pb, grupo &PTa, grupo &PTb, grupo &PAa, grupo &PAb, grupo &PMa, grupo &PMb, grupo &La, grupo &Lb, grupo &LTa, grupo &LTb, grupo &LAa, grupo &LAb, grupo &LMa, grupo &LMb, grupo &IAa, grupo &IAb, double cons1, double cons2, Crandom &ran, int index, trabajadores *altos, trabajadores *bajos){
  double num[4];
  num[0] = cons1*(Pa.size() + PTa.size() + PMa.size() + La.size() + LTa.size() + LMa.size())/(double)Na;
  num[1] = cons2*(Pb.size() + PTb.size() + PMb.size() + Lb.size() + LTb.size() + LMb.size())/(double)Nb;
  num[2] = (1-alpha)*cons1*(IAa.size() + PAa.size() + LAa.size())/(double)Na;
  num[3] = (1-alpha)*cons2*(IAb.size() + PAb.size() + LAb.size())/(double)Nb;
  grupo aux1, aux2, aux3;
  if(num[0]+num[1]+num[2]+num[3] > 0.0){
    double num2 = ran.r()*(num[0] + num[1] + num[2] + num[3]);
    if(num2 < num[0]){return selection_infectious(Pa, PTa, PMa, La, LTa, LMa, ran, index, altos);}
    else if(num2 < num[0] + num[1]){return selection_infectious(Pb, PTb, PMb, Lb, LTb, LMb, ran, index, bajos) + Na;}
    else if(num2 < num[0] + num[1] + num[2]){return selection_infectious(IAa, PAa, LAa, aux1, aux2, aux3, ran, index, altos);}
    else{return selection_infectious(IAb, PAb, LAb, aux1, aux2, aux3, ran, index, bajos) + Na;}
  }
  else{
    return -1;
  }
}


int selection_infectious(grupo &Ga, grupo &Gb, grupo &Gc, grupo &Gd, grupo &Ge, grupo &Gf, Crandom &ran, int index, trabajadores *family){
  double num = ran.r()*(Ga.size() + Gb.size() + Gc.size() + Gd.size() + Ge.size() + Gf.size());
  int ind, agent;
  if(num < Ga.size()){ind = (int)(ran.r()*Ga.size());    agent = Ga[ind];}
  else if(num < Ga.size() + Gb.size()){ind = (int)(ran.r()*Gb.size());    agent = Gb[ind];}
  else if(num < Ga.size() + Gb.size() + Gc.size()){ind = (int)(ran.r()*Gc.size());    agent = Gc[ind];}
  else if(num < Ga.size() + Gb.size() + Gc.size() + Gd.size()){ind = (int)(ran.r()*Gd.size());    agent = Gd[ind];}
  else if(num < Ga.size() + Gb.size() + Gc.size() + Gd.size() + Ge.size()){ind = (int)(ran.r()*Ge.size());    agent = Ge[ind];}
  else{ind = (int)(ran.r()*Gf.size());    agent = Gf[ind];}
  family[agent].my_inf.push_back(index);
  return agent;
}


/* Susceptible a exupuesto. Alto */
void reaction0(std::vector<grupo> &Val, std::vector<grupo> &Vba, Crandom &ran, trabajadores *altos, trabajadores *bajos, std::vector<lognormal_d> &dist){
  unsigned int index = (int)(ran.r()*(Val[0].size() + Val[1].size()));
  int agentS, value1 = Val[2].size(), value2;
  /* Susceptible a Expuesto */
  if(index < Val[0].size()){agentS = Val[0][index];    mother_reaction(Val[0], Val[2], index, altos, 0, 2);}
  /* Susceptible testeado a Expuesto testado */
  else{agentS = Val[1][index-Val[0].size()];    mother_reaction(Val[1], Val[3], index-Val[0].size(), altos, 1, 3);}
  value2 = Val[2].size();

  int agentI = who_infected(Val[5], Vba[5], Val[6], Vba[6], Val[7], Vba[7], Val[8], Vba[8], Val[9], Vba[9], Val[10], Vba[10], Val[11], Vba[11], Val[12], Vba[12], Val[13], Vba[13], phi1, mu, ran, agentS, altos, bajos);
  if(value2 > value1){altos[Val[2].back()].DF = agentI;}
  else{altos[Val[3].back()].DF = agentI;}
  //std::cout << agentI << '\t' << agentS << std::endl;
}


/* Susceptible a exupuesto. Bajo */
void reaction1(std::vector<grupo> &Val, std::vector<grupo> &Vba, Crandom &ran, trabajadores *altos, trabajadores *bajos, std::vector<lognormal_d> &dist){
  unsigned int index = (int)(ran.r()*(Vba[0].size() + Vba[1].size()));
  int agentS, value1 = Vba[2].size(), value2;
  /* Susceptible a Expuesto */
  if(index < Vba[0].size()){agentS = Vba[0][index];    mother_reaction(Vba[0], Vba[2], index, bajos, 0, 2);}
  /* Susceptible testeado a Expuesto testeado */
  else{agentS = Vba[1][index-Vba[0].size()];    mother_reaction(Vba[1], Vba[3], index-Vba[0].size(), bajos, 1, 3);}
  value2 = Vba[2].size();

  int agentI = who_infected(Val[5], Vba[5], Val[6], Vba[6], Val[7], Vba[7], Val[8], Vba[8], Val[9], Vba[9], Val[10], Vba[10], Val[11], Vba[11], Val[12], Vba[12], Val[13], Vba[13], mu, chi, ran, agentS + Na, altos, bajos);  
  if(value2 > value1){bajos[Vba[2].back()].DF = agentI;}
  else{bajos[Vba[3].back()].DF = agentI;}
  //std::cout << agentI << '\t' << agentS+Na << std::endl;
}


/* Expuesto a Pre-sintomático. Alto */
void reaction2(std::vector<grupo> &Val, std::vector<grupo> &Vba, Crandom &ran, trabajadores *altos, trabajadores *bajos, std::vector<lognormal_d> &dist){
  std::vector<int> aux;

  for(unsigned int i=0; i<Val[2].size(); i++){if(altos[Val[2][i]].tstate > TM){aux.push_back(Val[2][i]);}} //Expuesto
  for(unsigned int i=0; i<Val[3].size(); i++){if(altos[Val[3][i]].tstate > TM){aux.push_back(Val[3][i]);}} //Expuesto testeado
  for(unsigned int i=0; i<Val[4].size(); i++){if(altos[Val[4][i]].tstate > TM){aux.push_back(Val[4][i]);}} //Expuesto aislado

  if(aux.size() != 0){
    double value = (double)ran.r();
    unsigned int index = index_time(aux, altos, dist[0], value);
    unsigned int agent = aux[index];
    std::vector<int>::iterator it;
    unsigned int ind;
    if(Val[2].size() != 0){// Expuesto a Pre-sintomático
      it = std::find(Val[2].begin(), Val[2].end(), agent);      ind = std::distance(Val[2].begin(), it);
      if(ind < Val[2].size()){mother_reaction(Val[2], Val[5], ind, altos, 2, 5);}
    }
    if(Val[3].size() != 0){// Expuesto testeado a Pre-sintomático masivo
      it = std::find(Val[3].begin(), Val[3].end(), agent);      ind = std::distance(Val[3].begin(), it);
      if(ind < Val[3].size()){mother_reaction(Val[3], Val[8], ind, altos, 3, 8);}
    }
    if(Val[4].size() != 0){// Expuesto aislado a Pre-sintomático aislado
      it = std::find(Val[4].begin(), Val[4].end(), agent);      ind = std::distance(Val[4].begin(), it);
      if(ind < Val[4].size()){mother_reaction(Val[4], Val[7], ind, altos, 4, 7);}
    }
  }
}


/* Expuesto a Pre-sintomático. Bajo */
void reaction3(std::vector<grupo> &Val, std::vector<grupo> &Vba, Crandom &ran, trabajadores *altos, trabajadores *bajos, std::vector<lognormal_d> &dist){
  std::vector<int> aux;

  for(unsigned int i=0; i<Vba[2].size(); i++){if(bajos[Vba[2][i]].tstate > TM){aux.push_back(Vba[2][i]);}} //Expuesto
  for(unsigned int i=0; i<Vba[3].size(); i++){if(bajos[Vba[3][i]].tstate > TM){aux.push_back(Vba[3][i]);}} //Expuesto testeado
  for(unsigned int i=0; i<Vba[4].size(); i++){if(bajos[Vba[4][i]].tstate > TM){aux.push_back(Vba[4][i]);}} //Expuesto aislado

  if(aux.size() != 0){
    double value = (double)ran.r();
    unsigned int index = index_time(aux, bajos, dist[1], value);
    unsigned int agent = aux[index];
    std::vector<int>::iterator it;
    unsigned int ind;
    if(Vba[2].size() != 0){// Expuesto a Pre-sintomático
      it = std::find(Vba[2].begin(), Vba[2].end(), agent);      ind = std::distance(Vba[2].begin(), it);
      if(ind < Vba[2].size()){mother_reaction(Vba[2], Vba[5], ind, bajos, 2, 5);}
    }
    if(Vba[3].size() != 0){// Expuesto testeado a Pre-sintomático masivo
      it = std::find(Vba[3].begin(), Vba[3].end(), agent);      ind = std::distance(Vba[3].begin(), it);
      if(ind < Vba[3].size()){mother_reaction(Vba[3], Vba[8], ind, bajos, 3, 8);}
    }
    if(Vba[4].size() != 0){// Expuesto aislado a Pre-sintomático aislado
      it = std::find(Vba[4].begin(), Vba[4].end(), agent);      ind = std::distance(Vba[4].begin(), it);
      if(ind < Vba[4].size()){mother_reaction(Vba[4], Vba[7], ind, bajos, 4, 7);}
    }
  }
}


/* Pre-sintomático a Leve. Alto */
void reaction4(std::vector<grupo> &Val, std::vector<grupo> &Vba, Crandom &ran, trabajadores *altos, trabajadores *bajos, std::vector<lognormal_d> &dist){
  std::vector<int> aux;

  for(unsigned int i=0; i<Val[5].size(); i++){if(altos[Val[5][i]].tstate > TM){aux.push_back(Val[5][i]);}} //Presintomático
  for(unsigned int i=0; i<Val[6].size(); i++){if(altos[Val[6][i]].tstate > TM){aux.push_back(Val[6][i]);}} //Presintomático testeado
  for(unsigned int i=0; i<Val[7].size(); i++){if(altos[Val[7][i]].tstate > TM){aux.push_back(Val[7][i]);}} //Presintomático aislado
  for(unsigned int i=0; i<Val[8].size(); i++){if(altos[Val[8][i]].tstate > TM){aux.push_back(Val[8][i]);}} //Presintomático masivo

  if(aux.size() != 0){
    double value = (double)ran.r();
    unsigned int index = index_time(aux, altos, dist[2], value);
    unsigned int agent = aux[index];
    std::vector<int>::iterator it;
    unsigned int ind;
    if(Val[5].size() != 0){// Presintomático a Leve
      it = std::find(Val[5].begin(), Val[5].end(), agent);      ind = std::distance(Val[5].begin(), it);
      if(ind < Val[5].size()){
	mother_reaction(Val[5], Val[9], ind, altos, 5, 9);
	if(ran.r() < iota){// Leve a Leve aislado
	  tested_reaction(Val[9], Val[11], Val[9].size()-1, altos, 9, 11, 0.0);
	  aux_main(1, Val, Vba, altos, bajos, ran, phi1, mu, true, 11);
	}
      }
    }
    if(Val[6].size() != 0){// Presintomático testeado a Leve testeado
      it = std::find(Val[6].begin(), Val[6].end(), agent);      ind = std::distance(Val[6].begin(), it);
      if(ind < Val[6].size()){
	mother_reaction(Val[6], Val[10], ind, altos, 6, 10);
	if(ran.r() < iota){// Leve testeado a Leve aislado
	  tested_reaction(Val[10], Val[11], Val[10].size()-1, altos, 10, 11, 0.0);
	  aux_main(1, Val, Vba, altos, bajos, ran, phi1, mu, true, 11);
	}
      }
    }
    if(Val[7].size() != 0){// Presintomático aislado a Leve aislado
      it = std::find(Val[7].begin(), Val[7].end(), agent);      ind = std::distance(Val[7].begin(), it);
      if(ind < Val[7].size()){mother_reaction(Val[7], Val[10], ind, altos, 7, 10);}
    }
    if(Val[8].size() != 0){// Presintomático masivo a Leve masivo
      it = std::find(Val[8].begin(), Val[8].end(), agent);      ind = std::distance(Val[8].begin(), it);
      if(ind < Val[8].size()){
	mother_reaction(Val[8], Val[12], ind, altos, 8, 12);
	if(ran.r() < iota){// Leve masivo a Leve aislado
	  tested_reaction(Val[12], Val[11], Val[12].size()-1, altos, 12, 11, 0.0);
	  aux_main(1, Val, Vba, altos, bajos, ran, phi1, mu, true, 11);
	}
      }
    }
  }
}


/* Pre-sintomático a Leve. Bajo */
void reaction5(std::vector<grupo> &Val, std::vector<grupo> &Vba, Crandom &ran, trabajadores *altos, trabajadores *bajos, std::vector<lognormal_d> &dist){
  std::vector<int> aux;

  for(unsigned int i=0; i<Vba[5].size(); i++){if(bajos[Vba[5][i]].tstate > TM){aux.push_back(Vba[5][i]);}} //Presintomático
  for(unsigned int i=0; i<Vba[6].size(); i++){if(bajos[Vba[6][i]].tstate > TM){aux.push_back(Vba[6][i]);}} //Presintomático testeado
  for(unsigned int i=0; i<Vba[7].size(); i++){if(bajos[Vba[7][i]].tstate > TM){aux.push_back(Vba[7][i]);}} //Presintomático aislado
  for(unsigned int i=0; i<Vba[8].size(); i++){if(bajos[Vba[8][i]].tstate > TM){aux.push_back(Vba[8][i]);}} //Presintomático masivo

  if(aux.size() != 0){
    double value = (double)ran.r();
    unsigned int index = index_time(aux, bajos, dist[3], value);
    unsigned int agent = aux[index];
    std::vector<int>::iterator it;
    unsigned int ind;
    if(Vba[5].size() != 0){// Presintomático a Leve
      it = std::find(Vba[5].begin(), Vba[5].end(), agent);      ind = std::distance(Vba[5].begin(), it);
      if(ind < Vba[5].size()){
	mother_reaction(Vba[5], Vba[9], ind, bajos, 5, 9);
	if(ran.r() < iota){// Leve a Leve aislado
	  tested_reaction(Vba[9], Vba[11], Vba[9].size()-1, bajos, 9, 11, 0.0);
	  aux_main(1, Val, Vba, altos, bajos, ran, mu, chi, false, 11);
	}
      }
    }
    if(Vba[6].size() != 0){// Presintomático testeado a Leve testeado
      it = std::find(Vba[6].begin(), Vba[6].end(), agent);      ind = std::distance(Vba[6].begin(), it);
      if(ind < Vba[6].size()){
	mother_reaction(Vba[6], Vba[10], ind, bajos, 6, 10);
	if(ran.r() < iota){// Leve testeado a Leve aislado
	  tested_reaction(Vba[10], Vba[11], Vba[10].size()-1, bajos, 10, 11, 0.0);
	  aux_main(1, Val, Vba, altos, bajos, ran, mu, chi, false, 11);
	}
      }
    }
    if(Vba[7].size() != 0){// Presintomático aislado a Leve aislado
      it = std::find(Vba[7].begin(), Vba[7].end(), agent);      ind = std::distance(Vba[7].begin(), it);
      if(ind < Vba[7].size()){mother_reaction(Vba[7], Vba[11], ind, bajos, 7, 11);}
    }
    if(Vba[8].size() != 0){// Presintomático masivo a Leve masivo
      it = std::find(Vba[8].begin(), Vba[8].end(), agent);      ind = std::distance(Vba[8].begin(), it);
      if(ind < Vba[8].size()){
	mother_reaction(Vba[8], Vba[12], ind, bajos, 8, 12);
	if(ran.r() < iota){// Leve masivo a Leve aislado
	  tested_reaction(Vba[12], Vba[11], Vba[12].size()-1, bajos, 12, 11, 0.0);
	  aux_main(1, Val, Vba, altos, bajos, ran, mu, chi, false, 11);
	}
      }
    }
  }
}


/* Pre-sintomático a Grave. Alto */
void reaction6(std::vector<grupo> &Val, std::vector<grupo> &Vba, Crandom &ran, trabajadores *altos, trabajadores *bajos, std::vector<lognormal_d> &dist){
  std::vector<int> aux;

  for(unsigned int i=0; i<Val[5].size(); i++){if(altos[Val[5][i]].tstate > TM){aux.push_back(Val[5][i]);}} //Presintomático
  for(unsigned int i=0; i<Val[6].size(); i++){if(altos[Val[6][i]].tstate > TM){aux.push_back(Val[6][i]);}} //Presintomático testeado
  for(unsigned int i=0; i<Val[7].size(); i++){if(altos[Val[7][i]].tstate > TM){aux.push_back(Val[7][i]);}} //Presintomático aislado
  for(unsigned int i=0; i<Val[8].size(); i++){if(altos[Val[8][i]].tstate > TM){aux.push_back(Val[8][i]);}} //Presintomático masivo

  if(aux.size() != 0){
    double value = (double)ran.r();
    unsigned int index = index_time(aux, altos, dist[4], value);
    unsigned int agent = aux[index];
    std::vector<int>::iterator it;
    unsigned int ind;
    if(Val[5].size() != 0){// Presintomático a Grave
      it = std::find(Val[5].begin(), Val[5].end(), agent);      ind = std::distance(Val[5].begin(), it);
      if(ind < Val[5].size()){
	mother_reaction(Val[5], Val[13], ind, altos, 5, 13);
	aux_main(1, Val, Vba, altos, bajos, ran, phi1, mu, true, 13);
      }
    }
    if(Val[6].size() != 0){// Presintomático testeado a Grave
      it = std::find(Val[6].begin(), Val[6].end(), agent);      ind = std::distance(Val[6].begin(), it);
      if(ind < Val[6].size()){
	mother_reaction(Val[6], Val[13], ind, altos, 6, 13);
	aux_main(1, Val, Vba, altos, bajos, ran, phi1, mu, true, 13);
      }
    }
    if(Val[7].size() != 0){// Presintomático aislado a Grave
      it = std::find(Val[7].begin(), Val[7].end(), agent);      ind = std::distance(Val[7].begin(), it);
      if(ind < Val[7].size()){mother_reaction(Val[7], Val[13], ind, altos, 7, 13);}
    }
    if(Val[8].size() != 0){// Presintomático masivo a Grave
      it = std::find(Val[8].begin(), Val[8].end(), agent);      ind = std::distance(Val[8].begin(), it);
      if(ind < Val[8].size()){
	mother_reaction(Val[8], Val[13], ind, altos, 8, 13);
	aux_main(1, Val, Vba, altos, bajos, ran, phi1, mu, true, 13);
      }
    }
  }
}


/* Pre-sintomático a Grave. Bajo */
void reaction7(std::vector<grupo> &Val, std::vector<grupo> &Vba, Crandom &ran, trabajadores *altos, trabajadores *bajos, std::vector<lognormal_d> &dist){
  std::vector<int> aux;

  for(unsigned int i=0; i<Vba[5].size(); i++){if(bajos[Vba[5][i]].tstate > TM){aux.push_back(Vba[5][i]);}} //Presintomático
  for(unsigned int i=0; i<Vba[6].size(); i++){if(bajos[Vba[6][i]].tstate > TM){aux.push_back(Vba[6][i]);}} //Presintomático testeado
  for(unsigned int i=0; i<Vba[7].size(); i++){if(bajos[Vba[7][i]].tstate > TM){aux.push_back(Vba[7][i]);}} //Presitnomático aislado
  for(unsigned int i=0; i<Vba[8].size(); i++){if(bajos[Vba[8][i]].tstate > TM){aux.push_back(Vba[8][i]);}} //Presitnomático masivo

  if(aux.size() != 0){
    double value = (double)ran.r();
    unsigned int index = index_time(aux, bajos, dist[5], value);
    unsigned int agent = aux[index];
    std::vector<int>::iterator it;
    unsigned int ind;
    if(Vba[5].size() != 0){// Presintomático a Grave
      it = std::find(Vba[5].begin(), Vba[5].end(), agent);      ind = std::distance(Vba[5].begin(), it);
      if(ind < Vba[5].size()){
	mother_reaction(Vba[5], Vba[13], ind, bajos, 5, 13);
	aux_main(1, Val, Vba, altos, bajos, ran, mu, chi, false, 13);
      }
    }
    if(Vba[6].size() != 0){// Presintomático testeado a Grave
      it = std::find(Vba[6].begin(), Vba[6].end(), agent);      ind = std::distance(Vba[6].begin(), it);
      if(ind < Vba[6].size()){
	mother_reaction(Vba[6], Vba[13], ind, bajos, 6, 13);
	aux_main(1, Val, Vba, altos, bajos, ran, mu, chi, false, 13);
      }
    }
    if(Vba[7].size() != 0){// Presintomático aislado a Grave
      it = std::find(Vba[7].begin(), Vba[7].end(), agent);      ind = std::distance(Vba[7].begin(), it);
      if(ind < Vba[7].size()){mother_reaction(Vba[7], Vba[13], ind, bajos, 7, 13);}
    }
    if(Vba[8].size() != 0){// Presintomático masivo a Grave
      it = std::find(Vba[8].begin(), Vba[8].end(), agent);      ind = std::distance(Vba[8].begin(), it);
      if(ind < Vba[8].size()){
	mother_reaction(Vba[8], Vba[13], ind, bajos, 8, 13);
	aux_main(1, Val, Vba, altos, bajos, ran, mu, chi, false, 13);
      }
    }
  }
}


/* Pre-sintomático a Recuperado. Alto */
void reaction8(std::vector<grupo> &Val, std::vector<grupo> &Vba, Crandom &ran, trabajadores *altos, trabajadores *bajos, std::vector<lognormal_d> &dist){
  std::vector<int> aux;

  for(unsigned int i=0; i<Val[5].size(); i++){if(altos[Val[5][i]].tstate > TM){aux.push_back(Val[5][i]);}} //Presintomático
  for(unsigned int i=0; i<Val[6].size(); i++){if(altos[Val[6][i]].tstate > TM){aux.push_back(Val[6][i]);}} //Presintomático testeado
  for(unsigned int i=0; i<Val[7].size(); i++){if(altos[Val[7][i]].tstate > TM){aux.push_back(Val[7][i]);}} //Presintomático aislado
  for(unsigned int i=0; i<Val[8].size(); i++){if(altos[Val[8][i]].tstate > TM){aux.push_back(Val[8][i]);}} //Presintomático masivo

  if(aux.size() != 0){
    double value = (double)ran.r();
    unsigned int index = index_time(aux, altos, dist[6], value);
    unsigned int agent = aux[index];
    std::vector<int>::iterator it;
    unsigned int ind;
    if(Val[5].size() != 0){// Presintomático a Recuperado no-detectado
      it = std::find(Val[5].begin(), Val[5].end(), agent);      ind = std::distance(Val[5].begin(), it);
      if(ind < Val[5].size()){mother_reaction(Val[5], Val[14], ind, altos, 5, 14);}
    }
    if(Val[6].size() != 0){// Presintomático testeado a Recuperado detectado o testeado
      it = std::find(Val[6].begin(), Val[6].end(), agent);      ind = std::distance(Val[6].begin(), it);
      if(ind < Val[6].size()){
        if(ran.r() < xi){// A Recuperado detectado
          mother_reaction(Val[6], Val[16], ind, altos, 6, 16);
 	  aux_main(1, Val, Vba, altos, bajos, ran, phi1, mu, true, 16);
        }
        else{// A Recuperado testeado
          mother_reaction(Val[6], Val[15], ind, altos, 6, 15);
        }
      }
    }
    if(Val[7].size() != 0){// Presintomático aislado a Recuperado detectado
      it = std::find(Val[7].begin(), Val[7].end(), agent);      ind = std::distance(Val[7].begin(), it);
      if(ind < Val[7].size()){mother_reaction(Val[7], Val[16], ind, altos, 7, 16);}
    }
    if(Val[8].size() != 0){// Presintomático masivo a Recuperado testeado
      it = std::find(Val[8].begin(), Val[8].end(), agent);      ind = std::distance(Val[8].begin(), it);
      if(ind < Val[8].size()){mother_reaction(Val[8], Val[15], ind, altos, 8, 15);}
    }
  }
}


/* Pre-sintomático a Recuperado. Bajo */
void reaction9(std::vector<grupo> &Val, std::vector<grupo> &Vba, Crandom &ran, trabajadores *altos, trabajadores *bajos, std::vector<lognormal_d> &dist){
  std::vector<int> aux;

  for(unsigned int i=0; i<Vba[5].size(); i++){if(bajos[Vba[5][i]].tstate > TM){aux.push_back(Vba[5][i]);}} //Presintomático
  for(unsigned int i=0; i<Vba[6].size(); i++){if(bajos[Vba[6][i]].tstate > TM){aux.push_back(Vba[6][i]);}} //Presintomático testeado
  for(unsigned int i=0; i<Vba[7].size(); i++){if(bajos[Vba[7][i]].tstate > TM){aux.push_back(Vba[7][i]);}} //Presintomático aislado
  for(unsigned int i=0; i<Vba[8].size(); i++){if(bajos[Vba[8][i]].tstate > TM){aux.push_back(Vba[8][i]);}} //Presintomático masivo

  if(aux.size() != 0){
    double value = (double)ran.r();
    unsigned int index = index_time(aux, bajos, dist[7], value);
    unsigned int agent = aux[index];
    std::vector<int>::iterator it;
    unsigned int ind;
    if(Vba[5].size() != 0){// Presintomático a Recuperado no-detectado
      it = std::find(Vba[5].begin(), Vba[5].end(), agent);      ind = std::distance(Vba[5].begin(), it);
      if(ind < Vba[5].size()){mother_reaction(Vba[5], Vba[14], ind, bajos, 5, 14);}
    }
    if(Vba[6].size() != 0){// Presintomático testeado a Recuperado detectado o testeado
      it = std::find(Vba[6].begin(), Vba[6].end(), agent);      ind = std::distance(Vba[6].begin(), it);
      if(ind < Vba[6].size()){
        if(ran.r() < xi){// A Recuperado detectado
          mother_reaction(Vba[6], Vba[16], ind, bajos, 6, 16);
 	  aux_main(1, Val, Vba, altos, bajos, ran, mu, chi, false, 16);
        }
        else{// A Recuperado testeado
          mother_reaction(Vba[6], Vba[15], ind, bajos, 6, 15);
        }
      }
    }
    if(Vba[7].size() != 0){// Presintomático aislado a Recuperado detectado
      it = std::find(Vba[7].begin(), Vba[7].end(), agent);      ind = std::distance(Vba[7].begin(), it);
      if(ind < Vba[7].size()){mother_reaction(Vba[7], Vba[16], ind, bajos, 7, 16);}
    }
    if(Vba[8].size() != 0){// Presintomático masivo a Recuperado testeado
      it = std::find(Vba[8].begin(), Vba[8].end(), agent);      ind = std::distance(Vba[8].begin(), it);
      if(ind < Vba[8].size()){mother_reaction(Vba[8], Vba[15], ind, bajos, 8, 15);}
    }
  }
}


/* Leve a Recuperado. Alto */
void reaction10(std::vector<grupo> &Val, std::vector<grupo> &Vba, Crandom &ran, trabajadores *altos, trabajadores *bajos, std::vector<lognormal_d> &dist){
  std::vector<int> aux;

  for(unsigned int i=0; i<Val[9].size(); i++){if(altos[Val[9][i]].tstate > TM){aux.push_back(Val[9][i]);}} //Leve
  for(unsigned int i=0; i<Val[10].size(); i++){if(altos[Val[10][i]].tstate > TM){aux.push_back(Val[10][i]);}} //Leve testeado
  for(unsigned int i=0; i<Val[11].size(); i++){if(altos[Val[11][i]].tstate > TM){aux.push_back(Val[11][i]);}} //Leve aislado
  for(unsigned int i=0; i<Val[12].size(); i++){if(altos[Val[12][i]].tstate > TM){aux.push_back(Val[12][i]);}} //Leve masivo

  if(aux.size() != 0){
    double value = (double)ran.r();
    unsigned int index = index_time(aux, altos, dist[8], value);
    unsigned int agent = aux[index];
    std::vector<int>::iterator it;
    unsigned int ind;
    if(Val[9].size() != 0){// Leve a Recuperado no-detectado
      it = std::find(Val[9].begin(), Val[9].end(), agent);      ind = std::distance(Val[9].begin(), it);
      if(ind < Val[9].size()){mother_reaction(Val[9], Val[14], ind, altos, 9, 14);}
    }
    if(Val[10].size() != 0){// Leve testeado a Recuperado detectado o testeado
      it = std::find(Val[10].begin(), Val[10].end(), agent);      ind = std::distance(Val[10].begin(), it);
      if(ind < Val[10].size()){
        if(ran.r() < xi){// A Recuperado detectado
          mother_reaction(Val[10], Val[16], ind, altos, 10, 16);
 	  aux_main(1, Val, Vba, altos, bajos, ran, phi1, mu, true, 16);
        }
        else{// A Recuperado testeado
          mother_reaction(Val[10], Val[16], ind, altos, 10, 16);
        }
      }
    }
    if(Val[11].size() != 0){// Leve aislado a Recuperado detectado
      it = std::find(Val[11].begin(), Val[11].end(), agent);      ind = std::distance(Val[11].begin(), it);
      if(ind < Val[11].size()){mother_reaction(Val[11], Val[16], ind, altos, 11, 16);}
    }
    if(Val[12].size() != 0){// Leve masivo a Recuperado testeado
      it = std::find(Val[12].begin(), Val[12].end(), agent);      ind = std::distance(Val[12].begin(), it);
      if(ind < Val[12].size()){mother_reaction(Val[12], Val[15], ind, altos, 12, 15);}
    }
  }
}


/* Leve a Recuperado. Bajo */
void reaction11(std::vector<grupo> &Val, std::vector<grupo> &Vba, Crandom &ran, trabajadores *altos, trabajadores *bajos, std::vector<lognormal_d> &dist){
  std::vector<int> aux;

  for(unsigned int i=0; i<Vba[9].size(); i++){if(bajos[Vba[9][i]].tstate > TM){aux.push_back(Vba[9][i]);}} //Leve
  for(unsigned int i=0; i<Vba[10].size(); i++){if(bajos[Vba[10][i]].tstate > TM){aux.push_back(Vba[10][i]);}} //Leve testeado
  for(unsigned int i=0; i<Vba[11].size(); i++){if(bajos[Vba[11][i]].tstate > TM){aux.push_back(Vba[11][i]);}} //Leve aislado
  for(unsigned int i=0; i<Vba[12].size(); i++){if(bajos[Vba[12][i]].tstate > TM){aux.push_back(Vba[12][i]);}} //Leve masivo

  if(aux.size() != 0){
    double value = (double)ran.r();
    unsigned int index = index_time(aux, bajos, dist[9], value);
    unsigned int agent = aux[index];
    std::vector<int>::iterator it;
    unsigned int ind;
    if(Vba[9].size() != 0){// Leve a Recuperado no-detectado
      it = std::find(Vba[9].begin(), Vba[9].end(), agent);      ind = std::distance(Vba[9].begin(), it);
      if(ind < Vba[9].size()){mother_reaction(Vba[9], Vba[14], ind, bajos, 9, 14);}
    }
    if(Vba[10].size() != 0){// Leve testeado a Recuperado detectado o testeado
      it = std::find(Vba[10].begin(), Vba[10].end(), agent);      ind = std::distance(Vba[10].begin(), it);
      if(ind < Vba[10].size()){
        if(ran.r() < xi){// A Recuperado detectado
          mother_reaction(Vba[10], Vba[16], ind, bajos, 10, 16);
 	  aux_main(1, Val, Vba, altos, bajos, ran, mu, chi, false, 16);
        }
        else{// A Recuperado testeado
          mother_reaction(Vba[10], Vba[15], ind, bajos, 10, 15);
        }
      }
    }
    if(Vba[11].size() != 0){// Leve aislado a Recuperado detectado
      it = std::find(Vba[11].begin(), Vba[11].end(), agent);      ind = std::distance(Vba[11].begin(), it);
      if(ind < Vba[11].size()){mother_reaction(Vba[11], Vba[16], ind, bajos, 11, 16);}
    }
    if(Vba[12].size() != 0){// Leve masivo a Recuperado testeado
      it = std::find(Vba[12].begin(), Vba[12].end(), agent);      ind = std::distance(Vba[12].begin(), it);
      if(ind < Vba[12].size()){mother_reaction(Vba[12], Vba[15], ind, bajos, 12, 15);}
    }
  }
}


/* Grave a Recuperado. Alto */
void reaction12(std::vector<grupo> &Val, std::vector<grupo> &Vba, Crandom &ran, trabajadores *altos, trabajadores *bajos, std::vector<lognormal_d> &dist){
  std::vector<int> aux;

  for(unsigned int i=0; i<Val[13].size(); i++){if(altos[Val[13][i]].tstate > TM){aux.push_back(Val[13][i]);}} //Grave

  if(aux.size() != 0){// Grave a Recuperado detectado
    unsigned int index = index_time(aux, altos, dist[10], (double)ran.r());
    unsigned int agent = aux[index];
    std::vector<int>::iterator it = std::find(Val[13].begin(), Val[13].end(), agent);
    unsigned int ind = std::distance(Val[13].begin(), it);
    if(ind < Val[13].size()){mother_reaction(Val[13], Val[16], ind, altos, 13, 16);}
  }
}


/* Grave a Recuperado. Bajo */
void reaction13(std::vector<grupo> &Val, std::vector<grupo> &Vba, Crandom &ran, trabajadores *altos, trabajadores *bajos, std::vector<lognormal_d> &dist){
  std::vector<int> aux;

  for(unsigned int i=0; i<Vba[13].size(); i++){if(bajos[Vba[13][i]].tstate > TM){aux.push_back(Vba[13][i]);}} //Grave

  if(aux.size() != 0){// Grave a Recuperado detectado
    unsigned int index = index_time(aux, bajos, dist[11], (double)ran.r());
    unsigned int agent = aux[index];
    std::vector<int>::iterator it = std::find(Vba[13].begin(), Vba[13].end(), agent);
    unsigned int ind = std::distance(Vba[13].begin(), it);
    if(ind < Vba[13].size()){mother_reaction(Vba[13], Vba[16], ind, bajos, 13, 16);}
  }
}
