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


void reaction0(grupo &Sa, grupo &Sb, grupo &STa, grupo &STb, grupo &Ea, grupo &Eb, grupo &ETa, grupo &ETb, grupo &Pa, grupo &Pb, grupo &PTa, grupo &PTb, grupo &PTAa, grupo &PTAb, grupo &La, grupo &Lb, grupo &LTa, grupo &LTb, grupo &LTAa, grupo &LTAb, grupo &IAa, grupo &IAb, grupo &RTa, grupo&RTb, grupo &RIa, grupo &RIb, grupo &RAa, grupo &RAb, Crandom &ran, trabajadores *altos, trabajadores *bajos, std::vector<lognormal_d> &dist){
  unsigned int index = (int)(ran.r()*(Sa.size() + STa.size()));
  int agentS, value1 = Ea.size(), value2;
  if(index < Sa.size()){agentS = Sa[index];    mother_reaction(Sa, Ea, index, altos, 0, 2);}
  else{agentS = STa[index-Sa.size()];    mother_reaction(STa, ETa, index-Sa.size(), altos, 1, 3);}
  value2 = Ea.size();

  int agentI = who_infected(Pa, Pb, PTa, PTb, PTAa, PTAb, La, Lb, LTa, LTb, LTAa, LTAb, IAa, IAb, phi1, mu, ran, agentS, altos, bajos);
  if(value2 > value1){altos[Ea.back()].DF = agentI;}
  else{altos[ETa.back()].DF = agentI;}
}


void reaction1(grupo &Sa, grupo &Sb, grupo &STa, grupo &STb, grupo &Ea, grupo &Eb, grupo &ETa, grupo &ETb, grupo &Pa, grupo &Pb, grupo &PTa, grupo &PTb, grupo &PTAa, grupo &PTAb, grupo &La, grupo &Lb, grupo &LTa, grupo &LTb, grupo &LTAa, grupo &LTAb, grupo &IAa, grupo &IAb, grupo &RTa, grupo&RTb, grupo &RIa, grupo &RIb, grupo &RAa, grupo &RAb, Crandom &ran, trabajadores *altos, trabajadores *bajos, std::vector<lognormal_d> &dist){
  unsigned int index = (int)(ran.r()*(Sb.size() + STb.size()));
  int agentS, value1 = Eb.size(), value2;
  if(index < Sb.size()){agentS = Sb[index];    mother_reaction(Sb, Eb, index, bajos, 0, 2);}
  else{agentS = STb[index-Sb.size()];    mother_reaction(STb, ETb, index-Sb.size(), bajos, 1, 3);}
  value2 = Eb.size();

  int agentI = who_infected(Pa, Pb, PTa, PTb, PTAa, PTAb, La, Lb, LTa, LTb, LTAa, LTAb, IAa, IAb, mu, chi, ran, agentS + Na, altos, bajos);
  if(value2 > value1){bajos[Eb.back()].DF = agentI;}
  else{bajos[ETb.back()].DF = agentI;}
}


void reaction2(grupo &Sa, grupo &Sb, grupo &STa, grupo &STb, grupo &Ea, grupo &Eb, grupo &ETa, grupo &ETb, grupo &Pa, grupo &Pb, grupo &PTa, grupo &PTb, grupo &PTAa, grupo &PTAb, grupo &La, grupo &Lb, grupo &LTa, grupo &LTb, grupo &LTAa, grupo &LTAb, grupo &IAa, grupo &IAb, grupo &RTa, grupo&RTb, grupo &RIa, grupo &RIb, grupo &RAa, grupo &RAb, Crandom &ran, trabajadores *altos, trabajadores *bajos, std::vector<lognormal_d> &dist){
  std::vector<int> aux;

  for(unsigned int i=0; i<Ea.size(); i++){if(altos[Ea[i]].tstate > TM){aux.push_back(Ea[i]);}}
  for(unsigned int i=0; i<ETa.size(); i++){if(altos[ETa[i]].tstate > TM){aux.push_back(ETa[i]);}}

  if(aux.size() != 0){
    double value = (double)ran.r();
    unsigned int index = index_time(aux, altos, dist[0], value);
    unsigned int agent = aux[index];
    std::vector<int>::iterator it;
    unsigned int ind;
    if(Ea.size() != 0){
      it = std::find(Ea.begin(), Ea.end(), agent);      ind = std::distance(Ea.begin(), it);
      if(ind < Ea.size()){mother_reaction(Ea, Pa, ind, altos, 2, 4);}
    }
    if(ETa.size() != 0){
      it = std::find(ETa.begin(), ETa.end(), agent);      ind = std::distance(ETa.begin(), it);
      if(ind < ETa.size()){mother_reaction(ETa, Pa, ind, altos, 3, 4);}
    }
  }
}


void reaction3(grupo &Sa, grupo &Sb, grupo &STa, grupo &STb, grupo &Ea, grupo &Eb, grupo &ETa, grupo &ETb, grupo &Pa, grupo &Pb, grupo &PTa, grupo &PTb, grupo &PTAa, grupo &PTAb, grupo &La, grupo &Lb, grupo &LTa, grupo &LTb, grupo &LTAa, grupo &LTAb, grupo &IAa, grupo &IAb, grupo &RTa, grupo&RTb, grupo &RIa, grupo &RIb, grupo &RAa, grupo &RAb, Crandom &ran, trabajadores *altos, trabajadores *bajos, std::vector<lognormal_d> &dist){
  std::vector<int> aux;

  for(unsigned int i=0; i<Eb.size(); i++){if(bajos[Eb[i]].tstate > TM){aux.push_back(Eb[i]);}}
  for(unsigned int i=0; i<ETb.size(); i++){if(bajos[ETb[i]].tstate > TM){aux.push_back(ETb[i]);}}  

  if(aux.size() != 0){
    double value = (double)ran.r();
    unsigned int index = index_time(aux, bajos, dist[1], value);
    unsigned int agent = aux[index];
    std::vector<int>::iterator it;
    unsigned int ind;
    if(Eb.size() != 0){
      it = std::find(Eb.begin(), Eb.end(), agent);      ind = std::distance(Eb.begin(), it);
      if(ind < Eb.size()){mother_reaction(Eb, Pb, ind, bajos, 2, 4);}
    }
    if(ETb.size() != 0){
      it = std::find(ETb.begin(), ETb.end(), agent);      ind = std::distance(ETb.begin(), it);
      if(ind < ETb.size()){mother_reaction(ETb, Pb, ind, bajos, 3, 4);}
    }    
  }
}


void reaction4(grupo &Sa, grupo &Sb, grupo &STa, grupo &STb, grupo &Ea, grupo &Eb, grupo &ETa, grupo &ETb, grupo &Pa, grupo &Pb, grupo &PTa, grupo &PTb, grupo &PTAa, grupo &PTAb, grupo &La, grupo &Lb, grupo &LTa, grupo &LTb, grupo &LTAa, grupo &LTAb, grupo &IAa, grupo &IAb, grupo &RTa, grupo&RTb, grupo &RIa, grupo &RIb, grupo &RAa, grupo &RAb, Crandom &ran, trabajadores *altos, trabajadores *bajos, std::vector<lognormal_d> &dist){
  std::vector<int> aux;

  for(unsigned int i=0; i<Pa.size(); i++){if(altos[Pa[i]].tstate > TM){aux.push_back(Pa[i]);}}
  for(unsigned int i=0; i<PTa.size(); i++){if(altos[PTa[i]].tstate > TM){aux.push_back(PTa[i]);}}
  for(unsigned int i=0; i<PTAa.size(); i++){if(altos[PTAa[i]].tstate > TM){aux.push_back(PTAa[i]);}}

  if(aux.size() != 0){
    double value = (double)ran.r();
    unsigned int index = index_time(aux, altos, dist[2], value);
    unsigned int agent = aux[index];
    std::vector<int>::iterator it;
    unsigned int ind;
    if(Pa.size() != 0){
      it = std::find(Pa.begin(), Pa.end(), agent);      ind = std::distance(Pa.begin(), it);
      if(ind < Pa.size()){
	mother_reaction(Pa, La, ind, altos, 4, 7);
	if(ran.r() < iota){
	  tested_reaction(La, LTAa, La.size()-1, altos, 7, 9, 0.0);
	  aux_main(1, LTAa, Sa, Sb, STa, STb, Ea, Eb, ETa, ETb, Pa, Pb, PTa, PTb, La, Lb, LTa, LTb, RIa, RIb, RTa, RTb, altos, bajos, ran, phi1, mu, true);
	}
      }
    }
    if(PTa.size() != 0){
      it = std::find(PTa.begin(), PTa.end(), agent);      ind = std::distance(PTa.begin(), it);
      if(ind < PTa.size()){mother_reaction(PTa, LTa, ind, altos, 5, 8);}
    }
    if(PTAa.size() != 0){
      it = std::find(PTAa.begin(), PTAa.end(), agent);      ind = std::distance(PTAa.begin(), it);
      if(ind < PTAa.size()){mother_reaction(PTAa, LTAa, ind, altos, 6, 9);}
    }
  }
}


void reaction5(grupo &Sa, grupo &Sb, grupo &STa, grupo &STb, grupo &Ea, grupo &Eb, grupo &ETa, grupo &ETb, grupo &Pa, grupo &Pb, grupo &PTa, grupo &PTb, grupo &PTAa, grupo &PTAb, grupo &La, grupo &Lb, grupo &LTa, grupo &LTb, grupo &LTAa, grupo &LTAb, grupo &IAa, grupo &IAb, grupo &RTa, grupo&RTb, grupo &RIa, grupo &RIb, grupo &RAa, grupo &RAb, Crandom &ran, trabajadores *altos, trabajadores *bajos, std::vector<lognormal_d> &dist){
  std::vector<int> aux;

  for(unsigned int i=0; i<Pb.size(); i++){if(bajos[Pb[i]].tstate > TM){aux.push_back(Pb[i]);}}
  for(unsigned int i=0; i<PTb.size(); i++){if(bajos[PTb[i]].tstate > TM){aux.push_back(PTb[i]);}}
  for(unsigned int i=0; i<PTAb.size(); i++){if(bajos[PTAb[i]].tstate > TM){aux.push_back(PTAb[i]);}}

  if(aux.size() != 0){
    double value = (double)ran.r();
    unsigned int index = index_time(aux, bajos, dist[3], value);
    unsigned int agent = aux[index];
    std::vector<int>::iterator it;
    unsigned int ind;
    if(Pb.size() != 0){
      it = std::find(Pb.begin(), Pb.end(), agent);      ind = std::distance(Pb.begin(), it);
      if(ind < Pb.size()){
	mother_reaction(Pb, Lb, ind, bajos, 4, 7);
	if(ran.r() < iota){
	  tested_reaction(Lb, LTAb, Lb.size()-1, bajos, 7, 9, 0.0);
	  aux_main(1, LTAb, Sa, Sb, STa, STb, Ea, Eb, ETa, ETb, Pa, Pb, PTa, PTb, La, Lb, LTa, LTb, RIa, RIb, RTa, RTb, altos, bajos, ran, mu, chi, false);
	}
      }
    }
    if(PTb.size() != 0){
      it = std::find(PTb.begin(), PTb.end(), agent);      ind = std::distance(PTb.begin(), it);
      if(ind < PTb.size()){mother_reaction(PTb, LTb, ind, bajos, 5, 8);}
    }
    if(PTAb.size() != 0){
      it = std::find(PTAb.begin(), PTAb.end(), agent);      ind = std::distance(PTAb.begin(), it);
      if(ind < PTAb.size()){mother_reaction(PTAb, LTAb, ind, bajos, 6, 9);}
    }
  }
}


void reaction6(grupo &Sa, grupo &Sb, grupo &STa, grupo &STb, grupo &Ea, grupo &Eb, grupo &ETa, grupo &ETb, grupo &Pa, grupo &Pb, grupo &PTa, grupo &PTb, grupo &PTAa, grupo &PTAb, grupo &La, grupo &Lb, grupo &LTa, grupo &LTb, grupo &LTAa, grupo &LTAb, grupo &IAa, grupo &IAb, grupo &RTa, grupo&RTb, grupo &RIa, grupo &RIb, grupo &RAa, grupo &RAb, Crandom &ran, trabajadores *altos, trabajadores *bajos, std::vector<lognormal_d> &dist){
  std::vector<int> aux;

  for(unsigned int i=0; i<Pa.size(); i++){if(altos[Pa[i]].tstate > TM){aux.push_back(Pa[i]);}}
  for(unsigned int i=0; i<PTa.size(); i++){if(altos[PTa[i]].tstate > TM){aux.push_back(PTa[i]);}}
  for(unsigned int i=0; i<PTAa.size(); i++){if(altos[PTAa[i]].tstate > TM){aux.push_back(PTAa[i]);}}

  if(aux.size() != 0){
    double value = (double)ran.r();
    unsigned int index = index_time(aux, altos, dist[4], value);
    unsigned int agent = aux[index];
    std::vector<int>::iterator it;
    unsigned int ind;
    if(Pa.size() != 0){
      it = std::find(Pa.begin(), Pa.end(), agent);      ind = std::distance(Pa.begin(), it);
      if(ind < Pa.size()){mother_reaction(Pa, IAa, ind, altos, 4, 10);}
    }
    if(PTa.size() != 0){
      it = std::find(PTa.begin(), PTa.end(), agent);      ind = std::distance(PTa.begin(), it);
      if(ind < PTa.size()){mother_reaction(PTa, IAa, ind, altos, 5, 10);}
    }
    if(PTAa.size() != 0){
      it = std::find(PTAa.begin(), PTAa.end(), agent);      ind = std::distance(PTAa.begin(), it);
      if(ind < PTAa.size()){mother_reaction(PTAa, IAa, ind, altos, 6, 10);}
    }
    aux_main(1, IAa, Sa, Sb, STa, STb, Ea, Eb, ETa, ETb, Pa, Pb, PTa, PTb, La, Lb, LTa, LTb, RIa, RIb, RTa, RTb, altos, bajos, ran, phi1, mu, true);
  }
}


void reaction7(grupo &Sa, grupo &Sb, grupo &STa, grupo &STb, grupo &Ea, grupo &Eb, grupo &ETa, grupo &ETb, grupo &Pa, grupo &Pb, grupo &PTa, grupo &PTb, grupo &PTAa, grupo &PTAb, grupo &La, grupo &Lb, grupo &LTa, grupo &LTb, grupo &LTAa, grupo &LTAb, grupo &IAa, grupo &IAb, grupo &RTa, grupo&RTb, grupo &RIa, grupo &RIb, grupo &RAa, grupo &RAb, Crandom &ran, trabajadores *altos, trabajadores *bajos, std::vector<lognormal_d> &dist){
  std::vector<int> aux;

  for(unsigned int i=0; i<Pb.size(); i++){if(bajos[Pb[i]].tstate > TM){aux.push_back(Pb[i]);}}
  for(unsigned int i=0; i<PTb.size(); i++){if(bajos[PTb[i]].tstate > TM){aux.push_back(PTb[i]);}}
  for(unsigned int i=0; i<PTAb.size(); i++){if(bajos[PTAb[i]].tstate > TM){aux.push_back(PTAb[i]);}}

  if(aux.size() != 0){
    double value = (double)ran.r();
    unsigned int index = index_time(aux, bajos, dist[5], value);
    unsigned int agent = aux[index];
    std::vector<int>::iterator it;
    unsigned int ind;
    if(Pb.size() != 0){
      it = std::find(Pb.begin(), Pb.end(), agent);      ind = std::distance(Pb.begin(), it);
      if(ind < Pb.size()){mother_reaction(Pb, IAb, ind, bajos, 4, 10);}
    }
    if(PTb.size() != 0){
      it = std::find(PTb.begin(), PTb.end(), agent);      ind = std::distance(PTb.begin(), it);
      if(ind < PTb.size()){mother_reaction(PTb, IAb, ind, bajos, 5, 10);}
    }
    if(PTAb.size() != 0){
      it = std::find(PTAb.begin(), PTAb.end(), agent);      ind = std::distance(PTAb.begin(), it);
      if(ind < PTAb.size()){mother_reaction(PTAb, IAb, ind, bajos, 6, 10);}
    }
    aux_main(1, IAb, Sa, Sb, STa, STb, Ea, Eb, ETa, ETb, Pa, Pb, PTa, PTb, La, Lb, LTa, LTb, RIa, RIb, RTa, RTb, altos, bajos, ran, mu, chi, false);
  }
}


void reaction8(grupo &Sa, grupo &Sb, grupo &STa, grupo &STb, grupo &Ea, grupo &Eb, grupo &ETa, grupo &ETb, grupo &Pa, grupo &Pb, grupo &PTa, grupo &PTb, grupo &PTAa, grupo &PTAb, grupo &La, grupo &Lb, grupo &LTa, grupo &LTb, grupo &LTAa, grupo &LTAb, grupo &IAa, grupo &IAb, grupo &RTa, grupo&RTb, grupo &RIa, grupo &RIb, grupo &RAa, grupo &RAb, Crandom &ran, trabajadores *altos, trabajadores *bajos, std::vector<lognormal_d> &dist){
  std::vector<int> aux;

  for(unsigned int i=0; i<Pa.size(); i++){if(altos[Pa[i]].tstate > TM){aux.push_back(Pa[i]);}}
  for(unsigned int i=0; i<PTa.size(); i++){if(altos[PTa[i]].tstate > TM){aux.push_back(PTa[i]);}}
  for(unsigned int i=0; i<PTAa.size(); i++){if(altos[PTAa[i]].tstate > TM){aux.push_back(PTAa[i]);}}

  if(aux.size() != 0){
    double value = (double)ran.r();
    unsigned int index = index_time(aux, altos, dist[6], value);
    unsigned int agent = aux[index];
    std::vector<int>::iterator it;
    unsigned int ind;
    if(Pa.size() != 0){
      it = std::find(Pa.begin(), Pa.end(), agent);      ind = std::distance(Pa.begin(), it);
      if(ind < Pa.size()){mother_reaction(Pa, RIa, ind, altos, 4, 11);}
    }
    if(PTa.size() != 0){
      it = std::find(PTa.begin(), PTa.end(), agent);      ind = std::distance(PTa.begin(), it);
      if(ind < PTa.size()){mother_reaction(PTa, RAa, ind, altos, 5, 13);}
    }
    if(PTAa.size() != 0){
      it = std::find(PTAa.begin(), PTAa.end(), agent);      ind = std::distance(PTAa.begin(), it);
      if(ind < PTAa.size()){mother_reaction(PTAa, RAa, ind, altos, 6, 13);}
    }
  }
}


void reaction9(grupo &Sa, grupo &Sb, grupo &STa, grupo &STb, grupo &Ea, grupo &Eb, grupo &ETa, grupo &ETb, grupo &Pa, grupo &Pb, grupo &PTa, grupo &PTb, grupo &PTAa, grupo &PTAb, grupo &La, grupo &Lb, grupo &LTa, grupo &LTb, grupo &LTAa, grupo &LTAb, grupo &IAa, grupo &IAb, grupo &RTa, grupo&RTb, grupo &RIa, grupo &RIb, grupo &RAa, grupo &RAb, Crandom &ran, trabajadores *altos, trabajadores *bajos, std::vector<lognormal_d> &dist){
  std::vector<int> aux;

  for(unsigned int i=0; i<Pb.size(); i++){if(bajos[Pb[i]].tstate > TM){aux.push_back(Pb[i]);}}
  for(unsigned int i=0; i<PTb.size(); i++){if(bajos[PTb[i]].tstate > TM){aux.push_back(PTb[i]);}}
  for(unsigned int i=0; i<PTAb.size(); i++){if(bajos[PTAb[i]].tstate > TM){aux.push_back(PTAb[i]);}}

  if(aux.size() != 0){
    double value = (double)ran.r();
    unsigned int index = index_time(aux, bajos, dist[7], value);
    unsigned int agent = aux[index];
    std::vector<int>::iterator it;
    unsigned int ind;
    if(Pb.size() != 0){
      it = std::find(Pb.begin(), Pb.end(), agent);      ind = std::distance(Pb.begin(), it);
      if(ind < Pb.size()){mother_reaction(Pb, RIb, ind, bajos, 4, 11);}
    }
    if(PTb.size() != 0){
      it = std::find(PTb.begin(), PTb.end(), agent);      ind = std::distance(PTb.begin(), it);
      if(ind < PTb.size()){mother_reaction(PTb, RAb, ind, bajos, 5, 13);}
    }
    if(PTAb.size() != 0){
      it = std::find(PTAb.begin(), PTAb.end(), agent);      ind = std::distance(PTAb.begin(), it);
      if(ind < PTAb.size()){mother_reaction(PTAb, RAb, ind, bajos, 6, 13);}
    }
  }
}


void reaction10(grupo &Sa, grupo &Sb, grupo &STa, grupo &STb, grupo &Ea, grupo &Eb, grupo &ETa, grupo &ETb, grupo &Pa, grupo &Pb, grupo &PTa, grupo &PTb, grupo &PTAa, grupo &PTAb, grupo &La, grupo &Lb, grupo &LTa, grupo &LTb, grupo &LTAa, grupo &LTAb, grupo &IAa, grupo &IAb, grupo &RTa, grupo&RTb, grupo &RIa, grupo &RIb, grupo &RAa, grupo &RAb, Crandom &ran, trabajadores *altos, trabajadores *bajos, std::vector<lognormal_d> &dist){
  std::vector<int> aux;

  for(unsigned int i=0; i<La.size(); i++){if(altos[La[i]].tstate > TM){aux.push_back(La[i]);}}
  for(unsigned int i=0; i<LTa.size(); i++){if(altos[LTa[i]].tstate > TM){aux.push_back(LTa[i]);}}
  for(unsigned int i=0; i<LTAa.size(); i++){if(altos[LTAa[i]].tstate > TM){aux.push_back(LTAa[i]);}}

  if(aux.size() != 0){
    double value = (double)ran.r();
    unsigned int index = index_time(aux, altos, dist[8], value);
    unsigned int agent = aux[index];
    std::vector<int>::iterator it;
    unsigned int ind;
    if(La.size() != 0){
      it = std::find(La.begin(), La.end(), agent);      ind = std::distance(La.begin(), it);
      if(ind < La.size()){mother_reaction(La, RIa, ind, altos, 7, 11);}
    }
    if(LTa.size() != 0){
      it = std::find(LTa.begin(), LTa.end(), agent);      ind = std::distance(LTa.begin(), it);
      if(ind < LTa.size()){mother_reaction(LTa, RAa, ind, altos, 8, 13);}
    }
    if(LTAa.size() != 0){
      it = std::find(LTAa.begin(), LTAa.end(), agent);      ind = std::distance(LTAa.begin(), it);
      if(ind < LTAa.size()){mother_reaction(LTAa, RAa, ind, altos, 9, 13);}
    }
  }
}


void reaction11(grupo &Sa, grupo &Sb, grupo &STa, grupo &STb, grupo &Ea, grupo &Eb, grupo &ETa, grupo &ETb, grupo &Pa, grupo &Pb, grupo &PTa, grupo &PTb, grupo &PTAa, grupo &PTAb, grupo &La, grupo &Lb, grupo &LTa, grupo &LTb, grupo &LTAa, grupo &LTAb, grupo &IAa, grupo &IAb, grupo &RTa, grupo&RTb, grupo &RIa, grupo &RIb, grupo &RAa, grupo &RAb, Crandom &ran, trabajadores *altos, trabajadores *bajos, std::vector<lognormal_d> &dist){
  std::vector<int> aux;

  for(unsigned int i=0; i<Lb.size(); i++){if(bajos[Lb[i]].tstate > TM){aux.push_back(Lb[i]);}}
  for(unsigned int i=0; i<LTb.size(); i++){if(bajos[LTb[i]].tstate > TM){aux.push_back(LTb[i]);}}
  for(unsigned int i=0; i<LTAb.size(); i++){if(bajos[LTAb[i]].tstate > TM){aux.push_back(LTAb[i]);}}

  if(aux.size() != 0){
    double value = (double)ran.r();
    unsigned int index = index_time(aux, bajos, dist[9], value);
    unsigned int agent = aux[index];
    std::vector<int>::iterator it;
    unsigned int ind;
    if(Lb.size() != 0){
      it = std::find(Lb.begin(), Lb.end(), agent);      ind = std::distance(Lb.begin(), it);
      if(ind < Lb.size()){mother_reaction(Lb, RIb, ind, bajos, 7, 11);}
    }
    if(LTb.size() != 0){
      it = std::find(LTb.begin(), LTb.end(), agent);      ind = std::distance(LTb.begin(), it);
      if(ind < LTb.size()){mother_reaction(LTb, RAb, ind, bajos, 8, 13);}
    }
    if(LTAb.size() != 0){
      it = std::find(LTAb.begin(), LTAb.end(), agent);      ind = std::distance(LTAb.begin(), it);
      if(ind < LTAb.size()){mother_reaction(LTAb, RAb, ind, bajos, 9, 13);}
    }
  }
}


void reaction12(grupo &Sa, grupo &Sb, grupo &STa, grupo &STb, grupo &Ea, grupo &Eb, grupo &ETa, grupo &ETb, grupo &Pa, grupo &Pb, grupo &PTa, grupo &PTb, grupo &PTAa, grupo &PTAb, grupo &La, grupo &Lb, grupo &LTa, grupo &LTb, grupo &LTAa, grupo &LTAb, grupo &IAa, grupo &IAb, grupo &RTa, grupo&RTb, grupo &RIa, grupo &RIb, grupo &RAa, grupo &RAb, Crandom &ran, trabajadores *altos, trabajadores *bajos, std::vector<lognormal_d> &dist){
  std::vector<int> aux;
  
  for(unsigned int i=0; i<IAa.size(); i++){if(altos[IAa[i]].tstate > TM){aux.push_back(IAa[i]);}}

  if(aux.size() != 0){    
    unsigned int index = index_time(aux, altos, dist[10], (double)ran.r());
    unsigned int agent = aux[index];
    std::vector<int>::iterator it = std::find(IAa.begin(), IAa.end(), agent);
    unsigned int ind = std::distance(IAa.begin(), it);
    if(ind < IAa.size()){mother_reaction(IAa, RAa, ind, altos, 10, 13);}
  }
}


void reaction13(grupo &Sa, grupo &Sb, grupo &STa, grupo &STb, grupo &Ea, grupo &Eb, grupo &ETa, grupo &ETb, grupo &Pa, grupo &Pb, grupo &PTa, grupo &PTb, grupo &PTAa, grupo &PTAb, grupo &La, grupo &Lb, grupo &LTa, grupo &LTb, grupo &LTAa, grupo &LTAb, grupo &IAa, grupo &IAb, grupo &RTa, grupo&RTb, grupo &RIa, grupo &RIb, grupo &RAa, grupo &RAb, Crandom &ran, trabajadores *altos, trabajadores *bajos, std::vector<lognormal_d> &dist){
  std::vector<int> aux;
  
  for(unsigned int i=0; i<IAb.size(); i++){if(bajos[IAb[i]].tstate > TM){aux.push_back(IAb[i]);}}

  if(aux.size() != 0){
    unsigned int index = index_time(aux, bajos, dist[11], (double)ran.r());
    unsigned int agent = aux[index];
    std::vector<int>::iterator it = std::find(IAb.begin(), IAb.end(), agent);
    unsigned int ind = std::distance(IAb.begin(), it);
    if(ind < IAb.size()){mother_reaction(IAb, RAb, ind, bajos, 10, 13);}
  }  
}
