#include <reaction.h>
#include <trace.h>
#include <test.h>

void mother_reaction(group &Out, group &In, int index, Workers *family, Stages typeout, Stages typein){
  int agent = Out[index];
  Out.erase(Out.begin() + index);
  In.push_back(agent);
  family[agent].change(typein, typeout);
  family[agent].tstate = 0.0;
}


int index_time(group &Out, Workers *family, lognormal_d &dist, double value){
  //Jose
  unsigned int n = Out.size(), agent;
  double times[n];
  double param = value;

  //Hallo la diferencia del tiempo con el tiempo que lleve en el estado
  for(unsigned int i=0; i<n; i++){agent = Out[i];    times[i] = std::abs(cdf(dist, family[agent].tstate) - param);}

  //Retorno el índice donde está el tiempo mínimo
  return std::distance(times, std::min_element(times, times+n));
}


int infected(Workers *family, int max){
  Stages kind;
  for(int i=0; i<max; i++){
    kind = family[i].kind;
    if(kind == SUS || kind == SUST || kind == SUSA){
      return i;
    }
  }
  return -1;
}



/* Susceptible a expuesto. Alto */
void reaction0(std::map<std::string, group> &Val, std::map<std::string, group> &Vba, Crandom &ran, Workers *altos, Workers *bajos, std::vector<lognormal_d> &dist, int agentI){
  int agentS;
  unsigned int index;
  double auxind;
  if(agentI > -1){/* Si el que infecta no es externo */
    auxind = ran.r()*(Val["SUS"].size() + Val["SUST"].size() + (1-MyCons.alpha)*Val["SUSA"].size());

    /* Susceptible a Expuesto */
    if(auxind < Val["SUS"].size()){
      index = (int)(ran.r()*Val["SUS"].size());      agentS = Val["SUS"][index];
      mother_reaction(Val["SUS"], Val["EXP"], index, altos, SUS, EXP);
      altos[Val["EXP"].back()].DF = agentI;
    }
    /* Susceptible testeado a Expuesto testado */
    else if(auxind < Val["SUS"].size() + Val["SUST"].size()){
      index = (int)(ran.r()*Val["SUST"].size());      agentS = Val["SUST"][index];
      mother_reaction(Val["SUST"], Val["EXPT"], index, altos, SUST, EXPT);
      altos[Val["EXPT"].back()].DF = agentI;
    }
    /* Susceptible aislado a Expuesto aislado */
    else{
      index = (int)(ran.r()*Val["SUSA"].size());      agentS = Val["SUSA"][index];
      mother_reaction(Val["SUSA"], Val["EXPA"], index, altos, SUSA, EXPA);
      altos[Val["EXPA"].back()].DF = agentI;
    }

    /* Guardo a la persona que infecté */
    if(agentI < Na){altos[agentI].my_inf.push_back(agentS);}
    else{bajos[agentI-Na].my_inf.push_back(agentS);}
  }
  else{/* Si el que infecta es externo */
    agentS = infected(altos, Na);
    if(agentS != -1){/* Si la infección debido al externo se realizó */
      /* Susceptible a Expuesto */
      if(altos[agentS].kind == SUS){
	group::iterator it = std::find(Val["SUS"].begin(), Val["SUS"].end(), agentS);
	index = std::distance(Val["SUS"].begin(), it);
	mother_reaction(Val["SUS"], Val["EXP"], index, altos, SUS, EXP);
      }

      /* Susceptible testeado a Expuesto testado */
      else if(altos[agentS].kind == SUST){
	group::iterator it = std::find(Val["SUST"].begin(), Val["SUST"].end(), agentS);
	index = std::distance(Val["SUST"].begin(), it);
	mother_reaction(Val["SUST"], Val["EXPT"], index, altos, SUST, EXPT);
      }

      /* Susceptible aislado a Expuesto aislado */
      else{
	group::iterator it = std::find(Val["SUSA"].begin(), Val["SUSA"].end(), agentS);
	index = std::distance(Val["SUSA"].begin(), it);
	mother_reaction(Val["SUSA"], Val["EXPA"], index, altos, SUSA, EXPA);
      }
    }
  }
}


/* Susceptible a expuesto. Bajo */
void reaction1(std::map<std::string, group> &Val, std::map<std::string, group> &Vba, Crandom &ran, Workers *altos, Workers *bajos, std::vector<lognormal_d> &dist, int agentI){
  int agentS;
  unsigned int index;
  double auxind;
  if(agentI > -1){
    auxind = (ran.r()*(Vba["SUS"].size() + Vba["SUST"].size() + (1-MyCons.alpha)*Vba["SUSA"].size()));

    /* Susceptible a Expuesto */
    if(auxind < Vba["SUS"].size()){
      index = (int)(ran.r()*Vba["SUS"].size());      agentS = Vba["SUS"][index];
      mother_reaction(Vba["SUS"], Vba["EXP"], index, bajos, SUS, EXP);
      bajos[Vba["EXP"].back()].DF = agentI;
    }

    /* Susceptible testeado a Expuesto testeado */
    else if(auxind < Vba["SUS"].size() + Vba["SUST"].size()){
      index = (int)(ran.r()*Vba["SUST"].size());      agentS = Vba["SUST"][index];
      mother_reaction(Vba["SUST"], Vba["EXPT"], index, bajos, SUST, EXPT);
      bajos[Vba["EXPT"].back()].DF = agentI;
    }
    /* Susceptible aislado a Expuesto aislado */
    else{
      index = (int)(ran.r()*Vba["SUSA"].size());      agentS = Vba["SUSA"][index];
      mother_reaction(Vba["SUSA"], Vba["EXPA"], index, bajos, SUSA, EXPA);
      bajos[Vba["EXPA"].back()].DF = agentI;
    }

    /* Guardo a la persona que infecté */
    if(agentI < Na){altos[agentI].my_inf.push_back(agentS + Na);}
    else{bajos[agentI-Na].my_inf.push_back(agentS + Na);}
  }
  else{/* Si el que infecta es externo */
    agentS = infected(bajos, Nb);
    if(agentS != -1){/* Si la infección debido al externo se realizó */
      /* Susceptible a Expuesto */
      if(bajos[agentS].kind == SUS){
	group::iterator it = std::find(Vba["SUS"].begin(), Vba["SUS"].end(), agentS);
	index = std::distance(Vba["SUS"].begin(), it);
	mother_reaction(Vba["SUS"], Vba["EXP"], index, bajos, SUS, EXP);
      }
      /* Susceptible testeado a Expuesto testado */
      else if(bajos[agentS].kind == SUST){
	group::iterator it = std::find(Vba["SUST"].begin(), Vba["SUST"].end(), agentS);
	index = std::distance(Vba["SUST"].begin(), it);
	mother_reaction(Vba["SUST"], Vba["EXPT"], index, bajos, SUST, EXPT);
      }
      /* Susceptible aislado a Expuesto aislado */
      else{
	group::iterator it = std::find(Vba["SUSA"].begin(), Vba["SUSA"].end(), agentS);
	index = std::distance(Vba["SUSA"].begin(), it);
	mother_reaction(Vba["SUSA"], Vba["EXPA"], index, bajos, SUSA, EXPA);
      }
    }
  }
}


/* Expuesto a Pre-sintomático. Alto */
void reaction2(std::map<std::string, group> &Val, std::map<std::string, group> &Vba, Crandom &ran, Workers *altos, Workers *bajos, std::vector<lognormal_d> &dist, int agentI){
  group aux;

  for(unsigned int i=0; i<Val["EXP"].size(); i++){if(altos[Val["EXP"][i]].tstate > TM){aux.push_back(Val["EXP"][i]);}} //Expuesto
  for(unsigned int i=0; i<Val["EXPT"].size(); i++){if(altos[Val["EXPT"][i]].tstate > TM){aux.push_back(Val["EXPT"][i]);}} //Expuesto testeado
  for(unsigned int i=0; i<Val["EXPA"].size(); i++){if(altos[Val["EXPA"][i]].tstate > TM){aux.push_back(Val["EXPA"][i]);}} //Expuesto aislado

  if(aux.size() != 0){
    double value = (double)ran.r();
    unsigned int index = index_time(aux, altos, dist[0], value);
    unsigned int agent = aux[index];
    group::iterator it;
    unsigned int ind;
    if(Val["EXP"].size() != 0){// Expuesto a Pre-sintomático
      it = std::find(Val["EXP"].begin(), Val["EXP"].end(), agent);
      ind = std::distance(Val["EXP"].begin(), it);
      if(ind < Val["EXP"].size()){mother_reaction(Val["EXP"], Val["PRE"], ind, altos, EXP, PRE);}
    }
    if(Val["EXPT"].size() != 0){// Expuesto testeado a Pre-sintomático
      it = std::find(Val["EXPT"].begin(), Val["EXPT"].end(), agent);
      ind = std::distance(Val["EXPT"].begin(), it);
      if(ind < Val["EXPT"].size()){mother_reaction(Val["EXPT"], Val["PRE"], ind, altos, EXPT, PRE);}
    }
    if(Val["EXPA"].size() != 0){// Expuesto aislado a Pre-sintomático aislado
      it = std::find(Val["EXPA"].begin(), Val["EXPA"].end(), agent);
      ind = std::distance(Val["EXPA"].begin(), it);
      if(ind < Val["EXPA"].size()){mother_reaction(Val["EXPA"], Val["PREA"], ind, altos, EXPA, PREA);}
    }
  }
}


/* Expuesto a Pre-sintomático. Bajo */
void reaction3(std::map<std::string, group> &Val, std::map<std::string, group> &Vba, Crandom &ran, Workers *altos, Workers *bajos, std::vector<lognormal_d> &dist, int agentI){
  group aux;

  for(unsigned int i=0; i<Vba["EXP"].size(); i++){if(bajos[Vba["EXP"][i]].tstate > TM){aux.push_back(Vba["EXP"][i]);}} //Expuesto
  for(unsigned int i=0; i<Vba["EXPT"].size(); i++){if(bajos[Vba["EXPT"][i]].tstate > TM){aux.push_back(Vba["EXPT"][i]);}} //Expuesto testeado
  for(unsigned int i=0; i<Vba["EXPA"].size(); i++){if(bajos[Vba["EXPA"][i]].tstate > TM){aux.push_back(Vba["EXPA"][i]);}} //Expuesto aislado

  if(aux.size() != 0){
    double value = (double)ran.r();
    unsigned int index = index_time(aux, bajos, dist[1], value);
    unsigned int agent = aux[index];
    group::iterator it;
    unsigned int ind;
    if(Vba["EXP"].size() != 0){// Expuesto a Pre-sintomático
      it = std::find(Vba["EXP"].begin(), Vba["EXP"].end(), agent);
      ind = std::distance(Vba["EXP"].begin(), it);
      if(ind < Vba["EXP"].size()){mother_reaction(Vba["EXP"], Vba["PRE"], ind, bajos, EXP, PRE);}
    }
    if(Vba["EXPT"].size() != 0){// Expuesto testeado a Pre-sintomático
      it = std::find(Vba["EXPT"].begin(), Vba["EXPT"].end(), agent);
      ind = std::distance(Vba["EXPT"].begin(), it);
      if(ind < Vba["EXPT"].size()){mother_reaction(Vba["EXPT"], Vba["PRE"], ind, bajos, EXPT, PRE);}
    }
    if(Vba["EXPA"].size() != 0){// Expuesto aislado a Pre-sintomático aislado
      it = std::find(Vba["EXPA"].begin(), Vba["EXPA"].end(), agent);
      ind = std::distance(Vba["EXPA"].begin(), it);
      if(ind < Vba["EXPA"].size()){mother_reaction(Vba["EXPA"], Vba["PREA"], ind, bajos, EXPA, PREA);}
    }
  }
}


/* Pre-sintomático a Leve. Alto */
void reaction4(std::map<std::string, group> &Val, std::map<std::string, group> &Vba, Crandom &ran, Workers *altos, Workers *bajos, std::vector<lognormal_d> &dist, int agentI){
  group aux;

  for(unsigned int i=0; i<Val["PRE"].size(); i++){if(altos[Val["PRE"][i]].tstate > TM){aux.push_back(Val["PRE"][i]);}} //Presintomático
  for(unsigned int i=0; i<Val["PRET"].size(); i++){if(altos[Val["PRET"][i]].tstate > TM){aux.push_back(Val["PRET"][i]);}} //Presintomático testeado
  for(unsigned int i=0; i<Val["PREA"].size(); i++){if(altos[Val["PREA"][i]].tstate > TM){aux.push_back(Val["PREA"][i]);}} //Presintomático aislado

  if(aux.size() != 0){
    double value = (double)ran.r();
    unsigned int index = index_time(aux, altos, dist[2], value);
    unsigned int agent = aux[index];
    group::iterator it;
    unsigned int ind;
    if(Val["PRE"].size() != 0){// Presintomático a Leve
      it = std::find(Val["PRE"].begin(), Val["PRE"].end(), agent);
      ind = std::distance(Val["PRE"].begin(), it);
      if(ind < Val["PRE"].size()){
	mother_reaction(Val["PRE"], Val["MSYM"], ind, altos, PRE, MSYM);
	if(ran.r() < MyCons.iota){// Leve a Leve aislado
	  if(MyCons.AisLev){ // Aislamiento de leves
	    tested_reaction(Val["MSYM"], Val["MSYMA"], Val["MSYM"].size()-1, altos, MSYM, MSYMA, 0.0, false);
	    tested_lev_ais(Val["MSYMA"].back(), altos, ran.r(), true);
	  }
	  else{ // No aislamiento de leves
	    tested_reaction(Val["MSYM"], Val["MSYMT"], Val["MSYM"].size()-1, altos, MSYM, MSYMT, ran.r(), true);
	  }
	}
      }
    }

    if(Val["PRET"].size() != 0){// Presintomático testeado a Leve testeado
      it = std::find(Val["PRET"].begin(), Val["PRET"].end(), agent);
      ind = std::distance(Val["PRET"].begin(), it);
      if(ind < Val["PRET"].size()){
	mother_reaction(Val["PRET"], Val["MSYMT"], ind, altos, PRET, MSYMT);
	if(ran.r() < MyCons.iota){// Leve testeado a Leve aislado
	  if(MyCons.AisLev){ // Aislamiento de leves
	    tested_reaction(Val["MSYMT"], Val["MSYMA"], Val["MSYMT"].size()-1, altos, MSYMT, MSYMA, 0.0, false);
	    tested_lev_ais(Val["MSYMA"].back(), altos, ran.r(), true);
	  }
	}
      }
    }
    if(Val["PREA"].size() != 0){// Presintomático aislado a Leve aislado
      it = std::find(Val["PREA"].begin(), Val["PREA"].end(), agent);
      ind = std::distance(Val["PREA"].begin(), it);
      if(ind < Val["PREA"].size()){mother_reaction(Val["PREA"], Val["MSYMA"], ind, altos, PREA, MSYMA);}
    }
  }
}


/* Pre-sintomático a Leve. Bajo */
void reaction5(std::map<std::string, group> &Val, std::map<std::string, group> &Vba, Crandom &ran, Workers *altos, Workers *bajos, std::vector<lognormal_d> &dist, int agentI){
  group aux;

  for(unsigned int i=0; i<Vba["PRE"].size(); i++){if(bajos[Vba["PRE"][i]].tstate > TM){aux.push_back(Vba["PRE"][i]);}} //Presintomático
  for(unsigned int i=0; i<Vba["PRET"].size(); i++){if(bajos[Vba["PRET"][i]].tstate > TM){aux.push_back(Vba["PRET"][i]);}} //Presintomático testeado
  for(unsigned int i=0; i<Vba["PREA"].size(); i++){if(bajos[Vba["PREA"][i]].tstate > TM){aux.push_back(Vba["PREA"][i]);}} //Presintomático aislado

  if(aux.size() != 0){
    double value = (double)ran.r();
    unsigned int index = index_time(aux, bajos, dist[3], value);
    unsigned int agent = aux[index];
    group::iterator it;
    unsigned int ind;
    if(Vba["PRE"].size() != 0){// Presintomático a Leve
      it = std::find(Vba["PRE"].begin(), Vba["PRE"].end(), agent);
      ind = std::distance(Vba["PRE"].begin(), it);
      if(ind < Vba["PRE"].size()){
	mother_reaction(Vba["PRE"], Vba["MSYM"], ind, bajos, PRE, MSYM);
	if(ran.r() < MyCons.iota){// Leve a Leve aislado
	  if(MyCons.AisLev){ // Aislamiento de leves
	    tested_reaction(Vba["MSYM"], Vba["MSYMA"], Vba["MSYM"].size()-1, bajos, MSYM, MSYMA, 0.0, false);
	    tested_lev_ais(Vba["MSYMA"].back(), bajos, ran.r(), true);
	  }
	  else{
	    tested_reaction(Vba["MSYM"], Vba["MSYMT"], Vba["MSYM"].size()-1, bajos, MSYM, MSYMT, ran.r(), true);
	  }
	}
      }
    }
    if(Vba["PRET"].size() != 0){// Presintomático testeado a Leve testeado
      it = std::find(Vba["PRET"].begin(), Vba["PRET"].end(), agent);
      ind = std::distance(Vba["PRET"].begin(), it);
      if(ind < Vba["PRET"].size()){
	mother_reaction(Vba["PRET"], Vba["MSYMT"], ind, bajos, PRET, MSYMT);
	if(ran.r() < MyCons.iota){// Leve testeado a Leve aislado
	  if(MyCons.AisLev){ // Aislamiento de leves
	    tested_reaction(Vba["MSYMT"], Vba["MSYMA"], Vba["MSYMT"].size()-1, bajos, MSYMT, MSYMA, 0.0, false);
	    tested_lev_ais(Vba["MSYMA"].back(), bajos, ran.r(), true);
	  }
	}
      }
    }
    if(Vba["PREA"].size() != 0){// Presintomático aislado a Leve aislado
      it = std::find(Vba["PREA"].begin(), Vba["PREA"].end(), agent);
      ind = std::distance(Vba["PREA"].begin(), it);
      if(ind < Vba["PREA"].size()){mother_reaction(Vba["PREA"], Vba["MSYMA"], ind, bajos, PREA, MSYMA);}
    }
  }
}


/* Pre-sintomático a Grave. Alto */
void reaction6(std::map<std::string, group> &Val, std::map<std::string, group> &Vba, Crandom &ran, Workers *altos, Workers *bajos, std::vector<lognormal_d> &dist, int agentI){
  group aux;

  for(unsigned int i=0; i<Val["PRE"].size(); i++){if(altos[Val["PRE"][i]].tstate > TM){aux.push_back(Val["PRE"][i]);}} //Presintomático
  for(unsigned int i=0; i<Val["PRET"].size(); i++){if(altos[Val["PRET"][i]].tstate > TM){aux.push_back(Val["PRET"][i]);}} //Presintomático testeado
  for(unsigned int i=0; i<Val["PREA"].size(); i++){if(altos[Val["PREA"][i]].tstate > TM){aux.push_back(Val["PREA"][i]);}} //Presintomático aislado

  if(aux.size() != 0){
    double value = (double)ran.r();
    unsigned int index = index_time(aux, altos, dist[4], value);
    unsigned int agent = aux[index];
    group::iterator it;
    unsigned int ind;
    if(Val["PRE"].size() != 0){// Presintomático a Grave
      it = std::find(Val["PRE"].begin(), Val["PRE"].end(), agent);
      ind = std::distance(Val["PRE"].begin(), it);
      if(ind < Val["PRE"].size()){
	mother_reaction(Val["PRE"], Val["SSYMA"], ind, altos, PRE, SSYMA);
	aux_main(1, Val, Vba, altos, bajos, ran, MyCons.phi1, MyCons.mu, true, "SSYMA", true);
      }
    }
    if(Val["PRET"].size() != 0){// Presintomático testeado a Grave
      it = std::find(Val["PRET"].begin(), Val["PRET"].end(), agent);
      ind = std::distance(Val["PRET"].begin(), it);
      if(ind < Val["PRET"].size()){
	mother_reaction(Val["PRET"], Val["SSYMA"], ind, altos, PRET, SSYMA);
	aux_main(1, Val, Vba, altos, bajos, ran, MyCons.phi1, MyCons.mu, true, "SSYMA", true);
      }
    }
    if(Val["PREA"].size() != 0){// Presintomáticos aislado a Grave
      it = std::find(Val["PREA"].begin(), Val["PREA"].end(), agent);
      ind = std::distance(Val["PREA"].begin(), it);
      if(ind < Val["PREA"].size()){mother_reaction(Val["PREA"], Val["SSYMA"], ind, altos, PREA, SSYMA);}
    }
  }
}


/* Pre-sintomático a Grave. Bajo */
void reaction7(std::map<std::string, group> &Val, std::map<std::string, group> &Vba, Crandom &ran, Workers *altos, Workers *bajos, std::vector<lognormal_d> &dist, int agentI){
  group aux;

  for(unsigned int i=0; i<Vba["PRE"].size(); i++){if(bajos[Vba["PRE"][i]].tstate > TM){aux.push_back(Vba["PRE"][i]);}} //Presintomático
  for(unsigned int i=0; i<Vba["PRET"].size(); i++){if(bajos[Vba["PRET"][i]].tstate > TM){aux.push_back(Vba["PRET"][i]);}} //Presintomático testeado
  for(unsigned int i=0; i<Vba["PREA"].size(); i++){if(bajos[Vba["PREA"][i]].tstate > TM){aux.push_back(Vba["PREA"][i]);}} //Presitnomático aislado

  if(aux.size() != 0){
    double value = (double)ran.r();
    unsigned int index = index_time(aux, bajos, dist[5], value);
    unsigned int agent = aux[index];
    group::iterator it;
    unsigned int ind;
    if(Vba["PRE"].size() != 0){// Presintomático a Grave
      it = std::find(Vba["PRE"].begin(), Vba["PRE"].end(), agent);
      ind = std::distance(Vba["PRE"].begin(), it);
      if(ind < Vba["PRE"].size()){
	mother_reaction(Vba["PRE"], Vba["SSYMA"], ind, bajos, PRE, SSYMA);
	aux_main(1, Val, Vba, altos, bajos, ran, MyCons.mu, MyCons.chi, false, "SSYMA", true);
      }
    }
    if(Vba["PRET"].size() != 0){
      it = std::find(Vba["PRET"].begin(), Vba["PRET"].end(), agent);
      ind = std::distance(Vba["PRET"].begin(), it);
      if(ind < Vba["PRET"].size()){// Presintomático testeado a Grave
	mother_reaction(Vba["PRET"], Vba["SSYMA"], ind, bajos, PRET, SSYMA);
	aux_main(1, Val, Vba, altos, bajos, ran, MyCons.mu, MyCons.chi, false, "SSYMA", true);
      }
    }
    if(Vba["PREA"].size() != 0){// Presintomático aislado a Grave
      it = std::find(Vba["PREA"].begin(), Vba["PREA"].end(), agent);
      ind = std::distance(Vba["PREA"].begin(), it);
      if(ind < Vba["PREA"].size()){mother_reaction(Vba["PREA"], Vba["SSYMA"], ind, bajos, PREA, SSYMA);}
    }
  }
}


/* Pre-sintomático a Recuperado. Alto */
void reaction8(std::map<std::string, group> &Val, std::map<std::string, group> &Vba, Crandom &ran, Workers *altos, Workers *bajos, std::vector<lognormal_d> &dist, int agentI){
  group aux;

  for(unsigned int i=0; i<Val["PRE"].size(); i++){if(altos[Val["PRE"][i]].tstate > TM){aux.push_back(Val["PRE"][i]);}} //Presintomático
  for(unsigned int i=0; i<Val["PRET"].size(); i++){if(altos[Val["PRET"][i]].tstate > TM){aux.push_back(Val["PRET"][i]);}} //Presintomático testeado
  for(unsigned int i=0; i<Val["PREA"].size(); i++){if(altos[Val["PREA"][i]].tstate > TM){aux.push_back(Val["PREA"][i]);}} //Presintomático aislado

  if(aux.size() != 0){
    double value = (double)ran.r();
    unsigned int index = index_time(aux, altos, dist[6], value);
    unsigned int agent = aux[index];
    group::iterator it;
    unsigned int ind;
    if(Val["PRE"].size() != 0){// Presintomático a Recuperado no-detectado
      it = std::find(Val["PRE"].begin(), Val["PRE"].end(), agent);
      ind = std::distance(Val["PRE"].begin(), it);
      if(ind < Val["PRE"].size()){mother_reaction(Val["PRE"], Val["RECI"], ind, altos, PRE, RECI);}
    }
    if(Val["PRET"].size() != 0){// Presintomático testeado a Recuperado detectado o no-detectado
      it = std::find(Val["PRET"].begin(), Val["PRET"].end(), agent);
      ind = std::distance(Val["PRET"].begin(), it);
      if(ind < Val["PRET"].size()){
        if(ran.r() < MyCons.xi){//Presintomático testeado a Recuperado detectado
          mother_reaction(Val["PRET"], Val["RECA"], ind, altos, PRET, RECA);
 	  aux_main(1, Val, Vba, altos, bajos, ran, MyCons.phi1, MyCons.mu, true, "RECA", true);
        }
        else{//Presintomático testeado a Recuperado detectado
          mother_reaction(Val["PRET"], Val["RECT"], ind, altos, PRET, RECT);
        }
      }
    }
    if(Val["PREA"].size() != 0){// Presintomático aislado a Recuperado detectado
      it = std::find(Val["PREA"].begin(), Val["PREA"].end(), agent);
      ind = std::distance(Val["PREA"].begin(), it);
      if(ind < Val["PREA"].size()){mother_reaction(Val["PREA"], Val["RECA"], ind, altos, PREA, RECA);}
    }
  }
}


/* Pre-sintomático a Recuperado. Bajo */
void reaction9(std::map<std::string, group> &Val, std::map<std::string, group> &Vba, Crandom &ran, Workers *altos, Workers *bajos, std::vector<lognormal_d> &dist, int agentI){
  group aux;

  for(unsigned int i=0; i<Vba["PRE"].size(); i++){if(bajos[Vba["PRE"][i]].tstate > TM){aux.push_back(Vba["PRE"][i]);}} //Presintomático
  for(unsigned int i=0; i<Vba["PRET"].size(); i++){if(bajos[Vba["PRET"][i]].tstate > TM){aux.push_back(Vba["PRET"][i]);}} //Presintomático testeado
  for(unsigned int i=0; i<Vba["PREA"].size(); i++){if(bajos[Vba["PREA"][i]].tstate > TM){aux.push_back(Vba["PREA"][i]);}} //Presintomático aislado

  if(aux.size() != 0){
    double value = (double)ran.r();
    unsigned int index = index_time(aux, bajos, dist[7], value);
    unsigned int agent = aux[index];
    group::iterator it;
    unsigned int ind;
    if(Vba["PRE"].size() != 0){// Presintomático a Recuperado no-detectado
      it = std::find(Vba["PRE"].begin(), Vba["PRE"].end(), agent);
      ind = std::distance(Vba["PRE"].begin(), it);
      if(ind < Vba["PRE"].size()){mother_reaction(Vba["PRE"], Vba["RECI"], ind, bajos, PRE, RECI);}
    }
    if(Vba["PRET"].size() != 0){// Presintomático testeado a Recuperado detectado o no-detectado
      it = std::find(Vba["PRET"].begin(), Vba["PRET"].end(), agent);
      ind = std::distance(Vba["PRET"].begin(), it);
      if(ind < Vba["PRET"].size()){
        if(ran.r() < MyCons.xi){// Presintomático testeado a Recuperado detectado
          mother_reaction(Vba["PRET"], Vba["RECA"], ind, bajos, PRET, RECA);
 	  aux_main(1, Val, Vba, altos, bajos, ran, MyCons.mu, MyCons.chi, false, "RECA", true);
        }
        else{// Presintomático testeado a Recuperado no-detectado
          mother_reaction(Vba["PRET"], Vba["RECI"], ind, bajos, PRET, RECI);
        }
      }
    }
    if(Vba["PREA"].size() != 0){// Presintomático aislado a Recuperado detectado
      it = std::find(Vba["PREA"].begin(), Vba["PREA"].end(), agent);
      ind = std::distance(Vba["PREA"].begin(), it);
      if(ind < Vba["PREA"].size()){mother_reaction(Vba["PREA"], Vba["RECA"], ind, bajos, PREA, RECA);}
    }
  }
}


/* Leve a Recuperado. Alto */
void reaction10(std::map<std::string, group> &Val, std::map<std::string, group> &Vba, Crandom &ran, Workers *altos, Workers *bajos, std::vector<lognormal_d> &dist, int agentI){
  group aux;

  for(unsigned int i=0; i<Val["MSYM"].size(); i++){if(altos[Val["MSYM"][i]].tstate > TM){aux.push_back(Val["MSYM"][i]);}} //Leve
  for(unsigned int i=0; i<Val["MSYMT"].size(); i++){if(altos[Val["MSYMT"][i]].tstate > TM){aux.push_back(Val["MSYMT"][i]);}} //Leve testeado
  for(unsigned int i=0; i<Val["MSYMA"].size(); i++){if(altos[Val["MSYMA"][i]].tstate > TM){aux.push_back(Val["MSYMA"][i]);}} //Leve aislado

  if(aux.size() != 0){
    double value = (double)ran.r();
    unsigned int index = index_time(aux, altos, dist[8], value);
    unsigned int agent = aux[index];
    group::iterator it;
    unsigned int ind;
    if(Val["MSYM"].size() != 0){// Leve a Recuperado no-detectado
      it = std::find(Val["MSYM"].begin(), Val["MSYM"].end(), agent);
      ind = std::distance(Val["MSYM"].begin(), it);
      if(ind < Val["MSYM"].size()){mother_reaction(Val["MSYM"], Val["RECI"], ind, altos, MSYM, RECI);}
    }
    if(Val["MSYMT"].size() != 0){// Leve testeado a Recuperado detectado o no-detectado
      it = std::find(Val["MSYMT"].begin(), Val["MSYMT"].end(), agent);
      ind = std::distance(Val["MSYMT"].begin(), it);
      if(ind < Val["MSYMT"].size()){
        if(ran.r() < MyCons.xi){// Leve testeado a Recuperado detectado
          mother_reaction(Val["MSYMT"], Val["RECA"], ind, altos, MSYMT, RECA);
 	  aux_main(1, Val, Vba, altos, bajos, ran, MyCons.phi1, MyCons.mu, true, "RECA", false);
        }
        else{// Leve testeado a Recuperado no-detectado
          mother_reaction(Val["MSYMT"], Val["RECI"], ind, altos, MSYMT, RECI);
        }
      }
    }
    if(Val["MSYMA"].size() != 0){// Leve aislado a Recuperado detectado
      it = std::find(Val["MSYMA"].begin(), Val["MSYMA"].end(), agent);
      ind = std::distance(Val["MSYMA"].begin(), it);
      if(ind < Val["MSYMA"].size()){mother_reaction(Val["MSYMA"], Val["RECA"], ind, altos, MSYMA, RECA);}
    }
  }
}


/* Leve a Recuperado. Bajo */
void reaction11(std::map<std::string, group> &Val, std::map<std::string, group> &Vba, Crandom &ran, Workers *altos, Workers *bajos, std::vector<lognormal_d> &dist, int agentI){
  group aux;

  for(unsigned int i=0; i<Vba["MSYM"].size(); i++){if(bajos[Vba["MSYM"][i]].tstate > TM){aux.push_back(Vba["MSYM"][i]);}} //Leve
  for(unsigned int i=0; i<Vba["MSYMT"].size(); i++){if(bajos[Vba["MSYMT"][i]].tstate > TM){aux.push_back(Vba["MSYMT"][i]);}} //Leve testeado
  for(unsigned int i=0; i<Vba["MSYMA"].size(); i++){if(bajos[Vba["MSYMA"][i]].tstate > TM){aux.push_back(Vba["MSYMA"][i]);}} //Leve aislado

  if(aux.size() != 0){
    double value = (double)ran.r();
    unsigned int index = index_time(aux, bajos, dist[9], value);
    unsigned int agent = aux[index];
    group::iterator it;
    unsigned int ind;
    if(Vba["MSYM"].size() != 0){// Leve a Recuperado no-detectado
      it = std::find(Vba["MSYM"].begin(), Vba["MSYM"].end(), agent);
      ind = std::distance(Vba["MSYM"].begin(), it);
      if(ind < Vba["MSYM"].size()){mother_reaction(Vba["MSYM"], Vba["RECI"], ind, bajos, MSYM, RECI);}
    }
    if(Vba["MSYMT"].size() != 0){// Leve testeado a Recuperado detectado o no-detectado
      it = std::find(Vba["MSYMT"].begin(), Vba["MSYMT"].end(), agent);
      ind = std::distance(Vba["MSYMT"].begin(), it);
      if(ind < Vba["MSYMT"].size()){
        if(ran.r() < MyCons.xi){// Leve testeado a Recuperado detectado
          mother_reaction(Vba["MSYMT"], Vba["RECA"], ind, bajos, MSYMT, RECA);
 	  aux_main(1, Val, Vba, altos, bajos, ran, MyCons.mu, MyCons.chi, false, "RECA", false);
        }
        else{// Leve testeado a Recuperado no-detectado
          mother_reaction(Vba["MSYMT"], Vba["RECI"], ind, bajos, MSYMT, RECI);
        }
      }
    }
    if(Vba["MSYMA"].size() != 0){// Leve aislado a Recuperado detectado
      it = std::find(Vba["MSYMA"].begin(), Vba["MSYMA"].end(), agent);
      ind = std::distance(Vba["MSYMA"].begin(), it);
      if(ind < Vba["MSYMA"].size()){mother_reaction(Vba["MSYMA"], Vba["RECA"], ind, bajos, MSYMA, RECA);}
    }
  }
}


/* Grave a Recuperado. Alto */
void reaction12(std::map<std::string, group> &Val, std::map<std::string, group> &Vba, Crandom &ran, Workers *altos, Workers *bajos, std::vector<lognormal_d> &dist, int agentI){
  group aux;

  for(unsigned int i=0; i<Val["SSYMA"].size(); i++){if(altos[Val["SSYMA"][i]].tstate > TM){aux.push_back(Val["SSYMA"][i]);}} //Grave

  if(aux.size() != 0){// Grave a Recuperado detectado
    unsigned int index = index_time(aux, altos, dist[10], (double)ran.r());
    unsigned int agent = aux[index];
    group::iterator it = std::find(Val["SSYMA"].begin(), Val["SSYMA"].end(), agent);
    unsigned int ind = std::distance(Val["SSYMA"].begin(), it);
    if(ind < Val["SSYMA"].size()){mother_reaction(Val["SSYMA"], Val["RECA"], ind, altos, SSYMA, RECA);}
  }
}


/* Grave a Recuperado. Bajo */
void reaction13(std::map<std::string, group> &Val, std::map<std::string, group> &Vba, Crandom &ran, Workers *altos, Workers *bajos, std::vector<lognormal_d> &dist, int agentI){
  group aux;

  for(unsigned int i=0; i<Vba["SSYMA"].size(); i++){if(bajos[Vba["SSYMA"][i]].tstate > TM){aux.push_back(Vba["SSYMA"][i]);}} //Grave

  if(aux.size() != 0){// Grave a Recuperado detectado
    unsigned int index = index_time(aux, bajos, dist[11], (double)ran.r());
    unsigned int agent = aux[index];
    group::iterator it = std::find(Vba["SSYMA"].begin(), Vba["SSYMA"].end(), agent);
    unsigned int ind = std::distance(Vba["SSYMA"].begin(), it);
    if(ind < Vba["SSYMA"].size()){mother_reaction(Vba["SSYMA"], Vba["RECA"], ind, bajos, SSYMA, RECA);}
  }
}
