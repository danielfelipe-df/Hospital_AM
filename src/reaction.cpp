#include <iostream>
#include <reaction.h>
#include <trace.h>
#include <test.h>

void mother_reaction(grupo &Out, grupo &In, int index, Workers *family, int typeout, int typein){
  int agent = Out[index];
  Out.erase(Out.begin() + index);
  In.push_back(agent);
  family[agent].change(typein, typeout);
  family[agent].tstate = 0.0;
}


int index_time(grupo &Out, Workers *family, lognormal_d &dist, double value){
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
  int kind;
  for(int i=0; i<(max-1); i++){
    kind = family[i].kind;
    if(kind == 0 || kind == 1 || kind == 2){return i;}
  }
  return max-1;
}



/* Susceptible a expuesto. Alto */
void reaction0(std::vector<grupo> &Val, std::vector<grupo> &Vba, Crandom &ran, Workers *altos, Workers *bajos, std::vector<lognormal_d> &dist, int agentI){
  int agentS;
  unsigned int index;
  double auxind;
  if(agentI > -1){/* Si el que infecta no es externo */
    auxind = ran.r()*(Val[0].size() + Val[1].size() + (1-MyCons.alpha)*Val[2].size());

    /* Susceptible a Expuesto */
    if(auxind < Val[0].size()){index = (int)(ran.r()*Val[0].size());      agentS = Val[0][index];    mother_reaction(Val[0], Val[3], index, altos, 0, 3);      altos[Val[3].back()].DF = agentI;}
    /* Susceptible testeado a Expuesto testado */
    else if(auxind < Val[0].size() + Val[1].size()){index = (int)(ran.r()*Val[1].size());      agentS = Val[1][index];    mother_reaction(Val[1], Val[4], index, altos, 1, 4);      altos[Val[4].back()].DF = agentI;}
    /* Susceptible aislado a Expuesto aislado */
    else{index = (int)(ran.r()*Val[2].size());      agentS = Val[2][index];    mother_reaction(Val[2], Val[5], index, altos, 2, 5);      altos[Val[5].back()].DF = agentI;}

    /* Guardo a la persona que infecté */
    if(agentI < Na){altos[agentI].my_inf.push_back(agentS);}
    else{bajos[agentI-Na].my_inf.push_back(agentS);}
  }
  else{/* Si el que infecta es externo */
    agentS = infected(altos, Na);
    if(agentS != -1){/* Si la infección debido al externo se realizó */
      /* Susceptible a Expuesto */
      if(altos[agentS].kind == 0){
	std::vector<int>::iterator it = std::find(Val[0].begin(), Val[0].end(), agentS);
	index = std::distance(Val[0].begin(), it);
	mother_reaction(Val[0], Val[3], index, altos, 0, 3);}
      
      /* Susceptible testeado a Expuesto testado */
      else if(altos[agentS].kind == 1){
	std::vector<int>::iterator it = std::find(Val[1].begin(), Val[1].end(), agentS);
	index = std::distance(Val[1].begin(), it);
	mother_reaction(Val[1], Val[4], index, altos, 1, 4);
      }

      /* Susceptible aislado a Expuesto aislado */
      else{
	std::vector<int>::iterator it = std::find(Val[2].begin(), Val[2].end(), agentS);
	index = std::distance(Val[2].begin(), it);
	mother_reaction(Val[2], Val[5], index, altos, 2, 5);
      }
    }
  }
}


/* Susceptible a expuesto. Bajo */
void reaction1(std::vector<grupo> &Val, std::vector<grupo> &Vba, Crandom &ran, Workers *altos, Workers *bajos, std::vector<lognormal_d> &dist, int agentI){
  int agentS;
  unsigned int index;
  double auxind;
  if(agentI > -1){
    auxind = (ran.r()*(Vba[0].size() + Vba[1].size() + (1-MyCons.alpha)*Vba[2].size()));
    
    /* Susceptible a Expuesto */
    if(auxind < Vba[0].size()){index = (int)(ran.r()*Vba[0].size());      agentS = Vba[0][index];    mother_reaction(Vba[0], Vba[3], index, bajos, 0, 3);      bajos[Vba[3].back()].DF = agentI;}
    /* Susceptible testeado a Expuesto testeado */
    else if(auxind < Vba[0].size() + Vba[1].size()){index = (int)(ran.r()*Vba[1].size());      agentS = Vba[1][index];    mother_reaction(Vba[1], Vba[4], index, bajos, 1, 4);}
    /* Susceptible aislado a Expuesto aislado */
    else{index = (int)(ran.r()*Vba[2].size());      agentS = Vba[2][index];    mother_reaction(Vba[2], Vba[5], index, bajos, 2, 5);}
    
    /* Guardo a la persona que infecté */
    if(agentI < Na){altos[agentI].my_inf.push_back(agentS + Na);}
    else{bajos[agentI-Na].my_inf.push_back(agentS + Na);}
  }
  else{/* Si el que infecta es externo */
    agentS = infected(bajos, Nb);
    if(agentS != -1){/* Si la infección debido al externo se realizó */
      /* Susceptible a Expuesto */
      if(bajos[agentS].kind == 0){
	std::vector<int>::iterator it = std::find(Vba[0].begin(), Vba[0].end(), agentS);
	index = std::distance(Vba[0].begin(), it);
	mother_reaction(Vba[0], Vba[3], index, bajos, 0, 3);}
      /* Susceptible testeado a Expuesto testado */
      else if(bajos[agentS].kind == 1){
	std::vector<int>::iterator it = std::find(Vba[1].begin(), Vba[1].end(), agentS);
	index = std::distance(Vba[1].begin(), it);
	mother_reaction(Vba[1], Vba[4], index, bajos, 1, 4);
      }
      /* Susceptible aislado a Expuesto aislado */
      else{
	std::vector<int>::iterator it = std::find(Vba[2].begin(), Vba[2].end(), agentS);
	index = std::distance(Vba[2].begin(), it);
	mother_reaction(Vba[2], Vba[5], index, bajos, 2, 5);
      }
    }
  }
}


/* Expuesto a Pre-sintomático. Alto */
void reaction2(std::vector<grupo> &Val, std::vector<grupo> &Vba, Crandom &ran, Workers *altos, Workers *bajos, std::vector<lognormal_d> &dist, int agentI){
  std::vector<int> aux;

  for(unsigned int i=0; i<Val[3].size(); i++){if(altos[Val[3][i]].tstate > TM){aux.push_back(Val[3][i]);}} //Expuesto
  for(unsigned int i=0; i<Val[4].size(); i++){if(altos[Val[4][i]].tstate > TM){aux.push_back(Val[4][i]);}} //Expuesto testeado
  for(unsigned int i=0; i<Val[5].size(); i++){if(altos[Val[5][i]].tstate > TM){aux.push_back(Val[5][i]);}} //Expuesto aislado

  if(aux.size() != 0){
    double value = (double)ran.r();
    unsigned int index = index_time(aux, altos, dist[0], value);
    unsigned int agent = aux[index];
    std::vector<int>::iterator it;
    unsigned int ind;
    if(Val[3].size() != 0){// Expuesto a Pre-sintomático
      it = std::find(Val[3].begin(), Val[3].end(), agent);      ind = std::distance(Val[3].begin(), it);
      if(ind < Val[3].size()){mother_reaction(Val[3], Val[6], ind, altos, 3, 6);}
    }
    if(Val[4].size() != 0){// Expuesto testeado a Pre-sintomático
      it = std::find(Val[4].begin(), Val[4].end(), agent);      ind = std::distance(Val[4].begin(), it);
      if(ind < Val[4].size()){mother_reaction(Val[4], Val[7], ind, altos, 4, 7);}
    }
    if(Val[5].size() != 0){// Expuesto aislado a Pre-sintomático aislado
      it = std::find(Val[5].begin(), Val[5].end(), agent);      ind = std::distance(Val[5].begin(), it);
      if(ind < Val[5].size()){mother_reaction(Val[5], Val[8], ind, altos, 5, 8);}
    }
  }
}


/* Expuesto a Pre-sintomático. Bajo */
void reaction3(std::vector<grupo> &Val, std::vector<grupo> &Vba, Crandom &ran, Workers *altos, Workers *bajos, std::vector<lognormal_d> &dist, int agentI){
  std::vector<int> aux;

  for(unsigned int i=0; i<Vba[3].size(); i++){if(bajos[Vba[3][i]].tstate > TM){aux.push_back(Vba[3][i]);}} //Expuesto
  for(unsigned int i=0; i<Vba[4].size(); i++){if(bajos[Vba[4][i]].tstate > TM){aux.push_back(Vba[4][i]);}} //Expuesto testeado
  for(unsigned int i=0; i<Vba[5].size(); i++){if(bajos[Vba[5][i]].tstate > TM){aux.push_back(Vba[5][i]);}} //Expuesto aislado

  if(aux.size() != 0){
    double value = (double)ran.r();
    unsigned int index = index_time(aux, bajos, dist[1], value);
    unsigned int agent = aux[index];
    std::vector<int>::iterator it;
    unsigned int ind;
    if(Vba[3].size() != 0){// Expuesto a Pre-sintomático
      it = std::find(Vba[3].begin(), Vba[3].end(), agent);      ind = std::distance(Vba[3].begin(), it);
      if(ind < Vba[3].size()){mother_reaction(Vba[3], Vba[6], ind, bajos, 3, 6);}
    }
    if(Vba[4].size() != 0){// Expuesto testeado a Pre-sintomático
      it = std::find(Vba[4].begin(), Vba[4].end(), agent);      ind = std::distance(Vba[4].begin(), it);
      if(ind < Vba[4].size()){mother_reaction(Vba[4], Vba[6], ind, bajos, 4, 6);}
    }
    if(Vba[5].size() != 0){// Expuesto aislado a Pre-sintomático aislado
      it = std::find(Vba[5].begin(), Vba[5].end(), agent);      ind = std::distance(Vba[5].begin(), it);
      if(ind < Vba[5].size()){mother_reaction(Vba[5], Vba[8], ind, bajos, 5, 8);}
    }
  }
}


/* Pre-sintomático a Leve. Alto */
void reaction4(std::vector<grupo> &Val, std::vector<grupo> &Vba, Crandom &ran, Workers *altos, Workers *bajos, std::vector<lognormal_d> &dist, int agentI){
  std::vector<int> aux;

  for(unsigned int i=0; i<Val[6].size(); i++){if(altos[Val[6][i]].tstate > TM){aux.push_back(Val[6][i]);}} //Presintomático
  for(unsigned int i=0; i<Val[7].size(); i++){if(altos[Val[7][i]].tstate > TM){aux.push_back(Val[7][i]);}} //Presintomático testeado
  for(unsigned int i=0; i<Val[8].size(); i++){if(altos[Val[8][i]].tstate > TM){aux.push_back(Val[8][i]);}} //Presintomático aislado

  if(aux.size() != 0){
    double value = (double)ran.r();
    unsigned int index = index_time(aux, altos, dist[2], value);
    unsigned int agent = aux[index];
    std::vector<int>::iterator it;
    unsigned int ind;
    if(Val[6].size() != 0){// Presintomático a Leve
      it = std::find(Val[6].begin(), Val[6].end(), agent);      ind = std::distance(Val[6].begin(), it);
      if(ind < Val[6].size()){
	mother_reaction(Val[6], Val[9], ind, altos, 6, 9);
	if(ran.r() < MyCons.iota){// Leve a Leve aislado
	  if(MyCons.AisLev){ // Aislamiento de leves
	    tested_reaction(Val[9], Val[11], Val[9].size()-1, altos, 9, 11, 0.0, false);
	    tested_lev_ais(Val[11].back(), altos, ran.r(), true);
	  }
	  else{ // No aislamiento de leves
	    tested_reaction(Val[9], Val[10], Val[9].size()-1, altos, 9, 10, ran.r(), true);
	  }
	}
      }
    }
    if(Val[7].size() != 0){// Presintomático testeado a Leve testeado
      it = std::find(Val[7].begin(), Val[7].end(), agent);      ind = std::distance(Val[7].begin(), it);
      if(ind < Val[7].size()){
	mother_reaction(Val[7], Val[10], ind, altos, 7, 10);
	if(ran.r() < MyCons.iota){// Leve testeado a Leve aislado
	  if(MyCons.AisLev){ // Aislamiento de leves
	    tested_reaction(Val[10], Val[11], Val[10].size()-1, altos, 10, 11, 0.0, false);
	    tested_lev_ais(Val[11].back(), altos, ran.r(), true);
	  }
	}
      }
    }
    if(Val[8].size() != 0){// Presintomático aislado a Leve aislado
      it = std::find(Val[8].begin(), Val[8].end(), agent);      ind = std::distance(Val[8].begin(), it);
      if(ind < Val[8].size()){mother_reaction(Val[8], Val[11], ind, altos, 8, 11);}
    }
  }
}


/* Pre-sintomático a Leve. Bajo */
void reaction5(std::vector<grupo> &Val, std::vector<grupo> &Vba, Crandom &ran, Workers *altos, Workers *bajos, std::vector<lognormal_d> &dist, int agentI){
  std::vector<int> aux;

  for(unsigned int i=0; i<Vba[6].size(); i++){if(bajos[Vba[6][i]].tstate > TM){aux.push_back(Vba[6][i]);}} //Presintomático
  for(unsigned int i=0; i<Vba[7].size(); i++){if(bajos[Vba[7][i]].tstate > TM){aux.push_back(Vba[7][i]);}} //Presintomático testeado
  for(unsigned int i=0; i<Vba[8].size(); i++){if(bajos[Vba[8][i]].tstate > TM){aux.push_back(Vba[8][i]);}} //Presintomático aislado

  if(aux.size() != 0){
    double value = (double)ran.r();
    unsigned int index = index_time(aux, bajos, dist[3], value);
    unsigned int agent = aux[index];
    std::vector<int>::iterator it;
    unsigned int ind;
    if(Vba[6].size() != 0){// Presintomático a Leve
      it = std::find(Vba[6].begin(), Vba[6].end(), agent);      ind = std::distance(Vba[6].begin(), it);
      if(ind < Vba[6].size()){
	mother_reaction(Vba[6], Vba[9], ind, bajos, 6, 9);
	if(ran.r() < MyCons.iota){// Leve a Leve aislado
	  if(MyCons.AisLev){ // Aislamiento de leves
	    tested_reaction(Vba[9], Vba[11], Vba[9].size()-1, bajos, 9, 11, 0.0, false);
	    tested_lev_ais(Vba[11].back(), bajos, ran.r(), true);
	  }
	  else{
	    tested_reaction(Vba[9], Vba[10], Vba[9].size()-1, bajos, 9, 10, ran.r(), true);
	  }
	}
      }
    }
    if(Vba[7].size() != 0){// Presintomático testeado a Leve testeado
      it = std::find(Vba[7].begin(), Vba[7].end(), agent);      ind = std::distance(Vba[7].begin(), it);
      if(ind < Vba[7].size()){
	mother_reaction(Vba[7], Vba[10], ind, bajos, 7, 10);
	if(ran.r() < MyCons.iota){// Leve testeado a Leve aislado
	  if(MyCons.AisLev){ // Aislamiento de leves
	    tested_reaction(Vba[10], Vba[11], Vba[10].size()-1, bajos, 10, 11, 0.0, false);
	    tested_lev_ais(Vba[11].back(), bajos, ran.r(), true);
	  }
	}
      }
    }
    if(Vba[8].size() != 0){// Presintomático aislado a Leve aislado
      it = std::find(Vba[8].begin(), Vba[8].end(), agent);      ind = std::distance(Vba[8].begin(), it);
      if(ind < Vba[8].size()){mother_reaction(Vba[8], Vba[11], ind, bajos, 8, 11);}
    }
  }
}


/* Pre-sintomático a Grave. Alto */
void reaction6(std::vector<grupo> &Val, std::vector<grupo> &Vba, Crandom &ran, Workers *altos, Workers *bajos, std::vector<lognormal_d> &dist, int agentI){
  std::vector<int> aux;

  for(unsigned int i=0; i<Val[6].size(); i++){if(altos[Val[6][i]].tstate > TM){aux.push_back(Val[6][i]);}} //Presintomático
  for(unsigned int i=0; i<Val[7].size(); i++){if(altos[Val[7][i]].tstate > TM){aux.push_back(Val[7][i]);}} //Presintomático testeado
  for(unsigned int i=0; i<Val[8].size(); i++){if(altos[Val[8][i]].tstate > TM){aux.push_back(Val[8][i]);}} //Presintomático aislado

  if(aux.size() != 0){
    double value = (double)ran.r();
    unsigned int index = index_time(aux, altos, dist[4], value);
    unsigned int agent = aux[index];
    std::vector<int>::iterator it;
    unsigned int ind;
    if(Val[6].size() != 0){// Presintomático a Grave
      it = std::find(Val[6].begin(), Val[6].end(), agent);      ind = std::distance(Val[6].begin(), it);
      if(ind < Val[6].size()){
	mother_reaction(Val[6], Val[12], ind, altos, 6, 12);
	aux_main(1, Val, Vba, altos, bajos, ran, MyCons.phi1, MyCons.mu, true, 12, true);
      }
    }
    if(Val[7].size() != 0){// Presintomático testeado a Grave
      it = std::find(Val[7].begin(), Val[7].end(), agent);      ind = std::distance(Val[7].begin(), it);
      if(ind < Val[7].size()){
	mother_reaction(Val[7], Val[12], ind, altos, 7, 12);
	aux_main(1, Val, Vba, altos, bajos, ran, MyCons.phi1, MyCons.mu, true, 12, true);
      }
    }
    if(Val[8].size() != 0){// Presintomáticos aislado a Grave
      it = std::find(Val[8].begin(), Val[8].end(), agent);      ind = std::distance(Val[8].begin(), it);
      if(ind < Val[8].size()){mother_reaction(Val[8], Val[12], ind, altos, 8, 12);}
    }
  }
}


/* Pre-sintomático a Grave. Bajo */
void reaction7(std::vector<grupo> &Val, std::vector<grupo> &Vba, Crandom &ran, Workers *altos, Workers *bajos, std::vector<lognormal_d> &dist, int agentI){
  std::vector<int> aux;

  for(unsigned int i=0; i<Vba[6].size(); i++){if(bajos[Vba[6][i]].tstate > TM){aux.push_back(Vba[6][i]);}} //Presintomático
  for(unsigned int i=0; i<Vba[7].size(); i++){if(bajos[Vba[7][i]].tstate > TM){aux.push_back(Vba[7][i]);}} //Presintomático testeado
  for(unsigned int i=0; i<Vba[8].size(); i++){if(bajos[Vba[8][i]].tstate > TM){aux.push_back(Vba[8][i]);}} //Presitnomático aislado

  if(aux.size() != 0){
    double value = (double)ran.r();
    unsigned int index = index_time(aux, bajos, dist[5], value);
    unsigned int agent = aux[index];
    std::vector<int>::iterator it;
    unsigned int ind;
    if(Vba[6].size() != 0){// Presintomático a Grave
      it = std::find(Vba[6].begin(), Vba[6].end(), agent);      ind = std::distance(Vba[6].begin(), it);
      if(ind < Vba[6].size()){
	mother_reaction(Vba[6], Vba[12], ind, bajos, 6, 12);
	aux_main(1, Val, Vba, altos, bajos, ran, MyCons.mu, MyCons.chi, false, 12, true);
      }
    }
    if(Vba[7].size() != 0){
      it = std::find(Vba[7].begin(), Vba[7].end(), agent);      ind = std::distance(Vba[7].begin(), it);
      if(ind < Vba[7].size()){// Presintomático testeado a Grave
	mother_reaction(Vba[7], Vba[12], ind, bajos, 7, 12);
	aux_main(1, Val, Vba, altos, bajos, ran, MyCons.mu, MyCons.chi, false, 12, true);
      }
    }
    if(Vba[8].size() != 0){// Presintomático aislado a Grave
      it = std::find(Vba[8].begin(), Vba[8].end(), agent);      ind = std::distance(Vba[8].begin(), it);
      if(ind < Vba[8].size()){mother_reaction(Vba[8], Vba[12], ind, bajos, 8, 12);}
    }
  }
}


/* Pre-sintomático a Recuperado. Alto */
void reaction8(std::vector<grupo> &Val, std::vector<grupo> &Vba, Crandom &ran, Workers *altos, Workers *bajos, std::vector<lognormal_d> &dist, int agentI){
  std::vector<int> aux;

  for(unsigned int i=0; i<Val[6].size(); i++){if(altos[Val[6][i]].tstate > TM){aux.push_back(Val[6][i]);}} //Presintomático
  for(unsigned int i=0; i<Val[7].size(); i++){if(altos[Val[7][i]].tstate > TM){aux.push_back(Val[7][i]);}} //Presintomático testeado
  for(unsigned int i=0; i<Val[8].size(); i++){if(altos[Val[8][i]].tstate > TM){aux.push_back(Val[8][i]);}} //Presintomático aislado

  if(aux.size() != 0){
    double value = (double)ran.r();
    unsigned int index = index_time(aux, altos, dist[6], value);
    unsigned int agent = aux[index];
    std::vector<int>::iterator it;
    unsigned int ind;
    if(Val[6].size() != 0){// Presintomático a Recuperado no-detectado
      it = std::find(Val[6].begin(), Val[6].end(), agent);      ind = std::distance(Val[6].begin(), it);
      if(ind < Val[6].size()){mother_reaction(Val[6], Val[13], ind, altos, 6, 13);}
    }
    if(Val[7].size() != 0){// Presintomático testeado a Recuperado detectado o no-detectado
      it = std::find(Val[7].begin(), Val[7].end(), agent);      ind = std::distance(Val[7].begin(), it);
      if(ind < Val[7].size()){
        if(ran.r() < MyCons.xi){//Presintomático testeado a Recuperado detectado
          mother_reaction(Val[7], Val[15], ind, altos, 7, 15);
 	  aux_main(1, Val, Vba, altos, bajos, ran, MyCons.phi1, MyCons.mu, true, 15, true);
        }
        else{//Presintomático testeado a Recuperado detectado
          mother_reaction(Val[7], Val[13], ind, altos, 7, 13);
        }
      }
    }
    if(Val[8].size() != 0){// Presintomático aislado a Recuperado detectado
      it = std::find(Val[8].begin(), Val[8].end(), agent);      ind = std::distance(Val[8].begin(), it);
      if(ind < Val[8].size()){mother_reaction(Val[8], Val[15], ind, altos, 8, 15);}
    }
  }
}


/* Pre-sintomático a Recuperado. Bajo */
void reaction9(std::vector<grupo> &Val, std::vector<grupo> &Vba, Crandom &ran, Workers *altos, Workers *bajos, std::vector<lognormal_d> &dist, int agentI){
  std::vector<int> aux;

  for(unsigned int i=0; i<Vba[6].size(); i++){if(bajos[Vba[6][i]].tstate > TM){aux.push_back(Vba[6][i]);}} //Presintomático
  for(unsigned int i=0; i<Vba[7].size(); i++){if(bajos[Vba[7][i]].tstate > TM){aux.push_back(Vba[7][i]);}} //Presintomático testeado
  for(unsigned int i=0; i<Vba[8].size(); i++){if(bajos[Vba[8][i]].tstate > TM){aux.push_back(Vba[8][i]);}} //Presintomático aislado

  if(aux.size() != 0){
    double value = (double)ran.r();
    unsigned int index = index_time(aux, bajos, dist[7], value);
    unsigned int agent = aux[index];
    std::vector<int>::iterator it;
    unsigned int ind;
    if(Vba[6].size() != 0){// Presintomático a Recuperado no-detectado
      it = std::find(Vba[6].begin(), Vba[6].end(), agent);      ind = std::distance(Vba[6].begin(), it);
      if(ind < Vba[6].size()){mother_reaction(Vba[6], Vba[13], ind, bajos, 6, 13);}
    }
    if(Vba[7].size() != 0){// Presintomático testeado a Recuperado detectado o no-detectado
      it = std::find(Vba[7].begin(), Vba[7].end(), agent);      ind = std::distance(Vba[7].begin(), it);
      if(ind < Vba[7].size()){
        if(ran.r() < MyCons.xi){// Presintomático testeado a Recuperado detectado
          mother_reaction(Vba[7], Vba[15], ind, bajos, 7, 15);
 	  aux_main(1, Val, Vba, altos, bajos, ran, MyCons.mu, MyCons.chi, false, 15, true);
        }
        else{// Presintomático testeado a Recuperado no-detectado
          mother_reaction(Vba[7], Vba[13], ind, bajos, 7, 13);
        }
      }
    }
    if(Vba[8].size() != 0){// Presintomático aislado a Recuperado detectado
      it = std::find(Vba[8].begin(), Vba[8].end(), agent);      ind = std::distance(Vba[8].begin(), it);
      if(ind < Vba[8].size()){mother_reaction(Vba[8], Vba[15], ind, bajos, 8, 15);}
    }
  }
}


/* Leve a Recuperado. Alto */
void reaction10(std::vector<grupo> &Val, std::vector<grupo> &Vba, Crandom &ran, Workers *altos, Workers *bajos, std::vector<lognormal_d> &dist, int agentI){
  std::vector<int> aux;

  for(unsigned int i=0; i<Val[9].size(); i++){if(altos[Val[9][i]].tstate > TM){aux.push_back(Val[9][i]);}} //Leve
  for(unsigned int i=0; i<Val[10].size(); i++){if(altos[Val[10][i]].tstate > TM){aux.push_back(Val[10][i]);}} //Leve testeado
  for(unsigned int i=0; i<Val[11].size(); i++){if(altos[Val[11][i]].tstate > TM){aux.push_back(Val[11][i]);}} //Leve aislado

  if(aux.size() != 0){
    double value = (double)ran.r();
    unsigned int index = index_time(aux, altos, dist[8], value);
    unsigned int agent = aux[index];
    std::vector<int>::iterator it;
    unsigned int ind;
    if(Val[9].size() != 0){// Leve a Recuperado no-detectado
      it = std::find(Val[9].begin(), Val[9].end(), agent);      ind = std::distance(Val[9].begin(), it);
      if(ind < Val[9].size()){mother_reaction(Val[9], Val[13], ind, altos, 9, 13);}
    }
    if(Val[10].size() != 0){// Leve testeado a Recuperado detectado o no-detectado
      it = std::find(Val[10].begin(), Val[10].end(), agent);      ind = std::distance(Val[10].begin(), it);
      if(ind < Val[10].size()){
        if(ran.r() < MyCons.xi){// Leve testeado a Recuperado detectado
          mother_reaction(Val[10], Val[15], ind, altos, 10, 15);
 	  aux_main(1, Val, Vba, altos, bajos, ran, MyCons.phi1, MyCons.mu, true, 15, false);
        }
        else{// Leve testeado a Recuperado no-detectado
          mother_reaction(Val[10], Val[13], ind, altos, 10, 13);
        }
      }
    }
    if(Val[11].size() != 0){// Leve aislado a Recuperado detectado
      it = std::find(Val[11].begin(), Val[11].end(), agent);      ind = std::distance(Val[11].begin(), it);
      if(ind < Val[11].size()){mother_reaction(Val[11], Val[15], ind, altos, 11, 15);}
    }
  }
}


/* Leve a Recuperado. Bajo */
void reaction11(std::vector<grupo> &Val, std::vector<grupo> &Vba, Crandom &ran, Workers *altos, Workers *bajos, std::vector<lognormal_d> &dist, int agentI){
  std::vector<int> aux;

  for(unsigned int i=0; i<Vba[9].size(); i++){if(bajos[Vba[9][i]].tstate > TM){aux.push_back(Vba[9][i]);}} //Leve
  for(unsigned int i=0; i<Vba[10].size(); i++){if(bajos[Vba[10][i]].tstate > TM){aux.push_back(Vba[10][i]);}} //Leve testeado
  for(unsigned int i=0; i<Vba[11].size(); i++){if(bajos[Vba[11][i]].tstate > TM){aux.push_back(Vba[11][i]);}} //Leve aislado

  if(aux.size() != 0){
    double value = (double)ran.r();
    unsigned int index = index_time(aux, bajos, dist[9], value);
    unsigned int agent = aux[index];
    std::vector<int>::iterator it;
    unsigned int ind;
    if(Vba[9].size() != 0){// Leve a Recuperado no-detectado
      it = std::find(Vba[9].begin(), Vba[9].end(), agent);      ind = std::distance(Vba[9].begin(), it);
      if(ind < Vba[9].size()){mother_reaction(Vba[9], Vba[13], ind, bajos, 9, 13);}
    }
    if(Vba[10].size() != 0){// Leve testeado a Recuperado detectado o no-detectado
      it = std::find(Vba[10].begin(), Vba[10].end(), agent);      ind = std::distance(Vba[10].begin(), it);
      if(ind < Vba[10].size()){
        if(ran.r() < MyCons.xi){// Leve testeado a Recuperado detectado
          mother_reaction(Vba[10], Vba[15], ind, bajos, 10, 15);
 	  aux_main(1, Val, Vba, altos, bajos, ran, MyCons.mu, MyCons.chi, false, 15, false);
        }
        else{// Leve testeado a Recuperado no-detectado
          mother_reaction(Vba[10], Vba[13], ind, bajos, 10, 13);
        }
      }
    }
    if(Vba[11].size() != 0){// Leve aislado a Recuperado detectado
      it = std::find(Vba[11].begin(), Vba[11].end(), agent);      ind = std::distance(Vba[11].begin(), it);
      if(ind < Vba[11].size()){mother_reaction(Vba[11], Vba[15], ind, bajos, 11, 15);}
    }
  }
}


/* Grave a Recuperado. Alto */
void reaction12(std::vector<grupo> &Val, std::vector<grupo> &Vba, Crandom &ran, Workers *altos, Workers *bajos, std::vector<lognormal_d> &dist, int agentI){
  std::vector<int> aux;

  for(unsigned int i=0; i<Val[12].size(); i++){if(altos[Val[12][i]].tstate > TM){aux.push_back(Val[12][i]);}} //Grave

  if(aux.size() != 0){// Grave a Recuperado detectado
    unsigned int index = index_time(aux, altos, dist[10], (double)ran.r());
    unsigned int agent = aux[index];
    std::vector<int>::iterator it = std::find(Val[12].begin(), Val[12].end(), agent);
    unsigned int ind = std::distance(Val[12].begin(), it);
    if(ind < Val[12].size()){mother_reaction(Val[12], Val[15], ind, altos, 12, 15);}
  }
}


/* Grave a Recuperado. Bajo */
void reaction13(std::vector<grupo> &Val, std::vector<grupo> &Vba, Crandom &ran, Workers *altos, Workers *bajos, std::vector<lognormal_d> &dist, int agentI){
  std::vector<int> aux;

  for(unsigned int i=0; i<Vba[12].size(); i++){if(bajos[Vba[12][i]].tstate > TM){aux.push_back(Vba[12][i]);}} //Grave

  if(aux.size() != 0){// Grave a Recuperado detectado
    unsigned int index = index_time(aux, bajos, dist[11], (double)ran.r());
    unsigned int agent = aux[index];
    std::vector<int>::iterator it = std::find(Vba[12].begin(), Vba[12].end(), agent);
    unsigned int ind = std::distance(Vba[12].begin(), it);
    if(ind < Vba[12].size()){mother_reaction(Vba[12], Vba[15], ind, bajos, 12, 15);}
  }
}
