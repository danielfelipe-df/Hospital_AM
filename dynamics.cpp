#include <iostream>
#include <math.h>
#include <algorithm>
#include "dynamics.h"

std::vector<double> contagio(double Sa, double Sb, double STa, double STb, double Ea, double Eb, double ETa, double ETb, double Pa, double Pb, double PTa, double PTb, double PTAa, double PTAb, double La, double Lb, double LTa, double LTb, double LTAa, double LTAb, double IAa, double IAb, double prev, Crandom &ran, double t, double* tj){

  //Número de proponsidades
  const unsigned int n = 14;

  //Creo el arreglo de las propensidades
  double As[n];

  //Propensidades de exponerse
  As[0] = beta*(Sa+STa)*(phi1*(Pa+PTa+La+LTa)/(double)Na + mu*(Pb+PTb+Lb+LTb)/(double)Nb + (1-alpha)*phi1*(IAa+PTAa+LTAa)/(double)Na + (1-alpha)*mu*(IAb+PTAb+LTAb)/(double)Nb + eta*prev);
  As[1] = beta*(Sb+STb)*(mu*(Pa+PTa+La+LTa)/(double)Na + chi*(Pb+PTb+Lb+LTb)/(double)Nb + (1-alpha)*mu*(IAa+PTAa+LTAa)/(double)Na + (1-alpha)*chi*(IAb+PTAb+LTAb)/(double)Nb);

  //Propensidades de ser presintomático
  As[2] = USDe*(Ea+ETa);
  As[3] = USDe*(Eb+ETb);

  //Propensidades de ser leve
  As[4] = USDpl*(1-kappa)*psi*(Pa+PTa+PTAa);

  As[5] = USDpl*(1-kappa)*psi*(Pb+PTb+PTAb);

  //Propensidad de ser infeccioso grave y aislarse inmediatamente
  As[6] = USDpg*(1-kappa)*(1-psi)*(Pa+PTa+PTAa);
  As[7] = USDpg*(1-kappa)*(1-psi)*(Pb+PTb+PTAb);

  //Propensidad de recuperarse siendo pre-sintomático
  As[8] = USDplil*kappa*(Pa+PTa+PTAa);
  As[9] = USDplil*kappa*(Pb+PTb+PTAb);

  //Propensidad de recuperarse siendo leve
  As[10] = USDil*(La+LTa+LTAa);
  As[11] = USDil*(Lb+LTb+LTAb);

  //Propensidad de recuperarse siendo infeccioso grave
  As[12] = USDig*IAa;
  As[13] = USDig*IAb;

  //Las propensidades que sean cero las pongo con un mínimo
  //for(unsigned int i=2; i<n; i++){if(As[i] <= 0.0){As[i] = 1e-12;}}

  //Creo las distribuciones de cada propensidad
  std::vector<lognormal_d> dist;
  //for(unsigned int i=2; i<n; i++){lognormal_d my_dist(1.5, 1.0/As[i]);    dist.push_back(my_dist);}

  //Hallo el tiempo en el que va a pasar la siguiente reacción con el método NMGA
  //double prom = 230.0, sigma = 55.0, B = beta*Sa*eta*prev;
  double prom = 165.47, sigma = 29.24, B = beta*Sa*eta*prev;
  double tau = 0, index = 0;
  tau = biseccion(As, prom, sigma, t, B, std::log(ran.r()), tj, n, dist);

  if(tau < 1e6){//Si el tiempo en el que pasa la reacción es menor al máximo (1e6), entonces es porque la reacción si sucede
    //Hallo el vector que me guarda la propensidad acumulada en orden
    double cumulative[n];
    cumulative[0] = (As[0]-B) + B*std::exp(-(prom-t)*(prom-t)/(2*sigma*sigma));
    cumulative[1] = cumulative[0] + As[1];
    for(unsigned int i=2; i<n; i++){cumulative[i] = cumulative[i-1] + As[i];}//(pdf(dist[i-2], tj[i]+tau)/cdf(complement(dist[i-2], tj[i]+tau)));}

    //Escojo la reacción a escoger
    double lim = (double)(ran.r()*cumulative[n-1]);
    for(index = 0; index<(n-1); index++){
      if(lim < cumulative[(int)index]){break;}
    }
    if(lim >= cumulative[n-2]){index = n-1;}
  }

  //Le sumo el tau a todos los tiempos tj
  for(unsigned int i=0; i<n; i++){tj[i] += tau;}

  //Reinicio el tiempo tj de la reacción que se escogió
  tj[(int)index] = 0.0;

  //Creo el vector resultados
  std::vector<double> result(2);
  result[0] = tau; //Tiempo en el que sucede la reacción
  result[1] = index; //Número de la reacción que sucede

  return result;
}


void mother_reaction(grupo &Out, grupo &In, int index, trabajadores *family, int typeout, int typein){
  int agent = Out[index];
  Out.erase(Out.begin() + index);
  In.push_back(agent);
  family[agent].change(typein, typeout);
  family[agent].tstate = 0.0;
}


void tested_reaction(grupo &Out, grupo &In, int index, trabajadores *family, int typeout, int typein, double delta){
  int agent = Out[index];
  Out.erase(Out.begin() + index);
  In.push_back(agent);
  family[agent].change(typein, typeout);
  family[agent].time = 0.0;
  family[agent].tmax = delta;
}

void massive_reaction(grupo &Sa, grupo &Sb, grupo &SMa, grupo &SMb, grupo &Ea, grupo &Eb, grupo &EMa, grupo &EMb, grupo &Pa, grupo &Pb, grupo &PTa, grupo &PTb, grupo &La, grupo &Lb, grupo &LTa, grupo &LTb, grupo &RIa, grupo &RIb, grupo &RMa, grupo &RMb, Crandom &ran, trabajadores *altos, trabajadores *bajos){
  unsigned int M = Sa.size() + Sb.size() + Ea.size() + Eb.size() + Pa.size() + Pb.size() + La.size() + Lb.size() + RIa.size() + RIb.size();
  unsigned int num = (int)(ran.r()*M);
  
  if(num < Sa.size()){
    tested_reaction(Sa, SMa, (int)(ran.r()*Sa.size()), altos, 0, 1, Tt);
  }
  else if(num < Sa.size() + Sb.size()){
    tested_reaction(Sb, SMb, (int)(ran.r()*Sb.size()), bajos, 0, 1, Tt);
  }
  else if(num < Sa.size() + Sb.size() + Ea.size()){
    tested_reaction(Ea, EMa, (int)(ran.r()*Ea.size()), altos, 2, 3, Tt);
  }
  else if(num < Sa.size() + Sb.size() + Ea.size() + Eb.size()){
    tested_reaction(Eb, EMb, (int)(ran.r()*Eb.size()), bajos, 2, 3, Tt);
  }
  else if(num < Sa.size() + Sb.size() + Ea.size() + Eb.size() + Pa.size()){
    tested_reaction(Pa, PTa, (int)(ran.r()*Pa.size()), altos, 4, 5, Tt);
  }
  else if(num < Sa.size() + Sb.size() + Ea.size() + Eb.size() + Pa.size() + Pb.size()){
    tested_reaction(Pb, PTb, (int)(ran.r()*Pb.size()), bajos, 4, 5, Tt);
  }
  else if(num < Sa.size() + Sb.size() + Ea.size() + Eb.size() + Pa.size() + Pb.size() + La.size()){
    tested_reaction(La, LTa, (int)(ran.r()*La.size()), altos, 7, 8, Tt);
  }
  else if(num < Sa.size() + Sb.size() + Ea.size() + Eb.size() + Pa.size() + Pb.size() + La.size() + Lb.size()){
    tested_reaction(Lb, LTb, (int)(ran.r()*Lb.size()), bajos, 7, 8, Tt);
  }
  else if(num < Sa.size() + Sb.size() + Ea.size() + Eb.size() + Pa.size() + Pb.size() + La.size() + Lb.size() + RIa.size()){
    tested_reaction(RIa, RMa, (int)(ran.r()*RIa.size()), altos, 11, 12, Tt);
  }
  else{
    tested_reaction(RIb, RMb, (int)(ran.r()*RIb.size()), bajos, 11, 12, Tt);
  }
}


int tested_isolated_inf(grupo &T, grupo &TA, grupo &G, trabajadores *family, double time, int typeout, int typein1, int typein2, Crandom &ran){
  int index, contador = 0;
  for(unsigned int i=0; i<T.size(); i++){
    index = T[i];      family[index].time += time;
    if(family[index].time > family[index].tmax){
      if(ran.r() < xi){family[index].change(typein1, typeout);	  TA.push_back(index);	contador++;}
      else{family[index].change(typein2, typeout);	  G.push_back(index);}
      T.erase( T.begin() + i);
      family[index].time = 0.0;
      family[index].tmax = 0.0;
      i--;
    }
  }
  return contador;
}


void tested_massive(grupo &T, grupo &G, trabajadores *family, double time, int typeout, int typein){
  for(unsigned int i=0; i<T.size(); i++){
    family[T[i]].time += time;
    if(family[T[i]].time > family[T[i]].tmax){tested_reaction(T, G, i, family, typeout, typein, 0.0);      i--;}
  }
}


double biseccion(double* A, double prom, double sigma, double t, double B, double ranr, double* tj, int n, std::vector<lognormal_d> &dist){
  double m,fa,fm;
  double lim = 1e3, min = 0.0;
  double a = min, b = lim, eps = 1e-7, pmax = 100, p=0;
  fa = phi(A, tj, n, prom, sigma, B, a, t, dist) - ranr;

  while(b-a>eps && p<pmax){
    m = (a+b)/2;
    fm = phi(A, tj, n, prom, sigma, B, m, t, dist) - ranr;
    if(fa*fm<0){b = m;}
    else{a = m; fa = fm;}
    p++;
  }

  if(p==pmax || m > lim-(1e-4)){return 1e6;}
  else{return (a+b)/2;}
}


double function(double A, double prom, double sigma, double t, double tau){
  double A2 = std::sqrt(M_PI/2)*A*sigma, USsigma2 = 1.0/(std::sqrt(2)*sigma);
  return A2*(erf((prom-t)*USsigma2) - erf((prom-t-tau)*USsigma2));
}


double phi(double* A, double* tj, unsigned int n, double prom, double sigma, double B, double deltat, double t, std::vector<lognormal_d> &dist){
  double psi_num, psi_den;

  psi_num = -(A[0]-B)*(tj[0]+deltat) - function(B, prom, sigma, t, tj[0]+deltat);
  psi_den = -(A[0]-B)*tj[0] - function(B, prom, sigma, t, tj[0]);

  psi_num += -A[1]*(tj[1]+deltat);
  psi_den += -A[1]*tj[1];

  for(unsigned int i=2; i<n; i++){
    //psi_num += std::log(cdf(complement(dist[i-2], tj[i]+deltat)));
    //psi_den += std::log(cdf(complement(dist[i-2], tj[i])));
    psi_num += -A[i]*(tj[i] + deltat);
    psi_den += -A[i]*tj[i];
  }

  return psi_num-psi_den;
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


void update_times(grupo &G, trabajadores *family, double time){
  for(unsigned int i=0; i<G.size(); i++){family[G[i]].tstate += time;}
}


void update_massive(grupo &G, trabajadores *family, double time){
  for(unsigned int i=0; i<G.size(); i++){family[G[i]].time += time;}
}


void move_massive(grupo &T, grupo &G, trabajadores *family, unsigned int typeout, unsigned int typein){
  for(unsigned int i=0; i<T.size(); i++){
    if(family[T[i]].time > family[T[i]].tmax){tested_reaction(T, G, i, family, typeout, typein, 0.0);      i--;}
  }
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
    if(num2 < num[0]){
      return selection_infectious(Pa, PTa, La, LTa, ran, index, altos);
    }
    else if(num2 < num[0] + num[1]){
      return selection_infectious(Pb, PTb, Lb, LTb, ran, index, bajos) + Na;
    }
    else if(num2 < num[0] + num[1] + num[2]){
      return selection_infectious(IAa, PTAa, LTAa, aux, ran, index, altos);
    }
    else{
      return selection_infectious(IAb, PTAb, LTAb, aux, ran, index, bajos) + Na;
    }
  }
  else{return -1;}
}


int selection_infectious(grupo &Ga, grupo &Gb, grupo &Gc, grupo &Gd, Crandom &ran, int index, trabajadores *family){
  double num = ran.r()*(Ga.size() + Gb.size() + Gc.size() + Gd.size());
  int ind, agent;
  if(num < Ga.size()){
    ind = (int)(ran.r()*Ga.size());    agent = Ga[ind];    family[agent].my_inf.push_back(index);
  }
  else if(num < Ga.size() + Gb.size()){
    ind = (int)(ran.r()*Gb.size());    agent = Gb[ind];    family[agent].my_inf.push_back(index);
  }
  else if(num < Ga.size() + Gb.size() + Gc.size()){
    ind = (int)(ran.r()*Gc.size());    agent = Gc[ind];    family[agent].my_inf.push_back(index);
  }
  else{
    ind = (int)(ran.r()*Gd.size());    agent = Gd[ind];    family[agent].my_inf.push_back(index);
  }
  return agent;
}


void main_trace(grupo &Sa, grupo &Sb, grupo &STa, grupo &STb, grupo &Ea, grupo &Eb, grupo &ETa, grupo &ETb, grupo &Pa, grupo &Pb, grupo &PTa, grupo &PTb, grupo &PTAa, grupo &PTAb, grupo &La, grupo &Lb, grupo &LTa, grupo &LTb, grupo &LTAa, grupo &LTAb, grupo &IAa, grupo &IAb, grupo &RTa, grupo&RTb, grupo &RIa, grupo &RIb, trabajadores *altos, trabajadores *bajos, double time, Crandom &ran){
  int num;
  num = tested_isolated_inf(PTa, PTAa, Pa, altos, time, 5, 6, 4, ran);
  aux_main(num, PTAa, Sa, Sb, STa, STb, Ea, Eb, ETa, ETb, Pa, Pb, PTa, PTb, La, Lb, LTa, LTb, RIa, RIb, RTa, RTb, altos, bajos, ran, phi1, mu, true);

  num = tested_isolated_inf(PTb, PTAb, Pb, bajos, time, 5, 6, 4, ran);
  aux_main(num, PTAb, Sa, Sb, STa, STb, Ea, Eb, ETa, ETb, Pa, Pb, PTa, PTb, La, Lb, LTa, LTb, RIa, RIb, RTa, RTb, altos, bajos, ran, mu, chi, false);

  num = tested_isolated_inf(LTa, LTAa, La, altos, time, 8, 9, 7, ran);
  aux_main(num, LTAa, Sa, Sb, STa, STb, Ea, Eb, ETa, ETb, Pa, Pb, PTa, PTb, La, Lb, LTa, LTb, RIa, RIb, RTa, RTb, altos, bajos, ran, phi1, mu, true);

  num = tested_isolated_inf(LTb, LTAb, Lb, bajos, time, 8, 9, 7, ran);
  aux_main(num, LTAb, Sa, Sb, STa, STb, Ea, Eb, ETa, ETb, Pa, Pb, PTa, PTb, La, Lb, LTa, LTb, RIa, RIb, RTa, RTb, altos, bajos, ran, mu, chi, false);
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
