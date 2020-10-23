#include <iostream>
#include <math.h>
#include <algorithm>
#include "dynamics.h"

std::vector<double> contagio(double Sa, double Sb, double STa, double STb, double SMa, double SMb, double Ea, double Eb, double ETa, double ETb, double EMa, double EMb, double Pa, double Pb, double PTa, double PTb, double PTAa, double PTAb, double La, double Lb, double LTa, double LTb, double LTAa, double LTAb, double IAa, double IAb, double Na, double Nb, double prev, Crandom &ran, double t, double* tj){

  //Número de proponsidades
  const unsigned int n = 14;

  //Creo el arreglo de las propensidades
  double As[n];

  //Propensidades de exponerse
  As[0] = beta*(Sa+STa+SMa)*(phi1*(Pa+PTa+La+LTa)/Na + mu*(Pb+PTb+Lb+LTb)/Nb + (1-alpha)*phi1*(IAa+PTAa+LTAa)/Na + (1-alpha)*mu*(IAb+PTAb+LTAb)/Nb + eta*prev);
  As[1] = beta*(Sb+STb+SMb)*(mu*(Pa+PTa+La+LTa)/Na + chi*(Pb+PTb+Lb+LTb)/Nb + (1-alpha)*mu*(IAa+PTAa+LTAa)/Na + (1-alpha)*chi*(IAb+PTAb+LTAb)/Nb);

  //Propensidades de ser presintomático
  As[2] = USDe*(Ea+ETa+EMa);
  As[3] = USDe*(Eb+ETb+EMb);

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
  std::vector<weib_d> dist;
  //for(unsigned int i=2; i<n; i++){weib_d my_dist(1.5, 1.0/As[i]);    dist.push_back(my_dist);}

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
    tested_reaction(Sa, SMa, (int)(ran.r()*Sa.size()), altos, 5, 6, Tt);
  }
  else if(num < Sa.size() + Sb.size()){
    tested_reaction(Sb, SMb, (int)(ran.r()*Sb.size()), bajos, 5, 6, Tt);
  }
  else if(num < Sa.size() + Sb.size() + Ea.size()){
    tested_reaction(Ea, EMa, (int)(ran.r()*Ea.size()), altos, 5, 6, Tt);
  }
  else if(num < Sa.size() + Sb.size() + Ea.size() + Eb.size()){
    tested_reaction(Eb, EMb, (int)(ran.r()*Eb.size()), bajos, 5, 6, Tt);
  }
  else if(num < Sa.size() + Sb.size() + Ea.size() + Eb.size() + Pa.size()){
    tested_reaction(Pa, PTa, (int)(ran.r()*Pa.size()), altos, 5, 6, Tt);
  }
  else if(num < Sa.size() + Sb.size() + Ea.size() + Eb.size() + Pa.size() + Pb.size()){
    tested_reaction(Pb, PTb, (int)(ran.r()*Pb.size()), bajos, 5, 6, Tt);
  }
  else if(num < Sa.size() + Sb.size() + Ea.size() + Eb.size() + Pa.size() + Pb.size() + La.size()){
    tested_reaction(La, LTa, (int)(ran.r()*La.size()), altos, 5, 6, Tt);
  }
  else if(num < Sa.size() + Sb.size() + Ea.size() + Eb.size() + Pa.size() + Pb.size() + La.size() + Lb.size()){
    tested_reaction(Lb, LTb, (int)(ran.r()*Lb.size()), bajos, 5, 6, Tt);
  }
  else if(num < Sa.size() + Sb.size() + Ea.size() + Eb.size() + Pa.size() + Pb.size() + La.size() + Lb.size() + RIa.size()){
    tested_reaction(RIa, RMa, (int)(ran.r()*RIa.size()), altos, 5, 6, Tt);
  }
  else{
    tested_reaction(RIb, RMb, (int)(ran.r()*RIb.size()), bajos, 5, 6, Tt);
  }
}


void continue_reaction(grupo &L, grupo &LT, trabajadores *family, Crandom &ran){
  if(ran.r() < iota){tested_reaction(L, LT, L.size()-1, family, 5, 7, 0.0);}
}


void tested_isolated_inf(grupo &T, grupo &TA, grupo &G, trabajadores *family, double time, int typeout, int typein1, int typein2, Crandom &ran){
  int index;
  for(unsigned int i=0; i<T.size(); i++){
    index = T[i];      family[index].time += time;
    if(family[index].time > family[index].tmax){
      if(ran.r() < xi){family[index].change(typein1, typeout);	  TA.push_back(index);}
      else{family[index].change(typein2, typeout);	  G.push_back(index);}
      T.erase( T.begin() + i);
      family[index].time = 0.0;
      family[index].tmax = 0.0;
      i--;
    }
  }
}


void tested_isolated(grupo &T, grupo &G, trabajadores *family, double time, int typeout, int typein, Crandom &ran){
  int index;
  for(unsigned int i=0; i<T.size(); i++){
    index = T[i];      family[index].time += time;
    if(family[index].time > family[index].tmax){
      family[index].change(typein, typeout);
      family[index].time = 0.0;
      family[index].tmax = 0.0;
      T.erase( T.begin() + i);
      G.push_back(index);
      i--;
    }
  }
}


double biseccion(double* A, double prom, double sigma, double t, double B, double ranr, double* tj, int n, std::vector<weib_d> &dist){
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


double phi(double* A, double* tj, unsigned int n, double prom, double sigma, double B, double deltat, double t, std::vector<weib_d> &dist){
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


int index_time(grupo &Out, trabajadores *family, weib_d &dist, double value){
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


void move_massive(grupo &M, grupo &T, grupo &G, trabajadores *family){
  for(unsigned int i=0; i<M.size(); i++){
    if(family[M[i]].time > family[M[i]].tmax){G.push_back(M[i]);}
    else{T.push_back(M[i]);}
  }
  M.clear();
}


void reaction0(grupo &Sa, grupo &Sb, grupo &STa, grupo &STb, grupo &SMa, grupo &SMb, grupo &Ea, grupo &Eb, grupo &ETa, grupo &ETb, grupo &EMa, grupo &EMb, grupo &Pa, grupo &Pb, grupo &PTa, grupo &PTb, grupo &PTAa, grupo &PTAb, grupo &La, grupo &Lb, grupo &LTa, grupo &LTb, grupo &LTAa, grupo &LTAb, grupo &IAa, grupo &IAb, grupo &RTa, grupo&RTb, grupo &RIa, grupo &RIb, grupo &RMa, grupo &RMb, grupo &RAa, grupo &RAb, Crandom &ran, trabajadores *altos, trabajadores *bajos, std::vector<weib_d> &dist){
  unsigned int index = (int)(ran.r()*(Sa.size() + STa.size() + SMa.size()));
  if(index < Sa.size()){mother_reaction(Sa, Ea, index, altos, 0, 1);}
  else if(index < Sa.size() + STa.size()){mother_reaction(STa, ETa, index-Sa.size(), altos, 0, 1);}
  else{mother_reaction(SMa, EMa, index-Sa.size()-STa.size(), altos, 0, 1);}
}


void reaction1(grupo &Sa, grupo &Sb, grupo &STa, grupo &STb, grupo &SMa, grupo &SMb, grupo &Ea, grupo &Eb, grupo &ETa, grupo &ETb, grupo &EMa, grupo &EMb, grupo &Pa, grupo &Pb, grupo &PTa, grupo &PTb, grupo &PTAa, grupo &PTAb, grupo &La, grupo &Lb, grupo &LTa, grupo &LTb, grupo &LTAa, grupo &LTAb, grupo &IAa, grupo &IAb, grupo &RTa, grupo&RTb, grupo &RIa, grupo &RIb, grupo &RMa, grupo &RMb, grupo &RAa, grupo &RAb, Crandom &ran, trabajadores *altos, trabajadores *bajos, std::vector<weib_d> &dist){
  unsigned int index = (int)(ran.r()*(Sb.size() + STb.size() + SMb.size()));
  if(index < Sb.size()){mother_reaction(Sb, Eb, index, bajos, 0, 1);}
  else if(index < Sb.size() + STb.size()){mother_reaction(STb, ETb, index-Sb.size(), bajos, 0, 1);}
  else{mother_reaction(SMb, EMb, index-Sb.size()-STb.size(), bajos, 0, 1);}
}


void reaction2(grupo &Sa, grupo &Sb, grupo &STa, grupo &STb, grupo &SMa, grupo &SMb, grupo &Ea, grupo &Eb, grupo &ETa, grupo &ETb, grupo &EMa, grupo &EMb, grupo &Pa, grupo &Pb, grupo &PTa, grupo &PTb, grupo &PTAa, grupo &PTAb, grupo &La, grupo &Lb, grupo &LTa, grupo &LTb, grupo &LTAa, grupo &LTAb, grupo &IAa, grupo &IAb, grupo &RTa, grupo&RTb, grupo &RIa, grupo &RIb, grupo &RMa, grupo &RMb, grupo &RAa, grupo &RAb, Crandom &ran, trabajadores *altos, trabajadores *bajos, std::vector<weib_d> &dist){
  std::vector<int> aux;

  for(unsigned int i=0; i<Ea.size(); i++){if(altos[Ea[i]].tstate > TM){aux.push_back(Ea[i]);}}
  for(unsigned int i=0; i<ETa.size(); i++){if(altos[ETa[i]].tstate > TM){aux.push_back(ETa[i]);}}
  for(unsigned int i=0; i<EMa.size(); i++){if(altos[EMa[i]].tstate > TM){aux.push_back(EMa[i]);}}

  if(aux.size() != 0){
    double value = (double)ran.r();
    unsigned int index = index_time(aux, altos, dist[0], value);
    unsigned int agent = aux[index];
    std::vector<int>::iterator it;
    unsigned int ind;
    if(Ea.size() != 0){
      it = std::find(Ea.begin(), Ea.end(), agent);      ind = std::distance(Ea.begin(), it);
      if(ind < Ea.size()){mother_reaction(Ea, Pa, ind, altos, 1, 2);}
    }
    if(ETa.size() != 0){
      it = std::find(ETa.begin(), ETa.end(), agent);      ind = std::distance(ETa.begin(), it);
      if(ind < ETa.size()){mother_reaction(ETa, Pa, ind, altos, 1, 2);}
    }
    if(EMa.size() != 0){
      it = std::find(EMa.begin(), EMa.end(), agent);      ind = std::distance(EMa.begin(), it);
      if(ind < EMa.size()){mother_reaction(EMa, Pa, ind, altos, 1, 2);}
    }
  }
}


void reaction3(grupo &Sa, grupo &Sb, grupo &STa, grupo &STb, grupo &SMa, grupo &SMb, grupo &Ea, grupo &Eb, grupo &ETa, grupo &ETb, grupo &EMa, grupo &EMb, grupo &Pa, grupo &Pb, grupo &PTa, grupo &PTb, grupo &PTAa, grupo &PTAb, grupo &La, grupo &Lb, grupo &LTa, grupo &LTb, grupo &LTAa, grupo &LTAb, grupo &IAa, grupo &IAb, grupo &RTa, grupo&RTb, grupo &RIa, grupo &RIb, grupo &RMa, grupo &RMb, grupo &RAa, grupo &RAb, Crandom &ran, trabajadores *altos, trabajadores *bajos, std::vector<weib_d> &dist){
  std::vector<int> aux;

  for(unsigned int i=0; i<Eb.size(); i++){if(bajos[Eb[i]].tstate > TM){aux.push_back(Eb[i]);}}
  for(unsigned int i=0; i<ETb.size(); i++){if(bajos[ETb[i]].tstate > TM){aux.push_back(ETb[i]);}}
  for(unsigned int i=0; i<EMb.size(); i++){if(bajos[EMb[i]].tstate > TM){aux.push_back(EMb[i]);}}

  if(aux.size() != 0){
    double value = (double)ran.r();
    unsigned int index = index_time(aux, bajos, dist[1], value);
    unsigned int agent = aux[index];
    std::vector<int>::iterator it;
    unsigned int ind;
    if(Eb.size() != 0){
      it = std::find(Eb.begin(), Eb.end(), agent);      ind = std::distance(Eb.begin(), it);
      if(ind < Eb.size()){mother_reaction(Eb, Pb, ind, bajos, 1, 2);}
    }
    if(ETb.size() != 0){
      it = std::find(ETb.begin(), ETb.end(), agent);      ind = std::distance(ETb.begin(), it);
      if(ind < ETb.size()){mother_reaction(ETb, Pb, ind, bajos, 1, 2);}
    }
    if(EMb.size() != 0){
      it = std::find(EMb.begin(), EMb.end(), agent);      ind = std::distance(EMb.begin(), it);
      if(ind < EMb.size()){mother_reaction(EMb, Pb, ind, bajos, 1, 2);}
    }
  }
}


void reaction4(grupo &Sa, grupo &Sb, grupo &STa, grupo &STb, grupo &SMa, grupo &SMb, grupo &Ea, grupo &Eb, grupo &ETa, grupo &ETb, grupo &EMa, grupo &EMb, grupo &Pa, grupo &Pb, grupo &PTa, grupo &PTb, grupo &PTAa, grupo &PTAb, grupo &La, grupo &Lb, grupo &LTa, grupo &LTb, grupo &LTAa, grupo &LTAb, grupo &IAa, grupo &IAb, grupo &RTa, grupo&RTb, grupo &RIa, grupo &RIb, grupo &RMa, grupo &RMb, grupo &RAa, grupo &RAb, Crandom &ran, trabajadores *altos, trabajadores *bajos, std::vector<weib_d> &dist){
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
      if(ind < Pa.size()){mother_reaction(Pa, La, ind, altos, 2, 5);    continue_reaction(La, LTAa, altos, ran);}
    }
    if(PTa.size() != 0){
      it = std::find(PTa.begin(), PTa.end(), agent);      ind = std::distance(PTa.begin(), it);
      if(ind < PTa.size()){mother_reaction(PTa, LTa, ind, altos, 3, 6);}
    }
    if(PTAa.size() != 0){
      it = std::find(PTAa.begin(), PTAa.end(), agent);      ind = std::distance(PTAa.begin(), it);
      if(ind < PTAa.size()){mother_reaction(PTAa, LTAa, ind, altos, 4, 7);}
    }
  }
}


void reaction5(grupo &Sa, grupo &Sb, grupo &STa, grupo &STb, grupo &SMa, grupo &SMb, grupo &Ea, grupo &Eb, grupo &ETa, grupo &ETb, grupo &EMa, grupo &EMb, grupo &Pa, grupo &Pb, grupo &PTa, grupo &PTb, grupo &PTAa, grupo &PTAb, grupo &La, grupo &Lb, grupo &LTa, grupo &LTb, grupo &LTAa, grupo &LTAb, grupo &IAa, grupo &IAb, grupo &RTa, grupo&RTb, grupo &RIa, grupo &RIb, grupo &RMa, grupo &RMb, grupo &RAa, grupo &RAb, Crandom &ran, trabajadores *altos, trabajadores *bajos, std::vector<weib_d> &dist){
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
      if(ind < Pb.size()){mother_reaction(Pb, Lb, ind, bajos, 2, 5);    continue_reaction(Lb, LTAb, bajos, ran);}
    }
    if(PTb.size() != 0){
      it = std::find(PTb.begin(), PTb.end(), agent);      ind = std::distance(PTb.begin(), it);
      if(ind < PTb.size()){mother_reaction(PTb, LTb, ind, bajos, 3, 6);}
    }
    if(PTAb.size() != 0){
      it = std::find(PTAb.begin(), PTAb.end(), agent);      ind = std::distance(PTAb.begin(), it);
      if(ind < PTAb.size()){mother_reaction(PTAb, LTAb, ind, bajos, 4, 7);}
    }
  }
}


void reaction6(grupo &Sa, grupo &Sb, grupo &STa, grupo &STb, grupo &SMa, grupo &SMb, grupo &Ea, grupo &Eb, grupo &ETa, grupo &ETb, grupo &EMa, grupo &EMb, grupo &Pa, grupo &Pb, grupo &PTa, grupo &PTb, grupo &PTAa, grupo &PTAb, grupo &La, grupo &Lb, grupo &LTa, grupo &LTb, grupo &LTAa, grupo &LTAb, grupo &IAa, grupo &IAb, grupo &RTa, grupo&RTb, grupo &RIa, grupo &RIb, grupo &RMa, grupo &RMb, grupo &RAa, grupo &RAb, Crandom &ran, trabajadores *altos, trabajadores *bajos, std::vector<weib_d> &dist){
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
      if(ind < Pa.size()){mother_reaction(Pa, IAa, ind, altos, 2, 9);}
    }
    if(PTa.size() != 0){
      it = std::find(PTa.begin(), PTa.end(), agent);      ind = std::distance(PTa.begin(), it);
      if(ind < PTa.size()){mother_reaction(PTa, IAa, ind, altos, 3, 9);}
    }
    if(PTAa.size() != 0){
      it = std::find(PTAa.begin(), PTAa.end(), agent);      ind = std::distance(PTAa.begin(), it);
      if(ind < PTAa.size()){mother_reaction(PTAa, IAa, ind, altos, 4, 9);}
    }
  }
}


void reaction7(grupo &Sa, grupo &Sb, grupo &STa, grupo &STb, grupo &SMa, grupo &SMb, grupo &Ea, grupo &Eb, grupo &ETa, grupo &ETb, grupo &EMa, grupo &EMb, grupo &Pa, grupo &Pb, grupo &PTa, grupo &PTb, grupo &PTAa, grupo &PTAb, grupo &La, grupo &Lb, grupo &LTa, grupo &LTb, grupo &LTAa, grupo &LTAb, grupo &IAa, grupo &IAb, grupo &RTa, grupo&RTb, grupo &RIa, grupo &RIb, grupo &RMa, grupo &RMb, grupo &RAa, grupo &RAb, Crandom &ran, trabajadores *altos, trabajadores *bajos, std::vector<weib_d> &dist){
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
      if(ind < Pb.size()){mother_reaction(Pb, IAb, ind, bajos, 2, 9);}
    }
    if(PTb.size() != 0){
      it = std::find(PTb.begin(), PTb.end(), agent);      ind = std::distance(PTb.begin(), it);
      if(ind < PTb.size()){mother_reaction(PTb, IAb, ind, bajos, 3, 9);}
    }
    if(PTAb.size() != 0){
      it = std::find(PTAb.begin(), PTAb.end(), agent);      ind = std::distance(PTAb.begin(), it);
      if(ind < PTAb.size()){mother_reaction(PTAb, IAb, ind, bajos, 4, 9);}
    }
  }
}


void reaction8(grupo &Sa, grupo &Sb, grupo &STa, grupo &STb, grupo &SMa, grupo &SMb, grupo &Ea, grupo &Eb, grupo &ETa, grupo &ETb, grupo &EMa, grupo &EMb, grupo &Pa, grupo &Pb, grupo &PTa, grupo &PTb, grupo &PTAa, grupo &PTAb, grupo &La, grupo &Lb, grupo &LTa, grupo &LTb, grupo &LTAa, grupo &LTAb, grupo &IAa, grupo &IAb, grupo &RTa, grupo&RTb, grupo &RIa, grupo &RIb, grupo &RMa, grupo &RMb, grupo &RAa, grupo &RAb, Crandom &ran, trabajadores *altos, trabajadores *bajos, std::vector<weib_d> &dist){
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
      if(ind < Pa.size()){mother_reaction(Pa, RIa, ind, altos, 2, 10);}
    }
    if(PTa.size() != 0){
      it = std::find(PTa.begin(), PTa.end(), agent);      ind = std::distance(PTa.begin(), it);
      if(ind < PTa.size()){mother_reaction(PTa, RAa, ind, altos, 3, 10);}
    }
    if(PTAa.size() != 0){
      it = std::find(PTAa.begin(), PTAa.end(), agent);      ind = std::distance(PTAa.begin(), it);
      if(ind < PTAa.size()){mother_reaction(PTAa, RAa, ind, altos, 4, 10);}
    }
  }
}


void reaction9(grupo &Sa, grupo &Sb, grupo &STa, grupo &STb, grupo &SMa, grupo &SMb, grupo &Ea, grupo &Eb, grupo &ETa, grupo &ETb, grupo &EMa, grupo &EMb, grupo &Pa, grupo &Pb, grupo &PTa, grupo &PTb, grupo &PTAa, grupo &PTAb, grupo &La, grupo &Lb, grupo &LTa, grupo &LTb, grupo &LTAa, grupo &LTAb, grupo &IAa, grupo &IAb, grupo &RTa, grupo&RTb, grupo &RIa, grupo &RIb, grupo &RMa, grupo &RMb, grupo &RAa, grupo &RAb, Crandom &ran, trabajadores *altos, trabajadores *bajos, std::vector<weib_d> &dist){
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
      if(ind < Pb.size()){mother_reaction(Pb, RIb, ind, bajos, 2, 10);}
    }
    if(PTb.size() != 0){
      it = std::find(PTb.begin(), PTb.end(), agent);      ind = std::distance(PTb.begin(), it);
      if(ind < PTb.size()){mother_reaction(PTb, RAb, ind, bajos, 3, 10);}
    }
    if(PTAb.size() != 0){
      it = std::find(PTAb.begin(), PTAb.end(), agent);      ind = std::distance(PTAb.begin(), it);
      if(ind < PTAb.size()){mother_reaction(PTAb, RAb, ind, bajos, 4, 10);}
    }
  }
}


void reaction10(grupo &Sa, grupo &Sb, grupo &STa, grupo &STb, grupo &SMa, grupo &SMb, grupo &Ea, grupo &Eb, grupo &ETa, grupo &ETb, grupo &EMa, grupo &EMb, grupo &Pa, grupo &Pb, grupo &PTa, grupo &PTb, grupo &PTAa, grupo &PTAb, grupo &La, grupo &Lb, grupo &LTa, grupo &LTb, grupo &LTAa, grupo &LTAb, grupo &IAa, grupo &IAb, grupo &RTa, grupo&RTb, grupo &RIa, grupo &RIb, grupo &RMa, grupo &RMb, grupo &RAa, grupo &RAb, Crandom &ran, trabajadores *altos, trabajadores *bajos, std::vector<weib_d> &dist){
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
      if(ind < La.size()){mother_reaction(La, RIa, ind, altos, 5, 10);}
    }
    if(LTa.size() != 0){
      it = std::find(LTa.begin(), LTa.end(), agent);      ind = std::distance(LTa.begin(), it);
      if(ind < LTa.size()){mother_reaction(LTa, RAa, ind, altos, 6, 10);}
    }
    if(LTAa.size() != 0){
      it = std::find(LTAa.begin(), LTAa.end(), agent);      ind = std::distance(LTAa.begin(), it);
      if(ind < LTAa.size()){mother_reaction(LTAa, RAa, ind, altos, 7, 10);}
    }
  }
}


void reaction11(grupo &Sa, grupo &Sb, grupo &STa, grupo &STb, grupo &SMa, grupo &SMb, grupo &Ea, grupo &Eb, grupo &ETa, grupo &ETb, grupo &EMa, grupo &EMb, grupo &Pa, grupo &Pb, grupo &PTa, grupo &PTb, grupo &PTAa, grupo &PTAb, grupo &La, grupo &Lb, grupo &LTa, grupo &LTb, grupo &LTAa, grupo &LTAb, grupo &IAa, grupo &IAb, grupo &RTa, grupo&RTb, grupo &RIa, grupo &RIb, grupo &RMa, grupo &RMb, grupo &RAa, grupo &RAb, Crandom &ran, trabajadores *altos, trabajadores *bajos, std::vector<weib_d> &dist){
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
      if(ind < Lb.size()){mother_reaction(Lb, RIb, ind, bajos, 5, 10);}
    }
    if(LTb.size() != 0){
      it = std::find(LTb.begin(), LTb.end(), agent);      ind = std::distance(LTb.begin(), it);
      if(ind < LTb.size()){mother_reaction(LTb, RAb, ind, bajos, 6, 10);}
    }
    if(LTAb.size() != 0){
      it = std::find(LTAb.begin(), LTAb.end(), agent);      ind = std::distance(LTAb.begin(), it);
      if(ind < LTAb.size()){mother_reaction(LTAb, RAb, ind, bajos, 7, 10);}
    }
  }
}


void reaction12(grupo &Sa, grupo &Sb, grupo &STa, grupo &STb, grupo &SMa, grupo &SMb, grupo &Ea, grupo &Eb, grupo &ETa, grupo &ETb, grupo &EMa, grupo &EMb, grupo &Pa, grupo &Pb, grupo &PTa, grupo &PTb, grupo &PTAa, grupo &PTAb, grupo &La, grupo &Lb, grupo &LTa, grupo &LTb, grupo &LTAa, grupo &LTAb, grupo &IAa, grupo &IAb, grupo &RTa, grupo&RTb, grupo &RIa, grupo &RIb, grupo &RMa, grupo &RMb, grupo &RAa, grupo &RAb, Crandom &ran, trabajadores *altos, trabajadores *bajos, std::vector<weib_d> &dist){
  std::vector<int> aux;
  
  for(unsigned int i=0; i<IAa.size(); i++){if(altos[IAa[i]].tstate > TM){aux.push_back(IAa[i]);}}

  if(aux.size() != 0){    
    unsigned int index = index_time(aux, altos, dist[10], (double)ran.r());
    unsigned int agent = aux[index];
    std::vector<int>::iterator it = std::find(IAa.begin(), IAa.end(), agent);
    unsigned int ind = std::distance(IAa.begin(), it);
    if(ind < IAa.size()){mother_reaction(IAa, RAa, ind, altos, 9, 10);}
  }
}


void reaction13(grupo &Sa, grupo &Sb, grupo &STa, grupo &STb, grupo &SMa, grupo &SMb, grupo &Ea, grupo &Eb, grupo &ETa, grupo &ETb, grupo &EMa, grupo &EMb, grupo &Pa, grupo &Pb, grupo &PTa, grupo &PTb, grupo &PTAa, grupo &PTAb, grupo &La, grupo &Lb, grupo &LTa, grupo &LTb, grupo &LTAa, grupo &LTAb, grupo &IAa, grupo &IAb, grupo &RTa, grupo&RTb, grupo &RIa, grupo &RIb, grupo &RMa, grupo &RMb, grupo &RAa, grupo &RAb, Crandom &ran, trabajadores *altos, trabajadores *bajos, std::vector<weib_d> &dist){
  std::vector<int> aux;
  
  for(unsigned int i=0; i<IAb.size(); i++){if(bajos[IAb[i]].tstate > TM){aux.push_back(IAb[i]);}}

  if(aux.size() != 0){
    unsigned int index = index_time(aux, bajos, dist[11], (double)ran.r());
    unsigned int agent = aux[index];
    std::vector<int>::iterator it = std::find(IAb.begin(), IAb.end(), agent);
    unsigned int ind = std::distance(IAb.begin(), it);
    if(ind < IAb.size()){mother_reaction(IAb, RAb, ind, bajos, 9, 10);}
  }  
}
