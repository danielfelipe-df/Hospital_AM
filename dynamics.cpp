#include <iostream>
#include <math.h>
#include <algorithm>
#include "dynamics.h"

std::vector<double> contagio(double Sa, double Sb, double STa, double STb, double Ea, double Eb, double ETa, double ETb, double Pa, double Pb, double PTa, double PTb, double PTAa, double PTAb, double La, double Lb, double LTa, double LTb, double LTAa, double LTAb, double IAa, double IAb, double Na, double Nb, double prev, Crandom &ran, double t, double* tj){

  //Número de proponsidades
  const unsigned int n = 14;

  //Creo el arreglo de las propensidades
  double As[n];

  //Propensidades de exponerse
  As[0] = beta*(Sa+STa)*((Pa+PTa+La+LTa)/Na + mu*(Pb+PTb+Lb+LTb)/Nb + (1-alpha)*(IAa+PTAa+LTAa)/Na + (1-alpha)*mu*(IAb+PTAb+LTAb)/Nb + eta*prev);
  As[1] = beta*(Sb+STb)*(mu*(Pa+PTa+La+LTa)/Na + chi*(Pb+PTb+Lb+LTb)/Nb + (1-alpha)*mu*(IAa+PTAa+LTAa)/Na + (1-alpha)*chi*(IAb+PTAb+LTAb)/Nb);

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
  std::vector<weib_d> dist;
  //for(unsigned int i=2; i<n; i++){weib_d my_dist(1.5, 1.0/As[i]);    dist.push_back(my_dist);}

  //Hallo el tiempo en el que va a pasar la siguiente reacción con el método NMGA
  double prom = 230.0, sigma = 55.0, B = beta*Sa*eta*prev;
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


void massive_reaction(grupo &S, grupo &ST, grupo &E, grupo &ET, grupo &P, grupo &PT, grupo &L, grupo &LT, grupo &RI, grupo &RT, Crandom &ran, trabajadores *family){
  unsigned int M = S.size() + E.size() + P.size() + L.size() + RI.size();
  unsigned int num = (int)(ran.r()*M);  

  unsigned int agent, index;  
  if(num < P.size()){
    if(ran.r() < xi){
      //mother_reaction(P, PT, (int)(ran.r()*P.size()), family, 2, 3);
      index = (int)(ran.r()*P.size());      agent = P[index];
      P.erase(P.begin() + index);      PT.push_back(agent);
      family[agent].change(3, 2);
      family[agent].time = 0.0;
      family[agent].tmax = Tt;
    }
  }
  else if(num < P.size() + L.size()){
    if(ran.r() < xi){
      //mother_reaction(L, LT, (int)(ran.r()*L.size()), family, 5, 6);
      index = (int)(ran.r()*L.size());      agent = L[index];
      L.erase(L.begin() + index);      LT.push_back(agent);
      family[agent].change(6, 5);
      family[agent].time = 0.0;
      family[agent].tmax = Tt;
    }
  }
  else if(num < P.size() + L.size() + S.size()){
    index = (int)(ran.r()*S.size());      agent = S[index];
    S.erase(S.begin() + index);      ST.push_back(agent);
    family[agent].change(6, 5);
    family[agent].time = 0.0;
    family[agent].tmax = Tt;
  }
  else if(num < P.size() + L.size() + S.size() + E.size()){
    index = (int)(ran.r()*E.size());      agent = E[index];
    E.erase(E.begin() + index);      ET.push_back(agent);
    family[agent].change(6, 5);
    family[agent].time = 0.0;
    family[agent].tmax = Tt;
  }
  else{
    index = (int)(ran.r()*RI.size());      agent = RI[index];
    RI.erase(RI.begin() + index);      RT.push_back(agent);
    family[agent].change(6, 5);
    family[agent].time = 0.0;
    family[agent].tmax = Tt;
  }
}


void continue_reaction(grupo &L, grupo &LT, trabajadores *family, Crandom &ran){
  if(ran.r() < iota*xi){
    int agent = L.back();
    L.erase( L.end() - 1);
    LT.push_back(agent);
    family[agent].change(6,5);
    family[agent].time = 0.0;
    family[agent].tmax = Tt;
  }
}


void tested_isolated(grupo &T, grupo &TA, trabajadores *family, double time, int typeout, int typein){
  int index;
  for(unsigned int i=0; i<T.size(); i++){
    index = T[i];
    family[index].time += time;
    if(family[index].time > family[index].tmax){
      family[index].change(typein, typeout);
      family[index].time = 0.0;
      family[index].tmax = 0.0;
      T.erase( T.begin() + i);
      TA.push_back(index);
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


void reaction0(grupo &Sa, grupo &Sb, grupo &STa, grupo &STb, grupo &Ea, grupo &Eb, grupo &ETa, grupo &ETb, grupo &Pa, grupo &Pb, grupo &PTa, grupo &PTb, grupo &PTAa, grupo &PTAb, grupo &La, grupo &Lb, grupo &LTa, grupo &LTb, grupo &LTAa, grupo &LTAb, grupo &IAa, grupo &IAb, grupo &RTa, grupo&RTb, grupo &RIa, grupo &RIb, grupo &RAa, grupo &RAb, Crandom &ran, trabajadores *altos, trabajadores *bajos, std::vector<weib_d> &dist){
  unsigned int index = (int)(ran.r()*(Sa.size() + STa.size()));
  if(index < Sa.size()){mother_reaction(Sa, Ea, index, altos, 0, 1);}
  else{mother_reaction(STa, ETa, index-Sa.size(), altos, 0, 1);}
}


void reaction1(grupo &Sa, grupo &Sb, grupo &STa, grupo &STb, grupo &Ea, grupo &Eb, grupo &ETa, grupo &ETb, grupo &Pa, grupo &Pb, grupo &PTa, grupo &PTb, grupo &PTAa, grupo &PTAb, grupo &La, grupo &Lb, grupo &LTa, grupo &LTb, grupo &LTAa, grupo &LTAb, grupo &IAa, grupo &IAb, grupo &RTa, grupo&RTb, grupo &RIa, grupo &RIb, grupo &RAa, grupo &RAb, Crandom &ran, trabajadores *altos, trabajadores *bajos, std::vector<weib_d> &dist){
  unsigned int index = (int)(ran.r()*(Sb.size() + STb.size()));
  if(index < Sb.size()){mother_reaction(Sb, Eb, index, bajos, 0, 1);}
  else{mother_reaction(STb, ETb, index-Sb.size(), bajos, 0, 1);}
}


void reaction2(grupo &Sa, grupo &Sb, grupo &STa, grupo &STb, grupo &Ea, grupo &Eb, grupo &ETa, grupo &ETb, grupo &Pa, grupo &Pb, grupo &PTa, grupo &PTb, grupo &PTAa, grupo &PTAb, grupo &La, grupo &Lb, grupo &LTa, grupo &LTb, grupo &LTAa, grupo &LTAb, grupo &IAa, grupo &IAb, grupo &RTa, grupo&RTb, grupo &RIa, grupo &RIb, grupo &RAa, grupo &RAb, Crandom &ran, trabajadores *altos, trabajadores *bajos, std::vector<weib_d> &dist){
  int index[2] = {-10, -10};
  double time[2] = {1e6, 1e6};
  double value = (double)ran.r();
  double min = 1e6;
  if(Ea.size() != 0){
    index[0] = index_time(Ea, altos, dist[0], value);
    time[0] = altos[Ea[index[0]]].tstate;
    min = ((min > time[0]) ? time[0] : min);
  }
  if(ETa.size() != 0){
    index[1] = index_time(ETa, altos, dist[0], value);
    time[1] = altos[ETa[index[1]]].tstate;
    min = ((min > time[1]) ? time[1] : min);
  }  

  if(min == time[0]){mother_reaction(Ea, Pa, index[0], altos, 1, 2);}
  else{mother_reaction(ETa, Pa, index[1], altos, 1, 2);}
}


void reaction3(grupo &Sa, grupo &Sb, grupo &STa, grupo &STb, grupo &Ea, grupo &Eb, grupo &ETa, grupo &ETb, grupo &Pa, grupo &Pb, grupo &PTa, grupo &PTb, grupo &PTAa, grupo &PTAb, grupo &La, grupo &Lb, grupo &LTa, grupo &LTb, grupo &LTAa, grupo &LTAb, grupo &IAa, grupo &IAb, grupo &RTa, grupo&RTb, grupo &RIa, grupo &RIb, grupo &RAa, grupo &RAb, Crandom &ran, trabajadores *altos, trabajadores *bajos, std::vector<weib_d> &dist){
  int index[2] = {-10, -10};
  double time[2] = {1e6, 1e6};
  double value = (double)ran.r();
  double min = 1e6;
  if(Eb.size() != 0){
    index[0] = index_time(Eb, bajos, dist[1], value);
    time[0] = bajos[Eb[index[0]]].tstate;
    min = ((min > time[0]) ? time[0] : min);
  }
  if(ETb.size() != 0){
    index[1] = index_time(ETb, bajos, dist[1], value);
    time[1] = bajos[ETb[index[1]]].tstate;
    min = ((min > time[1]) ? time[1] : min);
  }  

  if(min == time[0]){mother_reaction(Eb, Pb, index[0], bajos, 1, 2);}
  else{mother_reaction(ETb, Pb, index[1], bajos, 1, 2);}
}


void reaction4(grupo &Sa, grupo &Sb, grupo &STa, grupo &STb, grupo &Ea, grupo &Eb, grupo &ETa, grupo &ETb, grupo &Pa, grupo &Pb, grupo &PTa, grupo &PTb, grupo &PTAa, grupo &PTAb, grupo &La, grupo &Lb, grupo &LTa, grupo &LTb, grupo &LTAa, grupo &LTAb, grupo &IAa, grupo &IAb, grupo &RTa, grupo&RTb, grupo &RIa, grupo &RIb, grupo &RAa, grupo &RAb, Crandom &ran, trabajadores *altos, trabajadores *bajos, std::vector<weib_d> &dist){
  int index[3] = {-10, -10, -10};
  double time[3] = {1e6, 1e6, 1e6};
  double value = (double)ran.r();
  double min = 1e6;
  if(Pa.size() != 0){
    index[0] = index_time(Pa, altos, dist[2], value);
    time[0] = altos[Pa[index[0]]].tstate;
    min = ((min > time[0]) ? time[0] : min);
  }
  if(PTa.size() != 0){
    index[1] = index_time(PTa, altos, dist[2], value);
    time[1] = altos[PTa[index[1]]].tstate;
    min = ((min > time[1]) ? time[1] : min);
  }
  if(PTAa.size() != 0){
    index[2] = index_time(PTAa, altos, dist[2], value);
    time[2] = altos[PTAa[index[2]]].tstate;
    min = ((min > time[2]) ? time[2] : min);
  }

  if(min == time[0]){mother_reaction(Pa, La, index[0], altos, 2, 5);    continue_reaction(La, LTa, altos, ran);}
  else if(min == time[1]){mother_reaction(PTa, LTa, index[1], altos, 3, 6);}
  else{mother_reaction(PTAa, LTAa, index[2], altos, 4, 7);}
}


void reaction5(grupo &Sa, grupo &Sb, grupo &STa, grupo &STb, grupo &Ea, grupo &Eb, grupo &ETa, grupo &ETb, grupo &Pa, grupo &Pb, grupo &PTa, grupo &PTb, grupo &PTAa, grupo &PTAb, grupo &La, grupo &Lb, grupo &LTa, grupo &LTb, grupo &LTAa, grupo &LTAb, grupo &IAa, grupo &IAb, grupo &RTa, grupo&RTb, grupo &RIa, grupo &RIb, grupo &RAa, grupo &RAb, Crandom &ran, trabajadores *altos, trabajadores *bajos, std::vector<weib_d> &dist){
  int index[3] = {-10, -10, -10};
  double time[3] = {1e6, 1e6, 1e6};
  double value = (double)ran.r();
  double min = 1e6;
  if(Pb.size() != 0){
    index[0] = index_time(Pb, bajos, dist[3], value);
    time[0] = bajos[Pb[index[0]]].tstate;
    min = ((min > time[0]) ? time[0] : min);
  }
  if(PTb.size() != 0){
    index[1] = index_time(PTb, bajos, dist[3], value);
    time[1] = bajos[PTb[index[1]]].tstate;
    min = ((min > time[1]) ? time[1] : min);
  }
  if(PTAb.size() != 0){
    index[2] = index_time(PTAb, bajos, dist[3], value);
    time[2] = bajos[PTAb[index[2]]].tstate;
    min = ((min > time[2]) ? time[2] : min);
  }

  if(min == time[0]){mother_reaction(Pb, Lb, index[0], bajos, 2, 5);    continue_reaction(Lb, LTb, bajos, ran);}
  else if(min == time[1]){mother_reaction(PTb, LTb, index[1], bajos, 3, 6);}
  else{mother_reaction(PTAb, LTAb, index[2], bajos, 4, 7);}
}


void reaction6(grupo &Sa, grupo &Sb, grupo &STa, grupo &STb, grupo &Ea, grupo &Eb, grupo &ETa, grupo &ETb, grupo &Pa, grupo &Pb, grupo &PTa, grupo &PTb, grupo &PTAa, grupo &PTAb, grupo &La, grupo &Lb, grupo &LTa, grupo &LTb, grupo &LTAa, grupo &LTAb, grupo &IAa, grupo &IAb, grupo &RTa, grupo&RTb, grupo &RIa, grupo &RIb, grupo &RAa, grupo &RAb, Crandom &ran, trabajadores *altos, trabajadores *bajos, std::vector<weib_d> &dist){  
  int index[3] = {-10, -10, -10};
  double time[3] = {1e6, 1e6, 1e6};
  double value = (double)ran.r();
  double min = 1e6;
  if(Pa.size() != 0){
    index[0] = index_time(Pa, altos, dist[4], value);
    time[0] = altos[Pa[index[0]]].tstate;
    min = ((min > time[0]) ? time[0] : min);
  }
  if(PTa.size() != 0){
    index[1] = index_time(PTa, altos, dist[4], value);
    time[1] = altos[PTa[index[1]]].tstate;
    min = ((min > time[1]) ? time[1] : min);
  }
  if(PTAa.size() != 0){
    index[2] = index_time(PTAa, altos, dist[4], value);
    time[2] = altos[PTAa[index[2]]].tstate;
    min = ((min > time[2]) ? time[2] : min);
  }

  if(min == time[0]){mother_reaction(Pa, IAa, index[0], altos, 2, 9);}
  else if(min == time[1]){mother_reaction(PTa, IAa, index[1], altos, 3, 9);}
  else{mother_reaction(PTAa, IAa, index[2], altos, 4, 9);}
}


void reaction7(grupo &Sa, grupo &Sb, grupo &STa, grupo &STb, grupo &Ea, grupo &Eb, grupo &ETa, grupo &ETb, grupo &Pa, grupo &Pb, grupo &PTa, grupo &PTb, grupo &PTAa, grupo &PTAb, grupo &La, grupo &Lb, grupo &LTa, grupo &LTb, grupo &LTAa, grupo &LTAb, grupo &IAa, grupo &IAb, grupo &RTa, grupo&RTb, grupo &RIa, grupo &RIb, grupo &RAa, grupo &RAb, Crandom &ran, trabajadores *altos, trabajadores *bajos, std::vector<weib_d> &dist){
  int index[3] = {-10, -10, -10};
  double time[3] = {1e6, 1e6, 1e6};
  double value = (double)ran.r();
  double min = 1e6;
  if(Pb.size() != 0){
    index[0] = index_time(Pb, bajos, dist[5], value);
    time[0] = bajos[Pb[index[0]]].tstate;
    min = ((min > time[0]) ? time[0] : min);
  }
  if(PTb.size() != 0){
    index[1] = index_time(PTb, bajos, dist[5], value);
    time[1] = bajos[PTb[index[1]]].tstate;
    min = ((min > time[1]) ? time[1] : min);
  }
  if(PTAb.size() != 0){
    index[2] = index_time(PTAb, bajos, dist[5], value);
    time[2] = bajos[PTAb[index[2]]].tstate;
    min = ((min > time[2]) ? time[2] : min);
  }

  if(min == time[0]){mother_reaction(Pb, IAb, index[0], bajos, 2, 9);}
  else if(min == time[1]){mother_reaction(PTb, IAb, index[1], bajos, 3, 9);}
  else{mother_reaction(PTAb, IAb, index[2], bajos, 4, 9);}
}


void reaction8(grupo &Sa, grupo &Sb, grupo &STa, grupo &STb, grupo &Ea, grupo &Eb, grupo &ETa, grupo &ETb, grupo &Pa, grupo &Pb, grupo &PTa, grupo &PTb, grupo &PTAa, grupo &PTAb, grupo &La, grupo &Lb, grupo &LTa, grupo &LTb, grupo &LTAa, grupo &LTAb, grupo &IAa, grupo &IAb, grupo &RTa, grupo&RTb, grupo &RIa, grupo &RIb, grupo &RAa, grupo &RAb, Crandom &ran, trabajadores *altos, trabajadores *bajos, std::vector<weib_d> &dist){
  int index[3] = {-10, -10, -10};
  double time[3] = {1e6, 1e6, 1e6};
  double value = (double)ran.r();
  double min = 1e6;
  if(Pa.size() != 0){
    index[0] = index_time(Pa, altos, dist[6], value);
    time[0] = altos[Pa[index[0]]].tstate;
    min = ((min > time[0]) ? time[0] : min);
  }
  if(PTa.size() != 0){
    index[1] = index_time(PTa, altos, dist[6], value);
    time[1] = altos[PTa[index[1]]].tstate;
    min = ((min > time[1]) ? time[1] : min);
  }
  if(PTAa.size() != 0){
    index[2] = index_time(PTAa, altos, dist[6], value);
    time[2] = altos[PTAa[index[2]]].tstate;
    min = ((min > time[2]) ? time[2] : min);
  }

  if(min == time[0]){mother_reaction(Pa, RIa, index[0], altos, 2, 10);}
  else if(min == time[1]){mother_reaction(PTa, RAa, index[1], altos, 3, 10);}
  else{mother_reaction(PTAa, RAa, index[2], altos, 4, 10);}
}


void reaction9(grupo &Sa, grupo &Sb, grupo &STa, grupo &STb, grupo &Ea, grupo &Eb, grupo &ETa, grupo &ETb, grupo &Pa, grupo &Pb, grupo &PTa, grupo &PTb, grupo &PTAa, grupo &PTAb, grupo &La, grupo &Lb, grupo &LTa, grupo &LTb, grupo &LTAa, grupo &LTAb, grupo &IAa, grupo &IAb, grupo &RTa, grupo&RTb, grupo &RIa, grupo &RIb, grupo &RAa, grupo &RAb, Crandom &ran, trabajadores *altos, trabajadores *bajos, std::vector<weib_d> &dist){
  int index[3] = {-10, -10, -10};
  double time[3] = {1e6, 1e6, 1e6};
  double value = (double)ran.r();
  double min = 1e6;
  if(Pb.size() != 0){
    index[0] = index_time(Pb, bajos, dist[7], value);
    time[0] = bajos[Pb[index[0]]].tstate;
    min = ((min > time[0]) ? time[0] : min);
  }
  if(PTb.size() != 0){
    index[1] = index_time(PTb, bajos, dist[7], value);
    time[1] = bajos[PTb[index[1]]].tstate;
    min = ((min > time[1]) ? time[1] : min);
  }
  if(PTAb.size() != 0){
    index[2] = index_time(PTAb, bajos, dist[7], value);
    time[2] = bajos[PTAb[index[2]]].tstate;
    min = ((min > time[2]) ? time[2] : min);
  }

  if(min == time[0]){mother_reaction(Pb, RIb, index[0], bajos, 2, 10);}
  else if(min == time[1]){mother_reaction(PTb, RAb, index[1], bajos, 3, 10);}
  else{mother_reaction(PTAb, RAb, index[2], bajos, 4, 10);}
}


void reaction10(grupo &Sa, grupo &Sb, grupo &STa, grupo &STb, grupo &Ea, grupo &Eb, grupo &ETa, grupo &ETb, grupo &Pa, grupo &Pb, grupo &PTa, grupo &PTb, grupo &PTAa, grupo &PTAb, grupo &La, grupo &Lb, grupo &LTa, grupo &LTb, grupo &LTAa, grupo &LTAb, grupo &IAa, grupo &IAb, grupo &RTa, grupo&RTb, grupo &RIa, grupo &RIb, grupo &RAa, grupo &RAb, Crandom &ran, trabajadores *altos, trabajadores *bajos, std::vector<weib_d> &dist){
  int index[3] = {-10, -10, -10};
  double time[3] = {1e6, 1e6, 1e6};
  double value = (double)ran.r();
  double min = 1e6;
  if(La.size() != 0){
    index[0] = index_time(La, altos, dist[8], value);
    time[0] = altos[La[index[0]]].tstate;
    min = ((min > time[0]) ? time[0] : min);
  }
  if(LTa.size() != 0){
    index[1] = index_time(LTa, altos, dist[8], value);
    time[1] = altos[LTa[index[1]]].tstate;
    min = ((min > time[1]) ? time[1] : min);
  }
  if(LTAa.size() != 0){
    index[2] = index_time(LTAa, altos, dist[8], value);
    time[2] = altos[LTAa[index[2]]].tstate;
    min = ((min > time[2]) ? time[2] : min);
  }

  if(min == time[0]){mother_reaction(La, RIa, index[0], altos, 5, 10);}
  else if(min == time[1]){mother_reaction(LTa, RAa, index[1], altos, 6, 10);}
  else{mother_reaction(LTAa, RAa, index[2], altos, 7, 10);}
}


void reaction11(grupo &Sa, grupo &Sb, grupo &STa, grupo &STb, grupo &Ea, grupo &Eb, grupo &ETa, grupo &ETb, grupo &Pa, grupo &Pb, grupo &PTa, grupo &PTb, grupo &PTAa, grupo &PTAb, grupo &La, grupo &Lb, grupo &LTa, grupo &LTb, grupo &LTAa, grupo &LTAb, grupo &IAa, grupo &IAb, grupo &RTa, grupo&RTb, grupo &RIa, grupo &RIb, grupo &RAa, grupo &RAb, Crandom &ran, trabajadores *altos, trabajadores *bajos, std::vector<weib_d> &dist){
  int index[3] = {-10, -10, -10};
  double time[3] = {1e6, 1e6, 1e6};
  double value = (double)ran.r();
  double min = 1e6;
  if(Lb.size() != 0){
    index[0] = index_time(Lb, bajos, dist[9], value);
    time[0] = bajos[Lb[index[0]]].tstate;
    min = ((min > time[0]) ? time[0] : min);
  }
  if(LTb.size() != 0){
    index[1] = index_time(LTb, bajos, dist[9], value);
    time[1] = bajos[LTb[index[1]]].tstate;
    min = ((min > time[1]) ? time[1] : min);
  }
  if(LTAb.size() != 0){
    index[2] = index_time(LTAb, bajos, dist[9], value);
    time[2] = bajos[LTAb[index[2]]].tstate;
    min = ((min > time[2]) ? time[2] : min);
  }

  if(min == time[0]){mother_reaction(Lb, RIb, index[0], bajos, 5, 10);}
  else if(min == time[1]){mother_reaction(LTb, RAb, index[1], bajos, 6, 10);}
  else{mother_reaction(LTAb, RAb, index[2], bajos, 7, 10);}
}


void reaction12(grupo &Sa, grupo &Sb, grupo &STa, grupo &STb, grupo &Ea, grupo &Eb, grupo &ETa, grupo &ETb, grupo &Pa, grupo &Pb, grupo &PTa, grupo &PTb, grupo &PTAa, grupo &PTAb, grupo &La, grupo &Lb, grupo &LTa, grupo &LTb, grupo &LTAa, grupo &LTAb, grupo &IAa, grupo &IAb, grupo &RTa, grupo&RTb, grupo &RIa, grupo &RIb, grupo &RAa, grupo &RAb, Crandom &ran, trabajadores *altos, trabajadores *bajos, std::vector<weib_d> &dist){
  int index = index_time(IAa, altos, dist[10], (double)ran.r());
  mother_reaction(IAa, RAa, index, altos, 9, 10);
}


void reaction13(grupo &Sa, grupo &Sb, grupo &STa, grupo &STb, grupo &Ea, grupo &Eb, grupo &ETa, grupo &ETb, grupo &Pa, grupo &Pb, grupo &PTa, grupo &PTb, grupo &PTAa, grupo &PTAb, grupo &La, grupo &Lb, grupo &LTa, grupo &LTb, grupo &LTAa, grupo &LTAb, grupo &IAa, grupo &IAb, grupo &RTa, grupo&RTb, grupo &RIa, grupo &RIb, grupo &RAa, grupo &RAb, Crandom &ran, trabajadores *altos, trabajadores *bajos, std::vector<weib_d> &dist){
  int index = index_time(IAb, bajos, dist[11], (double)ran.r());
  mother_reaction(IAb, RAb, index, bajos, 9, 10);
}
