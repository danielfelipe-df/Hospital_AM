#include <math.h>
#include "dynamics.h"


std::vector<double> contagio(std::vector<grupo> &Val, std::vector<grupo> &Vba, double prev, Crandom &ran, double t, double* tj){
  //Calculo los tamaños de cada vector
  double Sa, Sb, STa, STb, Ea, Eb, ETa, ETb, EAa, EAb, Pa, Pb, PTa, PTb, PTAa, PTAb, La, Lb, LTa, LTb, LTAa, LTAb, IAa, IAb;
  Sa = Val[0].size();  Sb = Vba[0].size();  STa = Val[1].size();  STb = Vba[1].size();
  Ea = Val[2].size();  Eb = Vba[2].size();  ETa = Val[3].size();  ETb = Vba[3].size();  EAa = Val[4].size();  EAb = Val[4].size();
  Pa = Val[5].size();  Pb = Vba[5].size();  PTa = Val[6].size();  PTb = Vba[6].size();  PTAa = Val[7].size();  PTAb = Vba[7].size();
  La = Val[8].size();  Lb = Vba[8].size();  LTa = Val[9].size();  LTb = Vba[9].size();  LTAa = Val[10].size();  LTAb = Vba[10].size();
  IAa = Val[11].size();  IAb = Vba[11].size();

  //Número de propensidades
  const unsigned int n = 14;

  //Creo el arreglo de las propensidades
  double As[n];

  //Propensidades de exponerse
  As[0] = beta*N95*(Sa+STa)*(phi1*(Pa+PTa+La+LTa)/(double)Na + mu*(Pb+PTb+Lb+LTb)/(double)Nb + (1-alpha)*phi1*(IAa+PTAa+LTAa)/(double)Na + (1-alpha)*mu*(IAb+PTAb+LTAb)/(double)Nb + eta*prev);
  As[1] = beta*(Sb+STb)*(mu*N95*(Pa+PTa+La+LTa)/(double)Na + chi*TBQ*(Pb+PTb+Lb+LTb)/(double)Nb + (1-alpha)*mu*N95*(IAa+PTAa+LTAa)/(double)Na + (1-alpha)*chi*TBQ*(IAb+PTAb+LTAb)/(double)Nb);

  //Propensidades de ser presintomático
  As[2] = USDe*(Ea+ETa+EAa);
  As[3] = USDe*(Eb+ETb+EAb);

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
  double prom = 165.47, sigma = 29.24, B = beta*N95*(Sa+STa)*eta*prev;
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

  //Escojo a la persona que contagia, si hay infección
  int conta = 0;
  if((int)index == 0){conta = who_infected(Val[5], Vba[5], Val[6], Vba[6], Val[7], Vba[7], Val[8], Vba[8], Val[9], Vba[9], Val[10], Vba[10], Val[11], Vba[11], phi1, mu, ran, 1, prev, t, N95, N95);}
  else if((int)index == 1){conta = who_infected(Val[5], Vba[5], Val[6], Vba[6], Val[7], Vba[7], Val[8], Vba[8], Val[9], Vba[9], Val[10], Vba[10], Val[11], Vba[11], phi1, mu, ran, 0, prev, t, N95, TBQ);}

  //Creo el vector resultados
  std::vector<double> result(3);
  result[0] = tau; //Tiempo en el que sucede la reacción
  result[1] = index; //Número de la reacción que sucede
  result[2] = conta; //Agente que hizo el contagio

  return result;
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


int who_infected(grupo &Pa, grupo &Pb, grupo &PTa, grupo &PTb, grupo &PTAa, grupo &PTAb, grupo &La, grupo &Lb, grupo &LTa, grupo &LTb, grupo &LTAa, grupo &LTAb, grupo &IAa, grupo &IAb, double cons1, double cons2, Crandom &ran, int alti, double prev, double t, double TBa, double TBb){
  //Parámetros Gaussiana
  double prom = 165.47, sigma = 29.24;

  //Calculo la propensidad de cada grupo
  double num[5];
  num[0] = cons1*TBa*(Pa.size() + PTa.size() + La.size() + LTa.size())/(double)Na;
  num[1] = cons2*TBb*(Pb.size() + PTb.size() + Lb.size() + LTb.size())/(double)Nb;
  num[2] = (1-alpha)*cons1*TBa*(IAa.size() + PTAa.size() + LTAa.size())/(double)Na;
  num[3] = (1-alpha)*cons2*TBb*(IAb.size() + PTAb.size() + LTAb.size())/(double)Nb;
  num[4] = alti*eta*N95*prev*std::exp(-(prom-t)*(prom-t)/(2*sigma*sigma));

  //Hallo el individuo que contagia
  grupo aux;
  double num2 = ran.r()*(num[0] + num[1] + num[2] + num[3] + num[4]);
  if(num2 < num[0]){return selection_infectious(Pa, PTa, La, LTa, ran);}
  else if(num2 < num[0] + num[1]){return selection_infectious(Pb, PTb, Lb, LTb, ran) + Na;}
  else if(num2 < num[0] + num[1] + num[2]){return selection_infectious(IAa, PTAa, LTAa, aux, ran);}
  else if(num2 < num[0] + num[1] + num[2] + num[3]){return selection_infectious(IAb, PTAb, LTAb, aux, ran) + Na;}
  else{return -1;}
}


int selection_infectious(grupo &Ga, grupo &Gb, grupo &Gc, grupo &Gd, Crandom &ran){
  double num = ran.r()*(Ga.size() + Gb.size() + Gc.size() + Gd.size());
  int ind, agent;
  if(num < Ga.size()){ind = (int)(ran.r()*Ga.size());    agent = Ga[ind];}
  else if(num < Ga.size() + Gb.size()){ind = (int)(ran.r()*Gb.size());    agent = Gb[ind];}
  else if(num < Ga.size() + Gb.size() + Gc.size()){ind = (int)(ran.r()*Gc.size());    agent = Gc[ind];}
  else{ind = (int)(ran.r()*Gd.size());    agent = Gd[ind];}
  return agent;
}
