#include <math.h>
#include <algorithm>
#include "dynamics.h"

std::vector<double> contagio(double Sa, double Sb, double Ea, double Eb, double Pa, double Pb, double PTAa, double PTAb, double La, double Lb, double LTAa, double LTAb, double LAa, double LAb, double IAa, double IAb, double Na, double Nb, double prev, Crandom &ran, double t){

  //Número de proponsidades
  const int n = 22;

  //Creo el arreglo de las propensidades
  double As[n];

  //Propensidades de exponerse
  As[0] = beta*Sa*((Pa+La)/Na + mu*(Pb+Lb)/Nb + (1-alpha)*(IAa+PTAa+LTAa+LAa)/Na + (1-alpha)*mu*(IAb+PTAb+LTAb+LAb)/Nb + eta*prev);
  As[1] = beta*Sb*(mu*(Pa+La)/Na + chi*(Pb+Lb)/Nb + (1-alpha)*mu*(IAa+PTAa+LTAa+LAa)/Na + (1-alpha)*chi*(IAb+PTAb+LTAb+LAb)/Nb);

  //Propensidades de ser presintomático
  As[2] = USDe*Ea;
  As[3] = USDe*Eb;

  //Propensidades de ser leve
  As[4] = USDpl*(1-kappa)*psi*Pa;
  As[5] = USDpl*(1-kappa)*psi*Pb;

  //Propensidad de ser infeccioso grave y aislarse inmediatamente
  As[6] = USDpg*(1-kappa)*(1-psi)*Pa;
  As[7] = USDpg*(1-kappa)*(1-psi)*Pb;

  //Propensidad de ser aislado siendo leve (continuo)
  As[8] = (iota*xi/Tt)*La;
  As[9] = (iota*xi/Tt)*Lb;
  
  //Propensidades de recuperarse siendo asintomático
  As[10] = USDplil*kappa*Pa;
  As[11] = USDplil*kappa*Pb;
  
  //Propensidades de recuperarse siendo leve
  As[12] = USDil*(1-iota*xi)*La;
  As[13] = USDil*(1-iota*xi)*Lb;

  //Propensidades de recuperarse siendo asintomático aislado
  As[14] = USrho*PTAa;
  As[15] = USrho*PTAb;

  //Propensidades de recuperarse siendo leve aislado
  As[16] = USlambda*LAa;
  As[17] = USlambda*LAb;

  //Propensidades de recuperarse siendo leve testeado y aislado
  As[18] = UStau*LTAa;
  As[19] = UStau*LTAb;
  
  //Propensidades de recuperarse siendo infeccioso grave aislado
  As[20] = USepsilon*IAa;
  As[21] = USepsilon*IAb;
  

  //Hallo el tiempo en el que va a pasar la siguiente reacción con el método MDM
  double prom = 230.0, sigma = 55.0, A = beta*Sa*eta*prev, B = 0;
  double tau = 0, index = 0;
  for(unsigned int i=0; i<n; i++){B += As[i];}
  tau = biseccion(A, prom, sigma, t, B-A, -std::log(ran.r()));

  if(tau < 1e6){//Si el tiempo en el que pasa la reacción es menor al máximo (1e6), entonces es porque la reacción si sucede
    //Hallo el vector que me guarda la propensidad acumulada en orden
    double cumulative[n];
    cumulative[0] = (As[0]-A) + A*std::exp(-(prom-t)*(prom-t)/(2*sigma*sigma));
    for(unsigned int i=1; i<n; i++){cumulative[i] = cumulative[i-1] + As[i];}

    //Escojo la reacción a escoger
    double lim = ran.r()*cumulative[n-1];
    for(index = 0; index<n; index++){
      if(lim < cumulative[(int)index]){break;}
    }
  }

  //Creo el vector resultados
  std::vector<double> result(2);
  result[0] = tau; //Tiempo en el que sucede la reacción
  result[1] = index; //Número de la reacción que sucede

  return result;
}


void mother_reaction(grupo &Out, grupo &In, Crandom &ran, trabajadores *family, int typeout, int typein){
  int index = (int)(ran.r()*Out.size());
  int agent = Out[index];
  In.push_back(agent);
  Out.erase(Out.begin() + index);
  family[agent].kind[typeout] = false;
  family[agent].kind[typein] = true;
}


void massive_reaction(grupo &S, grupo &E, grupo &P, grupo &PTA, grupo &L, grupo &LTA, grupo &R, Crandom &ran, trabajadores *family){
  int N = S.size() + E.size() + P.size() + L.size() + R.size();
  int num = (int)(ran.r()*N);

  if(ran.r() < xi){
    if(num < P.size()){mother_reaction(P, PTA, ran, family, 2, 3);}
    else if(num < P.size() + L.size()){mother_reaction(L, LTA, ran, family, 4, 5);}
  }
}


double biseccion(double A, double prom, double sigma, double t, double B, double ranr){
  double m,fa,fm;
  double a = 0, b = 1e3, eps = 1e-7, nmax = 100, n=0;
  fa = B*a + function(A, prom, sigma, t, a) - ranr;
  
  while(b-a>eps && n<nmax){
    m = (a+b)/2;
    fm = B*m + function(A, prom, sigma, t, m) - ranr;
    if(fa*fm<0){b = m;}
    else{a = m; fa = fm;}
    n++;
  }

  if(n==nmax || m > (1e3)-(1e-4)){return 1e6;}
  else{return (a+b)/2;}
}


double function(double A, double prom, double sigma, double t, double tau){
  double A2 = std::sqrt(M_PI/2)*A*sigma, USsigma2 = 1.0/(std::sqrt(2)*sigma);
  return A2*(erf((prom-t)*USsigma2) - erf((prom-t-tau)*USsigma2));
}
