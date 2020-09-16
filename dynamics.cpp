#include <math.h>
#include <algorithm>
#include "dynamics.h"

std::vector<double> contagio(double Sa, double Sb, double Ea, double Eb, double Pa, double Pb, double PTa, double PTb, double PTAa, double PTAb, double La, double Lb, double LTa, double LTb, double LTAa, double LTAb, double LAa, double LAb, double IAa, double IAb, double Na, double Nb, double prev, Crandom &ran, double t){

  //Número de proponsidades
  const int n = 14;

  //Creo el arreglo de las propensidades
  double As[n];

  //Propensidades de exponerse
  As[0] = beta*Sa*((Pa+PTa+La+LTa)/Na + mu*(Pb+PTb+Lb+LTb)/Nb + (1-alpha)*(IAa+PTAa+LTAa+LAa)/Na + (1-alpha)*mu*(IAb+PTAb+LTAb+LAb)/Nb + eta*prev);
  As[1] = beta*Sb*(mu*(Pa+PTa+La+LTa)/Na + chi*(Pb+PTb+Lb+LTb)/Nb + (1-alpha)*mu*(IAa+PTAa+LTAa+LAa)/Na + (1-alpha)*chi*(IAb+PTAb+LTAb+LAb)/Nb);

  //Propensidades de ser presintomático
  As[2] = USDe*Ea;
  As[3] = USDe*Eb;

  //Propensidades de ser leve
  As[4] = USDpl*(1-kappa)*psi*Pa;
  As[5] = USDpl*(1-kappa)*psi*Pb;

  //Propensidad de ser infeccioso grave y aislarse inmediatamente
  As[6] = USDpg*(1-kappa)*(1-psi)*Pa;
  As[7] = USDpg*(1-kappa)*(1-psi)*Pb;

  //Propensidad de recuperarse siendo pre-sintomático
  As[8] = USDplil*kappa*(Pa+PTa+PTAa);
  As[9] = USDplil*kappa*(Pb+PTb+PTAb);
  
  //Propensidad de recuperarse siendo leve
  As[10] = USDil*(La+LTa+LTAa+LAa);
  As[11] = USDil*(Lb+LTb+LTAb+LAb);

  //Propensidad de recuperarse siendo infeccioso grave
  As[12] = USDig*IAa;
  As[13] = USDig*IAa;


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
  family[agent].change(typein, typeout);
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


void reaction0(grupo &Sa, grupo &Sb, grupo &Ea, grupo &Eb, grupo &Pa, grupo &Pb, grupo &PTa, grupo &PTb, grupo &PTAa, grupo &PTAb, grupo &La, grupo &LTa, grupo &LTb, grupo &Lb, grupo &LTAa, grupo &LTAb, grupo &LAa, grupo &LAb, grupo &IAa, grupo &IAb, grupo &Ra, grupo &Rb, Crandom &ran, trabajadores *altos, trabajadores *bajos){
  mother_reaction(Sa, Ea, ran, altos, 0, 1);
}


void reaction1(grupo &Sa, grupo &Sb, grupo &Ea, grupo &Eb, grupo &Pa, grupo &Pb, grupo &PTa, grupo &PTb, grupo &PTAa, grupo &PTAb, grupo &La, grupo &LTa, grupo &LTb, grupo &Lb, grupo &LTAa, grupo &LTAb, grupo &LAa, grupo &LAb, grupo &IAa, grupo &IAb, grupo &Ra, grupo &Rb, Crandom &ran, trabajadores *altos, trabajadores *bajos){
  mother_reaction(Sb, Eb, ran, bajos, 0, 1);
}


void reaction2(grupo &Sa, grupo &Sb, grupo &Ea, grupo &Eb, grupo &Pa, grupo &Pb, grupo &PTa, grupo &PTb, grupo &PTAa, grupo &PTAb, grupo &La, grupo &LTa, grupo &LTb, grupo &Lb, grupo &LTAa, grupo &LTAb, grupo &LAa, grupo &LAb, grupo &IAa, grupo &IAb, grupo &Ra, grupo &Rb, Crandom &ran, trabajadores *altos, trabajadores *bajos){
  mother_reaction(Ea, Pa, ran, altos, 1, 2);
}


void reaction3(grupo &Sa, grupo &Sb, grupo &Ea, grupo &Eb, grupo &Pa, grupo &Pb, grupo &PTa, grupo &PTb, grupo &PTAa, grupo &PTAb, grupo &La, grupo &LTa, grupo &LTb, grupo &Lb, grupo &LTAa, grupo &LTAb, grupo &LAa, grupo &LAb, grupo &IAa, grupo &IAb, grupo &Ra, grupo &Rb, Crandom &ran, trabajadores *altos, trabajadores *bajos){
  mother_reaction(Eb, Pb, ran, bajos, 1, 2);
}


void reaction4(grupo &Sa, grupo &Sb, grupo &Ea, grupo &Eb, grupo &Pa, grupo &Pb, grupo &PTa, grupo &PTb, grupo &PTAa, grupo &PTAb, grupo &La, grupo &LTa, grupo &LTb, grupo &Lb, grupo &LTAa, grupo &LTAb, grupo &LAa, grupo &LAb, grupo &IAa, grupo &IAb, grupo &Ra, grupo &Rb, Crandom &ran, trabajadores *altos, trabajadores *bajos){
  mother_reaction(Pa, La, ran, altos, 2, 5);
}


void reaction5(grupo &Sa, grupo &Sb, grupo &Ea, grupo &Eb, grupo &Pa, grupo &Pb, grupo &PTa, grupo &PTb, grupo &PTAa, grupo &PTAb, grupo &La, grupo &LTa, grupo &LTb, grupo &Lb, grupo &LTAa, grupo &LTAb, grupo &LAa, grupo &LAb, grupo &IAa, grupo &IAb, grupo &Ra, grupo &Rb, Crandom &ran, trabajadores *altos, trabajadores *bajos){
  mother_reaction(Pb, Lb, ran, bajos, 2, 5);
}


void reaction6(grupo &Sa, grupo &Sb, grupo &Ea, grupo &Eb, grupo &Pa, grupo &Pb, grupo &PTa, grupo &PTb, grupo &PTAa, grupo &PTAb, grupo &La, grupo &LTa, grupo &LTb, grupo &Lb, grupo &LTAa, grupo &LTAb, grupo &LAa, grupo &LAb, grupo &IAa, grupo &IAb, grupo &Ra, grupo &Rb, Crandom &ran, trabajadores *altos, trabajadores *bajos){
  mother_reaction(Pa, IAa, ran, altos, 2, 9);
}


void reaction7(grupo &Sa, grupo &Sb, grupo &Ea, grupo &Eb, grupo &Pa, grupo &Pb, grupo &PTa, grupo &PTb, grupo &PTAa, grupo &PTAb, grupo &La, grupo &LTa, grupo &LTb, grupo &Lb, grupo &LTAa, grupo &LTAb, grupo &LAa, grupo &LAb, grupo &IAa, grupo &IAb, grupo &Ra, grupo &Rb, Crandom &ran, trabajadores *altos, trabajadores *bajos){
  mother_reaction(Pb, IAb, ran, bajos, 2, 9);
}


void reaction8(grupo &Sa, grupo &Sb, grupo &Ea, grupo &Eb, grupo &Pa, grupo &Pb, grupo &PTa, grupo &PTb, grupo &PTAa, grupo &PTAb, grupo &La, grupo &LTa, grupo &LTb, grupo &Lb, grupo &LTAa, grupo &LTAb, grupo &LAa, grupo &LAb, grupo &IAa, grupo &IAb, grupo &Ra, grupo &Rb, Crandom &ran, trabajadores *altos, trabajadores *bajos){
  double N = Pa.size() + PTa.size() + PTAa.size();
  double num = ran.r()*N;
  if(num < Pa.size()){mother_reaction(Pa, Ra, ran, altos, 2, 10);}
  else if(num < Pa.size() + PTa.size()){mother_reaction(PTa, Ra, ran, altos, 3, 10);}
  else{mother_reaction(PTAa, Ra, ran, altos, 4, 10);}
}


void reaction9(grupo &Sa, grupo &Sb, grupo &Ea, grupo &Eb, grupo &Pa, grupo &Pb, grupo &PTa, grupo &PTb, grupo &PTAa, grupo &PTAb, grupo &La, grupo &LTa, grupo &LTb, grupo &Lb, grupo &LTAa, grupo &LTAb, grupo &LAa, grupo &LAb, grupo &IAa, grupo &IAb, grupo &Ra, grupo &Rb, Crandom &ran, trabajadores *altos, trabajadores *bajos){
  double N = Pb.size() + PTb.size() + PTAb.size();
  double num = ran.r()*N;
  if(num < Pb.size()){mother_reaction(Pb, Rb, ran, bajos, 2, 10);}
  else if(num < Pb.size() + PTb.size()){mother_reaction(PTb, Rb, ran, bajos, 3, 10);}
  else{mother_reaction(PTAb, Rb, ran, bajos, 4, 10);}
}


void reaction10(grupo &Sa, grupo &Sb, grupo &Ea, grupo &Eb, grupo &Pa, grupo &Pb, grupo &PTa, grupo &PTb, grupo &PTAa, grupo &PTAb, grupo &La, grupo &LTa, grupo &LTb, grupo &Lb, grupo &LTAa, grupo &LTAb, grupo &LAa, grupo &LAb, grupo &IAa, grupo &IAb, grupo &Ra, grupo &Rb, Crandom &ran, trabajadores *altos, trabajadores *bajos){
  double N = La.size() + LTa.size() + LTAa.size() + LAa.size();
  double num = ran.r()*N;
  if(num < La.size()){mother_reaction(La, Ra, ran, altos, 5, 10);}
  else if(num < La.size() + LTa.size()){mother_reaction(LTa, Ra, ran, altos, 6, 10);}
  else if(num < La.size() + LTa.size() + LTAa.size()){mother_reaction(LTAa, Ra, ran, altos, 7, 10);}
  else{mother_reaction(LAa, Ra, ran, altos, 8, 10);}
}


void reaction11(grupo &Sa, grupo &Sb, grupo &Ea, grupo &Eb, grupo &Pa, grupo &Pb, grupo &PTa, grupo &PTb, grupo &PTAa, grupo &PTAb, grupo &La, grupo &LTa, grupo &LTb, grupo &Lb, grupo &LTAa, grupo &LTAb, grupo &LAa, grupo &LAb, grupo &IAa, grupo &IAb, grupo &Ra, grupo &Rb, Crandom &ran, trabajadores *altos, trabajadores *bajos){
  double N = Lb.size() + LTb.size() + LTAb.size() + LAb.size();
  double num = ran.r()*N;
  if(num < Lb.size()){mother_reaction(Lb, Rb, ran, bajos, 5, 10);}
  else if(num < Lb.size() + LTb.size()){mother_reaction(LTb, Rb, ran, bajos, 6, 10);}
  else if(num < Lb.size() + LTb.size() + LTAb.size()){mother_reaction(LTAb, Rb, ran, bajos, 7, 10);}
  else{mother_reaction(LAb, Rb, ran, bajos, 8, 10);}
}


void reaction12(grupo &Sa, grupo &Sb, grupo &Ea, grupo &Eb, grupo &Pa, grupo &Pb, grupo &PTa, grupo &PTb, grupo &PTAa, grupo &PTAb, grupo &La, grupo &LTa, grupo &LTb, grupo &Lb, grupo &LTAa, grupo &LTAb, grupo &LAa, grupo &LAb, grupo &IAa, grupo &IAb, grupo &Ra, grupo &Rb, Crandom &ran, trabajadores *altos, trabajadores *bajos){
  mother_reaction(IAa, Ra, ran, altos, 9, 10);
}


void reaction13(grupo &Sa, grupo &Sb, grupo &Ea, grupo &Eb, grupo &Pa, grupo &Pb, grupo &PTa, grupo &PTb, grupo &PTAa, grupo &PTAb, grupo &La, grupo &LTa, grupo &LTb, grupo &Lb, grupo &LTAa, grupo &LTAb, grupo &LAa, grupo &LAb, grupo &IAa, grupo &IAb, grupo &Ra, grupo &Rb, Crandom &ran, trabajadores *altos, trabajadores *bajos){
  mother_reaction(IAb, Rb, ran, bajos, 9, 10);
}
