#include <math.h>
#include <algorithm>
#include "dynamics.h"

std::vector<double> contagio(double Sa, double Sb, double Ea, double Eb, double Pa, double Pb, double PTAa, double PTAb, double La, double Lb, double LTAa, double LTAb, double LAa, double LAb, double IAa, double IAb, double Na, double Nb, double prev, Crandom &ran){

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


  //Arreglo de los tiempos
  double ts[n];

  //Hallo el tiempo de cada propensidad
  for(unsigned int i=0; i<n; i++){
    if(As[i] != 0){ts[i] = (-1/As[i])*std::log(ran.r());}
    else{ts[i] = 1e6;}
  }

  //Hallo el tiempo mínimo y la propensidad mínima
  double *pointer = std::min_element(ts, ts+n);
  std::vector<double> result(2);
  result[0] = *pointer; //Tiempo mínimo
  result[1] = std::distance(ts, pointer); //Índice de la reacción

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
