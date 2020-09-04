#include <math.h>
#include <algorithm>
#include "dynamics.h"

std::vector<double> contagio(double Sa, double Sb, double Ea, double Eb, double Pa, double Pb, double PTAa, double PTAb, double La, double Lb, double LTAa, double LTAb, double LAa, double LAb, double IAa, double IAb, double Na, double Nb, double prev, Crandom &ran, double r){

  //Número de proponsidades
  const int n = 26;

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

  //Propensidad de ser testeado y aislado siendo pre-sintomático (masivo)
  As[8] = (r*theta*xi/Tt)*kappa*Pa;
  As[9] = (r*theta*xi/Tt)*kappa*Pb;

  //Propensidad de ser aislado siendo leve (continuo)
  As[10] = (iota*xi/Tt)*La;
  As[11] = (iota*xi/Tt)*Lb;

  //Propensidad de ser testeado y aislado siendo leve (masivo)
  As[12] = (r*theta*xi*(1-iota*xi)/Tt)*La;
  As[13] = (r*theta*xi*(1-iota*xi)/Tt)*Lb;

  //Propensidades de recuperarse siendo asintomático
  As[14] = USDplil*(1-r*theta*xi)*kappa*Pa;
  As[15] = USDplil*(1-r*theta*xi)*kappa*Pb;

  //Propensidades de recuperarse siendo leve
  As[16] = USDil*(1-r*theta*xi)*(1-iota*xi)*La;
  As[17] = USDil*(1-r*theta*xi)*(1-iota*xi)*Lb;

  //Propensidades de recuperarse siendo asintomático aislado
  As[18] = USrho*PTAa;
  As[19] = USrho*PTAb;

  //Propensidades de recuperarse siendo leve aislado
  As[20] = USlambda*LAa;
  As[21] = USlambda*LAb;

  //Propensidades de recuperarse siendo leve testeado y aislado
  As[22] = UStau*LTAa;
  As[23] = UStau*LTAb;
  
  //Propensidades de recuperarse siendo infeccioso grave aislado
  As[24] = USepsilon*IAa;
  As[25] = USepsilon*IAb;

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


void mother_reaction(grupo &Out, grupo &In, Crandom &ran){
  int index = (int)(ran.r()*Out.size());
  In.push_back(Out[index]);
  Out.erase(Out.begin() + index);
}

