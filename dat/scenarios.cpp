#include <fstream>
#include "myconstants.h"
#include "scenarios.h"

void main_scenario(std::string name, double xi, double theta, double iota, double N95, double TBQ, double HW, double SDP, double alpha, double mu, double chi, double phi1, double eta, double lambda, unsigned int trace, unsigned int trace_net, unsigned int nu, unsigned int delta, size_t N_gauss, const double *A_gauss, const double *Mu_gauss, const double *Sigma_gauss, bool aislev){
  std::ofstream fout(name);

  fout << "Sensibilidad de la prueba" << '\t' << "xi" << '\t' << xi << std::endl;
  fout << "Cobertura Testeo Masivo" << '\t' << "theta" << '\t' << theta << std::endl;
  fout << "Proporción sintomáticos leves testeados continuamente" << '\t' << "iota" << '\t' << iota << std::endl;

  fout << "Tapabocas N95" << '\t' << "N95" << '\t' << N95 << std::endl;
  fout << "Tapabocas quirúrgico" << '\t' << "TBQ" << '\t' << TBQ << std::endl;
  fout << "Lavado de manos" << '\t' << "HW" << '\t' << HW << std::endl;
  fout << "Distanciamiento físico" << '\t' << "SDP" << '\t' << SDP << std::endl;

  fout << "Adherencia al aislamiento" << '\t' << "alpha" << '\t' << alpha << std::endl;
  fout << "Contacto cruzado" << '\t' << "mu" << '\t' << mu << std::endl;
  fout << "Contacto bajo-bajo" << '\t' << "chi" << '\t' << chi << std::endl;
  fout << "Contacto alto-alto" << '\t' << "phi1" << '\t' << phi1 << std::endl;
  fout << "Contacto hospitalizados" << '\t' << "eta" << '\t' << eta << std::endl;
  fout << "Contacto familiares" << '\t' << "lambda" << '\t' << lambda << std::endl;

  fout << "Personas rastreos" << '\t' << "trace" << '\t' << trace << std::endl;
  fout << "Red de posibles rastreos" << '\t' << "trace_net" << '\t' << trace_net << std::endl;

  fout << "Días con testeo masivo" << '\t' << "nu" << '\t' << nu << std::endl;
  fout << "Días sin testeo masivo" << '\t' << "delta" << '\t' << delta << std::endl;

  fout << "Número de gaussianas en la prevalencia externa" << '\t' << "N_gauss" << '\t' << N_gauss << std::endl;
  fout << "Amplitud prevalencia externa" << '\t' << "A_gauss" << '\t';
  for(size_t i=0; i<(N_gauss-1); i++){fout << A_gauss[i] << '\t';}
  fout << A_gauss[N_gauss-1] << std::endl;
  fout << "Mu prevalencia externa" << '\t' << "Mu_gauss" << '\t';
  for(size_t i=0; i<(N_gauss-1); i++){fout << Mu_gauss[i] << '\t';}
  fout << Mu_gauss[N_gauss-1] << std::endl;
  fout << "Sigma prevalencia externa" << '\t' << "Sigma_gauss" << '\t';
  for(size_t i=0; i<(N_gauss-1); i++){fout << Sigma_gauss[i] << '\t';}
  fout << Sigma_gauss[N_gauss-1] << std::endl;

  fout << "Secciones lineales del beta laboral" << '\t' << "N_betaL" << '\t' << N_betaL << std::endl;
  fout << "Pendiente de la ecuación lineal del betaL" << '\t' << "m_betaL" << '\t';
  for(size_t i=0; i<(N_betaL-1); i++){fout << m_betaL[i] << '\t';}
  fout << m_betaL[N_betaL-1] << std::endl;
  fout << "Corte en Y de la ecuación lineal del betaL" << '\t' << "b_betaL" << '\t';
  for(size_t i=0; i<(N_betaL-1); i++){fout << b_betaL[i] << '\t';}
  fout << b_betaL[N_betaL-1] << std::endl;
  fout << "Tiempos entre cada una de las regiones del betaL" << '\t' << "lim_betaL" << '\t';
  for(size_t i=0; i<N_betaL; i++){fout << lim_betaL[i] << '\t';}
  fout << lim_betaL[N_betaL] << std::endl;

  fout << "Secciones lineales del beta familiar" << '\t' << "N_betaF" << '\t' << N_betaF << std::endl;
  fout << "Pendiente de la ecuación lineal del betaF" << '\t' << "m_betaF" << '\t';
  for(size_t i=0; i<(N_betaF-1); i++){fout << m_betaF[i] << '\t';}
  fout << m_betaF[N_betaF-1] << std::endl;
  fout << "Corte en Y de la ecuación lineal del betaF" << '\t' << "b_betaF" << '\t';
  for(size_t i=0; i<(N_betaF-1); i++){fout << b_betaF[i] << '\t';}
  fout << b_betaF[N_betaF-1] << std::endl;
  fout << "Tiempos entre cada una de las regiones del betaF" << '\t' << "lim_betaF" << '\t';
  for(size_t i=0; i<N_betaF; i++){fout << lim_betaF[i] << '\t';}
  fout << lim_betaF[N_betaF] << std::endl;

  fout << "Aislamiento de leves testeados" << '\t' << "aislev" << '\t' << aislev << std::endl;

  fout.close();
}


void CNada(std::string name, double iota, double mu, double chi, double phi1, double eta, double lambda, unsigned int trace_net, unsigned int prevtype){

  double xi = 0.0, theta = 0.0;
  double N95 = 1.0, TBQ = 1.0, HW = 1.0, SDP = 1.0;
  double alpha = 0.0;
  unsigned int trace = 0;
  unsigned int nu = 0, delta = 7;
  bool aislev = false;

  if(prevtype == 0){
    main_scenario(name, xi, theta, iota, N95, TBQ, HW, SDP, alpha, mu, chi, phi1, eta, lambda, trace, trace_net, nu, delta, N_gaussA, A_gaussA, Mu_gaussA, Sigma_gaussA, aislev);
  }
  else if(prevtype == 1){
    main_scenario(name, xi, theta, iota, N95, TBQ, HW, SDP, alpha, mu, chi, phi1, eta, lambda, trace, trace_net, nu, delta, N_gaussM, A_gaussM, Mu_gaussM, Sigma_gaussM, aislev);
  }
  else if(prevtype == 2){
    main_scenario(name, xi, theta, iota, N95, TBQ, HW, SDP, alpha, mu, chi, phi1, eta, lambda, trace, trace_net, nu, delta, N_gaussB, A_gaussB, Mu_gaussB, Sigma_gaussB, aislev);
  }

}


void CMedio_Bajo(std::string name, double iota, double mu, double chi, double phi1, double eta, double lambda, unsigned int trace_net, unsigned int prevtype){

  double xi = 0.9, theta = 0.0;
  double N95 = 0.84, TBQ = 0.84, HW = 0.69, SDP = 0.89;
  double alpha = 0.5;
  unsigned int trace = 0;
  unsigned int nu = 0, delta = 7;
  bool aislev = false;

  if(prevtype == 0){
    main_scenario(name, xi, theta, iota, N95, TBQ, HW, SDP, alpha, mu, chi, phi1, eta, lambda, trace, trace_net, nu, delta, N_gaussA, A_gaussA, Mu_gaussA, Sigma_gaussA, aislev);
  }
  else if(prevtype == 1){
    main_scenario(name, xi, theta, iota, N95, TBQ, HW, SDP, alpha, mu, chi, phi1, eta, lambda, trace, trace_net, nu, delta, N_gaussM, A_gaussM, Mu_gaussM, Sigma_gaussM, aislev);
  }
  else if(prevtype == 2){
    main_scenario(name, xi, theta, iota, N95, TBQ, HW, SDP, alpha, mu, chi, phi1, eta, lambda, trace, trace_net, nu, delta, N_gaussB, A_gaussB, Mu_gaussB, Sigma_gaussB, aislev);
  }

}


void CMedio_Alto(std::string name, double iota, double mu, double chi, double phi1, double eta, double lambda, unsigned int trace_net, unsigned int prevtype){

  double xi = 0.9, theta = 0.0;
  double N95 = 0.84, TBQ = 0.84, HW = 0.69, SDP = 0.89;
  double alpha = 0.8;
  unsigned int trace = 0;
  unsigned int nu = 0, delta = 7;
  bool aislev = true;

  if(prevtype == 0){
    main_scenario(name, xi, theta, iota, N95, TBQ, HW, SDP, alpha, mu, chi, phi1, eta, lambda, trace, trace_net, nu, delta, N_gaussA, A_gaussA, Mu_gaussA, Sigma_gaussA, aislev);
  }
  else if(prevtype == 1){
    main_scenario(name, xi, theta, iota, N95, TBQ, HW, SDP, alpha, mu, chi, phi1, eta, lambda, trace, trace_net, nu, delta, N_gaussM, A_gaussM, Mu_gaussM, Sigma_gaussM, aislev);
  }
  else if(prevtype == 2){
    main_scenario(name, xi, theta, iota, N95, TBQ, HW, SDP, alpha, mu, chi, phi1, eta, lambda, trace, trace_net, nu, delta, N_gaussB, A_gaussB, Mu_gaussB, Sigma_gaussB, aislev);
  }

}


void CTotal(std::string name, double iota, double mu, double chi, double phi1, double eta, double lambda, unsigned int trace_net, unsigned int prevtype){

  double xi = 0.9, theta = 0.5;
  double N95 = 0.31, TBQ = 0.84, HW = 0.69, SDP = 0.89;
  double alpha = 0.9;
  unsigned int trace = 3;
  unsigned int nu = 2, delta = 7;
  bool aislev = true;

  if(prevtype == 0){
    main_scenario(name, xi, theta, iota, N95, TBQ, HW, SDP, alpha, mu, chi, phi1, eta, lambda, trace, trace_net, nu, delta, N_gaussA, A_gaussA, Mu_gaussA, Sigma_gaussA, aislev);
  }
  else if(prevtype == 1){
    main_scenario(name, xi, theta, iota, N95, TBQ, HW, SDP, alpha, mu, chi, phi1, eta, lambda, trace, trace_net, nu, delta, N_gaussM, A_gaussM, Mu_gaussM, Sigma_gaussM, aislev);
  }
  else if(prevtype == 2){
    main_scenario(name, xi, theta, iota, N95, TBQ, HW, SDP, alpha, mu, chi, phi1, eta, lambda, trace, trace_net, nu, delta, N_gaussB, A_gaussB, Mu_gaussB, Sigma_gaussB, aislev);
  }

}
