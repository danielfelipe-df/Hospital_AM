#include <iostream>
#include <fstream>
#include "scenarios.h"
#include "myconstants.h"


int main(void)
{
  double iota = 0.6;
  double mu = 0.1, chi = 0.7, phi1 = 0.6;
  double eta = 0.7, lambda = 0.7;
  unsigned int trace_net = 6;

  double xi[3] = {0.9, 0.8, 0.65}, theta = 0.0;
  double N95 = 0.84, TBQ = 0.84, HW = 0.69, SDP = 0.89;
  double alpha = 0.8;
  unsigned int trace = 3;
  unsigned int nu = 2, delta = 7;
  bool aislev = true;

  std::string name;
  for(size_t i=0; i<3; i++){
    name = "Heat-map_Si-rastreo_" + std::to_string((int)(xi[i]*100)) + ".csv";
    main_scenario(name, xi[i], theta, iota, N95, TBQ, HW, SDP, alpha, mu, chi, phi1, eta, lambda, trace, trace_net, nu, delta, N_gaussM, A_gaussM, Mu_gaussM, Sigma_gaussM, aislev);
  }
  
  return 0;
}

