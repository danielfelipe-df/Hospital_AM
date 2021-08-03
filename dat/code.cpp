#include <iostream>
#include <fstream>
#include "scenarios.h"
#include "myconstants.h"


int main(void)
{
  double iota = 0.6;
  double mu = 0.1, chi = 0.7, phi1 = 0.6;
  double eta = 0.7, lambda = 0.186;
  unsigned int trace_net = 6;

  double xi = 0.9, theta = 0.5;
  double N95 = 0.84, TBQ = 0.84, HW = 0.69, SDP = 0.89;
  double alpha = 0.8;
  unsigned int trace = 0;
  unsigned int nu = 0, delta = 5;
  bool aislev = true;

  std::string name = "template.csv";
  main_scenario(name, xi, theta, iota, N95, TBQ, HW, SDP, alpha, mu, chi, phi1, eta, lambda, trace, trace_net, nu, delta, N_gaussM, A_gaussM, Mu_gaussM, Sigma_gaussM, aislev);

  return 0;
}

