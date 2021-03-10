#ifndef SCENARIOS_H
#define SCENARIOS_H

#include <string>


// Main function for the scenarios
void main_scenario(std::string name, double xi, double theta, double iota, double N95, double TBQ, double HW, double SDP, double alpha, double mu, double chi, double phi1, double eta, double lambda, unsigned int trace, unsigned int trace_net, unsigned int nu, unsigned int delta, size_t N_gauss, const double *A_gauss, const double *Mu_gauss, const double *Sigma_gauss, bool aislev);

// Function for 'Control Nulo'
void CNada(std::string name, double iota, double mu, double chi, double phi1, double eta, double lambda, unsigned int trace_net, unsigned int prevtype);


// Function for 'Control Medio-Bajo'
void CMedio_Bajo(std::string name, double iota, double mu, double chi, double phi1, double eta, double lambda, unsigned int trace_net, unsigned int prevtype);


// Function for 'Control Medio-Alto'
void CMedio_Alto(std::string name, double iota, double mu, double chi, double phi1, double eta, double lambda, unsigned int trace_net, unsigned int prevtype);


// Function for 'Control Total'
void CTotal(std::string name, double iota, double mu, double chi, double phi1, double eta, double lambda, unsigned int trace_net, unsigned int prevtype);

#endif
