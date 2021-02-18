#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <string>

class constants{
public:

  /* Attributes */
  double xi; //Xi
  double theta; //Theta
  double iota; //Iota

  double N95; //N95
  double TBQ; //TBQ
  double HW; //HW
  double SDP; //SDP

  double alpha; //Alpha
  double mu; //Mu
  double chi; //Chi
  double phi1; //Phi1
  double eta; //Eta
  double lambda; //Lambda

  unsigned int trace; //Trace
  unsigned int trace_net; //Trace_net

  unsigned int nu; //Nu
  unsigned int delta; //Delta

  size_t N_gauss; //N_gauss
  double *A_gauss; //A_gauss
  double *Mu_gauss; //Mu_gauss
  double *Sigma_gauss; //Sigma_gauss

  size_t N_betaL; //N_betaL
  double *m_betaL; //m_betaL
  double *b_betaL; //b_betaL
  double *lim_betaL; //lim_betaL

  size_t N_betaF; //N_betaF
  double *m_betaF; //m_betaF
  double *b_betaF; //b_betaF
  double *lim_betaF; //lim_betaF

  bool AisLev; //AisLev


  /* Constructor */
  constants(std::string name);

};

#endif
