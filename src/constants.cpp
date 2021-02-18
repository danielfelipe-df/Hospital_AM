#include <vector>
#include <constants.h>
#include <other_CSV.h>

constants::constants(std::string name){
  std::vector<std::vector<std::wstring> > data;

  csv_ws(data, name);

  this->xi = std::stod(data[0][2]); //Xi
  this->theta = std::stod(data[1][2]); //Theta
  this->iota = std::stod(data[2][2]); //Iota
  
  this->N95 = std::stod(data[3][2]); //N95
  this->TBQ = std::stod(data[4][2]); //TBQ
  this->HW = std::stod(data[5][2]); //HW
  this->SDP = std::stod(data[6][2]); //SDP
  
  this->alpha = std::stod(data[7][2]); //Alpha
  this->mu = std::stod(data[8][2]); //Mu
  this->chi = std::stod(data[9][2]); //Chi
  this->phi1 = std::stod(data[10][2]); //Phi1
  this->eta = std::stod(data[11][2]); //Eta
  this->lambda = std::stod(data[12][2]); //Lambda
  
  this->trace = std::stoi(data[13][2]); //Trace
  this->trace_net = std::stoi(data[14][2]); //Trace_net

  this->nu = std::stoi(data[15][2]); //Nu
  this->delta = std::stoi(data[16][2]); //Delta

  this->N_gauss = std::stoi(data[17][2]); //N_gauss
  this->A_gauss = new double[N_gauss];
  this->Mu_gauss = new double[N_gauss];
  this->Sigma_gauss = new double[N_gauss];
  for(size_t i=0; i<N_gauss; i++){
    this->A_gauss[i] = std::stod(data[18][2+i]); //A_gauss
    this->Mu_gauss[i] = std::stod(data[19][2+i]); //Mu_gauss
    this->Sigma_gauss[i] = std::stod(data[20][2+i]); //Sigma_gauss
  }

  this->N_betaL = std::stod(data[21][2]);
  this->m_betaL = new double[N_betaL];
  this->b_betaL = new double[N_betaL];
  this->lim_betaL = new double[N_betaL+1];
  for(size_t i=0; i<N_betaL; i++){
    this->m_betaL[i] = std::stod(data[22][2+i]); //m_betaL
    this->b_betaL[i] = std::stod(data[23][2+i]); //b_betaL
    this->lim_betaL[i] = std::stod(data[24][2+i]); //lim_betaL
  }
  this->lim_betaL[N_betaL] = std::stod(data[24][2+N_betaL]); //lim_betaL

  this->N_betaF = std::stod(data[25][2]); //N_betaF
  this->m_betaF = new double[N_betaF];
  this->b_betaF = new double[N_betaF];
  this->lim_betaF = new double[N_betaF+1];
  for(size_t i=0; i<N_betaF; i++){
    this->m_betaF[i] = std::stod(data[26][2+i]); //m_betaF
    this->b_betaF[i] = std::stod(data[27][2+i]); //b_betaF
    this->lim_betaF[i] = std::stod(data[28][2+i]); //lim_betaF
  }
  this->lim_betaF[N_betaF] = std::stod(data[28][2+N_betaF]); //lim_betaF

  this->AisLev = ((std::stoi(data[29][2]) == 1) ? true : false); //AisLev
}
