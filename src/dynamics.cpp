#include <math.h>
#include <dynamics.h>


std::vector<double> contagio(std::map<std::string, grupo> &Val, std::map<std::string, grupo> &Vba, Crandom &ran, double t, double* tj){
  //Calculo los tamaños de cada vector
  double Sa, Sb, STa, STb, SAa, SAb, Ea, Eb, ETa, ETb, EAa, EAb, Pa, Pb, PTa, PTb, PTAa, PTAb, La, Lb, LTa, LTb, LTAa, LTAb, IAa, IAb;
  Sa = Val["SUS"].size();  STa = Val["SUST"].size();  SAa = Val["SUSA"].size();
  Sb = Vba["SUS"].size();  STb = Vba["SUST"].size();  SAb = Vba["SUSA"].size();
  Ea = Val["EXP"].size();  ETa = Val["EXPT"].size();  EAa = Val["EXPA"].size();
  Eb = Vba["EXP"].size();  ETb = Vba["EXPT"].size();  EAb = Val["EXPA"].size();
  Pa = Val["PRE"].size();  PTa = Val["PRET"].size();  PTAa = Val["PREA"].size();
  Pb = Vba["PRE"].size();  PTb = Vba["PRET"].size();  PTAb = Vba["PREA"].size();
  La = Val["MSYM"].size();  LTa = Val["MSYMT"].size();  LTAa = Val["MSYMA"].size();
  Lb = Vba["MSYM"].size();  LTb = Vba["MSYMT"].size();  LTAb = Vba["MSYMA"].size();
  IAa = Val["SSYMA"].size();
  IAb = Vba["SSYMA"].size();

  //Número de propensidades
  const unsigned int n = 14;

  //Creo el arreglo de las propensidades
  double As[n];

  //Propensidades de exponerse. (El alto no tiene la parte la gaussiana. Esa está implementada en la función bisección.)
  As[0] = MyCons.HW*MyCons.SDP*MyCons.N95*(Sa+STa+(1-MyCons.alpha)*SAa)*(MyCons.phi1*(Pa+PTa+La+LTa)/(double)Na + MyCons.mu*(Pb+PTb+Lb+LTb)/(double)Nb + (1-MyCons.alpha)*MyCons.phi1*(IAa+PTAa+LTAa)/(double)Na + (1-MyCons.alpha)*MyCons.mu*(IAb+PTAb+LTAb)/(double)Nb);
  As[1] = MyCons.HW*MyCons.SDP*(Sb+STb+(1-MyCons.alpha)*SAb)*(MyCons.mu*MyCons.N95*(Pa+PTa+La+LTa)/(double)Na + MyCons.chi*MyCons.TBQ*(Pb+PTb+Lb+LTb)/(double)Nb + (1-MyCons.alpha)*MyCons.mu*MyCons.N95*(IAa+PTAa+LTAa)/(double)Na + (1-MyCons.alpha)*MyCons.chi*MyCons.TBQ*(IAb+PTAb+LTAb)/(double)Nb);

  //Propensidades de ser presintomático
  As[2] = USDe*(Ea+ETa+EAa);
  As[3] = USDe*(Eb+ETb+EAb);

  //Propensidades de ser leve
  As[4] = USDpl*(1-kappa)*psi*(Pa+PTa+PTAa);
  As[5] = USDpl*(1-kappa)*psi*(Pb+PTb+PTAb);

  //Propensidad de ser infeccioso grave y aislarse inmediatamente
  As[6] = USDpg*(1-kappa)*(1-psi)*(Pa+PTa+PTAa);
  As[7] = USDpg*(1-kappa)*(1-psi)*(Pb+PTb+PTAb);

  //Propensidad de recuperarse siendo pre-sintomático
  As[8] = USDplil*kappa*(Pa+PTa+PTAa);
  As[9] = USDplil*kappa*(Pb+PTb+PTAb);

  //Propensidad de recuperarse siendo leve
  As[10] = USDil*(La+LTa+LTAa);
  As[11] = USDil*(Lb+LTb+LTAb);

  //Propensidad de recuperarse siendo infeccioso grave
  As[12] = USDig*IAa;
  As[13] = USDig*IAb;

  //Las propensidades que sean cero las pongo con un mínimo
  //for(unsigned int i=2; i<n; i++){if(As[i] <= 0.0){As[i] = 1e-12;}}

  //Creo las distribuciones de cada propensidad
  std::vector<lognormal_d> dist;
  //for(unsigned int i=2; i<n; i++){lognormal_d my_dist(1.5, 1.0/As[i]);    dist.push_back(my_dist);}

  //Hallo el tiempo en el que va a pasar la siguiente reacción con el método NMGA
  double Ba1 = MyCons.HW*MyCons.SDP*MyCons.N95*(Sa+STa+(1-MyCons.alpha)*SAa)*MyCons.eta, Ba2 = MyCons.HW*MyCons.SDP*(Sa+STa+(1-MyCons.alpha)*SAa)*MyCons.lambda, Bb = MyCons.HW*MyCons.SDP*(Sb+STb+(1-MyCons.alpha)*SAb)*MyCons.lambda;
  double tau = 0, index = 0;
  tau = biseccion(As, t, Ba1, Ba2, Bb, std::log(ran.r()), tj, n, dist); //Aquí ya está implementada las gaussianas

  if(tau < 1e6){//Si el tiempo en el que pasa la reacción es menor al máximo (1e6), entonces es porque la reacción si sucede
    //Hallo el vector que me guarda la propensidad acumulada en orden
    double cumulative[n];
    double value = 0;
    for(size_t i=0; i<MyCons.N_gauss; i++){value += function_gauss(t, MyCons.A_gauss[i], MyCons.Mu_gauss[i], MyCons.Sigma_gauss[i]);}
    cumulative[0] = function_beta(t, MyCons.lim_betaL, MyCons.m_betaL, MyCons.b_betaL, MyCons.N_betaL)*As[0] + Ba1*value + Ba2*value;
    cumulative[1] = cumulative[0] + function_beta(t, MyCons.lim_betaF, MyCons.m_betaF, MyCons.b_betaF, MyCons.N_betaF)*As[1] + Bb*value;
    for(unsigned int i=2; i<n; i++){cumulative[i] = cumulative[i-1] + As[i];}//(pdf(dist[i-2], tj[i]+tau)/cdf(complement(dist[i-2], tj[i]+tau)));}

    //Escojo la reacción a escoger
    double lim = (double)(ran.r()*cumulative[n-1]);
    for(index = 0; index<(n-1); index++){
      if(lim < cumulative[(int)index]){break;}
    }
    if(lim >= cumulative[n-2]){index = n-1;}
  }

  //Le sumo el tau a todos los tiempos tj
  for(unsigned int i=0; i<n; i++){tj[i] += tau;}

  //Reinicio el tiempo tj de la reacción que se escogió
  tj[(int)index] = 0.0;

  //Escojo a la persona que contagia, si hay infección
  int conta = 0;
  if((int)index == 0){conta = who_infected(Val["PRE"], Vba["PRE"], Val["PRET"], Vba["PRET"], Val["PREA"], Vba["PREA"], Val["MSYM"], Vba["MSYM"], Val["MSYMT"], Vba["MSYMT"], Val["MSYMA"], Vba["MSYMA"], Val["SSYMA"], Vba["SSYMA"], MyCons.phi1, MyCons.mu, ran, 1, t, MyCons.N95, MyCons.N95);}
  else if((int)index == 1){conta = who_infected(Val["PRE"], Vba["PRE"], Val["PRET"], Vba["PRET"], Val["PREA"], Vba["PREA"], Val["MSYM"], Vba["MSYM"], Val["MSYMT"], Vba["MSYMT"], Val["MSYMA"], Vba["MSYMA"], Val["SSYMA"], Vba["SSYMA"], MyCons.phi1, MyCons.mu, ran, 0, t, MyCons.N95, MyCons.TBQ);}

  //Creo el vector resultados
  std::vector<double> result(3);
  result[0] = tau; //Tiempo en el que sucede la reacción
  result[1] = index; //Número de la reacción que sucede
  result[2] = conta; //Agente que hizo el contagio

  return result;
}


double biseccion(double* A, double t, double Ba1, double Ba2, double Bb, double ranr, double* tj, int n, std::vector<lognormal_d> &dist){
  double m,fa,fm;
  double lim = 1e3, min = 0.0;
  double a = min, b = lim, eps = 1e-7, pmax = 100, p=0;
  fa = phi(A, tj, n, Ba1, Ba2, Bb, a, t, dist) - ranr;

  while(b-a>eps && p<pmax){
    m = (a+b)/2;
    fm = phi(A, tj, n, Ba1, Ba2, Bb, m, t, dist) - ranr;
    if(fa*fm<0){b = m;}
    else{a = m; fa = fm;}
    p++;
  }

  if(p==pmax || m > lim-(1e-4)){return 1e6;}
  else{return (a+b)/2;}
}


double function(double A, double prom, double sigma, double t, double tau){
  double A2 = std::sqrt(M_PI/2)*A*sigma, USsigma2 = 1.0/(std::sqrt(2)*sigma);
  return A2*(erf((prom-t)*USsigma2) - erf((prom-t-tau)*USsigma2));
}


double phi(double* A, double* tj, unsigned int n, double Ba1, double Ba2, double Bb, double deltat, double t, std::vector<lognormal_d> &dist){
  double psi_num, psi_den;

  double gnumL1 = 0, gdenL1 = 0, gnumL2 = 0, gdenL2 = 0;
  double bnumL1 = 0, bdenL1 = 0, bnumL2 = 0, bdenL2 = 0;
  double gnumF1 = 0, gdenF1 = 0, gnumF2 = 0, gdenF2 = 0;
  double bnumF1 = 0, bdenF1 = 0, bnumF2 = 0, bdenF2 = 0;
  double t0 = t - std::floor(t);
  double diff = std::floor(t);

  main_aux_phi_function(t0, diff, gnumL2, bnumL2, t, MyCons.lim_betaL, MyCons.m_betaL, MyCons.b_betaL, MyCons.N_betaL);
  main_aux_phi_function(t0, diff, gdenL2, bdenL2, t, MyCons.lim_betaL, MyCons.m_betaL, MyCons.b_betaL, MyCons.N_betaL);
  main_aux_phi_function(t0, diff, gnumF2, bnumF2, t, MyCons.lim_betaF, MyCons.m_betaF, MyCons.b_betaF, MyCons.N_betaF);
  main_aux_phi_function(t0, diff, gdenF2, bdenF2, t, MyCons.lim_betaF, MyCons.m_betaF, MyCons.b_betaF, MyCons.N_betaF);

  main_aux_phi_function(tj[0]+deltat, diff, gnumL1, bnumL1, t, MyCons.lim_betaL, MyCons.m_betaL, MyCons.b_betaL, MyCons.N_betaL);
  main_aux_phi_function(tj[0], diff, gdenL1, bdenL1, t, MyCons.lim_betaL, MyCons.m_betaL, MyCons.b_betaL, MyCons.N_betaL);
  main_aux_phi_function(tj[0]+deltat, diff, gnumF1, bnumF1, t, MyCons.lim_betaF, MyCons.m_betaF, MyCons.b_betaF, MyCons.N_betaF);
  main_aux_phi_function(tj[0], diff, gdenF1, bdenF1, t, MyCons.lim_betaF, MyCons.m_betaF, MyCons.b_betaF, MyCons.N_betaF);

  psi_num = -A[0]*(bnumL1 - bnumL2) - Ba1*(gnumL1 - gnumL2) - Ba2*(gnumF1 - gnumF2);
  psi_den = -A[0]*(bdenL1 - bdenL2) - Ba1*(gdenL1 - gdenL2) - Ba2*(gdenF1 - gdenF2);

  gnumL1 = 0;  gdenL1 = 0;  bnumL1 = 0;  bdenL1 = 0;
  gnumF1 = 0;  gdenF1 = 0;  bnumF1 = 0;  bdenF1 = 0;
  main_aux_phi_function(tj[1]+deltat, diff, gnumL1, bnumL1, t, MyCons.lim_betaL, MyCons.m_betaL, MyCons.b_betaL, MyCons.N_betaL);
  main_aux_phi_function(tj[1], diff, gdenL1, bdenL1, t, MyCons.lim_betaL, MyCons.m_betaL, MyCons.b_betaL, MyCons.N_betaL);
  main_aux_phi_function(tj[1]+deltat, diff, gnumF1, bnumF1, t, MyCons.lim_betaF, MyCons.m_betaF, MyCons.b_betaF, MyCons.N_betaF);
  main_aux_phi_function(tj[1], diff, gdenF1, bdenF1, t, MyCons.lim_betaF, MyCons.m_betaF, MyCons.b_betaF, MyCons.N_betaF);

  psi_num += -A[1]*(bnumL1 - bnumL2) - Bb*(gnumF1 - gnumF2);
  psi_den += -A[1]*(bdenL1 - bdenL2) - Bb*(gdenF1 - gdenF2);

  for(unsigned int i=2; i<n; i++){
    //psi_num += std::log(cdf(complement(dist[i-2], tj[i]+deltat)));
    //psi_den += std::log(cdf(complement(dist[i-2], tj[i])));
    psi_num += -A[i]*(tj[i] + deltat);
    psi_den += -A[i]*tj[i];
  }

  return psi_num-psi_den;
}


int who_infected(grupo &Pa, grupo &Pb, grupo &PTa, grupo &PTb, grupo &PTAa, grupo &PTAb, grupo &La, grupo &Lb, grupo &LTa, grupo &LTb, grupo &LTAa, grupo &LTAb, grupo &IAa, grupo &IAb, double cons1, double cons2, Crandom &ran, int alti, double t, double TBa, double TBb){  
  //Reviso si la infección se realizó dentro del hospital o afuera
  if(t - std::floor(t) < MyCons.lim_betaL[1]){// Si fue adentro, entonces acoto las propensidades a las que pertenecen a la dinámica interna
    //Parámetros Gaussiana
    double value = 0;
    for(size_t i=0; i<MyCons.N_gauss; i++){value += function_gauss(t, MyCons.A_gauss[i], MyCons.Mu_gauss[i], MyCons.Sigma_gauss[i]);}
    
    //Calculo la propensidad de cada grupo
    double num[5];
    num[0] = cons1*TBa*(Pa.size() + PTa.size() + La.size() + LTa.size())/(double)Na;
    num[1] = cons2*TBb*(Pb.size() + PTb.size() + Lb.size() + LTb.size())/(double)Nb;
    num[2] = (1-MyCons.alpha)*cons1*TBa*(IAa.size() + PTAa.size() + LTAa.size())/(double)Na;
    num[3] = (1-MyCons.alpha)*cons2*TBb*(IAb.size() + PTAb.size() + LTAb.size())/(double)Nb;
    num[4] = alti*MyCons.eta*MyCons.N95*value;
    
    //Hallo el individuo que contagia
    grupo aux;
    double num2 = ran.r()*(num[0] + num[1] + num[2] + num[3] + num[4]);
    if(num2 < num[0]){return selection_infectious(Pa, PTa, La, LTa, ran);}
    else if(num2 < num[0] + num[1]){return selection_infectious(Pb, PTb, Lb, LTb, ran) + Na;}
    else if(num2 < num[0] + num[1] + num[2]){return selection_infectious(IAa, PTAa, LTAa, aux, ran);}
    else if(num2 < num[0] + num[1] + num[2] + num[3]){return selection_infectious(IAb, PTAb, LTAb, aux, ran) + Na;}
    else{return -1;}
  }
  else{// Si fue afuera, entonces retorno el número asociado a la infección afuera
    return -2;
  }
}


int selection_infectious(grupo &Ga, grupo &Gb, grupo &Gc, grupo &Gd, Crandom &ran){
  double num = ran.r()*(Ga.size() + Gb.size() + Gc.size() + Gd.size());
  int ind, agent;
  if(num < Ga.size()){ind = (int)(ran.r()*Ga.size());    agent = Ga[ind];}
  else if(num < Ga.size() + Gb.size()){ind = (int)(ran.r()*Gb.size());    agent = Gb[ind];}
  else if(num < Ga.size() + Gb.size() + Gc.size()){ind = (int)(ran.r()*Gc.size());    agent = Gc[ind];}
  else{ind = (int)(ran.r()*Gd.size());    agent = Gd[ind];}
  return agent;
}



double function_gauss(double x, double A, double mu, double sigma){
  return A*std::exp(-(x-mu)*(x-mu)/(2*sigma*sigma));
}


void aux_phi_function(double t0, double diff, double &gnum, double &bnum, const double* lim_beta, const double* m_beta, const double* b_beta, const size_t N_beta){
  for(size_t i=1; i<(N_beta+1); i++){//Hacemos un ciclo sobre los tiempos de cambio del beta (se omite el primero porque es 0.0)
    if(t0 < lim_beta[i]){ //Revisamos si t0 está en el intervalo
      //Si lo está entonces usamos el algoritmo de integración
      
      for(size_t j=0; j<(i-1); j++){//Hacemos el algoritmo de integración acumulado durante los intervalos anteriores
	bnum += int_beta(lim_beta[j], lim_beta[j+1], m_beta[j], b_beta[j]);
	for(size_t k=0; k<MyCons.N_gauss; k++){
	  gnum += int_beta_gauss(diff+lim_beta[j], diff+lim_beta[j+1], MyCons.Mu_gauss[k], MyCons.Sigma_gauss[k], MyCons.A_gauss[k], m_beta[j], b_beta[j]);
	}
      }

      //Realizamos la integración en el último intervalo hasta t0
      bnum += int_beta(lim_beta[i-1], t0, m_beta[i-1], b_beta[i-1]);
      for(size_t k=0; k<MyCons.N_gauss; k++){
	gnum += int_beta_gauss(diff+lim_beta[i-1], diff+t0, MyCons.Mu_gauss[k], MyCons.Sigma_gauss[k], MyCons.A_gauss[k], m_beta[i-1], b_beta[i-1]);
      }

      //Rompemos el ciclo porque ya no es necesario seguir
      break;
    }
  }
}


void main_aux_phi_function(double t0, double diff, double &gnum, double &bnum, double t, const double* lim_beta, const double* m_beta, const double* b_beta, const size_t N_beta){
  //Si t0 está en el intervalo de los tiempos del beta, entonces hago la integración hasta ese intervalo
  if(t0 < 1.0){
    aux_phi_function(t0, diff, gnum, bnum, lim_beta, m_beta, b_beta, N_beta);
  }
  //Si no, hago la integración sobre los intervalos de beta necesarios
  else{
    //Hago la integración sobre todos los intervalos de tiempo de beta completos
    while(diff < (t+t0-1)){
      for(size_t j=0; j<N_beta; j++){
	for(size_t k=0; k<MyCons.N_gauss; k++){
	  gnum += int_beta_gauss(diff+lim_beta[j], diff+lim_beta[j+1], MyCons.Mu_gauss[k], MyCons.Sigma_gauss[k], MyCons.A_gauss[k], m_beta[j], b_beta[j]);
	}
	bnum += int_beta(lim_beta[j], lim_beta[j+1], m_beta[j], b_beta[j]);
      }
      diff++;
    }

    //Hago la integración sobre el intervalo incompleto
    t0 -= std::floor(t0);
    aux_phi_function(t0, diff, gnum, bnum, lim_beta, m_beta, b_beta, N_beta);
  }
}


double int_beta_gauss(double t0, double t1, double prom, double sigma, double A, double m, double b){
  double A1 = m*A*sigma*sigma, A2 = std::sqrt(M_PI_2)*sigma*A*(m*t0 + b), Un2S = 1/(M_SQRT2*sigma);
  double x1 = (t1 - prom)*Un2S, x0 = (t0 - prom)*Un2S;
  //Retorno la integral de una gaussiana multiplicada por x
  return A1*(std::exp(-x0*x0) - std::exp(-x1*x1)) + A2*(std::erf(x1) - std::erf(x0));
}


double int_beta(double t0, double t1, double m, double b){
  //Retorno la integral de una ecuación lineal
  return 0.5*m*(t1*t1 - t0*t0) + b*(t1 - t0);
}


double function_beta(double x, const double* lim_beta, const double* m_beta, const double* b_beta, const size_t N_beta){
  double myx = x - std::floor(x);
  size_t i;
  //Retorno la ecuación lineal de beta evaluada en x
  for(i=0; i<(N_beta-1); i++){if(myx <= lim_beta[i+1]){return m_beta[i]*myx + b_beta[i];}}
  return m_beta[i]*myx + b_beta[i];
}


size_t index_beta(double x, const double* lim_beta, const size_t N_beta){
  double myx = x - std::floor(x);
  size_t i;
  //Retorno el índice del intervalo en el que está x
  for(i=0; i<(N_beta-1); i++){if(myx <= lim_beta[i+1]){return i;}}
  return i;
}
