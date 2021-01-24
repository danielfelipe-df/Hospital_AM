#ifndef RANDOM64_H
#define RANDOM64_H

//Constantes del generador aleatorio
class Crandom{
  unsigned long long u,v,w;

public:
  Crandom(unsigned long long j);
  unsigned long long int64();
  double r() {return 5.42101086242752217E-20 * int64();}
  unsigned int int32(){return (unsigned int) int64();};
  double exponencial(float tau);
  double gauss(float mu,float sigma);
};

#endif
