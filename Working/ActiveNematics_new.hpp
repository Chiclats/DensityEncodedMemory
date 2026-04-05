/*--------------------
ver 260224
--------------------*/

#ifndef ActiveNematics_new_hpp
#define ActiveNematics_new_hpp

#include<omp.h>
#include"./BasicTools_new.hpp"



namespace ActiveNematics_2D{

  struct Particle{
    
    double x,y,vx,vy;
    
    unsigned int index;
    // to store the index of the particle in the original PG

    unsigned int box_index;
    // to store the index of the box to which the particle belongs, for updating MeshedIndexMat
  };

}

#endif//ActiveNematics_new_hpp
