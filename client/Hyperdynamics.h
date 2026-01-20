#ifndef HYPERDYNAMICS_H
#define HYPERDYNAMICS_H

#include "HelperFunctions.h"
#include "Parameters.h"
#include "Eigen.h"
#include "Matter.h"
 

class Hyperdynamics {
        public:
          Hyperdynamics(Matter *matt, Parameters *params);
          virtual ~Hyperdynamics();
          virtual double boost() = 0;
          virtual void initialize() = 0;
            
          static const char NONE[];
          static const char BOND_BOOST[];
          static const char RIDGE_BASED[];
          
};

#endif
