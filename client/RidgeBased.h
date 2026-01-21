#ifndef RIDGEBASED_H
#define RIDGEBASED_H

#include "HelperFunctions.h"
#include "Hyperdynamics.h"
#include "Matter.h"
#include "Parameters.h"
#include "ImprovedDimer.h"

#include "Eigen.h"

class RidgeBased: public Hyperdynamics{
public:
    RidgeBased(Matter *matt, Parameters *params);
    ~RidgeBased() override;
    bool contains_hyperF;
    bool reachridge;
    bool checknow;
    double T;
    double dR;
    double phi_tol;
    double rotationsMax;
    double EBoostMax;
    double EBoost;
    double EBoostOld;
    double maxStep;
    double curvature;
    double bisect_tol;
    double Eridge;
    double biasPot;
    int checksteps;
    int maxForceCalls;
    int forceCalls;
    int steps;
    int objectCalls;
    int opt_forceCalls;
    int dimerCalls;
    Matter R0;
    Matter R0old;
    Matter Rmin;
    Matter Rdimer;
    ImprovedDimer *dimer;
    AtomMatrix hyperF, mode, Nguess;
    VectorXd N, Ftrans, Nridge, Fridge, V, F0;
    

    void initialize() override;
    void bisect(std::vector<Matter> traj_dimer , double tol);
    void step();
    double boost() override;
    double get_biasPot();
    double search(double minForce, bool quite, int maxForceCalls, int interval);
    bool inbasin(Matter& Rmin0, Matter& Rcur);
    VectorXd getdimerforces();
    AtomMatrix rotateFridge(VectorXd Fridge,VectorXd Ncur,VectorXd Nridge);


private:
    Matter *matter; // Pointer to atom object \b outside the scope of the class.    
    Parameters *parameters; // Pointer to a structure outside the scope
                            // of the class containing runtime parameters. 
    double dt;
};

class RidgeDimerAtoms {
public:
    RidgeDimerAtoms(Matter *matt, Parameters *params);
    ~RidgeDimerAtoms();
    
    void initialize();

};
#endif
