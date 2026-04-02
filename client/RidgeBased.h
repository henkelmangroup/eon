#ifndef RIDGEBASED_H
#define RIDGEBASED_H

#include "HelperFunctions.h"
#include "Hyperdynamics.h"
#include "Matter.h"
#include "Parameters.h"
#include "ImprovedDimer.h"
#include <memory>
#include "Eigen.h"
#include "SaddleSearchJob.h"
#include "MinModeSaddleSearch.h"

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
    int inbasin_checks;
    int n_dof;
    Matter R0;
    Matter R0old;
    Matter Rmin;
    Matter Rdimer;
    std::unique_ptr<MinModeSaddleSearch> dim;
    std::unique_ptr<ImprovedDimer> idim;
    AtomMatrix hyperF, mode, Nguess;
    VectorXd V, dV, N, Nridge, Fridge, F0, Ftrans;
    std::vector<std::shared_ptr<Matter>> traj_dimer;
    

    void initialize() override;
    void bisect(const std::vector<std::shared_ptr<Matter>>& traj_dimer , double tol);
    void step();
    double boost() override;
    double get_biasPot();
    double search(double minForce, bool quite, int maxForceCalls, int interval);
    bool inbasin(Matter& Rmin0,const Matter& Rcur);
    VectorXd getdimerforces();
    AtomMatrix rotateFridge(VectorXd Fridge,VectorXd Ncur,VectorXd Nridge);
    std::vector<std::shared_ptr<Matter>> saddleTrajectory();
    bool identical(const Matter* m1, const Matter* m2, const double tolerance);


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
