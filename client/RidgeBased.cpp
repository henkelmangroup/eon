#include "RidgeBased.h"
#include "HelperFunctions.h"
#include "Matter.h"
#include "ImprovedDimer.h"
#include "Eigen.h"
#include <math.h>
#include <vector>

RidgeBased::RidgeBased(Matter *matt, Parameters *params)
  : Hyperdynamics(matt, params),
    R0(*matter), 
    R0old(*matter), 
    Rmin(*matter),
    Rdimer(*matter),
    dimer(&Rdimer, parameters)
{
  // Where I place parameters in the class variables
  matter = matt;    
  parameters = params;
  bool contains_hyperF = false;
  bool reachridge = false;
  bool checknow = false;
  double dt = 1/(parameters->timeUnit); //Timestep in Angstrom*sqrt(u/eV) where u is the atomic mass unit
  double T = parameters->temperature * parameters->kB;
  double dR = parameters->finiteDifference;
  double phi_tol = parameters->dimerConvergedAngle;
  double rotationsMax = parameters->dimerRotationsMax;
  double EBoostMax = parameters->ridgeBasedBMAX;
  double EBoost = EBoostMax;
  double EBoostOld = EBoostMax;
  double maxStep = parameters->ridgeBasedMAXSTEP;
  double curvature;
  double biasPot = 0.0;
  double bisect_tol = 0.02;
  double Eridge;
  int steps = 0;
  int checksteps = parameters->ridgeBasedCHECKSTEPS;
  int maxForceCalls = parameters->ridgeBasedMAXFCALLS;
  int forceCalls = 0;
  int objectCalls = 0; 
  int opt_forceCalls = 0;
  int dimerCalls = 0;
  AtomMatrix mode = Eigen::MatrixXd::Random(matter->getPositions().size()/3, 3);
  AtomMatrix Nguess = mode;
  AtomMatrix hyperF = R0.getForces();
  VectorXd N, Ftrans, Nridge, Fridge, V, F0;
  return;
}

RidgeBased::~RidgeBased(){ //destructor 
 // clean/free up memory/deallocate
 return;
}

void RidgeBased::initialize(){
  printf("Ridge-Based Hyperdynamics used\n");
}

double RidgeBased::boost(){
  EBoost = RidgeBased::get_biasPot(); //in Penghao's code this function is called get_biasPot
  return EBoost;
}

double RidgeBased::get_biasPot(){
  // Penghao's main code
  AtomMatrix mode;
  mode = matter->getPositions();
  mode = mode.setZero();
  mode(0, 0) = 1;

  // Run Dimer
  AtomMatrix r0_pos = R0.getPositions();
  AtomMatrix r0_old_pos = R0old.getPositions();
  if ((r0_pos - r0_old_pos).norm() == 0 && contains_hyperF == true) { 
    // Remember to chainge contains_hyperF to true whenever we assign hyper_F.
    return EBoostOld;
  }
  
  R0old = R0;
  steps = 0;
  dimerCalls = 0;
  int forceCalls = 0;
  opt_forceCalls = 0;
  checknow = false;
  V *= 0.0;
  objectCalls += 1;
  Nguess = Nguess / Nguess.norm();
 
  if (EBoost == EBoostMax && EBoostOld == EBoostMax){
    if (objectCalls %5 != 1){ //checking every 5 steps if we are in max Boost region.
      //Eigen::Matrix3f st = Eigen::Matrix3f::Zero(); //I think this matrix is in penghao code to accomodate for SSDimer. Uncomment if we need it for out dimer.
      hyperF = R0.getForces();
      contains_hyperF = true;
      return EBoostOld;
    }
  }
  else {
    double maxStep_copy = maxStep / 4;
    double phi_tol_copy = phi_tol / 3;
  }

  ImprovedDimer dimer = ImprovedDimer(&Rdimer, parameters);
  reachridge = false;
  bool quite = true;
  EBoostOld = EBoost;
  double minForce = 0.01;
  int interval = 1;
  RidgeBased::search(minForce, quite, maxForceCalls, interval);  
}

double RidgeBased::search(double minForce = 0.01, bool quite = false, int maxForceCalls = 500, int interval = 1){
  bool converged = false;
  double dE, dE_dimer;
  std::vector<Matter> traj_dimer;
  while (forceCalls < maxForceCalls){
    Matter pt = Rdimer;
    traj_dimer.push_back(pt);
    if (traj_dimer.size() > 1000){
      traj_dimer.erase(traj_dimer.begin());
      }
    dE = dimer.x0->getPotentialEnergy() - matter->getPotentialEnergy();
    try {
      dE_dimer = traj_dimer[-1].getPotentialEnergy() -traj_dimer[-2].getPotentialEnergy();
    }
    catch (double dE_dimer){
      dE_dimer = 100;
    }
    if ((dE >= EBoostMax && steps >= 1) || (dE_dimer < 0.01 && steps >= 3)) {
      printf("dE_dimer: %f", dE_dimer);
      double potdiff;
      for (Matter ptmp : traj_dimer){
        potdiff = ptmp.getPotentialEnergy() - R0.getPotentialEnergy();
        printf("%f",potdiff);
      }
      bool inbool = RidgeBased::inbasin(Rmin, *dimer.x0);
      if (inbool){
        if (dE >= EBoostMax){
          // exceeds the max energy
          hyperF = dimer.x0->getForces();
          return EBoostMax;
        }
        else {
          RidgeBased::bisect(traj_dimer, bisect_tol);
          // exceeds the max energy
          dE = Eridge - R0.getPotentialEnergy();
          if (dE  >= EBoostMax){
            hyperF = dimer.x0->getForces();
            return EBoostMax;
          }
          else {
            // reached the ridge
            hyperF = RidgeBased::rotateFridge(Fridge, N, Nridge);
            return std::max(dE, 0.0);
          }
        }
      }
    }
    else if (steps+2 == checksteps) {
      printf("reach self.checksteps:");
      for (Matter ptmp : traj_dimer){
        double potdiff = ptmp.getPotentialEnergy() - R0.getPotentialEnergy();
        printf("%f",potdiff);
      }
      bool inbool = RidgeBased::inbasin(Rmin, *dimer.x0);
      if (inbool == false){
        // reached the ridge
        RidgeBased::bisect(traj_dimer, bisect_tol);
        hyperF = rotateFridge(Fridge, N, Nridge);
        dE = Eridge - R0.getPotentialEnergy();
        return std::max(dE, 0.0);
      }
    }
    RidgeBased::step(); //Take a dimer step
    if (dimer.getEigenvalue() >= 0.0){
      hyperF = F0.reshaped(F0.size()/3,3);
      return EBoostMax;
    }
    if (quite == false){
      int ii = steps;
      int nf;
      double ff, cc, ee;
      ff = dimer.x0->maxForce();
      cc = dimer.getEigenvalue();
      ee = dimer.x0->getPotentialEnergy();
      nf = dimer.x0->getForceCalls() + opt_forceCalls;
      if (steps % 100 == 0 || steps == 1){
        printf("Iteration       Force       Curvature        Energy     ForceCalls");
        printf("-------------------------------------------------------------------------------");
        printf("%3i %13.6f %13.6f %13.6f %3i",ii,ff,cc,ee,nf);          
        }
      else{
        printf("%3i %13.6f %13.6f %13.6f %3i",ii,ff,cc,ee,nf);
        }
      }
    forceCalls = dimer.x0->getForceCalls() + opt_forceCalls;
  }
  if (dE < EBoostMax){
    // simple estimate
    // more accurate hyperforces, but has little affect on the rate
    hyperF = RidgeBased::rotateFridge(dimer.x0->getForces().reshaped<Eigen::RowMajor>(), N, dimer.tau);
    return std::max(dE, 0.0);
  }
  else {
    hyperF = dimer.x0->getForces();
    return EBoostMax;
  }
}

bool RidgeBased::inbasin(Matter& Rmin0, Matter& Rcur){
  // Minimize Rcur to see if minimize to Rmin0.
  // We start with minimizing using quickmin, then FIRE (this is how Penghao did it).
  Matter Rmin1(Rcur);
  parameters->optMethod = "qm";
  parameters->optConvergedForce = 0.001; //opt_ediffg
  parameters->optTimeStepInput = 0.2;
  parameters->optMaxIterations = 30; 
  bool first_check = Rmin1.relax();
  opt_forceCalls += Rmin1.getRelaxIterations();
  if (first_check == true){return first_check;}
  parameters->optMethod = "fire";
  parameters->optConvergedForce = 0.05; //opt_ediffg
  parameters->optTimeStepInput = 0.1; 
  parameters->optMaxMove = 0.1; 
  parameters->optMaxIterations = 200;
  bool second_check = Rmin1.relax();
  opt_forceCalls += Rmin1.getRelaxIterations(); 
  return second_check;
}

void RidgeBased::step(){
  Nguess = Nguess / Nguess.norm();
  steps += 1;
  double ENew = dimer.x0->getPotentialEnergy();
  //updating the forces parallel to the dimer.
  if (steps >= checksteps-1){
    Ftrans = RidgeBased::getdimerforces(); 
  }
  VectorXd temp_Ftrans = Ftrans;
  VectorXd dV = temp_Ftrans * dt; //Shouldn't this be dr??? Question for Penghao.
  VectorXd step;
  double inner_product; 
  if (((V.array(), Ftrans.array()).sum() > 0) || (steps == 0)){
    V = V / V.norm();
    step = maxStep * V;
    checknow = false;
  }
  else{
    V = dV;
    step = V * dt;
  }
  dimer.x0->setPositions(step.reshaped(V.size()/3, 3));
}

// This function returns the forces parallel to the dimer.
VectorXd RidgeBased::getdimerforces(){
  dimerCalls += 1;
  F0 = dimer.x0->getForces().reshaped<RowMajor>(); //dimer forces
  N = dimer.tau; // dimer mode
  return N;
}

// Finds the true point of the ridge along the dimer trajectory by linear inrterpolation between nearest two dimer points.
void RidgeBased::bisect(std::vector<Matter> traj_dimer, double tol = 0.02){
  bool inbool;
  Matter pout(*matter), pin(*matter), pt(*matter), pmid(*matter); 
  for (int i; i <= traj_dimer.size()-1; i++){
    Matter pout = traj_dimer[-1-i];
    Matter pt = traj_dimer[-2-i];
    inbool = RidgeBased::inbasin(Rmin, pt);
    if (inbool == true){
      Matter pin = pt;
      break;
    }
  }
  if ((traj_dimer.size() == 1) || (inbool == false)){
    printf("starting point of dimer is outside the basin");
    Fridge = F0.reshaped(F0.size()/3, 3);
    Eridge = R0.getPotentialEnergy();
    Nridge = dimer.tau;
    return;
  }
  AtomMatrix r0 = pin.getPositions();
  AtomMatrix r1 = pout.getPositions();
  AtomMatrix mode = r1 - r0;
  dR = tol + 1.0; 
  while (dR > tol){
    AtomMatrix rmid = 0.5*(r0 + r1);
    Matter pmid = pin;
    pmid.setPositions(rmid);
    inbool = RidgeBased::inbasin(Rmin, pmid);
    if (inbool = true){
      pin = pmid;
    }
    else{
      pout = pmid;
    }
    r0 = pin.getPositions();
    r1 = pout.getPositions();
    mode = r1 - r0;
    dR = mode.norm();
  }
  //std::vector<std::vector<double>> st(3, std::vector<double>(3, 0));
  Fridge = pmid.getForces().reshaped<Eigen::RowMajor>();
  Eridge = pmid.getPotentialEnergy();
  Nridge = mode.reshaped<Eigen::RowMajor>();
  reachridge = true;
}

AtomMatrix RidgeBased::rotateFridge(VectorXd Fridge, VectorXd Ncur, VectorXd Nridge){
  VectorXd N2 = Nridge;
  VectorXd N1 = Ncur;
  VectorXd N2perp = Ncur;
  VectorXd N1perp = Ncur;
  double N1N2 = N1.dot(N2);
  AtomMatrix hyper_F;
  if (1-abs(N1N2) < 1){
    hyper_F = Fridge.reshaped(Fridge.size()/3, 3);
  }
  else{
    N2perp = (N1 - N1N2 * N2) / (N1 - N1N2 * N2).norm();
    N1perp = (N2 - N1N2 * N1) / (N2 - N1N2 * N1).norm();
    double Fmag = (Fridge.array(), N2perp.array()).sum();
    hyper_F = (Fridge - Fmag * N2perp + Fmag * N1perp).reshaped((Fridge - Fmag * N2perp + Fmag * N1perp).size()/3,3);
  }
  return hyper_F;
}
