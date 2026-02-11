#include "RidgeBased.h"
#include "HelperFunctions.h"
#include "Matter.h"
#include "ImprovedDimer.h"
#include "Eigen.h"
#include <math.h>
#include <vector>
#include "Log.h"
#include <memory>
#include "SaddleSearchJob.h"
#include "MinModeSaddleSearch.h"

static const char LOG_PREFIX[] = "[RidgeBased]";

RidgeBased::RidgeBased(Matter *matt, Parameters *params)
  : Hyperdynamics(matt, params),
    R0(params), 
    R0old(params), 
    Rmin(params),
    Rdimer(params)
{
  // Where I place parameters in the class variables
  matter = matt;    
  parameters = params;
  R0 = *matt;     
  R0old = R0;
  Rmin = *matt;
  Rdimer = *matt;
  contains_hyperF = false;
  reachridge = false;
  checknow = false;
  dt = parameters->mdTimeStep; //Timestep in Angstrom*sqrt(u/eV) where u is the atomic mass unit
  T = parameters->temperature * parameters->kB;
  dR = parameters->finiteDifference;
  phi_tol = parameters->dimerConvergedAngle;
  rotationsMax = parameters->dimerRotationsMax;
  EBoostMax = parameters->ridgeBasedBMAX;
  EBoost = EBoostMax;
  log("%s EBoost value at initialization: %3.3f\n", LOG_PREFIX, EBoost);
  EBoostOld = EBoostMax;
  maxStep = parameters->ridgeBasedMAXSTEP;
  biasPot = 0.0;
  bisect_tol = 0.02;
  Eridge= 0.0;
  steps = 0;
  checksteps = parameters->ridgeBasedCHECKSTEPS;
  maxForceCalls = parameters->ridgeBasedMAXFCALLS;
  forceCalls = 0;
  objectCalls = 0; 
  opt_forceCalls = 0;
  dimerCalls = 0;
  inbasin_checks = 0;
  mode = Eigen::MatrixXd::Random(matter->getPositions().size()/3, 3);
  Nguess = mode;
  AtomMatrix hyperF = R0.getForces();
  Ftrans = R0.getForces();
  dim = std::make_unique<MinModeSaddleSearch>(&Rdimer, mode, R0.getPotentialEnergy(), parameters);
  int n_dof = matter->getPositions().size();
  V = Eigen::MatrixXd::Random(matter->getPositions().size()/3, 3);
  N = VectorXd::Zero(n_dof);
  F0 = VectorXd::Zero(n_dof);
  Fridge = VectorXd::Zero(n_dof);
  Nridge = VectorXd::Zero(n_dof);
  return;
}

RidgeBased::~RidgeBased(){ //destructor 
 // clean/free up memory/deallocate
 return;
}

void RidgeBased::initialize(){
}

double RidgeBased::boost(){
  R0.setPositions(matter->getPositions()); // gets new positions so that we can compare to old. (R0old)
  R0.getPotentialEnergy();
  EBoost = RidgeBased::get_biasPot(); //in Penghao's code this function is called get_biasPot
  log("%s EBoost value:		%3.7f\n", LOG_PREFIX, EBoost);
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
  if ((r0_pos - r0_old_pos).norm() < 1e-12 && contains_hyperF == true) { 
    // Remember to chainge contains_hyperF to true whenever we assign hyper_F.
    log("%s returning old right before boost() function\n", LOG_PREFIX); //I dont want this to return 
    return EBoostOld;
  }
  R0old = R0;
  steps = 0;
  dimerCalls = 0;
  opt_forceCalls = 0;
  forceCalls = 0;
  checknow = false;
  V *= 0.0;
  objectCalls += 1;
  inbasin_checks = 0;
  Nguess = Nguess / Nguess.norm();
 
  
  Rdimer = R0;
  double minForce = 0.01;
  int interval = 1;
  
  if (EBoost == EBoostMax && EBoostOld == EBoostMax){
    if (objectCalls %5 != 0){ //checking every 5 steps if we are in max Boost region.
      //Eigen::Matrix3f st = Eigen::Matrix3f::Zero(); //I think this matrix is in penghao code to accomodate for SSDimer. Uncomment if we need it for out dimer.
      hyperF = R0.getForces();
      contains_hyperF = true;
      return EBoostOld;
    }
    else{ // make dimer settings more liberal
      parameters->finiteDifference = 0.4;
      parameters->dimerConvergedAngle = 3.0;
      //return EBoostOld; //Not yet I think we need to allow for dimer to search
    }
  }
  //else {
  //  parameters->finiteDifference = parameters->finiteDifference / 4;
  //  parameters->dimerConvergedAngle = parameters->dimerConvergedAngle / 3;
  //}
  //dim->matter = &Rdimer; //Use the global dim so that both this function and search can see it and the modifications done to it.
  EBoostOld = EBoost;
  reachridge = false;
  bool quiet = false;
  //dim = std::make_unique<MinModeSaddleSearch>(&Rdimer, mode, R0.getPotentialEnergy(), parameters);
  EBoost = RidgeBased::search(minForce, quiet, maxForceCalls, interval); //Should be out of max boost region
  log("%s EBoost value from ridge-dimer search: %10.4f\n", LOG_PREFIX,EBoost);
  log("%s number of times inbasin was called for this EBoost value: %d\n", LOG_PREFIX, inbasin_checks);
  traj_dimer.clear(); //we don't want data to transfer over to next dimer search. 
  EBoostOld = EBoost;
  return EBoost;
}

double RidgeBased::search(double minForce = 0.01, bool quiet = false, int maxForceCalls = 500, int interval = 1){
  double dE, dE_dimer;
  while (forceCalls < maxForceCalls){
    Matter pt = Rdimer; //might segfault . The MinModeSaddleSearch::run modifies the Rdimer directly. At least basinhoppingsaddlesearch implies that.
    if (traj_dimer.size() > 1000){ // restaring the traj object 
      traj_dimer.erase(traj_dimer.begin());
      }
    dE = Rdimer.getPotentialEnergy() - R0.getPotentialEnergy(); //Comparing dimer image energy with energy of current MD step.
                                                                    //This dimer energy should be different because the dimer
                                                                    //Should be approaching the saddle.
    log("%s dE in forceCalls while loop:  %f  at forceCalls: %d\n", LOG_PREFIX, dE, forceCalls);
    log("%s Rdimer energy: %f\n",LOG_PREFIX, Rdimer.getPotentialEnergy());
    //log("%s traj_dimer last image energy: %f\n", LOG_PREFIX, traj_dimer.back()->getPotentialEnergy());
    if (traj_dimer.size() >= 2){
	    dE_dimer = traj_dimer.back()->getPotentialEnergy() - traj_dimer[traj_dimer.size()-2]->getPotentialEnergy(); //Comparing previous dimer image with current if it already ran.
    }
    else {
	    dE_dimer = 100.0;
    }

    if ((dE >= EBoostMax && steps >= 1) || (dE_dimer < 0.01 && steps >= 3)) {
      double potdiff;
      for (const std::shared_ptr<Matter>&  ptmp : traj_dimer){
        potdiff = abs(ptmp->getPotentialEnergy() - R0.getPotentialEnergy());
        log("%s Energy difference between Saddle Search and current MD step: %3.3f.\n",LOG_PREFIX, potdiff);
      }
      bool inbool = RidgeBased::inbasin(Rmin, Rdimer); //Minimize the climber to see if goes back to Rmin.
      inbasin_checks += 1;
      log("%sIn bool value: %d\n", LOG_PREFIX, inbool);
      if (inbool == true){
        if (dE >= EBoostMax){
          // Exceeds the max boost region energy so boost value is still max boost value.
          hyperF = Rdimer.getForces();
          log("%s In basin; Basin check did not minimize to antoher state.; In max boost region. \n", LOG_PREFIX);
          return EBoostMax;
        }
      }
      else {
        RidgeBased::bisect(traj_dimer, bisect_tol); // Dimer escaped ridge.
        dE = Eridge - R0.getPotentialEnergy(); // bisect updates Eridge and we comapare this value to R0.
        if (dE  >= EBoostMax){
          hyperF = R0.getForces();
          log("%s Climber found ridge region;System still in max boost region; Returning normal forces\n", LOG_PREFIX);
          return EBoostMax;
        }
        else {
          // reached the ridge
          hyperF = RidgeBased::rotateFridge(Fridge, N, Nridge);
          log("%s Climber found ridge region;System out of max boost region; Returning rotating ridge forces and returning hyper forces\n", LOG_PREFIX);
          return std::max(abs(dE), 0.0);
       }
      }
    }
    else if (steps+2 == checksteps) {
      for (const std::shared_ptr<Matter>& ptmp : traj_dimer){
        double potdiff = ptmp->getPotentialEnergy() - R0.getPotentialEnergy();
        printf("%f",potdiff);
      }
      bool inbool = RidgeBased::inbasin(Rmin, Rdimer);
      inbasin_checks += 1;
      if (inbool == false){
        // reached the ridge
        RidgeBased::bisect(traj_dimer, bisect_tol);
        hyperF = rotateFridge(Fridge, N, Nridge);
        dE = Eridge - R0.getPotentialEnergy();
        log("%s Climber found ridge region;Check because check+2 reached checksteps\n", LOG_PREFIX);
        return std::max(abs(dE), 0.0);
      }
    }
    else{
      log("%s skipped if statements. Here is traj_dimer size: %zu \n", LOG_PREFIX, traj_dimer.size());
    }
    dim->iteration = 0;

    RidgeBased::step(); //Take a dimer step
    traj_dimer.push_back(std::make_shared<Matter>(Rdimer));
    Rdimer.matter2con("Rdimer.con", true);
    log("%s Just took a RidgeBased step. Here is traj_dimer size: %zu \n", LOG_PREFIX, traj_dimer.size());
    if (dim->getEigenvalue() >= 0.0){
      hyperF = F0.reshaped(F0.size()/3,3);
      log("%s Positive curvature; Returning max boost value\n", LOG_PREFIX);
      return EBoostMax;
    }
    if (quiet == false){
      int ii = steps;
      int nf;
      double ff, cc, ee;
      ff = Rdimer.maxForce();
      cc = dim->getEigenvalue();
      ee = Rdimer.getPotentialEnergy();
      nf = Rdimer.getForceCalls() + opt_forceCalls;
      if (steps % 100 == 0 || steps == 1){
        log("%s Iteration       Force       Curvature        Energy     ForceCalls\n", LOG_PREFIX);
        log("-------------------------------------------------------------------------------\n");
        log("%s %3i %13.6f %13.6f %13.6f %3i\n",LOG_PREFIX,ii,ff,cc,ee,nf);          
        }
      else{
        log("%s %3i %13.6f %13.6f %13.6f %3i\n",LOG_PREFIX,ii,ff,cc,ee,nf);
        }
      }
    forceCalls = Rdimer.getForceCalls() + opt_forceCalls;
  }
  log("%s Rdimer.getForceCalls():   %d\n", LOG_PREFIX, Rdimer.getForceCalls());
  Rdimer.resetForceCalls();// zero out these force calls 

  if (dE < EBoostMax){
    // simple estimate
    // more accurate hyperforces, but has little affect on the rate
    log("%s dE (%f) should be lower than EBoostMax (%f)\n", LOG_PREFIX, dE, EBoostMax);
    hyperF = RidgeBased::rotateFridge(Rdimer.getForces().reshaped<Eigen::RowMajor>(), N, dim->getEigenvector().reshaped());
    log("%s abs(dE) before return statement: %3.3f\n", LOG_PREFIX, abs(dE));
    return std::max(abs(dE), 0.0);
  }
  else {
    hyperF = Rdimer.getForces();
    log("%s Final return statement if nothing else returns\n.", LOG_PREFIX);
    return EBoostMax;
  }
}

bool RidgeBased::inbasin(Matter& Rmin0,const Matter& Rcur){
  // Minimize Rcur to see if minimize to Rmin0.
  // We start with minimizing using quickmin, then FIRE (this is how Penghao did it).
  Matter Rmin1(parameters);
  Rmin1 = Rcur;
  //To check geometries at these high energy structures dimer is finding.
  parameters->optMethod = "qm";
  parameters->optConvergedForce = 0.001; //opt_ediffg
  parameters->optTimeStepInput = 0.2;
  parameters->optMaxIterations = 30; 
  Rmin1.relax();
  opt_forceCalls += Rmin1.getRelaxIterations(); 
  parameters->optMethod = "fire";
  parameters->optConvergedForce = 0.005; //opt_ediffg
  parameters->optTimeStepInput = 0.1; 
  parameters->optMaxMove = 0.1; 
  parameters->optMaxIterations = 1000;
  Rmin1.relax();
  opt_forceCalls += Rmin1.getRelaxIterations(); 
  bool same_basin = RidgeBased::identical(&Rmin1, &Rmin0, 1.0);
  if (same_basin ==false){
    const_cast<Matter&>(Rcur).matter2con("diff_basin.con", true);
  }
  log("%s same_basin ?? %b \n",LOG_PREFIX, same_basin);
  return same_basin;
}

void RidgeBased::step(){
  Nguess = Nguess / Nguess.norm();
  steps += 1;
  parameters->saddleMaxIterations=1;
  //updating the forces parallel to the dimer.
  dimerCalls += 1;
  //log("%s energy of Rdimer at this point should not be Rmin energy: %f\n", LOG_PREFIX, ENew);
  //log("%s energy of Rdimer should be R0, which here is : %f\n", LOG_PREFIX, R0.getPotentialEnergy());
  //log("%s The dimer step size and rotation tolerance, respectively: %f, %f \n",LOG_PREFIX, parameters->finiteDifference, parameters->dimerConvergedAngle);
  dim->run();//The position and forces will be updated here automatically.
  if (dimerCalls == 1){
    F0 = Rdimer.getForces().reshaped<RowMajor>(); //dimer forces
    N = dim->getEigenvector().reshaped(); // dimer mode
  }
    if (steps <= checksteps-1){
    //Ftrans = RidgeBased::getdimerforces();  //instead of getdimerforces, since we have MinModeSaddleSearch in eon, we can just take a step with that.
    log("%s making the dimer min mode search, dimerCalls: %d\n", LOG_PREFIX, dimerCalls);
  }
}

// Finds the true point of the ridge along the dimer trajectory by linear inrterpolation between nearest two dimer points.
void RidgeBased::bisect(const std::vector<std::shared_ptr<Matter>>& traj_dimer, double tol = 0.02){
  log("%s In bisect function\n", LOG_PREFIX);
  bool inbool;
  //const Matter& pout = *traj_dimer[traj_dimer.size()-1];
  //const Matter& pt = *traj_dimer[traj_dimer.size()-1] ;
  Matter pout(*matter), pin(*matter), pt(*matter), pmid(*matter); 
  log("%s traj_dimer.size(): %d\n", LOG_PREFIX, traj_dimer.size());
  if ((traj_dimer.size() == 1) || (inbool == false)){
    printf("starting point of dimer is outside the basin");
    Fridge = F0.reshaped(F0.size()/3, 3);
    Eridge = R0.getPotentialEnergy();
    Nridge = dim->getEigenvector().reshaped();
    return;
  }
  for (size_t i = traj_dimer.size()-1; i > 0; i--){
    log("%s i:  %d\n", LOG_PREFIX, i);
    pout = *traj_dimer[i];
    pt = *traj_dimer[i-1];
    log("%s Before in basin inside of bisect\n", LOG_PREFIX);
    inbool = RidgeBased::inbasin(Rmin, pt);
    inbasin_checks += 1;
    if (inbool == true){
      Matter pin = pt;
      break;
    }
  }
  AtomMatrix r0 = pin.getPositions();
  AtomMatrix r1 = pout.getPositions();
  AtomMatrix mode = r1 - r0;
  dR = tol + 1.0; 
  while (dR > tol){
    AtomMatrix rmid = 0.5*(r0 + r1);
    Matter pmid = pin;
    pmid.setPositions(rmid);
    log("%s Before in basin inside of bisect but the one in while loop\n", LOG_PREFIX);
    inbool = RidgeBased::inbasin(Rmin, pmid);
    if (inbool == true){
      pin = pmid;
    }
    else{
      pout = pmid;
    }
    r0 = pin.getPositions();
    r1 = pout.getPositions();
    mode = r1 - r0;
    dR = mode.norm();
    log("%s dR:   %10.6f  tol:  %10.6f\n", LOG_PREFIX, dR, tol);

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
  log("%s N1N2:  %3.3f\n",LOG_PREFIX,N1N2);
  if (1-abs(N1N2) < 1){
    hyperF = Fridge.reshaped(Fridge.size()/3, 3);
  }
  else{
    N2perp = (N1 - N1N2 * N2) / (N1 - N1N2 * N2).norm();
    N1perp = (N2 - N1N2 * N1) / (N2 - N1N2 * N1).norm();
    double Fmag = (Fridge.array(), N2perp.array()).sum();
    hyperF = (Fridge - Fmag * N2perp + Fmag * N1perp).reshaped((Fridge - Fmag * N2perp + Fmag * N1perp).size()/3,3);
  }
  matter->setForces(hyperF);
  log("%s passed hyper forces to matter object\n", LOG_PREFIX);
  contains_hyperF = true;
  return hyperF;
}

std::vector<std::shared_ptr<Matter>> RidgeBased::saddleTrajectory(){
  return traj_dimer;
}

bool RidgeBased::identical(const Matter* m1, const Matter* m2, const double tolerance)
{
    if (m1->numberOfAtoms() != m2->numberOfAtoms()) return false;
    
    int N = m1->numberOfAtoms();
    AtomMatrix r1 = m1->getPositions();
    AtomMatrix r2 = m2->getPositions();
    
    std::vector<bool> m2_claimed(N, false);
    int total_matched = 0;

    for (int i = 0; i < N; i++) {
        bool found_match = false;
        int type1 = m1->getAtomicNr(i);

        for (int j = 0; j < N; j++) {
            // Skip atoms in m2 already matched to something in m1
            if (m2_claimed[j]) continue;

            if (type1 == m2->getAtomicNr(j)) {
                double dist = (m1->pbc(r1.row(i) - r2.row(j))).norm();
                if (dist < tolerance) {
                    m2_claimed[j] = true;
                    total_matched++;
                    found_match = true;
                    break; 
                }
            }
        }

        // XXX Early Abort: If atom i couldn't find a partner, they aren't identical
        if (!found_match) return false;
    }

    return total_matched == N;
}
