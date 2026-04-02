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
  checksteps = 5;
  maxForceCalls = parameters->ridgeBasedMAXFCALLS;
  forceCalls = 0;
  objectCalls = 0; 
  opt_forceCalls = 0;
  dimerCalls = 0;
  inbasin_checks = 0;
  mode = Eigen::MatrixXd::Random(matter->getPositions().size()/3, 3);
  Nguess = mode;
  AtomMatrix hyperF = R0.getForces();
  //dim = std::make_unique<MinModeSaddleSearch>(&Rdimer, mode, R0.getPotentialEnergy(), parameters);
  idim = std::make_unique<ImprovedDimer>(&Rdimer, parameters); 
  n_dof = matter->getPositions().size();
  V = VectorXd::Zero(n_dof);
  dV = VectorXd::Zero(n_dof);
  N = VectorXd::Zero(n_dof);
  F0 = VectorXd::Zero(n_dof);
  Fridge = VectorXd::Zero(n_dof);
  Nridge = VectorXd::Zero(n_dof);
  Ftrans = VectorXd::Zero(n_dof);
  curvature = 0.0;
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
 
  
  Rdimer.setPositions(R0.getPositions());
  double minForce = 0.01;
  int interval = 1;
  
  if (EBoost == EBoostMax && EBoostOld == EBoostMax){
    maxStep = 0.4;
    parameters->dimerConvergedAngle = 3.0;
    if (objectCalls % 5 != 0){ //checking every 5 steps if we are in max Boost region.
      //Eigen::Matrix3f st = Eigen::Matrix3f::Zero(); //I think this matrix is in penghao code to accomodate for SSDimer. Uncomment if we need it for out dimer.
      hyperF = R0.getForces();
      contains_hyperF = true;
      return EBoostOld;
    }
  }
  else{ // No longer in a positive curvature region. Make dimer more conservative.
    maxStep = 0.1;
    parameters->dimerConvergedAngle = 1.0;
      //return EBoostOld; //Not yet I think we need to allow for dimer to search
  }
  EBoostOld = EBoost;
  reachridge = false;
  bool quiet = false;
  //dim = std::make_unique<MinModeSaddleSearch>(&Rdimer, mode, R0.getPotentialEnergy(), parameters);
  EBoost = RidgeBased::search(minForce, quiet, maxForceCalls, interval); //Should be out of max boost region
  //log("%s EBoost value from ridge-dimer search: %10.4f\n", LOG_PREFIX,EBoost);
  //log("%s number of times inbasin was called for this EBoost value: %d\n", LOG_PREFIX, inbasin_checks);
  // Assuming traj_dimer is a pointer to a vector, which the error message suggests
  //for (const auto& frame_ptr : traj_dimer) {
  //  if (frame_ptr) { // Always good to ensure the pointer isn't null
  //      frame_ptr->matter2con("traj_dimer.con", true);
  //  }
  //}
  traj_dimer.clear(); //we don't want data to transfer over to next dimer search. 
  EBoostOld = EBoost;
  return EBoost;
}

double RidgeBased::search(double minForce = 0.01, bool quiet = false, int maxForceCalls = 500, int interval = 1){
  double dE, dE_dimer;
  while (forceCalls < maxForceCalls){
    Matter pt(parameters);
    pt = Rdimer;
    if (traj_dimer.size() > 1000){ // restaring the traj object 
      traj_dimer.erase(traj_dimer.begin());
      }
    //log("%s Rdimer potential energy %f\n",LOG_PREFIX, Rdimer.getPotentialEnergy());
    //log("%s R0 potential energy %f\n",LOG_PREFIX, R0.getPotentialEnergy());

    dE = Rdimer.getPotentialEnergy() - R0.getPotentialEnergy(); //Comparing dimer image energy with energy of current MD step.
                                                                    //This dimer energy should be different because the dimer
                                                                    //Should be approaching the saddle.
    //log("%s dE in forceCalls while loop:  %f  at forceCalls: %d\n", LOG_PREFIX, dE, forceCalls);
    //log("%s Rdimer energy: %f\n",LOG_PREFIX, Rdimer.getPotentialEnergy());
    //log("%s traj_dimer last image energy: %f\n", LOG_PREFIX, traj_dimer.back()->getPotentialEnergy());
    if (traj_dimer.size() >= 2){
	    dE_dimer = traj_dimer.back()->getPotentialEnergy() - traj_dimer[traj_dimer.size()-2]->getPotentialEnergy(); //Comparing previous dimer image with current if it already ran.
    }
    else {
	    dE_dimer = 100.0;
    }

    if ((dE >= EBoostMax && steps >= 1) || (dE_dimer < 0.01 && steps >= 3)) {
      //double potdiff;
      //for (const std::shared_ptr<Matter>&  ptmp : traj_dimer){
      //  potdiff = abs(ptmp->getPotentialEnergy() - R0.getPotentialEnergy());
        //log("%s Energy difference between Saddle Search and current MD step: %3.3f.\n",LOG_PREFIX, potdiff);
      //}
      bool inbool = RidgeBased::inbasin(Rmin, Rdimer); //Minimize the climber to see if goes back to Rmin.
      inbasin_checks += 1;
      //log("%sIn bool value: %d\n", LOG_PREFIX, inbool);
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
      //log("%s skipped if statements. Here is traj_dimer size: %zu \n", LOG_PREFIX, traj_dimer.size());
    }
    //dim->iteration = 0;

    RidgeBased::step(); //Take a dimer step
    std::shared_ptr<Matter> new_image = std::make_shared<Matter>(parameters);
    *new_image = Rdimer;
    traj_dimer.push_back(new_image);
    //log("%s Just took a RidgeBased step. Here is its curvature: %f \n", LOG_PREFIX, curvature);
    //log("%s Some of the dimer settings\n", LOG_PREFIX);
    //log("%s Converged angle tolerance: %f \n", LOG_PREFIX, parameters->dimerConvergedAngle);
    //log("%s Finite Difference: %f \n", LOG_PREFIX, parameters->finiteDifference);
    if (curvature >= 0.0){
      hyperF = F0.reshaped(F0.size()/3,3);
      log("%s Positive curvature; Returning max boost value\n", LOG_PREFIX);
      return EBoostMax;
    }
    if (quiet == false){
      int ii = steps;
      int nf;
      double ff, cc, ee;
      ff = Rdimer.maxForce();
      cc = idim->getEigenvalue();
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
    hyperF = RidgeBased::rotateFridge(Rdimer.getForces().reshaped<Eigen::RowMajor>(), N, idim->getEigenvector().reshaped());
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
  std::string old_optMethod = parameters->optMethod;
  double old_optConvergedForce = parameters->optConvergedForce;
  double old_optTimeStepInput = parameters->optTimeStepInput;
  int old_optMaxIterations = parameters->optMaxIterations;
  double old_optMaxMove = parameters->optMaxMove;
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
  // 3. RESTORE global parameters
  parameters->optMethod = old_optMethod;
  parameters->optConvergedForce = old_optConvergedForce;
  parameters->optTimeStepInput = old_optTimeStepInput;
  parameters->optMaxIterations = old_optMaxIterations;
  parameters->optMaxMove = old_optMaxMove;

  //bool same_basin = identical(&Rmin1, &Rmin0, 0.5); Function that is defined only here for testing.
  bool same_basin = helper_functions::identical(&Rmin1, &Rmin0, 0.5);
  //bool same_basin = Rmin1.compare(&Rmin0);
  if (same_basin == false){
    const_cast<Matter&>(Rcur).matter2con("diff_basin.con", true);
  }
  //log("%s same_basin ?? %s \n", LOG_PREFIX, same_basin ? "true" : "false");
  return same_basin;
}

void RidgeBased::step(){
  //log("%s In step function\n", LOG_PREFIX);
  Nguess = Nguess / Nguess.norm();
  steps += 1;
  //updating the forces parallel to the dimer.
  //dim->run();//The position and forces will be updated here automatically??
  //if (dimerCalls == 1){
  //  F0 = Rdimer.getForces().reshaped<RowMajor>(); //dimer forces
  //  N = idim->getEigenvector().reshaped(); // dimer mode
  //}
  //log("%s steps inside of step function %d\n", LOG_PREFIX, steps);
  //log("%s checksteps inside of step function %d\n", LOG_PREFIX, checksteps);
  //log("%s rb_checksteps inside of step function %d\n", LOG_PREFIX, parameters->ridgeBasedCHECKSTEPS);
  //if (steps <= checksteps-1){
    //log("%s Right before running getdimerforces\n", LOG_PREFIX);
  Ftrans = RidgeBased::getdimerforces();
    //log("%s making the dimer min mode search, dimerCalls: %d\n", LOG_PREFIX, dimerCalls);
  //}
  VectorXd temp_Ftrans = Ftrans;
  dV = temp_Ftrans * dt;
  VectorXd dim_step;
  double VdotFtrans = V.dot(Ftrans);
  if (VdotFtrans > 0 || steps == 1){
    V = temp_Ftrans / temp_Ftrans.norm();
    //dim_step = maxStep * V;
    dim_step = maxStep * V;
    checknow = false;
  }
  else{
    V = dV;
    //dim_step = V * parameters->optTimeStepInput;
    dim_step = V * parameters->optTimeStep;
  }
  log("%s optTimeStep : %f\n", LOG_PREFIX, parameters->optTimeStep);
  log("%s maxStep : %f\n", LOG_PREFIX, maxStep);
  Rdimer.setPositionsV(Rdimer.getPositionsV() + dim_step);
}

// Finds the true point of the ridge along the dimer trajectory by linear inrterpolation between nearest two dimer points.
void RidgeBased::bisect(const std::vector<std::shared_ptr<Matter>>& traj_dimer, double tol=0.02){
  log("%s In bisect function\n", LOG_PREFIX);
  bool inbool = false;
  Matter pout(*matter), pin(*matter), pt(*matter);//, pmid(*matter);
  Matter pmid(parameters);
  pmid = *matter;
  if (traj_dimer.empty()) {
    log("%s traj_dimer is empty!\n", LOG_PREFIX);
    return;
  } 
  // The following loop parses thorugh climber trajecotry to see which index is the last to be in the reactant basin.
  for (size_t i = traj_dimer.size()-1; i > 0; i--){
    pout = *traj_dimer[i];
    pt = *traj_dimer[i-1];
    inbool = RidgeBased::inbasin(Rmin, pt);
    inbasin_checks += 1;
    if (inbool == true){
      Matter pin = pt;
      break;
    }
  }

  if ((traj_dimer.size() == 1) || (inbool == false)){
    log("%s starting point of dimer is outside the basin\n", LOG_PREFIX);
    //Fridge = F0.reshaped(F0.size()/3, 3);
    Fridge = F0;
    Eridge = R0.getPotentialEnergy();
    Nridge = idim->getEigenvector().reshaped();
    return;
  }
  
  // At this point we know where the climber left the basin. Now we need to find the exact ridge point.
  AtomMatrix r0 = pin.getPositions();
  AtomMatrix r1 = pout.getPositions();
  /*// Doing this for 2-D Potential
  r0(2) = 0.00;
  r1(2) = 0.00;
  */ 
  //log("%s In bisect function coordinates:\n", LOG_PREFIX);
  //for (int i = 0; i < r1.size(); ++i) {
  //    log("%s r0 [%d]: %f\n", LOG_PREFIX, i, r0(i));
  //    log("%s r1 [%d]: %f\n", LOG_PREFIX, i, r1(i));
  //}
  AtomMatrix mode = r1 - r0;
  dR = tol + 1.0; 
  while (dR > tol){
    AtomMatrix rmid = 0.5*(r0 + r1);
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
  /*
  log("%s Requesting forces from pmid...\n", LOG_PREFIX);
  Eigen::MatrixXd temp_forces = pmid.getForces();
  log("%s Forces retrieved! Reshaping...\n", LOG_PREFIX);
  Fridge = temp_forces.reshaped<Eigen::RowMajor>();
  log("%s Fridge assignment successful!\n", LOG_PREFIX);

  log("%s Before matter2con\n", LOG_PREFIX);
  const_cast<Matter&>(pmid).matter2con("ridge.con", true);

  log("%s After matter2con, getting forces\n", LOG_PREFIX);
  Fridge = pmid.getForces().reshaped<Eigen::RowMajor>();

  log("%s After forces, getting energy\n", LOG_PREFIX);
  Eridge = pmid.getPotentialEnergy();

  log("%s Finished block safely\n", LOG_PREFIX);
  */
  const_cast<Matter&>(pmid).matter2con("ridge.con", true);
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

// This fucntion is for testing 2-D potentials such as Voter97.
bool RidgeBased::identical(const Matter* m1, const Matter* m2, const double tolerance)
{
    if (m1->numberOfAtoms() != m2->numberOfAtoms()) return false;
    
    double distance;
    bool same;
    //int N = m1->numberOfAtoms();
    AtomMatrix r1 = m1->getPositions();
    AtomMatrix r2 = m2->getPositions();
    r2(2) = 0.00;
    r1(2) = 0.00;
    //for (int i = 0; i < r1.size(); ++i) {
    //  log("%s r1 [%d]: %f\n", LOG_PREFIX, i, r1(i));
    //  log("%s r2 [%d]: %f\n", LOG_PREFIX, i, r2(i));
    //}
    distance = (r1-r2).norm();
    if (distance < parameters->distanceDifference){
      same = true;
    }
    else{
      same = false;
      log("%s distance difference: %f same: %s\n", LOG_PREFIX, distance, same ? "true" : "false");
    }
  
    return same;
    /*std::vector<bool> m2_claimed(N, false);
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
    */
}

VectorXd RidgeBased::getdimerforces(){
  //log("%s In getdimerforces\n", LOG_PREFIX);
  dimerCalls += 1;
  if (dimerCalls == 1){
    idim->compute(&Rdimer, Nguess); //computing dimer forces. 
  }
  else{
    idim->compute(&Rdimer, AtomMatrix::Map(N.data(), n_dof, 3)); //computing dimer forces. 
  }
  F0 = Rdimer.getForcesV(); //get forces at the dimer step.
  N = Eigen::VectorXd::Map(idim->getEigenvector().data(), n_dof); // get mode computed from dimer.
  VectorXd Fpara = F0.dot(N) * N; // np.vdot(F0, N) * N;
  Ftrans = -Fpara;
  if (dimerCalls == 1){
     curvature = idim->getEigenvalue();
     log("%s curvature in getdimerforces: %f\n", LOG_PREFIX, curvature);
  }
  return Ftrans;
}
