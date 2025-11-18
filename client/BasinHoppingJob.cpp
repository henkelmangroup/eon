
#include <math.h>
#include <stdio.h>
#include <string>
#include <sstream>
#include <vector>
#include <algorithm>
#include <map>

#include "BasinHoppingJob.h"
#include "Dynamics.h"
#include "HelperFunctions.h"
#include "ObjectiveFunction.h"
#include "Optimizer.h"
#include "Potential.h"
#include "Log.h"

using namespace std;
using namespace helper_functions;

BasinHoppingJob::BasinHoppingJob(Parameters *params)
{
    parameters = params;
    current = new Matter(parameters);
    trial = new Matter(parameters);
    fcalls = Potential::fcalls;
    atomListsCached = false;
}

BasinHoppingJob::~BasinHoppingJob()
{
    delete current;
    delete trial;

    for (unsigned int i = 0; i < uniqueStructures.size(); i++)
    {
        delete uniqueStructures[i];
    }
}

std::vector<std::string> BasinHoppingJob::run(void)
{
    bool swapMove;
    double swap_accept = 0.0;
    jump_count = 0; // count of jump movies
    swap_count = 0; // count of swap moves
    disp_count = 0; // count of displacement moves
    int consecutive_rejected_trials = 0;
    double totalAccept = 0.0;
    Matter *minTrial = new Matter(parameters);
    Matter *swapTrial = new Matter(parameters);

    string conFilename = getRelevantFile(parameters->conFilename);
    current->con2matter(conFilename);

    // Sanity Check
    vector<long> Elements;
    Elements = getElements(current);
    if (parameters->basinHoppingSwapProbability > 0 && Elements.size() == 1)
    {
        char msg[] = "error: [Basin Hopping] swap move probability "
                     "must be zero if there is only one element type\n";
        log(msg);
        exit(1);
    }

    double randomProb = parameters->basinHoppingInitialRandomStructureProbability;
    if (randomProb > 0.0)
    {
        log("generating random structure with probability %.4f\n", randomProb);
    }
    double u = helper_functions::random();
    if (u < parameters->basinHoppingInitialRandomStructureProbability)
    {
        AtomMatrix randomPositions = current->getPositionsFree();
        for (int i = 0; i < current->numberOfFreeAtoms(); i++)
        {
            for (int j = 0; j < 3; j++)
            {
                randomPositions(i, j) = helper_functions::random();
            }
        }
        randomPositions *= current->getCell();
        current->setPositionsFree(randomPositions);

        pushApart(current, parameters->basinHoppingPushApartDistance);
    }

    *trial = *current;
    *minTrial = *current;

    current->relax(true);

    double currentEnergy = current->getPotentialEnergy();
    double minimumEnergy = currentEnergy;

    Matter *minimumEnergyStructure = new Matter(parameters);
    *minimumEnergyStructure = *current;
    int nsteps = parameters->basinHoppingSteps + parameters->basinHoppingQuenchingSteps;
    long totalfc;
    FILE *pFile;
    pFile = fopen("bh.dat", "w");

    log("[Basin Hopping] %4s %12s %12s %12s %4s %5s %5s\n",
        "step", "current", "trial", "global min", "fc", "ar", "md");
    log("[Basin Hopping] %4s %12s %12s %12s %4s %5s %5s\n",
        "----", "-------", "-----", "----------", "--", "--", "--");

    int recentAccept = 0;
    double curDisplacement = parameters->basinHoppingDisplacement;

    // Initialize atom selection lists if needed
    log("[Basin Hopping] DEBUG: Initializing atom selection\n");
    std::vector<int> initialSelection = selectAtomsForDisplacement();
    log("[Basin Hopping] DEBUG: Atom selection initialized: selectedAtoms.size()=%zu\n", initialSelection.size());
    if (!initialSelection.empty())
    {
        log("[Basin Hopping] DEBUG: Selected atoms: ");
        for (size_t i = 0; i < initialSelection.size(); i++)
        {
            log("%d ", initialSelection[i]);
        }
        log("\n");
    }

    for (int step = 0; step < nsteps; step++)
    {

        // Swap or displace
        if (randomDouble(1.0) < parameters->basinHoppingSwapProbability &&
            step < parameters->basinHoppingSteps)
        {
            *swapTrial = *current;
            randomSwap(swapTrial);
            swapMove = true;
            *minTrial = *swapTrial;
        }
        else
        {
            AtomMatrix displacement;
            std::vector<int> selectedAtoms = selectAtomsForDisplacement();

            log("[Basin Hopping] DEBUG: Step %d displacement: selectedAtoms.size()=%zu\n",
                step + 1, selectedAtoms.size());

            // Use atom selection if configured, otherwise use default
            if (!selectedAtoms.empty())
            {
                displacement = displaceRandomSelected(curDisplacement, selectedAtoms);
                log("[Basin Hopping] DEBUG: Using atom selection for displacement, epicenters (0-based): ");
                for (size_t i = 0; i < selectedAtoms.size(); i++)
                {
                    log("%d ", selectedAtoms[i]);
                }
                log("\n");
            }
            else
            {
                displacement = displaceRandom(curDisplacement);
                log("[Basin Hopping] DEBUG: Using default displacement (all atoms)\n");
            }

            trial->setPositions(current->getPositions() + displacement);
            swapMove = false;
            pushApart(trial, parameters->basinHoppingPushApartDistance);

            // Debug: Log coordinates of selected atoms before relax
            if (!selectedAtoms.empty())
            {
                log("[Basin Hopping] DEBUG: Coordinates of selected atoms BEFORE relax:\n");
                for (size_t i = 0; i < selectedAtoms.size(); i++)
                {
                    int idx = selectedAtoms[i];
                    if (idx >= 0 && idx < trial->numberOfAtoms())
                    {
                        log("[Basin Hopping] DEBUG: Atom #%d (idx0=%d): [%.6f, %.6f, %.6f]\n",
                            idx, idx, trial->getPosition(idx, 0), trial->getPosition(idx, 1), trial->getPosition(idx, 2));
                    }
                }
            }

            *minTrial = *trial;
        }

        if (parameters->writeMovies == true)
        {
            trial->matter2con("trials", true);
        }

        Potential::fcalls = 0;
        minTrial->relax(true);
        int minfcalls = Potential::fcalls;

        // Debug: Log coordinates of selected atoms after relax
        std::vector<int> selectedAtomsAfter = selectAtomsForDisplacement();
        if (!selectedAtomsAfter.empty())
        {
            log("[Basin Hopping] DEBUG: Coordinates of selected atoms AFTER relax:\n");
            for (size_t i = 0; i < selectedAtomsAfter.size(); i++)
            {
                int idx = selectedAtomsAfter[i];
                if (idx >= 0 && idx < minTrial->numberOfAtoms())
                {
                    log("[Basin Hopping] DEBUG: Atom #%d (idx0=%d): [%.6f, %.6f, %.6f]\n",
                        idx, idx, minTrial->getPosition(idx, 0), minTrial->getPosition(idx, 1), minTrial->getPosition(idx, 2));
                }
            }
        }

        double deltaE = minTrial->getPotentialEnergy() - currentEnergy;
        double p = 0.0;
        if (step >= parameters->basinHoppingSteps)
        {
            if (deltaE <= 0.0)
            {
                p = 1.0;
            }
        }
        else
        {
            if (deltaE <= 0.0)
            {
                p = 1.0;
            }
            else
            {
                p = exp(-deltaE / (parameters->temperature * 8.6173324e-5));
            }
        }

        bool accepted = false;
        if (randomDouble(1.0) < p)
        {
            accepted = true;
            if (parameters->basinHoppingSignificantStructure)
            {
                *current = *minTrial;
            }
            else
            {
                *current = *trial;
            }
            if (swapMove)
            {
                swap_accept += 1;
            }
            if (step < parameters->basinHoppingSteps)
            {
                totalAccept += 1;
                recentAccept += 1;
            }

            currentEnergy = minTrial->getPotentialEnergy();

            if (currentEnergy < minimumEnergy)
            {
                minimumEnergy = currentEnergy;
                *minimumEnergyStructure = *minTrial;
                minimumEnergyStructure->matter2con("min.con");
            }

            consecutive_rejected_trials = 0;

            if (parameters->basinHoppingWriteUnique)
            {
                bool newStructure = true;
                for (unsigned int i = 0; i < uniqueEnergies.size(); i++)
                {
                    // if minTrial has a different energy or a different structure
                    // it is new, otherwise it is old
                    if (fabs(currentEnergy - uniqueEnergies[i]) <
                        parameters->energyDifference)
                    {
                        if (current->compare(uniqueStructures[i],
                                             parameters->indistinguishableAtoms) == true)
                        {
                            newStructure = false;
                        }
                    }
                }

                if (newStructure)
                {
                    uniqueEnergies.push_back(currentEnergy);
                    Matter *currentCopy = new Matter(parameters);
                    *currentCopy = *current;
                    uniqueStructures.push_back(currentCopy);

                    char fname[128];
                    snprintf(fname, 128, "min_%.5i.con", step + 1);
                    current->matter2con(fname);
                    returnFiles.push_back(fname);

                    snprintf(fname, 128, "energy_%.5i.dat", step + 1);
                    returnFiles.push_back(fname);
                    FILE *fh = fopen(fname, "w");
                    fprintf(fh, "%.10e\n", currentEnergy);
                    fclose(fh);
                }
            }
        }
        else
        {
            consecutive_rejected_trials++;
        }

        if (parameters->writeMovies == true)
        {
            minTrial->matter2con("movie", true);
        }

        totalfc = Potential::fcallsTotal;
        char acceptReject[2];
        acceptReject[1] = '\0';
        if (accepted)
        {
            acceptReject[0] = 'A';
        }
        else
        {
            acceptReject[0] = 'R';
        }
        log("[Basin Hopping] %5i %12.3f %12.3f %12.3f %4i %5.3f %5.3f %1s\n",
            step + 1, currentEnergy, minTrial->getPotentialEnergy(), minimumEnergy,
            minfcalls, totalAccept / ((double)step + 1), curDisplacement, acceptReject);
        fprintf(pFile, "%6i %9ld %12.4e %12.4e\n", step + 1, totalfc, currentEnergy,
                minTrial->getPotentialEnergy());

        if (minimumEnergy < parameters->basinHoppingStopEnergy)
        {
            break;
        }

        if (consecutive_rejected_trials == parameters->basinHoppingJumpMax &&
            step < parameters->basinHoppingSteps)
        {
            consecutive_rejected_trials = 0;
            AtomMatrix jump;
            for (int j = 0; j < parameters->basinHoppingJumpSteps; j++)
            {
                jump_count++;
                jump = displaceRandom(curDisplacement);
                current->setPositions(current->getPositions() + jump);
                if (parameters->basinHoppingSignificantStructure)
                {
                    pushApart(current, parameters->basinHoppingPushApartDistance);
                    current->relax(true);
                }
                currentEnergy = current->getPotentialEnergy();
                if (currentEnergy < minimumEnergy)
                {
                    minimumEnergy = currentEnergy;
                    *minimumEnergyStructure = *current;
                }
            }
        }

        int nadjust = parameters->basinHoppingAdjustPeriod;
        double adjustFraction = parameters->basinHoppingAdjustFraction;
        if ((step + 1) % nadjust == 0 &&
            parameters->basinHoppingAdjustDisplacement == true)
        {
            double recentRatio = ((double)recentAccept) / ((double)nadjust);
            if (recentRatio > parameters->basinHoppingTargetRatio)
            {
                curDisplacement *= 1.0 + adjustFraction;
            }
            else
            {
                curDisplacement *= 1.0 - adjustFraction;
            }

            // log("recentRatio %.3f md: %.3f\n", recentRatio, curDisplacement);
            recentAccept = 0;
        }
    }
    fclose(pFile);

    /* Save Results */

    FILE *fileResults, *fileProduct;

    std::string resultsFilename("results.dat");
    returnFiles.push_back(resultsFilename);
    fileResults = fopen(resultsFilename.c_str(), "wb");

    if (parameters->writeMovies == true)
    {
        std::string movieFilename("movie.xyz");
        returnFiles.push_back(movieFilename);
    }

    fprintf(fileResults, "%d termination_reason\n", 0);
    fprintf(fileResults, "%.6f minimum_energy\n", minimumEnergy);
    fprintf(fileResults, "%ld random_seed\n", parameters->randomSeed);
    fprintf(fileResults, "%.3f acceptance_ratio\n",
            totalAccept / parameters->basinHoppingSteps);
    if (parameters->basinHoppingSwapProbability > 0)
    {
        fprintf(fileResults, "%.3f swap_acceptance_ratio\n",
                swap_accept / double(swap_count));
    }
    fprintf(fileResults, "%ld total_normal_displacement_steps\n",
            disp_count - jump_count - parameters->basinHoppingQuenchingSteps);
    fprintf(fileResults, "%d total_jump_steps\n", jump_count);
    fprintf(fileResults, "%d total_swap_steps\n", swap_count);
    fprintf(fileResults, "%d total_force_calls\n", Potential::fcallsTotal);
    fclose(fileResults);

    std::string productFilename("min.con");
    returnFiles.push_back(productFilename);
    fileProduct = fopen(productFilename.c_str(), "wb");
    minimumEnergyStructure->matter2con(fileProduct);
    fclose(fileProduct);

    std::string bhFilename("bh.dat");
    returnFiles.push_back(bhFilename);

    delete minTrial;
    delete minimumEnergyStructure;
    delete swapTrial;
    return returnFiles;
}

AtomMatrix BasinHoppingJob::displaceRandom(double curDisplacement)
{
    disp_count++;
    // Create a random displacement
    AtomMatrix displacement;
    displacement.resize(trial->numberOfAtoms(), 3);
    displacement.setZero();
    VectorXd distvec = calculateDistanceFromCenter(current);
    int num = trial->numberOfAtoms();
    int m = 0;
    if (parameters->basinHoppingSingleAtomDisplace)
    {
        m = randomInt(0, trial->numberOfAtoms() - 1);
        num = m + 1;
    }

    for (int i = m; i < num; i++)
    {
        double dist = distvec(i);
        double disp = 0.0; // displacement size, possibly scaled

        if (!trial->getFixed(i))
        {
            if (parameters->basinHoppingDisplacementAlgorithm == "standard")
            {
                disp = curDisplacement;
            }
            // scale displacement linearly with the particle radius
            else if (parameters->basinHoppingDisplacementAlgorithm == "linear")
            {
                double Cs = curDisplacement / distvec.maxCoeff();
                disp = Cs * dist;
            }
            // scale displacement quadratically with the particle radius
            else if (parameters->basinHoppingDisplacementAlgorithm == "quadratic")
            {
                double Cq = curDisplacement / (distvec.maxCoeff() * distvec.maxCoeff());
                disp = Cq * dist * dist;
            }
            else
            {
                log("Unknown displacement_algorithm\n");
                std::exit(1);
            }
            for (int j = 0; j < 3; j++)
            {
                if (parameters->basinHoppingDisplacementDistribution == "uniform")
                {
                    displacement(i, j) = randomDouble(2 * disp) - disp;
                }
                else if (parameters->basinHoppingDisplacementDistribution == "gaussian")
                {
                    displacement(i, j) = gaussRandom(0.0, disp);
                }
                else
                {
                    log("Unknown displacement_distribution\n");
                    std::exit(1);
                }
            }
        }
    }
    return displacement;
}

void BasinHoppingJob::randomSwap(Matter *matter)
{
    swap_count++;
    vector<long> Elements;
    Elements = getElements(matter);

    long ela;
    long elb;
    long ia = randomInt(0, Elements.size() - 1);
    ela = Elements.at(ia);
    Elements.erase(Elements.begin() + ia);

    long ib = randomInt(0, Elements.size() - 1);
    elb = Elements.at(ib);

    int changera = 0;
    int changerb = 0;

    changera = randomInt(0, matter->numberOfAtoms() - 1);
    while (matter->getAtomicNr(changera) != ela)
    {
        changera = randomInt(0, matter->numberOfAtoms() - 1);
    }

    changerb = randomInt(0, matter->numberOfAtoms() - 1);
    while (matter->getAtomicNr(changerb) != elb)
    {
        changerb = randomInt(0, matter->numberOfAtoms() - 1);
    }

    double posax = matter->getPosition(changera, 0);
    double posay = matter->getPosition(changera, 1);
    double posaz = matter->getPosition(changera, 2);

    matter->setPosition(changera, 0, matter->getPosition(changerb, 0));
    matter->setPosition(changera, 1, matter->getPosition(changerb, 1));
    matter->setPosition(changera, 2, matter->getPosition(changerb, 2));

    matter->setPosition(changerb, 0, posax);
    matter->setPosition(changerb, 1, posay);
    matter->setPosition(changerb, 2, posaz);
}

vector<long> BasinHoppingJob::getElements(Matter *matter)
{
    int allElements[118] = {0};
    vector<long> Elements;

    for (long y = 0; y < matter->numberOfAtoms(); y++)
    {
        if (!matter->getFixed(y))
        {
            int index = matter->getAtomicNr(y);
            allElements[index] = 1;
        }
    }

    for (int i = 0; i < 118; i++)
    {
        if (allElements[i] != 0)
        {
            Elements.push_back(i);
        }
    }

    return Elements;
}

VectorXd BasinHoppingJob::calculateDistanceFromCenter(Matter *matter)
{
    AtomMatrix pos = matter->getPositions();
    Vector3d cen(0, 0, 0);
    int num = matter->numberOfAtoms();

    cen = pos.colwise().sum() / (double)num;

    VectorXd dist(num);

    for (int n = 0; n < num; n++)
    {
        pos.row(n) -= cen;
        dist(n) = pos.row(n).norm();
    }

    return dist;
}

/**
 * Parse atom list string into vector of atom indices.
 *
 * IMPORTANT: Atom indexing convention
 * -----------------------------------
 * EON uses 0-based indexing internally (C/C++ convention):
 *   - First atom: index 0

// ============================================================================
// ATOM INDEXING CONVENTION
// ============================================================================
// EON uses 0-based indexing throughout (C/C++ array convention):
//   - Atom indices range from 0 to (nAtoms-1)
//   - First atom: index 0, Last atom: index (nAtoms-1)
//
// pos.con file displays 0-based numbering:
//   - First atom displayed as "0", Last atom displayed as "(nAtoms-1)"
//   - This matches the internal indexing
//
// Configuration (config.ini) uses 0-based indexing:
//   - displace_atom_list=0,1,2    selects atoms at indices 0, 1, 2 (first three atoms)
//   - displace_atom_list=0        selects first atom (displayed as "0" in pos.con)
//   - displace_atom_list=320      selects 321st atom (displayed as "320" in pos.con)
//
// displace_type_list uses atomic numbers (element-based, not index-based):
//   - displace_type_list=78       selects all Pt atoms
//   - displace_type_list=Pt       same as above (element symbol supported)
//   - displace_type_list=8,78     selects all O and Pt atoms
// ============================================================================

/**
 * Parse atom list string from config.ini.
 *
 * INDEXING: Uses 0-based indexing (C/C++ convention):
 *   - First atom: index 0
 *   - Last atom:  index (nAtoms - 1)
 *
 * pos.con file displays 0-based numbering:
 *   - First atom shown as: 0
 *   - Last atom shown as:  (nAtoms - 1)
 *
 * Configuration file (config.ini) uses 0-based indexing:
 *   - displace_atom_list=0,1,2  means atoms at indices 0, 1, 2 (first three atoms)
 *   - displace_atom_list=0      refers to first atom (shown as "0" in pos.con)
 *   - displace_atom_list=320    refers to 321st atom (shown as "320" in pos.con)
 *
 * Supports comma-separated indices and negative indices (-1 = last atom).
 * Example: "0,1,2" selects first 3 atoms, "10,20,-1" selects atoms at index 10, 20, and last.
 *
 * @param atomListStr Comma-separated string of atom indices
 * @param nAtoms Total number of atoms in the system
 * @return Vector of valid atom indices (sorted, duplicates removed)
 */
std::vector<int> BasinHoppingJob::parseAtomList(const std::string &atomListStr, int nAtoms)
{
    std::vector<int> atomIndices;
    if (atomListStr.empty())
    {
        return atomIndices;
    }

    std::stringstream ss(atomListStr);
    std::string item;

    while (std::getline(ss, item, ','))
    {
        // Trim whitespace (safe implementation: check npos before erase)
        size_t first = item.find_first_not_of(" \t");
        if (first != std::string::npos)
        {
            item.erase(0, first);
        }
        else
        {
            item.clear(); // All whitespace, clear the string
        }

        size_t last = item.find_last_not_of(" \t");
        if (last != std::string::npos)
        {
            item.erase(last + 1);
        }
        // If last == npos, string is empty or all whitespace, no need to erase

        if (item.empty())
            continue;

        int parsedIdx;
        try
        {
            size_t pos = 0;
            parsedIdx = std::stoi(item, &pos);
            if (pos != item.size())
                throw std::invalid_argument("trailing");
        }
        catch (...)
        {
            log("[Basin Hopping] WARNING: invalid atom index '%s', skipping\n", item.c_str());
            continue;
        }

        int idx = -1;

        if (parsedIdx < 0)
        {
            // Negative indices remain relative to the end: -1 -> last atom (nAtoms-1)
            idx = nAtoms + parsedIdx;
        }
        else
        {
            // 0-based indexing: use index as-is
            idx = parsedIdx;
        }

        // Validate index after conversion
        if (idx >= 0 && idx < nAtoms)
        {
            atomIndices.push_back(idx);
        }
        else
        {
            log("[Basin Hopping] WARNING: displace_atom_list entry '%s' (index %d) is out of range [0, %d), skipping\n",
                item.c_str(), idx, nAtoms);
        }
    }

    // Remove duplicates and sort
    std::sort(atomIndices.begin(), atomIndices.end());
    atomIndices.erase(std::unique(atomIndices.begin(), atomIndices.end()), atomIndices.end());

    return atomIndices;
}

/**
 * Get atom indices matching specified atomic numbers or element symbols.
 *
 * IMPORTANT: Returns 0-based indices
 * -----------------------------------
 * This function returns atom indices using 0-based indexing (C/C++ convention).
 * The returned indices can be directly used to access arrays in Matter class.
 *
 * Parses comma-separated atomic numbers/symbols and finds all matching free atoms.
 * Example: "29,28" or "Cu,Ni" or "Pt,8,Au" (mixed format supported)
 *
 * Only includes FREE atoms (atoms with getFixed(i) == false).
 * Fixed atoms are automatically excluded from selection.
 * Parses comma-separated atomic numbers and finds all matching free atoms.
 * Example: "29,28" finds all Cu (29) and Ni (28) atoms.
 *
 * @param typeListStr Comma-separated string of atomic numbers
 * @param matter Matter object containing atomic structure
 * @return Vector of atom indices matching the types (only free atoms)
 */
std::vector<int> BasinHoppingJob::getAtomsByType(const std::string &typeListStr, Matter *matter)
{
    std::vector<int> atomIndices;
    if (typeListStr.empty())
    {
        return atomIndices;
    }

    // Parse atomic numbers/symbols from comma-separated string
    std::vector<long> targetTypes;
    std::stringstream ss(typeListStr);
    std::string item;

    while (std::getline(ss, item, ','))
    {
        // Trim whitespace (safe implementation: check npos before erase)
        size_t first = item.find_first_not_of(" \t");
        if (first != std::string::npos)
        {
            item.erase(0, first);
        }
        else
        {
            item.clear(); // All whitespace, clear the string
        }

        size_t last = item.find_last_not_of(" \t");
        if (last != std::string::npos)
        {
            item.erase(last + 1);
        }
        // If last == npos, string is empty or all whitespace, no need to erase

        if (item.empty())
            continue;

        // Try to parse as number first
        try
        {
            long atomicNr = std::stol(item);
            targetTypes.push_back(atomicNr);
            log("[Basin Hopping] DEBUG: Parsed atomic number: %ld\n", atomicNr);
        }
        catch (const std::invalid_argument &)
        {
            // Not a number, try as element symbol - use Matter::symbol2atomicNumber()
            int atomicNr = Matter::symbol2atomicNumber(item.c_str());
            if (atomicNr >= 0)
            {
                targetTypes.push_back(atomicNr);
                log("[Basin Hopping] DEBUG: Parsed element symbol '%s' -> atomic number %d\n",
                    item.c_str(), atomicNr);
            }
            else
            {
                log("[Basin Hopping] WARNING: Unknown element symbol or invalid number: '%s', skipping\n",
                    item.c_str());
            }
        }
    }

    // Find all atoms matching the types
    for (int i = 0; i < matter->numberOfAtoms(); i++)
    { // Loop through all atoms using 0-based index
        long atomicNr = matter->getAtomicNr(i);
        if (std::find(targetTypes.begin(), targetTypes.end(), atomicNr) != targetTypes.end())
        {
            if (!matter->getFixed(i))
            {                             // Only include free atoms
                atomIndices.push_back(i); // Store 0-based index
            }
        }
    }

    return atomIndices;
}

/**
 * Select atoms for displacement based on configuration parameters.
 * Priority: displace_atom_list > displace_type_list > all atoms (default).
 * Caches the atom lists on first call for efficiency.
 *
 * @return Vector of atom indices to use as displacement epicenters
 */
std::vector<int> BasinHoppingJob::selectAtomsForDisplacement()
{
    // Cache atom lists on first call
    if (!atomListsCached)
    {
        int nAtoms = current->numberOfAtoms();
        log("[Basin Hopping] DEBUG: Caching atom lists (nAtoms=%d)\n", nAtoms);

        // Parse atom list if provided
        if (!parameters->basinHoppingDisplaceAtomList.empty())
        {
            log("[Basin Hopping] DEBUG: Parsing displace_atom_list: '%s'\n",
                parameters->basinHoppingDisplaceAtomList.c_str());
            cachedAtomList = parseAtomList(parameters->basinHoppingDisplaceAtomList, nAtoms);
            log("[Basin Hopping] DEBUG: Parsed %zu atoms from displace_atom_list\n", cachedAtomList.size());
        }

        // Parse type list if provided
        if (!parameters->basinHoppingDisplaceTypeList.empty())
        {
            log("[Basin Hopping] DEBUG: Parsing displace_type_list: '%s'\n",
                parameters->basinHoppingDisplaceTypeList.c_str());
            cachedTypeList = getAtomsByType(parameters->basinHoppingDisplaceTypeList, current);
            log("[Basin Hopping] DEBUG: Found %zu atoms matching displace_type_list\n", cachedTypeList.size());
        }

        atomListsCached = true;
    }

    std::vector<int> selectedAtoms;

    // Priority: atom_list > type_list > all atoms
    // When atom selection is specified, always use ALL selected atoms (not random selection)
    // This ensures only the specified atoms participate in perturbation
    if (!cachedAtomList.empty())
    {
        // Always use all listed atoms when displace_atom_list is specified
        selectedAtoms = cachedAtomList;
        log("[Basin Hopping] DEBUG: Using all %zu atoms from displace_atom_list: ", selectedAtoms.size());
        for (size_t i = 0; i < selectedAtoms.size(); i++)
        {
            log("%d ", selectedAtoms[i]);
        }
        log("\n");
    }
    else if (!cachedTypeList.empty())
    {
        // Always use all atoms of listed types when displace_type_list is specified
        selectedAtoms = cachedTypeList;
        log("[Basin Hopping] DEBUG: Using all %zu atoms from displace_type_list\n", selectedAtoms.size());
    }
    else
    {
        log("[Basin Hopping] DEBUG: No atom selection specified, will use default (all atoms)\n");
    }
    // If no selection specified, return empty (will use default behavior)

    return selectedAtoms;
}

/**
 * Apply random displacement to selected atoms only (no neighbors).
 * Only the atoms specified in selectedAtoms will be displaced.
 * Uses the same displacement algorithm and distribution as displaceRandom().
 *
 * @param maxDisplacement Maximum displacement magnitude
 * @param selectedAtoms Vector of atom indices to displace
 * @return Displacement matrix (zero for atoms not selected or fixed)
 */
AtomMatrix BasinHoppingJob::displaceRandomSelected(double maxDisplacement, const std::vector<int> &selectedAtoms)
{
    disp_count++;
    AtomMatrix displacement;
    displacement.resize(trial->numberOfAtoms(), 3);
    displacement.setZero();

    if (selectedAtoms.empty())
    {
        // Fall back to default behavior
        return displaceRandom(maxDisplacement);
    }

    VectorXd distvec = calculateDistanceFromCenter(current);

    // Mark only the selected atoms (no neighbors)
    std::vector<bool> atomsToDisplace(trial->numberOfAtoms(), false);

    for (size_t i = 0; i < selectedAtoms.size(); i++)
    {
        int atomIdx = selectedAtoms[i];
        if (atomIdx >= 0 && atomIdx < trial->numberOfAtoms())
        {
            atomsToDisplace[atomIdx] = true;
        }
    }

    // Apply displacement to selected atoms only
    for (int i = 0; i < trial->numberOfAtoms(); i++)
    {
        if (!atomsToDisplace[i] || trial->getFixed(i))
        {
            continue;
        }

        double dist = distvec(i);
        double disp = 0.0;

        // Calculate displacement magnitude (same as displaceRandom)
        if (parameters->basinHoppingDisplacementAlgorithm == "standard")
        {
            disp = maxDisplacement;
        }
        else if (parameters->basinHoppingDisplacementAlgorithm == "linear")
        {
            double Cs = maxDisplacement / distvec.maxCoeff();
            disp = Cs * dist;
        }
        else if (parameters->basinHoppingDisplacementAlgorithm == "quadratic")
        {
            double Cq = maxDisplacement / (distvec.maxCoeff() * distvec.maxCoeff());
            disp = Cq * dist * dist;
        }
        else
        {
            log("Unknown displacement_algorithm\n");
            std::exit(1);
        }

        // Apply random displacement
        for (int j = 0; j < 3; j++)
        {
            if (parameters->basinHoppingDisplacementDistribution == "uniform")
            {
                displacement(i, j) = randomDouble(2 * disp) - disp;
            }
            else if (parameters->basinHoppingDisplacementDistribution == "gaussian")
            {
                displacement(i, j) = gaussRandom(0.0, disp);
            }
            else
            {
                log("Unknown displacement_distribution\n");
                std::exit(1);
            }
        }
    }

    // Debug: Log displacement vectors for verification
    log("[Basin Hopping] DEBUG: Displacement vectors applied:\n");
    for (size_t i = 0; i < selectedAtoms.size(); i++)
    {
        int idx = selectedAtoms[i];
        if (idx >= 0 && idx < trial->numberOfAtoms())
        {
            double dx = displacement(idx, 0);
            double dy = displacement(idx, 1);
            double dz = displacement(idx, 2);
            double magnitude = sqrt(dx * dx + dy * dy + dz * dz);
            log("[Basin Hopping] DEBUG: Atom #%d (idx0=%d): displacement=[%.6f, %.6f, %.6f], magnitude=%.6f Å\n",
                idx, idx, dx, dy, dz, magnitude);
        }
    }

    // Verify: Check that non-selected atoms have zero displacement
    int nonZeroCount = 0;
    for (int i = 0; i < trial->numberOfAtoms(); i++)
    {
        // Check if atom i is in selectedAtoms
        bool isSelected = false;
        for (size_t j = 0; j < selectedAtoms.size(); j++)
        {
            if (selectedAtoms[j] == i)
            {
                isSelected = true;
                break;
            }
        }

        if (!isSelected && !trial->getFixed(i))
        {
            double dx = displacement(i, 0);
            double dy = displacement(i, 1);
            double dz = displacement(i, 2);
            double magnitude = sqrt(dx * dx + dy * dy + dz * dz);
            if (magnitude > 1e-10)
            {
                log("[Basin Hopping] DEBUG: WARNING: Atom #%d (idx0=%d, not selected) has non-zero displacement: %.10f Å\n",
                    i, i, magnitude);
                nonZeroCount++;
            }
        }
    }
    if (nonZeroCount == 0)
    {
        log("[Basin Hopping] DEBUG: Verification passed: All non-selected atoms have zero displacement\n");
    }
    else
    {
        log("[Basin Hopping] DEBUG: ERROR: %d non-selected atoms have non-zero displacement!\n", nonZeroCount);
    }

    return displacement;
}
