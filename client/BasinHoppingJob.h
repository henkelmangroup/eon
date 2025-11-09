#ifndef BASINHOPPINGJOB_H
#define BASINHOPPINGJOB_H

#include "Job.h"
#include "Matter.h"
#include "Parameters.h"

class BasinHoppingJob : public Job {
    public:
        BasinHoppingJob(Parameters *params);
        ~BasinHoppingJob(void);
        std::vector<std::string> run(void);

    private:
        VectorXd calculateDistanceFromCenter(Matter *matter);
        AtomMatrix displaceRandom(double maxDisplacement);
        AtomMatrix displaceRandomSelected(double maxDisplacement, const std::vector<int>& selectedAtoms);
        void randomSwap(Matter *matter);
        std::vector<int> parseAtomList(const std::string& atomListStr, int nAtoms);
        std::vector<int> getAtomsByType(const std::string& typeListStr, Matter *matter);
        std::vector<int> selectAtomsForDisplacement();
        Parameters *parameters;
        Matter *current;
        Matter *trial; // initial configuration
        vector<long> getElements(Matter *matter);
        std::vector<std::string> returnFiles;
        int jump_count; // count of jump moves
        int disp_count; // count of displacement moves
        int swap_count; // count of swap moves
        int fcalls;

        std::vector<Matter *> uniqueStructures;
        std::vector<double> uniqueEnergies;
        
        // Cached atom selection lists (computed once)
        std::vector<int> cachedAtomList;
        std::vector<int> cachedTypeList;
        bool atomListsCached;
};

#endif
