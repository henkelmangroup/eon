
#ifndef STRUCTURECOMPARISONJOB_H
#define STRUCTURECOMPARISONJOB_H

#include "Job.h"
#include "Parameters.h"
#include "Matter.h"

class StructureComparisonJob: public Job {
    public:
        StructureComparisonJob(Parameters *params);
        ~StructureComparisonJob(void);
        std::vector<std::string> run(void);
        bool ComparePositions(Matter *matt1, Matter *matt2);
    private:
        Parameters *parameters;
};

#endif
