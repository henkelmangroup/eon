#include "StructureComparisonJob.h"
#include "HelperFunctions.h"
#include "Log.h"
#include "Matter.h"
#include "Optimizer.h"

StructureComparisonJob::StructureComparisonJob(Parameters *params) {
    parameters = params;
}

StructureComparisonJob::~StructureComparisonJob() { }

std::vector<std::string> StructureComparisonJob::run(void) {
    std::vector<std::string> returnFiles;
    Matter *matter1 = new Matter(parameters);
    matter1->con2matter("matter1.con");
    return returnFiles;
}

bool StructureComparisonJob::ComparePositions(Matter *matt1, Matter *matt2){
    bool same = false;
    AtomMatrix pos1 = matt1->getPositions();
    AtomMatrix pos2 = matt2->getPositions();
    return same;
}
