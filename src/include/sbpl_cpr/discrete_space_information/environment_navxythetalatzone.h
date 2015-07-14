#ifndef __ENVIRONMENT_NAVXYTHETALATZONE_H_
#define __ENVIRONMENT_NAVXYTHETALATONE_H_

#include <sbpl_cpr/discrete_space_information/environment_navxythetalat.h>

class SBPL2DZoneGridSearch;


class EnvironmentNAVXYTHETALATZone : public EnvironmentNAVXYTHETALAT
{
public:
    EnvironmentNAVXYTHETALATZone();
    ~EnvironmentNAVXYTHETALATZone();

    /**
     * \brief see comments on the same function in the parent class
     */
    virtual int GetGoalHeuristic(int stateID);

    /**
     * \brief see comments on the same function in the parent class
     */
    virtual int GetStartHeuristic(int stateID);

    /**
     * \brief see comments on the same function in the parent class
     */
    virtual void EnsureHeuristicsUpdated(bool bGoalHeuristics);

    /**
     * @brief SetSearchCostFunction
     * @param cost_func pointer to function that generates the direction dependant costs
     * @param cost_data arbitrary data to support the cost decisions
     */
    virtual void SetSearchCostFunction(unsigned char (*cost_func)(size_t, size_t, unsigned char, void*), void* cost_data);

protected:
    SBPL2DZoneGridSearch* grid2Dsearchfromstart_zone; //computes h-values that provide distances from start x,y to all cells
    SBPL2DZoneGridSearch* grid2Dsearchfromgoal_zone; //computes h-values that provide distances to goal x,y from all cells

    virtual int GetActionCost(int SourceX, int SourceY, int SourceTheta, EnvNAVXYTHETALATAction_t* action);

    virtual void ComputeHeuristicValues();
    virtual void PrintHeuristicValues();

    unsigned char (*fwd_cost_func_)(size_t, size_t, unsigned char, void*);
    void* fwd_cost_data_;
};

#endif
