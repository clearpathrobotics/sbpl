#ifndef __ENVIRONMENT_NAVXYTHETADIRLAT_H_
#define __ENVIRONMENT_NAVXYTHETADIRLAT_H_

#include <sbpl_cpr/discrete_space_information/environment_navxythetalat.h>

class SBPL2DDirGridSearch;

/**
 * @brief The EnvironmentNAVXYTHETADIRLAT class
 * This extends EnvironmentNAVXYTHETALAT with a custom 2D search that
 * considers direction of motion and a custom GetActionCost that also
 * considers the direction of motion.
 */
class EnvironmentNAVXYTHETADIRLAT : public EnvironmentNAVXYTHETALAT
{
public:
    EnvironmentNAVXYTHETADIRLAT();
    ~EnvironmentNAVXYTHETADIRLAT();

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
    virtual void SetSearchCostFunction(unsigned char (*cost_func)(size_t, size_t, int, int, void*), void* cost_data);

protected:
    SBPL2DDirGridSearch* grid2Dsearchfromstart_dir; //computes h-values that provide distances from start x,y to all cells
    SBPL2DDirGridSearch* grid2Dsearchfromgoal_dir; //computes h-values that provide distances to goal x,y from all cells

    virtual int GetActionCost(int SourceX, int SourceY, int SourceTheta, EnvNAVXYTHETALATAction_t* action);

    virtual void ComputeHeuristicValues();
    virtual void PrintHeuristicValues();

    unsigned char (*fwd_cost_func_)(size_t, size_t, int, int, void*);
    void* fwd_cost_data_;
};

#endif
