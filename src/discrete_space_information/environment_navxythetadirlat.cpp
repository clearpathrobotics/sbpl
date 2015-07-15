/*
 * Copyright (c) 2008, Maxim Likhachev
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * Neither the name of the University of Pennsylvania nor the names of its
 *       contributors may be used to endorse or promote products derived from
 *       this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */

#include <cmath>
#include <cstring>
#include <ctime>
#include <sbpl_cpr/discrete_space_information/environment_navxythetadirlat.h>
#include <sbpl_cpr/utils/2Ddirgridsearch.h>
#include <sbpl_cpr/utils/key.h>
#include <sbpl_cpr/utils/mdp.h>
#include <sbpl_cpr/utils/mdpconfig.h>

using namespace std;

#if TIME_DEBUG
static clock_t time3_addallout = 0;
static clock_t time_gethash = 0;
static clock_t time_createhash = 0;
static clock_t time_getsuccs = 0;
#endif

static long int checks = 0;

#define XYTHETA2INDEX(X,Y,THETA) (THETA + X*EnvNAVXYTHETALATCfg.NumThetaDirs + \
                                  Y*EnvNAVXYTHETALATCfg.EnvWidth_c*EnvNAVXYTHETALATCfg.NumThetaDirs)

//-----------------constructors/destructors-------------------------------

EnvironmentNAVXYTHETADIRLAT::EnvironmentNAVXYTHETADIRLAT()
 : EnvironmentNAVXYTHETALAT()
{
  fwd_cost_func_ = NULL;
  fwd_cost_data_ = NULL;
}

EnvironmentNAVXYTHETADIRLAT::~EnvironmentNAVXYTHETADIRLAT()
{
    SBPL_PRINTF("destroying EnvironmentNAVXYTHETADIRLAT\n");
    if (grid2Dsearchfromstart_dir != NULL) delete grid2Dsearchfromstart_dir;
    grid2Dsearchfromstart_dir = NULL;

    if (grid2Dsearchfromgoal_dir != NULL) delete grid2Dsearchfromgoal_dir;
    grid2Dsearchfromgoal_dir = NULL;

}

// fall-back cost function if none is provided
static unsigned char default_cost_function(size_t x, size_t y, int /*dx*/, int /*dy*/, void* data)
{
    EnvironmentNAVXYTHETADIRLAT* env = (EnvironmentNAVXYTHETADIRLAT*)data;
    return env->GetMapCost(x, y);
}

//------------------------------Heuristic computation--------------------------
void EnvironmentNAVXYTHETADIRLAT::EnsureHeuristicsUpdated(bool bGoalHeuristics)
{
    if (bNeedtoRecomputeStartHeuristics && !bGoalHeuristics) {
        grid2Dsearchfromstart_dir->search(
                                      EnvNAVXYTHETALATCfg.StartX_c, EnvNAVXYTHETALATCfg.StartY_c,
                                      EnvNAVXYTHETALATCfg.EndX_c, EnvNAVXYTHETALATCfg.EndY_c,
                                      (fwd_cost_func_ && fwd_cost_data_) ? fwd_cost_func_ : default_cost_function,
                                      (fwd_cost_func_ && fwd_cost_data_) ? fwd_cost_data_ : (void*)this,
                                      EnvNAVXYTHETALATCfg.cost_inscribed_thresh,
                                      SBPL_2DGRIDSEARCH_TERM_CONDITION_TWOTIMESOPTPATH);
        bNeedtoRecomputeStartHeuristics = false;
        SBPL_PRINTF("2dsolcost_infullunits=%d\n",
                    (int)(grid2Dsearchfromstart_dir->getlowerboundoncostfromstart_inmm(EnvNAVXYTHETALATCfg.EndX_c,
                                                                                   EnvNAVXYTHETALATCfg.EndY_c) /
                          EnvNAVXYTHETALATCfg.nominalvel_mpersecs));

    }

    if (bNeedtoRecomputeGoalHeuristics && bGoalHeuristics) {
        grid2Dsearchfromgoal_dir->search(
                                     EnvNAVXYTHETALATCfg.EndX_c, EnvNAVXYTHETALATCfg.EndY_c,
                                     EnvNAVXYTHETALATCfg.StartX_c, EnvNAVXYTHETALATCfg.StartY_c,
                                     (fwd_cost_func_ && fwd_cost_data_) ? fwd_cost_func_ : default_cost_function,
                                     (fwd_cost_func_ && fwd_cost_data_) ? fwd_cost_data_ : (void*)this,
                                     EnvNAVXYTHETALATCfg.cost_inscribed_thresh,
                                     SBPL_2DGRIDSEARCH_TERM_CONDITION_TWOTIMESOPTPATH);
        bNeedtoRecomputeGoalHeuristics = false;
        SBPL_PRINTF("2dsolcost_infullunits=%d\n",
                    (int)(grid2Dsearchfromgoal_dir->getlowerboundoncostfromstart_inmm(EnvNAVXYTHETALATCfg.StartX_c,
                                                                                  EnvNAVXYTHETALATCfg.StartY_c) /
                          EnvNAVXYTHETALATCfg.nominalvel_mpersecs));
    }
}


void EnvironmentNAVXYTHETADIRLAT::SetSearchCostFunction(unsigned char (*cost_func)(size_t, size_t, int, int, void*), void* cost_data)
{
    fwd_cost_func_ = cost_func;
    fwd_cost_data_ = cost_data;
}

void EnvironmentNAVXYTHETADIRLAT::ComputeHeuristicValues()
{
    //whatever necessary pre-computation of heuristic values is done here
    SBPL_PRINTF("Precomputing heuristics...\n");

    const bool reversed_search = true;
    //allocated 2D grid searches
    grid2Dsearchfromstart_dir = new SBPL2DDirGridSearch(EnvNAVXYTHETALATCfg.EnvWidth_c, EnvNAVXYTHETALATCfg.EnvHeight_c,
                                                 (float)EnvNAVXYTHETALATCfg.cellsize_m, !reversed_search);
    grid2Dsearchfromgoal_dir = new SBPL2DDirGridSearch(EnvNAVXYTHETALATCfg.EnvWidth_c, EnvNAVXYTHETALATCfg.EnvHeight_c,
                                                 (float)EnvNAVXYTHETALATCfg.cellsize_m,  reversed_search);

    SBPL_PRINTF("done\n");
}

void EnvironmentNAVXYTHETADIRLAT::PrintHeuristicValues()
{
#ifndef ROS
    const char* heur = "heur.txt";
#endif
    FILE* fHeur = SBPL_FOPEN(heur, "w");
    if (fHeur == NULL) {
        SBPL_ERROR("ERROR: could not open debug file to write heuristic\n");
        throw new SBPL_Exception();
    }
    SBPL2DDirGridSearch* grid2Dsearch = NULL;

    for (int i = 0; i < 2; i++) {
        if (i == 0 && grid2Dsearchfromstart_dir != NULL) {
            grid2Dsearch = grid2Dsearchfromstart_dir;
            SBPL_FPRINTF(fHeur, "start heuristics:\n");
        }
        else if (i == 1 && grid2Dsearchfromgoal_dir != NULL) {
            grid2Dsearch = grid2Dsearchfromgoal_dir;
            SBPL_FPRINTF(fHeur, "goal heuristics:\n");
        }
        else
            continue;

        for (int y = 0; y < EnvNAVXYTHETALATCfg.EnvHeight_c; y++) {
            for (int x = 0; x < EnvNAVXYTHETALATCfg.EnvWidth_c; x++) {
                if (grid2Dsearch->getlowerboundoncostfromstart_inmm(x, y) < INFINITECOST)
                    SBPL_FPRINTF(fHeur, "%5d ", grid2Dsearch->getlowerboundoncostfromstart_inmm(x, y));
                else
                    SBPL_FPRINTF(fHeur, "XXXXX ");
            }
            SBPL_FPRINTF(fHeur, "\n");
        }
    }
    SBPL_FCLOSE(fHeur);
}

int EnvironmentNAVXYTHETADIRLAT::GetGoalHeuristic(int stateID)
{
#if USE_HEUR==0
    return 0;
#endif

#if DEBUG
    if (stateID >= (int)StateID2CoordTable.size()) {
        SBPL_ERROR("ERROR in EnvNAVXYTHETALAT... function: stateID illegal\n");
        throw new SBPL_Exception();
    }
#endif

    EnvNAVXYTHETALATHashEntry_t* HashEntry = StateID2CoordTable[stateID];
    //computes distances from start state that is grid2D, so it is EndX_c EndY_c
    int h2D = grid2Dsearchfromgoal_dir->getlowerboundoncostfromstart_inmm(HashEntry->X, HashEntry->Y);
    int hEuclid = (int)(NAVXYTHETALAT_COSTMULT_MTOMM * EuclideanDistance_m(HashEntry->X, HashEntry->Y,
                                                                           EnvNAVXYTHETALATCfg.EndX_c,
                                                                           EnvNAVXYTHETALATCfg.EndY_c));

    //define this function if it is used in the planner (heuristic backward search would use it)
    return (int)(((double)__max(h2D, hEuclid)) / EnvNAVXYTHETALATCfg.nominalvel_mpersecs);
}

int EnvironmentNAVXYTHETADIRLAT::GetStartHeuristic(int stateID)
{
#if USE_HEUR==0
    return 0;
#endif

#if DEBUG
    if (stateID >= (int)StateID2CoordTable.size()) {
        SBPL_ERROR("ERROR in EnvNAVXYTHETALAT... function: stateID illegal\n");
        throw new SBPL_Exception();
    }
#endif

    EnvNAVXYTHETALATHashEntry_t* HashEntry = StateID2CoordTable[stateID];
    int h2D = grid2Dsearchfromstart_dir->getlowerboundoncostfromstart_inmm(HashEntry->X, HashEntry->Y);
    int hEuclid = (int)(NAVXYTHETALAT_COSTMULT_MTOMM * EuclideanDistance_m(EnvNAVXYTHETALATCfg.StartX_c,
                                                                           EnvNAVXYTHETALATCfg.StartY_c, HashEntry->X,
                                                                           HashEntry->Y));

    //define this function if it is used in the planner (heuristic backward search would use it)
    return (int)(((double)__max(h2D, hEuclid)) / EnvNAVXYTHETALATCfg.nominalvel_mpersecs);
}

int EnvironmentNAVXYTHETADIRLAT::GetActionCost(int SourceX, int SourceY, int SourceTheta,
                                                EnvNAVXYTHETALATAction_t* action)
{
    sbpl_2Dcell_t cell;
    sbpl_xy_theta_cell_t interm3Dcell;
    int i;

    //TODO - go over bounding box (minpt and maxpt) to test validity and skip
    //testing boundaries below, also order intersect cells so that the four
    //farthest pts go first

    if (!IsValidCell(SourceX, SourceY)) return INFINITECOST;
    if (!IsValidCell(SourceX + action->dX, SourceY + action->dY)) return INFINITECOST;

    if (EnvNAVXYTHETALATCfg.Grid2D[SourceX + action->dX][SourceY + action->dY] >=
        EnvNAVXYTHETALATCfg.cost_inscribed_thresh)
    {
        return INFINITECOST;
    }

    //need to iterate over discretized center cells and compute cost based on them
    unsigned char maxcellcost = 0;
    const size_t size = action->interm3DcellsV.size();
    if(size)
    {
        int dx, dy; // deltas
        int px, py; // positions
        sbpl_xy_theta_cell_t* icells = &(action->interm3DcellsV.at(0));

        for (i = 0; i < size; i++)
        {
            if(i < size-1)
            {
              dx = icells[i+1].x - icells[i].x;
              dy = icells[i+1].y - icells[i].y;
            }
            else
            {
              if(size > 1)
              {
                dx = icells[i].x - icells[i-1].x;
                dy = icells[i].y - icells[i-1].y;
              }
              else
              {
                dx = 0;
                dy = 0;
              }
            }

            px = SourceX + icells[i].x;
            py = SourceY + icells[i].y;

            if (px < 0 || px >= EnvNAVXYTHETALATCfg.EnvWidth_c ||
                py < 0 || py >= EnvNAVXYTHETALATCfg.EnvHeight_c)
                return INFINITECOST;

            unsigned char this_cell_cost;

            if(fwd_cost_func_ && fwd_cost_data_)
            {
                this_cell_cost = fwd_cost_func_(px, py, dx, dy, fwd_cost_data_);
            }
            else
            {
                this_cell_cost = EnvNAVXYTHETALATCfg.Grid2D[px][py];
            }

            maxcellcost = __max(maxcellcost, this_cell_cost);

            //check that the robot is NOT in the cell at which there is no valid orientation
            if (maxcellcost >= EnvNAVXYTHETALATCfg.cost_inscribed_thresh) return INFINITECOST;
        }
    }

    //check collisions that for the particular footprint orientation along the action
    if (EnvNAVXYTHETALATCfg.FootprintPolygon.size() > 1 && (int)maxcellcost >=
        EnvNAVXYTHETALATCfg.cost_possibly_circumscribed_thresh)
    {
        checks++;

        //for (i = 0; i < (int)action->intersectingcellsV.size(); i++) {
        std::vector<sbpl_2Dcell_t>::const_iterator it;
        std::vector<sbpl_2Dcell_t>::const_iterator vend = action->intersectingcellsV.begin();
        for(it=action->intersectingcellsV.begin(); it != vend; ++it)
        {
            const int cell_x = (*it).x + SourceX;
            const int cell_y = (*it).y + SourceY;

            //check validity
            if (!IsValidCell(cell_x, cell_y)) return INFINITECOST;

            //if(EnvNAVXYTHETALATCfg.Grid2D[cell.x][cell.y] > currentmaxcost)
            ////cost computation changed: cost = max(cost of centers of the
            //robot along action)
            //	currentmaxcost = EnvNAVXYTHETALATCfg.Grid2D[cell.x][cell.y];
            //	//intersecting cells are only used for collision checking
        }
    }

    //to ensure consistency of h2D:
    maxcellcost = __max(maxcellcost, EnvNAVXYTHETALATCfg.Grid2D[SourceX][SourceY]);
    int currentmaxcost =
            (int)__max(maxcellcost, EnvNAVXYTHETALATCfg.Grid2D[SourceX + action->dX][SourceY + action->dY]);

    return action->cost * (currentmaxcost + 1); //use cell cost as multiplicative factor
}
