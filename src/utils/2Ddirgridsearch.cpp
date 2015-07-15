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

#include <cstdio>
#include <ctime>
#include <sbpl_cpr/utils/2Ddirgridsearch.h>
#include <sbpl_cpr/utils/heap.h>
#include <sbpl_cpr/utils/list.h>
#include <string.h>
#include <math.h>

// fitting direction open and closed flag in a single char
// still have room for 2 bits
#define DIRECTION_MASK 0x1F
#define CLOSED_FLAG    0x80

//    6     4
// 8  7  5  3   2
//    9  0  1
//10 11 13 15  16
//   12    14
static int dir_dx[17] = {0, 1, 2, 1, 1, 0, -1, -1, -2, -1, -2, -1, -1, 0, 1, 1, 2};
static int dir_dy[17] = {0, 0, 1, 1, 2, 1, 2, 1, 1, 0, -1, -1, -2, -1, -2, -1, -1};
#define SQ2 1.4142136
#define SQ5 2.2360680
static float dir_dist[17] = {0, 1, SQ5, SQ2, SQ5, 1, SQ5, SQ2, SQ5, 1, SQ5, SQ2, SQ5, 1, SQ5, SQ2, SQ5};

// reverse directions
static unsigned char rev_dir[17] = {0, 9, 10, 11, 12, 13, 14, 15, 16, 1, 2, 3, 4, 5, 6, 7, 8};

static float unit_dirs[17][2] = { {0,0},
{        1,         0}, { 0.894427,  0.447214}, { 0.707107,  0.707107}, { 0.447214,  0.894427},
{        0,         1}, {-0.447214,  0.894427}, {-0.707107,  0.707107}, {-0.894427,  0.447214},
{       -1,         0}, {-0.894427, -0.447214}, {-0.707107, -0.707107}, {-0.447214, -0.894427},
{        0,        -1}, { 0.447214, -0.894427}, { 0.707107, -0.707107}, { 0.894427, -0.447214} };

static unsigned char cache_friendly_dirs_no_origin[16] = {12, 14, 10, 11, 13, 15, 16, 9, 1, 8, 7, 5, 3, 2, 6, 4};

using namespace std;

// decide which of the 3 directions dot with the source vector to give the largest product
inline int max_dir(const int a, const int b, const int c, const int dx, const int dy)
{
  const float da = dx * unit_dirs[a][0] + dy * unit_dirs[a][1];
  const float db = dx * unit_dirs[b][0] + dy * unit_dirs[b][1];
  const float dc = dx * unit_dirs[c][0] + dy * unit_dirs[c][1];

  if(da > db)
  {
    if(da > dc)
      return a;
    return c;
  }

  if(dc > db)
  {
    if(dc > da)
      return c;
    return a;
  }

  return b;
}

namespace SBPL
{
// this function finds the best matching direction via dot products and
// a binary search
unsigned char dxdyToDir(int dx, int dy)
{
  if(dx == 0 && dy == 0)
    return 0;

  const float d1 = dx * unit_dirs[ 1][0] + dy * unit_dirs[ 1][1];
  const float d5 = dx * unit_dirs[ 5][0] + dy * unit_dirs[ 5][1];
  const float d9 = dx * unit_dirs[ 9][0] + dy * unit_dirs[ 9][1];
  const float d13= dx * unit_dirs[13][0] + dy * unit_dirs[13][1];

  if(d1 > d9)
  {
    if(d5 > d13)
    {
      if(d1 > d5)
      {
        // check 1 2 3
        return max_dir(1, 2, 3, dx, dy);
      }
      else
      {
        // check 3 4 5
        return max_dir(3, 4, 5, dx, dy);
      }
    }
    else
    {
      if(d1 > d13)
      {
        // 1 15 16
        return max_dir(1, 15, 16, dx, dy);
      }
      else
      {
        // 13 14 15
        return max_dir(13, 14, 15, dx, dy);
      }
    }
  }
  else // d9 > d1
  {
    if(d5 > d13)
    {
      if(d5 > d9)
      {
        // 5 6 7
        return max_dir(5, 6, 7, dx, dy);
      }
      else
      {
        // 7 8 9
        return max_dir(7, 8, 9, dx, dy);
      }
    }
    else
    {
      if(d9 > d13)
      {
        // 9 10 11
        return max_dir(9, 10, 11, dx, dy);
      }
      else
      {
        // 11 12 13
        return max_dir(11, 12, 13, dx, dy);
      }
    }
  }

  return 0;
}
}

//---------------------initialization and destruction routines--------------------------------------------------------
SBPL2DDirGridSearch::SBPL2DDirGridSearch(int width_x, int height_y, float cellsize_m, bool reverse_search)
{
    width_ = width_x;
    height_ = height_y;
    cellSize_m_ = cellsize_m;

    startX_ = -1;
    startY_ = -1;
    goalX_  = -1;
    goalY_  = -1;

    largestcomputedoptf_ = 0;

    term_condition_usedlast = SBPL_2DGRIDSEARCH_TERM_CONDITION_ALLCELLS;

    //allocate memory
    g_value = new size_t[width_ * height_];
    site_data = new unsigned char[width_ * height_];

    // cost of diagonal move in mm:
    int maxdistance = 2237 * 1000 * cellSize_m_;
    int numofbuckets = 255 * maxdistance;
    OPEN2DBLIST_ = new TSlidingBucket<unsigned char*>(numofbuckets);

    reverse_search_direction_ = reverse_search;
}

size_t SBPL2DDirGridSearch::xy2idx(size_t x, size_t y) const
{
    return x + width_ * y;
}

SBPL2DDirGridSearch::~SBPL2DDirGridSearch()
{
    // destroy the 2D states:
    if (g_value)
    {
      delete [] g_value;
      delete [] site_data;
    }
    g_value = 0;

    if (OPEN2DBLIST_ != NULL) {
        delete OPEN2DBLIST_;
        OPEN2DBLIST_ = NULL;
    }
}

void SBPL2DDirGridSearch::setInfGValues()
{
  for (int i=0; i<width_*height_; i++)
    g_value[i] = INFINITECOST;
}


void SBPL2DDirGridSearch::site_data_offset(const unsigned char* c, size_t& idx)
{
  std::ptrdiff_t diff = c - site_data;
  idx = static_cast<size_t>(diff);
}

void SBPL2DDirGridSearch::site_data_coordinates(const unsigned char* c, size_t& x, size_t& y)
{
    std::ptrdiff_t diff = c - site_data;
    const size_t offset = static_cast<size_t>(diff);

    x = offset % width_;
    y = (offset - x) / width_;
}

//-----------------------------------------main functions--------------------------------------------------------------

bool SBPL2DDirGridSearch::search(
            int startx_c, int starty_c, int goalx_c, int goaly_c,
            unsigned char (*cost_func)(size_t, size_t, int, int, void*),
            void* cost_data,
            unsigned char obsthresh,
            SBPL_2DGRIDSEARCH_TERM_CONDITION termination_condition)
{
    bzero(site_data, width_ * height_ * sizeof(unsigned char));
#if DEBUG
#ifndef ROS
    const char* f2dgriddebug = "2dgriddebug.txt";
#endif
    FILE* f2Dsearch = SBPL_FOPEN(f2dgriddebug, "w");
#endif

    //init start and goal coordinates
    startX_ = startx_c;
    startY_ = starty_c;
    goalX_ = goalx_c;
    goalY_ = goaly_c;

    int numofExpands = 0;

    //check the validity of start/goal
    if (!withinMap(startx_c, starty_c) || !withinMap(goalx_c, goaly_c)) {
        SBPL_ERROR("ERROR: grid2DZoneSearch is called on invalid start (%d %d) or goal(%d %d)\n", startx_c, starty_c,
                   goalx_c, goaly_c);
#if DEBUG
        SBPL_FCLOSE(f2Dsearch);
#endif
        return false;
    }

    //reset OPEN
    OPEN2DBLIST_->reset();

    //set the term. condition
    term_condition_usedlast = termination_condition;

    // initialize the start and goal states
    const size_t start_idx = xy2idx(startX_, startY_);
    const size_t  goal_idx = xy2idx( goalX_,  goalY_);

    setInfGValues(); // setting all g values to infinity.

    g_value[ start_idx ] = 0;
    OPEN2DBLIST_->insert(&site_data[ start_idx ],  g_value[ start_idx ]);

    //set the termination condition
    float term_factor = 0.0;
    switch (termination_condition) {
    case SBPL_2DGRIDSEARCH_TERM_CONDITION_OPTPATHFOUND:
        term_factor = 1;
        break;
    case SBPL_2DGRIDSEARCH_TERM_CONDITION_20PERCENTOVEROPTPATH:
        term_factor = (float)(1.0 / 1.2);
        break;
    case SBPL_2DGRIDSEARCH_TERM_CONDITION_TWOTIMESOPTPATH:
        term_factor = 0.5;
        break;
    case SBPL_2DGRIDSEARCH_TERM_CONDITION_THREETIMESOPTPATH:
        term_factor = (float)(1.0 / 3.0);
        break;
    case SBPL_2DGRIDSEARCH_TERM_CONDITION_ALLCELLS:
        term_factor = 0.0;
        break;
    default:
        SBPL_ERROR("ERROR: incorrect termination factor for grid2DZoneSearch\n");
        term_factor = 0.0;
    };

    //get the current time
    clock_t starttime = clock();

    //the main repetition of expansions
    SBPL_PRINTF("2D zone search with sliding buckets and term_factor=%.3f\n", term_factor);

    // --------------------------- Main Search Loop --------------------------
    while (!OPEN2DBLIST_->empty() && g_value[ goal_idx ] > term_factor * OPEN2DBLIST_->getminkey())
    {
#if DEBUG
        SBPL_FPRINTF(f2Dsearch, "currentminelement_priority before pop=%d\n", OPEN2DBLIST_->getminkey());
#endif

        //get the next state for expansion
        const unsigned char* searchExpState = OPEN2DBLIST_->popminelement();

        size_t exp_x, exp_y;
        site_data_coordinates(searchExpState, exp_x, exp_y);
        const size_t exp_idx = xy2idx(exp_x, exp_y);

        //close the state if it wasn't yet
        site_data[exp_idx] |= CLOSED_FLAG;

#if DEBUG
        SBPL_FPRINTF(f2Dsearch, "expanding state <%d %d> with g=%d "
                     "(currentminelement_priority=%d, currentfirstbucket_bindex=%d, currentfirstbucket_priority=%d)\n",
                     exp_x, exp_y, g_value[exp_idx], OPEN2DBLIST_->getminkey(),
                     OPEN2DBLIST_->currentfirstbucket_bindex, OPEN2DBLIST_->currentfirstbucket_priority);
#endif

        //expand
        numofExpands++;

        //iterate over successors
        for (size_t _d=0; _d<16; _d++)
        {
            size_t dir = cache_friendly_dirs_no_origin[_d];
            const int newx = exp_x + dir_dx[dir];
            const int newy = exp_y + dir_dy[dir];

            //make sure it is inside the map
            if (!withinMap(newx, newy)) continue;

            const size_t new_idx = xy2idx(newx, newy);

            //make sure it is not closed
            if(site_data[new_idx] & CLOSED_FLAG) continue;

            //compute the cost - using user supplied function
            const size_t mapcost = cost_func(newx, newy,
                                             reverse_search_direction_?-dir_dx[dir]:dir_dx[dir],
                                             reverse_search_direction_?-dir_dy[dir]:dir_dy[dir],
                                             cost_data);

            //check for obstacle
            if (mapcost >= obsthresh) continue;

            // interpretting cost as cost per mm
            const int cost = (mapcost + 1) * dir_dist[dir] * cellSize_m_ * 1000.0;

            // update if this move is better
            const size_t new_value = cost + g_value[ exp_idx ];
            if ( new_value < g_value[ new_idx ] )
            {
                g_value[ new_idx ] = new_value;

                // point to source (not needed - this is actually a flood fill algorithm)
                //site_data[ new_idx ] &= ~DIRECTION_MASK;
                //site_data[ new_idx ] |= (rev_dir[dir] & DIRECTION_MASK);
#if DEBUG
                SBPL_FPRINTF(f2Dsearch, "inserting state <%d %d> with g=%d\n", newx, newy, g_value[ new_idx ]);
#endif

                //put it into the list
                OPEN2DBLIST_->insert(&site_data[ new_idx ],  g_value[ new_idx ]);
            }
        } //over directions
    }//while

    //set lower bounds for the remaining states
    if (!OPEN2DBLIST_->empty())
    {
        unsigned char* m = OPEN2DBLIST_->getminelement();
        size_t idx;
        site_data_offset(m, idx);
        largestcomputedoptf_ = g_value[idx];
    }
    else
    {
        largestcomputedoptf_ = INFINITECOST;
    }

    SBPL_PRINTF( "# of expands during 2dzonegridsearch=%d time=%d msecs 2Dsolcost_inmm=%d "
                "largestoptfval=%d (start=%d %d goal=%d %d)\n",
                numofExpands, (int)(((clock() - starttime) / (double)CLOCKS_PER_SEC) * 1000),
                searchStates2D_[goalx_c][goaly_c].g, largestcomputedoptf_, startx_c, starty_c, goalx_c, goaly_c);

#if DEBUG
    SBPL_FCLOSE(f2Dsearch);
#endif
    return true;
}
