#ifndef __2DDIRGRIDSEARCH_H_
#define __2DDIRGRIDSEARCH_H_

#include <cstdlib>
#include <sbpl_cpr/planners/planner.h>
#include <sbpl_cpr/utils/key.h>
#include <sbpl_cpr/utils/utils.h>
#include <sbpl_cpr/utils/list.h>
#include <sbpl_cpr/utils/2Dgridsearch.h>

#define SBPL_2DGRIDSEARCH_HEUR2D(x,y)  ((int)(1000*cellSize_m_*__max(abs(x-goalX_),abs(y-goalY_))))

/**
 * \brief Custom 2D search with cost scaling based on direction
 */
class SBPL2DDirGridSearch
{
public:
    /**
     * @brief SBPL2DDirGridSearch Create a search space for 2D grids
     * @param width_x grid width
     * @param height_y grid height
     * @param cellsize_m resolution
     */
    SBPL2DDirGridSearch(int width_x, int height_y, float cellsize_m, bool reverse_search);
    ~SBPL2DDirGridSearch();

    /**
     * @brief search Flood fill from source to destination to be used as a distance heuristic
     * @param startx_c search source x cell coordinate
     * @param starty_c search source y cell coordinate
     * @param goalx_c  search destination x cell coordinate
     * @param goaly_c  search destination y cell coordinate
     * @param cost_func cost function used to scale base cose based on direction
     * @param cost_data arbitrary extra data to be supplied to cost function
     * @param obsthresh obstacle treshold value
     * @param termination_condition SBPL termination condition
     * @return
     */
    bool search(int startx_c, int starty_c, int goalx_c, int goaly_c,
                unsigned char (*cost_func)(size_t, size_t, int, int, void*),
                void* cost_data,
                unsigned char obsthresh,
                SBPL_2DGRIDSEARCH_TERM_CONDITION termination_condition);

    /**
     * \brief print all the values
     */
    void printvalues();

    /**
     * \brief returns the computed distance from the start to <x,y>. If not computed, then returns lower bound on it.
     */
    inline int getlowerboundoncostfromstart_inmm(int x, int y)
    {
        if (!withinMap(x,y))
          return largestcomputedoptf_;
        const size_t idx = xy2idx(x,y);

        if (term_condition_usedlast == SBPL_2DGRIDSEARCH_TERM_CONDITION_OPTPATHFOUND) {
            //heuristic search
            int h = SBPL_2DGRIDSEARCH_HEUR2D(x,y);
            //the logic is that if s wasn't expanded, then g(s) + h(s) >=
            //maxcomputed_fval => g(s) >= maxcomputed_fval - h(s)
            return ((
                    g_value[idx] + h <= largestcomputedoptf_) ? g_value[idx] :
                    largestcomputedoptf_ < INFINITECOST       ? largestcomputedoptf_ - h :
                                                                             INFINITECOST);
        }
        else {
            //Dijkstra's search
            //the logic is that if s wasn't expanded, then g(s) >= maxcomputed_fval => g(s) >= maxcomputed_fval - h(s)
            return std::min( g_value[idx], largestcomputedoptf_ );
        }
    }

    /**
     * \brief returns largest optimal g-value computed by search - a lower
     *        bound on the state values of unexpanded states
     */
    size_t getlargestcomputedoptimalf_inmm() { return largestcomputedoptf_; }

    size_t xy2idx(size_t x, size_t y) const;

private:
    inline bool withinMap(int x, int y)
    {
        return (x >= 0 && y >= 0 && x < width_ && y < height_);
    }

    void setInfGValues();

    //2D search data
    TSlidingBucket<unsigned char*>* OPEN2DBLIST_;

    // buffer to hold the costs from the source to all other visited points
    size_t* g_value;

    // flags for each site
    unsigned char* site_data;

    bool reverse_search_direction_;

    void site_data_coordinates(const unsigned char* c, size_t& x, size_t& y);
    void site_data_offset(const unsigned char* c, size_t& idx);

    //start and goal configurations
    int startX_, startY_;
    int goalX_, goalY_;

    //map parameters
    int width_, height_;
    float cellSize_m_;

    //largest optimal g-value computed by search
    size_t largestcomputedoptf_;

    //termination criterion used in the search
    SBPL_2DGRIDSEARCH_TERM_CONDITION term_condition_usedlast;
};

namespace SBPL
{
  unsigned char dxdyToDir(int dx, int dy);
}

#endif
