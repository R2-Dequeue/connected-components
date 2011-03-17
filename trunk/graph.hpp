/*!
 * \file
 * \author Chris de Pujo
 */

#ifndef __GRAPH__
#define __GRAPH__

#include <vector>
#include <algorithm>

#include <ginac/ginac.h>

template <typename T, typename Compare>
class Graph
{
public:

    typename unsigned int GraphIndex; // graphIndex

    Graph(std::vector<T> & vertexSet, Compare comp);

    //bool areAdjacent(const GraphIndex i, const GraphIndex j) const;
    bool areConnected(const GraphIndex i, const GraphIndex j) const;

private:

    std::vector< std::vector<GraphIndex> > partition;
};

////////////////////////////////////////////////////////////////////////////////

template <typename T, typename Compare>
inline Graph::Graph(std::vector<T> & set, Compare comp)
{
    typedef unsigned int uint;
    const uint n = vertexSet.size();

    GiNaC::matrix adjMatrix(n);

    for (uint j = 0; j < n; ++j)
        for (uint i = j; i < n; ++i)
            if (comp(i, j))
            {
                adjMatrix(i, j) = 1;
                adjMatrix(j, i) = 1;
            }

    adjMatrix = adjMatrix.pow(n);

    std::vector<GraphIndex> cells;
    cells.reserve(n);

    for (uint i = 0; i < n; ++i)
    {
        cells.clear();

        for (GraphIndex j = 0; j < n; ++j)
            if (adjMatrix(i, j) != 0)
                cells.push_back(j);

        uint k = 0;
        for (k = 0; k < partitions.size(); ++k)
            if (partitions[k][0] == cells[0])
                break;

        if (k >= partitions.size())
            partitions.push_back(cells);
    }
}

//bool Graph::areAdjacent(const GraphIndex i, const GraphIndex j) const

bool Graph::areConnected(GraphIndex i, GraphIndex j) const
{
    std::vector< std::vector<GraphIndex> >::iterator p, e, q;

    if (i > j)
        std::swap(i, j);

    for (p = partitions.begin(), e = partitions.end(); p != e; ++p)
    {
        q = std::find(p->begin(), p->end(), i);

        if (q != p->end())
        {
            if (std::find(q, p->end(), j) != p->end()) // should work if i == j.
                return true;
            else
                return false;
        }
    }

    assert(false);

    return false;
}

#endif // __GRAPH__
