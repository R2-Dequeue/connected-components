/*!
 * \file
 * \author Chris de Pujo
 */

#include <cassert>
#include <algorithm>

#include "algebraic.hpp"

template <class T, class Allocator = std::allocator<T>>
class UniqueSorter
{
private:

    struct ModifyingSetCompare {
        bool operator()(T::value_type * const lhs, T::value_type * const rhs)
        {
            return (*lhs < *rhs);
        }
    }

    //typename T::value_type * * ptrs;
    std::set<typename T::value_type *, ModifyingSetCompare, Allocator> ptrs;

public:

    typedef T::value_type       value_type;
    typedef Allocator           allocator_type;

    UniqueSorter();
    UniqueSorter(const T & container);
    ~UniqueSorter();
}; // std::vector<T*, boost::pool_allocator> ptrs;

inline UniqueSorter::UniqueSorter() : ptrs(NULL) {}

inline UniqueSorter::UniqueSorter(const T & container)
{
    unsigned int i = 0;

    ptrs = new (typename T::value_type) [container.size()];

    for (T::const_iterator
         e = container.begin(), end = container.end(); e != end; ++e)
        ptrs[i++] = &(*e);
}

inline UniqueSorter::~UniqueSorter()
{
    delete [] ptrs;
    ptrs = NULL;
}

// Sort, make unique
// want normal non-modifying sort (stl), plus a modifying sort for Algebraic
// objects.
