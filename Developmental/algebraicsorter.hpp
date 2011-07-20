/*!
 * \file
 * \author Chris de Pujo
 */

#include <cassert>
#include <algorithm>

#include "algebraic.hpp"

class AlgebraicSorter
{
public:

    typedef Algebraic value_type;
    // Need to change to a pool allocator.
    typedef std::allocator<Algebraic *> allocator_type;
    typedef std::iterator iterator;
    //typedef std::iterator const_iterator;
    typedef unsigned int size_type;

    //AlgebraicSorter();
    AlgebraicSorter(std::vector<Algebraic> & alphas);

    size_type size() const;
    bool empty() const;

    iterator begin();
    //const_iterator begin() const;
    iterator end();
    //const_iterator end() const;

private:

    struct ModifyingSetCompare;

    std::set<Algebraic *, ModifyingSetCompare, allocator_type> ptrs;
};

inline AlgebraicSorter::AlgebraicSorter(std::vector<Algebraic> & alphas)
{
    for (std::vector<Algebraic>::iterator
         alpha = alphas.begin(), e = alphas.end(); alpha != e; ++alpha)
        ptrs.insert( &(*alpha) );
}

inline AlgebraicSorter::size_type AlgebraicSorter::size() const
{
    return ptrs.size();
}

inline bool AlgebraicSorter::empty() const
{
    return ptrs.empty();
}

inline AlgebraicSorter::iterator AlgebraicSorter::begin()
{}

inline AlgebraicSorter::iterator AlgebraicSorter::end()
{}

struct AlgebraicSorter::ModifyingSetCompare
{
    bool operator()(Algebraic * const lhs, Algebraic * const rhs)
    {
        assert(lhs->Invariants());
        assert(rhs->Invariants());

        if (lhs->polynomial != rhs->polynomial)
            while (true)
            {
                if (lhs->rootinterval.lower() > rhs->rootinterval.upper())
                    return false;
                if (lhs->rootinterval.upper() < rhs->rootinterval.lower())
                    return true;

                lhs->tightenInterval();
                rhs->tightenInterval();
            }
        else
            while (true)
            {
                if (lhs->rootinterval.lower() >= rhs->rootinterval.lower() &&
                    lhs->rootinterval.upper() <= rhs->rootinterval.upper())
                    return false;
                if (lhs->rootinterval.lower() > rhs->rootinterval.upper())
                    return false;
                if (lhs->rootinterval.upper() < rhs->rootinterval.lower())
                    return true;

                lhs->tightenInterval();
            }

        assert(false);
        return true;
    }
}
