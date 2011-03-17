/*!
 * \file
 * \author Chris de Pujo
 */

#ifndef __TEMPLATEHELP__
#define __TEMPLATEHELP__

#include <algorithm>

class Algebraic;

namespace tmp {

template <typename E>
inline void ReserveHelper(std::vector<E> & container, unsigned int n)
{
    container.reserve(n);
}

template <typename C>
inline void ReserveHelper(C & container, unsigned int n)
{}

template <typename T1, typename T2>
inline void PushBack(T1 & container, const T2 & element)
{
    container.insert(container.end(), element);
}

template <typename container>
inline void Sort(container & group)
{
    std::sort(group.begin(), group.end());
}

template <typename container>
inline void Unique(container & group)
{
    typename container::iterator newEnd = std::unique(group.begin(), group.end());
    group.resize(newEnd - group.begin());
}

template <typename container>
inline void SortUnique(container & group)
{
    Sort(group);
    Unique(group);
}

template <typename container>
inline bool IsSorted(const container & group)
{
    if (group.size() <= 1)
        return true;

    typename container::const_iterator
        cur = group.begin(), next = group.begin(), end = group.end();
    ++next;

    while (next != end)
    {
        if (*cur > *next)
            return false;

        ++cur;
        ++next;
    }

    return true;
}

template <typename container>
inline bool IsUnique(const container & group)
{
    if (group.size() <= 1)
        return true;

    typename container::const_iterator i = group.begin(), j, end = group.end();

    do {
        j = i;

        while (++j != end)
            if (*i == *j)
                return false;
    } while (++i != end);

    return true;
}

template <typename container>
inline bool IsSortedUnique(const container & group)
{
    return (IsSorted(group) && IsUnique(group));
}

template <typename T>
struct DerefLess : public std::binary_function<const T *, const T *, bool>
{
    inline bool operator()(const T * lhs, const T * rhs)
    {
        return (*lhs < *rhs);
    }
};

}

#endif // __TEMPLATEHELP__
