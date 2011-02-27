/*!
 * \file
 * \author Chris de Pujo
 */

#ifndef __INTERVALQ__
#define __INTERVALQ__

#include <algorithm>
#include <functional>

#include <boost/array.hpp>

#include <ginac/ginac.h>

/*!
 * \brief A class representing an interval over the rational numbers.
 * \detail The class attempts to be syntax-compatible with the Boost library
 *         interval class (this class was quickly created because I was
 *         getting some odd results).
 */
class IntervalQ
{
private:

    GiNaC::numeric l;
    GiNaC::numeric u;

public:

    IntervalQ();
    IntervalQ(const GiNaC::numeric & a);
    IntervalQ(const GiNaC::numeric & a, const GiNaC::numeric & b);

    const GiNaC::numeric & lower() const;
    const GiNaC::numeric & upper() const;

    GiNaC::numeric median() const;

    void assign(const GiNaC::numeric & a, const GiNaC::numeric & b);
    void assignLower(const GiNaC::numeric & a);
    void assignUpper(const GiNaC::numeric & b);

    IntervalQ & operator+=(const IntervalQ & rhs);
    IntervalQ & operator*=(const IntervalQ & rhs);

    bool Invariants() const;
};

inline IntervalQ::IntervalQ() : l(0), u(0)
{}

inline IntervalQ::IntervalQ(const GiNaC::numeric & a) : l(a), u(a)
{
    assert(Invariants());
}

inline IntervalQ::IntervalQ(const GiNaC::numeric & a, const GiNaC::numeric & b)
    : l(a), u(b)
{
    assert(Invariants());
}

inline const GiNaC::numeric & IntervalQ::lower() const
{
    return l;
}

inline const GiNaC::numeric & IntervalQ::upper() const
{
    return u;
}

inline GiNaC::numeric IntervalQ::median() const
{
    return ( (l+u)/2 );
}

inline void IntervalQ::assign(const GiNaC::numeric & a,
                              const GiNaC::numeric & b)
{
    l = a;
    u = b;

    assert(Invariants());
}

inline void IntervalQ::assignLower(const GiNaC::numeric & a)
{
    l = a;

    assert(Invariants());
}

inline void IntervalQ::assignUpper(const GiNaC::numeric & b)
{
    u = b;

    assert(Invariants());
}

/*!
 * \detail Works correctly in cases such as: I += I.
 */
inline IntervalQ & IntervalQ::operator+=(const IntervalQ & rhs)
{
    assert(Invariants() && rhs.Invariants());

    l += rhs.l;
    u += rhs.u;

    return *this;
}

template <typename T>
struct DerefLess : public std::binary_function<const T *, const T *, bool>
{
    inline bool operator()(const T * lhs, const T * rhs)
    {
        return (*lhs < *rhs);
    }
};

/*!
 * \detail Works correctly in cases such as: I *= I.
 */
inline IntervalQ & IntervalQ::operator*=(const IntervalQ & rhs)
{
    assert(Invariants() && rhs.Invariants());

    boost::array<GiNaC::numeric, 4> a;
    a[0] = l*rhs.l;
    a[1] = l*rhs.u;
    a[2] = u*rhs.l;
    a[3] = u*rhs.u;
    boost::array<GiNaC::numeric *, 4> b;
    b[0] = &(a[0]);
    b[1] = &(a[1]);
    b[2] = &(a[2]);
    b[3] = &(a[3]);

    std::sort(b.begin(), b.end(), DerefLess<GiNaC::numeric>());
    l = *(b[0]);
    u = *(b[3]);

    assert(Invariants());

    return *this;
}

inline bool IntervalQ::Invariants() const
{
    if (!l.is_rational())
        return false;

    if (!u.is_rational())
        return false;

    if (l > u)
        return false;

    return true;
}

#endif // __INTERVALQ__
