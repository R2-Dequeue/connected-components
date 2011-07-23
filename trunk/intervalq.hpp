#ifndef __INTERVALQ__
#define __INTERVALQ__

#include <cassert>
#include <algorithm>
#include <functional>

#include <boost/array.hpp>

#include <ginac/ginac.h>

/*!
 * \brief A class representing an interval over the rational numbers.
 * \details The class attempts to be syntax-compatible with the Boost library
 *          interval class (this class was quickly created because I was
 *          getting some odd results).
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
 * \details Works correctly in cases such as: I += I.
 */
inline IntervalQ & IntervalQ::operator+=(const IntervalQ & rhs)
{
    assert(Invariants() && rhs.Invariants());

    l += rhs.l;
    u += rhs.u;

    return *this;
}

/*!
 * \details Works correctly in cases such as: I *= I.
 * \todo Switch to boost simultaneous min/max.
 */
inline IntervalQ & IntervalQ::operator*=(const IntervalQ & rhs)
{
    assert(Invariants() && rhs.Invariants());

    boost::array<GiNaC::numeric, 4> a;
    a[0] = this->l*rhs.l;
    a[1] = this->l*rhs.u;
    a[2] = this->u*rhs.l;
    a[3] = this->u*rhs.u;

    l = *std::min_element(a.begin(), a.end());
    u = *std::max_element(a.begin(), a.end());

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
