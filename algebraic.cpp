/*!
 * \file
 * \author Chris de Pujo
 */

#include "algebraic.hpp"

#include <cassert>

inline GiNaC::numeric Algebraic::lower() const
{
    assert(Invariant());

    return rootinterval.lower();
}

inline GiNaC::numeric Algebraic::upper() const
{
    assert(Invariant());

    return rootinterval.upper();
}

/*!
 * \param B Any Algebraic number, including *this.
 * \return -1 if *this in less than B.
 *          0 if *this is equal to B.
 *          1 if *this is greater than B.
 */
int Algebraic::Compare(const Algebraic & B) const
{
    assert(Invariant());

    Algebraic a(*this), b(B);

    if (a.polynomial != b.polynomial) // a.Compare(b) != 0, or !a.equalTo(b)
    {
        while (true)
        {
            if (a.lower() > b.upper())
                return 1;

            if (a.upper() < b.lower())
                return -1;

            a.tightenInterval();
            b.tightenInterval();
        }
    }
    else
    {
        while (true)
        {
            if (a.lower() >= b.lower() && a.upper() <= b.upper())
                return 0;

            if (a.lower() > b.upper())
                return 1;

            if (a.upper() < b.lower())
                return -1;
        }
    }

    assert(false);

    return -41;
}

Algebraic & Algebraic::tightenInterval()
{
    assert(Invariant());

    const GiNaC::numeric l = rootinterval.lower(),
                         u = rootinterval.upper(); // invoke copy cons?
    const GiNaC::numeric m = (l + u)/2;
    const GiNaC::numeric sample = polynomial.eval(m);

    if (sample == 0)
        rootinterval.assign((3*l + u)/4, (l + 3*u)/4);
    else
    {
        GiNaC::ex num = sample * polynomial.eval(rootinterval.upper());

        if (num <= 0)
            rootinterval.assign(m, u);
        else
            rootinterval.assign(l, m);
    }

    return *this;
}

/*!
 * \detail Uses Newton's method to find a floating point approximation of the
 *         algebraic number.
 */
float Algebraic::Approximate() const
{
    assert(Invariant());

    //GiNaC::ex der = polynomial.diff();
    //GiNaC::ex x0 = (rootinterval.lower() + rootinterval.upper()) / 2; // rootinterval.median();

    return 0; // temp, remove later
}

/*!
 * \todo Add check for polynomial irreducibility.
 */
bool Algebraic::Invariant() const
{
    if (this == NULL)
        return false;

    // Interval checks ////////////////////////////////////////////////////////

    if (rootinterval.lower() > rootinterval.upper())
        return false;

    if (!rootinterval.lower().is_rational())
        return false;

    if (!rootinterval.upper().is_rational())
        return false;

    // Polynomial checks //////////////////////////////////////////////////////

    if (!polynomial.Invariant())
        return false;

    return true;
}

/*!
 * \detail Modifies the parameters by resizing the intervals until they do not
 *         intersect.
 */
void SeparateIntervals(Algebraic & a, Algebraic & b)
{
    assert(a.Invariant());
    assert(b.Invariant());

    while ( !(a.upper() < b.lower() || b.upper() < a.lower()) )
    {
        a.tightenInterval();
        b.tightenInterval();
    }
}
