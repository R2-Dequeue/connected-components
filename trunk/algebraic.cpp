/*!
 * \file
 * \author Chris de Pujo
 */

#include "algebraic.hpp"

#include <cassert>

inline GiNaC::numeric Algebraic::lower() const
{
    assert(Invariants());

    return rootinterval.lower();
}

inline GiNaC::numeric Algebraic::upper() const
{
    assert(Invariants());

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
    assert(Invariants());

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
    assert(Invariants());

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
 * \todo Double check that the max degree of an irreducible polynomial over
 *		 the rationals is 2.
 */
double Algebraic::Approximate() const
{
    assert(Invariants());
}

bool Algebraic::Invariants() const
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

    if (!polynomial.Invariants())
        return false;
    
    // Can a non-monic polynomial be irreducible?
    
    if (!polynomial.isIrreducible())
    	return false;
    
    if (polynomial.isConstant())
    	return false;
    
    if (!polynomial.isMonic())
    	return false;

    return true;
}

/*!
 * \detail Modifies the parameters by resizing the intervals until they do not
 *         intersect.
 */
void Algebraic::SeparateIntervals(Algebraic & a, Algebraic & b)
{
    assert(a.Invariants());
    assert(b.Invariants());

    while ( !(a.upper() < b.lower() || b.upper() < a.lower()) )
    {
        a.tightenInterval();
        b.tightenInterval();
    }
}
