/*!
 * \file
 * \author Chris de Pujo
 */

#include "algebraic.hpp"

#include <cassert>
#include <sstream>

//const GiNaC::numeric Algebraic::delta = GiNaC::numeric(1, 64);

/*!
 * \param B Any Algebraic number, including *this.
 * \return -1 if *this in less than B.
 *          0 if *this is equal to B.
 *          1 if *this is greater than B.
 */
int Algebraic::compare(const Algebraic & B) const
{
    assert(Invariants());

    Algebraic a(*this), b(B);

    if (a.polynomial != b.polynomial) // a.compare(b) != 0, or !a.equalTo(b)
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

            a.tightenInterval();
        }
    }

    assert(false);

    return -41;
}

Algebraic & Algebraic::tightenInterval()
{
    assert(Invariants());

    const GiNaC::numeric & l = rootinterval.lower(),
                         & u = rootinterval.upper(); // invoke copy cons?
    const GiNaC::numeric m = (l + u)/2;
    const GiNaC::numeric sample = polynomial.eval(m);

    if (sample.is_zero())
        rootinterval.assign((l+m)/2, (u+m)/2);
    else
    {
        const GiNaC::numeric num = sample * polynomial.eval(u);

        if (!num.is_positive()) // (num <= 0)
            rootinterval.assignLower(m);
        else
            rootinterval.assignUpper(m);
    }

    return *this;
}

/*!
 * \todo Can a non-monic polynomial be irreducible?
 */
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

    if (!polynomial.isIrreducible())
    	return false;

    if (polynomial.isConstant())
    	return false;

    if (!polynomial.isMonic())
    	return false;

    return true;
}

std::string Algebraic::getString() const
{
    std::stringstream s;
    s << *this;
    return s.str();
}
