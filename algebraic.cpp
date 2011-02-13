/*!
 * \file
 * \author Chris de Pujo
 */

#include "algebraic.hpp"

#include <cassert>

Algebraic::Algebraic(const PolynomialQ & p, const IntervalQ & i)
    : polynomial(p), rootinterval(i)
{
    assert(Invariants());
}

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

    const GiNaC::numeric l = rootinterval.lower(),
                         u = rootinterval.upper(); // invoke copy cons?
    const GiNaC::numeric m = (l + u)/2;
    const GiNaC::numeric sample = polynomial.eval(m);

    if (sample == 0)
        //rootinterval.assign((3*l + u)/4, (l + 3*u)/4);
        rootinterval.assign((l+m)/2, (u+m)/2);
    else
    {
        GiNaC::ex num = sample * polynomial.eval(u);

        if (num <= 0)
            rootinterval.assign(m, u);
        else
            rootinterval.assign(l, m);
    }
/*local v,l,u,m,p,pm,rr;
 v := Var(r[2]);
 l := r[1][1];
 u := r[1][2];
 p := r[2];
 m := (u+l)/2;
 pm := eval(p,v=m);
 if pm = 0 then
    rr := [[(l+m)/2,(u+m)/2],p];
 elif pm * eval(p,v=u) <= 0 then
   rr := [[m,u],p];
 else
   rr := [[l,m],p];
 fi;
 return rr; */

    return *this;
}

/*!
 * \todo Double check that the max degree of an irreducible polynomial over
 *		 the rationals is 2 (its not, duh).
 */
GiNaC::numeric Algebraic::Approximate() const
{
    assert(Invariants());

    GiNaC::numeric value = GiNaC::fsolve(polynomial.getEx(),
                                         polynomial.getVariable(),
                                         rootinterval.lower(),
                                         rootinterval.upper());

    return value;
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

Algebraic Algebraic::MakeRational(const GiNaC::numeric & a)
{
    Algebraic alpha(PolynomialQ::GetVar() - a, IntervalQ(a, a));

    return alpha;
}

/*{
    Algebraic alpha(PolynomialQ::GetVar() - a,
                    IntervalQ(a - Algebraic::delta, a + Algebraic::delta));

    return alpha;
}*/

const GiNaC::numeric Algebraic::delta = GiNaC::numeric(1, 64);

std::ostream & operator<<(std::ostream & output, const Algebraic & alpha)
{
    output << "( [" << GiNaC::ex(alpha.rootinterval.lower()) << ", " << GiNaC::ex(alpha.rootinterval.upper())
    	   << "], " << alpha.polynomial.getEx() << " )";

    return output;
}
