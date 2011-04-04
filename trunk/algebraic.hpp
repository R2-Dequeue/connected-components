/*!
 * \file
 * \author Chris de Pujo
 */

#ifndef __ALGEBRAIC__
#define __ALGEBRAIC__

#include <cassert>
#include <vector>
#include <set>
#include <list>

#include <ginac/ginac.h>

#include "polynomialq.hpp"

/*!
 * \brief A class representing an algebraic number.
 *
 * \detail Internally, the number is represented by a polynomial and a closed
 *         interval on which the polynomial has one and only one root, which is
 *         the number. Note that the interval can be be a single point
 *         (i.e. [a, a]). Internally, the polynomial is represented by a
 *         GiNaC::ex object which has copy-on-write semantics.
 */
class Algebraic
{
private:

	/*!
     * \brief The polynomial that has as one of its roots this number.
     * \detail The polynomial is always monic, irreducible, and non-constant
     *		   with rational coefficients.
     */
    PolynomialQ polynomial;
    //! The interval that contains this number with rational endpoints.
    mutable IntervalQ rootinterval;

    struct ModifyingSetCompare;

public:

    typedef std::vector<Algebraic>                  vector;
    typedef std::set<Algebraic>                     set;
    typedef std::set<Algebraic,ModifyingSetCompare> modifying_set;
    typedef std::list<Algebraic>                    list;

    //! The default constructor.
    Algebraic();
    //! Basic constructor to assemble an algebraic number.
    Algebraic(const PolynomialQ & p, const IntervalQ & i);

    // can I make a constructor that will only take STATIC ints or longs?
    // (for initialization).

    //! Returns the lower bound of the interval.
    const GiNaC::numeric & lower() const;
    //!< Returns the upper bound of the interval.
    const GiNaC::numeric & upper() const;

    bool isZero() const;
    bool isNonZero() const;
    bool isPositive() const;
    bool isNegative() const;
    bool isNonNegative() const;
    bool isNonPositive() const;
    bool isRational() const;
    bool isInteger() const;
    bool isNatural() const;     //!< true iff alpha is in {1,2,3,4,...}.

    const IntervalQ & getInterval() const;
    const GiNaC::ex & getEx() const;
    const PolynomialQ & getPolynomial() const;
    Algebraic getAbs() const; //!< Return the absolute value of alpha.

    int sgn() const;
    int compare(const Algebraic & B) const; //!< Compare this number with another.
    void tightenInterval() const; //!< Shrinks the interval to a proper subset.
    Algebraic & takeAbs();
    Algebraic & mul(const GiNaC::numeric & a);
    GiNaC::numeric Approximate() const; //!< A floating point approximation of this number.
    int roundToInt();
    GiNaC::numeric approx(unsigned int d) const;

    std::string getString() const;

    //! Shrinks the internal intervals of a & b so that they don't intersect.
    static void SeparateIntervals(const Algebraic & a, const Algebraic & b);
    static void SeparateIntervals(Algebraic & a, Algebraic & b,
                                  const GiNaC::numeric & delta);
    template <typename Container>
    static void SeparateIntervals(Container & alphas);
    template <typename Container>
    static void SeparateIntervals(Container & alphas,
                                  const GiNaC::numeric & delta);
    //! Returns an object with a degree 1 polynomial and interval [a, a].
    static Algebraic MakeRational(const GiNaC::numeric & a);
    //! Returns an object with a degree 1 polynomial and interval [a-delta, a+delta].
    static Algebraic MakeWideRational(const GiNaC::numeric & a);
    //! Sets the radius for MakeWideRational.
    //static const GiNaC::numeric delta;

    /*!
     * \brief Helper method for internal 'assert' checks.
     * \detail This method is public but shouldn't really be published.
     */
	bool Invariants() const;

	friend class PolynomialQ;
	friend struct ModifyingSetCompare;
	//friend class AlgebraicSorter;
};

////////////////////////////////////////////////////////////////////////////////

inline bool operator==(const Algebraic & alpha, const Algebraic & beta);
inline bool operator!=(const Algebraic & alpha, const Algebraic & beta);
inline bool operator<(const Algebraic & alpha, const Algebraic & beta);
inline bool operator<=(const Algebraic & alpha, const Algebraic & beta);
inline bool operator>(const Algebraic & alpha, const Algebraic & beta);
inline bool operator>=(const Algebraic & alpha, const Algebraic & beta);
inline std::ostream & operator<<(std::ostream & output, const Algebraic & alpha);

/*!
 * \detail Used for comparison in sets containing roots.  Provides a speedup
 *         by not repeating \p tightenInterval calls.
 */
struct Algebraic::ModifyingSetCompare
{
bool operator()(const Algebraic & lhs, const Algebraic & rhs)
{
    Algebraic & l = const_cast<Algebraic &>(lhs);
    Algebraic & r = const_cast<Algebraic &>(rhs);

    assert(l.Invariants());
    assert(r.Invariants());

    if (l.polynomial != r.polynomial)
        while (true)
        {
            if (l.rootinterval.lower() > r.rootinterval.upper())
                return false;
            if (l.rootinterval.upper() < r.rootinterval.lower())
                return true;

            l.tightenInterval();
            r.tightenInterval();
        }
    else
        while (true)
        {
            if (l.rootinterval.lower() >= r.rootinterval.lower() &&
                l.rootinterval.upper() <= r.rootinterval.upper())
                return false;
            if (l.rootinterval.lower() > r.rootinterval.upper())
                return false;
            if (l.rootinterval.upper() < r.rootinterval.lower())
                return true;

            l.tightenInterval();
        }

    assert(false);
    return true;
}
};

#define __DELTA ( GiNaC::numeric(1, 64) )

#include "polynomialq.hpp"

inline Algebraic::Algebraic()
    : polynomial(GiNaC::ex(PolynomialQ::GetVar())),
      rootinterval(-__DELTA, __DELTA)
{
    assert(Invariants());
}

inline Algebraic::Algebraic(const PolynomialQ & p, const IntervalQ & i)
    : polynomial(p), rootinterval(i)
{
    assert(Invariants());
}

inline const GiNaC::numeric & Algebraic::lower() const
{
    assert(Invariants());

    return rootinterval.lower();
}

inline const GiNaC::numeric & Algebraic::upper() const
{
    assert(Invariants());

    return rootinterval.upper();
}

inline bool Algebraic::isZero() const
{
    assert(Invariants());

    return (polynomial.getEx() == GiNaC::ex(PolynomialQ::GetVar()));
}

inline bool Algebraic::isNonZero() const
{
    return (!this->isZero());
}

inline bool Algebraic::isPositive() const
{
    return (this->sgn() == 1);
}

inline bool Algebraic::isNegative() const
{
    return (this->sgn() == -1);
}

inline bool Algebraic::isNonNegative() const
{
    return (this->sgn() >= 0);
}

inline bool Algebraic::isNonPositive() const
{
    return (this->sgn() <= 0);
}

inline bool Algebraic::isRational() const
{
    assert(Invariants());

    return (polynomial.degree() == 1);
}

inline bool Algebraic::isInteger() const
{
    assert(Invariants());

    return (polynomial.degree() == 1 &&
            polynomial.getCoeff(0).is_integer());
}

inline bool Algebraic::isNatural() const
{
    assert(Invariants());

    return (polynomial.degree() == 1 &&
            polynomial.getCoeff(0).is_integer() &&
            polynomial.getCoeff(0).is_negative());// x+a=0 => x=-a
}

inline const IntervalQ & Algebraic::getInterval() const
{
    assert(Invariants());

    return rootinterval;
}

inline const GiNaC::ex & Algebraic::getEx() const
{
    assert(Invariants());

    return polynomial.getEx();
}

inline const PolynomialQ & Algebraic::getPolynomial() const
{
    assert(Invariants());

    return polynomial;
}

inline Algebraic Algebraic::getAbs() const
{
    assert(Invariants());

    Algebraic alpha(*this);

    while (alpha.lower().csgn() != alpha.upper().csgn())
        alpha.tightenInterval();

    if (alpha.lower().csgn() >= 0)
        return *this;

    alpha.polynomial.negate();
    const GiNaC::numeric a = this->lower(), b = this->upper();
    alpha.rootinterval.assign(-b, -a);

    assert(alpha.Invariants());

    return alpha;
}

inline int Algebraic::sgn() const
{
    assert(Invariants());

    if (this->isZero())
        return 0;

    Algebraic alpha(*this);

    while (alpha.lower().csgn() != alpha.upper().csgn())
        alpha.tightenInterval();

    return alpha.lower().csgn();
}

/*!
 * \todo Fix polynomial multiplication.
 */
inline Algebraic & Algebraic::takeAbs()
{
    if (this->isZero())
        return *this;

    while (this->lower().csgn() != this->upper().csgn())
        this->tightenInterval();

    if (this->lower().csgn() >= 0)
        return *this;

    const GiNaC::numeric a = this->lower(), b = this->upper();
    polynomial.negate();
    rootinterval.assign(-b, -a);

    unsigned int d = polynomial.degree();
    if ((d & 1) != 0)
        //this->polynomial = this->polynomial * GiNaC::numeric(-1);
        polynomial.makeMonic();

    return *this;
}

inline Algebraic & Algebraic::mul(const GiNaC::numeric & a)
{
    assert(a.is_rational());

    if (a.is_zero())
    {
        *this = Algebraic();
        return *this;
    }

    GiNaC::numeric l = this->lower(), u = this->upper();

    if (a.is_positive())
        rootinterval.assign(a*l, a*u);
    else
        rootinterval.assign(a*u, a*l);

    polynomial.subs(a.inverse()*PolynomialQ::GetVar());

    assert(lower() <= upper());
    assert(lower().is_rational());
    assert(upper().is_rational());

    return *this;
}

/*!
 * \todo Double check that the max degree of an irreducible polynomial over
 *		 the rationals is 2 (its not, duh).
 */
inline GiNaC::numeric Algebraic::Approximate() const
{
    assert(Invariants());

    Algebraic temp(*this);

    while (temp.upper() - temp.lower() > GiNaC::numeric(1,10000))
        temp.tightenInterval();
/*
    GiNaC::numeric value = GiNaC::fsolve(polynomial.getEx(),
                                         polynomial.getVariable(),
                                         temp.lower(), //rootinterval.lower(),
                                         temp.upper()); //rootinterval.upper());

    return value;*/
    return (temp.lower() + temp.upper()) / 2;
}

inline int Algebraic::roundToInt()
{
    assert(Invariants());

    //return GiNaC::fsolve(polynomial.getEx(),
    //                     PolynomialQ::GetVar(),
    //                     rootinterval.lower(),
    //                     rootinterval.upper());

    const GiNaC::numeric delta(1, 10);

    while (this->rootinterval.upper() - this->rootinterval.lower() > delta)
        this->tightenInterval();

    const GiNaC::numeric a =
        (this->rootinterval.lower()+this->rootinterval.upper())/2;

    // Round to nearest
    // Round half away from zero for tie-breaking

    assert(a.is_rational());

    const GiNaC::numeric n = a.numer(), d = a.denom();

    assert(n.is_integer());
    assert(d.is_nonneg_integer());
    assert(GiNaC::gcd(n, d) == 1);

    GiNaC::numeric rem;
    const GiNaC::numeric quo = GiNaC::iquo(n, d, rem);

    assert(quo.is_integer());//(quo.is_nonneg_integer());
    assert(rem.is_integer());
    assert(rem.is_zero() || rem.csgn() == n.csgn());

    if (a.is_positive())
    {
        if (rem < d/2)
            return quo.to_int();
        else
            return (quo.to_int()+1);
    }

    if (rem <= (d*(-1))/2)
        return (quo.to_int()*(-1)-1);
    else
        return (quo.to_int()*(-1));
}

inline GiNaC::numeric Algebraic::approx(unsigned int d) const
{
    assert(Invariants());

    Algebraic temp(*this);

    const GiNaC::numeric denom = GiNaC::pow(GiNaC::numeric(10),
                                            GiNaC::numeric(d));
    const GiNaC::numeric delta(denom.inverse());

    while (temp.upper() - temp.lower() > delta)
        temp.tightenInterval();

    return (temp.lower()+temp.upper())/2;
}

/*!
 * \detail Modifies the parameters by resizing the intervals until they do not
 *         intersect.
 */
inline void Algebraic::SeparateIntervals(const Algebraic & a, const Algebraic & b)
{
    assert(a.Invariants());
    assert(b.Invariants());

    if (a != b)
        while ( !(a.upper() < b.lower() || b.upper() < a.lower()) )
        {
            a.tightenInterval();
            b.tightenInterval();
        }
}

inline void Algebraic::SeparateIntervals(Algebraic & a, Algebraic & b,
                                         const GiNaC::numeric & delta)
{
    assert(a.Invariants());
    assert(b.Invariants());

    while (a.upper() - a.lower() > delta)
        a.tightenInterval();
    while (b.upper() - b.lower() > delta)
        b.tightenInterval();

    if (a != b)
    {
        while ( !(a.upper() < b.lower() || b.upper() < a.lower()) )
        {
            a.tightenInterval();
            b.tightenInterval();
        }
    }
}

/*!
 * \param alphas A container of SORTED Algebraic objects.
 */
template <typename Container>
void Algebraic::SeparateIntervals(Container & alphas)
{
    if (alphas.size() <= 1)
        return;

    typename Container::const_iterator alpha = alphas.begin(),
                                 endAlpha = alphas.end(),
                                 nextAlpha = alphas.begin();
    ++nextAlpha;

    while (nextAlpha != endAlpha)
    {
        assert(*alpha < *nextAlpha);

        Algebraic::SeparateIntervals(*alpha, *nextAlpha);

        ++alpha;
        ++nextAlpha;
    }
}

/*!
 * \param alphas A container of SORTED Algebraic objects.
 */
template <typename Container>
void Algebraic::SeparateIntervals(Container & alphas,
                                  const GiNaC::numeric & delta)
{
    if (alphas.size() <= 1)
        return;

    typename Container::iterator alpha = alphas.begin(),
                                 endAlpha = alphas.end(),
                                 nextAlpha = alphas.begin();
    ++nextAlpha;

    while (nextAlpha != endAlpha)
    {
        Algebraic::SeparateIntervals(*alpha, *nextAlpha);

        while (alpha->upper() - alpha->lower() > delta)
            alpha->tightenInterval();

        ++alpha;
        ++nextAlpha;
    }

    while (alpha->upper() - alpha->lower() > delta)
        alpha->tightenInterval();
}

/*!
 * \detail Creates a new \p Algebraic object with polynomial <p> x - a <\p>
 *         and interval [\p a, \p a].
 */
inline Algebraic Algebraic::MakeRational(const GiNaC::numeric & a)
{
    Algebraic alpha(PolynomialQ::GetVar() - a, IntervalQ(a, a));

    return alpha;
}

/*!
 * \detail Creates a new \p Algebraic object with a non-trivial interval
 *         [<p>x - __DELTA, x + __DELTA<\p>].
 */
inline Algebraic Algebraic::MakeWideRational(const GiNaC::numeric & a)
{
    Algebraic alpha(PolynomialQ::GetVar() - a,
                    IntervalQ(a-__DELTA, a+__DELTA));

    return alpha;
}

inline bool operator==(const Algebraic & alpha, const Algebraic & beta)
{
    assert(alpha.Invariants());
    assert(beta.Invariants());

    return (alpha.compare(beta) == 0);
}

inline bool operator!=(const Algebraic & alpha, const Algebraic & beta)
{
    assert(alpha.Invariants());
    assert(beta.Invariants());

    return (alpha.compare(beta) != 0);
}

inline bool operator<(const Algebraic & alpha, const Algebraic & beta)
{
    assert(alpha.Invariants());
    assert(beta.Invariants());

    return (alpha.compare(beta) == -1);
}

inline bool operator<=(const Algebraic & alpha, const Algebraic & beta)
{
    assert(alpha.Invariants());
    assert(beta.Invariants());

    return (alpha.compare(beta) <= 0);
}

inline bool operator>(const Algebraic & alpha, const Algebraic & beta)
{
    assert(alpha.Invariants());
    assert(beta.Invariants());

    return (alpha.compare(beta) == 1);
}

inline bool operator>=(const Algebraic & alpha, const Algebraic & beta)
{
    assert(alpha.Invariants());
    assert(beta.Invariants());

    return (alpha.compare(beta) >= 0);
}

/*!
 * \detail Outputs an exact representation of the number to \p output.
 */
inline std::ostream & operator<<(std::ostream & output, const Algebraic & alpha)
{
    output << "( [" << GiNaC::ex(alpha.lower()) << ", " << GiNaC::ex(alpha.upper())
    	   << "], " << alpha.getEx() << " )";

    return output;
}

#endif // __ALGEBRAIC__
