/*!
 * \file
 * \author Chris de Pujo
 */

#ifndef __ALGEBRAIC__
#define __ALGEBRAIC__

//#include <boost/numeric/interval.hpp>

#include <ginac/ginac.h>

#include "polynomialq.hpp"
//class PolynomialQ; // To allow a circular dependancy.

//! This type represents the interval in class Algebraic.
//typedef boost::numeric::interval<GiNaC::numeric> IntervalQ;

/*!
 * \brief A class representing an algebraic number.
 *
 * \detail Internally, the number is represented by a polynomial and an interval
 *         on which the polynomial has one and only one root, which is the
 *         number.
 */
class Algebraic
{
private:

    PolynomialQ polynomial; //!< Irreducible polynomial.
    IntervalQ rootinterval; //!< The interval that contains this number.

public:

    //! The default constructor.
    Algebraic() : polynomial(GiNaC::numeric(0)), rootinterval(0,0) {};
    //~Algebraic();

    //! Basic constructor to assemble an algebraic number.
    Algebraic(const PolynomialQ & p, const IntervalQ & i)
        : polynomial(p), rootinterval(i) { assert(Invariant()); };

    //Algebraic(const GiNaC::numeric & n)
    //    : polynomial(PolynomialQ::variable - n), interval(n)
    //    { assert(Invariant()); }

    // can I make a constructor that will only take STATIC ints or longs?
    // (for initialization).

    inline GiNaC::numeric lower() const; //!< Returns the lower bound of the interval.
    inline GiNaC::numeric upper() const; //!< Returns the upper bound of the interval.

    inline IntervalQ getInterval() const { return IntervalQ(rootinterval); }
    inline GiNaC::ex getEx() const { return polynomial.getEx(); }
    inline PolynomialQ getPolynomial() const { return polynomial; }

    int Compare(const Algebraic & B) const; //!< Compare this number with another.
    Algebraic & tightenInterval(); //!< Shrinks the interval to a proper subset.
    double Approximate() const; //!< A floating point approximation of this number.

    //! Shrinks the internal intervals of a & b so that they don't intersect.
    static void SeparateIntervals(Algebraic & a, Algebraic & b);

    /*!
     * \brief Helper method for internal 'assert' checks.
     * \detail This method is public but shouldn't really be published.
     */
	bool Invariant() const;

	inline bool operator==(const Algebraic & b) { return (Compare(b) == 0); }
    inline bool operator!=(const Algebraic & b) { return (Compare(b) != 0); }

	friend class PolynomialQ;
};

//Algebraic operator+(const Algebraic & lhs, const Algebraic & rhs);
//Algebraic operator-(const Algebraic & lhs, const Algebraic & rhs);
//Algebraic operator*(const Algebraic & lhs, const Algebraic & rhs);
//Algebraic operator/(const Algebraic & lhs, const Algebraic & rhs);
//Algebraic operator%(const Algebraic & lhs, const Algebraic & rhs);

#endif // __ALGEBRAIC__
