/*!
 * \file
 * \author Chris de Pujo
 */

#ifndef __ALGEBRAIC__
#define __ALGEBRAIC__

#include <ginac/ginac.h>

#include "polynomialq.hpp"

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

	/*!
     * \brief The polynomial that has as one of its roots this number.
     * \detail The polynomial is always monic, irreducible, and non-constant
     *		   with rational coefficients.
     */
    PolynomialQ polynomial;
    //! The interval that contains this number with rational endpoints.
    IntervalQ rootinterval;

public:

    //! The default constructor.
    Algebraic()
        : polynomial(GiNaC::ex(PolynomialQ::GetVar())), rootinterval(-1, 1) {};

    //! Basic constructor to assemble an algebraic number.
    Algebraic(const PolynomialQ & p, const IntervalQ & i);

    // can I make a constructor that will only take STATIC ints or longs?
    // (for initialization).

    GiNaC::numeric lower() const; //!< Returns the lower bound of the interval.
    GiNaC::numeric upper() const; //!< Returns the upper bound of the interval.

    inline IntervalQ getInterval() const { return IntervalQ(rootinterval); }
    inline GiNaC::ex getEx() const { return polynomial.getEx(); }
    inline PolynomialQ getPolynomial() const { return polynomial; }

    int compare(const Algebraic & B) const; //!< Compare this number with another.
    Algebraic & tightenInterval(); //!< Shrinks the interval to a proper subset.
    GiNaC::numeric Approximate() const; //!< A floating point approximation of this number.

    //! Shrinks the internal intervals of a & b so that they don't intersect.
    static void SeparateIntervals(Algebraic & a, Algebraic & b);

    static Algebraic MakeRational(const GiNaC::numeric & a);

    static const GiNaC::numeric delta;

    /*!
     * \brief Helper method for internal 'assert' checks.
     * \detail This method is public but shouldn't really be published.
     */
	bool Invariants() const;

	inline bool operator==(const Algebraic & b) { return (compare(b) == 0); }
    inline bool operator!=(const Algebraic & b) { return (compare(b) != 0); }

	friend class PolynomialQ;
	friend std::ostream & operator<<(std::ostream & output, const Algebraic & alpha);
};

inline bool operator<(const Algebraic & alpha, const Algebraic & beta)
    { return (alpha.compare(beta) == -1); }

//Algebraic operator+(const Algebraic & lhs, const Algebraic & rhs);
//Algebraic operator-(const Algebraic & lhs, const Algebraic & rhs);
//Algebraic operator*(const Algebraic & lhs, const Algebraic & rhs);
//Algebraic operator/(const Algebraic & lhs, const Algebraic & rhs);
//Algebraic operator%(const Algebraic & lhs, const Algebraic & rhs);

#endif // __ALGEBRAIC__
