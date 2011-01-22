/*!
 * \file
 * \author Chris de Pujo
 */

#ifndef __POLYNOMIALQ__
#define __POLYNOMIALQ__

#include <cassert>
#include <string>
#include <vector>

#include <ginac/ginac.h>

class Algebraic; // To allow a circular dependancy.
#include <boost/numeric/interval.hpp>
typedef boost::numeric::interval<GiNaC::numeric> IntervalQ;

/*!
 * \brief A class representing univariate polynomials over the rationals.
 *
 * \detail The rationals are represented by bignum rationals and the polynomials
 *         are only limited by memory.
 */
class PolynomialQ
{
private:

    static const GiNaC::symbol & variable;  //!< The variable of the polynomial.
    GiNaC::ex polynomial;                   //!< The internal representation of the polynomial.

public:

	PolynomialQ() : polynomial(0) {} //!< The default constructor.
	PolynomialQ(const std::string & s); //!< Parse a polynomial from a string.
	PolynomialQ(const char * const s);
	PolynomialQ(const GiNaC::ex & e);
	PolynomialQ(const GiNaC::numeric & n);

	inline int degree() const; //!< The degree of the polynomial.
	inline bool isMonic() const; //!< True iff the leading coefficient is zero.
	inline bool isZero() const; //!< True iff the polynomial is '0'.

    PolynomialQ getMonic();
	PolynomialQ & makeMonic();

    PolynomialQ getDerivative() const;
    PolynomialQ & differentiate();

    int signAt(const Algebraic & a) const;
	GiNaC::numeric eval(const GiNaC::numeric & value) const;
	//!< Returns the value of the polynomial at 'value'.

	GiNaC::ex getEx() const;

	PolynomialQ & operator+=(const PolynomialQ & rhs);
    PolynomialQ & operator-=(const PolynomialQ & rhs);
	PolynomialQ & operator*=(const PolynomialQ & rhs);
    PolynomialQ & operator/=(const PolynomialQ & rhs);
    PolynomialQ & operator%=(const PolynomialQ & rhs);

    static void TestClass();

    /*!
     * \brief Helper function for internal 'assert' checks.
     * \detail This member is public but shouldn't really be published.
     */
	bool Invariant() const;

    friend PolynomialQ operator+(const PolynomialQ & lhs, const PolynomialQ & rhs);
    friend PolynomialQ operator-(const PolynomialQ & lhs, const PolynomialQ & rhs);
    friend PolynomialQ operator*(const PolynomialQ & lhs, const PolynomialQ & rhs);
    friend PolynomialQ operator/(const PolynomialQ & lhs, const PolynomialQ & rhs);
    friend PolynomialQ operator%(const PolynomialQ & lhs, const PolynomialQ & rhs);

    friend PolynomialQ operator*(const PolynomialQ & lhs, const GiNaC::numeric & num);
    friend PolynomialQ operator/(const PolynomialQ & lhs, const GiNaC::numeric & num);

    friend inline bool operator==(const PolynomialQ & lhs, const PolynomialQ & rhs);
    friend inline bool operator!=(const PolynomialQ & lhs, const PolynomialQ & rhs);

protected:

    IntervalQ boundRange(const IntervalQ & interval) const;
    static void TestCompare(const PolynomialQ p,
                            const GiNaC::ex expected,
                            unsigned int & count);
};

std::ostream & operator<<(std::ostream & output, const PolynomialQ & p);

//PolynomialQ SubResultant(const PolynomialQ & f,
//                         const PolynomialQ & g,
//                         const int k,
//                         const GiNaC::symbol & x);
//std::vector<PolynomialQ> IrreducibleFactors(const std::vector<PolynomialQ> & F);
//std::list<Algebraic> FindRoots(std::list<PolynomialQ> P);
//PolynomialQ Resultant(const PolynomialQ & f, const PolynomialQ & g);

// These are the method versions of the binary arithmetic operators.
PolynomialQ operator+(const PolynomialQ & lhs, const PolynomialQ & rhs);
PolynomialQ operator-(const PolynomialQ & lhs, const PolynomialQ & rhs);
PolynomialQ operator*(const PolynomialQ & lhs, const PolynomialQ & rhs);
PolynomialQ operator/(const PolynomialQ & lhs, const PolynomialQ & rhs);
PolynomialQ operator%(const PolynomialQ & lhs, const PolynomialQ & rhs);

PolynomialQ operator*(const PolynomialQ & lhs, const GiNaC::numeric & num);
PolynomialQ operator/(const PolynomialQ & lhs, const GiNaC::numeric & num);

inline bool operator==(const PolynomialQ & lhs, const PolynomialQ & rhs);
inline bool operator!=(const PolynomialQ & lhs, const PolynomialQ & rhs);

#endif // __POLYNOMIALQ__
