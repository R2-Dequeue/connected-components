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

#include <boost/numeric/interval.hpp>
typedef boost::numeric::interval<GiNaC::numeric> IntervalQ;

class Algebraic; // To allow a circular dependancy.

/*!
 * \brief A class representing univariate polynomials over the rationals.
 *
 * \detail The rationals are represented by bignum rationals and the polynomials
 *         are only limited by memory.
 * \todo Add higher derivative function.
 * \todo Move test methods into proper test class and make it a friend.
 */
class PolynomialQ
{
private:

    static const GiNaC::symbol & variable;  //!< The variable of the polynomial.
    GiNaC::ex polynomial;                   //!< The internal representation of the polynomial.

public:

	PolynomialQ() : polynomial(0) {} //!< The default constructor.
	PolynomialQ(const std::string & s); //!< Parse a polynomial from a string.
	PolynomialQ(const char * const a);
	PolynomialQ(const GiNaC::ex & e);
	PolynomialQ(const GiNaC::numeric & n);

	inline int degree() const; //!< The degree of the polynomial.
	inline bool isMonic() const; //!< True iff the leading coefficient is 1.
	inline bool isZero() const; //!< True iff the polynomial is '0'.
	bool isIrreducible() const;
	inline bool isConstant() const;

	inline GiNaC::symbol getVariable() const { return variable; }
	inline GiNaC::ex getEx() const;

    GiNaC::numeric getCoeff(const unsigned int i) const;

    inline PolynomialQ getMonic() const;
	inline PolynomialQ & makeMonic();

    inline PolynomialQ getDerivative() const;
    inline PolynomialQ & differentiate();

    std::vector<PolynomialQ> getIrreducibleFactors() const;

	//! Returns the sign of f(a).
    int signAt(const Algebraic & a) const;
	GiNaC::numeric eval(const GiNaC::numeric & value) const;
	//!< Returns the value of the polynomial at 'value'.

	static std::vector<PolynomialQ>
        IrreducibleFactors(const std::vector<PolynomialQ> & F);

    static GiNaC::numeric
        Resultant(const PolynomialQ & f, const PolynomialQ & g);

    static std::vector<Algebraic> FindRoots(const std::vector<PolynomialQ> P);

	PolynomialQ & operator+=(const PolynomialQ & rhs);
    PolynomialQ & operator-=(const PolynomialQ & rhs);
	PolynomialQ & operator*=(const PolynomialQ & rhs);
    PolynomialQ & operator/=(const PolynomialQ & rhs);
    PolynomialQ & operator%=(const PolynomialQ & rhs);

    static void TestClass();

    /*!
     * \brief Helper method for internal 'assert' checks.
     * \detail This method is public but shouldn't really be published.
     */
	bool Invariants() const;

    friend PolynomialQ operator+(const PolynomialQ & lhs, const PolynomialQ & rhs);
    friend PolynomialQ operator-(const PolynomialQ & lhs, const PolynomialQ & rhs);
    friend PolynomialQ operator*(const PolynomialQ & lhs, const PolynomialQ & rhs);
    friend PolynomialQ operator/(const PolynomialQ & lhs, const PolynomialQ & rhs);
    friend PolynomialQ operator%(const PolynomialQ & lhs, const PolynomialQ & rhs);

    friend PolynomialQ operator*(const PolynomialQ & lhs, const GiNaC::numeric & num);
    friend PolynomialQ operator/(const PolynomialQ & lhs, const GiNaC::numeric & num);

    friend inline bool operator==(const PolynomialQ & lhs, const PolynomialQ & rhs) {
        assert(lhs.Invariants() && rhs.Invariants());
        return (lhs.polynomial == rhs.polynomial); }
    friend inline bool operator!=(const PolynomialQ & lhs, const PolynomialQ & rhs) {
        assert(lhs.Invariants() && rhs.Invariants());
        return (lhs.polynomial != rhs.polynomial); }

    friend class PolynomialQQ;

protected:

    IntervalQ boundRange(const IntervalQ & interval) const;
    static void TestCompare(const PolynomialQ p,
                            const GiNaC::ex expected,
                            unsigned int & count);
    static inline PolynomialQ ParseString(const std::string & s) const;
};

std::ostream & operator<<(std::ostream & output, const PolynomialQ & p);

//PolynomialQ SubResultant(const PolynomialQ & f,
//                         const PolynomialQ & g,
//                         const int k,
//                         const GiNaC::symbol & x);

#endif // __POLYNOMIALQ__
