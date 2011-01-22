/*!
 * \file
 * \author Chris de Pujo
 */

#ifndef __POLYNOMIALQQ__
#define __POLYNOMIALQQ__

#include <ginac/ginac.h>

#include "polynomialq.hpp"

#include <string>
#include <cassert>
#include <vector>

class PolynomialQQ
{
private:

    static const GiNaC::symbol & var1;
    static const GiNaC::symbol & var2;
    GiNaC::ex polynomial;

public:

    PolynomialQQ() : polynomial(0) {} //!< The default constructor.
    PolynomialQQ(const std::string & s);
    PolynomialQQ(const GiNaC::ex & e);
	PolynomialQQ(const GiNaC::numeric & n);

	PolynomialQQ & operator+=(const PolynomialQQ & rhs);
    PolynomialQQ & operator-=(const PolynomialQQ & rhs);
	PolynomialQQ & operator*=(const PolynomialQQ & rhs);
    //PolynomialQQ & operator/=(const PolynomialQQ & rhs);
    //PolynomialQQ & operator%=(const PolynomialQQ & rhs);

    // These are the method versions of the binary arithmetic operators.
	PolynomialQQ operator+(const PolynomialQQ & rhs) const; //!< Normal binary addition.
    PolynomialQQ operator-(const PolynomialQQ & rhs) const;
	PolynomialQQ operator*(const PolynomialQQ & rhs) const;
    //PolynomialQQ operator/(const PolynomialQQ & rhs) const;
    //PolynomialQQ operator%(const PolynomialQQ & rhs) const;

    inline bool operator==(const PolynomialQQ & rhs) const;
    inline bool operator!=(const PolynomialQQ & rhs) const;

    PolynomialQQ getDerivative(unsigned int variable) const;
    PolynomialQQ & differentiate(unsigned int variable);

    int signAt(const Algebraic & a, const Algebraic & b) const;

    /*!
     * \brief Helper function for internal 'assert' checks.
     * \detail This member is public but shouldn't really be published.
     */
    bool Invariant() const;
};

std::vector<PolynomialQQ> IrreducibleFactors(const std::vector<PolynomialQQ> & F);
PolynomialQ Resultant(const PolynomialQQ & f, const PolynomialQQ & g, unsigned int var);

#endif // __POLYNOMIALQQ__
