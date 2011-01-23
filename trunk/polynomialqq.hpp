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
    PolynomialQQ(const char * const s);
    PolynomialQQ(const GiNaC::ex & e);
	PolynomialQQ(const GiNaC::numeric & n);

    PolynomialQQ getDerivative(unsigned int variable) const;
    PolynomialQQ & differentiate(unsigned int variable);

    std::vector<PolynomialQQ> getIrreducibleFactors() const;

    int signAt(const Algebraic & a, const Algebraic & b) const;

    GiNaC::ex getEx() const;

    static std::vector<PolynomialQQ>
        IrreducibleFactors(const std::vector<PolynomialQQ> & F);

    static PolynomialQ Resultant(const PolynomialQQ & f,
                                 const PolynomialQQ & g,
                                 unsigned int var);

	PolynomialQQ & operator+=(const PolynomialQQ & rhs);
    PolynomialQQ & operator-=(const PolynomialQQ & rhs);
	PolynomialQQ & operator*=(const PolynomialQQ & rhs);
    //PolynomialQQ & operator/=(const PolynomialQQ & rhs);
    //PolynomialQQ & operator%=(const PolynomialQQ & rhs);

    /*!
     * \brief Helper method for internal 'assert' checks.
     * \detail This method is public but shouldn't really be published.
     */
    bool Invariant() const;

    // These are the method versions of the binary arithmetic operators.
	PolynomialQQ operator+(const PolynomialQQ & rhs) const;
    PolynomialQQ operator-(const PolynomialQQ & rhs) const;
	PolynomialQQ operator*(const PolynomialQQ & rhs) const;
    //PolynomialQQ operator/(const PolynomialQQ & rhs) const;
    //PolynomialQQ operator%(const PolynomialQQ & rhs) const;

    inline bool operator==(const PolynomialQQ & rhs) const;
    inline bool operator!=(const PolynomialQQ & rhs) const;
};

#endif // __POLYNOMIALQQ__
