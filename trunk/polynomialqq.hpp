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

#include <boost/tuple/tuple.hpp>

/*!
 * \todo Add support for different orderings (lots of functions to add...).
 */
class PolynomialQQ
{
private:

    static const GiNaC::symbol & var1;
    static const GiNaC::symbol & var2;
    GiNaC::ex polynomial;

public:

    PolynomialQQ() : polynomial(0) {} //!< The default constructor.
    PolynomialQQ(const std::string & s);
    PolynomialQQ(const char * const a);
    PolynomialQQ(const GiNaC::ex & e);
	PolynomialQQ(const GiNaC::numeric & n);

    int degree() const; //!< The degree of the polynomial.
    bool isMonic() const; //!< True iff the leading coefficient is 1.
    bool isZero() const; //!< True iff the polynomial is '0'.
	bool isIrreducible() const;
    bool isConstant() const;

    PolynomialQQ getDerivative(unsigned int variable) const;
    PolynomialQQ & differentiate(unsigned int variable);

    PolynomialQ subx(const GiNaC::numeric & a) const;
    PolynomialQ suby(const GiNaC::numeric & b) const;

    std::vector<PolynomialQQ> getIrreducibleFactors() const;

    int signAt(const Algebraic & alpha, const Algebraic & beta) const;

    inline GiNaC::ex getEx() const;

    static std::vector<PolynomialQQ>
        IrreducibleFactors(const std::vector<PolynomialQQ> & F);

    static PolynomialQ Resultant(const PolynomialQQ & f,
                                 const PolynomialQQ & g,
                                 unsigned int var);

    //static PolynomialQ Subresultant(const PolynomialQQ & f,
    //                                const PolynomialQQ & g,
    //                                const unsigned int k,
    //                                const unsigned int var);

	PolynomialQQ & operator+=(const PolynomialQQ & rhs);
    PolynomialQQ & operator-=(const PolynomialQQ & rhs);
	PolynomialQQ & operator*=(const PolynomialQQ & rhs);
    //PolynomialQQ & operator/=(const PolynomialQQ & rhs);
    //PolynomialQQ & operator%=(const PolynomialQQ & rhs);

    /*!
     * \brief Helper method for internal 'assert' checks.
     * \detail This method is public but shouldn't really be published.
     */
    bool Invariants() const;

    // These are the method versions of the binary arithmetic operators.
	PolynomialQQ operator+(const PolynomialQQ & rhs) const;
    PolynomialQQ operator-(const PolynomialQQ & rhs) const;
	PolynomialQQ operator*(const PolynomialQQ & rhs) const;
    //PolynomialQQ operator/(const PolynomialQQ & rhs) const;
    //PolynomialQQ operator%(const PolynomialQQ & rhs) const;

    inline bool operator==(const PolynomialQQ & rhs) const;
    inline bool operator!=(const PolynomialQQ & rhs) const;

//protected:

    static Algebraic ANComb(Algebraic alpha,
                            Algebraic beta,
                            const GiNaC::numeric & t);
    static boost::tuple<Algebraic, PolynomialQ, PolynomialQ>
        Simple(const Algebraic & alpha, const Algebraic & beta);

    static GiNaC::ex gcdex(GiNaC::ex f,
                    GiNaC::ex g,
                    const GiNaC::symbol & var,
                    GiNaC::ex & c1,
                    GiNaC::ex & c2);

    static GiNaC::ex sres(const GiNaC::ex & f,
                          const GiNaC::ex & g,
                          const int k,
                          const GiNaC::symbol & var);

    inline PolynomialQQ ParseString(const std::string & s) const;
};

std::ostream & operator<<(std::ostream & output, const PolynomialQQ & p);

#endif // __POLYNOMIALQQ__
