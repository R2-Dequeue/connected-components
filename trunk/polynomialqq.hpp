#ifndef __POLYNOMIALQQ__
#define __POLYNOMIALQQ__

#include <cassert>
#include <string>
#include <vector>
#include <set>
#include <list>

#include <boost/tuple/tuple.hpp>

#include <ginac/ginac.h>

#include "polynomialq.hpp"

/*!
 * \todo Add support for different orderings (lots of functions to add...).
 */
class PolynomialQQ
{
private:

    static const GiNaC::symbol & var1;
    static const GiNaC::symbol & var2;
    GiNaC::ex polynomial;

    struct less;

public:

    typedef std::vector<PolynomialQQ>       vector;
    typedef std::set<PolynomialQQ,less>     set;
    typedef std::list<PolynomialQQ>         list;
    typedef unsigned int                    size_type;

    PolynomialQQ();//!< The default constructor.
    PolynomialQQ(const std::string & s);
    PolynomialQQ(const char * const a);
    PolynomialQQ(const GiNaC::ex & e);
	PolynomialQQ(const GiNaC::numeric & n);

    int degreeX() const;
    int degreeY() const;
    //int degree() const; //!< The degree of the polynomial.
    //bool isMonic() const; //!< True iff the leading coefficient is 1.
    bool isZero() const; //!< True iff the polynomial is '0'.
	//bool isIrreducible() const;
    bool isConstant() const;

    const GiNaC::ex & getEx() const;
//    PolynomialQ getCoeffX(const size_type i) const;
//    PolynomialQ getCoeffY(const size_type i) const;
    PolynomialQQ getDerivative(unsigned int variable) const;

    void switchVariables();

    PolynomialQ subx(const GiNaC::numeric & a) const;
    PolynomialQ suby(const GiNaC::numeric & b) const;

    PolynomialQQ & linearSubs(const GiNaC::numeric & xa,
                              const GiNaC::numeric & xb,
                              const GiNaC::numeric & ya,
                              const GiNaC::numeric & yb);
    PolynomialQQ & differentiate(unsigned int variable);

    template <class T> std::auto_ptr<T> getIrreducibleFactors() const;
    template <class T> T & addIrreducibleFactorsTo(T & factors) const;

    int signAt(const Algebraic & alpha, const Algebraic & beta) const;
    bool signIsZeroAt(const Algebraic & alpha, const Algebraic & beta) const;
    int signAt2(const Algebraic & alpha, const Algebraic & beta) const;

    IntervalQ boundRange(const IntervalQ & iX, const IntervalQ & iY) const;

    static PolynomialQQ::vector
        IrreducibleFactors(const PolynomialQQ::vector & F);
    static PolynomialQ Resultant(const PolynomialQQ & f,
                                 const PolynomialQQ & g,
                                 unsigned int var);
    //static PolynomialQQ Subresultant(const PolynomialQQ & f,
    //                                 const PolynomialQQ & g,
    //                                 const unsigned int k,
    //                                 const unsigned int var);

	PolynomialQQ & operator+=(const PolynomialQQ & rhs);
    PolynomialQQ & operator-=(const PolynomialQQ & rhs);
	PolynomialQQ & operator*=(const PolynomialQQ & rhs);

	std::string getString() const;

    /*!
     * \brief Helper method for internal 'assert' checks.
     * \details This method is public but shouldn't really be published.
     */
    bool Invariants() const;

//protected:

    IntervalQ BoundRange1(const GiNaC::ex & poly,
                          const GiNaC::symbol & var,
                          const IntervalQ & interval) const;

    static Algebraic ANComb(Algebraic alpha,
                            Algebraic beta,
                            const GiNaC::numeric & t);
    static boost::tuple<Algebraic, PolynomialQ, PolynomialQ>
        Simple(const Algebraic & alpha, const Algebraic & beta);
    static boost::tuple<Algebraic, PolynomialQ, PolynomialQ>
        Simple2(const Algebraic & alpha, const Algebraic & beta);

    static GiNaC::ex gcdex(const GiNaC::ex & f,
                    const GiNaC::ex & g,
                    const GiNaC::symbol & var,
                    GiNaC::ex & c1,
                    GiNaC::ex & c2);

    static GiNaC::ex sres(const GiNaC::ex & f,
                          const GiNaC::ex & g,
                          const int k,
                          const GiNaC::symbol & var);

    inline PolynomialQQ ParseString(const std::string & s) const;

    template <typename T>
    inline void ReserveHelper(T & container, unsigned int n) const;
    template <class T>
    inline void InsertMonic(T & factors, const GiNaC::ex & p) const;
};

std::ostream & operator<<(std::ostream & output, const PolynomialQQ & p);

inline PolynomialQQ::PolynomialQQ()
{}

/*!
 * \throws parse_error Thrown by GiNaC if parsing fails (inherits
 *		   invalid_argument).
 */
inline PolynomialQQ::PolynomialQQ(const std::string & s)
{
	polynomial = PolynomialQQ::ParseString(s).polynomial;

    assert(Invariants());
}

/*!
 * \throws parse_error Thrown by GiNaC if parsing fails (inherits
 *		   invalid_argument).
 */
inline PolynomialQQ::PolynomialQQ(const char * const a)
{
	std::string s(a);

    polynomial = PolynomialQQ::ParseString(s).polynomial;

    assert(Invariants());
}

inline PolynomialQQ::PolynomialQQ(const GiNaC::ex & e) : polynomial(e)
{
    if (!Invariants())
        throw std::invalid_argument("PolynomialQQ Constructor: GiNaC expression "
        							"e is not a valid polynomial.");
}

inline PolynomialQQ::PolynomialQQ(const GiNaC::numeric & n) : polynomial(n)
{
    if (!Invariants())
        throw std::invalid_argument("PolynomialQQ Constructor: GiNaC numeric "
        							"n is not a valid rational number.");
}

inline int PolynomialQQ::degreeX() const
{
    assert(Invariants());

    return this->polynomial.degree(this->var1);
}

inline int PolynomialQQ::degreeY() const
{
    assert(Invariants());

    return this->polynomial.degree(this->var2);
}

/*inline int PolynomialQQ::degree() const
{
	assert(Invariants());

	return polynomial.degree(variable);
}

inline bool PolynomialQQ::isMonic() const
{
	assert(Invariants());

	return (polynomial.lcoeff(variable) == 1);
}*/

inline bool PolynomialQQ::isZero() const
{
	assert(Invariants());

	return (polynomial.is_zero()); // (polynomial == 0);
}
/*
inline bool PolynomialQQ::isIrreducible() const
{
	assert(Invariants());

	if (this->degree() > 2)
		return false;

	// Assuming getIrreducibleFactors ignores constant factors.
	if (this->getIrreducibleFactors().size() > 1)
		return false;

	return true;
}
*/
inline bool PolynomialQQ::isConstant() const
{
	// Ways to check if constant:
	// 1) this->degree() == 0 (maybe add a check for -1 in case degree function
	//						   changes behaviour in the future).
	// 2) If the set of symbols in 'polynomial' is empty.
	// 3) If lcoeff == tcoeff.
	// 4) If GiNaC::is_a<GiNaC::numeric>(polynomial) is true.

	assert(Invariants());

	//return (this->degree() == 0 || this->degree() == -1);
	return (GiNaC::is_a<GiNaC::numeric>(polynomial));
}

inline const GiNaC::ex & PolynomialQQ::getEx() const
{
    return polynomial;
}

inline PolynomialQQ PolynomialQQ::getDerivative(unsigned int variable) const
{
    assert(Invariants());
    assert(variable == 1 || variable == 2);

    if (variable != 1 && variable != 2)
        throw std::invalid_argument("Tried to take derivative wrt an invalid variable.");

    PolynomialQQ temp(*this);

    temp.differentiate(variable);

    return temp;
}

inline void PolynomialQQ::switchVariables()
{
    polynomial.subs(GiNaC::lst(var1 == var2, var2 == var1));
}

inline PolynomialQ PolynomialQQ::subx(const GiNaC::numeric & a) const
{
	assert(Invariants());
	assert(a.is_rational()); // throw?

	PolynomialQ p(polynomial.subs(var1 == a).subs(var2 == var1).expand());

	assert(p.Invariants());

	return p;
}

inline PolynomialQ PolynomialQQ::suby(const GiNaC::numeric & b) const
{
	assert(Invariants());
	assert(b.is_rational()); // throw?

	PolynomialQ p(polynomial.subs(var2 == b).expand());

	assert(p.Invariants());

	return p;
}

inline PolynomialQQ & PolynomialQQ::linearSubs(const GiNaC::numeric & xa,
                                       const GiNaC::numeric & xb,
                                       const GiNaC::numeric & ya,
                                       const GiNaC::numeric & yb)
{
    polynomial = polynomial.subs(GiNaC::lst(var1 == xa*var1 + xb,
                                            var2 == ya*var2 + yb)).expand();

    return *this;
}

inline PolynomialQQ & PolynomialQQ::differentiate(unsigned int variable)
{
    assert(Invariants());
    assert(variable == 1 || variable == 2);

    if (variable != 1 && variable != 2)
        throw std::invalid_argument("Tried to take derivative wrt an invalid variable.");

    polynomial = polynomial.diff((variable == 1) ? var1 : var2);

    return *this;
}

template <typename T>
inline void PolynomialQQ::ReserveHelper(T & container, unsigned int n) const
{}

template <>
inline void PolynomialQQ::ReserveHelper< std::vector<Algebraic> >
    (std::vector<Algebraic> & container, unsigned int n) const
{
    container.reserve(n);
}

template <class T>
std::auto_ptr<T> PolynomialQQ::getIrreducibleFactors() const
{
    std::auto_ptr<T> factors(new T);

    //ReserveHelper(*factors, this->degree());

    this->addIrreducibleFactorsTo(*factors);

    return factors;
}

template <class T>
inline void PolynomialQQ::InsertMonic(T & factors, const GiNaC::ex & p) const
{
    const PolynomialQQ temp(p);

    factors.insert(factors.end(), temp);
}


template <class T>
T & PolynomialQQ::addIrreducibleFactorsTo(T & factors) const
{
    assert(Invariants());

    GiNaC::ex p = GiNaC::factor(this->polynomial);

    if (GiNaC::is_a<GiNaC::power>(p))
    {
        InsertMonic(factors, p.op(0));
        return factors;
    }
    else if (!GiNaC::is_a<GiNaC::mul>(p))
    {
        if (!GiNaC::is_a<GiNaC::numeric>(p))
            InsertMonic(factors, p);
        return factors;
    }

    for (GiNaC::const_iterator i = p.begin(); i != p.end(); ++i)
    {
        if (GiNaC::is_a<GiNaC::numeric>(*i))
            continue;
        else if (GiNaC::is_a<GiNaC::power>(*i))
            InsertMonic(factors, (*i).op(0));
        else
            InsertMonic(factors, *i);
	}

    return factors;
}

/*!
 * \details Should work in cases such as 'p += p;'
 */
inline PolynomialQQ & PolynomialQQ::operator+=(const PolynomialQQ & rhs)
{
	assert(Invariants());
	assert(rhs.Invariants());

	polynomial += rhs.polynomial;

	assert(Invariants());

	return *this;
}

/*!
 * \details Should work in cases such as 'p -= p;'
 */
inline PolynomialQQ & PolynomialQQ::operator-=(const PolynomialQQ & rhs)
{
	assert(Invariants());
	assert(rhs.Invariants());

	polynomial -= rhs.polynomial;

	assert(Invariants());

	return *this;
}

/*!
 * \details Should work in cases such as 'p *= p;'
 */
inline PolynomialQQ & PolynomialQQ::operator*=(const PolynomialQQ & rhs)
{
	assert(Invariants());
	assert(rhs.Invariants());

	polynomial = expand(polynomial * rhs.polynomial); // need to figure out 'expand'

	assert(Invariants());

	return *this;
}
/*
inline PolynomialQQ & PolynomialQQ::operator/=(const PolynomialQQ & rhs)
{
	assert(Invariants());
	assert(rhs.Invariants());

	polynomial = expand(polynomial / rhs.polynomial);

	assert(Invariants());

	return *this;
}

inline PolynomialQQ & PolynomialQQ::operator%=(const PolynomialQQ & rhs)
{
	assert(Invariants());
	assert(rhs.Invariants());

	polynomial = expand(polynomial % rhs.polynomial);

	assert(Invariants());

	return *this;
}
*/
inline PolynomialQQ operator+(const PolynomialQQ & lhs, const PolynomialQQ & rhs)
{
	PolynomialQQ temp = lhs;
	temp += rhs;

	return temp;
}

inline PolynomialQQ operator-(const PolynomialQQ & lhs, const PolynomialQQ & rhs)
{
	PolynomialQQ temp = lhs;
	temp -= rhs;

	return temp;
}

inline PolynomialQQ operator*(const PolynomialQQ & lhs, const PolynomialQQ & rhs)
{
	PolynomialQQ temp = lhs;
	temp *= rhs;

	return temp;
}
/*
inline PolynomialQQ operator/(const PolynomialQQ & lhs, const PolynomialQQ & rhs)
{
	PolynomialQQ temp = lhs;
	temp /= rhs;

	return temp;
}

inline PolynomialQQ operator%(const PolynomialQQ & lhs, const PolynomialQQ & rhs)
{
	PolynomialQQ temp = lhs;
	temp %= rhs;

	return temp;
}
*/
inline bool operator==(const PolynomialQQ & lhs, const PolynomialQQ & rhs)
{
    assert(lhs.Invariants());
    assert(rhs.Invariants());

    return (lhs.getEx() == rhs.getEx());
}

inline bool operator!=(const PolynomialQQ & lhs, const PolynomialQQ & rhs)
{
    assert(lhs.Invariants());
    assert(rhs.Invariants());

    return (lhs.getEx() != rhs.getEx());
}

/*!
 * \details Parses the string wrt the internal variables.
 * \throws parse_error Thrown by GiNaC if parsing fails (inherits
 *		   invalid_argument).
 */
inline PolynomialQQ PolynomialQQ::ParseString(const std::string & s) const
{
    GiNaC::symtab table;

    table[var1.get_name()] = var1;
    table[var2.get_name()] = var2;

    GiNaC::parser reader(table, true);

    PolynomialQQ p(reader(s).expand()); // throws an exception if parsing fails.

    return p;
}

inline std::ostream & operator<<(std::ostream & output, const PolynomialQQ & p)
{
    output << p.getEx();

    return output;
}

struct PolynomialQQ::less
    : public std::binary_function<PolynomialQQ, PolynomialQQ, bool>
{
    bool operator()(const PolynomialQQ & lhs, const PolynomialQQ & rhs)
    {
        GiNaC::ex_is_less comp;
        return (comp(lhs.polynomial, rhs.polynomial) < 0);
    }
};

#endif // __POLYNOMIALQQ__
