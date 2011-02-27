/*!
 * \file
 * \author Chris de Pujo
 */

#ifndef __POLYNOMIALQ__
#define __POLYNOMIALQ__

#include <cassert>
#include <string>
#include <vector>
#include <memory>

//#include <boost/numeric/interval.hpp>
#include <boost/foreach.hpp>

#include <ginac/ginac.h>

// //! This type represents the interval in class Algebraic.
//typedef boost::numeric::interval<GiNaC::numeric> IntervalQ;

#include "intervalq.hpp"

class Algebraic; // To allow circular dependancy.

/*!
 * \brief A class representing univariate polynomials over the
 *        (infinite-precision) rationals.
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

    typedef std::vector<PolynomialQ> vector;
    typedef std::set<PolynomialQ>    set;
    typedef std::list<PolynomialQ>   list;
    typedef unsigned int size_type;

	PolynomialQ(); //!< The default constructor.
	PolynomialQ(const std::string & s); //!< Parse a polynomial from a string.
	PolynomialQ(const char * const a);
	PolynomialQ(const GiNaC::ex & e);
	PolynomialQ(const GiNaC::numeric & n);

    size_type degree() const;         //!< The degree of the polynomial.
    bool isMonic() const;       //!< True iff the leading coefficient is 1.
	bool isZero() const;        //!< True iff the polynomial is '0'.
	bool isIrreducible() const; //!< True iff irreducible over the rationals.
	bool isConstant() const;    //!< True iff this is a constant polynomial.

	GiNaC::symbol getVariable() const;
	GiNaC::ex getEx() const;
    GiNaC::numeric getCoeff(const unsigned int i) const;
    PolynomialQ getMonic() const;
    PolynomialQ getDerivative() const;

	PolynomialQ & makeMonic();
    PolynomialQ & differentiate();
    PolynomialQ & negate(); //!< Makes the substituition: f=f(-x).

    //! Returns the number of roots of f between a and b.
    unsigned int sturm(const GiNaC::numeric & a,
    				   const GiNaC::numeric & b) const;
    std::auto_ptr<PolynomialQ::vector> sturmseq() const;

    template <class T> std::auto_ptr<T> getIrreducibleFactors() const;
    template <class T> T & addIrreducibleFactorsTo(T & factors) const;
    template <class T> std::auto_ptr<T> getRoots() const;
    template <class T> T & addRootsTo(T & roots) const;

	//! Calculates the sign of f(a).
    int signAt(const Algebraic & a) const;
    //! Returns the value of the polynomial at 'value'.
	GiNaC::numeric eval(const GiNaC::numeric & value) const;
	PolynomialQ & subs(const GiNaC::ex & e);

    static GiNaC::symbol GetVar() { return variable; }
    static std::auto_ptr<PolynomialQ::vector> sturmseq(const PolynomialQ & p);
    static unsigned int sturm(const PolynomialQ::vector & F,
                              const GiNaC::numeric & a,
                              const GiNaC::numeric & b);
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

    void out();
    std::string toString();

    static void TestClass();

    /*!
     * \brief Helper method for internal 'assert' checks.
     * \detail This method is public but shouldn't really be published.
     */
	bool Invariants() const;

    friend class PolynomialQQ;

//protected:

    template <class T> std::auto_ptr<T> getRootsOfIrreducible() const;
    static std::vector<Algebraic> FindRootsOfIrreducible(const PolynomialQ & p);
    IntervalQ boundRange(const IntervalQ & interval) const;
    static void TestCompare(const PolynomialQ p,
                            const GiNaC::ex expected,
                            unsigned int & count);
    PolynomialQ ParseString(const std::string & s) const;
};

//PolynomialQ SubResultant(const PolynomialQ & f,
//                         const PolynomialQ & g,
//                         const int k,
//                         const GiNaC::symbol & x);

PolynomialQ operator+(const PolynomialQ & lhs, const PolynomialQ & rhs);
PolynomialQ operator-(const PolynomialQ & lhs, const PolynomialQ & rhs);
PolynomialQ operator*(const PolynomialQ & lhs, const PolynomialQ & rhs);
PolynomialQ operator/(const PolynomialQ & lhs, const PolynomialQ & rhs);
PolynomialQ operator%(const PolynomialQ & lhs, const PolynomialQ & rhs);
bool operator==(const PolynomialQ & lhs, const PolynomialQ & rhs);
bool operator!=(const PolynomialQ & lhs, const PolynomialQ & rhs);
PolynomialQ operator*(const PolynomialQ & lhs, const GiNaC::numeric & num);
PolynomialQ operator*(const GiNaC::numeric & num, const PolynomialQ & lhs);
PolynomialQ operator/(const PolynomialQ & lhs, const GiNaC::numeric & num);
PolynomialQ operator/(const GiNaC::numeric & num, const PolynomialQ & lhs);
std::ostream & operator<<(std::ostream & output, const PolynomialQ & p);

#include "algebraic.hpp"

inline PolynomialQ::PolynomialQ() : polynomial(0) {}

/*!
 * \throws parse_error Thrown if string s is an invalid polynomial.
 */
inline PolynomialQ::PolynomialQ(const std::string & s)
{
    polynomial = PolynomialQ::ParseString(s).polynomial;
    assert(Invariants());
}

/*!
 * \throws parse_error Thrown if string s is an invalid polynomial.
 */
inline PolynomialQ::PolynomialQ(const char * const a)
{
	std::string s(a);
	polynomial = PolynomialQ::ParseString(s).polynomial;
    assert(Invariants());
}

/*!
 * \throws invalid_argument Thrown if e is not a valid univariate rational
 *							polynomial in 'variable'.
 */
inline PolynomialQ::PolynomialQ(const GiNaC::ex & e) : polynomial(e)
{
    if (!Invariants())
        throw std::invalid_argument("PolynomialQ Constructor: GiNaC expression "
        							"e is not a valid polynomial.");
}

/*!
 * \throws invalid_argument Thrown if n is not a rational number.
 */
inline PolynomialQ::PolynomialQ(const GiNaC::numeric & n) : polynomial(n)
{
    if (!Invariants())
        throw std::invalid_argument("PolynomialQ Constructor: GiNaC numeric "
        							"n is not a valid rational number.");
}

/*!
 * \detail Just returns the degree from the underlying 'ex' class. This implies
 *         that the degree of the polynomial '0' is 0 and not -1 as it is
 *         sometimes defined.
 */
inline PolynomialQ::size_type PolynomialQ::degree() const
{
	assert(Invariants());

	return polynomial.degree(variable);
}

/*!
 * \detail The zero polynomial is considered monic.
 */
inline bool PolynomialQ::isMonic() const
{
	assert(Invariants());

	return (polynomial.lcoeff(variable) == 1);
}

inline bool PolynomialQ::isZero() const
{
	assert(Invariants());

	return (polynomial.is_zero()); // (polynomial == 0);
}

/*!
 * \todo Are constants irreducible?
 * \todo This check is horribly inefficient; fix it.
 */
inline bool PolynomialQ::isIrreducible() const
{
	assert(Invariants());

	// Assuming getIrreducibleFactors ignores constant factors.
	if (this->getIrreducibleFactors<PolynomialQ::vector>()->size() > 1)
		return false;

	return true;
}

inline bool PolynomialQ::isConstant() const
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

inline GiNaC::symbol PolynomialQ::getVariable() const
{
    assert(Invariants());

    return variable;
}

inline GiNaC::ex PolynomialQ::getEx() const
{
    assert(Invariants());
    return polynomial;
}

/*!
 * \param i i is valid for all unsigned int values.  0 will be returned for
 *			i > this->degree().
 */
inline GiNaC::numeric PolynomialQ::getCoeff(const unsigned int i) const
{
    assert(Invariants());

    GiNaC::ex c = polynomial.coeff(variable, i);

    assert(GiNaC::is_a<GiNaC::numeric>(c));

    GiNaC::numeric num = GiNaC::ex_to<GiNaC::numeric>(c);

    assert(num.is_rational());

    return num;
}

/*!
 * \todo Figure out how to handle the zero polynomial.
 */
inline PolynomialQ PolynomialQ::getMonic() const
{
    assert(Invariants());

    if (!this->isZero())
    {
        PolynomialQ temp((polynomial / polynomial.lcoeff(variable)).expand());
        assert(temp.Invariants());
        return temp;
    }

    PolynomialQ temp(*this);

    return temp;
}

/*!
 * \detail In the case that this is the zero polynomial, then this method
 *		   modifies nothing.
 */
inline PolynomialQ & PolynomialQ::makeMonic()
{
    assert(Invariants());

    if (!this->isZero()) // ==> Leading coefficient is non-zero.
    {
    	polynomial = (polynomial / polynomial.lcoeff(variable)).expand();

    	assert(Invariants());
    }

    return *this;
}

inline PolynomialQ PolynomialQ::getDerivative() const
{
    assert(Invariants());

    PolynomialQ temp(*this);
    temp.differentiate();

    return temp;
}

/*!
 * \detail Note that this modifies the polynomial object; this doesn't just
 *         return a new polynomial that is the derivative.
 */
inline PolynomialQ & PolynomialQ::differentiate()
{
	assert(Invariants());

    polynomial = polynomial.diff(variable);

	return *this;
}

inline PolynomialQ & PolynomialQ::negate()
{
    assert(Invariants());

    polynomial = polynomial.subs(variable == -variable);

    assert(Invariants());

    return *this;
}

inline std::auto_ptr<PolynomialQ::vector> PolynomialQ::sturmseq() const
{
    return PolynomialQ::sturmseq(*this);
}
// at, in, to, over, on, within, inside, into, with, amid, among, amongst,
// inside, into
//
// get, return
//
// store, put, save, add, set
//
// this->StoreIrreducibleFactorsIn(std::vector<Algebraic> & factors);
// this->PutIrreducibleFactorsIn(std::vector<Algebraic> & factors);
// this->SaveIrreducibleFactorsIn(std::vector<Algebraic> & factors);
// this->AddIrreducibleFactorsTo(std::vector<Algebraic> & factors);
//
// this->storeIrreducibleFactorsIn(std::vector<Algebraic> & factors);
// this->putIrreducibleFactorsIn(std::vector<Algebraic> & factors);
// this->saveIrreducibleFactorsIn(std::vector<Algebraic> & factors);
// this->addIrreducibleFactorsTo(std::vector<Algebraic> & factors);
//
//
// this->StoreIrreducibleFactorsAt(std::vector<Algebraic> & factors);
// this->StoreIrreducibleFactorsIn(std::vector<Algebraic> & factors);
// this->PutIrreducibleFactorsAt(std::vector<Algebraic> & factors);
// this->PutIrreducibleFactorsIn(std::vector<Algebraic> & factors);
// this->SaveIrreducibleFactorsAt(std::vector<Algebraic> & factors);
// this->SaveIrreducibleFactorsIn(std::vector<Algebraic> & factors);


// this->addIrreducibleFactorsTo(std::vector<Algebraic> & factors);
// this->getIrreducibleFactors();
// Class::IrreducibleFactors();

template <typename T>
inline void ReserveHelper(T & container, unsigned int n)
{}

template <>
inline void ReserveHelper< std::vector<Algebraic> >
    (std::vector<Algebraic> & container, unsigned int n)
{
    container.reserve(n);
}

/*!
 * \detail Returns the factors of the polynomial.
 * \return Each factor appears only once (multiplicities are ignored). Each
 *		   factor is monic, irreducible, and not a constant.
 */
template <class T>
inline std::auto_ptr<T> PolynomialQ::getIrreducibleFactors() const
{
    std::auto_ptr<T> factors(new T);

    ReserveHelper(*factors, this->degree());

    this->addIrreducibleFactorsTo(*factors);

    return factors;
}

template <class T>
inline void InsertMonic(T & factors, const GiNaC::ex & p)
{
    const PolynomialQ temp(p);

    if (temp.isMonic())
        factors.insert(factors.end(), temp);
    else
        factors.insert(factors.end(), temp.getMonic());
}

template <class T>
T & PolynomialQ::addIrreducibleFactorsTo(T & factors) const
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

template <class T> T & PolynomialQ::addRootsTo(T & roots) const
{
    //std::auto_ptr<PolynomialQ::vector> factors = this->getIrreducibleFactors();
    PolynomialQ::vector factors;
    factors.reserve(this->degree());
    this->addIrreducibleFactorsTo(factors);

    // Assuming all factors are monic (thanks to IrreducibleFactors).

    BOOST_FOREACH(const PolynomialQ & f, factors)
    {
        if (f.degree() <= 0) // just in case
            continue;
        if (f.degree() == 1)
        {
        	const GiNaC::numeric root = -f.getMonic().getCoeff(0);
        	const GiNaC::numeric delta(1, 64);
            //ptr->insert(ptr->end(), Algebraic::MakeWideRational(root));
            roots.insert(roots.end(),
                         Algebraic(PolynomialQ::GetVar() - root,
                                   IntervalQ(root-delta, root+delta)));
        }
        else
        {
            std::auto_ptr<T> ptr2 = f.getRootsOfIrreducible<T>();
            roots.insert(roots.end(), ptr2->begin(), ptr2->end());
        }
    }

    return roots;
}

template <class T>
inline std::auto_ptr<T> PolynomialQ::getRoots() const
{
    std::auto_ptr<T> ptr(new T);

    ReserveHelper(*ptr, this->degree());

    this->addRootsTo(*ptr);

    return ptr;
}

/*!
 * \detail *this must be irreducible.
 */
template <class T>
std::auto_ptr<T> PolynomialQ::getRootsOfIrreducible() const
{
    std::auto_ptr<T> ptr(new T);

    if (this->degree() <= 0)
        return ptr;

    std::auto_ptr<PolynomialQ::vector> F = this->sturmseq();

    GiNaC::numeric R = 0;

    for (unsigned int i = 0; i < this->degree(); i++)
        if (GiNaC::abs(this->getCoeff(i)) > R)
            R = GiNaC::abs(this->getCoeff(i));

    R /= GiNaC::abs(this->getCoeff(this->degree()));
    R += 1;

    GiNaC::numeric L = -R, U = R;
    GiNaC::numeric numroots = PolynomialQ::sturm(*F, L, U);
    GiNaC::numeric roots, L1, U1, NU;

    if (numroots == 0)
        return ptr;

    //numbers.reserve(numroots.to_int());
    ReserveHelper(*ptr, numroots.to_int());

    while (true)
    {
        numroots = PolynomialQ::sturm(*F, L, U);

        if (numroots == 1)
        {
            //numbers.push_back(Algebraic(p, IntervalQ(L, U)));
            ptr->insert(ptr->end(), Algebraic(*this, IntervalQ(L, U)));
            break;
        }

        L1 = L;
        U1 = U;

        while (true)
        {
            NU = (L1+U1)/2;
            roots = PolynomialQ::sturm(*F,L1,NU);

            if (roots < numroots)
            {
                if (roots == 1)
                {
                    //numbers.push_back(Algebraic(p, IntervalQ(L1, NU)));
                    ptr->insert(ptr->end(), Algebraic(*this, IntervalQ(L1, NU)));
                    U1 = NU;
                    break;
                }
                else if (roots == 0)
                    L1 = NU;
                else
                    U1 = NU;
            }
            else
                U1 = NU;
        }

        L = U1;
    }

    return ptr;
}

inline PolynomialQ & PolynomialQ::subs(const GiNaC::ex & e)
{
    assert(Invariants());

    polynomial = polynomial.subs(variable == e).expand();
    this->makeMonic();

    assert(Invariants());

    return *this;
}

/*!
 * \detail Should work in cases such as 'p += p;'
 */
inline PolynomialQ & PolynomialQ::operator+=(const PolynomialQ & rhs)
{
	assert(Invariants());
	assert(rhs.Invariants());

	polynomial += rhs.polynomial;

	assert(Invariants());

	return *this;
}

/*!
 * \detail Should work in cases such as 'p -= p;'
 */
inline PolynomialQ & PolynomialQ::operator-=(const PolynomialQ & rhs)
{
	assert(Invariants());
	assert(rhs.Invariants());

	polynomial -= rhs.polynomial;

	assert(Invariants());

	return *this;
}

/*!
 * \detail Should work in cases such as 'p *= p;'
 */
inline PolynomialQ & PolynomialQ::operator*=(const PolynomialQ & rhs)
{
	assert(Invariants());
	assert(rhs.Invariants());

	polynomial = expand(polynomial * rhs.polynomial); // need to figure out 'expand'

	assert(Invariants());

	return *this;
}

inline PolynomialQ & PolynomialQ::operator/=(const PolynomialQ & rhs)
{
    assert(Invariants());
    assert(rhs.Invariants());

    polynomial = quo(polynomial, rhs.polynomial, variable);

    assert(Invariants());

    return *this;
}

inline PolynomialQ & PolynomialQ::operator%=(const PolynomialQ & rhs)
{
    assert(Invariants());
    assert(rhs.Invariants());

    polynomial = rem(polynomial, rhs.polynomial, variable);

    assert(Invariants());

    return *this;
}

inline PolynomialQ operator+(const PolynomialQ & lhs, const PolynomialQ & rhs)
{
	PolynomialQ temp = lhs;
	temp += rhs;

	return temp;
}

inline PolynomialQ operator-(const PolynomialQ & lhs, const PolynomialQ & rhs)
{
	PolynomialQ temp = lhs;
	temp -= rhs;

	return temp;
}

inline PolynomialQ operator*(const PolynomialQ & lhs, const PolynomialQ & rhs)
{
	PolynomialQ temp = lhs;
	temp *= rhs;

	return temp;
}

inline PolynomialQ operator/(const PolynomialQ & lhs, const PolynomialQ & rhs)
{
	PolynomialQ temp = lhs;
	temp /= rhs;

	return temp;
}

inline PolynomialQ operator%(const PolynomialQ & lhs, const PolynomialQ & rhs)
{
	PolynomialQ temp = lhs;
	temp %= rhs;

	return temp;
}

inline bool operator==(const PolynomialQ & lhs, const PolynomialQ & rhs)
{
    assert(lhs.Invariants());
    assert(rhs.Invariants());

    return (lhs.getEx() == rhs.getEx());
}

inline bool operator!=(const PolynomialQ & lhs, const PolynomialQ & rhs)
{
    assert(lhs.Invariants());
    assert(rhs.Invariants());

    return (lhs.getEx() != rhs.getEx());
}

inline PolynomialQ operator*(const PolynomialQ & lhs, const GiNaC::numeric & num)
{
    assert(lhs.Invariants());
    assert(num.is_rational());

    PolynomialQ temp(lhs.getEx()*GiNaC::ex(num));

    return temp;
}

inline PolynomialQ operator*(const GiNaC::numeric & num, const PolynomialQ & lhs)
{
    return lhs*num;
}

inline PolynomialQ operator/(const PolynomialQ & lhs, const GiNaC::numeric & num)
{
    assert(lhs.Invariants());
    assert(num.is_rational());

    PolynomialQ temp(lhs.getEx()/GiNaC::ex(num));

    return temp; // rely on GiNaC throwing an exception if num == 0.
}

inline PolynomialQ operator/(const GiNaC::numeric & num, const PolynomialQ & lhs)
{
    return lhs/num; // let other method throw
}

inline std::ostream & operator<<(std::ostream & output, const PolynomialQ & p)
{
    output << p.getEx();

    return output;
}

#endif // __POLYNOMIALQ__
