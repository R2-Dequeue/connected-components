/*!
 * \file
 * \author Chris de Pujo
 */

#include "polynomialq.hpp"

#include <string>
#include <cassert>
#include <numeric>
#include <functional>
#include <stdexcept>

#include "polynomialbase.hpp"

#include <boost/foreach.hpp>

#include "algebraic.hpp"

// Field to use
// Number of variables
//
// Complex      C,i
// Reals        R
// Algebraic    A
// Rationals    Q,F,R
// Integers     Z,I

const GiNaC::symbol & PolynomialQ::variable = PolynomialBase::var1;

PolynomialQ::PolynomialQ(const std::string & s)
{
    using namespace GiNaC;

    symtab table;

    table[variable.get_name()] = variable;

    parser reader(table);

    // parser.strict = true; // must be a typo in the tutorial
    reader.strict = true;

    // reader(s).expand(); ?
    polynomial = reader(s); // throws an exception if parsing fails

    assert(Invariant());

    if (!Invariant())
        throw std::invalid_argument("PolynomialQ Constructor: Parsing"
                                    "univariate polynomial over the rationals"
                                    "failed.");
}

PolynomialQ::PolynomialQ(const char * const s)
{
    using namespace GiNaC;

    symtab table;

    table[variable.get_name()] = variable;

    parser reader(table);

    // parser.strict = true; // must be a typo in the tutorial
    reader.strict = true;

    std::string a(s);

    // reader(s).expand(); ?
    polynomial = reader(a); // throws an exception if parsing fails

    assert(Invariant());

    if (!Invariant())
        throw std::invalid_argument("PolynomialQ Constructor: Parsing"
                                    "univariate polynomial over the rationals"
                                    "failed.");
}

PolynomialQ::PolynomialQ(const GiNaC::ex & e)
    : polynomial(e)
{
    if (!Invariant())
        throw std::invalid_argument("Parsing of polynomial succeded, but the"
                                    "result is not canonical.");
}

PolynomialQ::PolynomialQ(const GiNaC::numeric & n)
    : polynomial(n)
{
    if (!Invariant())
        throw std::invalid_argument("Parsing of polynomial succeded, but the"
                                    "result is not canonical.");
}

/*!
 * \detail This method checks that the internal state is free of errors. It is
 *         intended to be used only with 'assert' at the begining of its
 *         methods ( assert(Invariant()); ).
 * \todo Add check for irreducibility.
 */
bool PolynomialQ::Invariant() const
{
    using namespace GiNaC;

	if (this == NULL)
		return false;

    // is this false if there are additional variables?
    if (!polynomial.is_polynomial(variable))
        return false;

    if (!polynomial.info(info_flags::rational_polynomial))
        return false;

    //if (polynomial == ex()) // should always be at least 'ex(0)'.
    //    return false;

    int deg = polynomial.degree(variable);

    for (int i = 0; i <= deg; i++)
    {
        if (is_a<numeric>(polynomial.coeff(variable, i)))
        {
            numeric c = ex_to<numeric>(polynomial.coeff(variable, i));

            if (!c.is_rational())
                return false;
        }
        else
            return false;
    }

    return true;
}

/*!
 * \detail Just returns the degree from the underlying 'ex' class. This implies
 *         that the degree of the polynomial '0' is 0 and not -1 as it is
 *         sometimes defined.
 */
inline int PolynomialQ::degree() const
{
	assert(Invariant());

	return polynomial.degree(variable);
}

inline bool PolynomialQ::isMonic() const
{
	assert(Invariant());

	return (polynomial.lcoeff(variable) == 1);
}

inline bool PolynomialQ::isZero() const
{
	assert(Invariant());

	return (polynomial == 0);
}

bool PolynomialQ::isIrreducible() const
{
	assert(Invariant());
	
	if (this->degree() > 2)
		return false;
	
	// Assuming getIrreducibleFactors ignores constant factors.
	if (this->getIrreducibleFactors().size() > 1)
		return false;
	
	return true;
}

GiNaC::numeric PolynomialQ::getCoeff(const unsigned int i) const
{
    assert(Invariant());

    GiNaC::ex c = polynomial.coeff(variable, i);

    assert(GiNaC::is_a<GiNaC::numeric>(c));

    GiNaC::numeric num = GiNaC::ex_to<GiNaC::numeric>(c);

    assert(num.is_rational());

    return num;
}

PolynomialQ PolynomialQ::getMonic()
{
    assert(Invariant());

    PolynomialQ temp((polynomial / polynomial.lcoeff(variable)).expand());

    assert(temp.Invariant());

    return temp;
}

PolynomialQ & PolynomialQ::makeMonic()
{
    assert(Invariant());

    polynomial = (polynomial / polynomial.lcoeff(variable)).expand();

    assert(Invariant());

    return *this;
}

PolynomialQ PolynomialQ::getDerivative() const
{
    assert(Invariant());

    PolynomialQ temp(*this);

    temp.differentiate();

    return temp;
}

/*!
 * \detail Note that this modifies the polynomial object; this doesn't just
 *         return a new polynomial that is the derivative.
 */
PolynomialQ & PolynomialQ::differentiate()
{
	assert(Invariant());

    polynomial = polynomial.diff(variable);

	return *this;
}

/*!
 * \detail Returns a vector containing the unique irreducible factors of the
 *         polynomial.
 * \return Each element should only appear once.
 * \todo Add more error checking and experiment with 'factor(GiNaC::ex)'.
 *		 How does 'factor' handle constant factors?
 */
std::vector<PolynomialQ> PolynomialQ::getIrreducibleFactors() const
{
    assert(Invariant());

    GiNaC::ex p = factor(polynomial);

    std::vector<PolynomialQ> factors;

    for (GiNaC::const_iterator i = p.begin(); i != p.end(); i++)
        factors.push_back(PolynomialQ(GiNaC::is_a<GiNaC::power>(*i) ?
                                      (*i).op(0) :
                                      *i));
        // The PolyQ constructor will throw if the ex is not a valid poly.

    return factors;
}

GiNaC::ex PolynomialQ::getEx() const
{
    GiNaC::ex temp(polynomial);

    return temp;
}

int PolynomialQ::signAt(const Algebraic & a) const
{
    assert(Invariant());

    PolynomialQ remainder(*this % a.polynomial);

    if (remainder == 0)
        return 0;

    Algebraic alpha(a);

    while (true)
    {
        IntervalQ Y = this->boundRange(alpha.getInterval());

        if (Y.upper() < 0) return -1;
        if (Y.lower() > 0) return  1;

        alpha.tightenInterval();
    }

    assert(false);

    return -41;
}

/*!
 * \param value Must be a rational number.
 * \throws invalid_argument If value is not rational.
 * \return Will be a rational number.
 */
GiNaC::numeric PolynomialQ::eval(const GiNaC::numeric & value) const
{
    assert(Invariant());
    
    if (!GiNaC::is_a<GiNaC::numeric>(value))
    	throw std::invalid_argument("PolynomialQ::eval: value passed is not"
    								"rational.");

    GiNaC::ex temp(polynomial.subs(variable == value));

    assert(GiNaC::is_a<GiNaC::numeric>(temp));

    return GiNaC::ex_to<GiNaC::numeric>(temp); // Will this work? ex_to returns
    										   // a reference.  Hopefully a new
    										   // object is created.
}

/*!
 * \return A vector of the roots of the polynomials in P.  The elements in the
 *		   returned vector are unique; multiple roots or roots that appear in
 *		   more than one polynomial are represented only once.
 * \todo Find resource (preferably consice) for exceptions thrown by the STL.
 */
std::vector<Algebraic> PolynomialQ::FindRoots(const std::vector<PolynomialQ> P)
{
    /*
    How should I handle 0 or other constant polynomials?
    */
    std::vector<PolynomialQ> factors = PolynomialQ::IrreducibleFactors(P);
    std::set<Algebraic> numberSet; // Use this to sort alpebraic numbers.

    if (factors.size() == 1)
        if (factors[0].isConstant())
        {
            std::vector<Algebraic> temp; // Empty; constant poly has no roots.
            return temp;
        }

    BOOST_FOREACH(const PolynomialQ & f, factors)
    {
        if (f.degree() == 1)
        {
            GiNaC::ex numberEx = GiNaC::lsolve(f.getEx() == 0, f.getVar());

            if (GiNaC::is_a<GiNaC::numeric>(numberEx))
            {
                GiNaC::numeric number = GiNaC::ex_to<GiNaC::numeric>(numberEx);

                if (number.is_rational())
                {
                    Algebraic a(f, IntervalQ(number));
                    numberSet.insert(a);
                    continue;
                }
            }

            throw std::runtime_error();

            /*
            GiNaC::ex numberEx = GiNaC::lsolve(f.getEx() == 0, f.getVar());
            GiNaC::numeric number = GiNaC::ex_to<GiNaC::numeric>(numberEx);

            if (!number.is_rational())
                throw std::runtime_error();

            Algebraic a(f, IntervalQ(number));
            numberSet.insert(a);
            continue;
            */
        }
        else (f.degree() == 2)
        {
            // Check for complex roots
            // \sqrt{b^2-4ac}; or if a=1, \sqrt{b^2-4c}
            GiNaC::numeric discriminant = f.coeff(1)^2 - 4*f.coeff(2)*f.coeff(0);

            if (discriminant >= 0) // ==> We have real roots.
            {
                // Know that f.coeff(2) != 0.
                GiNaC::numeric left = -f.coeff(1)/(f.coeff(2)*2);
                GiNaC::numeric right = GiNaC::sqrt(discriminant)/(f.coeff(2)*2);
            }
        }
        else
            throw std::runtime_error("PolynomialQ::FindRoots: Got a factor not"
                                     "of degree 1 or 2.");
    }

    std::vector<Algebraic> numbers(numberSet.begin(), numberSet.end());
    // numbers.reserve(numberSet.size());

    for (int i = 0; i < numbers.size()-1; i++)
        Algebraic::SeparateIntervals(numbers[i], numbers[i+1]);

    return numbers;
}

/*!
 * \detail Returns a vector of unique irreducible factors of the polynomials
 *         in F.
 * \return Each element should only appear once.
 * \todo Add more error checking and handle zero-polynomials.
 *		 Also, decide how to handle constant factors of polynomials.
 */
std::vector<PolynomialQ>
    PolynomialQ::IrreducibleFactors(const std::vector<PolynomialQ> & F)
{
    // remove '0' polynomials? Monic check?
    PolynomialQ fp = accumulate(F.begin(),
                                F.end(),
                                PolynomialQ((GiNaC::numeric)1),
                                std::multiplies<PolynomialQ>());
    // implement PolynomialQ.mul(set)?

    return p.getIrreducibleFactors();
}

/*!
 * \return A rational number.
 */
GiNaC::numeric
    PolynomialQ::Resultant(const PolynomialQ & f, const PolynomialQ & g)
{
    assert(f.Invariant());
    assert(g.Invariant());

    GiNaC::ex res = resultant(f.polynomial, g.polynomial, variable);
    
    assert(GiNaC::is_a<GiNaC::numeric>(res));
    
    GiNaC::numeric num = GiNaC::ex_to<GiNaC::numeric>(res);
    
    assert(num.is_rational());

    return num;
}

IntervalQ PolynomialQ::boundRange(const IntervalQ & interval) const
{
    assert(Invariant());

    unsigned int d = degree();
    IntervalQ range(getCoeff(d));

    for (int i = d-1; i >= 0; i--)
    {
        range *= interval;
        range += IntervalQ(getCoeff(i));
    }

    assert(range.lower().is_rational());
    assert(range.upper().is_rational());

    return range;
}

void PolynomialQ::TestClass()
{
    unsigned int errCount = 0;

    PolynomialQ x1("x-1"), x2("x-2"),
                x3("x-3"), pzero;

    std::cout << "Testing class PolynomialQ." << std::endl;

    TestCompare(x1, variable - 1, errCount);
    TestCompare(x2, variable - 2, errCount);
    TestCompare(x3, variable - 3, errCount);

    TestCompare(x1 + x2, variable*2 - 3, errCount);
    TestCompare(x2 + x3, variable*2 - 5, errCount);
    TestCompare(x1 + x3, variable*2 - 4, errCount);

    TestCompare(x1 + x1, variable*2 - 2, errCount);
    TestCompare(x2 + x2, variable*2 - 4, errCount);
    TestCompare(x3 + x3, variable*2 - 6, errCount);

    TestCompare(x1 * x1, variable*variable + variable*(-2) + 1, errCount);
    TestCompare(x2 * x2, variable*variable + variable*(-4) + 4, errCount);
    TestCompare(x3 * x3, variable*variable + variable*(-6) + 9, errCount);

    TestCompare(x1 * (GiNaC::numeric)0, 0, errCount);
    TestCompare(x2 * (GiNaC::numeric)0, 0, errCount);
    TestCompare(x3 * (GiNaC::numeric)0, 0, errCount);

    TestCompare(x1 * (GiNaC::numeric)1, variable - 1, errCount);
    TestCompare(x2 * (GiNaC::numeric)1, variable - 2, errCount);
    TestCompare(x3 * (GiNaC::numeric)1, variable - 3, errCount);

    TestCompare(x1 - x2, 1, errCount);
    TestCompare(x2 - x3, 1, errCount);
    TestCompare(x1 - x3, 2, errCount);

    TestCompare(x1 - x1, 0, errCount);
    TestCompare(x2 - x2, 0, errCount);
    TestCompare(x3 - x3, 0, errCount);

    TestCompare(x1.getDerivative(), 1, errCount);
    TestCompare(x2.getDerivative(), 1, errCount);
    TestCompare(x3.getDerivative(), 1, errCount);

    TestCompare(pzero, 0, errCount);
    TestCompare(pzero.getDerivative(), 0, errCount);

    PolynomialQ p1("5*x^3-x^2+7");

    TestCompare(p1, pow(variable,3)*5-pow(variable,2)+7, errCount);
    TestCompare(p1.getMonic(),
                pow(variable,3)-pow(variable,2)/5+GiNaC::numeric(7,5),
                errCount);
    TestCompare(p1.differentiate(),
                15*pow(variable,2)-2*variable,
                errCount);
    TestCompare(p1.differentiate(),
                30*variable-2,
                errCount);
    TestCompare(p1.differentiate(),
                30,
                errCount);

    TestCompare(x1.getMonic(), variable - 1, errCount);
    TestCompare(x2.getMonic(), variable - 2, errCount);
    TestCompare(x3.getMonic(), variable - 3, errCount);

    TestCompare(PolynomialQ("x-5")*PolynomialQ("x-7"),
                pow(variable,2)-12*variable+35,
                errCount);

    if (errCount == 0)
        std::cout << "All tests on class PolynomialQ passed." << std::endl;
    else
        std::cout << errCount << " tests on class PolynomialQ failed."
                  << std::endl;
}

void PolynomialQ::TestCompare(const PolynomialQ p,
                              const GiNaC::ex expected,
                              unsigned int & count)
{
    if (p.polynomial != expected)
    {
        std::cout << "Test failed. Expected \"" << expected << "\" but got: \""
                  << p.polynomial << "\"." << std::endl;

        count++;
    }
}

std::ostream & operator<<(std::ostream & output, const PolynomialQ & p)
{
    output << p.getEx() << std::endl;

    return output;
}

/*!
 * \detail Should work in cases such as 'p += p;'
 */
PolynomialQ & PolynomialQ::operator+=(const PolynomialQ & rhs)
{
	assert(Invariant());
	assert(rhs.Invariant());

	polynomial += rhs.polynomial;

	assert(Invariant());

	return *this;
}

/*!
 * \detail Should work in cases such as 'p -= p;'
 */
PolynomialQ & PolynomialQ::operator-=(const PolynomialQ & rhs)
{
	assert(Invariant());
	assert(rhs.Invariant());

	polynomial -= rhs.polynomial;

	assert(Invariant());

	return *this;
}

/*!
 * \detail Should work in cases such as 'p *= p;'
 */
PolynomialQ & PolynomialQ::operator*=(const PolynomialQ & rhs)
{
	assert(Invariant());
	assert(rhs.Invariant());

	polynomial = expand(polynomial * rhs.polynomial); // need to figure out 'expand'

	assert(Invariant());

	return *this;
}

PolynomialQ & PolynomialQ::operator/=(const PolynomialQ & rhs)
{
    assert(Invariant());
    assert(rhs.Invariant());

    polynomial = quo(polynomial, rhs.polynomial, variable);

    assert(Invariant());

    return *this;
}

PolynomialQ operator+(const PolynomialQ & lhs, const PolynomialQ & rhs)
{
	assert(lhs.Invariant());
	assert(rhs.Invariant());

	return PolynomialQ(lhs.polynomial + rhs.polynomial);
}

PolynomialQ operator-(const PolynomialQ & lhs, const PolynomialQ & rhs)
{
	assert(lhs.Invariant());
	assert(rhs.Invariant());

	return PolynomialQ(lhs.polynomial - rhs.polynomial);
}

PolynomialQ operator*(const PolynomialQ & lhs, const PolynomialQ & rhs)
{
	assert(lhs.Invariant());
	assert(rhs.Invariant());

    // need to figure out 'expand'
	return PolynomialQ(expand(lhs.polynomial * rhs.polynomial));
}

PolynomialQ operator/(const PolynomialQ & lhs, const PolynomialQ & rhs)
{
    assert(lhs.Invariant());
    assert(rhs.Invariant());

    return PolynomialQ(quo(lhs.polynomial,
                           rhs.polynomial,
                           PolynomialQ::variable));
}

PolynomialQ operator%(const PolynomialQ & lhs, const PolynomialQ & rhs)
{
    assert(lhs.Invariant());
    assert(rhs.Invariant());

    return PolynomialQ(rem(lhs.polynomial,
                           rhs.polynomial,
                           PolynomialQ::variable));
}

PolynomialQ operator*(const PolynomialQ & lhs, const GiNaC::numeric & num)
{
    assert(lhs.Invariant());
    assert(num.is_rational());

    return lhs.polynomial*(GiNaC::ex)num; // will convert
}

PolynomialQ operator/(const PolynomialQ & lhs, const GiNaC::numeric & num)
{
    assert(lhs.Invariant());
    assert(num.is_rational());

    return lhs.polynomial/(GiNaC::ex)num; // rely on GiNaC throwing an exception if num==0.
}
