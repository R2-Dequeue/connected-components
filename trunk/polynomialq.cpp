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

/*!
 * \throws parse_error Thrown if string s is an invalid polynomial.
 */
PolynomialQ::PolynomialQ(const std::string & s)
{
    polynomial = PolynomialQ::ParseString(s).polynomial;

    assert(Invariants());
}

/*!
 * \throws parse_error Thrown if string s is an invalid polynomial.
 */
PolynomialQ::PolynomialQ(const char * const a)
{
	std::string s(a);

	polynomial = PolynomialQ::ParseString(s).polynomial;

    assert(Invariants());
}

/*!
 * \throws invalid_argument Thrown if e is not a valid univariate rational
 *							polynomial in 'variable'.
 */
PolynomialQ::PolynomialQ(const GiNaC::ex & e)
    : polynomial(e)
{
    if (!Invariants())
        throw std::invalid_argument("PolynomialQ Constructor: GiNaC expression"
        							"e is not a valid polynomial.");
}

/*!
 * \throws invalid_argument Thrown if n is not a rational number.
 */
PolynomialQ::PolynomialQ(const GiNaC::numeric & n)
    : polynomial(n)
{
    if (!Invariants())
        throw std::invalid_argument("PolynomialQ Constructor: GiNaC numeric"
        							"n is not a valid rational number.");
}

/*!
 * \detail Just returns the degree from the underlying 'ex' class. This implies
 *         that the degree of the polynomial '0' is 0 and not -1 as it is
 *         sometimes defined.
 */
int PolynomialQ::degree() const
{
	assert(Invariants());

	return polynomial.degree(variable);
}

/*!
 * \detail The zero polynomial is considered monic.
 */
bool PolynomialQ::isMonic() const
{
	assert(Invariants());

	return (polynomial.lcoeff(variable) == 1);
}

bool PolynomialQ::isZero() const
{
	assert(Invariants());

	return (polynomial.is_zero()); // (polynomial == 0);
}

bool PolynomialQ::isIrreducible() const
{
	assert(Invariants());

	if (this->degree() > 2)
		return false;

	// Assuming getIrreducibleFactors ignores constant factors.
	if (this->getIrreducibleFactors().size() > 1)
		return false;

	return true;
}

bool PolynomialQ::isConstant() const
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

/*!
 * \param i i is valid for all unsigned int values.  0 will be returned for
 *			i > this->degree().
 */
GiNaC::numeric PolynomialQ::getCoeff(const unsigned int i) const
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
PolynomialQ PolynomialQ::getMonic() const
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
PolynomialQ & PolynomialQ::makeMonic()
{
    assert(Invariants());

    if (!this->isZero()) // ==> Leading coefficient is non-zero.
    {
    	polynomial = (polynomial / polynomial.lcoeff(variable)).expand();

    	assert(Invariants());
    }

    return *this;
}

PolynomialQ PolynomialQ::getDerivative() const
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
PolynomialQ & PolynomialQ::differentiate()
{
	assert(Invariants());

    polynomial = polynomial.diff(variable);

	return *this;
}

/*!
 * \detail a must be less than b.
 */
unsigned int PolynomialQ::sturm(const GiNaC::numeric & a,
                                const GiNaC::numeric & b) const
{
    PolynomialQ p2(*this);
    PolynomialQ p1(this->getDerivative());

    GiNaC::numeric f2 = p2.eval(a), g2 = p2.eval(b);
    GiNaC::numeric fa = p1.eval(a), gb = p1.eval(b);
    GiNaC::numeric va = 0,			vb = 0;

    if (f2 != 0) // if == 0 then do nothing
    {
    	if (fa == 0)
    		fa == f2;
    	else
    		if (csgn(fa) != csgn(f2))
    			va = 1;
    }

    if (g2 != 0) // if == 0 then do nothing
    {
    	if (gb == 0)
    		gb == g2;
    	else
    		if (csgn(gb) != csgn(g2))
    			vb = 1;
    }

    while (!p1.isZero())
    {
    	PolynomialQ rem = (p1 % p2)*(-1);
    	GiNaC::numeric alpha = rem.eval(a), beta = rem.eval(b);

    	if (alpha != 0)
    	{
    		if (fa != 0)
    			if (csgn(alpha) != csgn(fa))
    				va++;
    		fa = alpha;
    	}

    	if (beta != 0)
    	{
    		if (gb != 0)
    			if (csgn(beta) != csgn(gb))
    				vb++;
    		gb = beta;
    	}

    	p2 = p1;
    	p1 = rem;
    }

    assert((va-vb).is_nonneg_integer());

    return (va - vb).to_int();
}

/*!
 * \detail Returns the factors of the polynomial in a vector.
 * \return Each factor appears only once (multiplicities are ignored). Each
 *		   factor is monic, irreducible, and not a constant.
 * \todo Add more error checking and experiment with 'factor(GiNaC::ex)'.
 *		 How does 'factor' handle constant factors?
 * \todo Use a set and a comparison function.
 */
std::vector<PolynomialQ> PolynomialQ::getIrreducibleFactors() const
{
    assert(Invariants());

    GiNaC::ex p = factor(polynomial);

    std::vector<PolynomialQ> factors; // Make a comparison function and change
    								  // this to use a set.

    for (GiNaC::const_iterator i = p.begin(); i != p.end(); i++)
    {
        // The PolyQ constructor will throw if the ex is not a valid poly.

        if (GiNaC::is_a<GiNaC::numeric>(*i))
        {
        	// do nothing
        	continue;
        }
        else if (GiNaC::is_a<GiNaC::power>(*i))
        {
        	PolynomialQ p((*i).op(0));

        	factors.push_back(p.getMonic());
        }
        else
        {
        	PolynomialQ p(*i);

        	factors.push_back(p.getMonic());
        }
	}

    return factors;
}

int PolynomialQ::signAt(const Algebraic & a) const
{
    assert(Invariants());

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
    assert(Invariants());

    if (!value.is_rational())
    	throw std::invalid_argument("PolynomialQ::eval: value passed is not"
    								"rational.");

    GiNaC::ex temp(polynomial.subs(variable == value));

    assert(GiNaC::is_a<GiNaC::numeric>(temp));

    return GiNaC::ex_to<GiNaC::numeric>(temp); // Will this work? ex_to returns
    										   // a reference.  Hopefully a new
    										   // object is created.
}

/*!
 * \param P A vector of non-zero polynomials. Constant polynomials are valid,
 *			but will not contribute anything to the output.
 * \return A vector of the roots of the polynomials in P.  The elements in the
 *		   returned vector are unique; multiple roots or roots that appear in
 *		   more than one polynomial are represented only once.
 * \todo Find resource (preferably consice) for exceptions thrown by the STL.
 * \todo Optimize by replacing the call to SeparateIntervals with custom code.
 * \todo Check if GiNaC::numeric uses copy-on-write like GiNaC::ex.
 */
std::vector<Algebraic> PolynomialQ::FindRoots(const std::vector<PolynomialQ> P)
{
    std::vector<PolynomialQ> factors = PolynomialQ::IrreducibleFactors(P);
    std::set<Algebraic> numberSet; // Use this to sort the alpebraic numbers.

    // Assume all factors are monic (thanks to IrreducibleFactors).

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
        	const GiNaC::numeric root = -f.getCoeff(0);

            Algebraic alpha(f, IntervalQ(root, root));

            numberSet.insert(alpha);
        }
        else if (f.degree() == 2)
        {
            const GiNaC::numeric b = f.getCoeff(1);
            const GiNaC::numeric c = f.getCoeff(0);
            const GiNaC::numeric discriminant = pow(b,(GiNaC::numeric)2) - 4*c;

            if (discriminant > 0) // ==> We have real roots.
            {
                const GiNaC::numeric left = -b/2;
                const GiNaC::numeric right = GiNaC::sqrt(discriminant) / 2;
                // Note: the above is always > 0.

                Algebraic alpha(f, IntervalQ( left 		, left+right	));
                Algebraic beta (f, IntervalQ( left-right	, left			));

                Algebraic::SeparateIntervals(alpha, beta);

                numberSet.insert(alpha);
                numberSet.insert(beta);
            }
            else if (discriminant == 0) // ==> Just one (real) root.
            {
            	const GiNaC::numeric root = -b/2;

            	Algebraic alpha(f.getVariable() - root, IntervalQ(root , root));

            	numberSet.insert(alpha);
            }
        }
        else
            throw std::runtime_error("PolynomialQ::FindRoots: Got a factor not"
                                     "of degree 1 or 2.");
    }

    std::vector<Algebraic> numbers(numberSet.begin(), numberSet.end());

    for (std::vector<Algebraic>::size_type i = 0; i < numbers.size()-1; i++)
        Algebraic::SeparateIntervals(numbers[i], numbers[i+1]);

    return numbers;
}

/*!
 * \detail Returns the factors of the polynomials in F in a vector.
 * \param F A vector of non-zero polynomials. Constant polynomials are valid,
 *			but will not contribute anything to the output.
 * \return Each factor appears only once (multiplicities and common roots are
 *		   ignored). Each factor is monic, irreducible, and not a constant.
 */
std::vector<PolynomialQ>
    PolynomialQ::IrreducibleFactors(const std::vector<PolynomialQ> & F)
{
    // Probably faster to get factors separately and merge.
    PolynomialQ fp = accumulate(F.begin(),
                                F.end(),
                                PolynomialQ((GiNaC::ex)1),
                                std::multiplies<PolynomialQ>());
    // implement PolynomialQ::mul(set)?

    assert(!fp.isZero());

    return fp.getIrreducibleFactors();
}

/*!
 * \return A rational number.
 */
GiNaC::numeric
    PolynomialQ::Resultant(const PolynomialQ & f, const PolynomialQ & g)
{
    assert(f.Invariants());
    assert(g.Invariants());

    GiNaC::ex res = resultant(f.polynomial, g.polynomial, variable);

    assert(GiNaC::is_a<GiNaC::numeric>(res));

    GiNaC::numeric num = GiNaC::ex_to<GiNaC::numeric>(res);

    assert(num.is_rational());

    return num;
}

IntervalQ PolynomialQ::boundRange(const IntervalQ & interval) const
{
    assert(Invariants());

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

/*!
 * \detail This method checks that the internal state is free of errors. It is
 *         intended to be used only with 'assert' at the begining of its
 *         methods ( assert(Invariants()); ).
 * \todo Lookup how to get the set of symbols in an expression.
 */
bool PolynomialQ::Invariants() const
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
 * \detail Parses the string wrt the internal variable.
 * \throws parse_error Thrown by GiNaC if parsing fails (inherits
 *		   invalid_argument).
 * \todo Report 'parser.strict' typo in the tutorial.
 */
PolynomialQ PolynomialQ::ParseString(const std::string & s) const
{
	GiNaC::symtab table;

    table[variable.get_name()] = variable;

    GiNaC::parser reader(table); // reader(table, true);
    reader.strict = true; // Tells reader to throw if variables besides
    					  // 'variable' appear in the passed string.

    PolynomialQ p(reader(s).expand()); // throws an exception if parsing fails

    return p;
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
	assert(Invariants());
	assert(rhs.Invariants());

	polynomial += rhs.polynomial;

	assert(Invariants());

	return *this;
}

/*!
 * \detail Should work in cases such as 'p -= p;'
 */
PolynomialQ & PolynomialQ::operator-=(const PolynomialQ & rhs)
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
PolynomialQ & PolynomialQ::operator*=(const PolynomialQ & rhs)
{
	assert(Invariants());
	assert(rhs.Invariants());

	polynomial = expand(polynomial * rhs.polynomial); // need to figure out 'expand'

	assert(Invariants());

	return *this;
}

PolynomialQ & PolynomialQ::operator/=(const PolynomialQ & rhs)
{
    assert(Invariants());
    assert(rhs.Invariants());

    polynomial = quo(polynomial, rhs.polynomial, variable);

    assert(Invariants());

    return *this;
}

PolynomialQ operator+(const PolynomialQ & lhs, const PolynomialQ & rhs)
{
	assert(lhs.Invariants());
	assert(rhs.Invariants());

	return PolynomialQ(lhs.polynomial + rhs.polynomial);
}

PolynomialQ operator-(const PolynomialQ & lhs, const PolynomialQ & rhs)
{
	assert(lhs.Invariants());
	assert(rhs.Invariants());

	return PolynomialQ(lhs.polynomial - rhs.polynomial);
}

PolynomialQ operator*(const PolynomialQ & lhs, const PolynomialQ & rhs)
{
	assert(lhs.Invariants());
	assert(rhs.Invariants());

    // need to figure out 'expand'
	return PolynomialQ(expand(lhs.polynomial * rhs.polynomial));
}

PolynomialQ operator/(const PolynomialQ & lhs, const PolynomialQ & rhs)
{
    assert(lhs.Invariants());
    assert(rhs.Invariants());

    return PolynomialQ(quo(lhs.polynomial,
                           rhs.polynomial,
                           PolynomialQ::variable));
}

PolynomialQ operator%(const PolynomialQ & lhs, const PolynomialQ & rhs)
{
    assert(lhs.Invariants());
    assert(rhs.Invariants());

    return PolynomialQ(rem(lhs.polynomial,
                           rhs.polynomial,
                           PolynomialQ::variable));
}

PolynomialQ operator*(const PolynomialQ & lhs, const GiNaC::numeric & num)
{
    assert(lhs.Invariants());
    assert(num.is_rational());

    return lhs.polynomial*(GiNaC::ex)num; // will convert
}

PolynomialQ operator/(const PolynomialQ & lhs, const GiNaC::numeric & num)
{
    assert(lhs.Invariants());
    assert(num.is_rational());

    return lhs.polynomial/(GiNaC::ex)num; // rely on GiNaC throwing an exception
    									  // if num == 0.
}
