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
#include <memory>

#include <boost/foreach.hpp>

#include "polynomialbase.hpp"
#include "algebraic.hpp"

// Field to use
// Number of variables
//
// Complex      C,i
// Reals        R
// Algebraic    A
// Rationals    Q,F,R
// Integers     Z,I

// 0
// 1
// -1
// 11
// -11
// x - 1
// -x - 1
// 2*x + 2
// x^2 + 1
// 2*x^2 + 1
// x^2 - 5*x + 6

const GiNaC::symbol & PolynomialQ::variable = PolynomialBase::var1;

/*!
 * \detail a must be less than b.
 * \todo Fix old code and uncomment.
 */
unsigned int PolynomialQ::sturm(const GiNaC::numeric & a,
                                const GiNaC::numeric & b) const
{
    PolynomialQ::vector polynomials;
    polynomials.reserve(this->degree()+2);
    this->sturmseq(polynomials);
    return PolynomialQ::sturm(polynomials, a, b);

    /*PolynomialQ p2(*this);
    PolynomialQ p1(this->getDerivative());

    GiNaC::numeric  f2 = p2.eval(a),    g2 = p2.eval(b);
    GiNaC::numeric  fa = p1.eval(a),    gb = p1.eval(b);
    unsigned int    va = 0,             vb = 0;

    if (f2 != 0) // if == 0 then do nothing
    {
    	if (fa == 0)
    		fa = f2;
    	else
    		if (GiNaC::csgn(fa) != GiNaC::csgn(f2))
    			va = 1;
    }

    if (g2 != 0) // if == 0 then do nothing
    {
    	if (gb == 0)
    		gb = g2;
    	else
    		if (GiNaC::csgn(gb) != GiNaC::csgn(g2))
    			vb = 1;
    }

    while (!p1.isZero())
    {
    	PolynomialQ rem = (p1 % p2)*(-1);
    	GiNaC::numeric alpha = rem.eval(a), beta = rem.eval(b);

    	if (alpha != 0)
    	{
    		if (fa != 0)
    			if (GiNaC::csgn(alpha) != GiNaC::csgn(fa))
    				va++;
    		fa = alpha;
    	}

    	if (beta != 0)
    	{
    		if (gb != 0)
    			if (GiNaC::csgn(beta) != GiNaC::csgn(gb))
    				vb++;
    		gb = beta;
    	}

    	p2 = p1;
    	p1 = rem;
    }

    assert(int(va)-int(vb) >= 0);

    return (va - vb);*/
}

int PolynomialQ::signAt(const Algebraic & a) const
{
    assert(Invariants());

    PolynomialQ remainder(*this % a.polynomial);

    if (remainder.isZero())
        return 0;

    Algebraic alpha(a);

    while (true)
    {
        IntervalQ Y = remainder.boundRange(alpha.getInterval());

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

    //if (!value.is_rational())
    //	throw std::invalid_argument("PolynomialQ::eval: value passed is not"
    //								"rational.");

    GiNaC::ex temp(polynomial.subs(variable == value));

    assert(GiNaC::is_a<GiNaC::numeric>(temp));

    return GiNaC::numeric(GiNaC::ex_to<GiNaC::numeric>(temp));
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
        if (f.degree() <= 0) // just in case
            continue;
        if (f.degree() == 1)
        {
        	const GiNaC::numeric root = -f.getMonic().getCoeff(0);

            Algebraic alpha(Algebraic::MakeWideRational(root));

            numberSet.insert(alpha);
        }
        else
        {
            //std::vector<Algebraic> nums = PolynomialQ::FindRootsOfIrreducible(f);
            std::auto_ptr<Algebraic::vector> ptr =
                f.getRootsOfIrreducible<Algebraic::vector>();

            numberSet.insert(ptr->begin(), ptr->end());
        }
    }

    std::vector<Algebraic> numbers(numberSet.begin(), numberSet.end());

    if (numbers.size() > 1)
        for (std::vector<Algebraic>::size_type i = 0; i < numbers.size()-1; i++)
            Algebraic::SeparateIntervals(numbers[i], numbers[i+1]);

    return numbers;
}

/*
 * \param p Must be irreducible.
 * \todo Make this static.
 * \todo Does irreducible include 0? 1? any constant?
 */
//std::vector<IntervalQ> PolynomialQ::BoundRealRoots(const PolynomialQ & p)
//std::vector<Algebraic>
//    PolynomialQ::FindRootsOfIrreducible(const PolynomialQ & p);

/*!
 * \param p Any polynomial (including zero and other constants).
 */
std::auto_ptr<PolynomialQ::vector> PolynomialQ::sturmseq(const PolynomialQ & p)
{
    std::auto_ptr<PolynomialQ::vector> ptr(new PolynomialQ::vector);
    ptr->reserve(p.degree()+2);
    ptr->push_back(p);
    ptr->push_back(p.getDerivative());

    while (!ptr->back().isZero())
        ptr->push_back((ptr->at(ptr->size()-2) % ptr->back())*(-1));
    ptr->pop_back();

    return ptr;
}

/*!
 * \param p Any polynomial (including zero and other constants).
 */
PolynomialQ::vector & PolynomialQ::sturmseq(PolynomialQ::vector & polys) const
{
    //polys.reserve(p.degree()+2);
    polys.push_back(*this);
    polys.push_back(this->getDerivative());

    while (!polys.back().isZero())
        polys.push_back((polys[polys.size()-2] % polys.back())*(-1));
    polys.pop_back();

    return polys;
}

/*!
 * \detail a < b.
 */
unsigned int PolynomialQ::sturm(const std::vector<PolynomialQ> & F,
                                const GiNaC::numeric & a,
                                const GiNaC::numeric & b)
{
    assert(a <= b);

    if (F.size() == 1)
        return 0;

    GiNaC::numeric f2 = F[0].eval(a), g2 = F[0].eval(b);
    GiNaC::numeric fa = F[1].eval(a), gb = F[1].eval(b);
    unsigned int   va = 0,            vb = 0;

    if (f2 != 0) // if == 0 then do nothing
    {
    	if (fa == 0)
    		fa = f2;
    	else
    		if (fa.csgn() != f2.csgn())
    			va = 1;
    }

    if (g2 != 0) // if == 0 then do nothing
    {
    	if (gb == 0)
    		gb = g2;
    	else
    		if (gb.csgn() != g2.csgn())
    			vb = 1;
    }

    for (unsigned int i = 2; i < F.size(); i++)
    {
    	GiNaC::numeric alpha = F[i].eval(a), beta = F[i].eval(b);

    	if (alpha != 0)
    	{
    		if (fa != 0)
    			if (alpha.csgn() != fa.csgn())
    				++va;
    		fa = alpha;
    	}

    	if (beta != 0)
    	{
    		if (gb != 0)
    			if (beta.csgn() != gb.csgn())
    				++vb;
    		gb = beta;
    	}
    }

    assert(int(va) - int(vb) >= 0);

    return (va - vb);
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

    PolynomialQ::vector temp;
    fp.addIrreducibleFactorsTo(temp);

    return temp;
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

/*!
 * \detail Returns a bounding interval for the image of 'interval' in f. Uses
 *         Horner's scheme.
 * \todo Does the polynomial need to be monic?
 */
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
    assert(range.lower() <= range.upper());

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

    GiNaC::parser reader(table, true);

    PolynomialQ p(reader(s).expand()); // throws an exception if parsing fails

    return p;
}

void PolynomialQ::out()
{
    std::cout << *this;
}
#include <sstream>
std::string PolynomialQ::toString()
{
    std::stringstream s;
    s << *this;
    return s.str();
}
