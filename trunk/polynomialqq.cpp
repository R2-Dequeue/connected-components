/*!
 * \file
 * \author Chris de Pujo
 */

#include "polynomialqq.hpp"

#include "polynomialbase.hpp"
#include "algebraic.hpp"

#include <stdexcept>
#include <numeric>
#include <functional>

#include <boost/foreach.hpp>

//#include <boost/numeric/ublas/matrix.hpp>

const GiNaC::symbol & PolynomialQQ::var1 = PolynomialBase::var1;
const GiNaC::symbol & PolynomialQQ::var2 = PolynomialBase::var2;

PolynomialQQ::PolynomialQQ(const std::string & s)
{
	polynomial = PolynomialQQ::ParseString(s);
    
    assert(Invariants());
}

PolynomialQQ::PolynomialQQ(const char * const a)
{
	std::string s(a);
	
    polynomial = PolynomialQQ::ParseString(s);
    
    assert(Invariants());
}

PolynomialQQ::PolynomialQQ(const GiNaC::ex & e)
    : polynomial(e)
{
    if (!Invariants())
        throw std::invalid_argument("PolynomialQQ Constructor: GiNaC expression"
        							"e is not a valid polynomial.");
}

PolynomialQQ::PolynomialQQ(const GiNaC::numeric & n)
    : polynomial(n)
{
    if (!Invariants())
        throw std::invalid_argument("PolynomialQQ Constructor: GiNaC numeric"
        							"n is not a valid rational number.");
}

inline int PolynomialQQ::degree() const
{
	assert(Invariants());

	return polynomial.degree(variable);
}

inline bool PolynomialQQ::isMonic() const
{
	assert(Invariants());

	return (polynomial.lcoeff(variable) == 1);
}

inline bool PolynomialQQ::isZero() const
{
	assert(Invariants());

	return (polynomial.is_zero()); // (polynomial == 0);
}

bool PolynomialQQ::isIrreducible() const
{
	assert(Invariants());
	
	if (this->degree() > 2)
		return false;
	
	// Assuming getIrreducibleFactors ignores constant factors.
	if (this->getIrreducibleFactors().size() > 1)
		return false;
	
	return true;
}

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

/*!
 * \detail Should work in cases such as 'p += p;'
 */
PolynomialQQ & PolynomialQQ::operator+=(const PolynomialQQ & rhs)
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
PolynomialQQ & PolynomialQQ::operator-=(const PolynomialQQ & rhs)
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
PolynomialQQ & PolynomialQQ::operator*=(const PolynomialQQ & rhs)
{
	assert(Invariants());
	assert(rhs.Invariants());

	polynomial = expand(polynomial * rhs.polynomial); // need to figure out 'expand'

	assert(Invariants());

	return *this;
}

PolynomialQQ PolynomialQQ::getDerivative(unsigned int variable) const
{
    assert(Invariants());
    assert(variable == 1 || variable == 2);

    if (variable != 1 && variable != 2)
        throw std::invalid_argument("Tried to take derivative wrt an invalid variable.");

    PolynomialQQ temp(*this);

    temp.differentiate(variable);

    return temp;
}

PolynomialQQ & PolynomialQQ::differentiate(unsigned int variable)
{
    assert(Invariants());
    assert(variable == 1 || variable == 2);

    if (variable != 1 && variable != 2)
        throw std::invalid_argument("Tried to take derivative wrt an invalid variable.");

    polynomial = polynomial.diff((variable == 1) ? var1 : var2);

    return *this;
}

std::vector<PolynomialQQ> PolynomialQQ::getIrreducibleFactors() const
{
    assert(Invariants());

    GiNaC::ex p = factor(polynomial);

    std::vector<PolynomialQQ> factors;

    for (GiNaC::const_iterator i = p.begin(); i != p.end(); i++)
        factors.push_back(PolynomialQQ(GiNaC::is_a<GiNaC::power>(*i) ?
                                      (*i).op(0) :
                                      *i));
        // The PolyQ constructor will throw if the ex is not a valid poly.

    return factors;
}

int PolynomialQQ::signAt(const Algebraic & alpha, const Algebraic & beta) const
{
    assert(Invariants());

    boost::tuple<Algebraic, PolynomialQ, PolynomialQ> t =
        PolynomialQQ::Simple(alpha, beta);
    Algebraic & gamma   = t.get<0>();
    PolynomialQ & s1    = t.get<1>();
    PolynomialQ & s2    = t.get<2>();

    PolynomialQ h(polynomial.subs(
                        GiNaC::lst(var1 == s1.getEx(), var2 == s2.getEx())));

    return h.signAt(gamma);

/*
 local gamma, Ss, _z, h;
 gamma, Ss := Simple(alphas,_z);
 h := eval(f,{seq(Var(alphas[i])=Ss[i],i=1..nops(Ss))});
 return Sign1(gamma, h, _z);
*/
}

inline GiNaC::ex PolynomialQQ::getEx() const
{
    GiNaC::ex temp(polynomial);

    return temp;
}

std::vector<PolynomialQQ>
    PolynomialQQ::IrreducibleFactors(const std::vector<PolynomialQQ> & F)
{
    // remove '0' polynomials?
    PolynomialQQ fp = std::accumulate(F.begin(),
                                 F.end(),
                                 PolynomialQQ(GiNaC::numeric(1)),
                                 std::multiplies<PolynomialQQ>());

    assert(fp.Invariants()); // isn't this redundant?

    GiNaC::ex p = factor(fp.polynomial);

    std::vector<PolynomialQQ> factors;

    for (GiNaC::const_iterator i = p.begin(); i != p.end(); i++)
        factors.push_back(PolynomialQQ(GiNaC::is_a<GiNaC::power>(*i) ?
                                      (*i).op(0) :
                                      *i));
        // The PolyQ constructor will throw if the ex is not a valid poly.

    return factors;
}

PolynomialQ PolynomialQQ::Resultant(const PolynomialQQ & f,
                                    const PolynomialQQ & g,
                                    unsigned int var)
{
    assert(f.Invariants());
    assert(g.Invariants());
    assert(var == 1 || var == 2);

    if (var != 1 && var != 2)
        throw std::invalid_argument("Tried to calculate Resultant wrt an invalid variable.");

    GiNaC::ex res = resultant(f.polynomial,
                              g.polynomial,
                              (var == 1) ? PolynomialQQ::var1 : PolynomialQQ::var2);

    if (var == 1)
        res = res.subs(var2 == var1);

    return PolynomialQ(res); // Hoping constructor will be called and check
                             // Invariants().
}

/*!
 * \detail Returns the kth subresultant of f and g with respect to the variable
 *         x, including constants like 1 and 0.
 *
 * \param f Any univariate polynomial in x.
 * \param g Same as f.
 * \param k 0 <= k <= min(deg(f), deg(g)).
 *//*
PolynomialQ Subresultant(const PolynomialQQ & f,
                         const PolynomialQQ & g,
                         const unsigned int k,
                         const unsigned int var)
{
    int m = f.degree();
    int n = g.degree();

    assert(f.Invariants());
    assert(g.Invariants());

    if (k > min(m, n))
        throw std::invalid_argument("Subresultant: k is out of bounds.");

    if (var != 1 && var != 2)
        throw std::invalid_argument("Tried to calculate Subresultant wrt an invalid variable.");

    if (m == k && n == k)
        return ex(g);

    boost::numeric::ublas::matrix<GiNaC::numeric> M(n+m-2*k, n+m-2*k);

    for (int i = 0; i < n-k; i++)
    {
        for (int j = i; j <= m+i; j++)
            M(i, j) = f.coeff(x, m-(j-i));

        M(i, n+m-2*k-1) = f*pow(x, n-k-(i+1));
    }

    for (int i = n-k+1; i < n+m-2*k; i++)
    {
        for (int j = i-(n-k+1); j <= n+(i-(n-k+1)); j++)
            M(i, j) = g.coeff(x, n-(j-i));

        M(i-(n-k+1), n+m-2*k-1) = f*pow(x, n-k-(i+1));
    }

    return M.determinant();
}*/

bool PolynomialQQ::Invariants() const
{
    using namespace GiNaC;

    if (this == NULL)
        return false;

    if (!polynomial.is_polynomial(lst(var1, var2)))
        return false;

    if (!polynomial.info(info_flags::rational_polynomial))
        return false;

    if (polynomial == ex()) // should always be at least 'ex(0)'.
        return false;

    int deg1 = polynomial.degree(var1);

    for (int i = 0; i <= deg1; i++)
    {
        ex p2 = polynomial.coeff(var1, i);

        if (!p2.is_polynomial(var2))
            return false;

        if (!p2.info(info_flags::rational_polynomial))
            return false;

        int deg2 = p2.degree(var2);

        for (int j = 0; j <= deg2; j++)
        {
            ex co = p2.coeff(var2, j);

            if (is_a<numeric>(co))
            {
                numeric c = ex_to<numeric>(co);

                if (!c.is_rational())
                    return false;
            }
            else
                return false;
        }
    }

    return true;
}

PolynomialQQ PolynomialQQ::operator+(const PolynomialQQ & rhs) const {
	assert(Invariants() && rhs.Invariants());
	return PolynomialQQ(polynomial + rhs.polynomial);
}
PolynomialQQ PolynomialQQ::operator-(const PolynomialQQ & rhs) const {
	assert(Invariants() && rhs.Invariants());
	return PolynomialQQ(polynomial - rhs.polynomial);
}
PolynomialQQ PolynomialQQ::operator*(const PolynomialQQ & rhs) const {
	assert(Invariants() && rhs.Invariants());
	return PolynomialQQ(expand(polynomial * rhs.polynomial));
}

Algebraic
    PolynomialQQ::ANComb(const Algebraic & alpha, const Algebraic & beta, int t)
{
    Algebraic a(alpha), b(beta);
    GiNaC::ex polya = a.getEx();
    GiNaC::symbol tmp, var(a.getPolynomial().getVariable());;

    GiNaC::ex res = resultant(polya.subs(var == tmp - var*t),
                              b.getEx(), var);
    res = res.subs(tmp == var);

    PolynomialQ r(res); // throws on error

    std::vector<Algebraic> gammas =
    	PolynomialQ::FindRoots(r.getIrreducibleFactors());

    // maybe declare a & b here

    while (true)
    {
        IntervalQ IJ = a.getInterval() + b.getInterval()*GiNaC::numeric(t);

        BOOST_FOREACH(const Algebraic & gamma, gammas)
        {
            IntervalQ K = gamma.getInterval();

            if (K.lower() <= IJ.lower() && IJ.upper() <= K.upper())
                return gamma;
        }

        a.tightenInterval();
        b.tightenInterval();
    }

    assert(false);
    throw std::runtime_error("ANComb error. This point should not have been"
    						 "reached");

    return Algebraic();

/*
 A := alpha[2]; // A = [ [l, u], poly ]
 B := beta[2];
 u := Var(A);
 v := Var(B);
 r := resultant(eval(A,u=_z-t*v),B,v);
 r := eval(r,_z=w);
 gammas := FindRoots(IrreducibleFactors([r]));
 alphap := alpha;
 betap  := beta;
 while true do
   IJ := alphap[1] + t * betap[1];
   for gamma in gammas do
     K := gamma[1];
     if K[1] <= IJ[1] and IJ[2] <= K[2] then
       return gamma;
     fi;
   od;
   alphap := TightenInterval(alphap);
   betap  := TightenInterval(betap);
 od; */
}

boost::tuple<Algebraic, PolynomialQ, PolynomialQ>
    PolynomialQQ::Simple(const Algebraic & alpha, const Algebraic & beta)
{
    GiNaC::symbol var = alpha.getPolynomial().getVariable();
    GiNaC::symbol tmp;
    GiNaC::ex s1;
    GiNaC::ex g;
    Algebraic gamma;
    int t = 1;

    while (true)
    {
        gamma = PolynomialQQ::ANComb(alpha, beta, t);
        s1 = PolynomialQQ::sres(
                            alpha.getEx().subs(var == tmp - var*t),
                            beta.getEx(),
                            1,
                            var);

        // s1 is now a polynomial in 'tmp'.

        g = gcdex(gamma.getEx().subs(gamma.getVariable() ==, s1.coeff(tmp, 1))

        if (g.degree() == 0)
            break;

        t++;
    }

    PolynomialQ T(s1.coeff(alpha.getPolynomial().getVariable(), 0));
    PolynomialQ S(T.getVariable() - T.getEx()*t);

    return boost::make_tuple(gamma,S, T);

/* A := alpha[2];
 B := beta[2];
 u := Var(A);
 v := Var(B);
 t := 1;
 while true do
   gamma := ANComb(alpha,beta,_z,t);
   C     := gamma[2];
   s1    := Sres(eval(A,u=_z-t*v),B,1,v);
   g     := gcdex(C, coeff(s1,v,1), _z, 'c1','c2');
   if degree(g,_z) = 0 then
       break
   fi;
   t := t + 1;
 od;
 T  := rem((-coeff(s1,v,0)) * c2 / g, C, _z);
 S  := _z - t*T;
 gamma := eval(gamma,_z=w);
 S     := eval(S,    _z=w);
 T     := eval(T,    _z=w);
 return gamma, S, T; */
}

GiNaC::ex PolynomialQQ::sres(const GiNaC::ex & f,
                              const GiNaC::ex & g,
                              const int k,
                              const GiNaC::symbol & var)
{
    const int m = f.degree(var);
    const int n = g.degree(var);

    if (m == k && n == k)
        return g;

    GiNaC::matrix M(n+m-2*k, n+m-2*k); // k == 1

    for (       int i = 1;      i <= n-k;           i++)
        for (   int j = 1;      j <= n+m-2*k-1;     j++)

            M(i, j) = f.coeff(var, m-j+i);

    for (       int i = n-k+1;  i <= n+m-2*k;       i++)
        for (   int j = 1;      j <= n+m-2*k-1;     j++)

            M(i, j) = g.coeff(var, n-j+(i-(n-k)));

    for (int i = 1;         i <= n-k;       i++)
        M(i, n+m-2*k) = pow(var, n-k-i)             *f;

    for (int i = n-k+1;     i <= n+m-2*k;   i++)
        M(i, n+m-2*k) = pow(var, m-k-(i-(n-k)))     *g;

    return M.determinant();

/*
 local m,n,i,j,Mf,Mfl,Mg,Mgl,M;
 m   := degree(f,v);
 n   := degree(g,v);
 if m=k and n=k then return g fi;
 Mf  := Matrix(n-k,n-k+m-k-1, (i,j) -> coeff(f,v,m-j+i));
 Mfl := Vector(n-k, i -> v^(n-k-i) * f);
 Mg  := Matrix(m-k,n-k+m-k-1, (i,j) -> coeff(g,v,n-j+i));
 Mgl := Vector(m-k, i -> v^(m-k-i) * g);
 M   := <<Mf|Mfl>,<Mg|Mgl>>;
 return Determinant(M);
*/
}

/*!
 * \detail Parses the string wrt the internal variables.
 * \throws parse_error Thrown by GiNaC if parsing fails (inherits
 *		   invalid_argument).
 */
inline PolynomialQQ PolynomialQQ::ParseString(const std::string & s) const
{
    GiNaC::symtab table;

    table[var1.get_name()] = var1;
    table[var2.get_name()] = var2;

    GiNaC::parser reader(table); // reader(table, true);
    reader.strict = true;

    // reader(s).expand(); ?
    PolynomialQQ p(reader(s)); // throws an exception if parsing fails.
    
    return p;
}
