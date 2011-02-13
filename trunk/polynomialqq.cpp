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

/*!
 * \throws parse_error Thrown by GiNaC if parsing fails (inherits
 *		   invalid_argument).
 */
PolynomialQQ::PolynomialQQ(const std::string & s)
{
	polynomial = PolynomialQQ::ParseString(s).polynomial;

    assert(Invariants());
}

/*!
 * \throws parse_error Thrown by GiNaC if parsing fails (inherits
 *		   invalid_argument).
 */
PolynomialQQ::PolynomialQQ(const char * const a)
{
	std::string s(a);

    polynomial = PolynomialQQ::ParseString(s).polynomial;

    assert(Invariants());
}

PolynomialQQ::PolynomialQQ(const GiNaC::ex & e)
    : polynomial(e)
{
    if (!Invariants())
        throw std::invalid_argument("PolynomialQQ Constructor: GiNaC expression "
        							"e is not a valid polynomial.");
}

PolynomialQQ::PolynomialQQ(const GiNaC::numeric & n)
    : polynomial(n)
{
    if (!Invariants())
        throw std::invalid_argument("PolynomialQQ Constructor: GiNaC numeric "
        							"n is not a valid rational number.");
}
/*
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
*/
bool PolynomialQQ::isZero() const
{
	assert(Invariants());

	return (polynomial.is_zero()); // (polynomial == 0);
}
/*
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

PolynomialQ PolynomialQQ::subx(const GiNaC::numeric & a) const
{
	assert(Invariants());
	assert(a.is_rational()); // throw?

	PolynomialQ p(polynomial.subs(var1 == a).subs(var2 == var1).expand());

	assert(p.Invariants());

	return p;
}

PolynomialQ PolynomialQQ::suby(const GiNaC::numeric & b) const
{
	assert(Invariants());
	assert(b.is_rational()); // throw?

	PolynomialQ p(polynomial.subs(var2 == b).expand());

	assert(p.Invariants());

	return p;
}

std::vector<PolynomialQQ> PolynomialQQ::getIrreducibleFactors() const
{
    assert(Invariants());

    std::vector<PolynomialQQ> factors; // Make a comparison function and change
                                       // this to use a set.

    GiNaC::ex p = GiNaC::factor(this->polynomial);

    if (GiNaC::is_a<GiNaC::power>(p))
    {
        PolynomialQQ temp(p.op(0));
        factors.push_back(temp);
    }
    else if (!GiNaC::is_a<GiNaC::mul>(p)) // what about '(x-1)^3' ?
    {
        if (!GiNaC::is_a<GiNaC::numeric>(p))
            factors.push_back(PolynomialQQ(p)); // (this->getMonic());
        return factors;
    }

    for (GiNaC::const_iterator i = p.begin(); i != p.end(); i++)
    {
        if (GiNaC::is_a<GiNaC::numeric>(*i))
            continue;
        else if (GiNaC::is_a<GiNaC::power>(*i))
        {
        	PolynomialQQ temp((*i).op(0));

        	factors.push_back(temp); // (p.getMonic());
        }
        else
        {
        	PolynomialQQ temp(*i);

        	factors.push_back(temp); // (p.getMonic());
        }
	}

    return factors;
}

int PolynomialQQ::signAt(const Algebraic & alpha, const Algebraic & beta) const
{
    assert(Invariants());

    boost::tuple<Algebraic, PolynomialQ, PolynomialQ> t =
        PolynomialQQ::Simple(alpha, beta);
    Algebraic & gamma   = t.get<0>();
    PolynomialQ & S     = t.get<1>();
    PolynomialQ & T     = t.get<2>();
std::cout << gamma << "    ";
std::cout << "( " << S << ", " << T << " )" << std::endl;
    PolynomialQ h(polynomial.subs(GiNaC::lst(var1 == S.getEx(),
                                             var2 == T.getEx())).expand());
std::cout << h << std::endl;
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
                                 PolynomialQQ(GiNaC::ex(1)),
                                 std::multiplies<PolynomialQQ>());

    //assert(!fp.isZero());

    return fp.getIrreducibleFactors();
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

    return PolynomialQ(res);
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

    //if (polynomial == ex()) // should always be at least 'ex(0)'.
    //    return false;

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

/*!
 * \todo Does not currently support single-point intervals in algebraic
 *       numbers.
 */
Algebraic PolynomialQQ::ANComb(Algebraic alpha,
                               Algebraic beta,
                               const GiNaC::numeric & t)
{
/**///std::cout << "**** ANComb **** A *******************" << std::endl;
/**///std::cout << "Alpha: " << alpha << std::endl;
/**///std::cout << "Beta: " << beta << std::endl;
/**///std::cout << "t: " << GiNaC::ex(t) << std::endl;
    GiNaC::ex polya = alpha.getEx();
    GiNaC::symbol tmp, var(alpha.getPolynomial().getVariable());;

    GiNaC::ex res = resultant(polya.subs(var == tmp - var*t).expand(),
                              beta.getEx(), var);
/**///std::cout << "res: " << res << std::endl;
    PolynomialQ r(res.subs(tmp == var)); // throws on error
/**///std::cout << "r: " << r << std::endl;
    std::vector<Algebraic> gammas =
    	PolynomialQ::FindRoots(r.getIrreducibleFactors());
/**///std::cout << "Size of gammas: " << gammas.size() << std::endl;
/**///BOOST_FOREACH(const Algebraic & gamma, gammas)
/**///    std::cout << "    " << gamma << std::endl;
    while (true)
    {
        IntervalQ IJ(alpha.lower() + t*beta.lower(),
                     alpha.upper() + t*beta.upper());

        BOOST_FOREACH(const Algebraic & gamma, gammas)
        {
            IntervalQ K = gamma.getInterval();

            if (K.lower() <= IJ.lower() && IJ.upper() <= K.upper())
            {
/**///std::cout << "**** ANComb **** B *******************" << std::endl;
                return gamma;
            }
        }

        alpha.tightenInterval();
        beta.tightenInterval();
    }
/*A := alpha[2]; // A = [ [l, u], poly ]
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
    assert(false);
    throw std::runtime_error("ANComb error. This point should not have been"
    						 "reached");

    return Algebraic();
}

boost::tuple<Algebraic, PolynomialQ, PolynomialQ>
    PolynomialQQ::Simple(const Algebraic & alpha, const Algebraic & beta)
{
    boost::tuple<Algebraic, PolynomialQ, PolynomialQ> t =
        Simple2(beta, alpha);

    Algebraic &     gamma   = t.get<0>();
    PolynomialQ &   S       = t.get<1>();
    PolynomialQ &   T       = t.get<2>();

    return boost::make_tuple(gamma, T % gamma.getPolynomial(), S);
// Simple(alpha, beta)
//      gamma, S,T := Simple2(beta,alpha);
//      return gamma, T % gamma.getPolynomial(), S;
}

void ANCombCheck(const Algebraic & gamma,
                 const Algebraic & alpha,
                 const Algebraic & beta,
                 const GiNaC::numeric & t)
{
    std::cout << GiNaC::ex(gamma.Approximate()) << " ?= "
              << GiNaC::ex(alpha.Approximate() + beta.Approximate()*t)
              << std::endl;
}

void gcdexCheck(const GiNaC::ex & gcd, const GiNaC::ex & f, const GiNaC::ex g,
                const GiNaC::symbol & var, const GiNaC::ex & c1,
                const GiNaC::ex & c2)
{
    std::cout << "** " << gcd << " ?= " << (f*c1 + g*c2).expand() << std::endl;
}

boost::tuple<Algebraic, PolynomialQ, PolynomialQ>
    PolynomialQQ::Simple2(const Algebraic & alpha, const Algebraic & beta)
{
    GiNaC::ex A(alpha.getEx());
    GiNaC::ex B;

    GiNaC::symbol u = alpha.getPolynomial().getVariable();
    GiNaC::symbol v; // makes a new, unique symbol
    GiNaC::symbol w = u; // The return variable.

    B = beta.getEx().subs(u == v);

    GiNaC::numeric t(1);

    GiNaC::symbol _z; // The temporary variable used.
    GiNaC::ex C;
    GiNaC::ex s1;
    GiNaC::ex c1, c2;
    GiNaC::ex g;

    Algebraic gamma = Algebraic();

    while (true)
    {
        gamma = PolynomialQQ::ANComb(alpha, beta, t);
        //ANCombCheck(gamma, alpha, beta, t);
        C = gamma.getEx().subs(w == _z);
        s1 = PolynomialQQ::sres(A.subs(u == _z - t*v).expand(), B, 1, v);
        g = PolynomialQQ::gcdex(C, s1.coeff(v, 1), _z, c1, c2);
        //gcdexCheck(g,C,s1.coeff(v,1),_z,c1,c2);

        if (g.degree(_z) == 0)
            break;

        t = t+1;
    }

    GiNaC::ex T = GiNaC::rem((-s1.coeff(v, 0)) * c2 / g, C, _z);
    GiNaC::ex S = _z - T*t;

    S = S.subs(_z == w);
    T = T.subs(_z == w);

    return boost::make_tuple(gamma, PolynomialQ(S), PolynomialQ(T));

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

GiNaC::ex PolynomialQQ::gcdex(GiNaC::ex f,
							  GiNaC::ex g,
							  const GiNaC::symbol & var,
							  GiNaC::ex & c1,
							  GiNaC::ex & c2)
{
	GiNaC::ex x = 0;	c1 = 1;
	GiNaC::ex y = 1;	c2 = 0;
	GiNaC::ex t1, t2, quotient;

	while (g != 0)
	{
		quotient = GiNaC::quo(f, g, var);

		t1 = f;
		t2 = g;
		f = t2;
		g = GiNaC::rem(t1, t2, var);

		t1 = x;
		t2 = c1;
		x = t2 - quotient*t1;
		c1 = t1;

		t1 = y;
		t2 = c2;
		y = t2 - quotient*t1;
		c2 = t1;
	}

	GiNaC::ex lead = f.lcoeff(var);

	f /= lead;
	c1 /= lead;
	c2 /= lead;

	return f;

//function extended_gcd(a, b)
//    x := 0    lastx := 1
//    y := 1    lasty := 0
//    while b != 0
//        quotient := a div b
//
//        {a, b} = {b, a mod b}
//        {x, lastx} = {lastx - quotient*x, x}
//        {y, lasty} = {lasty - quotient*y, y}
//    return {lastx, lasty, a}
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

            M(i-1, j-1) = f.coeff(var, m-j+i);

    for (       int i = n-k+1;  i <= n+m-2*k;       i++)
        for (   int j = 1;      j <= n+m-2*k-1;     j++)

            M(i-1, j-1) = g.coeff(var, n-j+(i-(n-k)));

    for (int i = 1;         i <= n-k;       i++)
        M(i-1, n+m-2*k-1) = pow(var, n-k-i)             *f;

    for (int i = n-k+1;     i <= n+m-2*k;   i++)
        M(i-1, n+m-2*k-1) = pow(var, m-k-(i-(n-k)))     *g;

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

    PolynomialQQ p(reader(s).expand()); // throws an exception if parsing fails.

    return p;
}

std::ostream & operator<<(std::ostream & output, const PolynomialQQ & p)
{
    output << p.getEx();

    return output;
}
