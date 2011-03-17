/*!
 * \file
 * \author Chris de Pujo
 */

#include "polynomialqq.hpp"

#include <stdexcept>
#include <numeric>
#include <functional>
#include <sstream> // for getString()

#include <boost/foreach.hpp>

#include "polynomialbase.hpp"
#include "algebraic.hpp"
#include "templatehelp.hpp"

const GiNaC::symbol & PolynomialQQ::var1 = PolynomialBase::var1;
const GiNaC::symbol & PolynomialQQ::var2 = PolynomialBase::var2;

int PolynomialQQ::signAt(const Algebraic & alpha, const Algebraic & beta) const
{
    assert(Invariants());

    boost::tuple<Algebraic, PolynomialQ, PolynomialQ> t =
        PolynomialQQ::Simple(alpha, beta);

    Algebraic & gamma   = t.get<0>();
    PolynomialQ & S     = t.get<1>();
    PolynomialQ & T     = t.get<2>();

    PolynomialQ h(polynomial.subs(GiNaC::lst(var1 == S.getEx(),
                                             var2 == T.getEx())).expand());

    return h.signAt(gamma);

/*
 local gamma, Ss, _z, h;
 gamma, Ss := Simple(alphas,_z);
 h := eval(f,{seq(Var(alphas[i])=Ss[i],i=1..nops(Ss))});
 return Sign1(gamma, h, _z);
*/
}

PolynomialQQ::vector
    PolynomialQQ::IrreducibleFactors(const PolynomialQQ::vector & F)
{
    // remove '0' polynomials?
    //assert(!fp.isZero());

    /*PolynomialQQ::vector factors;
    BOOST_FOREACH(const PolynomialQQ & f, F)
        f.addIrreducibleFactorsTo(factors);*/

    /*PolynomialQQ::set factors;

    for (PolynomialQQ::vector::const_iterator i = F.begin(), e = F.end();
         i != e; ++i)
        i->addIrreducibleFactorsTo(factors);

    PolynomialQQ::vector temp(factors.begin(), factors.end());*/

    PolynomialQQ::vector factors, temp;
    factors.reserve(F.size());
    temp.reserve(F.size());

    for (PolynomialQQ::vector::const_iterator f = F.begin(), e1 = F.end();
         f != e1; ++f)
    {
        f->addIrreducibleFactorsTo(temp);

        for (PolynomialQQ::vector::iterator g = temp.begin(), e2 = temp.end();
             g != e2; ++g)
            if (find(factors.begin(), factors.end(), *g) == factors.end())
                factors.push_back(*g);

        temp.clear();
    }

    assert(tmp::IsUnique(factors));

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

/*!
 * \todo Does not currently support single-point intervals in algebraic
 *       numbers.
 */
Algebraic PolynomialQQ::ANComb(Algebraic alpha,
                               Algebraic beta,
                               const GiNaC::numeric & t)
{
    assert(t.is_nonneg_integer() && !t.is_zero());

    Algebraic::SeparateIntervals(alpha, beta);//, GiNaC::numeric(1, 128));

    GiNaC::ex A = alpha.getEx();
    GiNaC::symbol _z("_z"), v("v"), u = alpha.getPolynomial().getVariable();

    GiNaC::ex B = beta.getEx().subs(u == v);

    GiNaC::ex temp = A.subs(u == _z - v*t).expand();
    GiNaC::ex res = //PolynomialQQ::sres(temp, B, 0, v);
                    GiNaC::resultant(temp, B, v);

    PolynomialQ r(res.subs(_z == u)); // throws on error

    Algebraic::vector gammas; /**/ gammas.reserve(r.degree());
    r.addRootsTo(gammas);
    Algebraic::SeparateIntervals(gammas);//, GiNaC::numeric(1,128));

    while (true)
    {
        IntervalQ IJ(alpha.lower() + t*beta.lower(),
                     alpha.upper() + t*beta.upper());

        for (Algebraic::vector::const_iterator
             gamma = gammas.begin(), e = gammas.end(); gamma != e; ++gamma)
            if (gamma->lower() <= IJ.lower() && IJ.upper() <= gamma->upper())
            {
                GiNaC::numeric a = alpha.approx(10), b = beta.approx(10),
                               g = gamma->approx(10);
                if (GiNaC::abs((a + t*b) - g) >= GiNaC::numeric(1,1000))
                    std::cout << " ";
                return *gamma;
            }

        alpha.tightenInterval();
        beta.tightenInterval();
    }

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

bool validSimple2(Algebraic alpha, Algebraic beta,
                  Algebraic gamma, PolynomialQ S, PolynomialQ T, PolynomialQ n)
{
    GiNaC::numeric delta(1, 100);
    GiNaC::numeric g = gamma.approx(6);
    GiNaC::ex s = S.getEx().subs(S.getVariable() == g);
    GiNaC::ex t = T.getEx().subs(T.getVariable() == g);
    GiNaC::ex a = alpha.approx(6);
    GiNaC::ex b = beta.approx(6);

    if (GiNaC::abs(a-s) < delta && GiNaC::abs(b-t) < delta)
        return true;
    else
        return false;
}

static std::string Print(const GiNaC::ex & a)
{
    std::stringstream ss;
    ss << a;
    return ss.str();
}

static std::string Print(const GiNaC::numeric & a)
{
    std::stringstream ss;
    ss << (GiNaC::ex)a;
    return ss.str();
}

static std::string Print(const GiNaC::symbol & a)
{
    std::stringstream ss;
    ss << a;
    return ss.str();
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
        s1 = PolynomialQQ::sres(A.subs(u == _z - t*v).expand(), B, 1, v);//expand may not be necessary
        g = PolynomialQQ::gcdex(C, s1.coeff(v, 1), _z, c1, c2);
        //gcdexCheck(g,C,s1.coeff(v,1),_z,c1,c2);

        if (g.degree(_z) == 0)
            break;

        ++t;
    }

    GiNaC::ex T = GiNaC::rem(((-s1.coeff(v, 0)) * c2 / g).expand(), C, _z);
    GiNaC::ex S = (_z - (GiNaC::ex)t*T).expand();// expands may not be necessary

    S = S.subs(_z == w);
    T = T.subs(_z == w);

    if (!validSimple2(alpha, beta, gamma, PolynomialQ(S), PolynomialQ(T),t))
    {
        std::stringstream s;
        s << s1;
        int i = 0;
        ++i;
    }

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

/*GiNaC::ex PolynomialQQ::gcdex(GiNaC::ex f,
							  GiNaC::ex g,
							  const GiNaC::symbol & var,
							  GiNaC::ex & c1,
							  GiNaC::ex & c2)
{
    const GiNaC::ex _f(f), _g(g);

	GiNaC::ex x = 0;	c1 = 1;
	GiNaC::ex y = 1;	c2 = 0;
	GiNaC::ex t, quotient;

	while (g != 0)
	{
		quotient = GiNaC::quo(f, g, var);

		t = g;
		g = GiNaC::rem(f, g, var);
		f = t;

		t = x;
		x = c1 - quotient*x;
		c1 = t;

		t = y;
		y = c2 - quotient*y;
		c2 = t;
	}

	assert(f == (_f*c1 + _g*c2).expand());

	return f;
}*/

GiNaC::ex PolynomialQQ::gcdex(const GiNaC::ex & f,
							  const GiNaC::ex & g,
							  const GiNaC::symbol & var,
							  GiNaC::ex & c1,
							  GiNaC::ex & c2)
{
    GiNaC::ex D(0);

	if (g.is_zero())
	{
	    if (f.is_zero())
	    {
	        c1 = 0;
	        c2 = 0;
	        return D;
	    }

	    D = f;
	    GiNaC::ex lead = f.lcoeff(var);
	    D = (D/lead).expand();
	    c1 = GiNaC::ex(1) / lead;
	    c2 = 0;
	    return D;
	}

	GiNaC::ex q = GiNaC::quo(f,g,var), r = GiNaC::rem(f,g,var), ss, tt;
	D = PolynomialQQ::gcdex(g, r, var, ss, tt);
	c1 = tt;
	c2 = ss - q*tt;

	return D;
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

    GiNaC::matrix M(n+m-2*k, n+m-2*k);

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

std::string PolynomialQQ::getString() const
{
    std::stringstream s;
    s << *this;
    return s.str();
}
