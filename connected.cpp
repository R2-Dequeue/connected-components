#include <cassert>

#include <ginac/ginac.h>
#include <boost/numeric/interval.hpp>

#include <numeric>
#include <functional>

// sres
// sresultant
// subresultant
// SubResultant

/*!
 * \brief Returns the kth subresultant of f and g with respect to the variable
 *        x, including constants like 1 and 0.
 *
 * \param f Any univariate polynomial in x.
 * \param g Same as f.
 * \param k 0 <= k <= min(deg(f), deg(g)).
 */
ex SubResultant(const ex & f, const ex & g, const int k, const symbol & x)
{
// local m,n,i,j,Mf,Mfl,Mg,Mgl,M;
// m   := degree(f,v);
// n   := degree(g,v);
// if m=k and n=k then return g fi;

// Mf  := Matrix(n-k,n-k+m-k-1, (i,j) -> coeff(f,v,m-j+i));
// Mfl := Vector(n-k, i -> v^(n-k-i) * f);
// Mg  := Matrix(m-k,n-k+m-k-1, (i,j) -> coeff(g,v,n-j+i));
// Mgl := Vector(m-k, i -> v^(m-k-i) * g);

// M   := <<Mf|Mfl>,<Mg|Mgl>>;
// return Determinant(M);

    assert(k >= 0);

    assert(f.is_polynomial(x));
    assert(g.is_polynomial(x));

    int m = f.degree(x);
    int n = g.degree(x);

    assert(m >= 0 && n >= 0);
    assert(k <= min(m, n));)

    if (m == k && n == k)
        return ex(g);

    matrix M(n+m-2*k, n+m-2*k);

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
}
