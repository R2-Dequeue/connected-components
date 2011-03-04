/*!
 * \file
 * \author Chris de Pujo
 */

#include "cad.hpp"

#include <cassert>
#include <numeric>
//#include <functional>

#include <boost/foreach.hpp>
#include <boost/numeric/ublas/io.hpp>

/*!
 * \throws parse_error Thrown by GiNaC if parsing of polynomials fails (inherits
 *		   invalid_argument).
 * \todo Check if doxygen collates all 'throws' statements from called
 *		 functions.
 */
CAD::CAD(const std::vector<std::string> & F)
{
    // Be careful with 'F' and 'this->F'.

    for (std::vector<std::string>::const_iterator f = F.begin(), e = F.end();
         f != e; ++f)
        this->F.push_back(PolynomialQQ(*f)); // throws on error

    PolynomialQQ::vector temp;
    for (PolynomialQQ::vector::const_iterator f = this->F.begin(),
                                              e = this->F.end();
         f != e; ++f)
        f->addIrreducibleFactorsTo(temp);

    for (PolynomialQQ::vector::iterator i = temp.begin(), e = temp.end();
        i != e; ++i)
        if (find(this->irreducibles.begin(), this->irreducibles.end(), *i) ==
            this->irreducibles.end())
            this->irreducibles.push_back(*i);

    alphas = CAD::SamplePoints(PolynomialQ::FindRoots(CAD::Project(this->F)));

    stacks.reserve(alphas.size());

    for (Algebraic::vector::iterator alpha = alphas.begin(), e = alphas.end();
         alpha != e;
         ++alpha)
    {
        std::vector<Algebraic> betas =
            CAD::SamplePoints(CAD::FindRoots2(*alpha, this->F));

        stacks.push_back(std::vector<Sample>()); // Should I name this parameter?
        std::vector<Sample> & T = stacks.back();

        T.reserve(betas.size());

        BOOST_FOREACH(const Algebraic & beta, betas)
        {
            std::vector<char> signs(this->F.size()); // preallocate

            for (std::vector<char>::size_type i = 0, e = this->F.size();
                 i < e;
                 ++i)
                signs[i] = this->F[i].signAt(*alpha, beta);

            T.push_back(Sample(beta, signs));
        }
    }

    MakeConnectivityMatrix();

    Partition();

    // maybe this should be an assert?
    if (!Invariants())
        throw std::logic_error("CAD Constructor: CAD creation completed, but"
                               "the result is not canonical.");
}

bool CAD::Connectivity(const Point & p1, const Point & p2) const
{
    assert(Invariants());

    return (cMatrix(CellNumber(Cell(p1)),
                    CellNumber(Cell(p2))) == 1);
}

CellIndex CAD::Cell(const Point & p) const
{
    assert(Invariants());

    typedef unsigned int uint;

	uint ix = 2*alphas.size();

	for (uint i = 0; i < alphas.size(); ++i) // i from 1 to nops(alphas)
	{
		int t = p.first.compare(alphas[i]);

		if (t <= 0)
		{
			if (t == 0)
				ix = 2*i + 1;
			else
				ix = 2*i;

			break;
		}
	}

	std::vector<Algebraic> betas = CAD::FindRoots2(p.first, F);

	uint iy = 2*betas.size();

	for (uint i = 0; i < betas.size(); ++i) // i from 1 to nops(betas)
	{
		int t = p.second.compare(betas[i]);

		if (t <= 0)
		{
			if (t == 0)
				iy = 2*i + 1;
			else
				iy = 2*i;

			break;
		}
	}

	return CellIndex(ix, iy);
}

void MatrixOut(GiNaC::matrix & M)
{
    using std::cout;
    using std::endl;
    typedef unsigned int uint;

    cout << M << endl;
}

unsigned int CAD::CellNumber(const CellIndex & i) const
{
    // Indices start at 0, not 1.
    assert(Invariants());
    assert(i.first < stacks.size());
    assert(i.second < stacks[i.first].size());

    unsigned int num = 0;

    for (unsigned int j = 0; j < i.first; j++)
        num += stacks[j].size();

    num += i.second;

    return num;
}

void CAD::out() const
{
	assert(Invariants());

	using namespace std;

	typedef unsigned int uint;

	cout << "Number of functions: " << F.size() << endl;
	cout << "Functions: ";

	if (F.size() >= 1)
		cout << F[0] << endl;

	for (uint i = 1; i < F.size(); i++)
		cout << "           " << F[i] << endl;

	cout << "Number of alphas: " << alphas.size() << endl;
	cout << "Number of stacks: " << stacks.size() << endl;
	cout << "Individual stack sizes: " << endl;

	for (uint i = 0; i < stacks.size(); i++)
		cout << "    " << "Size of stack " << i << ": " << stacks[i].size() << endl;

	cout << "Samples: " << endl;

	for (uint i = 0; i < alphas.size(); i++)
		cout << "    " << alphas[i] << endl;

	cout << "Stacks: " << endl;

	for (uint i = 0; i < stacks.size(); i++)
	{
		cout << "    Stack " << i << ": " << endl;

		for (uint j = 0; j < stacks[i].size(); j++)
		{
			cout << "        " << stacks[i][j].y << " , Signs: ";

			for (uint k = 0; k < stacks[i][j].signs.size(); k++)
				cout << int(stacks[i][j].signs[k]) << " ";

			cout << endl;
		}
	}

	cout << "Connectivity matrix: " << endl;

	for (uint i = 0; i < cMatrix.rows(); ++i)
	{
        for (uint j = 0; j < cMatrix.cols(); ++j)
            cout << cMatrix(i, j) << " ";
        cout << endl;
	}

	cout << "Partitions: " << endl;

	for (uint i = 0; i < partitions.size(); ++i)
    {
        for (uint j = 0; j < partitions[i].size(); ++j)
            cout << partitions[i][j] << " ";
        cout << endl;
    }
}

bool CAD::Invariants() const
{
    if (this == NULL)
        return false;

    if (stacks.empty())
        return false;

    if (alphas.size() != stacks.size())
        return false;

    BOOST_FOREACH(const std::vector<Sample> & stack, stacks)
        if (stack.empty())
            return false;

    BOOST_FOREACH(const Algebraic & alpha, alphas)
        if (!alpha.Invariants())
            return false;

    BOOST_FOREACH(const std::vector<Sample> & stack, stacks)
        BOOST_FOREACH(const Sample & sample, stack)
        {
            if (!sample.y.Invariants())
                return false;

            BOOST_FOREACH(const char & sign, sample.signs)
                if (sign != -1 && sign !=  0 && sign !=  1)
                    return false;
        }

    return true;
}

void CAD::MakeConnectivityMatrix()
{
    typedef unsigned int uint;

	uint n = 0;

	BOOST_FOREACH(const std::vector<Sample> & s, stacks)
		n += s.size();

    cMatrix = GiNaC::matrix(n, n);

    for (uint i = 0; i < n; ++i)
        cMatrix(i, i) = 1;

	for (unsigned int k = 1; k < stacks.size(); k += 2)
	{
	    std::vector< std::pair<CellIndex, CellIndex> > A =
            this->AdjacencyLeft(k);
        std::vector< std::pair<CellIndex, CellIndex> >::iterator c, e;
        for (c = A.begin(), e = A.end(); c != e; ++c)
        {
            uint i = CellNumber(c->first), j = CellNumber(c->second);
            cMatrix(i, j) = 1;
            cMatrix(j, i) = 1;
        }

        A = this->AdjacencyRight(k);
        for (c = A.begin(), e = A.end(); c != e; ++c)
        {
            uint i = CellNumber(c->first), j = CellNumber(c->second);
            cMatrix(i, j) = 1;
            cMatrix(j, i) = 1;
        }
	}

	cMatrix = cMatrix.pow(n);

	for (uint i = 0; i < n; ++i)
        for (uint j = 0; j < n; ++j)
            if (cMatrix(i, j) != 0)
                cMatrix(i, j) = 1;
}

void CAD::Partition()
{
    typedef unsigned int uint;

    partitions.clear();

    // std::vector< std::vector<unsigned int> > partitions;
    std::vector<uint> cells;
    cells.reserve(cMatrix.rows());

    for (uint i = 0; i < cMatrix.rows(); ++i)
    {
        cells.clear();

        for (uint j = 0; j < cMatrix.cols(); ++j)
            if (cMatrix(i, j) != 0)
                cells.push_back(j);

        uint k = 0;
        for (k = 0; k < partitions.size(); ++k)
            if (partitions[k][0] == cells[0])
                break;

        if (k >= partitions.size())
            partitions.push_back(cells);

        // if (cells in partition)
        //     partitions.push_back(cells);
    }
}

CellIndex CAD::BranchCount(const CellIndex & ci)
{
    Algebraic & x = alphas[ci.first];
    Algebraic & y = stacks[ci.first][ci.second].y;

    unsigned int nroots = 0, L = 0, R = 0;

    while (true)
    {
        nroots = 0;

        for (PolynomialQQ::vector::iterator f = F.begin(), e = F.end();
             f != e; ++f)
            nroots += f->suby(y.lower()).sturm(x.lower(), x.upper()) +
            		  f->suby(y.upper()).sturm(x.lower(), x.upper());

        if (nroots == 0)
            break;

        x.tightenInterval();
    }

    for (PolynomialQQ::vector::iterator f = F.begin(), e = F.end(); f != e; ++f)
    {
        L += f->subx(x.lower()).sturm(y.lower(), y.upper());
        R += f->subx(x.upper()).sturm(y.lower(), y.upper());
    }

    return CellIndex(L, R);
}

bool ContainersEqual(const std::vector<char> & S1, const std::vector<char> & S2)
{
    assert(S1.size() == S2.size());

    for (std::vector<char>::size_type i = 0, e = S1.size(); i < e; ++i)
        if (S1[i] != S2[i])
            return false;

    return true;
}

std::vector< std::pair<CellIndex, CellIndex> >
    CAD::AdjacencyLeft(const unsigned int k)
{
    assert(Invariants());
    assert(k < stacks.size());
    assert(isOdd(k)); // Uses a private helper method.

    typedef unsigned int uint;

    std::vector<uint> r; // Find a cap so I can use 'reserve'.
    CellIndex bcount;

    for (uint i = 0; i < stacks[k].size(); ++i)
    {
    	if (isOdd(i))
    	{
    		bcount = BranchCount(CellIndex(k,i));
    		r.insert(r.end(), bcount.first, i);
    		// void insert ( iterator position, size_type n, const T& x );
    		// r := [op(r), seq(i,j=1..bcount[1])];
    	}
    	else
    		r.push_back(i);
    }

    std::vector<uint> l; // 1, 1); // l = { 1 };
    l.push_back(0);

    for (uint i = 1; i < r.size(); ++i)
    {
    	if (isOdd(r[i]) && isOdd(r[i-1]))
    		l.push_back(l.back()+2);
    	else if (isEven(r[i]) && isEven(r[i-1]))
    		l.push_back(l.back());
    	else
    		l.push_back(l.back()+1);
    }

    std::vector< std::pair<CellIndex, CellIndex> > c;

    for (uint i = 0, e = r.size(); i < e; ++i)
        if (ContainersEqual(stacks[k-1][l[i]].signs, stacks[k][r[i]].signs))
            c.push_back(make_pair(CellIndex(k-1,    l[i]),
                                  CellIndex(k,      r[i]) ));

    return c;
}

std::vector< std::pair<CellIndex, CellIndex> >
    CAD::AdjacencyRight(const unsigned int k)
{
    assert(Invariants());
    assert(k < stacks.size());
    assert(isOdd(k)); // Uses a private helper method.

    typedef unsigned int uint;

    std::vector<uint> l;
    CellIndex bcount;

    for (uint i = 0, e = stacks[k].size(); i < e; ++i)
    {
        if (isOdd(i))
        {
            bcount = BranchCount(CellIndex(k, i));
            l.insert(l.end(), bcount.second, i);
        }
        else
            l.push_back(i);
    }

    std::vector<uint> r; //(1, 1); // r = { 1 };
    r.push_back(0);

    for (uint i = 1, e = l.size(); i < e; ++i)
    {
        if (isOdd(l[i]) && isOdd(l[i-1]))
            r.push_back(r.back()+2);
        else if (isEven(l[i]) && isEven(l[i-1]))
            r.push_back(r.back());
        else
            r.push_back(r.back()+1);
    }

    std::vector< std::pair<CellIndex, CellIndex> > c;

    for (uint i = 0, e = l.size(); i < e; ++i)
        if (ContainersEqual(stacks[k][l[i]].signs, stacks[k+1][r[i]].signs))
            c.push_back(make_pair(CellIndex(k,      l[i]),
                                  CellIndex(k+1,    r[i]) ));

    return c;
}

/*!
 * \todo Reserve the space for P upfront.
 */
std::vector<PolynomialQ> CAD::Project(const std::vector<PolynomialQQ> & F)
{
    std::vector<PolynomialQQ> G = PolynomialQQ::IrreducibleFactors(F);
    std::vector<PolynomialQ> P;

    for (PolynomialQQ::vector::size_type i = 0;        i < G.size()-1; i++)
        for (PolynomialQQ::vector::size_type j = i+1;  j < G.size();   j++)
            P.push_back(PolynomialQQ::Resultant(G[i], G[j], 2));

    for (PolynomialQQ::vector::size_type i = 0;        i < G.size(); i++)
        P.push_back(PolynomialQQ::Resultant(G[i], G[i].getDerivative(2), 2));

    return PolynomialQ::IrreducibleFactors(P);
}

/*!
 * \param roots Must be an ordered vector of Algebraic numbers.
 * \return A vector of Algebraic numbers with new numbers inserted between
 *		   and on each end of the numbers in roots.
 */
std::vector<Algebraic>
    CAD::SamplePoints(const std::vector<Algebraic> & roots)
{
    std::vector<Algebraic> S;
    S.reserve(2*roots.size() + 1);

    if (roots.size() == 0)
    {
        S.push_back(Algebraic());
        return S;
    }

    // Assume roots.size() >= 1.

    GiNaC::numeric s = roots.front().lower() - 1;
    //S.push_back(Algebraic(CAD::MakePoly(s), IntervalQ(s)));
    S.push_back(Algebraic::MakeRational(s));

    S.push_back(roots[0]);

    for (unsigned int i = 1; i < roots.size(); i++)
    {
    	s = (roots[i-1].upper() + roots[i].lower())/2;
    	S.push_back(Algebraic::MakeRational(s));

    	S.push_back(roots[i]);
    }

    s = roots.back().upper() + 1;
    S.push_back(Algebraic::MakeRational(s));

	return S;
}

/*!
 * \param F A vector of non-zero polynomials.
 */
std::vector<Algebraic>
    CAD::FindRoots2(const Algebraic & alpha,
    				const std::vector<PolynomialQQ> & F)
{
    PolynomialQQ fs = accumulate(F.begin(),
                                 F.end(),
                                 PolynomialQQ((GiNaC::ex)1),
                                 std::multiplies<PolynomialQQ>());

    assert(!fs.isZero());

    // make a conversion operator from PolyQ to PolyQQ for alpha.
    // maybe PolynomialQ::getPolynomialQQ?
    PolynomialQ r(PolynomialQQ::Resultant(PolynomialQQ(alpha.getPolynomial().getEx()),
                                          fs,
                                          1));

    Algebraic::vector P = PolynomialQ::FindRoots(*r.getIrreducibleFactors<PolynomialQ::vector>());
    Algebraic::vector R;
    R.reserve(P.size());

    for (Algebraic::vector::const_iterator p = P.begin(), e = P.end();
         p != e; ++p)
        if (fs.signAt(alpha, *p) == 0)
            R.push_back(*p);

    return R;
}
