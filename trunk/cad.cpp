/*!
 * \file
 * \author Chris de Pujo
 */

#include "cad.hpp"

#include <cassert>
#include <numeric>
#include <iterator>

#include <boost/foreach.hpp>

#include "templatehelp.hpp"

bool CAD::Connectivity(const Point & p1, const Point & p2) const
{
    assert(Invariants());

    return (cMatrix(CellNumber(Cell(p1)),
                    CellNumber(Cell(p2))) == 1);
}

StacksIndex CAD::Cell(const Point & p) const
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

	std::vector<Algebraic> betas; /**/ betas.reserve(stacks[ix].size());
	CAD::addRoots2To(betas, p.first, this->irreducibles);

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

	return StacksIndex(ix, iy);
}

unsigned int CAD::CellNumber(const StacksIndex & i) const
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

bool CAD::SignsNonzero(const std::vector<char> & s)
{
    for (unsigned int i = 0, e = s.size(); i < e; ++i)
        if (s[i] == 0)
            return false;

    return true; // ( std::find(s.begin(), s.end(), 0) == s.end() )
}

bool CAD::SignsEqual(const std::vector<char> & S1, const std::vector<char> & S2)
{
    if (S1.size() != S2.size())
    {
        int i;
        i = 5;
        ++i;
    }
    assert(S1.size() == S2.size());

    for (std::vector<char>::size_type i = 0, e = S1.size(); i < e; ++i)
        if (S1[i] != S2[i])
            return false;

    return true; // std::equal(S1.begin(), S1.end(), S2.begin());
}

bool CAD::SignsDiffByOne(const std::vector<char> & s1,
                         const std::vector<char> & s2)
{
    assert(s1.size() == s2.size());

    unsigned int count = 0;

    for (unsigned int i = 0, e = s1.size(); i < e; ++i)
        if (s1[i] != s2[i])
            ++count;

    return (count == 1);
}

const CAD::Sample & CAD::RefCell(CellIndex i) const
{
    for (std::vector< std::vector<Sample> >::const_iterator s = stacks.begin(),
         e = stacks.end(); s != e; ++s)
    {
        if (i < s->size())
            return s->at(i);

        i -= s->size();
    }

    assert(false);

    return *(stacks.back().end());
}

bool CAD::areComponentsAdjacent(ComponentIndex i, ComponentIndex j) const
{
    assert(i < components.size());
    assert(j < components.size());

    return SignsDiffByOne(RefCell(partitions[components[i]][0]).signs,
                          RefCell(partitions[components[j]][0]).signs);
}

void MatrixOut(GiNaC::matrix & M)
{
    using std::cout;
    using std::endl;
    typedef unsigned int uint;

    cout << M << endl;
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

	for (uint i = 1; i < F.size(); ++i)
		cout << "           " << F[i] << endl;

    cout << "Number of irreducibles: " << irreducibles.size() << endl;
	cout << "Irreducibles: ";

	if (irreducibles.size() >= 1)
		cout << irreducibles[0] << endl;

	for (uint i = 1; i < irreducibles.size(); ++i)
		cout << "              " << irreducibles[i] << endl;

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
            cout << cMatrix(i, j);
        cout << endl;
	}

	cout << "Partitions: " << endl;

	for (uint i = 0; i < partitions.size(); ++i)
    {
        for (uint j = 0; j < partitions[i].size(); ++j)
            cout << partitions[i][j] << " ";
        cout << endl;
    }

    cout << "Components(" << components.size() << "): " << endl;

    for (uint i = 0; i < components.size(); ++i)
    {
        cout << components[i] << ": ";
        for (uint j = 0; j < partitions[components[i]].size(); ++j)
            cout << partitions[components[i]][j] << " ";
        cout << endl;
    }
}

/*!
 * \todo Update to reflect recent changes.
 */
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

template <typename T>
void PrintContainer(const T & c)
{
    std::copy(c.begin(), c.end(),
              std::ostream_iterator<typename T::value_type>(std::cout, "\n"));
    std::cout << std::endl;
}

template <typename T>
void PrintContainerAlgebraicf(const T & c)
{
    for (typename T::const_iterator i = c.begin(), e = c.end(); i != e; ++i)
        std::cout << *i << "\t\t" << i->approx(3).to_double() << std::endl;
    std::cout << std::endl;
}

void CAD::ConstructorPart2()
{
    for (unsigned int i = 0, e = this->F.size(); i < e; ++i)
        this->F[i].addIrreducibleFactorsTo(this->irreducibles);
    //tmp::Unique(this->irreducibles);
    assert(tmp::IsUnique(this->irreducibles));

    PolynomialQ::vector tempQ;
    CAD::addIrreducibleProjectionTo(tempQ, this->irreducibles);
    Algebraic::set setA;

    for (PolynomialQ::vector::const_iterator
         f = tempQ.begin(), e = tempQ.end(); f != e; ++f)
        f->addRootsTo(setA);

    Algebraic::SeparateIntervals(setA);

    alphas.reserve(2*setA.size() + 1);
    CAD::addSamplePointsTo(alphas, setA);

    assert(tmp::IsSortedUnique(alphas));

    stacks = std::vector< std::vector<Sample> >(alphas.size());
    std::vector< std::vector<Sample> >::iterator T = stacks.begin();
    setA.clear();

    for (Algebraic::vector::iterator
         alpha = alphas.begin(), e = alphas.end(); alpha != e; ++alpha)
    {
        CAD::addRoots2To(setA, *alpha, this->irreducibles);
        Algebraic::SeparateIntervals(setA);
        Algebraic::vector betas;
        CAD::addSamplePointsTo(betas, setA);

        T->reserve(betas.size());

        for (Algebraic::vector::const_iterator
             beta = betas.begin(), e = betas.end(); beta != e; ++beta)
        {
            T->push_back(this->F.size());

            for (unsigned int i = 0, e = this->F.size(); i < e; ++i)
                T->back().signs[i] = this->F[i].signAt2(*alpha, *beta);

            T->back().y = *beta;
        }

        ++T;
        setA.clear();
    }

    MakeConnectivityMatrix();

    Partition();

    // maybe this should be an assert?
    if (!Invariants())
        throw std::logic_error("CAD Constructor: CAD creation completed, but"
                               "the result is not canonical.");
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
	    std::vector< std::pair<StacksIndex, StacksIndex> > A =
            this->AdjacencyLeft(k);
        std::vector< std::pair<StacksIndex, StacksIndex> >::iterator c, e;
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
    }

    for (PartitionIndex i = 0; i < partitions.size(); ++i)
    {
        if (SignsNonzero(RefCell(partitions[i][0]).signs))
            components.push_back(i);
    }
}

StacksIndex CAD::BranchCount(const StacksIndex & ci)
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

    return StacksIndex(L, R);
}

std::vector< std::pair<StacksIndex, StacksIndex> >
    CAD::AdjacencyLeft(const unsigned int k)
{
    assert(Invariants());
    assert(k < stacks.size());
    assert(isOdd(k)); // Uses a private helper method.

    typedef unsigned int uint;

    std::vector<uint> r; // Find a cap so I can use 'reserve'.
    StacksIndex bcount;

    for (uint i = 0; i < stacks[k].size(); ++i)
    {
    	if (isOdd(i))
    	{
    		bcount = BranchCount(StacksIndex(k,i));
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

    assert(l.size() == r.size());

    std::vector< std::pair<StacksIndex, StacksIndex> > c;

    for (uint i = 0, e = r.size(); i < e; ++i)
        if (SignsEqual(stacks[k-1][l[i]].signs, stacks[k][r[i]].signs))
            c.push_back(make_pair(StacksIndex(k-1,    l[i]),
                                  StacksIndex(k,      r[i]) ));

    return c;
}

std::vector< std::pair<StacksIndex, StacksIndex> >
    CAD::AdjacencyRight(const unsigned int k)
{
    assert(Invariants());
    assert(k < stacks.size());
    assert(isOdd(k)); // Uses a private helper method.

    typedef unsigned int uint;

    std::vector<uint> l;
    StacksIndex bcount;

    for (uint i = 0, e = stacks[k].size(); i < e; ++i)
    {
        if (isOdd(i))
        {
            bcount = BranchCount(StacksIndex(k, i));
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

    assert(l.size() == r.size());

    std::vector< std::pair<StacksIndex, StacksIndex> > c;

    for (uint i = 0, e = l.size(); i < e; ++i)
        if (SignsEqual(stacks[k][l[i]].signs, stacks[k+1][r[i]].signs))
            c.push_back(make_pair(StacksIndex(k,      l[i]),
                                  StacksIndex(k+1,    r[i]) ));

    return c;
}

/*!
 * \todo Reserve the space for P upfront.
 */
PolynomialQ::vector CAD::Project(const PolynomialQQ::vector & F)
{
    PolynomialQQ::vector G = PolynomialQQ::IrreducibleFactors(F);
    PolynomialQ::vector P;

    for (PolynomialQQ::vector::size_type i = 0;        i < G.size()-1; i++)
        for (PolynomialQQ::vector::size_type j = i+1;  j < G.size();   j++)
            P.push_back(PolynomialQQ::Resultant(G[i], G[j], 2));

    for (PolynomialQQ::vector::size_type i = 0;        i < G.size(); i++)
        P.push_back(PolynomialQQ::Resultant(G[i], G[i].getDerivative(2), 2));

    return PolynomialQ::IrreducibleFactors(P);
}

//ProjectIrreducibles
void CAD::addIrreducibleProjectionTo(PolynomialQ::vector & P,
                                     const PolynomialQQ::vector & G)
{
    P.reserve(P.size() + (G.size()*(G.size() + 1))/2);

    for (PolynomialQQ::vector::size_type i = 0;        i < G.size()-1; i++)
        for (PolynomialQQ::vector::size_type j = i+1;  j < G.size();   j++)
            P.push_back(PolynomialQQ::Resultant(G[i], G[j], 2));

    for (PolynomialQQ::vector::size_type i = 0;        i < G.size(); i++)
        P.push_back(PolynomialQQ::Resultant(G[i], G[i].getDerivative(2), 2));
}

/*!
 * \param roots Must be an ordered vector of Algebraic numbers.
 * \return A vector of Algebraic numbers with new numbers inserted between
 *		   and on each end of the numbers in roots.
 */
std::vector<Algebraic> CAD::SamplePoints(const std::vector<Algebraic> & roots)
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
    S.push_back(Algebraic::MakeWideRational(s));

    S.push_back(roots[0]);

    for (unsigned int i = 1; i < roots.size(); i++)
    {
    	s = (roots[i-1].upper() + roots[i].lower())/2;
    	S.push_back(Algebraic::MakeWideRational(s));

    	S.push_back(roots[i]);
    }

    s = roots.back().upper() + 1;
    S.push_back(Algebraic::MakeWideRational(s));

	return S;
}
