/*!
 * \file
 * \author Chris de Pujo
 */

#include "cad.hpp"

#include <cassert>

#include <boost/foreach.hpp>

/*!
 * \throws parse_error Thrown by GiNaC if parsing of polynomials fails (inherits
 *		   invalid_argument).
 * \todo Check if doxygen collates all 'throws' statements from called
 *		 functions.
 */
CAD::CAD(const std::vector<std::string> & F)
{
    // can constructor be called on already instantiated object?
    // Be careful with 'F' and 'this->F'.
std::cout << "Functions: " << std::endl;
    BOOST_FOREACH(const std::string & f, F)
{        this->F.push_back(PolynomialQQ(f)); // throws on error
std::cout << this->F.back() << std::endl;}
    samples = CAD::SamplePoints(PolynomialQ::FindRoots(CAD::Project(this->F)));
    std::vector<Algebraic> & alphas = samples;
std::cout << "Samples: " << std::endl;
BOOST_FOREACH(const Algebraic & alpha, alphas)
std::cout << alpha << std::endl;
    stacks.reserve(alphas.size());

    BOOST_FOREACH(const Algebraic & alpha, alphas)
    {
        std::vector<Algebraic> betas =
            CAD::SamplePoints(CAD::FindRoots2(alpha, this->F));

        stacks.push_back(std::vector<Sample>()); // Should I name this parameter?
        std::vector<Sample> & T = stacks.back();

        T.reserve(betas.size());
std::cout << "Betas: " << std::endl;
        BOOST_FOREACH(const Algebraic & beta, betas)
        {
            std::vector<char> signs(this->F.size()); // allocates 'size' of
            										 // these up front.
std::cout << beta << ": ";
            for (std::vector<char>::size_type i = 0;
                 i < this->F.size();
                 i++)
{                signs[i] = this->F[i].signAt(alpha, beta);
std::cout << int(signs[i]) << " ";}
std::cout << std::endl;
            T.push_back(Sample(beta, signs));
        }
    }

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

	int ix = samples.size();

	for (unsigned int i = 1; i < samples.size(); i += 2) // i from 1 to nops(alphas)
	{
		int t = p.first.compare(samples[i]);

		if (t <= 0)
		{
			if (t == 0)
				ix = i;
			else
				ix = i - 1;

			break;
		}
	}

	std::vector<Algebraic> betas = CAD::FindRoots2(p.first, F);

	int iy = 2*betas.size();

	for (unsigned int i = 0; i < betas.size(); i++) // i from 1 to nops(betas)
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

	cout << "Number of samples: " << samples.size() << endl;
	cout << "Number of stacks: " << stacks.size() << endl;
	cout << "Individual stack sizes: " << endl;

	for (uint i = 0; i < stacks.size(); i++)
		cout << "    " << "Size of stack " << i << ": " << stacks[i].size() << endl;

	cout << "Samples: " << endl;

	for (uint i = 0; i < samples.size(); i++)
		cout << "    " << samples[i] << endl;

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
}

bool CAD::Invariants() const
{
    if (this == NULL)
        return false;

    if (stacks.empty())
        return false;

    if (samples.size() != stacks.size())
        return false;

    BOOST_FOREACH(const std::vector<Sample> & stack, stacks)
        if (stack.empty())
            return false;

    BOOST_FOREACH(const Algebraic & alpha, samples)
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

void CAD::ConnectivityMatrix()
{
	const std::vector< std::vector<Sample> > & S = stacks;

	unsigned int n = 0;

	BOOST_FOREACH(const std::vector<Sample> & s, S)
		n += s.size();

	// matrix mat(n,n,identity);

	for (unsigned int k = 1; k < S.size(); k += 2)
	{
		cMatrix(k,k/*i, j*/) = 1; // The diagonal already set and symmetry is
						   // maintained internally by the matrix class.
	}
/*
 local n,i,j,k,adjM,a,c,s;
 n := add(nops(s), s in S);
 adjM := Matrix(n, (i,j) -> `if`(i=j,1,0));
 for k from 2 to nops(S) by 2 do
   a := [op(AdjacencyLeft(F,S,k)), op(AdjacencyRight(F,S,k))];
   for c in a do
     i := CellNumber(S, c[1]);
     j := CellNumber(S, c[2]);
     adjM[i,j] := 1;
     adjM[j,i] := 1;
   od;
 od;
 return Closure(adjM);
 */
}

CellIndex CAD::BranchCount(const CellIndex & ci)
{
    Algebraic & xcoord = samples[ci.first];
    Algebraic & ycoord = stacks[ci.first][ci.second].y;

    unsigned int nroots = 0, L = 0, R = 0;

    while (true)
    {
        nroots = 0;

        BOOST_FOREACH(const PolynomialQQ & f, F)
            nroots += f.suby(ycoord.lower()).sturm(xcoord.lower(),
            									   xcoord.upper()) +
            		  f.suby(ycoord.upper()).sturm(xcoord.lower(),
            									   xcoord.upper());

        if (nroots == 0)
            break;

        xcoord.tightenInterval();
    }

    BOOST_FOREACH(const PolynomialQQ & f, F)
    {
        L += f.subx(xcoord.lower()).sturm(ycoord.lower(), ycoord.upper());
        R += f.subx(xcoord.upper()).sturm(ycoord.lower(), ycoord.upper());
    }

    return CellIndex(L, R);
}

void CAD::AdjacencyLeft(const unsigned int k)
{
    assert(Invariants());
    assert(k < stacks.size());
    assert(isOdd(k)); // Uses a private helper method.

    typedef unsigned int uint;

    std::vector<uint> r; // Find a cap so I can use 'reserve'.
    CellIndex bcount;

    for (uint i = 0; i < stacks[k].size(); i++)
    {
    	if (isOdd(i))
    	{
    		bcount = BranchCount(CellIndex(k,i));
    		r.insert(r.end(), bcount.first, i);
    	}
    	else
    		r.push_back(i);
    }

    std::vector<uint> l(1, 1); // l = { 1 };

    for (uint i = 1; i < r.size(); i++)
    {
    	if (isOdd(r[i]) && isOdd(r[i-1]))
    		l.push_back(l[i-1]+2);
    	else if (isEven(r[i]) && isEven(r[i-1]))
    		l.push_back(l[i-1]);
    	else
    		l.push_back(l[i-1]+1);
    }

    std::vector<CellIndex> c;
/*
    for (uint i = 0; i < r.size(); i++)
    	if (stacks[k-1][...].signs == stacks[k][...].signs) // Do pairwise comparison.
    		c.push_back(std::pair<CellIndex, CellIndex>(CellIndex(,),
    													CellIndex(,)));
*/
//    return c;

/*
 for i from 1 to nops(r) do
   if evalb(S[k-1][l[i]][3]=S[k][r[i]][3]) then
     c := [op(c), [[k-1,l[i]],[k,r[i]]]];
   fi;
 od;
 return c; */
}

/*!
 * \todo Reserve the space for P upfront.
 */
std::vector<PolynomialQ> CAD::Project(const std::vector<PolynomialQQ> & F)
{
    std::vector<PolynomialQQ> G = PolynomialQQ::IrreducibleFactors(F);
    std::vector<PolynomialQ> P;

    for (std::vector<PolynomialQQ>::size_type i = 0;        i < G.size()-1; i++)
        for (std::vector<PolynomialQQ>::size_type j = i+1;  j < G.size();   j++)
            P.push_back(PolynomialQQ::Resultant(G[i], G[j], 2));

    for (std::vector<PolynomialQQ>::size_type i = 0;        i < G.size(); i++)
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
    S.push_back(Algebraic(CAD::MakePoly(s), IntervalQ(s)));

    S.push_back(roots[0]);

    for (unsigned int i = 1; i < roots.size(); i++)
    {
    	s = (roots[i-1].upper() + roots[i].lower())/2;
    	S.push_back(Algebraic(CAD::MakePoly(s), IntervalQ(s)));

    	S.push_back(roots[i]);
    }

    s = roots.back().upper() + 1;
    S.push_back(Algebraic(CAD::MakePoly(s), IntervalQ(s)));

	return S;
/*local n,S,s,w,i,lb,ub;
 n := nops(R);
 if n = 0 then
   return [[[0,0],v-0]];
 fi;
 s := R[1][1][1]-1;
 S  := [[[s,s],v-s],R[1]];
 for i from 2 to n do
   lb := R[i-1][1][2];
   ub := R[i  ][1][1];
   s  := (ub + lb)/2;
   S  := [op(S),[[s,s],v-s],R[i]];
 od;
 s := R[n][1][2]+1;
 S  := [op(S),[[s,s],v-s]];
 return S;	*/
}

inline PolynomialQ CAD::MakePoly(const GiNaC::numeric & num)
{
	PolynomialQ p;

	return PolynomialQ(GiNaC::ex(p.getVariable() - num));
}

/*!
 * \param F A vector of non-zero polynomials.
 */
std::vector<Algebraic>
    CAD::FindRoots2(const Algebraic & alpha,
    				const std::vector<PolynomialQQ> & F)
{
    PolynomialQQ fs((GiNaC::ex)1);

    // Change this to use 'accumulate',
    // or implement & use some mul(vector) method.
    BOOST_FOREACH(const PolynomialQQ & f, F)
        fs *= f;

    assert(!fs.isZero());
/**/std::cout << "fs: " << fs << std::endl;
    // make a conversion operator from PolyQ to PolyQQ for alpha.
    // maybe PolynomialQ::getPolynomialQQ?
    PolynomialQ r(PolynomialQQ::Resultant(PolynomialQQ(alpha.getPolynomial().getEx()),
                                          fs,
                                          1));
/**/std::cout << "r: " << r << std::endl;
    std::vector<Algebraic> P = PolynomialQ::FindRoots(r.getIrreducibleFactors());
    std::vector<Algebraic> R;
    R.reserve(P.size());
std::cout << "Size of P: " << P.size() << std::endl;
BOOST_FOREACH(const Algebraic & p, P)
std::cout << p << std::endl;
std::cout << "R's: " << std::endl;
    // Why did I add this step? Can't remember, but removing it breaks the code.
    BOOST_FOREACH(const Algebraic & p, P)
    {
        if (fs.signAt(alpha, p) == 0)
        {
            std::cout << p << std::endl;
            R.push_back(p);
        }
    }

    return R;

/*  local Fs,r,G,P,R,i;
 Fs := product(F[i],i=1..nops(F));
 r  := resultant(alpha[2],Fs,x);
 G  := IrreducibleFactors([r]);
 P  := FindRoots(G);
 R  := select(p->Sign([alpha,p],Fs)=0,P);
 return R; */
}
