/*!
 * \file
 * \author Chris de Pujo
 */

#include "cad.hpp"

#include <cassert>

#include <boost/foreach.hpp>

CAD::CAD(const std::list<std::string> & F)
{
    // can constructor be called on already instantiated object?
    // Be careful with 'F' and 'this->F'.

    BOOST_FOREACH(const std::string & f, F)
        this->F.push_back(PolynomialQQ(f)); // throws on error

    samples = CAD::SamplePoints(PolynomialQ::FindRoots(CAD::Project(this->F)));
    std::vector<Algebraic> & alphas = samples;

    stacks.reserve(alphas.size());

    BOOST_FOREACH(const Algebraic & alpha, alphas)
    {
        std::vector<Algebraic> betas =
            CAD::SamplePoints(CAD::FindRoots2(alpha, this->F));

        stacks.push_back(std::vector<Sample>()); // Should I name this parameter?
        std::vector<Sample> & T = stacks.back;

        T.reserve(betas.size());

        BOOST_FOREACH(const Algebraic & beta, betas)
        {
            std::vector<Algebraic> numbers;
            numbers.reserve(2);
            numbers.push_back(alpha);
            numbers.push_back(beta);

            std::vector<char> signs(this->F.size()); // allocates 'size' of
            										 // these up front.

            for (int i = 0; i < this->F.size(); i++)
                signs[i] = Sign(numbers, (this->F)[i]);

            T.push_back(Sample(beta, signs));
        }
    }

    // maybe this should be an assert?
    if (!Invariants())
        throw logic_error("CAD Constructor: CAD creation completed, but the"
                          "result is not canonical.");
}

bool CAD::Connectivity(const Point & p1, const Point & p2) const
{
    assert(Invariants());

    return (M[CellNumber(Cell(p1))][CellNumber(Cell(p2))] == 1);
}

CellIndex CAD::Cell(const Point & p) const
{
    assert(Invariants());
}

unsigned int CAD::CellNumber(const CellIndex & i) const
{
    // Indices start at 0, not 1.
    assert(Invariants());
    assert(i.first < stacks.size());
    assert(i.second < stacks[i.first].size());

    unsigned int num = 0;

    for (int j = 0; j < i.first; j++)
        num += stacks[j].size();

    num += i.second;

    return num;
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

    BOOST_FOREACH(const std::vector<Sample> & stack, stacks)
        BOOST_FOREACH(const Sample & sample, stack)
        {
            if (!sample.x.Invariants() || !sample.y.Invariants())
                return false;

            if (sample.signs != -1 &&
                sample.signs !=  0 &&
                sample.signs !=  1)
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
		cMatrix(i, j) = 1; // The diagonal already set and symmetry is
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
            nroots += sturm(sturmseq(eval(f,y=ycoord.lower()),x),x,op(xcoord[1])) +
                      sturm(sturmseq(eval(f,y=ycoord.upper()),x),x,op(xcoord[1]));

        if (nroots == 0) then
            break;

        xcoord.tightenInterval();
    }

    BOOST_FOREACH(const PolynomialQQ & f, F)
    {
        L += sturm(sturmseq(eval(f,x=xcoord.lower()),y),y,op(ycoord[1]));
        R += sturm(sturmseq(eval(f,x=xcoord.upper()),y),y,op(ycoord[1]));
    }

    return CellIndex(L, R);
}

void CAD::AdjacencyLeft(const unsigned int k)
{
    assert(Invariants());
    assert(k < stacks.size());
    assert(isOdd(k)); // Uses a private helper method.

    std::list<unsigned int> r;

    for (unsigned int i = 0; i < S[k].size(); i++) // i from 1 to nops(S[k])
    {
        if (isOdd(i))
        {
            CellIndex bcount = BranchCount(k,i);
            // r := [op(r), seq(i,j=1..bcount.first)];
            r.insert(r.end(), bcount.first, i);
        }
        else
            r.push_back(i); // r := [op(r), i];
    }

    std::list<unsigned int> l(1, 1); // l := [1]

    for (unsigned int i = 1; i < r.size(); i++) // i from 2 to nops(r)
    {
        if (isOdd(r[i]) && isOdd(r[i-1]))
            l := [op(l), l[i-1]+2];
        else if (isEven(r[i]) && isEven(r[i-1]))
            l := [op(l), l[i-1]];
        else
            l := [op(l), l[i-1]+1];
    }

    c := [];

    for i from 1 to nops(r)
        if evalb(S[k-1][l[i]][3]=S[k][r[i]][3])
            c := [op(c), [[k-1,l[i]],[k,r[i]]]];

    return c;
}

std::vector<PolynomialQ> CAD::Project(const std::vector<PolynomialQQ> & F) const
{
    std::vector<PolynomialQQ> G = PolynomialQQ::IrreducibleFactors(F);
    std::vector<PolynomialQ> P;

    for (int i = 0; i < G.size()-1; i++)
        for (int j = i+1; j < G.size(); j++)
            P.push_back(PolynomialQQ::Resultant(G[i], G[j], 2));

    for (int i = 0; i < G.size(); i++)
        P.push_back(PolynomialQQ::Resultant(G[i], G[i].getDerivative(), 2));

    return PolynomialQ::IrreducibleFactors(P);
}

/*!
 * \param roots Must be an ordered vector of Algebraic numbers.
 * \return A vector of Algebraic numbers with new numbers inserted between
 *		   and on each end of the numbers in roots.
 */
std::vector<Algebraic>
    CAD::SamplePoints(const std::vector<Algebraic> & roots) const;
{
    std::vector<Algebraic> S;
    S.reserve(2*roots.size() + 1);

    if (roots.size() == 0)
    {
        S.push_back(Algebraic());
        return points;
    }

    // Assume roots.size() >= 1.
    
    GiNaC::numeric s = roots.front().lower() - 1;
    S.push_back(Algebraic(CAD::MakePoly(s), IntervalQ(s)));
    
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
	
	return PolynomialQ(GiNaC::ex(p.getVar() - num));
}

/*!
 * \param F Should not contain 0 polynomials; this 
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

    // make a conversion operator from PolyQ to PolyQQ for alpha.
    // maybe PolynomialQ::getPolynomialQQ?
    PolynomialQ r = PolynomialQQ::Resultant(alpha.getPolynomial(), fs, 1);

    std::vector<Algebraic> P = PolynomialQ::FindRoots(r.getIrreducibleFactors());
    std::vector<Algebraic> R;
    R.reserve(P.size());

    BOOST_FOREACH(const Algebraic & p, P)
        if (fs.signAt(alpha, p) == 0)
            R.push_back(p);

    return R;

/*  local Fs,r,G,P,R,i;
 Fs := product(F[i],i=1..nops(F));
 r  := resultant(alpha[2],Fs,x);
 G  := IrreducibleFactors([r]);
 P  := FindRoots(G);
 R  := select(p->Sign([alpha,p],Fs)=0,P);
 return R; */
}
