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
        this->F.push_back(PolynomialQQ(f)); // throws invalid_argument on error

    samples = SamplePoints(FindRoots(Project()));
    std::vector<Algebraic> & alphas = samples;

    stacks.reserve(alphas.size());

    BOOST_FOREACH(const Algebraic & alpha, alphas)
    {
        std::vector<Algebraic> betas =
            SamplePoints(CAD::FindRoots2(alpha, this->F));

        stacks.push_back(std::vector<Sample>());
        std::vector<Sample> & T = stacks.back;

        T.reserve(betas.size());

        BOOST_FOREACH(const Algebraic & beta, betas)
        {
            std::vector<Algebraic> numbers;
            numbers.reserve(2);
            numbers.push_back(alpha);
            numbers.push_back(beta);

            std::vector<char> signs(this->F.size());

            for (int i = 0; i < this->F.size(); i++)
                signs[i] = Sign(numbers, (this->F)[i]);

            T.push_back(Sample(beta, signs));
        }
    }

    // maybe this should be an assert?
    if (!Invariant())
        throw logic_error("CAD Constructor: CAD creation completed, but the"
                          "result is not canonical.");
}

bool CAD::Connectivity(const Point & p1, const Point & p2) const
{
    assert(Invariant());

    return (M[CellNumber(Cell(p1))][CellNumber(Cell(p2))] == 1);
}

CellIndex CAD::Cell(const Point & p) const
{
    assert(Invariant());
}

unsigned int CAD::CellNumber(const CellIndex & i) const
{
    // Indices start at 0, not 1.
    assert(Invariant());
    assert(i.first < stacks.size());
    assert(i.second < stacks[i.first].size());

    unsigned int num = 0;

    for (int j = 0; j < i.first; j++)
        num += stacks[j].size();

    num += i.second;

    return num;
}

bool CAD::Invariant() const
{
    // make sure stacks is not empty.
    // make sure there are no empty stacks in 'stacks'.
    // make sure all signs are -1, 0, or 1.
    // check all algebraics are invariant.

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
            if (!sample.x.Invariant() || !sample.y.Invariant())
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
}

CellIndex CAD::BranchCount(const CellIndex & ci)
{
    Algebraic & xcoord = samples[ci.first]; // stacks[ci.first][ci.second].x;
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

        xcoord.TightenInterval();
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
    assert(Invariant());
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

std::vector<Algebraic>
    CAD::SamplePoints(const std::vector<Algebraic> & roots) const;
{
    std::vector<Algebraic> points;
    points.reserve(2*roots.size() + 1);

    if (roots.size() == 0)
    {
        points.push_back(Algebraic());
        return points;
    }

    // Assume roots.size() >= 1.

    BOOST_FOREACH(const Algebraic & a, roots)
}
