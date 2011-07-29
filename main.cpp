#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include <boost/foreach.hpp>
#include <boost/tuple/tuple.hpp>

#include <ginac/ginac.h>

#include "cad.hpp"
#include "polynomialq.hpp"
#include "polynomialqq.hpp"

std::ostream & operator<<(std::ostream & output, const GiNaC::numeric & num)
{
    output << GiNaC::ex(num);

    return output;
}

void TestGCD(const PolynomialQ & f, const PolynomialQ & g)
{
    using namespace std;
    using namespace GiNaC;

    ex s,t;
    ex gcd = PolynomialQQ::gcdex
        (f.getEx(), g.getEx(), PolynomialQ::GetVar(), s, t);
    PolynomialQ gg(gcd);
    cout << "**********************************************************************" << endl;
    cout << "f: " << f << ". g: " << g << "." << endl;
    cout << "gcdex: " << gcd << ". s: " << s << ". t: " << t << "." << endl;
    cout << "f*s + g*t: " << (f.getEx()*s+g.getEx()*t).expand() << endl;
    cout << "f/gcd: " <<  f/gg << ". f%gcd: " << f%gg << "." << endl;
    cout << "f/gcd: " <<  g/gg << ". f%gcd: " << g%gg << "." << endl;
    cout << "**********************************************************************" << endl;
}

void TestSimple2(const Algebraic & alpha, const Algebraic & beta)
{
    using namespace std;
    using namespace GiNaC;

    boost::tuple<Algebraic, PolynomialQ, PolynomialQ> t =
        PolynomialQQ::Simple2(alpha, beta);

    Algebraic &     gamma   = t.get<0>();
    PolynomialQ &   S       = t.get<1>();
    PolynomialQ &   T       = t.get<2>();

    cout << "Simple2: alpha = " << alpha << ", beta = " << beta << endl;
    cout << "\tgamma: " << gamma << endl;
    cout << "\tS: "     << S << endl;
    cout << "\tT: "     << T << endl;

    numeric a = alpha.Approximate();
    numeric b = beta.Approximate();
    numeric g = gamma.Approximate();

    cout << "\t" << ex(S.eval(g)) << " ?= " << ex(a) << endl;
    cout << "\t" << ex(T.eval(g)) << " ?= " << ex(b) << endl;
}

void TestANComb(const Algebraic & alpha,
                const Algebraic & beta,
                const GiNaC::numeric & t)
{
    using namespace std;
    using namespace GiNaC;

    Algebraic gamma = PolynomialQQ::ANComb(alpha, beta, t);

    cout << "ANComb( " << alpha << " , " << beta << " , " << t << " ) =" << endl;
    cout << "\t" << gamma << " ~=" << endl;
    cout << "\t" << gamma.approx(6).to_double() << " ?= "
         << (alpha.approx(6)+t*beta.approx(6)).to_double() << endl;
}

int main()
{
    using namespace std;
    using namespace GiNaC;

    // Code for CAD stuff /////////////////////////////////////////////////////

    try
    {
        {
            vector<string> F;
            //F.push_back("x^2+y^2-4");
            //F.push_back("x+y");
            //F.push_back("x-y");
            //F.push_back("x+y");
            //F.push_back("x^2+y^2-1");
            //F.push_back("2*x+3*y^2-2");
            //F.push_back("-x+y^2");
            //F.push_back("x-y-5");
            F.push_back("(x^2 + y^2)^4 - (x^2 - y^2)^2");
            F.push_back("x + y + 1/2");

            CAD mainCAD(F);

            mainCAD.out();
        }
        /*{
            PolynomialQQ Fs("(-x+y^2)*(x-y-5)");
            Algebraic alpha(PolynomialQ("25 - 11*x + x^2"),
                            IntervalQ(numeric(247,32),numeric(143,16)));

            Algebraic a1(PolynomialQ("x^2+x-5"),
                         IntervalQ(numeric(-3),numeric(-21,8)));
            Algebraic a2(PolynomialQ("x^2-x-5"),
                         IntervalQ(numeric(-15,8),numeric(-3,2)));
            Algebraic a3(PolynomialQ("x^2+x-5"),
                         IntervalQ(numeric(3,2),numeric(15,8)));
            Algebraic a4(PolynomialQ("x^2-x-5"),
                         IntervalQ(numeric(21,8),numeric(3)));

            cout << Fs.signAt(alpha, a1) << endl;
            cout << Fs.signAt(alpha, a2) << endl;
            cout << Fs.signAt(alpha, a3) << endl;
            cout << Fs.signAt(alpha, a4) << endl;

            cout << "************************************" << endl;
        }*/
        /*{
            // [[[1,2],x^2-2], [[-2,-1],x^2-2]]
            Algebraic alpha(PolynomialQ("x^2-2"), IntervalQ(1,2));
            Algebraic beta(PolynomialQ("x^2-2"), IntervalQ(-2,-1));
            boost::tuple<Algebraic, PolynomialQ, PolynomialQ> t =
                PolynomialQQ::Simple(alpha, beta);

            cout << "Simple:" << endl;
            cout << "gamma: " << t.get<0>() << endl;
            cout << "S: " << t.get<1>() << endl;
            cout << "T: " << t.get<2>() << endl;

            alpha = Algebraic(PolynomialQ("x^2-2"), IntervalQ(1,2));
            beta = Algebraic(PolynomialQ("x^2-3"), IntervalQ(1,2));

            TestSimple2(alpha, beta);

            alpha = Algebraic(PolynomialQ("x^2-3"), IntervalQ(1,2));
            beta = Algebraic(PolynomialQ("x^2-5"), IntervalQ(2,3));

            TestSimple2(alpha, beta);

            alpha = Algebraic(PolynomialQ("x^2-2"), IntervalQ(1,2));
            beta = Algebraic(PolynomialQ("x^2-2"), IntervalQ(-2,-1));

            TestSimple2(alpha, beta);

            // ------------------------------------------------------------

            alpha = Algebraic(PolynomialQ("x^2-2"), IntervalQ(1,2));
            beta = Algebraic(PolynomialQ("x^2-3"), IntervalQ(1,2));
            TestANComb(alpha, beta, 1);
            TestANComb(alpha, beta, 2);

            alpha = Algebraic(PolynomialQ("x^2-3"), IntervalQ(1,2));
            beta  = Algebraic(PolynomialQ("x^2-5"), IntervalQ(2,3));
            TestANComb(alpha, beta, 1);
            TestANComb(alpha, beta, 2);

            cout << "sres *************************************" << endl;
            PolynomialQ f("(x-1)*(x-2)*(x-3)*(x-4)*(x-5)");
            PolynomialQ g("(x+1)*(x+2)*(x+3)*(x+4)");
            for (int k = 0; k <= 4; k++)
                cout << "Sres" << k << ": " << PolynomialQQ::sres
                    (f.getEx(), g.getEx(), k, PolynomialQ::GetVar()) << endl;
            g = PolynomialQ("(x-1)*(x+2)*(x+3)*(x+4)");
            for (int k = 0; k <= 4; k++)
                cout << "Sres" << k << ": " << PolynomialQQ::sres
                    (f.getEx(), g.getEx(), k, PolynomialQ::GetVar()) << endl;
            g = PolynomialQ("(x-1)*(x-2)*(x+3)*(x+4)");
            for (int k = 0; k <= 4; k++)
                cout << "Sres" << k << ": " << PolynomialQQ::sres
                    (f.getEx(), g.getEx(), k, PolynomialQ::GetVar()) << endl;
            g = PolynomialQ("(x-1)*(x-2)*(x-3)*(x+4)");
            for (int k = 0; k <= 4; k++)
                cout << "Sres" << k << ": " << PolynomialQQ::sres
                    (f.getEx(), g.getEx(), k, PolynomialQ::GetVar()) << endl;

            cout << endl;

            PolynomialQ p("672050176/29225227377515625-8044544/29225227377515625*x+16384/29225227377515625*x^2");
            cout << p << endl;

            PolynomialQ::vector v;
            p.addIrreducibleFactorsTo(v);
        }*/
        /*{
            cout << "gcdex ************************************" << endl;
            ex s,t,gcd;
            PolynomialQ f = PolynomialQ("x^3-1");
            PolynomialQ g = PolynomialQ("x^2-1");
            TestGCD(f,g);
            f = PolynomialQ("x^3+1");
            g = PolynomialQ("x^2 + 2*x + 1");
            TestGCD(f,g);
        }*/
/*
        Algebraic alpha(PolynomialQ("x^2-2"), IntervalQ(1, 2));
        Algebraic beta(PolynomialQ("x^2-3"), IntervalQ(1, 2));

        cout << alpha << endl;
        cout << beta << endl;
GiNaC::symbol v = alpha.getPolynomial().getVariable();
        cout << GiNaC::factor(v*v*v*v - 10*v*v +1) << endl;
        {
            cout << "****************" << endl;
            ex f = v*v - 5*v + 6;
            ex fac = factor(f);
            cout << ((f == fac) ? "same" : "diff") << endl;
            cout << "Is f an add?" << (is_a<add>(f) ? "Yes" : "No") << endl;
            cout << "Is f a mul?" << (is_a<mul>(f) ? "Yes" : "No") << endl;
            cout << "Is fac a add?" << (is_a<add>(fac) ? "Yes" : "No") << endl;
            cout << "Is fac a mul?" << (is_a<mul>(fac) ? "Yes" : "No") << endl;
            cout << f << endl;
            cout << fac << endl;
            ex g = ex((v-2)*(v-3)*(v-5)).expand();
            ex gfac = factor(g);
            cout << "g: " << g << endl;
            cout << "g factored: " << gfac << endl;
            cout << "Is g an add?" << (is_a<add>(g) ? "Yes" : "No") << endl;
            cout << "Is g a mul?" << (is_a<mul>(g) ? "Yes" : "No") << endl;
            cout << "Is gfac a add?" << (is_a<add>(gfac) ? "Yes" : "No") << endl;
            cout << "Is gfac a mul?" << (is_a<mul>(gfac) ? "Yes" : "No") << endl;
            cout << "****************" << endl;
        }*/
        /*{
            PolynomialQ p("(x-3)*(x-5)*(x-7)");
            cout << p << endl;
            vector<PolynomialQ> seq = PolynomialQ::sturmseq(p);
            cout << "Sturm sequence:" << endl;
            BOOST_FOREACH(const PolynomialQ & f, seq)
                cout << "    " << f << endl;
            cout << "sturm(-100,100): " << ex(PolynomialQ::sturm(seq,-100,100)) << endl;
            cout << "sturm(-6,6): " << ex(PolynomialQ::sturm(seq,-6,6)) << endl;
            cout << "sturm(-1,1): " << ex(PolynomialQ::sturm(seq,-1,1)) << endl;
            cout << "sturm(6,8): " << ex(PolynomialQ::sturm(seq,6,8)) << endl;
            cout << endl << endl;
            vector<Algebraic> alphas =
                PolynomialQ::FindRootsOfIrreducible(PolynomialQ("x*x-2"));
            cout << "Number of alphas: " << alphas.size() << endl;
            BOOST_FOREACH(const Algebraic & alpha, alphas)
                cout << alpha << endl;
            alphas =
                PolynomialQ::FindRootsOfIrreducible(PolynomialQ("x*x-3"));
            cout << "Number of alphas: " << alphas.size() << endl;
            BOOST_FOREACH(const Algebraic & alpha, alphas)
                cout << alpha << endl;
            alphas =
                PolynomialQ::FindRootsOfIrreducible(PolynomialQ("x*x-5"));
            cout << "Number of alphas: " << alphas.size() << endl;
            BOOST_FOREACH(const Algebraic & alpha, alphas)
                cout << alpha << endl;
            alphas =
                PolynomialQ::FindRootsOfIrreducible(PolynomialQ("x^2+x+10"));
            cout << "Number of alphas: " << alphas.size() << endl;
            BOOST_FOREACH(const Algebraic & alpha, alphas)
                cout << alpha << endl;
            cout << "**** FindRoots Test **************************" << endl;
            vector<PolynomialQ> plist;
            plist.push_back(PolynomialQ("(x-3)*(x-7)*(x-11)"));
            alphas = PolynomialQ::FindRoots(plist);
            cout << "Number of alphas: " << alphas.size() << endl;
            BOOST_FOREACH(const Algebraic & alpha, alphas)
                cout << alpha << endl;
            plist.push_back(PolynomialQ("x^2-2"));
            alphas = PolynomialQ::FindRoots(plist);
            cout << "Number of alphas: " << alphas.size() << endl;
            BOOST_FOREACH(const Algebraic & alpha, alphas)
                cout << alpha << endl;
        }*/
/*
        boost::tuple<Algebraic, PolynomialQ, PolynomialQ> t =
            PolynomialQQ::Simple(alpha, beta);

        cout << t.get<0>() << endl;
        cout << t.get<1>() << endl;
        cout << t.get<2>() << endl;
*/
//        Algebraic gamma = PolynomialQQ::ANComb(alpha, beta, 1);

//        cout << gamma << endl;
    }
    catch (exception & e)
    {
        cout << "Exception thrown in CAD code: " << e.what() << endl;
    }
/*
    PolynomialQ p(std::string("1+x+x^2+7/2*x^5")),
                f(std::string("x+2")),g(std::string("x+3"));

    std::cout << p;
    std::cout << f;
    std::cout << g;
    std::cout << f*g;
    std::cout << (f*g)/PolynomialQ(std::string("x+3"));

    PolynomialQ::TestClass();

    PolynomialQ //p("(x-1)(x-2)"),
                q("x^2-3*x+2");

    std::cout << std::endl << q;

    std::vector<PolynomialQ> factors = q.getIrreducibleFactors();
    std::cout << "Number of factors: " << factors.size() << std::endl;

    BOOST_FOREACH(const PolynomialQ & p, factors)
        std::cout << p;

    PolynomialQ p(PolynomialQ("x-5")*PolynomialQ("x-7")*PolynomialQ("x-11"));

    std::cout << std::endl << p;

    factors = p.getIrreducibleFactors();
    std::cout << "Number of factors: " << factors.size() << std::endl;

    BOOST_FOREACH(const PolynomialQ & poly, factors)
        std::cout << poly;

    PolynomialQ p1(PolynomialQ("x-5")*PolynomialQ("x-5")*PolynomialQ("x-11"));

    std::cout << std::endl << p1;

    factors = p1.getIrreducibleFactors();
    std::cout << "Number of factors: " << factors.size() << std::endl;

    BOOST_FOREACH(const PolynomialQ & poly, factors)
        std::cout << poly;

    PolynomialQ p2(PolynomialQ("x^2-2")*PolynomialQ("x-5")*PolynomialQ("x-11"));

    std::cout << std::endl << p2;

    factors = p2.getIrreducibleFactors();
    std::cout << "Number of factors: " << factors.size() << std::endl;

    BOOST_FOREACH(const PolynomialQ & poly, factors)
        std::cout << poly;
*/
    //std::cin.get();

    return 0;
}
