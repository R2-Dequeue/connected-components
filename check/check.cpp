#include <ginac/ginac.h>

#include "../polynomialq.hpp"
#include "../polynomialqq.hpp"

//#define BOOST_TEST_NO_MAIN
#define BOOST_TEST_MODULE Connected Components
#include <boost/test/included/unit_test.hpp>

BOOST_AUTO_TEST_SUITE(Class_PolynomialQ)

BOOST_AUTO_TEST_CASE(Constructor_Tests)
{
    PolynomialQ pnull;
    BOOST_CHECK_EQUAL(pnull.degree(), 0);
    //BOOST_CHECK(pnull.isMonic());
    BOOST_CHECK(pnull.isZero());
    BOOST_CHECK(pnull.isIrreducible());
    BOOST_CHECK(pnull.isConstant());
    BOOST_CHECK_EQUAL(pnull.getEx(), (GiNaC::ex)0);
    //BOOST_CHECK_EQUAL(pnull, pnull.getMonic());
    BOOST_CHECK_EQUAL(pnull, pnull.getDerivative());
    BOOST_CHECK(pnull.getIrreducibleFactors<PolynomialQ::vector>()->empty());
    BOOST_CHECK_EQUAL(pnull.getCoeff(0), 0);
    BOOST_CHECK_EQUAL(pnull.getCoeff(1), 0);

    BOOST_CHECK_EQUAL(pnull.differentiate(), PolynomialQ(GiNaC::numeric(0)));
    BOOST_CHECK_EQUAL(pnull.getEx(), GiNaC::ex(0));

    PolynomialQ nmonic("5*x^2-25*x+30");
    BOOST_CHECK_EQUAL(nmonic.degree(), 2);
    BOOST_CHECK(!nmonic.isMonic());
    BOOST_CHECK(!nmonic.isZero());
    BOOST_CHECK(!nmonic.isIrreducible());
    BOOST_CHECK(!nmonic.isConstant());
    BOOST_CHECK_EQUAL(nmonic.getEx(), 5*pow(nmonic.getVariable(), GiNaC::ex(2))-25*nmonic.getVariable()+30);
    BOOST_CHECK_EQUAL(nmonic.getDerivative(), PolynomialQ("10*x-25"));
    BOOST_CHECK(nmonic.getIrreducibleFactors<PolynomialQ::vector>()->size() == 2);

    PolynomialQ nm2("5*(x-2)*(x-3)");
    BOOST_CHECK_EQUAL(nmonic, nm2);
}

BOOST_AUTO_TEST_CASE(Derivative_Tests)
{
    PolynomialQ p("x^2+1");
    BOOST_CHECK_EQUAL(p.getDerivative().getEx(), GiNaC::ex(p.getVariable()*2));
    BOOST_CHECK_EQUAL(p.getDerivative().getDerivative().getEx(), GiNaC::ex(2));
    BOOST_CHECK_EQUAL(p.getDerivative().getDerivative().getDerivative().getEx(), GiNaC::ex(0));
    BOOST_CHECK_EQUAL(p.getDerivative().getDerivative().getDerivative().getDerivative().getEx(), GiNaC::ex(0));
    BOOST_CHECK_EQUAL(p.getDerivative().getDerivative().getDerivative().getDerivative().getDerivative().getEx(), GiNaC::ex(0));
    BOOST_CHECK_EQUAL(p.getDerivative().getDerivative().getDerivative().getDerivative().getDerivative().getDerivative().getEx(), GiNaC::ex(0));
    p.differentiate();
    BOOST_CHECK_EQUAL(p, GiNaC::ex(p.getVariable()*2));
    p.differentiate();
    BOOST_CHECK_EQUAL(p, GiNaC::ex(2));
    p.differentiate();
    BOOST_CHECK_EQUAL(p, GiNaC::ex(0));
    p.differentiate();
    BOOST_CHECK_EQUAL(p, GiNaC::ex(0));
    p.differentiate();
    BOOST_CHECK_EQUAL(p, GiNaC::ex(0));
    p.differentiate();
    BOOST_CHECK_EQUAL(p, GiNaC::ex(0));
}

BOOST_AUTO_TEST_CASE(old_tests)
{
    BOOST_CHECK_EQUAL(2, 2);
}

BOOST_AUTO_TEST_CASE(gcd_tests)
{
    GiNaC::ex c1, c2;

    {PolynomialQ p("x^3-1"), q("x^2-1"), ans("-1+x"), ansc1("1"), ansc2("-x");
    BOOST_CHECK_EQUAL(PolynomialQQ::gcdex(p.getEx(), q.getEx(), p.GetVar(),
                                          c1, c2), ans.getEx());
    BOOST_CHECK_EQUAL(c1, ansc1.getEx());
    BOOST_CHECK_EQUAL(c2, ansc2.getEx());}

    {PolynomialQ p("x^2 + 7*x + 6"), q("x^2 - 5*x - 6"), ans = PolynomialQ("x+1");
    BOOST_CHECK_EQUAL(PolynomialQQ::gcdex(p.getEx(), q.getEx(), p.GetVar(),
                                          c1, c2), ans.getEx());}

    {PolynomialQ p("2*x^5 - 2*x"), q("(x^2 - 1)^2"), ans("x^2 - 1"),
                 ansc1("x / 4"), ansc2("(-4 - 2*x^2)/4");
    BOOST_CHECK_EQUAL(PolynomialQQ::gcdex(p.getEx(), q.getEx(), p.GetVar(),
                                          c1, c2), ans.getEx());
    BOOST_CHECK_EQUAL(c1, ansc1.getEx());
    BOOST_CHECK_EQUAL(c2, ansc2.getEx());}

    {PolynomialQ p("x^3 + 1"), q("x^2 + 2*x + 1"), ans("x + 1"),
                 ansc1("1/3"), ansc2("2/3 - x/3");
    BOOST_CHECK_EQUAL(PolynomialQQ::gcdex(p.getEx(), q.getEx(), p.GetVar(),
                                          c1, c2), ans.getEx());
    BOOST_CHECK_EQUAL(c1, ansc1.getEx());
    BOOST_CHECK_EQUAL(c2, ansc2.getEx());} // x^3 + 1, x^2 + 2*x + 1
}

BOOST_AUTO_TEST_SUITE_END()
