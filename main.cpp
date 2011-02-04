/*!
 * \file
 * \author Chris de Pujo
 */

#include <iostream>
#include <stdexcept>
#include <string>

#include <ginac/ginac.h>

//#include "cad.hpp"
#include "polynomialq.hpp"

#include <boost/foreach.hpp>

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
    BOOST_CHECK(pnull.getIrreducibleFactors().empty());
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
    BOOST_CHECK(nmonic.getIrreducibleFactors().size() == 2);

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

BOOST_AUTO_TEST_SUITE_END()
/*
int main()
{
    using namespace std;
    using namespace GiNaC;

    cout << "Hello world!" << endl;

    symbol x("x"), y("y");

    ex poly = 1*pow(x, 0) + 2*pow(x, 1) + 3*pow(x, 2) + 4*pow(x, 3);
    ex poly2 = (x - 1)*(x - 2)*(x - 3)*(x - 4)*(x - 5);

    poly.expand();
    poly2.expand();

    cout << poly << ((poly.is_polynomial(x)) ? " is" : " is not") << " a polynomial." << endl;
    cout << poly2 << ((poly2.is_polynomial(x)) ? " is" : " is not") << " a polynomial." << endl;

    cout << poly.expand() << endl;
    cout << poly2.expand() << endl;

    poly = expand(poly*poly2);

    cout << "poly*poly2 = " << poly << endl << endl;

    // Code for CAD stuff /////////////////////////////////////////////////////

    try
    {
        list<string> F;
        F.push_back("x^2+y^2-4");
        F.push_back("x+y");

        CAD mainCAD(F);
    }
    catch (exception & e)
    {
        cout << "Exception thrown in CAD code: " << e.what() << endl;
    }

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

    std::cin.get();

    return 0;
}*/
