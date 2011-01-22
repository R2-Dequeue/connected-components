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

int main()
{
    /*
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
    */

    // Code for CAD stuff /////////////////////////////////////////////////////

    /*
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
    */

    /*
    PolynomialQ p(std::string("1+x+x^2+7/2*x^5")),
                f(std::string("x+2")),g(std::string("x+3"));

    std::cout << p;
    std::cout << f;
    std::cout << g;
    std::cout << f*g;
    std::cout << (f*g)/PolynomialQ(std::string("x+3"));
    */

    PolynomialQ::TestClass();

    std::cin.get();

    return 0;
}
