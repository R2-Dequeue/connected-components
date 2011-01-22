/*!
 * \file
 * @author Chris de Pujo
 * @version 0.1
 *
 * \brief Polynomial class header file.
 */

#ifndef __QX__
#define __QX__

#include <list>
#include <string>

#include <ginac/ginac.h>

	// assign
	// readin
	// parse
	// parseString
	// stringin
	// fromstring

/*!
 * \brief The main class for univariate polynomials over the rationals.
 *
 * \detail The coefficients are all stored as unlimited precision rational
 *         numbers and the polynomials themselves should only be limited by
 *         available memory.
 */
class qx
{
private:

    static const symbol variable("x"); //!< The variable of the polynomial.
    ex polynomial; //!< The internal representation of the polynomial.

public:

	qx(); //!< The default constructor.
	//~qx(); //!< The destructor.

	//qx(qx const& p);
	//qx(string const& p);

    //! Parses a polynomial from a 'string'
	bool parseString(const std::string & input);

    // These are the method versions of the binary arithmetic operators.
	qx operator+(const qx & rhs) const; //!< Normal binary addition.
    qx operator-(const qx & rhs) const;
	qx operator*(const qx & rhs) const;
//	qx operator/(const qx & rhs) const;
//	qx operator%(const qx & rhs) const;

	qx & operator+=(const qx & rhs);
    qx & operator-=(const qx & rhs);
	qx & operator*=(const qx & rhs);
//	qx & operator/=(const qx & rhs);
//	qx & operator%=(const qx & rhs);

//  ^ operator

//	qx operator*(mpq_class const& rhs) const;
//	qx operator/(mpq_class const& rhs) const;

//	qx & operator*=(mpq_class const& rhs);
//	qx & operator/=(mpq_class const& rhs);

	//operator std::string() const;
	std::string toString() const;

	// qx & makeMonic();

	mpz_class degree() const; //!< The degree of the polynomial.
	bool monic() const; //!< True iff the leading coefficient is zero.
	bool zero() const; //!< True iff the polynomial is '0'.

	qx & diff(); //!< Differentiates the polynomial.
	//qx & diff(mpz_class const& n);

protected:
	bool Invariant() const;
};

#endif // __QX__
