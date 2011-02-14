/*!
 * \file
 * \author Chris de Pujo
 */

#ifndef __CAD__
#define __CAD__

#include <vector>
#include <utility>
#include <string>

#include <boost/numeric/ublas/symmetric.hpp>

#include <ginac/ginac.h>

#include "polynomialq.hpp"
#include "polynomialqq.hpp"
#include "algebraic.hpp"

//! The type used for the connectivity matrix.
typedef boost::numeric::ublas::symmetric_matrix<unsigned char> uBitMatrix;
//! Represents a point in the plane.
typedef std::pair<Algebraic, Algebraic> Point;
/*!
 * \brief Stores pairs of indices into the CAD.
 * \detail Indices start at 0, not at 1.
 */
typedef std::pair<unsigned int, unsigned int> CellIndex;

/*!
 * \brief The cyclical algebraic decomposition of a set of polynomials.
 *
 * \detail This object presents an interface for determining the partitions
 *         of the plane made by two polynomials. For example [x+y, x-y]
 *         creates four partitions, and [x^2+y^2-1, x^2+y^2-4] creates three
 *         partitions.
 * \todo Maybe change to lists for internal private functions for fast
 *		 modifications.
 * \todo Maybe template methods so they can take/return various containers?
 */
class CAD
{
private:

    class Sample
    {
    public:
        Sample(const Algebraic & b, const std::vector<char> & s)
            : y(b), signs(s) {};

        Algebraic y;
        /*!
         * \todo Change this to a compact type like a pair of bitsets or a
         *       3-value equivalent of a bitset.
         */
        std::vector<char> signs; // -1, 0, or 1
    };

    std::vector<PolynomialQQ> F; //!< The set of polynomials dividing the plane.
    std::vector<Algebraic> samples; //!< x-coords of vertical dividing lines.
    std::vector< std::vector<Sample> > stacks; //!< The CAD proper.
    /*!
     * \brief Stores the connectivity information.
     * \detail The indices run from 0 to n-1. Symmetric, transitive, and
     *         reflexive (m_i,i == 1). Entries are either 0 or 1.
     * \todo Use bitset-like class for better storage.
     */
    uBitMatrix cMatrix;

public:

    CAD(const std::vector<std::string> & F);

	//! True if p1 and p2 are in the same component.
    bool Connectivity(const Point & p1, const Point & p2) const;

    //! Returns the index of the cell containing p.
    CellIndex Cell(const Point & p) const;
    //! Returns the cell's number according to the lexicographical ordering.
    unsigned int CellNumber(const CellIndex & i) const;

    void out() const;

    /*!
     * \brief Helper method for internal 'assert' checks.
     * \detail This method is public but shouldn't really be published.
     */
    bool Invariants() const;

    /*
     * 1) toString() method
     * 2) print() method to cout
     * 3) print(stream) method
     * 4) operator>> overload for use with streams
     */

private:

    //! Calculates the connectivity matrix based on the private CAD.
    void ConnectivityMatrix();

    CellIndex BranchCount(const CellIndex & ci);

    void AdjacencyLeft(const unsigned int k);
    void AdjacencyRight(const unsigned int k);

    static std::vector<PolynomialQ>
        Project(const std::vector<PolynomialQQ> & F);

	//! Inserts points between the passed vector.
    static std::vector<Algebraic>
        SamplePoints(const std::vector<Algebraic> & roots);

	//! Returns the roots of (f1*...*fn)(alpha, y).
    static std::vector<Algebraic>
        FindRoots2(const Algebraic & alpha, const std::vector<PolynomialQQ> & F);

	//! Internal helper method.
	inline static bool isEven(const unsigned int i) { return ((i & 1) == 0); }
    //! Internal helper method.
    inline static bool isOdd(const unsigned int i) { return ((i & 1) == 1); }
};

#endif // __CAD__
