/*!
 * \file
 * \author Chris de Pujo
 */

#ifndef __CAD__
#define __CAD__

#include <ginac/ginac.h>

#include "polynomialqq.h"
#include "algebraic.h"

#include <vector>
#include <utility>
#include <string>

#include <boost/numeric/ublas/symmetric.hpp>

//! The type used for the connectivity matrix.
typedef boost::numeric::ublas::symmetric_matrix<unsigned char> uBitMatrix;

typedef pair<Algebraic, Algebraic> Point; //!< Represents a point in the plane.
/*!
 * \brief Stores pairs of indices into the CAD.
 *
 * \detail Indices start at 0, not at 1.
 */
typedef pair<unsigned int, unsigned int> CellIndex;

/*!
 * \brief The cyclical algebraic decomposition of a set of polynomials.
 *
 * \detail
 */
class CAD
{
private:

    class Sample
    {
    public:
        Sample(const Algebraic & b, const std::vector<unsigned char> & s)
            : y(b), signs(s) {};

        Algebraic y;
        std::vector<char> signs; // -1, 0, or 1
    }

    std::vector<PolynomialQQ> F; //!< The set of polynomials dividing the plane.
    std::vector<Algebraic> samples;
    std::vector<std::vector<Sample>> stacks; //!< The CAD proper.
    /*!
     * \brief Stores the connectivity information.
     * \detail The indices run from 0 to n-1. Symmetric, transitive, and
     *         reflexive (m_i,i == 1). Entries are either 0 or 1.
     * \todo Use bitset-like class for better storage.
     */
    uBitMatrix cmatrix;

public:

    //CAD();
    CAD(const std::list<std::string> & F);

    bool Connectivity(const Point & p1, const Point & p2) const;

    //! Returns the index of the cell containing p.
    CellIndex Cell(const Point & p) const;
    unsigned int CellNumber(const CellIndex & i) const;

    /*!
     * \detail This member is public but shouldn't really be published.
     */
    bool Invariant() const;

private:

    //! Calculates the connectivity matrix based on the private CAD.
    void ConnectivityMatrix();
    CellIndex BranchCount(const CellIndex & ci);
    void AdjacencyLeft(const unsigned int k);
    void AdjacencyRight(const unsigned int k);

    std::vector<PolynomialQ> Project(const std::vector<PolynomialQQ> & F)
};

#endif // __CAD__
