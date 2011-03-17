/*!
 * \file
 * \author Chris de Pujo
 */

#ifndef __CAD__
#define __CAD__

#include <cassert>
#include <vector>
#include <utility>
#include <string>

#include <ginac/ginac.h>

#include "polynomialq.hpp"
#include "polynomialqq.hpp"
#include "algebraic.hpp"
#include "templatehelp.hpp"

//! Represents a point in the plane.
typedef std::pair<Algebraic, Algebraic> Point;
/*!
 * \brief Stores pairs of indices into the CAD.
 * \detail Indices start at 0, not at 1.
 */
typedef std::pair<unsigned int, unsigned int>   StacksIndex;
typedef unsigned int                            CellIndex;
typedef unsigned int                            PartitionIndex;
typedef unsigned int                            ComponentIndex;

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

    PolynomialQQ::vector F; //!< The set of polynomials dividing the plane.
    PolynomialQQ::vector irreducibles;
    Algebraic::vector alphas; //!< x-coords of vertical dividing lines.
    std::vector< std::vector<Sample> > stacks; //!< The CAD proper.
    std::vector< std::vector<CellIndex> > partitions;
    std::vector<PartitionIndex> components;
    /*!
     * \brief Stores the connectivity information.
     * \detail The indices run from 0 to n-1. Symmetric, transitive, and
     *         reflexive (m_i,i == 1). Entries are either 0 or 1.
     * \todo Use bitset-like class for better storage.
     */
    GiNaC::matrix cMatrix;

public:

    CAD(const std::vector<std::string> & F);
    CAD(const PolynomialQQ::vector & G);

    const PolynomialQQ::vector &                        GetF() const;
    const PolynomialQQ::vector &                        GetIrreducibles() const;
    const Algebraic::vector &                           GetAlphas() const;
    const std::vector< std::vector<Sample> > &          GetStacks() const;
    const std::vector< std::vector<unsigned int> > &    GetPartitions() const;

	//! True if p1 and p2 are in the same component.
    bool Connectivity(const Point & p1, const Point & p2) const;

    //! Returns the index of the cell containing p.
    StacksIndex Cell(const Point & p) const;
    //! Returns the cell's number according to the lexicographical ordering.
    unsigned int CellNumber(const StacksIndex & i) const;
    const Sample & RefCell(CellIndex i) const;

    bool areComponentsAdjacent(ComponentIndex i, ComponentIndex j) const;

    static void addIrreducibleProjectionTo(PolynomialQ::vector & P,
                                           const PolynomialQQ::vector & G);
    static PolynomialQ::vector Project(const PolynomialQQ::vector & F);

    /*!
     * \brief Prints a summary of *this to cout.
     */
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

    void ConstructorPart2();

    static bool SignsNonzero(const std::vector<char> & s);
    static bool SignsEqual(const std::vector<char> & S1,
                           const std::vector<char> & S2);
    static bool SignsDiffByOne(const std::vector<char> & s1,
                               const std::vector<char> & s2);

    //! Calculates the connectivity matrix based on the private CAD.
    void MakeConnectivityMatrix();
    void Partition();

    StacksIndex BranchCount(const StacksIndex & ci);
    std::vector< std::pair<StacksIndex, StacksIndex> >
        AdjacencyLeft(const unsigned int k);
    std::vector< std::pair<StacksIndex, StacksIndex> >
        AdjacencyRight(const unsigned int k);

    template <typename sample_type, typename root_type>
    void addSamplePointsTo(sample_type & S, const root_type & roots);
	//! Inserts points between the passed vector.
    static Algebraic::vector SamplePoints(const Algebraic::vector & roots);

	//! Returns the roots of (f1*...*fn)(alpha, y).
    static Algebraic::vector
        FindRoots2(const Algebraic & alpha, const PolynomialQQ::vector & F);

	//! Internal helper method.
	inline static bool isEven(const unsigned int i) { return ((i & 1) == 0); }
    //! Internal helper method.
    inline static bool isOdd(const unsigned int i) { return ((i & 1) == 1); }
};

///////////////////////////////////////////////////////////////////////////////

/*!
 * \throws parse_error Thrown by GiNaC if parsing of polynomials fails (inherits
 *		   invalid_argument).
 * \todo Check if doxygen collates all 'throws' statements from called
 *		 functions.
 */
inline CAD::CAD(const std::vector<std::string> & F)
{
    // Be careful with 'F' and 'this->F'.

    for (std::vector<std::string>::const_iterator f = F.begin(), e = F.end();
         f != e; ++f)
        this->F.push_back(PolynomialQQ(*f)); // throws on error

    this->ConstructorPart2();
}

inline CAD::CAD(const PolynomialQQ::vector & G) : F(G)
{
    this->ConstructorPart2();
}

inline const PolynomialQQ::vector & CAD::GetF() const
{
    return this->F;
}

inline const PolynomialQQ::vector & CAD::GetIrreducibles() const
{
    return this->irreducibles;
}

inline const Algebraic::vector & CAD::GetAlphas() const
{
    return this->alphas;
}

inline const std::vector< std::vector<CAD::Sample> > & CAD::GetStacks() const
{
    return this->stacks;
}

inline const std::vector< std::vector<unsigned int> > & CAD::GetPartitions() const
{
    return this->partitions;
}

/*!
 * \param roots Must have disjoint intervals.
 */
template <typename sample_type, typename root_type>
void CAD::addSamplePointsTo(sample_type & S, const root_type & roots)
{
    if (roots.size() == 0)
    {
        tmp::PushBack(S, Algebraic());
        return;
    }

    typename root_type::const_iterator root = roots.begin(),
                                       nextRoot,
                                       endRoot = roots.end();

    tmp::ReserveHelper(S, 2*roots.size() + 1);

    tmp::PushBack(S, Algebraic::MakeWideRational( root->lower() - 1 ));

    for ( ; ; ++root)
    {
        nextRoot = root; // This block is necessary b/c root is not
        ++nextRoot;      // necessarily a RandomAccessIterator.
        if (nextRoot == endRoot)
            break;

        tmp::PushBack(S, *root);

    	GiNaC::numeric s( (root->upper() + nextRoot->lower())/2 );
    	tmp::PushBack(S, Algebraic::MakeWideRational(s));

    	assert(*root < *nextRoot);
    }

    tmp::PushBack(S, *root);
    tmp::PushBack(S, Algebraic::MakeWideRational( root->upper() + 1 ));
}

#endif // __CAD__
