/*!
 * \file
 * \author Chris de Pujo
 */

#if HAVE_ALLOCA_H
#    include <alloca.h>
#elif defined _MSC_VER || defined __BORLANDC__ || defined __MINGW32__
#    include <malloc.h>
#elif defined __GNUC__
#    define alloca __builtin_alloca
#elif defined _AIX
#    define alloca __alloca
#else
#    include <stddef.h>
#    ifdef  __cplusplus
extern "C"
#    endif
void * alloca(size_t);
#endif

#include <cassert>      // assert

#include <boost/tuple/tuple.hpp>

#include "pool.hpp"

namespace aBaBaBaB {
// formacroonly, formacrosonly, macro, macrohelp, macroreserved,
// reservedformacro, macroonly, formacrouse, formacrouseonly

size_t count = 0; // so so ghetto...
size_t width = 0;

template <size_t d>
inline size_t ____A(size_t n)
{
    static_assert(d != 0);
    static_assert( (d & (d-1)) == 0 ); // ==> is a power of 2.

    assert(n > 0);      // throw?

    width = d;
    count = n;

    return (d * (n + 1));
}

inline boost::tuple<char *, size_t, size_t> ____B(void * ptr)
{
    return boost::make_tuple(reinterpret_cast<char *>(ptr), width, count);
}

}

#define salloc(C, N) \
    aBaBaBaB::____B( alloca( aBaBaBaB::____A<(C)>( (N) ) ) )

template <typename T>
class PoolAllocator;

// Different names:
// TempPool - emphasizes the transient nature of stack memory;
//            the name implies the what its used for (temporary memory)
//            and not where it is stored (the stack).

class StackPool : public SimplePool
{
public:

    /*!
     * \detail \c t should always be the return value from \c salloc.
     */
    //                                    data,   bsize,  count
    explicit StackPool(const boost::tuple<char *, size_t, size_t> & t) :
        SimplePool(t.get<0>(), false, t.get<1>(), t.get<2>()) {}

private:

	void * operator new(size_t n);

	// Make '&' operator private?
};

// Namespaces: pool, pl, aloc, allo, alloc, mem, memory

////////////////////////////////////////////////////////////////////////////////

/*#define salloc2(B, N)

#define MakePool(B, N) \
    StackPool pool; \
    do { \
        size_t count = aBaBaBaB::____C<(B)>( pool, (N) ); \
        void * ptr = alloc(count); \
        aBaBaBaB::____D(pool); \
    } while (false)*/

// stack_malloc, smalloc, sallocate, stack_allocate, stack_new, snew, salloc

// StackBuffer buffer(salloc(32, r.degree()));
// StackPool pool(salloc(32, r.degree()));

// ____StackPoolPointerHelp
// ____StackPoolAlign
// ____Align

// 1)
// StackPool<PolynomialQ> polys(stack_malloc(PolynomialQ, 32));

// buffer, array, vector, tuple, sequence, pool

////////////////////
// Could make a macro to be called only once per function (really, per-scope)
// that would call stack_malloc and create a StackPool object with a standard
// name that would clash if the macro was repeatedly invoked.
////////////////////

/*
#include <limits>
static_assert(std::numeric_limits<N>::is_integer == true);
static_assert(std::numeric_limits<N>::is_signed == false);
static_assert(align_size != 0);
static_assert(align_size & (align_size-1) == 0); // ==> is a power of 2.


long
unsigned
unsigned int
size_t
std::size_t

max_align_t // C++0x, in <stddef> (or maybe <cstddef>)

typeof
typeid
decltype

Ways to keep the macro as one statement:
	({ ...... })								// GCC extension
	do { ...... } while (false)
	comma operator
*/
