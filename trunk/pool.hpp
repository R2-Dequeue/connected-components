/*!
 * \file
 * \author Chris de Pujo
 */

#ifndef __SIMPLEPOOL__
#define __SIMPLEPOOL__

#include <new>          // std::bad_alloc
#include <limits>       // std::numeric_limits
#include <cassert>      // assert
#include <stdexcept>    // std::invalid_argument

#include <boost/noncopyable.hpp>
#include <boost/scoped_array.hpp>

#define static_assert(x) \
    do { typedef char ____static_assert_check[ (x) ? 1 : -1 ]; } while (false)

class SimplePool : public boost::noncopyable
{
public:

    static char * Align(const char * ptr, size_t width);

    template <typename T>
	friend class PoolAllocator; // Allows for the most minimal public interface.

protected:

    SimplePool(char * mem, bool takeOwnership, size_t bsize, size_t n) :
        block_size(bsize),
        data(takeOwnership ? mem : NULL),
        head_element(SimplePool::Align(mem, bsize)),
        end_element(head_element + (bsize*n)) {}

    void * allocate();
    void deallocate(void * mem);

	const size_t                block_size;
	boost::scoped_array<char>   data;           // const?

	char *                      head_element;
	const char *                end_element;

	// Really need pointer to the begining for error checking in deallocate.
};

inline char * SimplePool::Align(const char * ptr, size_t width)
{
    typedef size_t integral_type;

	static_assert(std::numeric_limits<integral_type>::is_integer == true);
	static_assert(std::numeric_limits<integral_type>::is_signed == false);
	static_assert(sizeof(integral_type) == sizeof(void *));
	static_assert(sizeof(integral_type) == sizeof(char *));
	static_assert(sizeof(char) == 1);

	assert(ptr != NULL); // throw ?

	integral_type p = reinterpret_cast<integral_type>(ptr);

	if (p % width != 0)
        p += width - (p % width);

    return reinterpret_cast<char *>(p);
}

inline void * SimplePool::allocate()
{
    if (head_element == end_element)
        throw std::bad_alloc();

    void * ret = this->head_element;

    this->head_element += this->block_size;

    return ret;
}

inline void SimplePool::deallocate(void * mem)
{
    if (reinterpret_cast<char *>(mem) >= end_element)
        throw std::invalid_argument("Attempt to deallocate memory not in "
                                    "this pool.");
}

#endif // __SIMPLEPOOL__
