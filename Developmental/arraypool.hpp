/*!
 * \file
 * \author Chris de Pujo
 */

#ifndef __ARRAYPOOL__
#define __ARRAYPOOL__

#include "pool.hpp"

template <size_t BlockSize, size_t StackCount = 32>
class ArrayPool : public SimplePool
{
public:

    explicit ArrayPool(size_t n) :
        SimplePool((n > StackCount) ? new char[BlockSize * (n + 1)] : carray,
                   (n > StackCount),
                   BlockSize,
                   n) {}

private:

    // Most likely should use HeapPool in this case.
    void * operator new(size_t n);

    char                        carray[BlockSize * (StackCount + 1)];
};

#endif // __ARRAYPOOL__
