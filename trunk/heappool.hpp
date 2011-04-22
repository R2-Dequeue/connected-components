/*!
 * \file
 * \author Chris de Pujo
 */

#ifndef __HEAPPOOL__
#define __HEAPPOOL__

#include "pool.hpp"

// 1) Does it make sense to allow pools to be copied?
// 2) Does it make sense to have a default constructor?
// 3) Does assignment make sense?

class HeapPool : public SimplePool
{
public:

    // new throws? What happens if 'new' is passed 0?
    explicit HeapPool(size_t bsize, size_t n) :
        SimplePool(new char[bsize*(n + 1)], true, bsize, n) {}
};

#endif // __HEAPPOOL__
