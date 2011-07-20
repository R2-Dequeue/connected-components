/*!
 * \file
 * \author Chris de Pujo
 */

#ifndef __POOLALLOCATOR__
#define __POOLALLOCATOR__

#include <cstddef>      // size_t, ptrdiff_t
#include <limits>       // numeric_limits
#include <new>          // std::bad_alloc
#include <exception>    // std::exception

#include "stackmemory.hpp"

template <typename T>
class PoolAllocator
{
public:

    typedef     T           value_type;
    typedef     T *         pointer;
    typedef     const T *   const_pointer;
    typedef     T &         reference;
    typedef     const T &   const_reference;
    typedef     size_t      size_type;
    typedef     ptrdiff_t   difference_type;

    template <typename U>
    struct rebind { typedef PoolAllocator<U> other; };

    void check() const { assert(sizeof(T) <= pool->block_size); }

    PoolAllocator(SimplePool & p) : pool(&p) { check(); }
    template <typename U>
    PoolAllocator(const PoolAllocator<U> & a) : pool(a.pool) { check(); }

    pointer address(reference r) const;
    const_pointer address(const_reference r) const;

    size_type max_size() const;

    pointer allocate(size_type n);                  // const?
    void deallocate(pointer p, size_type n);        // const?

    void construct(pointer p, const_reference r);
    void destroy(pointer p);

    template <typename U>
    friend class PoolAllocator;

private:

    SimplePool * const pool;
};

template <typename T, typename U>
bool operator==(const PoolAllocator<T> & lhs, const PoolAllocator<U> & rhs);
template <typename T, typename U>
bool operator!=(const PoolAllocator<T> & lhs, const PoolAllocator<U> & rhs);

////////////////////////////////////////////////////////////////////////////////

template <typename T>
inline typename PoolAllocator<T>::pointer
    PoolAllocator<T>::address(reference r) const
{ return &r; }

template <typename T>
inline typename PoolAllocator<T>::const_pointer
    PoolAllocator<T>::address(const_reference r) const
{ return &r; }

template <typename T>
inline typename PoolAllocator<T>::size_type
    PoolAllocator<T>::max_size() const
{ return std::numeric_limits<size_type>::max() / sizeof(T); }

template <typename T>
inline typename PoolAllocator<T>::pointer
    PoolAllocator<T>::allocate(size_type n)
{
    if (n != 1)
        throw std::invalid_argument("allocated more than 1 element");
    else
        return reinterpret_cast<pointer>(pool->allocate());
}

template <typename T>
inline void PoolAllocator<T>::deallocate(pointer p, size_type n)
{
    if (n != 1)
        throw std::out_of_range("deallocated more than 1 element");
    else
        pool->deallocate(p);
}

template <typename T>
inline void PoolAllocator<T>::construct(pointer p, const_reference r)
{
    new (p) T(r);
}

template <typename T>
inline void PoolAllocator<T>::destroy(pointer p)
{
    p->~T();
}

template <typename T, typename U>
inline bool operator==(const PoolAllocator<T> & lhs, const PoolAllocator<U> & rhs)
{ return true; }

template <typename T, typename U>
inline bool operator!=(const PoolAllocator<T> & lhs, const PoolAllocator<U> & rhs)
{ return ( !(lhs == rhs) ); }

#endif // __POOLALLOCATOR__
