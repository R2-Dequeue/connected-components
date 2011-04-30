/*!
 * \file
 * \author Chris de Pujo
 */

#include <cstddef>      // size_t, ptrdiff_t
#include <memory>

/*!
 * \details Useful mostly for classes with trivial default constructors as
 *          standard C arrays are used.
 */
template <typename T, size_t N = 32> // std::size_t ?
class StackArray
{
public:

    typedef     T           value_type;
    typedef     T *         iterator;
    typedef     const T *   const_iterator;
    typedef     T &         reference;
    typedef     const T &   const_reference;
    typedef     size_t      size_type;
    typedef     ptrdiff_t   difference_type;

    const size_t static_size = N;

    StackArray(size_type n);

    reference operator[](size_type pos);
    const_reference operator[](size_type pos) const;

    size_type size() const;

private:

    value_type          data[N];
    std::auto_ptr<T>    heap_data;            // scoped_array ?
    pointer             ptr;
    size_type           count;
};

StackArray::StackArray(size_type n) : count(n)
{
    // upper bound check for n?

    if (n > N)
    {
        heap_data = new T[n];
        ptr = heap_data.get();
    }
    else
        ptr = data;
}

inline reference StackArray::operator[](size_type pos)
{
    assert(pos < count && "out of range");
    return ptr[pos];
}

inline const_reference StackArray::operator[](size_type pos) const
{
    assert(pos < count && "out of range");
    return ptr[pos];
}

inline size_type StackArray::size() const
{
    return count;
}
