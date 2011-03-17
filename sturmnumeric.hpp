/*!
 * \file
 * \author Chris de Pujo
 */

#ifndef __STURMNUMERIC__
#define __STURMNUMERIC__

#include <ginac/ginac.h>

class sturmNumeric
{
public:

    GiNaC::numeric  number;
    unsigned int    signChanges;
    bool            isDataValid;
    bool            isRoot;

    sturmNumeric() : isDataValid(false) {}
    sturmNumeric(const GiNaC::numeric & a) : number(a), isDataValid(false) {}

    void assignMedian(const sturmNumeric & l, const sturmNumeric & u);
};

inline void sturmNumeric::assignMedian(const sturmNumeric & l,
                                       const sturmNumeric & u)
{
    number = (l.number + u.number) / 2;
    isDataValid = false;
}

#endif // __STURMNUMERIC__
