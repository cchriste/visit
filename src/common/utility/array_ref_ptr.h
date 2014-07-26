/*****************************************************************************
*
* Copyright (c) 2000 - 2014, Lawrence Livermore National Security, LLC
* Produced at the Lawrence Livermore National Laboratory
* LLNL-CODE-442911
* All rights reserved.
*
* This file is  part of VisIt. For  details, see https://visit.llnl.gov/.  The
* full copyright notice is contained in the file COPYRIGHT located at the root
* of the VisIt distribution or at http://www.llnl.gov/visit/copyright.html.
*
* Redistribution  and  use  in  source  and  binary  forms,  with  or  without
* modification, are permitted provided that the following conditions are met:
*
*  - Redistributions of  source code must  retain the above  copyright notice,
*    this list of conditions and the disclaimer below.
*  - Redistributions in binary form must reproduce the above copyright notice,
*    this  list of  conditions  and  the  disclaimer (as noted below)  in  the
*    documentation and/or other materials provided with the distribution.
*  - Neither the name of  the LLNS/LLNL nor the names of  its contributors may
*    be used to endorse or promote products derived from this software without
*    specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT  HOLDERS AND CONTRIBUTORS "AS IS"
* AND ANY EXPRESS OR  IMPLIED WARRANTIES, INCLUDING,  BUT NOT  LIMITED TO, THE
* IMPLIED WARRANTIES OF MERCHANTABILITY AND  FITNESS FOR A PARTICULAR  PURPOSE
* ARE  DISCLAIMED. IN  NO EVENT  SHALL LAWRENCE  LIVERMORE NATIONAL  SECURITY,
* LLC, THE  U.S.  DEPARTMENT OF  ENERGY  OR  CONTRIBUTORS BE  LIABLE  FOR  ANY
* DIRECT,  INDIRECT,   INCIDENTAL,   SPECIAL,   EXEMPLARY,  OR   CONSEQUENTIAL
* DAMAGES (INCLUDING, BUT NOT  LIMITED TO, PROCUREMENT OF  SUBSTITUTE GOODS OR
* SERVICES; LOSS OF  USE, DATA, OR PROFITS; OR  BUSINESS INTERRUPTION) HOWEVER
* CAUSED  AND  ON  ANY  THEORY  OF  LIABILITY,  WHETHER  IN  CONTRACT,  STRICT
* LIABILITY, OR TORT  (INCLUDING NEGLIGENCE OR OTHERWISE)  ARISING IN ANY  WAY
* OUT OF THE  USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
* DAMAGE.
*
*****************************************************************************/

#ifndef ARRAY_REF_PTR_H
#define ARRAY_REF_PTR_H

#include <visitstream.h>
#include "VisItAtomic.h"

// ****************************************************************************
//  Class: array_ref_ptr
//
//  Purpose:
//      Reference counted pointer for array variables.
//      This is similar to the STL auto_ptr, but it is actually useful.
//
//  Notes:      Copied from ref_ptr class.
//
//  Programmer: Hank Childs
//  Creation:   September 17, 2001
//
//  Modifications:
//
//    Tom Fogal, Mon Jun  9 14:30:36 EDT 2008
//    Added [] operator, for direct array access.
//
// ****************************************************************************

template <class T>
class array_ref_ptr
{
 public:
    array_ref_ptr();
    array_ref_ptr(T *p_, VISIT_ATOMIC_TYPE * = NULL);
    array_ref_ptr(const array_ref_ptr<T> &rhs);
    ~array_ref_ptr();
    void operator=(T *rhs);
    void operator=(const array_ref_ptr<T> &rhs);
    T *operator*();
    const T *operator*() const;
    T *operator->();
    const T *operator->() const;
    T &operator[](size_t idx);
    const T &operator[](size_t idx) const;
    template <class S>  void CopyTo(array_ref_ptr<S> &);
    void Print(ostream&);
 private:
    T                 *p;
    VISIT_ATOMIC_TYPE *n;

    void   AddReference(T *, VISIT_ATOMIC_TYPE *);
    void   RemoveReference(void);
};

// ****************************************************************************
//  Inline Methods
// ****************************************************************************
template <class T>
inline
array_ref_ptr<T>::array_ref_ptr()
{
    AddReference(NULL, NULL);
}

template <class T>
inline
array_ref_ptr<T>::array_ref_ptr(T *p_, VISIT_ATOMIC_TYPE *n_)
{
    AddReference(p_, n_);
}

template <class T>
inline
array_ref_ptr<T>::array_ref_ptr(const array_ref_ptr<T> &rhs)
{
    AddReference(rhs.p, rhs.n);
}

template <class T>
inline
array_ref_ptr<T>::~array_ref_ptr()
{
    RemoveReference();
}

template <class T>
inline void
array_ref_ptr<T>::operator=(T *rhs)
{
    RemoveReference();
    AddReference(rhs, NULL);
}

template <class T>
inline void
array_ref_ptr<T>::operator=(const array_ref_ptr<T> &rhs)
{
    if (p == rhs.p && n == rhs.n)
        return;
    RemoveReference();
    AddReference(rhs.p, rhs.n);
}

template <class T>
inline T *
array_ref_ptr<T>::operator*()
{
    return p;
}

template <class T>
inline const T *
array_ref_ptr<T>::operator*() const
{
    return p;
}

template <class T>
inline T *
array_ref_ptr<T>::operator->()
{
    return p;
}

template <class T>
inline const T *
array_ref_ptr<T>::operator->() const
{
    return p;
}

template <class T>
inline T &
array_ref_ptr<T>::operator[](size_t idx)
{
    return p[idx];
}

template <class T>
const T &
array_ref_ptr<T>::operator[](size_t idx) const
{
    return p[idx];
}

template <class T>
void
array_ref_ptr<T>::AddReference(T *p_, VISIT_ATOMIC_TYPE *n_)
{
    p = p_;
    if (p != NULL)
    {
        if (n_ == NULL)
        {
            n = new VISIT_ATOMIC_TYPE(1); 
        }
        else
        {
            n = n_;
            AtomicInc(*n);
        }
    }
    else
    {
        n = NULL;
    }
}

template <class T>
void
array_ref_ptr<T>::RemoveReference(void)
{
    if (p)
    {
        VISIT_ATOMIC_TYPE a = AtomicDec(*n);
        if (a == 0)
        {
            delete [] p;
            delete n;
        }
    }
}

#if defined(_MSC_VER) && (_MSC_VER <= 1200) // MSVC 6
template <class T, class S>
#else
template <class T>
template <class S>
#endif
void
array_ref_ptr<T>::CopyTo(array_ref_ptr<S> &rhs)
{
    array_ref_ptr<S>  newone((S *) p, n);
    rhs = newone;
}

template <class T>
inline void
array_ref_ptr<T>::Print(ostream &out)
{
    out << "<" << *n << ", " << (void*) p << ">";
}

typedef array_ref_ptr<char> CharStrRef;

#endif

