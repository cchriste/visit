/*****************************************************************************
*
* Copyright (c) 2000 - 2006, The Regents of the University of California
* Produced at the Lawrence Livermore National Laboratory
* All rights reserved.
*
* This file is part of VisIt. For details, see http://www.llnl.gov/visit/. The
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
*    documentation and/or materials provided with the distribution.
*  - Neither the name of the UC/LLNL nor  the names of its contributors may be
*    used to  endorse or  promote products derived from  this software without
*    specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT  HOLDERS AND CONTRIBUTORS "AS IS"
* AND ANY EXPRESS OR  IMPLIED WARRANTIES, INCLUDING,  BUT NOT  LIMITED TO, THE
* IMPLIED WARRANTIES OF MERCHANTABILITY AND  FITNESS FOR A PARTICULAR  PURPOSE
* ARE  DISCLAIMED.  IN  NO  EVENT  SHALL  THE  REGENTS  OF  THE  UNIVERSITY OF
* CALIFORNIA, THE U.S.  DEPARTMENT  OF  ENERGY OR CONTRIBUTORS BE  LIABLE  FOR
* ANY  DIRECT,  INDIRECT,  INCIDENTAL,  SPECIAL,  EXEMPLARY,  OR CONSEQUENTIAL
* DAMAGES (INCLUDING, BUT NOT  LIMITED TO, PROCUREMENT OF  SUBSTITUTE GOODS OR
* SERVICES; LOSS OF  USE, DATA, OR PROFITS; OR  BUSINESS INTERRUPTION) HOWEVER
* CAUSED  AND  ON  ANY  THEORY  OF  LIABILITY,  WHETHER  IN  CONTRACT,  STRICT
* LIABILITY, OR TORT  (INCLUDING NEGLIGENCE OR OTHERWISE)  ARISING IN ANY  WAY
* OUT OF THE  USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
* DAMAGE.
*
*****************************************************************************/

// ----------------------------------------------------------------------------
// File:  Matrix.h
//
// Programmer: Jeremy Meredith
// Date:       August 11, 2003
// ----------------------------------------------------------------------------

#ifndef MATRIX_H
#define MATRIX_H

class Vector;
#include <visitstream.h>

// ****************************************************************************
//  Class:  Matrix
//
//  Purpose:
//    Encapsulation of a 4x4 matrix.
//
//  Programmer:  Jeremy Meredith
//  Creation:    April 10, 2001
//
// ****************************************************************************
class Matrix
{
  public:
    float m[4][4];
    float openglm[16];
  public:
    Matrix();
    Matrix(const Matrix&);

    // assignment operator
    void   operator=(const Matrix&);

    // multiply matrix*matrix
    Matrix operator*(const Matrix&) const;
    // transform point
    Vector operator*(const Vector&) const;
    // transform vector
    Vector operator^(const Vector&) const;

    void   Inverse();
    void   Transpose();

    // utility
    void   CreateIdentity();
    void   CreateZero();
    void   CreateTrackball(float,float, float,float);
    void   CreateTranslate(float, float, float);
    void   CreateRBT(const Vector&, const Vector&, const Vector&);
    void   CreateScale(float,float,float);
    void   CreateScale(float);
    void   CreatePerspectiveProjection(float,float, float, float);
    void   CreateOrthographicProjection(float, float,float, float);
    void   CreateView(const Vector&, const Vector&, const Vector&);

    float* GetOpenGLMatrix();

    // friends
    friend ostream& operator<<(ostream&,const Matrix&);
};

#endif
