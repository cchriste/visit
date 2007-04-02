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

// ************************************************************************* //
//                                 avtExtents.h                              //
// ************************************************************************* //

#ifndef AVT_EXTENTS_H
#define AVT_EXTENTS_H

#include <pipeline_exports.h>

#include <visitstream.h>

class   avtDataObjectString;
class   avtDataObjectWriter;
class   vtkMatrix4x4;


// ****************************************************************************
//  Class: avtExtents
//
//  Purpose:
//      Manages a set of extents.
//
//      The extents are an array of 'dimension' doubles of min/max pairs.
//      Even numbered indices are for min values while odd numbered indices
//      are for max values.
//
//  Programmer: Hank Childs
//  Creation:   September 4, 2001
//
//  Modifications:
//
//    Jeremy Meredith, Mon Sep 24 14:18:03 PDT 2001
//    Added Transform function.
//
//    Hank Childs, Wed Feb 25 09:00:17 PST 2004
//    Added Print method.
//
// ****************************************************************************

class PIPELINE_API avtExtents
{
  public:
                          avtExtents(int);
                          avtExtents(const avtExtents &);
    virtual              ~avtExtents();

    void                  Print(ostream &);

    bool                  HasExtents(void);
    int                   GetDimension(void)  { return dimension; };

    void                  CopyTo(double *);
    void                  Set(const double *);
    void                  Clear(void);

    const avtExtents     &operator=(const avtExtents &);

    void                  Merge(const avtExtents &);
    void                  Merge(const double *);

    void                  Write(avtDataObjectString &,
                                const avtDataObjectWriter *);
    int                   Read(char *);

    void                  Transform(vtkMatrix4x4 *);

  protected:
    int                   dimension;
    double               *extents;
};


#endif


