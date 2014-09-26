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

// ************************************************************************* //
//                     avtExpressionDataTreeIterator.h                       //
// ************************************************************************* //

#ifndef AVT_EXPRESSION_DATA_TREE_ITERATOR_H
#define AVT_EXPRESSION_DATA_TREE_ITERATOR_H

#include <expression_exports.h>

#include <avtExpressionFilter.h>
#include <avtDataTreeIterator.h>

#include <string>


// ****************************************************************************
//  Class: avtExpressionDataTreeIterator
//
//  Purpose:
//      This is an abstract type that allows derived types to create 
//      expressions one VTK dataset at a time.
//
//  Notes:      The streaming functionality used to be part of 
//              avtExpressionFilter.  The creation date corresponds to when the
//              class was split.
//
//  Programmer: Hank Childs
//  Creation:   December 27, 2004
//
//  Modifications:
//
//    David Camp, Tue May 21 13:56:12 PDT 2013
//    Remove these variables because the worker threads will call
//    the ExecuteData function and they all can not share these variables.
//    I will change the code to pass the currentDomainsIndex into the
//    DeriveVariable function. I have found that no one uses the 
//    currentDomainsLabel variable, so I am not going to pass it.
//     std::string              currentDomainsLabel;
//     int                      currentDomainsIndex;
//
//    Eric Brugger, Wed Aug 20 16:24:51 PDT 2014
//    Modified the class to work with avtDataRepresentation.
//
// ****************************************************************************

class EXPRESSION_API avtExpressionDataTreeIterator 
                                    : virtual public avtDataTreeIterator, 
                                      virtual public avtExpressionFilter
{
  public:
                             avtExpressionDataTreeIterator();
    virtual                 ~avtExpressionDataTreeIterator();

  protected:
    virtual avtDataRepresentation *ExecuteData(avtDataRepresentation *);
    virtual vtkDataArray    *DeriveVariable(vtkDataSet *, int currentDomainsIndex) = 0;
};

#endif

