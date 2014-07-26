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
//                           ExpressionException.h                           //
// ************************************************************************* //

#ifndef EXPRESSION_EXCEPTION_H
#define EXPRESSION_EXCEPTION_H
#include <avtexception_exports.h>

#include <PipelineException.h>


// ****************************************************************************
//  Class: ExpressionException
//
//  Purpose:
//      Thrown when an expression fails to evaluate or an expression filter
//      fails.  Example: taking the log of a negative number or trying to
//      add variables that don't conform to the same mesh.
//
//
//  Programmer: Sean Ahern
//  Creation:   Fri Mar 22 13:21:20 PST 2002
//
//  Modifications:
//    Brad Whitlock, Fri Jun 28 13:26:48 PST 2002
//    Added windows api.
//
//    Sean Ahern, Mon Dec 10 09:53:05 EST 2007
//    Required the expression name.
//
//    Jeremy Meredith, December 18, 2007
//    Added a const char *name constructor to handle the case where the name
//    is NULL.
//
// **************************************************************************** 

class AVTEXCEPTION_API ExpressionException : public PipelineException
{
  public:
                          ExpressionException(std::string name, std::string reason);
                          ExpressionException(const char *name, std::string reason);
    virtual              ~ExpressionException() VISIT_THROW_NOTHING {;};
};



// ****************************************************************************
//  Class: ExpressionParseException
//
//  Purpose:
//      Thrown when an expression fails to parse.
//
//  Programmer: Sean Ahern
//  Creation:   Mon Dec 10 13:13:08 EST 2007
//
//  Modifications:
// **************************************************************************** 

class AVTEXCEPTION_API ExpressionParseException : public PipelineException
{
  public:
                          ExpressionParseException(std::string reason);
    virtual              ~ExpressionParseException() VISIT_THROW_NOTHING {;};
};


#endif


