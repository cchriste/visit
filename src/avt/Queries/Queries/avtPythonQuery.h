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
//                             avtPythonQuery.h                              //
// ************************************************************************* //

#ifndef AVT_PYTHON_QUERY_H
#define AVT_PYTHON_QUERY_H

#include <query_exports.h>

#include <avtDataObjectQuery.h>
#include <avtDatasetSink.h>
#include <QueryAttributes.h>

#include <string>
#include <vector>

class avtPythonFilterEnvironment;


// ****************************************************************************
//  Class: avtPythonQuery
//
//  Purpose:
//      Interface to python queries.
//
//  Programmer: Cyrus Harrison
//  Creation:   Tue Feb  9 15:15:41 PST 2010
//
//  Modifications:
//   Cyrus Harrison, Tue Sep 21 11:14:21 PDT 2010
//   Added SetPythonArgs()
//
//   Cyrus Harrison, Wed Jan 12 11:32:42 PST 2011
//   Added queryType & queryDescription members.
//
//   Cyrus Harrison, Fri Mar 30 13:51:24 PDT 2012
//   Convert python query filter to use new query params infrastructure.
//
// ****************************************************************************

class QUERY_API avtPythonQuery :  public avtDataObjectQuery,
                                  public avtDatasetSink
{
  public:
                                avtPythonQuery();
    virtual                    ~avtPythonQuery();
    void                        CleanUp();
    virtual void                SetInputParams(const MapNode &);
    void                        SetPythonScript(const std::string &py_script);
    void                        SetPythonArgs(const std::string &py_args);

    virtual const char         *GetType(void);
    virtual const char         *GetDescription(void);
    virtual void                PerformQuery(QueryAttributes *qa);

    void                        SetVariableNames(const std::vector<std::string> &vars)
                                    { varNames = vars;}

    virtual void                GetSecondaryVariables(std::vector<std::string> &res);

    void                        SetResultMessage(const std::string &msg)
                                 { resultMessage = msg; };
    void                        SetResultValues(const doubleVector &d)
                                 { resultValues = d;}
    void                        SetXmlResult(const std::string &xml)
                                 { resultXml= xml; };

    virtual std::string         GetResultMessage(void)
                                 { return resultMessage;}

 private:

    virtual void                PreExecute(void);
    virtual void                PostExecute(void);
    virtual void                UpdateContract(void);

    virtual void                Execute();

    QueryAttributes             queryAtts;

    avtPythonFilterEnvironment *pyEnv;
    std::string                 pyScript;
    std::string                 pyArgs;

    stringVector                varNames;

    std::string                 resultMessage;
    doubleVector                resultValues;
    std::string                 resultXml;

    std::string                 queryType;
    std::string                 queryDescription;
};


#endif


