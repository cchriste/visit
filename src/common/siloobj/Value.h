/*****************************************************************************
*
* Copyright (c) 2000 - 2008, Lawrence Livermore National Security, LLC
* Produced at the Lawrence Livermore National Laboratory
* LLNL-CODE-400142
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
//                                 Value.h                                   //
// ************************************************************************* //

#ifndef VALUE_H
#define VALUE_H
#include <siloobj_exports.h>

#include <visitstream.h>
#include <silo.h>


// ****************************************************************************
//  Class: Value
// 
//  Purpose:
//      Keeps information about a value.
//
//  Data Members:
//      nDomains       -  The number of domains.
//      nVals          -  The number of variables in the vector this Value 
//                        represents.
//      totalEntries   -  The total number of elements for some arrays.  This 
//                        is nDomains*nVals.
//      entryNames     -  List of array names for the values for each domain. 
//                        The array is of size nDomains.
//      offsets        -  List of offsets into the values arrays for each 
//                        domain.  The array is of size totalEntries.
//      lengths        -  List of lengths of values for each domain.  The array
//                        is of size totalEntries.
//      files          -  The file number that each domain should go into.
//                        The array is of size nDomains.
//
//  Programmer: Hank Childs
//  Creation:   December 2, 1999
//
// ****************************************************************************

class SILOOBJ_API Value
{
  public:
                      Value();
    virtual          ~Value();
    
    char             *GetName() { return name; };
   
    void              PrintSelf(ostream &);

    void              Read(DBobject *, DBfile *);
    virtual void      Write(DBfile *);

  protected:
    int               nDomains;
    int               nVals;
    int               totalEntries;

    char            **entryNames;
    int              *offsets;
    int              *lengths;
 
    char             *name;
   
    // Class-Scoped Constants
  public:
    static char * const NAME;
    static char * const SILO_TYPE;
  protected:
    static char * const ARRAY_STRING;
    static int    const SILO_NUM_COMPONENTS;
    static char * const SILO_OBJ_NAME;

    static char * const SILO_N_DOMAINS_NAME;
    static char * const SILO_N_VALS_NAME;
    static char * const SILO_LENGTHS_NAME;
    static char * const SILO_ENTRY_NAMES_NAME;
    static char * const SILO_OFFSETS_NAME;
};


#endif


