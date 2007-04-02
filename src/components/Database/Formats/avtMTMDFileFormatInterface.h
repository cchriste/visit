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
//                        avtMTMDFileFormatInterface.h                       //
// ************************************************************************* //

#ifndef AVT_MTMD_FILE_FORMAT_INTERFACE_H
#define AVT_MTMD_FILE_FORMAT_INTERFACE_H

#include <database_exports.h>

#include <avtFileFormatInterface.h>
#include <avtMTMDFileFormat.h>

class avtIOInformation;

// ****************************************************************************
//  Class: avtMTMDFileFormatInterface
//
//  Purpose:
//      This is an implementation of avtFileFormatInterface for file formats
//      that have Multiple Timesteps and a Multiple Domain.  This is the basic
//      abstraction used by the generic database.  However, this class is
//      needed to enable those classes that don't use the basic abstraction
//      provided by the generic database to also fit in.
//
//  Programmer: Hank Childs
//  Creation:   April 4, 2003
//
//  Modifications:
//    Brad Whitlock, Wed May 14 09:24:22 PDT 2003
//    Added an optional timeState argument to SetDatabaseMetaData.
//
//    Mark C. Miller, Mon Feb 23 20:38:47 PST 2004
//    Added ActivateTimestep method
//
//    Mark C. Miller, Tue May 17 18:48:38 PDT 2005
//    Added bool arg, forceReadAllCyclesTimes, to SetDatabaseMetaData
//
//    Mark C. Miller, Tue May 31 20:12:42 PDT 2005
//    Added method SetCycleTimeInDatabaseMetaData
// ****************************************************************************

class DATABASE_API avtMTMDFileFormatInterface : public avtFileFormatInterface
{
  public:
                    avtMTMDFileFormatInterface(avtMTMDFileFormat *);
    virtual        ~avtMTMDFileFormatInterface();

    virtual vtkDataSet     *GetMesh(int, int, const char *);
    virtual vtkDataArray   *GetVar(int, int, const char *);
    virtual vtkDataArray   *GetVectorVar(int, int, const char *);
    virtual void           *GetAuxiliaryData(const char *var, int, int,
                                             const char *type, void *args,
                                             DestructorFunction &);

    virtual const char     *GetFilename(int);
    virtual void            SetDatabaseMetaData(avtDatabaseMetaData *md,
                                int timeState = 0,
                                bool forceReadAllCyclesTimes = false);
    virtual void            SetCycleTimeInDatabaseMetaData(
                                avtDatabaseMetaData *md, int ts);
    virtual void            FreeUpResources(int, int);

    virtual void            ActivateTimestep(int ts);

    virtual void            PopulateIOInformation(int ts, avtIOInformation& ioInfo);

  protected:
    avtMTMDFileFormat      *format;

    virtual int             GetNumberOfFileFormats(void)
                              { return 1; };
    virtual avtFileFormat  *GetFormat(int n) const { return format; };
};


#endif


