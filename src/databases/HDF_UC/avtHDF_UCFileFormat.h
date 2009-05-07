/*****************************************************************************
*
* Copyright (c) 2000 - 2007, The Regents of the University of California
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
//                           avtHDF_UCFileFormat.h                           //
// ************************************************************************* //

#ifndef AVT_HDF_UC_FILE_FORMAT_H
#define AVT_HDF_UC_FILE_FORMAT_H

#include <avtMTSDFileFormat.h>
#include <avtHistogramSpecification.h>
#include <FileFunctions.h>

#include "hdf5_fastquery.h"

#include "histogramCache.h"


class DBOptionsAttributes;
class avtDataRangeSelection;
class avtIdentifierSelection;


// ****************************************************************************
//  Class: avtHDF_UCFileFormat
//
//  Purpose:
//      Reads in HDF_UC files as a plugin to VisIt.
//
//  Programmer: prabhat -- generated by xml2avt
//  Creation:   Wed Dec 19 09:24:48 PDT 2007
//
//  Modifications:
//
//    Hank Childs, Thu Mar  6 09:09:32 PST 2008
//    Add method ConstructIdentifiersFromDataRangeSelection.
//
//    Hank Childs, Fri Mar  7 11:16:39 PST 2008
//    Re-inherit from MTSD.
//
// ****************************************************************************

class avtHDF_UCFileFormat : public avtMTSDFileFormat
{
  public:
                       avtHDF_UCFileFormat(const char *filename, DBOptionsAttributes *);
    virtual           ~avtHDF_UCFileFormat() {
          delete histoCache;
    };

    //
    // This is used to return unconvention data -- ranging from material
    // information to information about block connectivity.
    //
    virtual void      *GetAuxiliaryData(const char *var, 
                                        int timestep,
                                        const char *type,
                                        void *args, DestructorFunction &);

    // PRABHAT
    virtual int            GetNTimesteps(void);

    //
    // These are used to declare what the current time and cycle are for the
    // file.  These should only be defined if the file format knows what the
    // time and/or cycle is.
    //
    // PRABHAT (low priority)
    //
    // virtual bool      ReturnsValidCycle() const { return true; };
    // virtual int       GetCycle(void);
    // virtual bool      ReturnsValidTime() const { return true; };
    // virtual double    GetTime(void);
    //
    
    virtual const char    *GetType(void)   { return "HDF_UC"; };
    virtual void           FreeUpResources(void); 

    virtual bool           CanCacheVariable(const char *) { return !fbIsAllowed; };
    virtual void           RegisterDataSelections(const std::vector<avtDataSelection_p>&, std::vector<bool> *);

    virtual vtkDataSet    *GetMesh(int ts, const char *);
    virtual vtkDataArray  *GetVar(int ts, const char *);
    virtual vtkDataArray  *GetVectorVar(int ts, const char *);
    void                   AddFileInThisDirectory(const std::string &filenameWithDir);

    virtual void           ActivateTimestep(int ts);

 private:
    std::string            stringify(double x);
    void                   ConstructHistogram(avtHistogramSpecification *spec);

    avtIdentifierSelection *ConstructIdentifiersFromDataRangeSelection(
                                       std::vector<avtDataSelection *> &);
    
    void                   findAllTimesteps(const char* filename);
    void                   updateReader(int ts, bool runQuery);
    void                   runQuery();

  protected:
    HDF5_FQ                           reader;
    bool                              is2D;
    bool                              fbIsAllowed;
    bool                              readAllData;
    std::vector<avtDataSelection_p>   selList;
    std::vector<bool>                *selsApplied;
    std::string                       queryString;
    std::vector<hsize_t>              queryResults;

    int                               currentTimestep;
    int                               readerTimestep;
    std::vector<string>               fileNames;
    bool                              processAllTimesteps;
    bool                              querySpecified;
    bool                              userMadeSelection;
    
    int                get_string_from_identifiers(const vector<double>& Identifiers, string& id_string);
    
    virtual void           PopulateDatabaseMetaData(avtDatabaseMetaData *);
    
    //Chache histograms
    histogramCache*                    histoCache;

    void                               printOffsets();
    void                               printMesh(int numParticles, float*ptr, string msg);

};


#endif


