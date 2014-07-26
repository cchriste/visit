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

#ifndef AVTFILTERS_EXPORTS_H
#define AVTFILTERS_EXPORTS_H

#if defined(_WIN32)
#if defined(AVTFILTERS_EXPORTS) || defined(avtfilters_ser_EXPORTS) || defined(avtfilters_par_EXPORTS)
#define AVTFILTERS_API __declspec(dllexport)
#else
#define AVTFILTERS_API __declspec(dllimport)
#endif
#if defined(_MSC_VER)
// Turn off warning about lack of DLL interface
#pragma warning(disable:4251)
// Turn off warning non-dll class is base for dll-interface class.
#pragma warning(disable:4275)
// Turn off warning about identifier truncation
#pragma warning(disable:4786)
#endif
#else
# if __GNUC__ >= 4 && (defined(AVTFILTERS_EXPORTS) || defined(avtfilters_ser_EXPORTS) || defined(avtfilters_par_EXPORTS))
#   define AVTFILTERS_API __attribute__ ((visibility("default")))
# else
#   define AVTFILTERS_API /* hidden by default */
# endif
#endif

#include <stdio.h>
#include <string>
#include <iostream>

// ****************************************************************************
//  Struct:  imgMetaData
//
//  Purpose:
//    Holds information about patches but not the image 
//
//  Programmer:  
//  Creation:   
//
// ****************************************************************************
struct imgMetaData{
  int procId;       // processor that produced the patch
  int patchNumber;  // id of the patch on that processor - with procId, acts as a key

  int destProcId;   // destination proc where this patch gets composited

  int inUse;   // whether the patch is composed locally or not
  int dims[2];      // height, width
  int screen_ll[2]; // position in the final image
  int screen_ur[2];
  float avg_z;      // depth of the patch
};


// ****************************************************************************
//  Struct:  imgData
//
//  Purpose:
//    Holds the image data generated
//
//  Programmer:  
//  Creation:    
//
// ****************************************************************************
struct imgData{
  int procId;         // processor that produced the patch
  int patchNumber;    // id of the patch on that processor  - with procId, acts as a key

  float *imagePatch;  // the image data - RGBA

  bool operator==(const imgData &a){
    return (patchNumber == a.patchNumber);// && (procId == a.procId);
 }   
};


// ****************************************************************************
//  Struct:  iotaMeta
//
//  Purpose:
//    Holds the image data generated
//
//  Programmer:  
//  Creation:    
//
// ****************************************************************************
struct iotaMeta{
  int procId;  
  int patchNumber; 
  
  int dims[2];      // height, width
  int screen_ll[2]; // position in the final image
  float avg_z;   

  bool operator==(const iotaMeta &a){
    return (patchNumber == a.patchNumber) && (procId == a.procId);
 }    
};

#endif
