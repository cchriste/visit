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
//                            avtImageCompositer.h                           //
// ************************************************************************* //

#ifndef AVT_IMAGE_COMPOSITER_H
#define AVT_IMAGE_COMPOSITER_H
#include <filters_exports.h>
#ifdef PARALLEL
#include <mpi.h>
#endif

#include <avtImage.h>
#include <avtImageToImageFilter.h>

// ****************************************************************************
//  Class: avtImageCompositer
//
//  Purpose:
//      This is a base class to support various kinds of image compositing
//      algorithms. This class is intended to provide an asbtract
//      interface for compositing images. When used within a parallel client,
//      this object can support composition of a single image rendered on
//      different processors into a complete image. Alternatively, it can also
//      handle sequential compositing of images from multiple sources on the
//      same processor, for example when two different engines are running and
//      serving up images to the viewer.
//
//  Programmer: Mark C. Miller 
//  Creation:   February 12, 2003
//
//  Modifications:
//
//    Hank Childs, Thu Feb  5 17:11:06 PST 2004
//    Moved inlined destructor definition to .C file because certain compilers
//    have problems with them.
//
//    Hank Childs, Sun Dec  4 17:07:52 PST 2005
//    Added description method.
//
// ****************************************************************************

class AVTFILTERS_API avtImageCompositer : public avtImageToImageFilter
{
   public:
                              avtImageCompositer();
      virtual                ~avtImageCompositer();

      virtual const char     *GetType(void) { return "avtImageCompositer"; };
      virtual const char     *GetDescription(void)
                                   { return "Compositing images."; };

      void                    SetOutputImageSize(const int numRows,
                                                 const int numCols);
      void                    GetOutputImageSize(int *numRows,
                                                 int *numCols) const;

      void                    AddImageInput(avtImage_p subImage,
                                 int rowOrigin, int colOrigin);

      void                    SetShouldOutputZBuffer(bool outputZBuffer);
      bool                    GetShouldOutputZBuffer();

      int                     SetRoot(const int mpiRank);
      int                     GetRoot() const;
      int                     GetRank() const;

#ifdef PARALLEL
      MPI_Comm                SetMPICommunicator(const MPI_Comm comm);
      MPI_Comm                GetMPICommunicator() const;
#endif

   protected:
      int                     outRows, outCols;
      bool                    shouldOutputZBuffer;
      std::vector<avtImage_p> inputImages;
      int                     mpiRoot;
      int                     mpiRank;
#ifdef PARALLEL
      MPI_Comm                mpiComm;
#endif

};

inline void avtImageCompositer::SetOutputImageSize(int numRows, int numCols)
{ outRows = numRows; outCols = numCols; }

inline void avtImageCompositer::SetShouldOutputZBuffer(bool outputZBuffer)
{ shouldOutputZBuffer = outputZBuffer; }

inline bool avtImageCompositer::GetShouldOutputZBuffer()
{ return shouldOutputZBuffer; }

inline int avtImageCompositer::GetRoot() const
{ return mpiRoot; }

inline int avtImageCompositer::GetRank() const
{ return mpiRank; }

#ifdef PARALLEL
inline MPI_Comm avtImageCompositer::GetMPICommunicator() const
{ return mpiComm; }
#endif

#endif
