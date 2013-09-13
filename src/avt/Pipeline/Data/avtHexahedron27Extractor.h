/*****************************************************************************
*
* Copyright (c) 2000 - 2013, Lawrence Livermore National Security, LLC
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
//                          avtHexahedron27Extractor.h                       //
// ************************************************************************* //

#ifndef AVT_HEXAHEDRON27_EXTRACTOR_H
#define AVT_HEXAHEDRON27_EXTRACTOR_H

#include <pipeline_exports.h>

#include <avtHexahedronExtractor.h>

// ****************************************************************************
//  Class: avtHexahedron27Extractor
//
//  Purpose:
//      Extracts sample points from a higher order hexahedron.  It assumes 
//      that the hexahedrons it has been given are in camera space and 
//      does not try to populate points that are not in the cube 
//      [-1, 1], [-1, 1], [-1, 1].
//
//  Programmer: Timo Bremer
//  Creation:   September 13, 2007
//
// ****************************************************************************

class PIPELINE_API avtHexahedron27Extractor : public avtHexahedronExtractor
{
  public:
    //
    // Typedef for the different method of treating the domain. The
    // domain can be 
    //
    // avtHex27Constant:  Approximate the domain using a single linear hex 
    //                    with vertices 0-7 as corners
    // avtHex27Lienar:    Approximate the domain using eight linear hexes 
    //                    where the face vertices are the average of the 
    //                    appropriate edges and the center the average of the 
    //                    size face centers
    // avtHex27Quadratic: Compute the correct quadratic domain for sampling
    //
    enum DomainAppStyle 
    {
        avtHex27Constant = 0,
        avtHex27Linear,
        avtHex27Quadratic,
    };

  
    static const DomainAppStyle DEFAULT_DOMAIN_APPROXIMATION;
  


                     avtHexahedron27Extractor(int, int, int, avtVolume *,
                                             avtCellList *);
    virtual         ~avtHexahedron27Extractor();

    void             Extract(const avtHexahedron27 &);
    void             ConstantHexExtract(const avtHexahedron27 &);
    void             LinearHexExtract(const avtHexahedron27 &);
    void             QuadraticHexExtract(const avtHexahedron27 &);


    virtual void      StoreRay(int, int, int, int,
                               const double (*)[AVT_VARIABLE_LIMIT]);
    double             TriLinearWeight(int, double, double, double);
    double             QuadraticWeight(int, double, double, double);
    double             ClosestParametricNeighbor(double, double, double);

  protected:
    static const int              sControlPointParam[27][3];
    static const unsigned char    sSubHexCorners[8][8];
    const avtHexahedron27        *currentHex;
    DomainAppStyle                domainApproximation;
};


#endif


