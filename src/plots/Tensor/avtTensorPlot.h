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
//                                 avtTensorPlot.h                           //
// ************************************************************************* //

#ifndef AVT_Tensor_PLOT_H
#define AVT_Tensor_PLOT_H


#include <avtLegend.h>
#include <avtPlot.h>

#include <TensorAttributes.h>

class     avtTensorFilter;
class     avtTensorGlyphMapper;
class     avtGhostZoneFilter;
class     avtVariableLegend;
class     avtLookupTable;


// ****************************************************************************
//  Class:  avtTensorPlot
//
//  Purpose:
//      A concrete type of avtPlot, this is the Tensor plot.
//
//  Programmer: childs -- generated by xml2info
//  Creation:   Tue Sep 23 20:57:03 PST 2003
//
//  Modifications:
//
//    Mark C. Miller, Wed Aug 11 23:42:18 PDT 2004
//    Added GetCellCountMultiplierForSRThreshold
//
//    Mark C. Miller, Mon Aug 23 20:24:31 PDT 2004
//    Changed GetCellCountMultiplierForSRThreshold to Set...
//
//    Hank Childs, Wed Aug 13 13:45:33 PDT 2008
//    Turn on NeedZBufferToCompositeEvenIn2D, as tensor glyphs can bleed
//    into other processor's portion of image space.
//
// ****************************************************************************

class avtTensorPlot : public avtPointDataPlot
{
  public:
                                avtTensorPlot();
    virtual                    ~avtTensorPlot();

    virtual const char         *GetName(void) { return "TensorPlot"; };
    virtual void                ReleaseData(void);

    static avtPlot             *Create();

    virtual void                SetAtts(const AttributeGroup*);
    virtual bool                NeedZBufferToCompositeEvenIn2D(void)
                                                          { return true; };


  protected:
    TensorAttributes              atts;
    bool                          colorsInitialized;

    avtTensorGlyphMapper         *tensorMapper;
    avtVariableLegend            *varLegend;
    avtLegend_p                   varLegendRefPtr;
    avtTensorFilter              *TensorFilter;
    avtGhostZoneFilter           *ghostFilter;
    avtLookupTable               *avtLUT;


    virtual avtMapper          *GetMapper(void);
    virtual avtDataObject_p     ApplyOperators(avtDataObject_p);
    virtual avtDataObject_p     ApplyRenderingTransformation(avtDataObject_p);
    virtual void                CustomizeBehavior(void);
    virtual void                CustomizeMapper(avtDataObjectInformation &);
    bool                        SetColorTable(const char *);

    virtual avtLegend_p         GetLegend(void) { return varLegendRefPtr; };
    void                        SetLegend(bool);
    void                        SetLegendRanges();

    virtual void                SetCellCountMultiplierForSRThreshold(const avtDataObject_p)
                                    { cellCountMultiplierForSRThreshold = 96.0; }; 
};


#endif


