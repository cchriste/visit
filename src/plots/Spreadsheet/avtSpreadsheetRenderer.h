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
#ifndef AVT_SPREADSHEET_RENDERER_H
#define AVT_SPREADSHEET_RENDERER_H

#include <avtCustomRenderer.h>
#include <SpreadsheetAttributes.h>

class avtSpreadsheetTraceRenderer;

// ****************************************************************************
//  Class: avtSpreadsheetRenderer
//
//  Purpose:
//      An implementation of an avtCustomRenderer for a spreadsheet plot.
//
//  Programmer: Brad Whitlock
//  Creation:   Tue Feb 6 16:01:07 PST 2007
//
//  Modifications:
//
// ****************************************************************************

class avtSpreadsheetRenderer : public avtCustomRenderer
{
public:
                            avtSpreadsheetRenderer();
    virtual                ~avtSpreadsheetRenderer();
    static avtSpreadsheetRenderer *New(void);

    virtual bool            OperatesOnScalars(void) { return true; };
    virtual void            Render(vtkDataSet *);
    virtual void            ReleaseGraphicsResources();
    virtual void            SetAlternateDisplay(void *);

    void                    SetAtts(const AttributeGroup *);
    bool                    SetColorTable(const char *);
    bool                    SetForegroundColor(const double *);
private:
    void RenderTracePlane(vtkDataSet *);

    void                        *plotDisplay;
    SpreadsheetAttributes        atts;
    double                       fgColor[3];
    avtSpreadsheetTraceRenderer *rendererImplementation;
};


typedef ref_ptr<avtSpreadsheetRenderer> avtSpreadsheetRenderer_p;


#endif
