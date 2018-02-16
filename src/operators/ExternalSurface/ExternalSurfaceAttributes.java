// ***************************************************************************
//
// Copyright (c) 2000 - 2018, Lawrence Livermore National Security, LLC
// Produced at the Lawrence Livermore National Laboratory
// LLNL-CODE-442911
// All rights reserved.
//
// This file is  part of VisIt. For  details, see https://visit.llnl.gov/.  The
// full copyright notice is contained in the file COPYRIGHT located at the root
// of the VisIt distribution or at http://www.llnl.gov/visit/copyright.html.
//
// Redistribution  and  use  in  source  and  binary  forms,  with  or  without
// modification, are permitted provided that the following conditions are met:
//
//  - Redistributions of  source code must  retain the above  copyright notice,
//    this list of conditions and the disclaimer below.
//  - Redistributions in binary form must reproduce the above copyright notice,
//    this  list of  conditions  and  the  disclaimer (as noted below)  in  the
//    documentation and/or other materials provided with the distribution.
//  - Neither the name of  the LLNS/LLNL nor the names of  its contributors may
//    be used to endorse or promote products derived from this software without
//    specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT  HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR  IMPLIED WARRANTIES, INCLUDING,  BUT NOT  LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND  FITNESS FOR A PARTICULAR  PURPOSE
// ARE  DISCLAIMED. IN  NO EVENT  SHALL LAWRENCE  LIVERMORE NATIONAL  SECURITY,
// LLC, THE  U.S.  DEPARTMENT OF  ENERGY  OR  CONTRIBUTORS BE  LIABLE  FOR  ANY
// DIRECT,  INDIRECT,   INCIDENTAL,   SPECIAL,   EXEMPLARY,  OR   CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT  LIMITED TO, PROCUREMENT OF  SUBSTITUTE GOODS OR
// SERVICES; LOSS OF  USE, DATA, OR PROFITS; OR  BUSINESS INTERRUPTION) HOWEVER
// CAUSED  AND  ON  ANY  THEORY  OF  LIABILITY,  WHETHER  IN  CONTRACT,  STRICT
// LIABILITY, OR TORT  (INCLUDING NEGLIGENCE OR OTHERWISE)  ARISING IN ANY  WAY
// OUT OF THE  USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ***************************************************************************

package llnl.visit.operators;

import llnl.visit.AttributeSubject;
import llnl.visit.CommunicationBuffer;
import llnl.visit.Plugin;

// ****************************************************************************
// Class: ExternalSurfaceAttributes
//
// Purpose:
//    This class contains attributes for the external surface operator.
//
// Notes:      Autogenerated by xml2java.
//
// Programmer: xml2java
// Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

public class ExternalSurfaceAttributes extends AttributeSubject implements Plugin
{
    private static int ExternalSurfaceAttributes_numAdditionalAtts = 2;

    public ExternalSurfaceAttributes()
    {
        super(ExternalSurfaceAttributes_numAdditionalAtts);

        removeGhosts = false;
        edgesIn2D = true;
    }

    public ExternalSurfaceAttributes(int nMoreFields)
    {
        super(ExternalSurfaceAttributes_numAdditionalAtts + nMoreFields);

        removeGhosts = false;
        edgesIn2D = true;
    }

    public ExternalSurfaceAttributes(ExternalSurfaceAttributes obj)
    {
        super(obj);

        removeGhosts = obj.removeGhosts;
        edgesIn2D = obj.edgesIn2D;

        SelectAll();
    }

    public int Offset()
    {
        return super.Offset() + super.GetNumAdditionalAttributes();
    }

    public int GetNumAdditionalAttributes()
    {
        return ExternalSurfaceAttributes_numAdditionalAtts;
    }

    public boolean equals(ExternalSurfaceAttributes obj)
    {
        // Create the return value
        return ((removeGhosts == obj.removeGhosts) &&
                (edgesIn2D == obj.edgesIn2D));
    }

    public String GetName() { return "ExternalSurface"; }
    public String GetVersion() { return "1.0"; }

    // Property setting methods
    public void SetRemoveGhosts(boolean removeGhosts_)
    {
        removeGhosts = removeGhosts_;
        Select(0);
    }

    public void SetEdgesIn2D(boolean edgesIn2D_)
    {
        edgesIn2D = edgesIn2D_;
        Select(1);
    }

    // Property getting methods
    public boolean GetRemoveGhosts() { return removeGhosts; }
    public boolean GetEdgesIn2D() { return edgesIn2D; }

    // Write and read methods.
    public void WriteAtts(CommunicationBuffer buf)
    {
        if(WriteSelect(0, buf))
            buf.WriteBool(removeGhosts);
        if(WriteSelect(1, buf))
            buf.WriteBool(edgesIn2D);
    }

    public void ReadAtts(int index, CommunicationBuffer buf)
    {
        switch(index)
        {
        case 0:
            SetRemoveGhosts(buf.ReadBool());
            break;
        case 1:
            SetEdgesIn2D(buf.ReadBool());
            break;
        }
    }

    public String toString(String indent)
    {
        String str = new String();
        str = str + boolToString("removeGhosts", removeGhosts, indent) + "\n";
        str = str + boolToString("edgesIn2D", edgesIn2D, indent) + "\n";
        return str;
    }


    // Attributes
    private boolean removeGhosts;
    private boolean edgesIn2D;
}

