// ***************************************************************************
//
// Copyright (c) 2000 - 2011, Lawrence Livermore National Security, LLC
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
// Class: SiloDumpAttributes
//
// Purpose:
//    This class contains attributes for the silo dump operator.
//
// Notes:      Autogenerated by xml2java.
//
// Programmer: xml2java
// Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

public class SiloDumpAttributes extends AttributeSubject implements Plugin
{
    private static int SiloDumpAttributes_numAdditionalAtts = 2;

    public SiloDumpAttributes()
    {
        super(SiloDumpAttributes_numAdditionalAtts);

        filename = new String("dump");
        display = false;
    }

    public SiloDumpAttributes(int nMoreFields)
    {
        super(SiloDumpAttributes_numAdditionalAtts + nMoreFields);

        filename = new String("dump");
        display = false;
    }

    public SiloDumpAttributes(SiloDumpAttributes obj)
    {
        super(SiloDumpAttributes_numAdditionalAtts);

        filename = new String(obj.filename);
        display = obj.display;

        SelectAll();
    }

    public int Offset()
    {
        return super.Offset() + super.GetNumAdditionalAttributes();
    }

    public int GetNumAdditionalAttributes()
    {
        return SiloDumpAttributes_numAdditionalAtts;
    }

    public boolean equals(SiloDumpAttributes obj)
    {
        // Create the return value
        return ((filename.equals(obj.filename)) &&
                (display == obj.display));
    }

    public String GetName() { return "SiloDump"; }
    public String GetVersion() { return "1.0"; }

    // Property setting methods
    public void SetFilename(String filename_)
    {
        filename = filename_;
        Select(0);
    }

    public void SetDisplay(boolean display_)
    {
        display = display_;
        Select(1);
    }

    // Property getting methods
    public String  GetFilename() { return filename; }
    public boolean GetDisplay() { return display; }

    // Write and read methods.
    public void WriteAtts(CommunicationBuffer buf)
    {
        if(WriteSelect(0, buf))
            buf.WriteString(filename);
        if(WriteSelect(1, buf))
            buf.WriteBool(display);
    }

    public void ReadAtts(int index, CommunicationBuffer buf)
    {
        switch(index)
        {
        case 0:
            SetFilename(buf.ReadString());
            break;
        case 1:
            SetDisplay(buf.ReadBool());
            break;
        }
    }

    public String toString(String indent)
    {
        String str = new String();
        str = str + stringToString("filename", filename, indent) + "\n";
        str = str + boolToString("display", display, indent) + "\n";
        return str;
    }


    // Attributes
    private String  filename;
    private boolean display;
}

