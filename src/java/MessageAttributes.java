// ***************************************************************************
//
// Copyright (c) 2000 - 2008, Lawrence Livermore National Security, LLC
// Produced at the Lawrence Livermore National Laboratory
// LLNL-CODE-400142
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

package llnl.visit;


// ****************************************************************************
// Class: MessageAttributes
//
// Purpose:
//    This class contains attributes for sending messages.
//
// Notes:      Autogenerated by xml2java.
//
// Programmer: xml2java
// Creation:   Wed Mar 14 17:55:00 PST 2007
//
// Modifications:
//   
// ****************************************************************************

public class MessageAttributes extends AttributeSubject
{
    // Enum values
    public final static int SEVERITY_ERROR = 0;
    public final static int SEVERITY_WARNING = 1;
    public final static int SEVERITY_MESSAGE = 2;
    public final static int SEVERITY_ERRORCLEAR = 3;


    public MessageAttributes()
    {
        super(2);

        text = new String("");
        severity = SEVERITY_MESSAGE;
    }

    public MessageAttributes(MessageAttributes obj)
    {
        super(2);

        text = new String(obj.text);
        severity = obj.severity;

        SelectAll();
    }

    public boolean equals(MessageAttributes obj)
    {
        // Create the return value
        return ((text == obj.text) &&
                (severity == obj.severity));
    }

    // Property setting methods
    public void SetText(String text_)
    {
        text = text_;
        Select(0);
    }

    public void SetSeverity(int severity_)
    {
        severity = severity_;
        Select(1);
    }

    // Property getting methods
    public String GetText() { return text; }
    public int    GetSeverity() { return severity; }

    // Write and read methods.
    public void WriteAtts(CommunicationBuffer buf)
    {
        if(WriteSelect(0, buf))
            buf.WriteString(text);
        if(WriteSelect(1, buf))
            buf.WriteInt(severity);
    }

    public void ReadAtts(int n, CommunicationBuffer buf)
    {
        for(int i = 0; i < n; ++i)
        {
            int index = (int)buf.ReadByte();
            switch(index)
            {
            case 0:
                SetText(buf.ReadString());
                break;
            case 1:
                SetSeverity(buf.ReadInt());
                break;
            }
        }
    }


    // Attributes
    private String text;
    private int    severity;
}

