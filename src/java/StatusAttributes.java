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
// Class: StatusAttributes
//
// Purpose:
//    This class contains the status that is displayed in the GUI's status bar.
//
// Notes:      Autogenerated by xml2java.
//
// Programmer: xml2java
// Creation:   Wed Mar 14 17:55:15 PST 2007
//
// Modifications:
//   
// ****************************************************************************

public class StatusAttributes extends AttributeSubject
{
    // Constants
    public final static int DEFAULT_DURATION = 5000;


    public StatusAttributes()
    {
        super(9);

        sender = new String("viewer");
        clearStatus = false;
        statusMessage = new String("");
        percent = 0;
        currentStage = 1;
        currentStageName = new String("stage1");
        maxStage = 1;
        messageType = 0;
        duration = 5000;
    }

    public StatusAttributes(StatusAttributes obj)
    {
        super(9);

        sender = new String(obj.sender);
        clearStatus = obj.clearStatus;
        statusMessage = new String(obj.statusMessage);
        percent = obj.percent;
        currentStage = obj.currentStage;
        currentStageName = new String(obj.currentStageName);
        maxStage = obj.maxStage;
        messageType = obj.messageType;
        duration = obj.duration;

        SelectAll();
    }

    public boolean equals(StatusAttributes obj)
    {
        // Create the return value
        return ((sender == obj.sender) &&
                (clearStatus == obj.clearStatus) &&
                (statusMessage == obj.statusMessage) &&
                (percent == obj.percent) &&
                (currentStage == obj.currentStage) &&
                (currentStageName == obj.currentStageName) &&
                (maxStage == obj.maxStage) &&
                (messageType == obj.messageType) &&
                (duration == obj.duration));
    }

    // Property setting methods
    public void SetSender(String sender_)
    {
        sender = sender_;
        Select(0);
    }

    public void SetClearStatus(boolean clearStatus_)
    {
        clearStatus = clearStatus_;
        Select(1);
    }

    public void SetStatusMessage(String statusMessage_)
    {
        statusMessage = statusMessage_;
        Select(2);
    }

    public void SetPercent(int percent_)
    {
        percent = percent_;
        Select(3);
    }

    public void SetCurrentStage(int currentStage_)
    {
        currentStage = currentStage_;
        Select(4);
    }

    public void SetCurrentStageName(String currentStageName_)
    {
        currentStageName = currentStageName_;
        Select(5);
    }

    public void SetMaxStage(int maxStage_)
    {
        maxStage = maxStage_;
        Select(6);
    }

    public void SetMessageType(int messageType_)
    {
        messageType = messageType_;
        Select(7);
    }

    public void SetDuration(int duration_)
    {
        duration = duration_;
        Select(8);
    }

    // Property getting methods
    public String  GetSender() { return sender; }
    public boolean GetClearStatus() { return clearStatus; }
    public String  GetStatusMessage() { return statusMessage; }
    public int     GetPercent() { return percent; }
    public int     GetCurrentStage() { return currentStage; }
    public String  GetCurrentStageName() { return currentStageName; }
    public int     GetMaxStage() { return maxStage; }
    public int     GetMessageType() { return messageType; }
    public int     GetDuration() { return duration; }

    // Write and read methods.
    public void WriteAtts(CommunicationBuffer buf)
    {
        if(WriteSelect(0, buf))
            buf.WriteString(sender);
        if(WriteSelect(1, buf))
            buf.WriteBool(clearStatus);
        if(WriteSelect(2, buf))
            buf.WriteString(statusMessage);
        if(WriteSelect(3, buf))
            buf.WriteInt(percent);
        if(WriteSelect(4, buf))
            buf.WriteInt(currentStage);
        if(WriteSelect(5, buf))
            buf.WriteString(currentStageName);
        if(WriteSelect(6, buf))
            buf.WriteInt(maxStage);
        if(WriteSelect(7, buf))
            buf.WriteInt(messageType);
        if(WriteSelect(8, buf))
            buf.WriteInt(duration);
    }

    public void ReadAtts(int n, CommunicationBuffer buf)
    {
        for(int i = 0; i < n; ++i)
        {
            int index = (int)buf.ReadByte();
            switch(index)
            {
            case 0:
                SetSender(buf.ReadString());
                break;
            case 1:
                SetClearStatus(buf.ReadBool());
                break;
            case 2:
                SetStatusMessage(buf.ReadString());
                break;
            case 3:
                SetPercent(buf.ReadInt());
                break;
            case 4:
                SetCurrentStage(buf.ReadInt());
                break;
            case 5:
                SetCurrentStageName(buf.ReadString());
                break;
            case 6:
                SetMaxStage(buf.ReadInt());
                break;
            case 7:
                SetMessageType(buf.ReadInt());
                break;
            case 8:
                SetDuration(buf.ReadInt());
                break;
            }
        }
    }


    // Attributes
    private String  sender;
    private boolean clearStatus;
    private String  statusMessage;
    private int     percent;
    private int     currentStage;
    private String  currentStageName;
    private int     maxStage;
    private int     messageType;
    private int     duration;
}

