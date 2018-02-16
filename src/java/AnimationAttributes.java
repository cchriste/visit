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

package llnl.visit;


// ****************************************************************************
// Class: AnimationAttributes
//
// Purpose:
//    This class contains the animation attributes.
//
// Notes:      Autogenerated by xml2java.
//
// Programmer: xml2java
// Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

public class AnimationAttributes extends AttributeSubject
{
    private static int AnimationAttributes_numAdditionalAtts = 5;

    // Enum values
    public final static int ANIMATIONMODE_REVERSEPLAYMODE = 0;
    public final static int ANIMATIONMODE_STOPMODE = 1;
    public final static int ANIMATIONMODE_PLAYMODE = 2;

    public final static int PLAYBACKMODE_LOOPING = 0;
    public final static int PLAYBACKMODE_PLAYONCE = 1;
    public final static int PLAYBACKMODE_SWING = 2;


    public AnimationAttributes()
    {
        super(AnimationAttributes_numAdditionalAtts);

        animationMode = ANIMATIONMODE_STOPMODE;
        pipelineCachingMode = false;
        frameIncrement = 1;
        timeout = 1;
        playbackMode = PLAYBACKMODE_LOOPING;
    }

    public AnimationAttributes(int nMoreFields)
    {
        super(AnimationAttributes_numAdditionalAtts + nMoreFields);

        animationMode = ANIMATIONMODE_STOPMODE;
        pipelineCachingMode = false;
        frameIncrement = 1;
        timeout = 1;
        playbackMode = PLAYBACKMODE_LOOPING;
    }

    public AnimationAttributes(AnimationAttributes obj)
    {
        super(obj);

        animationMode = obj.animationMode;
        pipelineCachingMode = obj.pipelineCachingMode;
        frameIncrement = obj.frameIncrement;
        timeout = obj.timeout;
        playbackMode = obj.playbackMode;

        SelectAll();
    }

    public int Offset()
    {
        return super.Offset() + super.GetNumAdditionalAttributes();
    }

    public int GetNumAdditionalAttributes()
    {
        return AnimationAttributes_numAdditionalAtts;
    }

    public boolean equals(AnimationAttributes obj)
    {
        // Create the return value
        return ((animationMode == obj.animationMode) &&
                (pipelineCachingMode == obj.pipelineCachingMode) &&
                (frameIncrement == obj.frameIncrement) &&
                (timeout == obj.timeout) &&
                (playbackMode == obj.playbackMode));
    }

    // Property setting methods
    public void SetAnimationMode(int animationMode_)
    {
        animationMode = animationMode_;
        Select(0);
    }

    public void SetPipelineCachingMode(boolean pipelineCachingMode_)
    {
        pipelineCachingMode = pipelineCachingMode_;
        Select(1);
    }

    public void SetFrameIncrement(int frameIncrement_)
    {
        frameIncrement = frameIncrement_;
        Select(2);
    }

    public void SetTimeout(int timeout_)
    {
        timeout = timeout_;
        Select(3);
    }

    public void SetPlaybackMode(int playbackMode_)
    {
        playbackMode = playbackMode_;
        Select(4);
    }

    // Property getting methods
    public int     GetAnimationMode() { return animationMode; }
    public boolean GetPipelineCachingMode() { return pipelineCachingMode; }
    public int     GetFrameIncrement() { return frameIncrement; }
    public int     GetTimeout() { return timeout; }
    public int     GetPlaybackMode() { return playbackMode; }

    // Write and read methods.
    public void WriteAtts(CommunicationBuffer buf)
    {
        if(WriteSelect(0, buf))
            buf.WriteInt(animationMode);
        if(WriteSelect(1, buf))
            buf.WriteBool(pipelineCachingMode);
        if(WriteSelect(2, buf))
            buf.WriteInt(frameIncrement);
        if(WriteSelect(3, buf))
            buf.WriteInt(timeout);
        if(WriteSelect(4, buf))
            buf.WriteInt(playbackMode);
    }

    public void ReadAtts(int index, CommunicationBuffer buf)
    {
        switch(index)
        {
        case 0:
            SetAnimationMode(buf.ReadInt());
            break;
        case 1:
            SetPipelineCachingMode(buf.ReadBool());
            break;
        case 2:
            SetFrameIncrement(buf.ReadInt());
            break;
        case 3:
            SetTimeout(buf.ReadInt());
            break;
        case 4:
            SetPlaybackMode(buf.ReadInt());
            break;
        }
    }

    public String toString(String indent)
    {
        String str = new String();
        str = str + indent + "animationMode = ";
        if(animationMode == ANIMATIONMODE_REVERSEPLAYMODE)
            str = str + "ANIMATIONMODE_REVERSEPLAYMODE";
        if(animationMode == ANIMATIONMODE_STOPMODE)
            str = str + "ANIMATIONMODE_STOPMODE";
        if(animationMode == ANIMATIONMODE_PLAYMODE)
            str = str + "ANIMATIONMODE_PLAYMODE";
        str = str + "\n";
        str = str + boolToString("pipelineCachingMode", pipelineCachingMode, indent) + "\n";
        str = str + intToString("frameIncrement", frameIncrement, indent) + "\n";
        str = str + intToString("timeout", timeout, indent) + "\n";
        str = str + indent + "playbackMode = ";
        if(playbackMode == PLAYBACKMODE_LOOPING)
            str = str + "PLAYBACKMODE_LOOPING";
        if(playbackMode == PLAYBACKMODE_PLAYONCE)
            str = str + "PLAYBACKMODE_PLAYONCE";
        if(playbackMode == PLAYBACKMODE_SWING)
            str = str + "PLAYBACKMODE_SWING";
        str = str + "\n";
        return str;
    }


    // Attributes
    private int     animationMode;
    private boolean pipelineCachingMode;
    private int     frameIncrement;
    private int     timeout;
    private int     playbackMode;
}

