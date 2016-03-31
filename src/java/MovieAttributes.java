// ***************************************************************************
//
// Copyright (c) 2000 - 2016, Lawrence Livermore National Security, LLC
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

import java.util.Vector;
import java.lang.Byte;
import java.lang.Integer;
import java.lang.Double;

// ****************************************************************************
// Class: MovieAttributes
//
// Purpose:
//    This class contains the attributes used for saving movies.
//
// Notes:      Autogenerated by xml2java.
//
// Programmer: xml2java
// Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

public class MovieAttributes extends AttributeSubject
{
    private static int MovieAttributes_numAdditionalAtts = 19;

    // Enum values
    public final static int MOVIETYPEENUM_SIMPLE = 0;
    public final static int MOVIETYPEENUM_USINGTEMPLATE = 1;

    public final static int GENERATIONMETHODENUM_NOWCURRENTINSTANCE = 0;
    public final static int GENERATIONMETHODENUM_NOWNEWINSTANCE = 1;
    public final static int GENERATIONMETHODENUM_LATER = 2;


    public MovieAttributes()
    {
        super(MovieAttributes_numAdditionalAtts);

        generationMethod = GENERATIONMETHODENUM_NOWCURRENTINSTANCE;
        movieType = MOVIETYPEENUM_SIMPLE;
        outputDirectory = new String(".");
        outputName = new String("movie");
        fileFormats = new Vector();
        useCurrentSize = new Vector();
        widths = new Vector();
        heights = new Vector();
        scales = new Vector();
        stereoFlags = new Vector();
        templateFile = new String("");
        sendEmailNotification = false;
        useScreenCapture = false;
        emailAddress = new String("");
        fps = 10;
        startIndex = 0;
        endIndex = 1000000000;
        stride = 1;
        initialFrameValue = 0;
    }

    public MovieAttributes(int nMoreFields)
    {
        super(MovieAttributes_numAdditionalAtts + nMoreFields);

        generationMethod = GENERATIONMETHODENUM_NOWCURRENTINSTANCE;
        movieType = MOVIETYPEENUM_SIMPLE;
        outputDirectory = new String(".");
        outputName = new String("movie");
        fileFormats = new Vector();
        useCurrentSize = new Vector();
        widths = new Vector();
        heights = new Vector();
        scales = new Vector();
        stereoFlags = new Vector();
        templateFile = new String("");
        sendEmailNotification = false;
        useScreenCapture = false;
        emailAddress = new String("");
        fps = 10;
        startIndex = 0;
        endIndex = 1000000000;
        stride = 1;
        initialFrameValue = 0;
    }

    public MovieAttributes(MovieAttributes obj)
    {
        super(MovieAttributes_numAdditionalAtts);

        int i;

        generationMethod = obj.generationMethod;
        movieType = obj.movieType;
        outputDirectory = new String(obj.outputDirectory);
        outputName = new String(obj.outputName);
        fileFormats = new Vector(obj.fileFormats.size());
        for(i = 0; i < obj.fileFormats.size(); ++i)
            fileFormats.addElement(new String((String)obj.fileFormats.elementAt(i)));

        useCurrentSize = new Vector(obj.useCurrentSize.size());
        for(i = 0; i < obj.useCurrentSize.size(); ++i)
        {
            Byte bv = (Byte)obj.useCurrentSize.elementAt(i);
            useCurrentSize.addElement(new Byte(bv.byteValue()));
        }

        widths = new Vector();
        for(i = 0; i < obj.widths.size(); ++i)
        {
            Integer iv = (Integer)obj.widths.elementAt(i);
            widths.addElement(new Integer(iv.intValue()));
        }
        heights = new Vector();
        for(i = 0; i < obj.heights.size(); ++i)
        {
            Integer iv = (Integer)obj.heights.elementAt(i);
            heights.addElement(new Integer(iv.intValue()));
        }
        scales = new Vector(obj.scales.size());
        for(i = 0; i < obj.scales.size(); ++i)
        {
            Double dv = (Double)obj.scales.elementAt(i);
            scales.addElement(new Double(dv.doubleValue()));
        }

        stereoFlags = new Vector();
        for(i = 0; i < obj.stereoFlags.size(); ++i)
        {
            Integer iv = (Integer)obj.stereoFlags.elementAt(i);
            stereoFlags.addElement(new Integer(iv.intValue()));
        }
        templateFile = new String(obj.templateFile);
        sendEmailNotification = obj.sendEmailNotification;
        useScreenCapture = obj.useScreenCapture;
        emailAddress = new String(obj.emailAddress);
        fps = obj.fps;
        startIndex = obj.startIndex;
        endIndex = obj.endIndex;
        stride = obj.stride;
        initialFrameValue = obj.initialFrameValue;

        SelectAll();
    }

    public int Offset()
    {
        return super.Offset() + super.GetNumAdditionalAttributes();
    }

    public int GetNumAdditionalAttributes()
    {
        return MovieAttributes_numAdditionalAtts;
    }

    public boolean equals(MovieAttributes obj)
    {
        int i;

        // Compare the elements in the fileFormats vector.
        boolean fileFormats_equal = (obj.fileFormats.size() == fileFormats.size());
        for(i = 0; (i < fileFormats.size()) && fileFormats_equal; ++i)
        {
            // Make references to String from Object.
            String fileFormats1 = (String)fileFormats.elementAt(i);
            String fileFormats2 = (String)obj.fileFormats.elementAt(i);
            fileFormats_equal = fileFormats1.equals(fileFormats2);
        }
        // Compare the elements in the useCurrentSize vector.
        boolean useCurrentSize_equal = (obj.useCurrentSize.size() == useCurrentSize.size());
        for(i = 0; (i < useCurrentSize.size()) && useCurrentSize_equal; ++i)
        {
            // Make references to Byte from Object.
            Byte useCurrentSize1 = (Byte)useCurrentSize.elementAt(i);
            Byte useCurrentSize2 = (Byte)obj.useCurrentSize.elementAt(i);
            useCurrentSize_equal = useCurrentSize1.equals(useCurrentSize2);
        }
        // Compare the elements in the widths vector.
        boolean widths_equal = (obj.widths.size() == widths.size());
        for(i = 0; (i < widths.size()) && widths_equal; ++i)
        {
            // Make references to Integer from Object.
            Integer widths1 = (Integer)widths.elementAt(i);
            Integer widths2 = (Integer)obj.widths.elementAt(i);
            widths_equal = widths1.equals(widths2);
        }
        // Compare the elements in the heights vector.
        boolean heights_equal = (obj.heights.size() == heights.size());
        for(i = 0; (i < heights.size()) && heights_equal; ++i)
        {
            // Make references to Integer from Object.
            Integer heights1 = (Integer)heights.elementAt(i);
            Integer heights2 = (Integer)obj.heights.elementAt(i);
            heights_equal = heights1.equals(heights2);
        }
        // Compare the elements in the scales vector.
        boolean scales_equal = (obj.scales.size() == scales.size());
        for(i = 0; (i < scales.size()) && scales_equal; ++i)
        {
            // Make references to Double from Object.
            Double scales1 = (Double)scales.elementAt(i);
            Double scales2 = (Double)obj.scales.elementAt(i);
            scales_equal = scales1.equals(scales2);
        }
        // Compare the elements in the stereoFlags vector.
        boolean stereoFlags_equal = (obj.stereoFlags.size() == stereoFlags.size());
        for(i = 0; (i < stereoFlags.size()) && stereoFlags_equal; ++i)
        {
            // Make references to Integer from Object.
            Integer stereoFlags1 = (Integer)stereoFlags.elementAt(i);
            Integer stereoFlags2 = (Integer)obj.stereoFlags.elementAt(i);
            stereoFlags_equal = stereoFlags1.equals(stereoFlags2);
        }
        // Create the return value
        return ((generationMethod == obj.generationMethod) &&
                (movieType == obj.movieType) &&
                (outputDirectory.equals(obj.outputDirectory)) &&
                (outputName.equals(obj.outputName)) &&
                fileFormats_equal &&
                useCurrentSize_equal &&
                widths_equal &&
                heights_equal &&
                scales_equal &&
                stereoFlags_equal &&
                (templateFile.equals(obj.templateFile)) &&
                (sendEmailNotification == obj.sendEmailNotification) &&
                (useScreenCapture == obj.useScreenCapture) &&
                (emailAddress.equals(obj.emailAddress)) &&
                (fps == obj.fps) &&
                (startIndex == obj.startIndex) &&
                (endIndex == obj.endIndex) &&
                (stride == obj.stride) &&
                (initialFrameValue == obj.initialFrameValue));
    }

    // Property setting methods
    public void SetGenerationMethod(int generationMethod_)
    {
        generationMethod = generationMethod_;
        Select(0);
    }

    public void SetMovieType(int movieType_)
    {
        movieType = movieType_;
        Select(1);
    }

    public void SetOutputDirectory(String outputDirectory_)
    {
        outputDirectory = outputDirectory_;
        Select(2);
    }

    public void SetOutputName(String outputName_)
    {
        outputName = outputName_;
        Select(3);
    }

    public void SetFileFormats(Vector fileFormats_)
    {
        fileFormats = fileFormats_;
        Select(4);
    }

    public void SetUseCurrentSize(Vector useCurrentSize_)
    {
        useCurrentSize = useCurrentSize_;
        Select(5);
    }

    public void SetWidths(Vector widths_)
    {
        widths = widths_;
        Select(6);
    }

    public void SetHeights(Vector heights_)
    {
        heights = heights_;
        Select(7);
    }

    public void SetScales(Vector scales_)
    {
        scales = scales_;
        Select(8);
    }

    public void SetStereoFlags(Vector stereoFlags_)
    {
        stereoFlags = stereoFlags_;
        Select(9);
    }

    public void SetTemplateFile(String templateFile_)
    {
        templateFile = templateFile_;
        Select(10);
    }

    public void SetSendEmailNotification(boolean sendEmailNotification_)
    {
        sendEmailNotification = sendEmailNotification_;
        Select(11);
    }

    public void SetUseScreenCapture(boolean useScreenCapture_)
    {
        useScreenCapture = useScreenCapture_;
        Select(12);
    }

    public void SetEmailAddress(String emailAddress_)
    {
        emailAddress = emailAddress_;
        Select(13);
    }

    public void SetFps(int fps_)
    {
        fps = fps_;
        Select(14);
    }

    public void SetStartIndex(int startIndex_)
    {
        startIndex = startIndex_;
        Select(15);
    }

    public void SetEndIndex(int endIndex_)
    {
        endIndex = endIndex_;
        Select(16);
    }

    public void SetStride(int stride_)
    {
        stride = stride_;
        Select(17);
    }

    public void SetInitialFrameValue(int initialFrameValue_)
    {
        initialFrameValue = initialFrameValue_;
        Select(18);
    }

    // Property getting methods
    public int     GetGenerationMethod() { return generationMethod; }
    public int     GetMovieType() { return movieType; }
    public String  GetOutputDirectory() { return outputDirectory; }
    public String  GetOutputName() { return outputName; }
    public Vector  GetFileFormats() { return fileFormats; }
    public Vector  GetUseCurrentSize() { return useCurrentSize; }
    public Vector  GetWidths() { return widths; }
    public Vector  GetHeights() { return heights; }
    public Vector  GetScales() { return scales; }
    public Vector  GetStereoFlags() { return stereoFlags; }
    public String  GetTemplateFile() { return templateFile; }
    public boolean GetSendEmailNotification() { return sendEmailNotification; }
    public boolean GetUseScreenCapture() { return useScreenCapture; }
    public String  GetEmailAddress() { return emailAddress; }
    public int     GetFps() { return fps; }
    public int     GetStartIndex() { return startIndex; }
    public int     GetEndIndex() { return endIndex; }
    public int     GetStride() { return stride; }
    public int     GetInitialFrameValue() { return initialFrameValue; }

    // Write and read methods.
    public void WriteAtts(CommunicationBuffer buf)
    {
        if(WriteSelect(0, buf))
            buf.WriteInt(generationMethod);
        if(WriteSelect(1, buf))
            buf.WriteInt(movieType);
        if(WriteSelect(2, buf))
            buf.WriteString(outputDirectory);
        if(WriteSelect(3, buf))
            buf.WriteString(outputName);
        if(WriteSelect(4, buf))
            buf.WriteStringVector(fileFormats);
        if(WriteSelect(5, buf))
            buf.WriteByteVector(useCurrentSize);
        if(WriteSelect(6, buf))
            buf.WriteIntVector(widths);
        if(WriteSelect(7, buf))
            buf.WriteIntVector(heights);
        if(WriteSelect(8, buf))
            buf.WriteDoubleVector(scales);
        if(WriteSelect(9, buf))
            buf.WriteIntVector(stereoFlags);
        if(WriteSelect(10, buf))
            buf.WriteString(templateFile);
        if(WriteSelect(11, buf))
            buf.WriteBool(sendEmailNotification);
        if(WriteSelect(12, buf))
            buf.WriteBool(useScreenCapture);
        if(WriteSelect(13, buf))
            buf.WriteString(emailAddress);
        if(WriteSelect(14, buf))
            buf.WriteInt(fps);
        if(WriteSelect(15, buf))
            buf.WriteInt(startIndex);
        if(WriteSelect(16, buf))
            buf.WriteInt(endIndex);
        if(WriteSelect(17, buf))
            buf.WriteInt(stride);
        if(WriteSelect(18, buf))
            buf.WriteInt(initialFrameValue);
    }

    public void ReadAtts(int index, CommunicationBuffer buf)
    {
        switch(index)
        {
        case 0:
            SetGenerationMethod(buf.ReadInt());
            break;
        case 1:
            SetMovieType(buf.ReadInt());
            break;
        case 2:
            SetOutputDirectory(buf.ReadString());
            break;
        case 3:
            SetOutputName(buf.ReadString());
            break;
        case 4:
            SetFileFormats(buf.ReadStringVector());
            break;
        case 5:
            SetUseCurrentSize(buf.ReadByteVector());
            break;
        case 6:
            SetWidths(buf.ReadIntVector());
            break;
        case 7:
            SetHeights(buf.ReadIntVector());
            break;
        case 8:
            SetScales(buf.ReadDoubleVector());
            break;
        case 9:
            SetStereoFlags(buf.ReadIntVector());
            break;
        case 10:
            SetTemplateFile(buf.ReadString());
            break;
        case 11:
            SetSendEmailNotification(buf.ReadBool());
            break;
        case 12:
            SetUseScreenCapture(buf.ReadBool());
            break;
        case 13:
            SetEmailAddress(buf.ReadString());
            break;
        case 14:
            SetFps(buf.ReadInt());
            break;
        case 15:
            SetStartIndex(buf.ReadInt());
            break;
        case 16:
            SetEndIndex(buf.ReadInt());
            break;
        case 17:
            SetStride(buf.ReadInt());
            break;
        case 18:
            SetInitialFrameValue(buf.ReadInt());
            break;
        }
    }

    public String toString(String indent)
    {
        String str = new String();
        str = str + indent + "generationMethod = ";
        if(generationMethod == GENERATIONMETHODENUM_NOWCURRENTINSTANCE)
            str = str + "GENERATIONMETHODENUM_NOWCURRENTINSTANCE";
        if(generationMethod == GENERATIONMETHODENUM_NOWNEWINSTANCE)
            str = str + "GENERATIONMETHODENUM_NOWNEWINSTANCE";
        if(generationMethod == GENERATIONMETHODENUM_LATER)
            str = str + "GENERATIONMETHODENUM_LATER";
        str = str + "\n";
        str = str + indent + "movieType = ";
        if(movieType == MOVIETYPEENUM_SIMPLE)
            str = str + "MOVIETYPEENUM_SIMPLE";
        if(movieType == MOVIETYPEENUM_USINGTEMPLATE)
            str = str + "MOVIETYPEENUM_USINGTEMPLATE";
        str = str + "\n";
        str = str + stringToString("outputDirectory", outputDirectory, indent) + "\n";
        str = str + stringToString("outputName", outputName, indent) + "\n";
        str = str + stringVectorToString("fileFormats", fileFormats, indent) + "\n";
        str = str + ucharVectorToString("useCurrentSize", useCurrentSize, indent) + "\n";
        str = str + intVectorToString("widths", widths, indent) + "\n";
        str = str + intVectorToString("heights", heights, indent) + "\n";
        str = str + doubleVectorToString("scales", scales, indent) + "\n";
        str = str + intVectorToString("stereoFlags", stereoFlags, indent) + "\n";
        str = str + stringToString("templateFile", templateFile, indent) + "\n";
        str = str + boolToString("sendEmailNotification", sendEmailNotification, indent) + "\n";
        str = str + boolToString("useScreenCapture", useScreenCapture, indent) + "\n";
        str = str + stringToString("emailAddress", emailAddress, indent) + "\n";
        str = str + intToString("fps", fps, indent) + "\n";
        str = str + intToString("startIndex", startIndex, indent) + "\n";
        str = str + intToString("endIndex", endIndex, indent) + "\n";
        str = str + intToString("stride", stride, indent) + "\n";
        str = str + intToString("initialFrameValue", initialFrameValue, indent) + "\n";
        return str;
    }


    // Attributes
    private int     generationMethod;
    private int     movieType;
    private String  outputDirectory;
    private String  outputName;
    private Vector  fileFormats; // vector of String objects
    private Vector  useCurrentSize; // vector of Byte objects
    private Vector  widths; // vector of Integer objects
    private Vector  heights; // vector of Integer objects
    private Vector  scales; // vector of Double objects
    private Vector  stereoFlags; // vector of Integer objects
    private String  templateFile;
    private boolean sendEmailNotification;
    private boolean useScreenCapture;
    private String  emailAddress;
    private int     fps;
    private int     startIndex;
    private int     endIndex;
    private int     stride;
    private int     initialFrameValue;
}

