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

package llnl.visit;


// ****************************************************************************
// Class: SaveWindowAttributes
//
// Purpose:
//    This class contains the attributes used for saving windows.
//
// Notes:      Autogenerated by xml2java.
//
// Programmer: xml2java
// Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

public class SaveWindowAttributes extends AttributeSubject
{
    private static int SaveWindowAttributes_numAdditionalAtts = 19;

    // Enum values
    public final static int FILEFORMAT_BMP = 0;
    public final static int FILEFORMAT_CURVE = 1;
    public final static int FILEFORMAT_JPEG = 2;
    public final static int FILEFORMAT_OBJ = 3;
    public final static int FILEFORMAT_PNG = 4;
    public final static int FILEFORMAT_POSTSCRIPT = 5;
    public final static int FILEFORMAT_POVRAY = 6;
    public final static int FILEFORMAT_PPM = 7;
    public final static int FILEFORMAT_RGB = 8;
    public final static int FILEFORMAT_STL = 9;
    public final static int FILEFORMAT_TIFF = 10;
    public final static int FILEFORMAT_ULTRA = 11;
    public final static int FILEFORMAT_VTK = 12;
    public final static int FILEFORMAT_PLY = 13;

    public final static int COMPRESSIONTYPE_NONE = 0;
    public final static int COMPRESSIONTYPE_PACKBITS = 1;
    public final static int COMPRESSIONTYPE_JPEG = 2;
    public final static int COMPRESSIONTYPE_DEFLATE = 3;

    public final static int RESCONSTRAINT_NOCONSTRAINT = 0;
    public final static int RESCONSTRAINT_EQUALWIDTHHEIGHT = 1;
    public final static int RESCONSTRAINT_SCREENPROPORTIONS = 2;


    public SaveWindowAttributes()
    {
        super(SaveWindowAttributes_numAdditionalAtts);

        outputToCurrentDirectory = true;
        outputDirectory = new String(".");
        fileName = new String("visit");
        family = true;
        format = FILEFORMAT_PNG;
        width = 1024;
        height = 1024;
        screenCapture = false;
        saveTiled = false;
        quality = 80;
        progressive = false;
        binary = false;
        lastRealFilename = new String("");
        stereo = false;
        compression = COMPRESSIONTYPE_PACKBITS;
        forceMerge = false;
        resConstraint = RESCONSTRAINT_SCREENPROPORTIONS;
        advancedMultiWindowSave = false;
        subWindowAtts = new SaveSubWindowsAttributes();
    }

    public SaveWindowAttributes(int nMoreFields)
    {
        super(SaveWindowAttributes_numAdditionalAtts + nMoreFields);

        outputToCurrentDirectory = true;
        outputDirectory = new String(".");
        fileName = new String("visit");
        family = true;
        format = FILEFORMAT_PNG;
        width = 1024;
        height = 1024;
        screenCapture = false;
        saveTiled = false;
        quality = 80;
        progressive = false;
        binary = false;
        lastRealFilename = new String("");
        stereo = false;
        compression = COMPRESSIONTYPE_PACKBITS;
        forceMerge = false;
        resConstraint = RESCONSTRAINT_SCREENPROPORTIONS;
        advancedMultiWindowSave = false;
        subWindowAtts = new SaveSubWindowsAttributes();
    }

    public SaveWindowAttributes(SaveWindowAttributes obj)
    {
        super(SaveWindowAttributes_numAdditionalAtts);

        outputToCurrentDirectory = obj.outputToCurrentDirectory;
        outputDirectory = new String(obj.outputDirectory);
        fileName = new String(obj.fileName);
        family = obj.family;
        format = obj.format;
        width = obj.width;
        height = obj.height;
        screenCapture = obj.screenCapture;
        saveTiled = obj.saveTiled;
        quality = obj.quality;
        progressive = obj.progressive;
        binary = obj.binary;
        lastRealFilename = new String(obj.lastRealFilename);
        stereo = obj.stereo;
        compression = obj.compression;
        forceMerge = obj.forceMerge;
        resConstraint = obj.resConstraint;
        advancedMultiWindowSave = obj.advancedMultiWindowSave;
        subWindowAtts = new SaveSubWindowsAttributes(obj.subWindowAtts);

        SelectAll();
    }

    public int Offset()
    {
        return super.Offset() + super.GetNumAdditionalAttributes();
    }

    public int GetNumAdditionalAttributes()
    {
        return SaveWindowAttributes_numAdditionalAtts;
    }

    public boolean equals(SaveWindowAttributes obj)
    {
        // Create the return value
        return ((outputToCurrentDirectory == obj.outputToCurrentDirectory) &&
                (outputDirectory.equals(obj.outputDirectory)) &&
                (fileName.equals(obj.fileName)) &&
                (family == obj.family) &&
                (format == obj.format) &&
                (width == obj.width) &&
                (height == obj.height) &&
                (screenCapture == obj.screenCapture) &&
                (saveTiled == obj.saveTiled) &&
                (quality == obj.quality) &&
                (progressive == obj.progressive) &&
                (binary == obj.binary) &&
                (lastRealFilename.equals(obj.lastRealFilename)) &&
                (stereo == obj.stereo) &&
                (compression == obj.compression) &&
                (forceMerge == obj.forceMerge) &&
                (resConstraint == obj.resConstraint) &&
                (advancedMultiWindowSave == obj.advancedMultiWindowSave) &&
                (subWindowAtts.equals(obj.subWindowAtts)));
    }

    // Property setting methods
    public void SetOutputToCurrentDirectory(boolean outputToCurrentDirectory_)
    {
        outputToCurrentDirectory = outputToCurrentDirectory_;
        Select(0);
    }

    public void SetOutputDirectory(String outputDirectory_)
    {
        outputDirectory = outputDirectory_;
        Select(1);
    }

    public void SetFileName(String fileName_)
    {
        fileName = fileName_;
        Select(2);
    }

    public void SetFamily(boolean family_)
    {
        family = family_;
        Select(3);
    }

    public void SetFormat(int format_)
    {
        format = format_;
        Select(4);
    }

    public void SetWidth(int width_)
    {
        width = width_;
        Select(5);
    }

    public void SetHeight(int height_)
    {
        height = height_;
        Select(6);
    }

    public void SetScreenCapture(boolean screenCapture_)
    {
        screenCapture = screenCapture_;
        Select(7);
    }

    public void SetSaveTiled(boolean saveTiled_)
    {
        saveTiled = saveTiled_;
        Select(8);
    }

    public void SetQuality(int quality_)
    {
        quality = quality_;
        Select(9);
    }

    public void SetProgressive(boolean progressive_)
    {
        progressive = progressive_;
        Select(10);
    }

    public void SetBinary(boolean binary_)
    {
        binary = binary_;
        Select(11);
    }

    public void SetLastRealFilename(String lastRealFilename_)
    {
        lastRealFilename = lastRealFilename_;
        Select(12);
    }

    public void SetStereo(boolean stereo_)
    {
        stereo = stereo_;
        Select(13);
    }

    public void SetCompression(int compression_)
    {
        compression = compression_;
        Select(14);
    }

    public void SetForceMerge(boolean forceMerge_)
    {
        forceMerge = forceMerge_;
        Select(15);
    }

    public void SetResConstraint(int resConstraint_)
    {
        resConstraint = resConstraint_;
        Select(16);
    }

    public void SetAdvancedMultiWindowSave(boolean advancedMultiWindowSave_)
    {
        advancedMultiWindowSave = advancedMultiWindowSave_;
        Select(17);
    }

    public void SetSubWindowAtts(SaveSubWindowsAttributes subWindowAtts_)
    {
        subWindowAtts = subWindowAtts_;
        Select(18);
    }

    // Property getting methods
    public boolean                  GetOutputToCurrentDirectory() { return outputToCurrentDirectory; }
    public String                   GetOutputDirectory() { return outputDirectory; }
    public String                   GetFileName() { return fileName; }
    public boolean                  GetFamily() { return family; }
    public int                      GetFormat() { return format; }
    public int                      GetWidth() { return width; }
    public int                      GetHeight() { return height; }
    public boolean                  GetScreenCapture() { return screenCapture; }
    public boolean                  GetSaveTiled() { return saveTiled; }
    public int                      GetQuality() { return quality; }
    public boolean                  GetProgressive() { return progressive; }
    public boolean                  GetBinary() { return binary; }
    public String                   GetLastRealFilename() { return lastRealFilename; }
    public boolean                  GetStereo() { return stereo; }
    public int                      GetCompression() { return compression; }
    public boolean                  GetForceMerge() { return forceMerge; }
    public int                      GetResConstraint() { return resConstraint; }
    public boolean                  GetAdvancedMultiWindowSave() { return advancedMultiWindowSave; }
    public SaveSubWindowsAttributes GetSubWindowAtts() { return subWindowAtts; }

    // Write and read methods.
    public void WriteAtts(CommunicationBuffer buf)
    {
        if(WriteSelect(0, buf))
            buf.WriteBool(outputToCurrentDirectory);
        if(WriteSelect(1, buf))
            buf.WriteString(outputDirectory);
        if(WriteSelect(2, buf))
            buf.WriteString(fileName);
        if(WriteSelect(3, buf))
            buf.WriteBool(family);
        if(WriteSelect(4, buf))
            buf.WriteInt(format);
        if(WriteSelect(5, buf))
            buf.WriteInt(width);
        if(WriteSelect(6, buf))
            buf.WriteInt(height);
        if(WriteSelect(7, buf))
            buf.WriteBool(screenCapture);
        if(WriteSelect(8, buf))
            buf.WriteBool(saveTiled);
        if(WriteSelect(9, buf))
            buf.WriteInt(quality);
        if(WriteSelect(10, buf))
            buf.WriteBool(progressive);
        if(WriteSelect(11, buf))
            buf.WriteBool(binary);
        if(WriteSelect(12, buf))
            buf.WriteString(lastRealFilename);
        if(WriteSelect(13, buf))
            buf.WriteBool(stereo);
        if(WriteSelect(14, buf))
            buf.WriteInt(compression);
        if(WriteSelect(15, buf))
            buf.WriteBool(forceMerge);
        if(WriteSelect(16, buf))
            buf.WriteInt(resConstraint);
        if(WriteSelect(17, buf))
            buf.WriteBool(advancedMultiWindowSave);
        if(WriteSelect(18, buf))
            subWindowAtts.Write(buf);
    }

    public void ReadAtts(int index, CommunicationBuffer buf)
    {
        switch(index)
        {
        case 0:
            SetOutputToCurrentDirectory(buf.ReadBool());
            break;
        case 1:
            SetOutputDirectory(buf.ReadString());
            break;
        case 2:
            SetFileName(buf.ReadString());
            break;
        case 3:
            SetFamily(buf.ReadBool());
            break;
        case 4:
            SetFormat(buf.ReadInt());
            break;
        case 5:
            SetWidth(buf.ReadInt());
            break;
        case 6:
            SetHeight(buf.ReadInt());
            break;
        case 7:
            SetScreenCapture(buf.ReadBool());
            break;
        case 8:
            SetSaveTiled(buf.ReadBool());
            break;
        case 9:
            SetQuality(buf.ReadInt());
            break;
        case 10:
            SetProgressive(buf.ReadBool());
            break;
        case 11:
            SetBinary(buf.ReadBool());
            break;
        case 12:
            SetLastRealFilename(buf.ReadString());
            break;
        case 13:
            SetStereo(buf.ReadBool());
            break;
        case 14:
            SetCompression(buf.ReadInt());
            break;
        case 15:
            SetForceMerge(buf.ReadBool());
            break;
        case 16:
            SetResConstraint(buf.ReadInt());
            break;
        case 17:
            SetAdvancedMultiWindowSave(buf.ReadBool());
            break;
        case 18:
            subWindowAtts.Read(buf);
            Select(18);
            break;
        }
    }

    public String toString(String indent)
    {
        String str = new String();
        str = str + boolToString("outputToCurrentDirectory", outputToCurrentDirectory, indent) + "\n";
        str = str + stringToString("outputDirectory", outputDirectory, indent) + "\n";
        str = str + stringToString("fileName", fileName, indent) + "\n";
        str = str + boolToString("family", family, indent) + "\n";
        str = str + indent + "format = ";
        if(format == FILEFORMAT_BMP)
            str = str + "FILEFORMAT_BMP";
        if(format == FILEFORMAT_CURVE)
            str = str + "FILEFORMAT_CURVE";
        if(format == FILEFORMAT_JPEG)
            str = str + "FILEFORMAT_JPEG";
        if(format == FILEFORMAT_OBJ)
            str = str + "FILEFORMAT_OBJ";
        if(format == FILEFORMAT_PNG)
            str = str + "FILEFORMAT_PNG";
        if(format == FILEFORMAT_POSTSCRIPT)
            str = str + "FILEFORMAT_POSTSCRIPT";
        if(format == FILEFORMAT_POVRAY)
            str = str + "FILEFORMAT_POVRAY";
        if(format == FILEFORMAT_PPM)
            str = str + "FILEFORMAT_PPM";
        if(format == FILEFORMAT_RGB)
            str = str + "FILEFORMAT_RGB";
        if(format == FILEFORMAT_STL)
            str = str + "FILEFORMAT_STL";
        if(format == FILEFORMAT_TIFF)
            str = str + "FILEFORMAT_TIFF";
        if(format == FILEFORMAT_ULTRA)
            str = str + "FILEFORMAT_ULTRA";
        if(format == FILEFORMAT_VTK)
            str = str + "FILEFORMAT_VTK";
        if(format == FILEFORMAT_PLY)
            str = str + "FILEFORMAT_PLY";
        str = str + "\n";
        str = str + intToString("width", width, indent) + "\n";
        str = str + intToString("height", height, indent) + "\n";
        str = str + boolToString("screenCapture", screenCapture, indent) + "\n";
        str = str + boolToString("saveTiled", saveTiled, indent) + "\n";
        str = str + intToString("quality", quality, indent) + "\n";
        str = str + boolToString("progressive", progressive, indent) + "\n";
        str = str + boolToString("binary", binary, indent) + "\n";
        str = str + stringToString("lastRealFilename", lastRealFilename, indent) + "\n";
        str = str + boolToString("stereo", stereo, indent) + "\n";
        str = str + indent + "compression = ";
        if(compression == COMPRESSIONTYPE_NONE)
            str = str + "COMPRESSIONTYPE_NONE";
        if(compression == COMPRESSIONTYPE_PACKBITS)
            str = str + "COMPRESSIONTYPE_PACKBITS";
        if(compression == COMPRESSIONTYPE_JPEG)
            str = str + "COMPRESSIONTYPE_JPEG";
        if(compression == COMPRESSIONTYPE_DEFLATE)
            str = str + "COMPRESSIONTYPE_DEFLATE";
        str = str + "\n";
        str = str + boolToString("forceMerge", forceMerge, indent) + "\n";
        str = str + indent + "resConstraint = ";
        if(resConstraint == RESCONSTRAINT_NOCONSTRAINT)
            str = str + "RESCONSTRAINT_NOCONSTRAINT";
        if(resConstraint == RESCONSTRAINT_EQUALWIDTHHEIGHT)
            str = str + "RESCONSTRAINT_EQUALWIDTHHEIGHT";
        if(resConstraint == RESCONSTRAINT_SCREENPROPORTIONS)
            str = str + "RESCONSTRAINT_SCREENPROPORTIONS";
        str = str + "\n";
        str = str + boolToString("advancedMultiWindowSave", advancedMultiWindowSave, indent) + "\n";
        str = str + indent + "subWindowAtts = {\n" + subWindowAtts.toString(indent + "    ") + indent + "}\n";
        return str;
    }


    // Attributes
    private boolean                  outputToCurrentDirectory;
    private String                   outputDirectory;
    private String                   fileName;
    private boolean                  family;
    private int                      format;
    private int                      width;
    private int                      height;
    private boolean                  screenCapture;
    private boolean                  saveTiled;
    private int                      quality;
    private boolean                  progressive;
    private boolean                  binary;
    private String                   lastRealFilename;
    private boolean                  stereo;
    private int                      compression;
    private boolean                  forceMerge;
    private int                      resConstraint;
    private boolean                  advancedMultiWindowSave;
    private SaveSubWindowsAttributes subWindowAtts;
}

