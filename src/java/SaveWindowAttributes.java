// ***************************************************************************
//
// Copyright (c) 2000 - 2007, The Regents of the University of California
// Produced at the Lawrence Livermore National Laboratory
// All rights reserved.
//
// This file is part of VisIt. For details, see http://www.llnl.gov/visit/. The
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
//    documentation and/or materials provided with the distribution.
//  - Neither the name of the UC/LLNL nor  the names of its contributors may be
//    used to  endorse or  promote products derived from  this software without
//    specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT  HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR  IMPLIED WARRANTIES, INCLUDING,  BUT NOT  LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND  FITNESS FOR A PARTICULAR  PURPOSE
// ARE  DISCLAIMED.  IN  NO  EVENT  SHALL  THE  REGENTS  OF  THE  UNIVERSITY OF
// CALIFORNIA, THE U.S.  DEPARTMENT  OF  ENERGY OR CONTRIBUTORS BE  LIABLE  FOR
// ANY  DIRECT,  INDIRECT,  INCIDENTAL,  SPECIAL,  EXEMPLARY,  OR CONSEQUENTIAL
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
// Creation:   Wed Mar 14 17:55:13 PST 2007
//
// Modifications:
//   
// ****************************************************************************

public class SaveWindowAttributes extends AttributeSubject
{
    // Enum values
    public final static int FILEFORMAT_BMP = 0;
    public final static int FILEFORMAT_CURVE = 1;
    public final static int FILEFORMAT_JPEG = 2;
    public final static int FILEFORMAT_OBJ = 3;
    public final static int FILEFORMAT_PNG = 4;
    public final static int FILEFORMAT_POSTSCRIPT = 5;
    public final static int FILEFORMAT_PPM = 6;
    public final static int FILEFORMAT_RGB = 7;
    public final static int FILEFORMAT_STL = 8;
    public final static int FILEFORMAT_TIFF = 9;
    public final static int FILEFORMAT_ULTRA = 10;
    public final static int FILEFORMAT_VTK = 11;

    public final static int COMPRESSIONTYPE_NONE = 0;
    public final static int COMPRESSIONTYPE_PACKBITS = 1;
    public final static int COMPRESSIONTYPE_JPEG = 2;
    public final static int COMPRESSIONTYPE_DEFLATE = 3;


    public SaveWindowAttributes()
    {
        super(16);

        outputToCurrentDirectory = true;
        outputDirectory = new String(".");
        fileName = new String("visit");
        family = true;
        format = FILEFORMAT_TIFF;
        maintainAspect = true;
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
    }

    public SaveWindowAttributes(SaveWindowAttributes obj)
    {
        super(16);

        outputToCurrentDirectory = obj.outputToCurrentDirectory;
        outputDirectory = new String(obj.outputDirectory);
        fileName = new String(obj.fileName);
        family = obj.family;
        format = obj.format;
        maintainAspect = obj.maintainAspect;
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

        SelectAll();
    }

    public boolean equals(SaveWindowAttributes obj)
    {
        // Create the return value
        return ((outputToCurrentDirectory == obj.outputToCurrentDirectory) &&
                (outputDirectory == obj.outputDirectory) &&
                (fileName == obj.fileName) &&
                (family == obj.family) &&
                (format == obj.format) &&
                (maintainAspect == obj.maintainAspect) &&
                (width == obj.width) &&
                (height == obj.height) &&
                (screenCapture == obj.screenCapture) &&
                (saveTiled == obj.saveTiled) &&
                (quality == obj.quality) &&
                (progressive == obj.progressive) &&
                (binary == obj.binary) &&
                (lastRealFilename == obj.lastRealFilename) &&
                (stereo == obj.stereo) &&
                (compression == obj.compression));
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

    public void SetMaintainAspect(boolean maintainAspect_)
    {
        maintainAspect = maintainAspect_;
        Select(5);
    }

    public void SetWidth(int width_)
    {
        width = width_;
        Select(6);
    }

    public void SetHeight(int height_)
    {
        height = height_;
        Select(7);
    }

    public void SetScreenCapture(boolean screenCapture_)
    {
        screenCapture = screenCapture_;
        Select(8);
    }

    public void SetSaveTiled(boolean saveTiled_)
    {
        saveTiled = saveTiled_;
        Select(9);
    }

    public void SetQuality(int quality_)
    {
        quality = quality_;
        Select(10);
    }

    public void SetProgressive(boolean progressive_)
    {
        progressive = progressive_;
        Select(11);
    }

    public void SetBinary(boolean binary_)
    {
        binary = binary_;
        Select(12);
    }

    public void SetLastRealFilename(String lastRealFilename_)
    {
        lastRealFilename = lastRealFilename_;
        Select(13);
    }

    public void SetStereo(boolean stereo_)
    {
        stereo = stereo_;
        Select(14);
    }

    public void SetCompression(int compression_)
    {
        compression = compression_;
        Select(15);
    }

    // Property getting methods
    public boolean GetOutputToCurrentDirectory() { return outputToCurrentDirectory; }
    public String  GetOutputDirectory() { return outputDirectory; }
    public String  GetFileName() { return fileName; }
    public boolean GetFamily() { return family; }
    public int     GetFormat() { return format; }
    public boolean GetMaintainAspect() { return maintainAspect; }
    public int     GetWidth() { return width; }
    public int     GetHeight() { return height; }
    public boolean GetScreenCapture() { return screenCapture; }
    public boolean GetSaveTiled() { return saveTiled; }
    public int     GetQuality() { return quality; }
    public boolean GetProgressive() { return progressive; }
    public boolean GetBinary() { return binary; }
    public String  GetLastRealFilename() { return lastRealFilename; }
    public boolean GetStereo() { return stereo; }
    public int     GetCompression() { return compression; }

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
            buf.WriteBool(maintainAspect);
        if(WriteSelect(6, buf))
            buf.WriteInt(width);
        if(WriteSelect(7, buf))
            buf.WriteInt(height);
        if(WriteSelect(8, buf))
            buf.WriteBool(screenCapture);
        if(WriteSelect(9, buf))
            buf.WriteBool(saveTiled);
        if(WriteSelect(10, buf))
            buf.WriteInt(quality);
        if(WriteSelect(11, buf))
            buf.WriteBool(progressive);
        if(WriteSelect(12, buf))
            buf.WriteBool(binary);
        if(WriteSelect(13, buf))
            buf.WriteString(lastRealFilename);
        if(WriteSelect(14, buf))
            buf.WriteBool(stereo);
        if(WriteSelect(15, buf))
            buf.WriteInt(compression);
    }

    public void ReadAtts(int n, CommunicationBuffer buf)
    {
        for(int i = 0; i < n; ++i)
        {
            int index = (int)buf.ReadByte();
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
                SetMaintainAspect(buf.ReadBool());
                break;
            case 6:
                SetWidth(buf.ReadInt());
                break;
            case 7:
                SetHeight(buf.ReadInt());
                break;
            case 8:
                SetScreenCapture(buf.ReadBool());
                break;
            case 9:
                SetSaveTiled(buf.ReadBool());
                break;
            case 10:
                SetQuality(buf.ReadInt());
                break;
            case 11:
                SetProgressive(buf.ReadBool());
                break;
            case 12:
                SetBinary(buf.ReadBool());
                break;
            case 13:
                SetLastRealFilename(buf.ReadString());
                break;
            case 14:
                SetStereo(buf.ReadBool());
                break;
            case 15:
                SetCompression(buf.ReadInt());
                break;
            }
        }
    }


    // Attributes
    private boolean outputToCurrentDirectory;
    private String  outputDirectory;
    private String  fileName;
    private boolean family;
    private int     format;
    private boolean maintainAspect;
    private int     width;
    private int     height;
    private boolean screenCapture;
    private boolean saveTiled;
    private int     quality;
    private boolean progressive;
    private boolean binary;
    private String  lastRealFilename;
    private boolean stereo;
    private int     compression;
}

