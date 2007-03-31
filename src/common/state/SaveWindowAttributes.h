#ifndef SAVEWINDOWATTRIBUTES_H
#define SAVEWINDOWATTRIBUTES_H
#include <state_exports.h>
#include <string>
#include <AttributeSubject.h>

// ****************************************************************************
// Class: SaveWindowAttributes
//
// Purpose:
//    This class contains the attributes used for saving windows.
//
// Notes:      Autogenerated by xml2atts.
//
// Programmer: xml2atts
// Creation:   Tue Jul 15 13:47:56 PST 2003
//
// Modifications:
//   
// ****************************************************************************

class STATE_API SaveWindowAttributes : public AttributeSubject
{
public:
    enum FileFormat
    {
        BMP,
        CURVE,
        JPEG,
        OBJ,
        PNG,
        POSTSCRIPT,
        PPM,
        RGB,
        STL,
        TIFF,
        ULTRA,
        VTK
    };

    SaveWindowAttributes();
    SaveWindowAttributes(const SaveWindowAttributes &obj);
    virtual ~SaveWindowAttributes();

    virtual void operator = (const SaveWindowAttributes &obj);
    virtual bool operator == (const SaveWindowAttributes &obj) const;
    virtual bool operator != (const SaveWindowAttributes &obj) const;

    virtual const std::string TypeName() const;
    virtual bool CopyAttributes(const AttributeGroup *);
    virtual AttributeSubject *CreateCompatible(const std::string &) const;

    // Property selection methods
    virtual void SelectAll();
    void SelectHostName();
    void SelectFileName();
    void SelectLastRealFilename();

    // Property setting methods
    void SetHostName(const std::string &hostName_);
    void SetFileName(const std::string &fileName_);
    void SetFamily(bool family_);
    void SetFormat(FileFormat format_);
    void SetMaintainAspect(bool maintainAspect_);
    void SetWidth(int width_);
    void SetHeight(int height_);
    void SetScreenCapture(bool screenCapture_);
    void SetSaveTiled(bool saveTiled_);
    void SetQuality(int quality_);
    void SetProgressive(bool progressive_);
    void SetBinary(bool binary_);
    void SetLastRealFilename(const std::string &lastRealFilename_);

    // Property getting methods
    const std::string &GetHostName() const;
          std::string &GetHostName();
    const std::string &GetFileName() const;
          std::string &GetFileName();
    bool              GetFamily() const;
    FileFormat        GetFormat() const;
    bool              GetMaintainAspect() const;
    int               GetWidth() const;
    int               GetHeight() const;
    bool              GetScreenCapture() const;
    bool              GetSaveTiled() const;
    int               GetQuality() const;
    bool              GetProgressive() const;
    bool              GetBinary() const;
    const std::string &GetLastRealFilename() const;
          std::string &GetLastRealFilename();

    // Persistence methods
    virtual bool CreateNode(DataNode *node, bool forceAdd);
    virtual void SetFromNode(DataNode *node);

    // Enum conversion functions
    static std::string FileFormat_ToString(FileFormat);
    static bool FileFormat_FromString(const std::string &, FileFormat &);
protected:
    static std::string FileFormat_ToString(int);
public:

    // Keyframing methods
    virtual std::string               GetFieldName(int index) const;
    virtual AttributeGroup::FieldType GetFieldType(int index) const;
    virtual std::string               GetFieldTypeName(int index) const;
    virtual bool                      FieldsEqual(int index, const AttributeGroup *rhs) const;

private:
    std::string hostName;
    std::string fileName;
    bool        family;
    int         format;
    bool        maintainAspect;
    int         width;
    int         height;
    bool        screenCapture;
    bool        saveTiled;
    int         quality;
    bool        progressive;
    bool        binary;
    std::string lastRealFilename;
};

#endif
