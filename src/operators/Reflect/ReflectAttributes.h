#ifndef REFLECTATTRIBUTES_H
#define REFLECTATTRIBUTES_H
#include <string>
#include <AttributeSubject.h>

// ****************************************************************************
// Class: ReflectAttributes
//
// Purpose:
//    This class contains attributes for the reflect operator.
//
// Notes:      Autogenerated by xml2atts.
//
// Programmer: xml2atts
// Creation:   Tue May 20 14:50:21 PST 2003
//
// Modifications:
//   
// ****************************************************************************

class ReflectAttributes : public AttributeSubject
{
public:
    enum Octant
    {
        PXPYPZ,
        NXPYPZ,
        PXNYPZ,
        NXNYPZ,
        PXPYNZ,
        NXPYNZ,
        PXNYNZ,
        NXNYNZ
    };

    ReflectAttributes();
    ReflectAttributes(const ReflectAttributes &obj);
    virtual ~ReflectAttributes();

    virtual void operator = (const ReflectAttributes &obj);
    virtual bool operator == (const ReflectAttributes &obj) const;
    virtual bool operator != (const ReflectAttributes &obj) const;

    virtual const std::string TypeName() const;
    virtual bool CopyAttributes(const AttributeGroup *);
    virtual AttributeSubject *CreateCompatible(const std::string &) const;

    // Property selection methods
    virtual void SelectAll();
    void SelectReflections();

    // Property setting methods
    void SetOctant(Octant octant_);
    void SetUseXBoundary(bool useXBoundary_);
    void SetSpecifiedX(double specifiedX_);
    void SetUseYBoundary(bool useYBoundary_);
    void SetSpecifiedY(double specifiedY_);
    void SetUseZBoundary(bool useZBoundary_);
    void SetSpecifiedZ(double specifiedZ_);
    void SetReflections(const int *reflections_);

    // Property getting methods
    Octant    GetOctant() const;
    bool      GetUseXBoundary() const;
    double    GetSpecifiedX() const;
    bool      GetUseYBoundary() const;
    double    GetSpecifiedY() const;
    bool      GetUseZBoundary() const;
    double    GetSpecifiedZ() const;
    const int *GetReflections() const;
          int *GetReflections();

    // Persistence methods
    virtual bool CreateNode(DataNode *node, bool forceAdd);
    virtual void SetFromNode(DataNode *node);

    // Enum conversion functions
    static std::string Octant_ToString(Octant);
    static bool Octant_FromString(const std::string &, Octant &);
protected:
    static std::string Octant_ToString(int);
public:

    // Keyframing methods
    virtual std::string               GetFieldName(int index) const;
    virtual AttributeGroup::FieldType GetFieldType(int index) const;
    virtual std::string               GetFieldTypeName(int index) const;
    virtual bool                      FieldsEqual(int index, const AttributeGroup *rhs) const;

    // User-defined methods
    virtual void ProcessOldVersions(DataNode *node, const char *configVersion);
private:
    int    octant;
    bool   useXBoundary;
    double specifiedX;
    bool   useYBoundary;
    double specifiedY;
    bool   useZBoundary;
    double specifiedZ;
    int    reflections[8];
};

#endif
