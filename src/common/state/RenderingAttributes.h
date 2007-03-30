#ifndef RENDERINGATTRIBUTES_H
#define RENDERINGATTRIBUTES_H
#include <state_exports.h>
#include <string>
#include <AttributeSubject.h>

// ****************************************************************************
// Class: RenderingAttributes
//
// Purpose:
//    This class contains special rendering attributes like antialiasing and stero settings.
//
// Notes:      Autogenerated by xml2atts.
//
// Programmer: xml2atts
// Creation:   Tue May 20 13:40:07 PST 2003
//
// Modifications:
//   
// ****************************************************************************

class STATE_API RenderingAttributes : public AttributeSubject
{
public:
    enum GeometryRepresentation
    {
        Surfaces,
        Wireframe,
        Points
    };
    enum StereoTypes
    {
        RedBlue,
        Interlaced,
        CrystalEyes
    };

    RenderingAttributes();
    RenderingAttributes(const RenderingAttributes &obj);
    virtual ~RenderingAttributes();

    virtual void operator = (const RenderingAttributes &obj);
    virtual bool operator == (const RenderingAttributes &obj) const;
    virtual bool operator != (const RenderingAttributes &obj) const;

    virtual const std::string TypeName() const;
    virtual bool CopyAttributes(const AttributeGroup *);
    virtual AttributeSubject *CreateCompatible(const std::string &) const;

    // Property selection methods
    virtual void SelectAll();

    // Property setting methods
    void SetAntialiasing(bool antialiasing_);
    void SetGeometryRepresentation(GeometryRepresentation geometryRepresentation_);
    void SetDisplayLists(bool displayLists_);
    void SetStereoRendering(bool stereoRendering_);
    void SetStereoType(StereoTypes stereoType_);
    void SetNotifyForEachRender(bool notifyForEachRender_);
    void SetScalableRendering(bool scalableRendering_);
    void SetScalableThreshold(int scalableThreshold_);

    // Property getting methods
    bool GetAntialiasing() const;
    GeometryRepresentation GetGeometryRepresentation() const;
    bool GetDisplayLists() const;
    bool GetStereoRendering() const;
    StereoTypes GetStereoType() const;
    bool GetNotifyForEachRender() const;
    bool GetScalableRendering() const;
    int  GetScalableThreshold() const;

    // Persistence methods
    virtual bool CreateNode(DataNode *node, bool forceAdd);
    virtual void SetFromNode(DataNode *node);

    // Enum conversion functions
    static std::string GeometryRepresentation_ToString(GeometryRepresentation);
    static bool GeometryRepresentation_FromString(const std::string &, GeometryRepresentation &);
protected:
    static std::string GeometryRepresentation_ToString(int);
public:
    static std::string StereoTypes_ToString(StereoTypes);
    static bool StereoTypes_FromString(const std::string &, StereoTypes &);
protected:
    static std::string StereoTypes_ToString(int);
public:

    // Keyframing methods
    virtual std::string               GetFieldName(int index) const;
    virtual AttributeGroup::FieldType GetFieldType(int index) const;
    virtual std::string               GetFieldTypeName(int index) const;
    virtual bool                      FieldsEqual(int index, const AttributeGroup *rhs) const;

private:
    bool antialiasing;
    int  geometryRepresentation;
    bool displayLists;
    bool stereoRendering;
    int  stereoType;
    bool notifyForEachRender;
    bool scalableRendering;
    int  scalableThreshold;
};

#endif
