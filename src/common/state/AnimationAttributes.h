#ifndef ANIMATIONATTRIBUTES_H
#define ANIMATIONATTRIBUTES_H
#include <state_exports.h>
#include <string>
#include <AttributeSubject.h>

// ****************************************************************************
// Class: AnimationAttributes
//
// Purpose:
//    This class contains the animation attributes.
//
// Notes:      Autogenerated by xml2atts.
//
// Programmer: xml2atts
// Creation:   Thu Dec 18 11:23:54 PDT 2003
//
// Modifications:
//   
// ****************************************************************************

class STATE_API AnimationAttributes : public AttributeSubject
{
public:
    enum PlaybackMode
    {
        Looping,
        PlayOnce,
        Swing
    };

    AnimationAttributes();
    AnimationAttributes(const AnimationAttributes &obj);
    virtual ~AnimationAttributes();

    virtual void operator = (const AnimationAttributes &obj);
    virtual bool operator == (const AnimationAttributes &obj) const;
    virtual bool operator != (const AnimationAttributes &obj) const;

    virtual const std::string TypeName() const;
    virtual bool CopyAttributes(const AttributeGroup *);
    virtual AttributeSubject *CreateCompatible(const std::string &) const;
    virtual AttributeSubject *NewInstance(bool) const;

    // Property selection methods
    virtual void SelectAll();

    // Property setting methods
    void SetPipelineCachingMode(bool pipelineCachingMode_);
    void SetTimeout(int timeout_);
    void SetPlaybackMode(PlaybackMode playbackMode_);

    // Property getting methods
    bool GetPipelineCachingMode() const;
    int  GetTimeout() const;
    PlaybackMode GetPlaybackMode() const;

    // Persistence methods
    virtual bool CreateNode(DataNode *node, bool completeSave, bool forceAdd);
    virtual void SetFromNode(DataNode *node);

    // Enum conversion functions
    static std::string PlaybackMode_ToString(PlaybackMode);
    static bool PlaybackMode_FromString(const std::string &, PlaybackMode &);
protected:
    static std::string PlaybackMode_ToString(int);
public:

    // Keyframing methods
    virtual std::string               GetFieldName(int index) const;
    virtual AttributeGroup::FieldType GetFieldType(int index) const;
    virtual std::string               GetFieldTypeName(int index) const;
    virtual bool                      FieldsEqual(int index, const AttributeGroup *rhs) const;

private:
    bool pipelineCachingMode;
    int  timeout;
    int  playbackMode;
};

#endif
