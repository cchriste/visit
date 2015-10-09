/*****************************************************************************
*
* Copyright (c) 2000 - 2015, Lawrence Livermore National Security, LLC
* Produced at the Lawrence Livermore National Laboratory
* LLNL-CODE-442911
* All rights reserved.
*
* This file is  part of VisIt. For  details, see https://visit.llnl.gov/.  The
* full copyright notice is contained in the file COPYRIGHT located at the root
* of the VisIt distribution or at http://www.llnl.gov/visit/copyright.html.
*
* Redistribution  and  use  in  source  and  binary  forms,  with  or  without
* modification, are permitted provided that the following conditions are met:
*
*  - Redistributions of  source code must  retain the above  copyright notice,
*    this list of conditions and the disclaimer below.
*  - Redistributions in binary form must reproduce the above copyright notice,
*    this  list of  conditions  and  the  disclaimer (as noted below)  in  the
*    documentation and/or other materials provided with the distribution.
*  - Neither the name of  the LLNS/LLNL nor the names of  its contributors may
*    be used to endorse or promote products derived from this software without
*    specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT  HOLDERS AND CONTRIBUTORS "AS IS"
* AND ANY EXPRESS OR  IMPLIED WARRANTIES, INCLUDING,  BUT NOT  LIMITED TO, THE
* IMPLIED WARRANTIES OF MERCHANTABILITY AND  FITNESS FOR A PARTICULAR  PURPOSE
* ARE  DISCLAIMED. IN  NO EVENT  SHALL LAWRENCE  LIVERMORE NATIONAL  SECURITY,
* LLC, THE  U.S.  DEPARTMENT OF  ENERGY  OR  CONTRIBUTORS BE  LIABLE  FOR  ANY
* DIRECT,  INDIRECT,   INCIDENTAL,   SPECIAL,   EXEMPLARY,  OR   CONSEQUENTIAL
* DAMAGES (INCLUDING, BUT NOT  LIMITED TO, PROCUREMENT OF  SUBSTITUTE GOODS OR
* SERVICES; LOSS OF  USE, DATA, OR PROFITS; OR  BUSINESS INTERRUPTION) HOWEVER
* CAUSED  AND  ON  ANY  THEORY  OF  LIABILITY,  WHETHER  IN  CONTRACT,  STRICT
* LIABILITY, OR TORT  (INCLUDING NEGLIGENCE OR OTHERWISE)  ARISING IN ANY  WAY
* OUT OF THE  USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
* DAMAGE.
*
*****************************************************************************/

#ifndef VIEW3DATTRIBUTES_H
#define VIEW3DATTRIBUTES_H
#include <state_exports.h>
#include <AttributeSubject.h>


// ****************************************************************************
// Class: View3DAttributes
//
// Purpose:
//    This class contains the 3d view attributes.
//
// Notes:      Autogenerated by xml2atts.
//
// Programmer: xml2atts
// Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

class STATE_API View3DAttributes : public AttributeSubject
{
public:
    // These constructors are for objects of this class
    View3DAttributes();
    View3DAttributes(const View3DAttributes &obj);
protected:
    // These constructors are for objects derived from this class
    View3DAttributes(private_tmfs_t tmfs);
    View3DAttributes(const View3DAttributes &obj, private_tmfs_t tmfs);
public:
    virtual ~View3DAttributes();

    virtual View3DAttributes& operator = (const View3DAttributes &obj);
    virtual bool operator == (const View3DAttributes &obj) const;
    virtual bool operator != (const View3DAttributes &obj) const;
private:
    void Init();
    void Copy(const View3DAttributes &obj);
public:

    virtual const std::string TypeName() const;
    virtual bool CopyAttributes(const AttributeGroup *);
    virtual AttributeSubject *CreateCompatible(const std::string &) const;
    virtual AttributeSubject *NewInstance(bool) const;

    // Property selection methods
    virtual void SelectAll();
    void SelectViewNormal();
    void SelectFocus();
    void SelectViewUp();
    void SelectImagePan();
    void SelectCenterOfRotation();
    void SelectAxis3DScales();
    void SelectShear();

    // Property setting methods
    void SetViewNormal(const double *viewNormal_);
    void SetFocus(const double *focus_);
    void SetViewUp(const double *viewUp_);
    void SetViewAngle(double viewAngle_);
    void SetParallelScale(double parallelScale_);
    void SetNearPlane(double nearPlane_);
    void SetFarPlane(double farPlane_);
    void SetImagePan(const double *imagePan_);
    void SetImageZoom(double imageZoom_);
    void SetPerspective(bool perspective_);
    void SetEyeAngle(double eyeAngle_);
    void SetCenterOfRotationSet(bool centerOfRotationSet_);
    void SetCenterOfRotation(const double *centerOfRotation_);
    void SetAxis3DScaleFlag(bool axis3DScaleFlag_);
    void SetAxis3DScales(const double *axis3DScales_);
    void SetShear(const double *shear_);
    void SetWindowValid(bool windowValid_);

    // Property getting methods
    const double *GetViewNormal() const;
          double *GetViewNormal();
    const double *GetFocus() const;
          double *GetFocus();
    const double *GetViewUp() const;
          double *GetViewUp();
    double       GetViewAngle() const;
    double       GetParallelScale() const;
    double       GetNearPlane() const;
    double       GetFarPlane() const;
    const double *GetImagePan() const;
          double *GetImagePan();
    double       GetImageZoom() const;
    bool         GetPerspective() const;
    double       GetEyeAngle() const;
    bool         GetCenterOfRotationSet() const;
    const double *GetCenterOfRotation() const;
          double *GetCenterOfRotation();
    bool         GetAxis3DScaleFlag() const;
    const double *GetAxis3DScales() const;
          double *GetAxis3DScales();
    const double *GetShear() const;
          double *GetShear();
    bool         GetWindowValid() const;

    // Persistence methods
    virtual bool CreateNode(DataNode *node, bool completeSave, bool forceAdd);
    virtual void SetFromNode(DataNode *node);


    // Keyframing methods
    virtual std::string               GetFieldName(int index) const;
    virtual AttributeGroup::FieldType GetFieldType(int index) const;
    virtual std::string               GetFieldTypeName(int index) const;
    virtual bool                      FieldsEqual(int index, const AttributeGroup *rhs) const;

    // User-defined methods
    void RotateAxis(int axis, double angle);
    void ResetView(const double *bbox);

    // IDs that can be used to identify fields in case statements
    enum {
        ID_viewNormal = 0,
        ID_focus,
        ID_viewUp,
        ID_viewAngle,
        ID_parallelScale,
        ID_nearPlane,
        ID_farPlane,
        ID_imagePan,
        ID_imageZoom,
        ID_perspective,
        ID_eyeAngle,
        ID_centerOfRotationSet,
        ID_centerOfRotation,
        ID_axis3DScaleFlag,
        ID_axis3DScales,
        ID_shear,
        ID_windowValid,
        ID__LAST
    };

private:
    double viewNormal[3];
    double focus[3];
    double viewUp[3];
    double viewAngle;
    double parallelScale;
    double nearPlane;
    double farPlane;
    double imagePan[2];
    double imageZoom;
    bool   perspective;
    double eyeAngle;
    bool   centerOfRotationSet;
    double centerOfRotation[3];
    bool   axis3DScaleFlag;
    double axis3DScales[3];
    double shear[3];
    bool   windowValid;

    // Static class format string for type map.
    static const char *TypeMapFormatString;
    static const private_tmfs_t TmfsStruct;
};
#define VIEW3DATTRIBUTES_TMFS "DDDddddDdbdbDbDDb"

#endif
