/*
   For more information, please see: http://software.sci.utah.edu

   The MIT License

   Copyright (c) 2008 Scientific Computing and Imaging Institute,
   University of Utah.

   
   Permission is hereby granted, free of charge, to any person obtaining a
   copy of this software and associated documentation files (the "Software"),
   to deal in the Software without restriction, including without limitation
   the rights to use, copy, modify, merge, publish, distribute, sublicense,
   and/or sell copies of the Software, and to permit persons to whom the
   Software is furnished to do so, subject to the following conditions:

   The above copyright notice and this permission notice shall be included
   in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
   OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
   THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
   FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
   DEALINGS IN THE SOFTWARE.
*/

/**
  \file    GLRaycaster.h
  \author  Jens Krueger
           SCI Institute
           University of Utah
  \version 1.0
  \date    November 2008
*/


#pragma once

#ifndef GLRAYCASTER_H
#define GLRAYCASTER_H

#include "../Renderer/GLRenderer.h"

/** \class GLRaycaster
 * GPU Rayster.
 *
 * GLRaycaster is a GPU based raycaster for volumetric scalar data which uses GLSL. */
class GLRaycaster : public GLRenderer {
  public:
    /** Constructs a VRer with immediate redraw, and
     * wireframe mode off.
     * \param pMasterController message routing object */
    GLRaycaster(MasterController* pMasterController, bool bUseOnlyPowerOfTwo);
    virtual ~GLRaycaster();

    /** Loads GLSL vertex and fragment shaders. */
    virtual bool Initialize();

    /** Deallocates GPU memory allocated during the rendering process. */
    virtual void Cleanup();

    virtual bool SupportsClearView() {return true;}

  protected:
    GLFBOTex*       m_pFBORayEntry;
    GLSLProgram*    m_pProgramRenderFrontFaces;
    GLSLProgram*    m_pProgramIso2;

    void SetBrickDepShaderVars(size_t iCurrentBrick);

    virtual void CreateOffscreenBuffers();
    void RenderBox(const FLOATVECTOR3& vCenter, const FLOATVECTOR3& vExtend, const FLOATVECTOR3& vMinCoords, const FLOATVECTOR3& vMaxCoords, bool bCullBack);

    virtual void Render3DPreLoop();
    virtual void Render3DInLoop(size_t iCurrentBrick);
    virtual void Render3DPostLoop();

    virtual void StartFrame();
    virtual void SetDataDepShaderVars();

    FLOATMATRIX4 ComputeEyeToTextureMatrix(FLOATVECTOR3 p1, FLOATVECTOR3 t1, 
                                           FLOATVECTOR3 p2, FLOATVECTOR3 t2);
};

#endif // GLRAYCASTER_H
