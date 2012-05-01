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
  \file    MultiplexOut.h
  \author    Jens Krueger
        SCI Institute
        University of Utah
  \version  1.0
  \date    September 2008
*/


#pragma once

#ifndef MULTIPLEXOUT_H
#define MULTIPLEXOUT_H

#include <vector>
#include "AbstrDebugOut.h"

class MultiplexOut : public AbstrDebugOut {
  public:
    MultiplexOut() {}
    ~MultiplexOut();

    void AddDebugOut(AbstrDebugOut* pDebugger);
    void RemoveDebugOut(AbstrDebugOut* pDebugger);

    virtual void printf(const char* format, ...) const;
    virtual void Message(const char* source, const char* format, ...);
    virtual void Warning(const char* source, const char* format, ...);
    virtual void Error(const char* source, const char* format, ...);

    virtual void SetShowMessages(bool bShowMessages);
    virtual void SetShowWarnings(bool bShowWarnings);
    virtual void SetShowErrors(bool bShowErrors);
    virtual void SetShowOther(bool bShowOther);

    size_t size() const { return m_vpDebugger.size(); }
    void clear();

  private:
    std::vector<AbstrDebugOut*> m_vpDebugger;
};
#endif // MULTIPLEXOUT_H
