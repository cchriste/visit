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
  \file    TextfileOut.h
  \author    Jens Krueger
        SCI Institute
        University of Utah
  \version  1.0
  \date    August 2008
*/


#pragma once

#ifndef TEXTFILEOUT_H
#define TEXTFILEOUT_H

#include <string>
#include "AbstrDebugOut.h"

class TextfileOut : public AbstrDebugOut{
  public:
    TextfileOut(std::string strFilename="logfile.txt");
    ~TextfileOut();
    virtual void printf(const char* format, ...) const;
    virtual void Message(const char* source, const char* format, ...);
    virtual void Warning(const char* source, const char* format, ...);
    virtual void Error(const char* source, const char* format, ...);

    const std::string& GetFileName() const {return m_strFilename;}

  private:
    TextfileOut(const TextfileOut &); ///< unimplemented.

  private:
    std::string m_strFilename;

    /// same as printf above but does regard m_bShowOther
    void _printf(const char* format, ...) const;
};

#endif // TEXTFILEOUT_H
