/*****************************************************************************
*
* Copyright (c) 2000 - 2013, Lawrence Livermore National Security, LLC
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

// ************************************************************************* //
//                              StringHelpers.h                              //
// ************************************************************************* //

#ifndef STRINGHELPERS_H
#define STRINGHELPERS_H
#include <utility_exports.h>

#include <iostream>
#include <set>
#include <string>
#include <vector>
#include <cerrno>
#include <cstdio>
#include <cstdlib>
#include <snprintf.h>
#include <visit-config.h> // FOR VISIT_SLASH_STRING

#if __GNUC__ >= 3
#   define MUST_CHECK __attribute__ ((warn_unused_result))
#else
#   define MUST_CHECK /*nothing*/
#endif



namespace StringHelpers
{
    const std::string NON_RELEVANT_CHARS =
      "`~!@#$%^&*()-_=+{[}]|\\:;\"'<,>.?/0123456789";

    enum FindResult {FindNone = -1, FindError = -2};

    void UTILITY_API GroupStrings(std::vector<std::string> stringList,
                       std::vector<std::vector<std::string> > &stringGroups,
                       std::vector<std::string> &groupNames,
                       int numLeadingVals = 3,
                       std::string nonRelevantChars = NON_RELEVANT_CHARS);
    void UTILITY_API GroupStringsAsPaths(std::vector<std::string> stringList,
                       std::vector<std::vector<std::string> > &stringGroups,
                       std::vector<std::string> &groupNames);
    void UTILITY_API GroupStringsFixedAlpha(std::vector<std::string> stringList,
                       int numGroups,
                       std::vector<std::vector<std::string> > &stringGroups);
    void UTILITY_API GroupStringsFixedAlpha(
                       const std::set<std::string> &stringList,
                       int numGroups,
                       std::vector<std::set<std::string> > &stringGroups);

    int UTILITY_API FindRE(const std::string &s, const std::string &re);
    int UTILITY_API FindRE(const char *stringToSearch, const char *re);
    bool UTILITY_API ReplaceRE(std::string &s, const std::string &re,
                               const std::string &repl);
    std::string UTILITY_API ExtractRESubstr(const char *stringToSearch,
                                            const char *re);

    bool UTILITY_API ValidatePrintfFormatString(const char *fmtStr,
                                                const char *arg1Type, ...);

    const char UTILITY_API *Basename(const char *path);
    std::string UTILITY_API Basename(std::string const path);
    const char UTILITY_API *Dirname(const char *path);
    std::string UTILITY_API  Dirname(std::string const path);
    const char UTILITY_API *Absname(const char *cwd_context, 
                                const char *path,
                                const char *pathSep = VISIT_SLASH_STRING);
    std::string UTILITY_API Absname(std::string const cwd_context, 
                                std::string const path,
                                std::string const pathSep = VISIT_SLASH_STRING);
    const char UTILITY_API *Normalize(const char *path,
                                const char *pathSep = VISIT_SLASH_STRING);
    std::string UTILITY_API Normalize(const std::string& path,
                                std::string const pathSep = VISIT_SLASH_STRING);

    std::string UTILITY_API car(const std::string, const char separator);
    std::string UTILITY_API cdr(const std::string, const char separator);
    void UTILITY_API append(std::vector<std::string> &,
                            std::vector<std::string>);
    std::vector<std::string> UTILITY_API split(const std::string,
                                               const char separator);
    void UTILITY_API rtrim(std::string &var);
    void UTILITY_API ltrim(std::string &var);
    void UTILITY_API  trim(std::string &var);

    std::string UTILITY_API Replace(const std::string &source,
                                    const std::string &before,
                                    const std::string &after);
    std::string UTILITY_API Plural(const std::string &noun);
    std::string UTILITY_API Plural(int, const std::string &noun);
    std::string UTILITY_API HumanReadableList(const std::vector<std::string>&);
    bool UTILITY_API IsPureASCII(const std::string &txt);
    bool UTILITY_API IsPureASCII(const char *const txt, size_t length);
    bool UTILITY_API CaseInsenstiveEqual(const std::string &str_a,
                                         const std::string &str_b);

// ****************************************************************************
//  Function: str_to_u_numeric
//
//  Purpose: Converts a string value into an unsigned numeric type as given by
//           the template parameter.
//           WARNING: This is likely to compile and silently fail if given a
//                    signed type in the template parameter.
//
//  Programmer: Tom Fogal
//  Creation:   August 11, 2008
//
//  Modifications:
//
//    Tom Fogal, Fri Aug 29 16:15:17 EDT 2008
//    Reorganized to propagate error upward.
//
//    Tom Fogal, Tue Sep 23 11:08:02 MDT 2008
//    Removed a statically-false branch which was causing an annoying warning.
//
// ****************************************************************************
    template<typename UT>
    MUST_CHECK bool str_to_u_numeric(const char * const s, UT *retval)
    {
        // strtoul() will happily convert a negative string into an unsigned
        // integer.  Do a simple check and bail out if we're given a negative
        // number.
        if(s[0] == '-')
        {
            return false;
        }

        const char *str = s;
        // get rid of leading 0's; they confuse strtoul.
        if(str[0] == '0' && str[1] != 'x')
        {
            while(*str == '0') { ++str; }
        }

        // One might want to think about switching this to an `unsigned long
        // long' and using `strtoull' below.  That will catch more cases, but
        // this is more portable.
        unsigned long ret;
        char *end;
        errno = 0;
        ret = strtoul(str, &end, 0);
        *retval = static_cast<UT>(ret);
        switch(errno)
        {
            case 0: /* success */ break;
            case ERANGE:
                // Constant does not fit in sizeof(unsigned long) bytes.
                return false;
                break;
            case EINVAL:
                // Bad base (3rd arg) given; this should be impossible.
                return false;
                break;
            default:
                // Unknown error.
                return false;
                break;
        }
        if(end == s) {
            // junk characters start the string .. is this a number?
            return false;
        }
        return true;
    }
}
#endif
