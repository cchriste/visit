/*****************************************************************************
*
* Copyright (c) 2000 - 2007, The Regents of the University of California
* Produced at the Lawrence Livermore National Laboratory
* All rights reserved.
*
* This file is part of VisIt. For details, see http://www.llnl.gov/visit/. The
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
*    documentation and/or materials provided with the distribution.
*  - Neither the name of the UC/LLNL nor  the names of its contributors may be
*    used to  endorse or  promote products derived from  this software without
*    specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT  HOLDERS AND CONTRIBUTORS "AS IS"
* AND ANY EXPRESS OR  IMPLIED WARRANTIES, INCLUDING,  BUT NOT  LIMITED TO, THE
* IMPLIED WARRANTIES OF MERCHANTABILITY AND  FITNESS FOR A PARTICULAR  PURPOSE
* ARE  DISCLAIMED.  IN  NO  EVENT  SHALL  THE  REGENTS  OF  THE  UNIVERSITY OF
* CALIFORNIA, THE U.S.  DEPARTMENT  OF  ENERGY OR CONTRIBUTORS BE  LIABLE  FOR
* ANY  DIRECT,  INDIRECT,  INCIDENTAL,  SPECIAL,  EXEMPLARY,  OR CONSEQUENTIAL
* DAMAGES (INCLUDING, BUT NOT  LIMITED TO, PROCUREMENT OF  SUBSTITUTE GOODS OR
* SERVICES; LOSS OF  USE, DATA, OR PROFITS; OR  BUSINESS INTERRUPTION) HOWEVER
* CAUSED  AND  ON  ANY  THEORY  OF  LIABILITY,  WHETHER  IN  CONTRACT,  STRICT
* LIABILITY, OR TORT  (INCLUDING NEGLIGENCE OR OTHERWISE)  ARISING IN ANY  WAY
* OUT OF THE  USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
* DAMAGE.
*
*****************************************************************************/

#ifndef ATOMIC_PROPERTIES_H
#define ATOMIC_PROPERTIES_H
#include <utility_exports.h>

#define MAX_ELEMENT_NUMBER 109
#define KNOWN_AMINO_ACIDS  23

extern UTILITY_API unsigned char aminocolors[KNOWN_AMINO_ACIDS][3];
extern UTILITY_API unsigned char shapelycolors[KNOWN_AMINO_ACIDS][3];
extern UTILITY_API unsigned char jmolcolors[MAX_ELEMENT_NUMBER][3];
extern UTILITY_API unsigned char rasmolcolors[MAX_ELEMENT_NUMBER][3];
extern UTILITY_API float         atomic_radius[MAX_ELEMENT_NUMBER];
extern UTILITY_API float         covalent_radius[MAX_ELEMENT_NUMBER];

extern UTILITY_API const char   *element_names[MAX_ELEMENT_NUMBER];
extern UTILITY_API const char   *residue_names[KNOWN_AMINO_ACIDS];

int  UTILITY_API ElementNameToAtomicNumber(const char *element);
int  UTILITY_API ResiduenameToNumber(const char *name);
UTILITY_API const char *NumberToResiduename(int num);
int UTILITY_API NumberOfKnownResidues();
int  UTILITY_API ResiduenameToAbbr(const char *name);
void UTILITY_API AddResiduename(const char *name, const char *longname);
UTILITY_API const char *ResiduenameToLongName(const char *name);
int UTILITY_API ResidueLongnameMaxlen();

#endif
