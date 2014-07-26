#ifndef MOVIE_SEQUENCE_TRANSITION_H
#define MOVIE_SEQUENCE_TRANSITION_H
/*****************************************************************************
*
* Copyright (c) 2000 - 2014, Lawrence Livermore National Security, LLC
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
#include <MovieSequence.h>

// ****************************************************************************
// Class: MovieSequenceTransition
//
// Purpose:
//   Base class for movie transitions.
//
// Notes:      
//
// Programmer: Brad Whitlock
// Creation:   Thu Oct 12 14:25:12 PST 2006
//
// Modifications:
//   Brad Whitlock, Tue Apr  8 09:21:55 PDT 2008
//   Moved SequenceMenuName methods to the .C file.
//
// ****************************************************************************

class MovieSequenceTransition : public MovieSequence
{
public:
    MovieSequenceTransition();
    virtual ~MovieSequenceTransition();

    virtual int SubMenuIndex() const { return 0; }

    virtual QWidget *CreateUI();
    virtual void ReadUIValues(QWidget *ui, DataNode *node);
    virtual void WriteUIValues(QWidget *ui, DataNode *node);
};

class MovieSequenceFade : public MovieSequenceTransition
{
public:
    MovieSequenceFade() : MovieSequenceTransition() { }
    virtual ~MovieSequenceFade()                    { }

    virtual MovieSequence *NewInstance() const      { return new MovieSequenceFade; }
    virtual int SequenceId() const                  { return 100; }
    virtual std::string SequenceName() const        { return "Fade"; }
    virtual QString SequenceMenuName() const;
    virtual QPixmap Pixmap();
};

//
// Types of wipe screen transitions.
//

class MovieSequenceLRWipe : public MovieSequenceTransition
{
public:
    MovieSequenceLRWipe() : MovieSequenceTransition() { }
    virtual ~MovieSequenceLRWipe()                    { }

    virtual MovieSequence *NewInstance() const        { return new MovieSequenceLRWipe; }
    virtual int SequenceId() const                    { return 101; }
    virtual std::string SequenceName() const          { return "LRWipe"; }
    virtual QString SequenceMenuName() const;
    virtual QPixmap Pixmap();
};

class MovieSequenceRLWipe : public MovieSequenceTransition
{
public:
    MovieSequenceRLWipe() : MovieSequenceTransition() { }
    virtual ~MovieSequenceRLWipe()                    { }

    virtual MovieSequence *NewInstance() const        { return new MovieSequenceRLWipe; }
    virtual int SequenceId() const                    { return 102; }
    virtual std::string SequenceName() const          { return "RLWipe"; }
    virtual QString SequenceMenuName() const;
    virtual QPixmap Pixmap();
};

class MovieSequenceTBWipe : public MovieSequenceTransition
{
public:
    MovieSequenceTBWipe() : MovieSequenceTransition() { }
    virtual ~MovieSequenceTBWipe()                    { }

    virtual MovieSequence *NewInstance() const        { return new MovieSequenceTBWipe; }
    virtual int SequenceId() const                    { return 103; }
    virtual std::string SequenceName() const          { return "TBWipe"; }
    virtual QString SequenceMenuName() const;
    virtual QPixmap Pixmap();
};

class MovieSequenceBTWipe : public MovieSequenceTransition
{
public:
    MovieSequenceBTWipe() : MovieSequenceTransition() { }
    virtual ~MovieSequenceBTWipe()                    { }

    virtual MovieSequence *NewInstance() const        { return new MovieSequenceBTWipe; }
    virtual int SequenceId() const                    { return 104; }
    virtual std::string SequenceName() const          { return "BTWipe"; }
    virtual QString SequenceMenuName() const;
    virtual QPixmap Pixmap();
};

//
// Types of slide screen transitions.
//

class MovieSequenceLRSlide : public MovieSequenceTransition
{
public:
    MovieSequenceLRSlide() : MovieSequenceTransition() { }
    virtual ~MovieSequenceLRSlide()                    { }

    virtual MovieSequence *NewInstance() const        { return new MovieSequenceLRSlide; }
    virtual int SequenceId() const                    { return 105; }
    virtual std::string SequenceName() const          { return "LRSlide"; }
    virtual QString SequenceMenuName() const;
    virtual QPixmap Pixmap();
};

class MovieSequenceRLSlide : public MovieSequenceTransition
{
public:
    MovieSequenceRLSlide() : MovieSequenceTransition() { }
    virtual ~MovieSequenceRLSlide()                    { }

    virtual MovieSequence *NewInstance() const        { return new MovieSequenceRLSlide; }
    virtual int SequenceId() const                    { return 106; }
    virtual std::string SequenceName() const          { return "RLSlide"; }
    virtual QString SequenceMenuName() const;
    virtual QPixmap Pixmap();
};

class MovieSequenceTBSlide : public MovieSequenceTransition
{
public:
    MovieSequenceTBSlide() : MovieSequenceTransition() { }
    virtual ~MovieSequenceTBSlide()                    { }

    virtual MovieSequence *NewInstance() const        { return new MovieSequenceTBSlide; }
    virtual int SequenceId() const                    { return 107; }
    virtual std::string SequenceName() const          { return "TBSlide"; }
    virtual QString SequenceMenuName() const;
    virtual QPixmap Pixmap();
};

class MovieSequenceBTSlide : public MovieSequenceTransition
{
public:
    MovieSequenceBTSlide() : MovieSequenceTransition() { }
    virtual ~MovieSequenceBTSlide()                    { }

    virtual MovieSequence *NewInstance() const        { return new MovieSequenceBTSlide; }
    virtual int SequenceId() const                    { return 108; }
    virtual std::string SequenceName() const          { return "BTSlide"; }
    virtual QString SequenceMenuName() const;
    virtual QPixmap Pixmap();
};

//
// Other transitions
//

class MovieSequenceCircle : public MovieSequenceTransition
{
public:
    MovieSequenceCircle() : MovieSequenceTransition() { }
    virtual ~MovieSequenceCircle()                    { }

    virtual MovieSequence *NewInstance() const        { return new MovieSequenceCircle; }
    virtual int SequenceId() const                    { return 109; }
    virtual std::string SequenceName() const          { return "CircularReveal"; }
    virtual QString SequenceMenuName() const;
    virtual QPixmap Pixmap();
};


class MovieSequenceRipple: public MovieSequenceTransition
{
public:
    MovieSequenceRipple() : MovieSequenceTransition() { }
    virtual ~MovieSequenceRipple()                    { }

    virtual MovieSequence *NewInstance() const        { return new MovieSequenceRipple; }
    virtual int SequenceId() const                    { return 110; }
    virtual std::string SequenceName() const          { return "RippleFade"; }
    virtual QString SequenceMenuName() const;
    virtual QPixmap Pixmap();
};


class MovieSequenceBlock: public MovieSequenceTransition
{
public:
    MovieSequenceBlock() : MovieSequenceTransition() { }
    virtual ~MovieSequenceBlock()                    { }

    virtual MovieSequence *NewInstance() const        { return new MovieSequenceBlock; }
    virtual int SequenceId() const                    { return 111; }
    virtual std::string SequenceName() const          { return "BlockFade"; }
    virtual QString SequenceMenuName() const;
    virtual QPixmap Pixmap();
};


#endif

