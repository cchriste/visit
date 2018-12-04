#!/bin/sh
#
# Purpose: Maintain common, global variables and functions
#  used by hook scripts
#
# Programmer: Mark C. Miller
#
# Modifications:
#
#   Mark C. Miller, Tue Dec  9 23:05:19 PST 2008
#   Replaced explicit enumeration of various needed commands as
#   exported variables with just setting PATH variable.
#   Added log() and HandleCommonSkipCases() functions and
#   exported them to subshells.
#
#   Tom Fogal, Sat May 16 14:55:06 MDT 2009
#   Don't do any checks for third_party files.


export PATH=/bin:/usr/bin
export VISIT_GROUP_NAME=visitdev
export hadError=0

function log()
{
    hadError=1
    echo "$@" 1>&2
}
export -f log

#
# Determine if a given file can be skipped by hooks
# based on its status as given by svnlook, whether
# it contains a trailing slash ('/') character
# indicating its a directory and whether it is a text
# file or not.
#
function HandleCommonSkipCases()
{
    # $1 is file status string from svnlook
    # $2 is file name
    # return of 0 means TRUE (e.g. yes, skip it)
    # return of 1 means FALSE (e.g. DO NOT SKIP)

    #
    # Skip deletions (D) or just changes in props (_U)
    #
    if test $1 = D || test $1 = _U; then
        return 0
    fi

    #
    # third_party is a free for all.
    #
    echo "$2" | /bin/grep "third_party" >/dev/null
    if test $? -eq 0; then
      return 0
    fi

    #
    # Skip entries which are the names of directories.
    # This logic checks if the last character of the entry
    # is a slash.
    #
    if test ${2:${#2}-1:1} = /; then
        return 0
    fi

    #
    # Only do this check for files svn thinks are 'text' files
    #
    hasMimeTypeProp=`svnlook proplist -t $TXN $REPOS $2 | grep mime-type`
    if test -n "$hasMimeTypeProp"; then
        mimeTypeProp=`svnlook propget -t $TXN $REPOS svn:mime-type $2`
        if test -n "$mimeTypeProp"; then
            if test -z "`echo $mimeTypeProp | grep ^text/`"; then
                return 0
            fi
        fi
    fi

    return 1
}
export -f HandleCommonSkipCases
