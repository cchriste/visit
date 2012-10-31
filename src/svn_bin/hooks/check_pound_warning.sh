#!/bin/sh
##############################################################################
#
# Purpose: Check for existence of #warning compiler directives. 
#
# Programmer: Mark C. Miller
# Created:    May 20, 2008 
#
# Modifications:
#   Mark C. Miller, Tue Dec  9 00:19:04 PST 2008
#   Obtain list of changed files via FLIST ($3) argument and loop
#   over them via 'read' sh builtin method.
#
#   Mark C. Miller, Tue Dec  9 23:11:58 PST 2008
#   Re-factored most of skipping logic to CanHandleCommonSkipCases.
#   Adjusted main file loop to accomodate fact that FLIST file now contains
#   file status chars as well as file names.
#
#   Mark C. Miller, Mon Feb  7 20:33:45 PST 2011
#   Permit #warning in branches.
#
#   Kathleen Biagas, Tue Nov 29 22:30:09 PST 2011
#   Skip windowsbuild directory.
#
##############################################################################
REPOS="$1"
TXN="$2"
FLIST="$3"

while read fline; do

    #
    # Get file 'svnlook' status and name
    #
    fstat=`echo $fline | tr -s ' ' | cut -d' ' -f1`
    fname=`echo $fline | tr -s ' ' | cut -d' ' -f2`

    #
    # Skip common cases of deletions, dirs, non-text files
    #
    if `HandleCommonSkipCases $fstat $fname`; then
        continue
    fi

    case $fname in
        trunk/windowsbuild/*)
            # Skip anything in windowsbuild
            continue
            ;;
        trunk/*)
            # Do not skip (continue) these
            ;;
        *)
            # Skip anything not on trunk
            continue
            ;;
    esac

    svnlook cat -t $TXN $REPOS $fname | grep -q '^#warning' 1>/dev/null 2>&1
    commitFileHasPoundWarnings=$?

    # If the file we're committing has #warnings, reject it
    if test $commitFileHasPoundWarnings -eq 0; then
        log "File \"$fname\" appears to contain '#warning' compilation directives."
        log "Please remove them before committing."
        exit 1
    fi

done < $FLIST

# all is well!
exit 0
