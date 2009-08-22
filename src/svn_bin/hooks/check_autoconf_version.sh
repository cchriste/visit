#!/bin/sh
##############################################################################
#
# Purpose: Check that configure.in/configure always get updated together.
#          Also, check that correct version of autoconf was used to
#          generate configure.
#
# Programmer: Mark C. Miller
# Created:    Mon Apr  7 18:16:51 PDT 2008
#
# Modifications:
#
#   Tom Fogal, Tue Jun 24 11:55:51 EDT 2008
#   Added a bypass/allow for modifying configure if you've modified any m4
#   files.  Also, 's/\t/\ \ \ \ /g'
#
#   Mark C. Miller, Tue Nov 18 18:22:40 PST 2008
#   Fixed typo in logic involving || for $m4_file
#
#   Mark C. Miller, Mon Dec  8 12:51:21 PST 2008
#   Fixed typo in logic involving && for $m4_file
#
#   Mark C. Miller, Tue Dec  9 00:19:04 PST 2008
#   Obtain list of changed files via FLIST ($3) argument and loop
#   over them via 'read' sh builtin method.
#
#   Mark C. Miller, Tue Dec  9 23:08:23 PST 2008
#   Adjusted file loop to deal with fact that FLIST file now contains
#   both status chars and filename. Replaced ${<VARNAME>} command refs
#   with just the commands themselves.
#
#   Mark C. Miller, Tue Mar 24 22:02:44 PDT 2009
#   Changed logic to simply record the configure file in one variable
#   and all other files that serve as input to configure in another
#   variable.
#
#   Tom Fogal, Sat Aug 22 16:21:30 MDT 2009
#   Ignore the check on any of my branches.
#
##############################################################################
REPOS="$1"
TXN="$2"
FLIST="$3"

theConfigureFile=""
configureInputFiles=""
while read fline; do

    #
    # Get file 'svnlook' status and name
    #
    fstat=`echo $fline | tr -s ' ' | cut -d' ' -f1`
    fname=`echo $fline | tr -s ' ' | cut -d' ' -f2`

    #
    # We're looking for only these specific files
    # for this hook.
    case $fname in
        *fogal1*)
            # nothing.
            ;;
        *src/configure)
            theConfigureFile=$fname
            ;;
        *src/configure.in|*src/*m4|*src/ac/*)
            configureInputFiles="$configureInputFiles $fname"
            ;;
    esac

done < $FLIST

# Check that if any files that serve as input to configure are getting
# updated, then configure itself is also getting updated. 
if test -n "$configureInputFiles" && test -z "$theConfigureFile"; then
    log "Attempt to commit changes to $configureInputFiles"
    log "without also committing changes to \"configure\""
    exit 1
fi

# Check that if configure itself is getting updated, then at least
# one of the files it depends on is also getting updated.
if test -z "$configureInputFiles" && test -n "$theConfigureFile"; then
    log "Attempt to commit changes to $theConfigureFile"
    log "without also committing changes to any of its input files."
    exit 1
fi

# get and check autoconf version number
if test -n "$theConfigureFile"; then
    acVno=`svnlook cat -t $TXN $REPOS $theConfigureFile | grep 'Generated by GNU Autoconf' | tr -s ' ' | cut -d' ' -f6 | cut -d'.' -f1,2`
    if test "$acVno" != "$AUTOCONF_VERSION_NUMBER"; then
        log "You must use autoconf version $AUTOCONF_VERSION_NUMBER to re-generate configure."
        log "You have used $acVno"
        exit 1
    fi
fi

# all is well!
exit 0
