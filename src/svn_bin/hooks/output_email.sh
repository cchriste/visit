#!/bin/sh
# This is a subversion hook script intended to identify the email address of a
# user who committed a change.  It relies on a shell script existing in your
# repository which maps usernames to email addresses.  When run, this uses
# svnlook to grab your script out of the repository and store it in a
# temporary file.  Then it runs the script, giving the username of the
# committer as the sole argument.  The output of the script should be the
# committers email, and nothing else.
# The script needs the repository path and revision identifier as arguments, in
# that order.
### For example:
# We have /trunk/some_dir/emails.sh which looks something like this:
#   #!/bin/sh
#   case $1 in
#      "zaphod")
#         echo "zaphod@beeblebrox.com"
#         ;;
#      "paul")
#         echo "paul@atreides.org"
#         ;;
#      "gently")
#         echo "dirk@holisticdetectives.com"
#         ;;
#      *)
#         echo "mailinglist@yourproject.org"
#         ;;
#   esac
# So we set the below SCRIPTFILE variable to be 'trunk/some_dir/emails.sh'.
#####
#
# Script file location: location of the script relative to the repository root.
SCRIPTFILE="trunk/src/clearcase_bin/nersc_username_to_email"
# What should we use for the temporary file?
TMPFILE="/tmp/email.sh"
#### Program locations.
SVNLOOK=/usr/bin/svnlook
RM=/bin/rm
SH=/bin/sh

###########################################################
# You should not need to modify anything below this line. #
###########################################################
REPOS=$1
REV=$2

function error()
{
    echo $@ 1>&2
    exit 1
}

if test "x${REPOS}" = "x"; then
    error "Repository path not set in $0."
fi
if test "x${REV}" = "x"; then
    error "Revision number not set in $0."
fi

author=`${SVNLOOK} author -r${REV} ${REPOS}`
${SVNLOOK} cat ${REPOS} ${SCRIPTFILE} > ${TMPFILE}
if test "x$?" = "x1"; then
    error "email script does not exist in repository."
fi
address=`${SH} ${TMPFILE} $author`
${RM} -f ${TMPFILE}
/bin/echo "${address}"
