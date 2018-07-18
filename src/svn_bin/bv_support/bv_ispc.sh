function bv_ispc_initialize
{
    export DO_ISPC="no"
    export USE_SYSTEM_ISPC="no"
    export ISPC_INSTALL_DIR=""
    add_extra_commandline_args "ispc" "alt-ispc-dir" 1 "Use alternative directory for ispc"
}

function bv_ispc_enable
{
    DO_ISPC="yes"
}

function bv_ispc_disable
{
    DO_ISPC="no"
}

function bv_ispc_alt_ispc_dir
{
    echo "Using alternate ispc directory"
    bv_ispc_enable
    USE_SYSTEM_ISPC="yes"
    ISPC_INSTALL_DIR="$1"
}

function bv_ispc_depends_on
{
    echo ""
}

function bv_ispc_initialize_vars
{
    info "initializing ispc vars"
    if [[ "$DO_ISPC" == "yes" ]] ; then
        if [[ "$USE_SYSTEM_ISPC" == "no" ]]; then
            ISPC_INSTALL_DIR=$VISITDIR/ispc/$ISPC_VERSION/$VISITARCH
        fi
    fi
}

function bv_ispc_info
{
    export ISPC_VERSION=${ISPC_VERSION:-"1.9.2"}
    if [[ "$OPSYS" == "Darwin" ]] ; then
        export ISPC_FILE=${ISPC_FILE:-"ispc-v${ISPC_VERSION}-osx.tar.gz"}
        export ISPC_URL=${ISPC_URL:-"http://sdvis.org/ospray/download/dependencies/osx/"}
    else
        export ISPC_FILE=${ISPC_FILE:-"ispc-v${ISPC_VERSION}-linux.tar.gz"}
        export ISPC_URL=${ISPC_URL:-"http://sdvis.org/ospray/download/dependencies/linux/"}
    fi
    export ISPC_COMPATIBILITY_VERSION=${ISPC_COMPATIBILITY_VERSION:-"${ISPC_VERSION}"}
    export ISPC_BUILD_DIR=${ISPC_BUILD_DIR:-"${ISPC_VERSION}"}
    export ISPC_INSTALL_DIR_NAME=ispc-v$ISPC_VERSION-linux
    export ISPC_MD5_CHECKSUM=""
    export ISPC_SHA256_CHECKSUM=""
}

function bv_ispc_print
{
    printf "%s%s\n" "ISPC_FILE=" "${ISPC_FILE}"
    printf "%s%s\n" "ISPC_VERSION=" "${ISPC_VERSION}"
    printf "%s%s\n" "ISPC_COMPATIBILITY_VERSION=" "${ISPC_COMPATIBILITY_VERSION}"
    printf "%s%s\n" "ISPC_BUILD_DIR=" "${ISPC_BUILD_DIR}"
}

function bv_ispc_host_profile
{
    if [[ "$DO_ISPC" == "yes" ]]; then
        echo >> $HOSTCONF
        echo "##" >> $HOSTCONF
        echo "## ISPC" >> $HOSTCONF
        echo "##" >> $HOSTCONF
        if [[ "$USE_SYSTEM_ISPC" == "no" ]]; then
            echo "VISIT_OPTION_DEFAULT(VISIT_ISPC_ROOT \${VISITHOME}/ispc/$ISPC_VERSION/\${VISITARCH})" >> $HOSTCONF
        else
            echo "VISIT_OPTION_DEFAULT(VISIT_ISPC_ROOT ${ISPC_INSTALL_DIR})" >> $HOSTCONF
        fi
    fi
}

function bv_ispc_print_usage
{
    #ispc does not have an option, it is only dependent on ispc.
    printf "%-15s %s [%s]\n" "--ispc" "Build ISPC" "$DO_ISPC"
}

function bv_ispc_ensure
{
    if [[ "$DO_ISPC" == "yes" && "$USE_SYSTEM_ISPC" == "no" ]] ; then
        ensure_built_or_ready "ispc" $ISPC_VERSION $ISPC_BUILD_DIR $ISPC_FILE $ISPC_URL
        if [[ $? != 0 ]] ; then
            ANY_ERRORS="yes"
            DO_ISPC="no"
            error "Unable to build ISPC.  ${ISPC_FILE} not found."
        fi
    fi
}

function bv_ispc_dry_run
{
    if [[ "$DO_ISPC" == "yes" ]] ; then
        echo "Dry run option not set for ISPC."
    fi
}

# ***************************************************************************
# build_ispc
#
# Modifications:
#
# ***************************************************************************

function build_ispc
{
    # Unzip the ISPC tarball and copy it to the VisIt installation.
    info "Installing prebuilt ISPC"    
    tar zxvf $ISPC_FILE
    mkdir -p $VISITDIR/ispc/$ISPC_VERSION/$VISITARCH
    cp -R $ISPC_INSTALL_DIR_NAME/* "$VISITDIR/ispc/$ISPC_VERSION/$VISITARCH"
    rm -rf $ISPC_INSTALL_DIR_NAME
    if [[ "$DO_GROUP" == "yes" ]] ; then
        chmod -R ug+w,a+rX "$VISITDIR/ispc/$ISPC_VERSION/$VISITARCH"
        chgrp -R ${GROUP} "$VISITDIR/ispc/$ISPC_VERSION/$VISITARCH"
    fi
    cd "$START_DIR"
    info "Done with ISPC"
    return 0
}

function bv_ispc_is_enabled
{
    if [[ $DO_ISPC == "yes" ]]; then
        return 1    
    fi
    return 0
}

function bv_ispc_is_installed
{
    if [[ "$USE_SYSTEM_ISPC" == "yes" ]]; then   
        return 1
    fi

    check_if_installed "ispc" $ISPC_VERSION
    if [[ $? == 0 ]] ; then
        return 1
    fi
    return 0
}

function bv_ispc_build
{
    if [[ "$DO_ISPC" == "yes" && "$USE_SYSTEM_ISPC" == "no" ]] ; then
        check_if_installed "ispc" $ISPC_VERSION
        if [[ $? == 0 ]] ; then
            info "Skipping build of ISPC"
        else
            build_ispc
            if [[ $? != 0 ]] ; then
                error "Unable to build or install ISPC.  Bailing out."
            fi
            info "Done building ISPC"
        fi
    fi
}

