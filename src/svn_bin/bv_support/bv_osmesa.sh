function bv_osmesa_initialize
{
    export DO_OSMESA="no"
}

function bv_osmesa_enable
{
    DO_OSMESA="yes"
}

function bv_osmesa_disable
{
    DO_OSMESA="no"
}

function bv_osmesa_depends_on
{
    depends_on="llvm"

    echo ${depends_on}
}

function bv_osmesa_info
{
    export OSMESA_VERSION=${OSMESA_VERSION:-"17.1.9"}
    export OSMESA_FILE=${OSMESA_FILE:-"mesa-$OSMESA_VERSION.tar.gz"}
    export OSMESA_BUILD_DIR=${OSMESA_BUILD_DIR:-"mesa-$OSMESA_VERSION"}
    export OSMESA_URL=${OSMESA_URL:-"ftp://ftp.freedesktop.org/pub/mesa"}
    export OSMESA_MD5_CHECKSUM=""
    export OSMESA_SHA256_CHECKSUM=""
}

function bv_osmesa_print
{
    printf "%s%s\n" "OSMESA_FILE=" "${OSMESA_FILE}"
    printf "%s%s\n" "OSMESA_VERSION=" "${OSMESA_VERSION}"
    printf "%s%s\n" "OSMESA_TARGET=" "${OSMESA_TARGET}"
    printf "%s%s\n" "OSMESA_BUILD_DIR=" "${OSMESA_BUILD_DIR}"
}

function bv_osmesa_print_usage
{
    printf "%-15s %s [%s]\n" "--osmesa" "Build OSMesa" "$DO_OSMESA"
}

function bv_osmesa_host_profile
{
    # If we are using osmesa as the GL for VTK in a static build, we'll tell
    # VisIt about osmesa using a different mechanism.
    addhp="yes"
    if [[ "$DO_STATIC_BUILD" == "yes" ]] ; then
        if [[ "$DO_SERVER_COMPONENTS_ONLY" == "yes" || "$DO_ENGINE_ONLY" == "yes" ]] ; then
            addhp="no"
        fi
    fi

    if [[ "$DO_OSMESA" == "yes" && "$addhp" == "yes" ]] ; then
        echo >> $HOSTCONF
        echo "##" >> $HOSTCONF
        echo "## OSMesa" >> $HOSTCONF
        echo "##" >> $HOSTCONF
        echo "VISIT_OPTION_DEFAULT(VISIT_OSMESA_DIR \${VISITHOME}/osmesa/$OSMESA_VERSION/\${VISITARCH})" >> $HOSTCONF
    fi
}

function bv_osmesa_selected
{
    args=$@
    if [[ $args == "--osmesa" ]]; then
        DO_OSMESA="yes"
        return 1
    fi

    return 0
}

function bv_osmesa_initialize_vars
{
    info "initalizing osmesa vars"
    if [[ "$DO_OSMESA" == "yes" ]]; then
        OSMESA_INSTALL_DIR="${VISITDIR}/osmesa/${OSMESA_VERSION}/${VISITARCH}"
        OSMESA_INCLUDE_DIR="${OSMESA_INSTALL_DIR}/include"
        OSMESA_LIB_DIR="${OSMESA_INSTALL_DIR}/lib"
        if [[ "$DO_STATIC_BUILD" == "yes" ]]; then
            OSMESA_LIB="${OSMESA_LIB_DIR}/libOSMesa.a"
        else
            OSMESA_LIB="${OSMESA_LIB_DIR}/libOSMesa.${SO_EXT}"
        fi
    fi
}

function bv_osmesa_ensure
{
    if [[ "$DO_DBIO_ONLY" != "yes" ]]; then
        if [[ "$DO_OSMESA" == "yes" ]] ; then
            ensure_built_or_ready "osmesa"   $OSMESA_VERSION   $OSMESA_BUILD_DIR   $OSMESA_FILE $OSMESA_URL
            if [[ $? != 0 ]] ; then
                return 1
            fi
        fi
    fi
}

function bv_osmesa_dry_run
{
    if [[ "$DO_OSMESA" == "yes" ]] ; then
        echo "Dry run option not set for osmesa."
    fi
}

function build_osmesa
{
    #
    # prepare build dir
    #
    prepare_build_dir $OSMESA_BUILD_DIR $OSMESA_FILE
    untarred_osmesa=$?
    if [[ $untarred_osmesa == -1 ]] ; then
        warn "Unable to prepare Mesa build directory. Giving Up!"
        return 1
    fi

    #
    # Build OSMESA.
    #
    cd $OSMESA_BUILD_DIR || error "Couldn't cd to osmesa build dir."

    if [[ "$DO_STATIC_BUILD" == "yes" ]]; then
        OSMESA_STATIC_DYNAMIC="--disable-shared --disable-shared-glapi --enable-static --enable-static-glapi"
    fi

    info "Configuring OSMesa . . ."
    echo CXXFLAGS="${CXXFLAGS} ${CXX_OPT_FLAGS}" \
        CXX=${CXX_COMPILER} \
        CFLAGS="${CFLAGS} ${C_OPT_FLAGS}" \
        CC=${C_COMPILER} \
        ./autogen.sh \
        --prefix=${VISITDIR}/osmesa/${OSMESA_VERSION}/${VISITARCH} \
        --disable-gles1 \
        --disable-gles2 \
        --disable-dri \
        --disable-dri3 \
        --disable-glx \
        --disable-glx-tls \
        --disable-egl \
        --disable-gbm \
        --disable-xvmc \
        --disable-vdpau \
        --disable-omx \
        --disable-va \
        --with-platforms= \
        --with-gallium-drivers=swrast,swr \
        --enable-gallium-osmesa $OSMESA_STATIC_DYNAMIC \
        --with-llvm-prefix=${VISIT_LLVM_DIR}
    env CXXFLAGS="${CXXFLAGS} ${CXX_OPT_FLAGS}" \
        CXX=${CXX_COMPILER} \
        CFLAGS="${CFLAGS} ${C_OPT_FLAGS}" \
        CC=${C_COMPILER} \
        ./autogen.sh \
        --prefix=${VISITDIR}/osmesa/${OSMESA_VERSION}/${VISITARCH} \
        --disable-gles1 \
        --disable-gles2 \
        --disable-dri \
        --disable-dri3 \
        --disable-glx \
        --disable-glx-tls \
        --disable-egl \
        --disable-gbm \
        --disable-xvmc \
        --disable-vdpau \
        --disable-omx \
        --disable-va \
        --with-platforms= \
        --with-gallium-drivers=swrast,swr \
        --enable-gallium-osmesa $OSMESA_STATIC_DYNAMIC \
        --with-llvm-prefix=${VISIT_LLVM_DIR}

    if [[ $? != 0 ]] ; then
        warn "OSMesa configure failed.  Giving up"
        return 1
    fi

    info "Building OSMesa . . ."
    ${MAKE} ${MAKE_OPT_FLAGS}
    if [[ $? != 0 ]] ; then
        warn "OSMesa build failed.  Giving up"
        return 1
    fi

    info "Installing OSMesa ..."
    ${MAKE} ${MAKE_OPT_FLAGS} install
    if [[ $? != 0 ]] ; then
        warn "OSMesa install failed.  Giving up"
        return 1
    fi

    if [[ "$DO_GROUP" == "yes" ]] ; then
        chmod -R ug+w,a+rX "$VISITDIR/osmesa"
        chgrp -R ${GROUP} "$VISITDIR/osmesa"
    fi
    cd "$START_DIR"
    info "Done with OSMesa"
    return 0
}

function bv_osmesa_is_enabled
{
    if [[ $DO_OSMESA == "yes" ]]; then
        return 1    
    fi
    return 0
}

function bv_osmesa_is_installed
{
    check_if_installed "osmesa" $OSMESA_VERSION
    if [[ $? == 0 ]] ; then
        return 1
    fi
    return 0
}

function bv_osmesa_build
{
    #
    # Build OSMesa
    #
    cd "$START_DIR"
    if [[ "$DO_OSMESA" == "yes" ]] ; then
        check_if_installed "osmesa" $OSMESA_VERSION
        if [[ $? == 0 ]] ; then
            info "Skipping OSMesa build.  OSMesa is already installed."
        else
            info "Building OSMesa (~20 minutes)"
            build_osmesa
            if [[ $? != 0 ]] ; then
                error "Unable to build or install OSMesa.  Bailing out."
            fi
            info "Done building OSMesa"
        fi
    fi
}
