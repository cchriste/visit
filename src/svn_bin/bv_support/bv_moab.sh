function bv_moab_initialize
{
    export DO_MOAB="no"
    export ON_MOAB="off"
}

function bv_moab_enable
{
    DO_MOAB="yes"
    ON_MOAB="on"
}

function bv_moab_disable
{
    DO_MOAB="no"
    ON_MOAB="off"
}

function bv_moab_depends_on
{
    local depends_on="hdf5"

    if [[ "$DO_SZIP" == "yes" ]] ; then
        depends_on="$depends_on szip"
    fi

    if [[ "$DO_ZLIB" == "yes" ]] ; then
        depends_on="$depends_on zlib"
    fi

    echo $depends_on
}

function bv_moab_info
{
    export MOAB_VERSION=${MOAB_VERSION:-"4.9.2-RC0"}
    export MOAB_FILE=${MOAB_FILE:-"moab-${MOAB_VERSION}.tar.gz"}
    export MOAB_BUILD_DIR=${MOAB_BUILD_DIR:-"moab-4.9.2"}
    export MOAB_URL=${MOAB_URL:-"ftp://ftp.mcs.anl.gov/pub/fathom"}
}

function bv_moab_print
{
    printf "%s%s\n" "MOAB_FILE=" "${MOAB_FILE}"
    printf "%s%s\n" "MOAB_VERSION=" "${MOAB_VERSION}"
    printf "%s%s\n" "MOAB_BUILD_DIR=" "${MOAB_BUILD_DIR}"
}

function bv_moab_print_usage
{
    printf "%-15s %s [%s]\n" "--moab" "Build moab support" "$DO_MOAB"
}

function bv_moab_graphical
{
    local graphical_out="moab     $MOAB_VERSION($MOAB_FILE)      $ON_MOAB"
    echo "$graphical_out"
}

function bv_moab_host_profile
{
    if [[ "$DO_MOAB" == "yes" ]] ; then
        echo >> $HOSTCONF
        echo "##" >> $HOSTCONF
        echo "## MOAB " >> $HOSTCONF
        echo "##" >> $HOSTCONF
        echo \
            "VISIT_OPTION_DEFAULT(VISIT_MOAB_DIR \${VISITHOME}/moab/$MOAB_VERSION/\${VISITARCH})" \
            >> $HOSTCONF
        echo \
            "VISIT_OPTION_DEFAULT(VISIT_MOAB_LIBDEP HDF5_LIBRARY_DIR hdf5 \${VISIT_HDF5_LIBDEP} TYPE STRING)" \
            >> $HOSTCONF
        if [[ -n "$PAR_COMPILER" ]]; then
            echo \
                "VISIT_OPTION_DEFAULT(VISIT_MOAB_MIPPAR_DIR \${VISITHOME}/moab/$MOAB_VERSION/mpipar/\${VISITARCH})" \
                >> $HOSTCONF
            echo \
                "VISIT_OPTION_DEFAULT(VISIT_MOAB_MPIPAR_LIBDEP HDF5_MPIPAR_LIBRARY_DIR hdf5 \${VISIT_HDF5_MPIPAR_LIBDEP} TYPE STRING)" \
                >> $HOSTCONF
        fi
    fi
}

function bv_moab_ensure
{
    if [[ "$DO_MOAB" == "yes" ]] ; then
        ensure_built_or_ready "moab" $MOAB_VERSION $MOAB_BUILD_DIR $MOAB_FILE $MOAB_URL
        if [[ $? != 0 ]] ; then
            ANY_ERRORS="yes"
            DO_MOAB="no"
            error "Unable to build moab.  ${MOAB_FILE} not found."
        fi
    fi
}

function bv_moab_dry_run
{
    if [[ "$DO_MOAB" == "yes" ]] ; then
        echo "Dry run option not set for moab."
    fi
}

# *************************************************************************** #
#                            Function 8, build_moab
# *************************************************************************** #
function build_moab
{
    #
    # Prepare build dir
    #
    prepare_build_dir $MOAB_BUILD_DIR $MOAB_FILE
    untarred_moab=$?
    if [[ $untarred_moab == -1 ]] ; then
        warn "Unable to prepare moab build directory. Giving Up!"
        return 1
    fi

    cd $MOAB_BUILD_DIR || error "Can't cd to moab build dir."
    rm -f src/moab/MOABConfig.h # work around a potential issue in MOAB tarball

    par_build_types="serial"
    if [[ -n "$PAR_COMPILER_CXX" ]]; then
        par_build_types="$par_build_types parallel"
    fi

    for bt in $par_build_types; do 

        mkdir build_$bt
        pushd build_$bt

        cf_mpi_arg=""
        cf_par_prefix=""
        if [[ "$bt" == "serial" ]]; then
            cf_c_compiler="$C_COMPILER"
            cf_cxx_compiler="$CXX_COMPILER"
        elif [[ "$bt" == "parallel" ]]; then
            cf_mpi_arg="--with-mpi"
            cf_par_prefix="mpipar/"
            cf_c_compiler="$PAR_COMPILER"
            cf_cxx_compiler="$PAR_COMPILER_CXX"
        fi

        cf_prefix_arg="--prefix=$VISITDIR/moab/$MOAB_VERSION/${par_prefix}$VISITARCH"
        cf_common_args="--with-pic --disable-fortran"

        if [[ "DO_STATIC_BUILD" == "yes" ]]; then
            cf_static_args="--enable-static --disable-shared"
        else
            cf_static_args="--disable-static --enable-shared"
        fi

        cf_hdf5_ldflags_arg=""
        cf_szip_arg=""
        cf_zlib_arg=""
        cf_hdf5_arg="--with-hdf5=$VISITDIR/hdf5/$HDF5_VERSION/${par_prefix}$VISITARCH"
        if [[ "$DO_SZIP" == "yes" ]] ; then
            cf_szip_arg="--with-szip=$VISITDIR/szip/$SZIP_VERSION/$VISITARCH"
            cf_hdf5_ldflags_arg="-lsz"
        fi
        if [[ "$DO_ZLIB" == "yes" ]] ; then
            cf_zlib_arg="--with-zlib=$VISITDIR/zlib/$ZLIB_VERSION/$VISITARCH"
            cf_hdf5_ldflags_arg="$cf_hdf5_ldflags_arg -lz"
        fi
        if [[ -n "$cf_hdf5_ldflags_arg" ]]; then
            cf_hdf5_ldflags_arg="--with-hdf5-ldflags=\"$cf_hdf5_ldflags_arg\""
        fi

        info "Configuring $bt moab . . ."
        info ../configure CXX=\"$cf_cxx_compiler\" CXXFLAGS=\"$CXXFLAGS $CXX_OPT_FLAGS\" \
            CC=\"$cf_c_compiler\" CFLAGS=\"$CFLAGS $C_OPT_FLAGS\" \
            ${cf_prefix_arg} ${cf_mpi_arg} ${cf_common_args} ${cf_static_args} \
            ${cf_hdf5_arg} ${cf_hdf5_ldflags_arg} \
            ${cf_szip_arg} ${cf_zlib_arg}

        sh -c "../configure \
            CXX=\"$cf_cxx_compiler\" CXXFLAGS=\"$CXXFLAGS $CXX_OPT_FLAGS\" \
            CC=\"$cf_c_compiler\" CFLAGS=\"$CFLAGS $C_OPT_FLAGS\" \
            ${cf_prefix_arg} ${cf_mpi_arg} ${cf_common_args} ${cf_static_args} \
            ${cf_hdf5_arg} ${cf_hdf5_ldflags_arg} \
            ${cf_szip_arg} ${cf_zlib_arg}"

        if [[ $? != 0 ]] ; then
            warn "$bt MOAB configure failed.  Giving up"
            return 1
        fi

        #
        # Build moab
        #

        info "Building $bt moab . . . (~2 minutes)"
        $MAKE $MAKE_OPT_FLAGS
        if [[ $? != 0 ]] ; then
            warn "$bt moab build failed.  Giving up"
            return 1
        fi

        #
        # Install into the VisIt third party location.
        #
        info "Installing $bt moab"
        $MAKE install

        if [[ "$DO_GROUP" == "yes" ]] ; then
            chmod -R ug+w,a+rX "$VISITDIR/moab"
            chgrp -R ${GROUP} "$VISITDIR/moab"
        fi

        popd
    done

    cd "$START_DIR"
    info "Done with moab"
    return 0
}


function bv_moab_is_enabled
{
    if [[ $DO_MOAB == "yes" ]]; then
        return 1    
    fi
    return 0
}

function bv_moab_is_installed
{
    check_if_installed "moab" $MOAB_VERSION
    if [[ $? == 0 ]] ; then
        return 1
    fi
    return 0
}

function bv_moab_build
{
    cd "$START_DIR"
    if [[ "$DO_MOAB" == "yes" ]] ; then
        check_if_installed "moab" $MOAB_VERSION
        if [[ $? == 0 ]] ; then
            info "Skipping moab build.  moab is already installed."
        else
            info "Building moab (~2 minutes)"
            build_moab
            if [[ $? != 0 ]] ; then
                error "Unable to build or install moab.  Bailing out."
            fi
            info "Done building moab"
        fi
    fi
}
