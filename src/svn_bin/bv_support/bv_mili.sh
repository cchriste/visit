function bv_mili_initialize
{
export DO_MILI="no"
export ON_MILI="off"
}

function bv_mili_enable
{
DO_MILI="yes"
ON_MILI="on"
}

function bv_mili_disable
{
DO_MILI="no"
ON_MILI="off"
}

function bv_mili_depends_on
{
return ""
}

function bv_mili_info
{
export MILI_FILE=${MILI_FILE:-"mili-1.10a.tar.gz"}
export MILI_VERSION=${MILI_VERSION:-"1.10.0"}
export MILI_COMPATIBILITY_VERSION=${MILI_COMPATIBILITY_VERSION:-"1.0"}
export MILI_BUILD_DIR=${MILI_BUILD_DIR:-"mili"}
}

function bv_mili_print
{
  printf "%s%s\n" "MILI_FILE=" "${MILI_FILE}"
  printf "%s%s\n" "MILI_VERSION=" "${MILI_VERSION}"
  printf "%s%s\n" "MILI_COMPATIBILITY_VERSION=" "${MILI_COMPATIBILITY_VERSION}"
  printf "%s%s\n" "MILI_BUILD_DIR=" "${MILI_BUILD_DIR}"
}

function bv_mili_print_usage
{
printf "\t\t%15s\n" "NOTE: not available for download from web"
printf "%-15s %s [%s]\n" "--mili" "Build Mili" "$DO_MILI"
}

function bv_mili_graphical
{
local graphical_out="MILI    $MILI_VERSION($MILI_FILE)     $ON_MILI"
echo "$graphical_out"
}

function bv_mili_host_profile
{
    if [[ "$DO_MILI" == "yes" ]] ; then
        echo >> $HOSTCONF
        echo "##" >> $HOSTCONF
        echo "## Mili" >> $HOSTCONF
        echo "##" >> $HOSTCONF
        echo \
        "VISIT_OPTION_DEFAULT(VISIT_MILI_DIR \${VISITHOME}/mili/$MILI_VERSION/\${VISITARCH})" \
        >> $HOSTCONF
    fi
}

function bv_mili_ensure
{
    if [[ "$DO_MILI" == "yes" ]] ; then
        ensure_built_or_ready "mili" $MILI_VERSION $MILI_BUILD_DIR $MILI_FILE
        if [[ $? != 0 ]] ; then
            warn "Unable to build Mili.  ${MILI_FILE} not found."
            ANY_ERRORS="yes"
            DO_MILI="no"
            if [[ "$DO_SVN" != "yes" ]] ; then
                warn "Note: You have requested to build the Mili library." 
                warn "Mili is not available for public download and" 
                warn "is only available through Subversion access." 
            fi
            error
        fi
    fi
}
# *************************************************************************** #
#                          Function 8.2, build_mili                           #
# *************************************************************************** #

function apply_mili_110_darwin_patch
{
   patch -p0 << \EOF
diff -c a/src/mili_internal.h mili/src/mili_internal.h
*** a/src/mili_internal.h
--- mili/src/mili_internal.h
***************
*** 54,59 ****
--- 54,60 ----
  #include <stdio.h>
  #include <stdlib.h>
  #include <dirent.h>
+ #include <sys/types.h>
  #include "list.h"
  #include "misc.h"
  #include "mili.h"
EOF
   if [[ $? != 0 ]] ; then
      warn "Mili patch for V1.10 on Darwin failed."
      return 1
   fi

   return 0
}

function apply_mili_110_patch
{
   if [[ "$OPSYS" == "Darwin" ]]; then
       apply_mili_110_darwin_patch
       if [[ $? != 0 ]] ; then
           return 1
       fi
   fi

   return 0
}

function apply_mili_patch
{
   if [[ ${MILI_VERSION} == 1.10.0 ]] ; then
       apply_mili_110_patch
       if [[ $? != 0 ]] ; then
           return 1
       fi
   fi

   return 0
}

function build_mili
{
    #
    # Prepare build dir
    #
    prepare_build_dir $MILI_BUILD_DIR $MILI_FILE
    untarred_mili=$?
    if [[ $untarred_mili == -1 ]] ; then
       warn "Unable to prepare Mili Build Directory. Giving Up"
       return 1
    fi

    #
    # Apply patches
    #
    info "Patching mili . . ."
    apply_mili_patch
    if [[ $? != 0 ]] ; then
       if [[ $untarred_mili == 1 ]] ; then
          warn "Giving up on Mili build because the patch failed."
          return 1
       else
          warn "Patch failed, but continuing.  I believe that this script " 
          warn "tried to apply a patch to an existing directory which had " 
          warn "already been patched ... that is, that the patch is "
          warn "failing harmlessly on a second application."
       fi
    fi

    info "Configuring Mili . . ."
    cd $MILI_BUILD_DIR || error "Can't cd to mili build dir."

    info "Invoking command to configure Mili"
    ./configure CXX="$CXX_COMPILER" CC="$C_COMPILER" \
        CFLAGS="$CFLAGS $C_OPT_FLAGS" CXXFLAGS="$CXXFLAGS $CXX_OPT_FLAGS" \
        --prefix="$VISITDIR/mili/$MILI_VERSION/$VISITARCH"
    if [[ $? != 0 ]] ; then
       warn "Mili configure failed.  Giving up"
       return 1
    fi

    #
    # Build Mili
    #
    info "Building Mili . . . (~2 minutes)"

    cd MILI-$OPSYS-*
    cd src
    $C_COMPILER $CFLAGS $C_OPT_FLAGS -D_LARGEFILE64_SOURCE -c \
        mili.c direc.c param.c io.c util.c dep.c svar.c \
        srec.c mesh_u.c wrap_c.c io_mem.c eprtf.c \
        sarray.c gahl.c util.c partition.c ti.c tidirc.c
    if [[ $? != 0 ]] ; then
        warn "Mili build failed.  Giving up"
        return 1
    fi
    #
    # Install into the VisIt third party location.
    #
    info "Installing Mili . . ." 

    mkdir "$VISITDIR/mili"
    mkdir "$VISITDIR/mili/$MILI_VERSION"
    mkdir "$VISITDIR/mili/$MILI_VERSION/$VISITARCH"
    mkdir "$VISITDIR/mili/$MILI_VERSION/$VISITARCH/lib"
    mkdir "$VISITDIR/mili/$MILI_VERSION/$VISITARCH/include"
    cp mili.h mili_enum.h  "$VISITDIR/mili/$MILI_VERSION/$VISITARCH/include"
    if [[ "$DO_STATIC_BUILD" == "no" && "$OPSYS" == "Darwin" ]]; then
        if [[ $ABS_PATH == "yes" ]]; then
           INSTALLNAMEPATH="$VISITDIR/mili/${MILI_VERSION}/$VISITARCH/lib"
        else
           INSTALLNAMEPATH="@executable_path/../lib"
        fi
        $C_COMPILER -dynamiclib -o libmili.$SO_EXT *.o \
          -Wl,-headerpad_max_install_names \
          -Wl,-install_name,$INSTALLNAMEPATH/libmili.${SO_EXT} \
          -Wl,-compatibility_version,$MILI_COMPATIBILITY_VERSION \
          -Wl,-current_version,$MILI_VERSION
        if [[ $? != 0 ]] ; then
          warn "Mili dynamic library build failed.  Giving up"
          return 1
        fi
        cp libmili.$SO_EXT "$VISITDIR/mili/$MILI_VERSION/$VISITARCH/lib"
    else
        ar -rc libmili.a *.o 
        if [[ $? != 0 ]] ; then
          warn "Mili install failed.  Giving up"
          return 1
        fi
        cp libmili.a "$VISITDIR/mili/$MILI_VERSION/$VISITARCH/lib"
    fi

    if [[ "$DO_GROUP" == "yes" ]] ; then
       chmod -R ug+w,a+rX "$VISITDIR/mili"
       chgrp -R ${GROUP} "$VISITDIR/mili"
    fi
    cd "$START_DIR"
    info "Done with Mili"
    return 0
}

function bv_mili_build
{
cd "$START_DIR"
if [[ "$DO_MILI" == "yes" ]] ; then
    check_if_installed "mili" $MILI_VERSION
    if [[ $? == 0 ]] ; then
        info "Skipping Mili build.  Mili is already installed."
    else
        info "Building Mili (~2 minutes)"
        build_mili
        if [[ $? != 0 ]] ; then
            error "Unable to build or install Mili.  Bailing out."
        fi
        info "Done building Mili"
    fi
fi
}

