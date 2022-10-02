AC_PREREQ([2.63])

define([fastmesh_ver_major], 0)
define([fastmesh_ver_minor], 1)
define([fastmesh_ver_patch], 0)
define([ts], esyscmd([sh -c "date +%Y%m%d%H%M%S"]))

# This is the API version (see libtool library versioning)
# http://www.gnu.org/software/libtool/manual/html_node/Updating-version-info.html
# current:rev:age
define([libfastmesh_so_version], 0:1:0)

AC_INIT([fastmesh], [fastmesh_ver_major.fastmesh_ver_minor])
AC_LANG(C++)
AC_CONFIG_MACRO_DIRS([m4])
AC_CHECK_PROG(GITBIN, git, yes)
AS_IF([test x"${GITBIN}" = x"yes"],
      [AC_SUBST(SCM_VERSION, esyscmd([sh -c 'git rev-parse --short=7 HEAD']))],
      [AC_SUBST(SCM_VERSION, "0000000")])

MAJOR_VERSION=fastmesh_ver_major
MINOR_VERSION=fastmesh_ver_minor
PATCH_VERSION=fastmesh_ver_patch
VERSION=$MAJOR_VERSION.$MINOR_VERSION.$PATCH_VERSION
SOVERSION=libfastmesh_so_version
AC_SUBST(MAJOR_VERSION)
AC_SUBST(MINOR_VERSION)
AC_SUBST(PATCH_VERSION)
AC_SUBST(SCM_VERSION)
AC_SUBST(SOVERSION)

LT_INIT
AM_INIT_AUTOMAKE([subdir-objects])
AC_CONFIG_FILES([Makefile
                 api/version.h])
if test -z $CXXFLAGS; then
    CXXFLAGS='-O3'
fi
AC_PROG_CXX

AC_ARG_ENABLE([mpi],
              [AS_HELP_STRING([--enable-mpi],
                              [build the parallel version of libFastMesh with MPI])])
AS_IF([test "x$enable_mpi" = "xyes"],
      [AX_MPI([AC_SUBST(CXX, "$MPICXX")],
              [AC_MSG_ERROR([Unable to locate MPI C++ compiler wrapper])])])
AM_CONDITIONAL([WITH_MPI], [test "x$enable_mpi" = "xyes"])



AC_ARG_ENABLE(debug,
AS_HELP_STRING([--enable-debug],
               [enable debugging, default: no]),
[case "${enableval}" in
             yes) debug=true ;;
             no)  debug=false ;;
             *)   AC_MSG_ERROR([bad value ${enableval} for --enable-debug]) ;;
esac],
[debug=false])

AM_CONDITIONAL(DEBUG, test x"$debug" = x"true")

AC_OUTPUT
