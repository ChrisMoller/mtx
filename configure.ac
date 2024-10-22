AC_INIT([APL matrix native functions], [2.1.0], [henrik@henrimoller.me],[mtx])
AC_PREFIX_DEFAULT( /usr/local )
AC_USE_SYSTEM_EXTENSIONS
AC_CONFIG_HEADERS([mtx_config.h])
AC_CONFIG_MACRO_DIR([m4])
AM_INIT_AUTOMAKE([-Wall -Werror foreign])
AM_MAINTAINER_MODE([enable])

AC_PROG_CXX
AM_PROG_AR
#AM_PROG_LIBTOOL
LT_INIT

libdir=$libdir/apl

APL_SOURCES=`apl --show_src_dir`/..


AC_CHECK_FILE($APL_SOURCES, [foundit=yes], [foundit=no])
if test x$foundit = xno; then
AC_ARG_WITH([apl-source-dir],
  AS_HELP_STRING([--apl-source-dir=arg],
  [Location of the GNU APL source code.]))
  APL_SOURCES=$with_apl_source_dir
fi

PKG_CHECK_MODULES([GSL], [gsl >= 2.8])
AC_SUBST(GSL_CFLAGS)
AC_SUBST(GSL_LIBS)

#AC_LANG_PUSH([C++])
#AC_CHECK_HEADER(eigen3/Eigen/Dense,
#[AC_DEFINE([HAVE_EIGEN], [1], [Has Eigen3])],
#[
#  echo "";
#  echo "Eigen3 package not found.";
#  echo "If you need eigenvectors or eigenvalues, you need to install the package.";
#  echo "It can be found at https://eigen.tuxfamily.org/index.php?title=Main_Page";
#  echo "";
#])
#AC_LANG_POP([C++])

AC_ARG_WITH([apl], AS_HELP_STRING([--with-apl],
[Specify the APL source directory]))
if test "x$with_apl" != "x"; then
  APL_SOURCES=$with_apl
fi

AC_SUBST(APL_SOURCES)

AC_CONFIG_FILES([
  Makefile
  src/Makefile
])
AC_OUTPUT
