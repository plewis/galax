AC_INIT([galax], [1.2.0], [bug-automake@gnu.org])
AM_INIT_AUTOMAKE([-Wall -Werror foreign])
LT_INIT
AC_CONFIG_MACRO_DIR([m4])
AX_BOOST_BASE([1.55],, [AC_MSG_ERROR([galax needs Boost, but it was not found in your system])])
AX_BOOST_REGEX
AX_BOOST_PROGRAM_OPTIONS
AC_PROG_CXX
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([
 Makefile
 src/Makefile
])
AC_OUTPUT
