#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([cmap], [0.2], [])
AC_CONFIG_SRCDIR([src/pdb.c])

# Checks for programs.
AC_PROG_CC
AC_PROG_CC_STDC

# Checks for libraries.
AC_CHECK_LIB([m], [log10])
AC_CHECK_LIB([ncurses], [waddwstr],, [
                    AC_CHECK_LIB([ncursesw], [waddwstr],, [
                                            AC_MSG_ERROR([Error: requires ncurses library with wide character support.])
                                                          ])
             ])

# If the ncursesw5-config utility is installed, use it to get a list of
# necessary compiler/linker flags:
command -v ncursesw5-config
if test $? -eq 0 ; then 
        echo Found ncursesw5-config utility; Adding suggested compiler flags.
        AC_SUBST([CURSESFLAGS], [$(ncursesw5-config --cflags)])
        AC_SUBST([CURSESLIBS], [$(ncursesw5-config --libs)])
fi

# Checks for header files.
AC_CHECK_HEADERS([locale.h stdlib.h string.h wchar.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_CHECK_HEADER_STDBOOL
AC_TYPE_SIZE_T

AC_FUNC_MALLOC
AC_CHECK_FUNCS([floor pow setlocale sqrt])

AM_INIT_AUTOMAKE([-Wall foreign])
AC_CONFIG_FILES([Makefile])
AC_OUTPUT
AC_CONFIG_FILES([src/Makefile])
AC_OUTPUT
