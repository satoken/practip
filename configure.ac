#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([practip], [0.0.1], [satoken@bio.keio.ac.jp])
AM_INIT_AUTOMAKE([-Wall foreign])
AM_CONFIG_HEADER(config.h)
AC_CONFIG_SRCDIR([src/practip.cpp])

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
AC_LANG(C++)
AC_CXX_HAVE_STL
AC_PROG_RANLIB
AX_CXX_COMPILE_STDCXX_11
#AX_CXX_HEADER_UNORDERED_MAP

# debug option
AC_ARG_ENABLE([debug], AS_HELP_STRING([--enable-debug], [turn on debugging (default: off)]))
if test "$enable_debug" == "yes"; then
  CFLAGS="-DDEBUG $CFLAGS"
  CXXFLAGS="-DDEBUG $CXXFLAGS"
else
  CFLAGS="-DNDEBUG $CFLAGS"
  CXXFLAGS="-DNDEBUG $CXXFLAGS"
fi

# Checks for libraries.
# FIXME: Replace `main' with a function in `-lglpk':
# glpk
AC_ARG_WITH([glpk], AS_HELP_STRING([--with-glpk], [support GNU linear programming kit]),
  [
    if test "$with_glpk" != "no"; then
      CPPFLAGS="-DWITH_GLPK $CPPFLAGS"
      AC_CHECK_HEADERS([glpk.h])
      AC_CHECK_LIB([z], [main])
      AC_CHECK_LIB([gmp], [main])
      AC_CHECK_LIB([dl], [main])
      AC_CHECK_LIB([ltdl], [main])
      AC_CHECK_LIB([glpk], [glp_create_prob])
    fi
  ]
)

# ILOG CPLEX
AC_ARG_WITH([cplex], AS_HELP_STRING([--with-cplex], [support CPLEX optimizer]),
  [
    if test "$with_cplex" != "no"; then
      CPPFLAGS="-DWITH_CPLEX -DIL_STD $CPPFLAGS"
      AC_CHECK_HEADERS([ilconcert/ilomodel.h])
      AC_CHECK_LIB([concert], [main])
      AC_CHECK_HEADERS([ilcplex/ilocplex.h])
      AC_CHECK_LIB([cplex], [main])
      AC_CHECK_LIB([ilocplex], [main])
    fi
  ]
)

# gurobi optimizer
AC_ARG_WITH([gurobi], AS_HELP_STRING([--with-gurobi], [support Gurobi optimizer]),
  [
    if test "$with_gurobi" != "no"; then
      CPPFLAGS="-DWITH_GUROBI $CPPFLAGS"
      if test "$with_gurobi" != "yes"; then
        CPPFLAGS="-I${with_gurobi}/include $CPPFLAGS"
      	LDFLAGS="-L${with_gurobi}/lib $LDFLAGS"
      fi
      AC_CHECK_HEADERS([gurobi_c++.h])
      AC_CHECK_LIB([gurobi], [GRBloadenv], LIBS="-lgurobi_c++ -lgurobi $LIBS")
      AC_CHECK_LIB([gurobi20], [GRBloadenv], LIBS="-lgurobi_c++ -lgurobi20 $LIBS")
      AC_CHECK_LIB([gurobi30], [GRBloadenv], LIBS="-lgurobi_c++ -lgurobi30 $LIBS")
      AC_CHECK_LIB([gurobi40], [GRBloadenv], LIBS="-lgurobi_c++ -lgurobi40 $LIBS")
      AC_CHECK_LIB([gurobi46], [GRBloadenv], LIBS="-lgurobi_c++ -lgurobi46 $LIBS")
    fi
  ]
)


# Checks for header files.
AC_CHECK_HEADERS([stdlib.h string.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_CHECK_HEADER_STDBOOL
AC_C_INLINE

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_STRTOD
AC_CHECK_FUNCS([sqrt strchr strtol])

AC_CONFIG_FILES([Makefile src/Makefile])
AC_OUTPUT
