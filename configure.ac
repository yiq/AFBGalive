# Perlude.
AC_INIT([AFBGalive], [1.0], [yi.qiao@genetics.utah.edu])
AC_CONFIG_AUX_DIR([build-aux])
AC_CONFIG_MACRO_DIR([m4])
AM_INIT_AUTOMAKE([foreign -Wall -Werror])

# Fixing file creation time
find ${srcdir} -name aclocal.m4 -exec touch {} \;
find ${srcdir} -name \*.in -exec touch {} \;
find ${srcdir} -name configure -exec touch {} \;

# Checks for programs
AC_PROG_CC
AC_PROG_CXX
AM_PROG_AR
LT_PREREQ([2.2])
LT_INIT()
AC_LANG([C++])
AX_CXX_COMPILE_STDCXX_11

AX_PTHREAD

AC_CONFIG_HEADERS([config.h])

## START: DEBUG / RELEASE ##
AC_ARG_ENABLE([debug],
			  [AS_HELP_STRING([--enable-debug],
							  [enable debugging, default: no])],
			  [case "${enableval}" in
			   yes) debug=true ;;
		       no)  debug=false ;;
		       *)   AC_MSG_ERROR([bad value ${enableval} for --enable-debug]) ;;
		       esac],
			   [debug=false])
## END: DEBUG / RELEASE ##

## Options
AC_ARG_WITH([redis-hostname],
			[AS_HELP_STRING([--with-redis-hostname=<hostname>],
							[Specify a hostname for the redis server. Default: redis])],
							[redis_host=${withval}],
							[redis_host=redis])
AC_DEFINE_UNQUOTED([REDIS_HOST], ["${redis_host}"], [The hostname of the redis server])

AC_ARG_WITH([redis-port],
			[AS_HELP_STRING([--with-redis-port=<port>],
							[Specify a port to connect to the redis server. Default: 6379])],
							[redis_port=${withval}],
							[redis_port=6379])
AC_DEFINE_UNQUOTED([REDIS_PORT], [${redis_port}], [The port of the redis server])


AS_IF([test "x$debug" = xtrue], [AC_DEFINE([LOGLVL], [LOGLV_DEBUG], [Define log level to be debug for debug mode])])

AC_CONFIG_SUBDIRS([vendor/sqlite3 vendor/htslib-1.2.1])

AC_CONFIG_FILES([Makefile src/Makefile])
AC_OUTPUT

echo \
"----------------------------------------------
${PACKAGE_NAME} Version ${PACKAGE_VERSION}

Prefix: '${prefix}'
Compiler: '${CC}'
CFLAGS: '${CFLAGS}'
CPPFLAGS: '${CPPFLAGS}'
CXXFLAGS: '${CXXFLAGS}'

Jansson path: '${jansson_path}'

Extra variables:
LDFLAGS='${LDFLAGS}'
LDADD='${LDADD}'
LIBS='${LIBS}'

PTHREAD_LIBS='${PTHREAD_LIBS}'
PTHREAD_CFLAGS='${PTHREAD_CFLAGS}'
PTHREAD_CC='${PTHREAD_CC}'

DATADIR=${DATADIR}
debug=${debug}
redis_host=${redis_host}
"
