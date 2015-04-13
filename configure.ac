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

## START: CHECKING FOR JANSSON ##
AC_ARG_WITH([jansson],
			[AS_HELP_STRING([--with-jansson=DIR], 
							[specify the path to search for jansson])],
			[jansson_path=${withval}])

if test "x${jansson_path}" != x; then
	if test -e ${jansson_path}/lib; then LDFLAGS+="-L${jansson_path}/lib"; fi
	if test -e ${jansson_path}/include; then CPPFLAGS+="-I${jansson_path}/include"; fi
fi

have_jansson=no

AC_SEARCH_LIBS([json_object], [jansson], [have_jansson=yes], [have_jansson=no])

if test "x${have_jansson}" = xyes; then
	AC_CHECK_HEADERS([jansson.h], [], [have_jansson=no])
fi

if test "x${have_jansson}" = xno; then
	AC_MSG_ERROR([Unable to locate jansson])
fi

## END: CHECKINF FOR JANSSON ##

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

AS_IF([test "x$debug" = xtrue], [AC_DEFINE([LOGLVL], [LOGLV_DEBUG], [Define log level to be debug for debug mode])])

AC_CONFIG_SUBDIRS([vendor/YiCppLib/RedisKVStore vendor/htslib-1.2.1])

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
