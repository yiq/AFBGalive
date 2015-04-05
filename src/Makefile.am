AM_CPPFLAGS=-I$(top_srcdir)/vendor/YiCppLib/RedisKVStore/src \
			-I$(top_srcdir)/vendor/htslib-1.2.1

AM_CFLAGS=${PTHREAD_CFLAGS}
AM_CXXFLAGS=${AM_CFLAGS}
AM_LIBS=${PTHREAD_LIBS}

AM_LDFLAGS=-all-static

#noinst_LIBRARIES=libtaxangler.a
#libtaxangler_a_SOURCES=TaxDB.h \
#					   TaxDB.cc \
#					   TaxDBManager.h\
#					   TaxDBManager.cc \
#					   TaxResult.h \
#					   TaxResult.cc \
#					   common.h \
#					   log.h \
#					   util_result_formatter.cc
#
#bin_PROGRAMS=taxangler
#taxangler_SOURCES=main.cc
#taxangler_LDADD=libtaxangler.a $(top_builddir)/vendor/YiCppLib/RedisKVStore/src/libyi_rediskvstore.la

bin_PROGRAMS=afbgalive vcfload
afbgalive_SOURCES=main.cc
afbgalive_LDADD=$(top_builddir)/vendor/YiCppLib/RedisKVStore/src/libyi_rediskvstore.la \
				$(top_builddir)/vendor/htslib-1.2.1/libhts.la -lz

vcfload_SOURCES=vcfload.cc
vcfload_LDADD=$(top_builddir)/vendor/YiCppLib/RedisKVStore/src/libyi_rediskvstore.la \
			  $(top_builddir)/vendor/htslib-1.2.1/libhts.la -lz


