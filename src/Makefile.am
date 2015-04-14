AM_CPPFLAGS=-I$(top_srcdir)/vendor/YiCppLib/RedisKVStore/src \
			-I$(top_srcdir)/vendor/htslib-1.2.1

AM_CFLAGS=${PTHREAD_CFLAGS} -std=gnu99
AM_CXXFLAGS=${PTHREAD_CFLAGS}
AM_LIBS=${PTHREAD_LIBS}

#AM_LDFLAGS=-all-static

bin_PROGRAMS=afbgalive vcfload
afbgalive_SOURCES=main.cc h2ip.c
afbgalive_LDADD=$(top_builddir)/vendor/YiCppLib/RedisKVStore/src/libyi_rediskvstore.la \
				$(top_builddir)/vendor/htslib-1.2.1/libhts.la -lz

vcfload_SOURCES=vcfload.cc h2ip.c
vcfload_LDADD=$(top_builddir)/vendor/YiCppLib/RedisKVStore/src/libyi_rediskvstore.la \
			  $(top_builddir)/vendor/htslib-1.2.1/libhts.la -lz


