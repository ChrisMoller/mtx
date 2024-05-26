
lib_LTLIBRARIES = libmtx.la

libmtx_la_SOURCES = mtx.cc eigens.cc Matrix.cc
libmtx_la_CPPFLAGS = -I$(APL_SOURCES) -I$(APL_SOURCES)/src

noinst_LTLIBRARIES =

#BUILT_SOURCES = gitversion.h

.FORCE:

gitversion.h : .FORCE
	git log  -1 --pretty='format:#define GIT_VERSION "%h %ci"' > xxx
	echo "" >> xxx
	cmp -s xxx $@; \
	rc=$$?; \
	if [ $$rc -eq 0 ]; then rm xxx; else mv xxx $@; fi
