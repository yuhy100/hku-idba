AM_CXXFLAGS = -Wall -O2 -fopenmp
AM_LDFLAGS = -lboost_program_options -lboost_filesystem
INCLUDES = -I$(top_srcdir)/src/common/include -I$(top_srcdir)/src/assembly/include
#LIBS = $(top_srcdir)/lib/libcommon.a 
LIBS = $(top_srcdir)/lib/libassembly.a

bin_PROGRAMS = idba
noinst_PROGRAMS = calcCoverBlat \
		  calcCover \
		  calcCoverMummer \
		  compare \
		  compareKmer \
		  fa2fq \
		  fq2fa \
		  genTaxonDatabase \
		  groupDistance \
		  mergeReads \
		  metaCluster \
		  normReads \
		  rawN50 \
		  simReads \
		  test \
		  metaClusterHi \
		  printClusterHi 



idba_SOURCES = $(top_srcdir)/src/tools/idba.cpp

calcCoverBlat_SOURCES = $(top_srcdir)/src/tools/calcCoverBlat.cpp
calcCover_SOURCES = $(top_srcdir)/src/tools/calcCover.cpp
calcCoverMummer_SOURCES = $(top_srcdir)/src/tools/calcCoverMummer.cpp
compare_SOURCES = $(top_srcdir)/src/tools/compare.cpp
compareKmer_SOURCES = $(top_srcdir)/src/tools/compareKmer.cpp
fa2fq_SOURCES = $(top_srcdir)/src/tools/fa2fq.cpp
fq2fa_SOURCES = $(top_srcdir)/src/tools/fq2fa.cpp
genTaxonDatabase_SOURCES = $(top_srcdir)/src/tools/genTaxonDatabase.cpp
groupDistance_SOURCES = $(top_srcdir)/src/tools/groupDistance.cpp
mergeReads_SOURCES = $(top_srcdir)/src/tools/mergeReads.cpp
metaCluster_SOURCES = $(top_srcdir)/src/tools/metaCluster.cpp
normReads_SOURCES = $(top_srcdir)/src/tools/normReads.cpp
rawN50_SOURCES = $(top_srcdir)/src/tools/rawN50.cpp
simReads_SOURCES = $(top_srcdir)/src/tools/simReads.cpp
test_SOURCES = $(top_srcdir)/src/tools/test.cpp
metaClusterHi_SOURCES = $(top_srcdir)/src/tools/metaClusterHi.cpp
printClusterHi_SOURCES = $(top_srcdir)/src/tools/printClusterHi.cpp

