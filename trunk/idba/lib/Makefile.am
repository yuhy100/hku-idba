AM_CXXFLAGS = -Wall -O2 -fopenmp
AM_LDFLAGS = -fopenmp -lboost_program_options-mt -lboost_filesystem-mt
INCLUDES = -I$(top_srcdir)/src/common/include -I$(top_srcdir)/src/assembly/include

noinst_LIBRARIES = libcommon.a libassembly.a
libcommon_a_SOURCES = \
	$(top_srcdir)/src/common/include/AbstractNode.h \
	$(top_srcdir)/src/common/include/BitOperation.h \
	$(top_srcdir)/src/common/include/BranchGroup.h \
	$(top_srcdir)/src/common/include/ContigBranchGroup.h \
	$(top_srcdir)/src/common/include/Contig.h \
	$(top_srcdir)/src/common/include/ContigGraph.h \
	$(top_srcdir)/src/common/include/globals.h \
	$(top_srcdir)/src/common/include/HashAlign.h \
	$(top_srcdir)/src/common/include/HashGraph.h \
	$(top_srcdir)/src/common/include/HashNode.h \
	$(top_srcdir)/src/common/include/HashTable.h \
	$(top_srcdir)/src/common/include/Kmer.h \
	$(top_srcdir)/src/common/include/KmerNode.h \
	$(top_srcdir)/src/common/include/KmerVector.h \
	$(top_srcdir)/src/common/include/Log.h \
	$(top_srcdir)/src/common/include/MapGraph.h \
	$(top_srcdir)/src/common/include/Read.h \
	$(top_srcdir)/src/common/include/Reader.h \
	$(top_srcdir)/src/common/include/ScaffoldGraph.h \
	$(top_srcdir)/src/common/include/Sequence.h \
	$(top_srcdir)/src/common/include/Utils.h \
	$(top_srcdir)/src/common/include/Writer.h \
	$(top_srcdir)/src/common/include/ContigNode.h \
	$(top_srcdir)/src/common/AbstractNode.cpp \
	$(top_srcdir)/src/common/BitOperation.cpp \
	$(top_srcdir)/src/common/BranchGroup.cpp \
	$(top_srcdir)/src/common/ContigBranchGroup.cpp \
	$(top_srcdir)/src/common/Contig.cpp \
	$(top_srcdir)/src/common/ContigGraph.cpp \
	$(top_srcdir)/src/common/globals.cpp \
	$(top_srcdir)/src/common/HashAlign.cpp \
	$(top_srcdir)/src/common/HashGraph.cpp \
	$(top_srcdir)/src/common/HashNode.cpp \
	$(top_srcdir)/src/common/HashTable.cpp \
	$(top_srcdir)/src/common/Kmer.cpp \
	$(top_srcdir)/src/common/KmerNode.cpp \
	$(top_srcdir)/src/common/KmerVector.cpp \
	$(top_srcdir)/src/common/Log.cpp \
	$(top_srcdir)/src/common/MapGraph.cpp \
	$(top_srcdir)/src/common/Read.cpp \
	$(top_srcdir)/src/common/Reader.cpp \
	$(top_srcdir)/src/common/ScaffoldGraph.cpp \
	$(top_srcdir)/src/common/Sequence.cpp \
	$(top_srcdir)/src/common/Utils.cpp \
	$(top_srcdir)/src/common/Writer.cpp \
	$(top_srcdir)/src/common/ContigNode.cpp

#libassembly_a_LIBADD = $(top_srcdir)/lib/libcommon.a
#libassembly_a_LDFLAGS = $(top_srcdir)/lib/libcommon.a
libassembly_a_SOURCES = $(libcommon_a_SOURCES) \
	$(top_srcdir)/src/assembly/include/AbstractAssembler.h \
	$(top_srcdir)/src/assembly/include/AssemblyAlgorithm.h \
	$(top_srcdir)/src/assembly/include/GappedScaffoldAssembler.h \
	$(top_srcdir)/src/assembly/include/IterativeAlgorithm.h \
	$(top_srcdir)/src/assembly/include/IterativeAssembler.h \
	$(top_srcdir)/src/assembly/include/ScaffoldAssembler.h \
	$(top_srcdir)/src/assembly/include/SimpleAlgorithm.h \
	$(top_srcdir)/src/assembly/include/SimpleAssembler.h \
	$(top_srcdir)/src/assembly/include/Option.h \
	$(top_srcdir)/src/assembly/AbstractAssembler.cpp \
	$(top_srcdir)/src/assembly/AssemblyAlgorithm.cpp \
	$(top_srcdir)/src/assembly/GappedScaffoldAssembler.cpp \
	$(top_srcdir)/src/assembly/IterativeAlgorithm.cpp \
	$(top_srcdir)/src/assembly/IterativeAssembler.cpp \
	$(top_srcdir)/src/assembly/ScaffoldAssembler.cpp \
	$(top_srcdir)/src/assembly/SimpleAlgorithm.cpp \
	$(top_srcdir)/src/assembly/SimpleAssembler.cpp \
	$(top_srcdir)/src/assembly/Option.cpp
#libassembly_a_LIBADD = $(top_srcdir)/lib/libcommon.a
