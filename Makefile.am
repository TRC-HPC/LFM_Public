#
# Copyright (C) Huawei Technologies Co., Ltd. 2020. ALL RIGHTS RESERVED.
#
# See file LICENSE for terms.
#


AM_LDFLAGS                   = -lm

ACLOCAL_AMFLAGS = -I m4

lib_LTLIBRARIES              = libfastmesh.la
libfastmesh_la_LDFLAGS       = -ldl -version-info $(SOVERSION)
libfastmesh_la_SOURCES       =  src/mesh_reader.cpp src/mesh_solver.cpp			\
								src/mpi_env.cpp 				\
								src/cfd_v0.cpp src/inputReader.cpp
libfastmesh_ladir            = $(includedir)
dist_libfastmesh_la_HEADERS  = api/fastmesh.h api/mpi_env.h						\
			       api/cfdv0_elements.h				\
			       api/cfdv0_solver.h \
			       api/polyMeshReaderOF.h api/dictReaderOF.h api/inputReader.h api/runTimeManager.h

bin_PROGRAMS                 = lfm_solve


# Solver binary:
lfm_solve_LDADD = libfastmesh.la -L$(FOAM_USER_LIBBIN) -lpolyMeshReaderOF -ldictReaderOF -lrunTimeManagerOF
#lfm_solve_LDADD = libfastmesh.la -L$(FOAM_USER_LIBBIN) -lpolyMeshReaderOF -ldictReaderOF -lrunTimeManagerOF -L$(EXTRAE_HOME)/lib -lmpitrace
lfm_solve_SOURCES = info/lfm_solve.cpp

# Examples
examplesdir = examples

install-data-local:
	test -d "$(prefix)/$(examplesdir)" || cp -r $(examplesdir) "$(prefix)/"

uninstall-local:
	rm -rf "$(prefix)/$(examplesdir)"

EXTRA_DIST = examples
