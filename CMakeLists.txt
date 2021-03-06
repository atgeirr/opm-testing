# -*- mode: cmake; tab-width: 2; indent-tabs-mode: t; truncate-lines: t; compile-command: "cmake -Wdev" -*-
# vim: set filetype=cmake autoindent tabstop=2 shiftwidth=2 noexpandtab softtabstop=2 nowrap:

###########################################################################
#                                                                         #
# Note: The bulk of the build system is located in the cmake/ directory.  #
#       This file only contains the specializations for this particular   #
#       project. Most likely you are interested in editing one of these   #
#       files instead:                                                    #
#                                                                         #
#       dune.module                              Name and version number  #
#       CMakeLists_files.cmake                   Path of source files     #
#       cmake/Modules/${project}-prereqs.cmake   Dependencies             #
#                                                                         #
###########################################################################

# Mandatory call to project
project(opm-testing C CXX)

cmake_minimum_required (VERSION 2.8)

option(SIBLING_SEARCH "Search for other modules in sibling directories?" ON)
set( USE_OPENMP_DEFAULT OFF ) # Use of OpenMP is considered experimental

if(SIBLING_SEARCH AND NOT opm-common_DIR)
  # guess the sibling dir
  get_filename_component(_leaf_dir_name ${PROJECT_BINARY_DIR} NAME)
  get_filename_component(_parent_full_dir ${PROJECT_BINARY_DIR} DIRECTORY)
  get_filename_component(_parent_dir_name ${_parent_full_dir} NAME)
  #Try if <module-name>/<build-dir> is used
  get_filename_component(_modules_dir ${_parent_full_dir} DIRECTORY)
  if(IS_DIRECTORY ${_modules_dir}/opm-common/${_leaf_dir_name})
    set(opm-common_DIR ${_modules_dir}/opm-common/${_leaf_dir_name})
  else()
    string(REPLACE ${PROJECT_NAME} opm-common _opm_common_leaf ${_leaf_dir_name})
    if(NOT _leaf_dir_name STREQUAL _opm_common_leaf
        AND IS_DIRECTORY ${_parent_full_dir}/${_opm_common_leaf})
      # We are using build directories named <prefix><module-name><postfix>
      set(opm-common_DIR ${_parent_full_dir}/${_opm_common_leaf})
    elseif(IS_DIRECTORY ${_parent_full_dir}/opm-common)
      # All modules are in a common build dir
      set(opm-common_DIR "${_parent_full_dir}/opm-common}")
    endif()
  endif()
endif()
if(opm-common_DIR AND NOT IS_DIRECTORY ${opm-common_DIR})
  message(WARNING "Value ${opm-common_DIR} passed to variable"
    " opm-common_DIR is not a directory")
endif()

find_package(opm-common REQUIRED)

include(OpmInit)

# not the same location as most of the other projects? this hook overrides
macro (dir_hook)
endmacro (dir_hook)

# project information is in dune.module. Read this file and set variables.
# we cannot generate dune.module since it is read by dunecontrol before
# the build starts, so it makes sense to keep the data there then.
include (OpmInit)

# Look for the opm-data repository; if found the variable
# HAVE_OPM_DATA will be set to true.

# list of prerequisites for this particular project; this is in a
# separate file (in cmake/Modules sub-directory) because it is shared
# with the find module
include ("${project}-prereqs")

# read the list of components from this file (in the project directory);
# it should set various lists with the names of the files to include
include (CMakeLists_files.cmake)

macro (config_hook)
	opm_need_version_of ("dune-common")
	opm_need_version_of ("dune-istl")
	opm_need_version_of ("ewoms")
endmacro (config_hook)

macro (prereqs_hook)
endmacro (prereqs_hook)

macro (sources_hook)
endmacro (sources_hook)

macro (fortran_hook)
endmacro (fortran_hook)

macro (files_hook)
endmacro (files_hook)

macro (tests_hook)
endmacro (tests_hook)

# all setup common to the OPM library modules is done here
include (OpmLibMain)


# opm_add_test(flow_simple
#   ONLY_COMPILE
#   ALWAYS_ENABLE
#   DEPENDS "opmsimulators"
#   LIBRARIES "opmsimulators"
#   SOURCES
#     examples/flow_simple.cpp)

opm_add_test(fluid_test
  ONLY_COMPILE
  ALWAYS_ENABLE
  DEPENDS "opmsimulators"
  LIBRARIES "opmsimulators"
  SOURCES
    examples/fluid_test)	

opm_add_test(flow_blackoil_cpr_tester
  ONLY_COMPILE
  ALWAYS_ENABLE
  DEPENDS "opmsimulators"
  LIBRARIES "opmsimulators"
  SOURCES
    examples/flow_blackoil_cpr_tester.cpp)	


opm_add_test(flow_onephase
  ONLY_COMPILE
  ALWAYS_ENABLE
  DEPENDS "opmsimulators"
  LIBRARIES "opmsimulators"
  SOURCES
  examples/flow_onephase.cpp)	
	

opm_add_test(linearsolver_test
  ONLY_COMPILE
  ALWAYS_ENABLE
  DEPENDS "opmsimulators"
  LIBRARIES "opmsimulators"
  SOURCES
    examples/linearsolver_test.cpp)	

opm_add_test(linearsolver_test_mpi
  ONLY_COMPILE
  ALWAYS_ENABLE
  DEPENDS "opmsimulators"
  LIBRARIES "opmsimulators"
  SOURCES
    examples/linearsolver_test_mpi.cpp)	

	
target_link_libraries(linearsolver_test boost_program_options)
target_link_libraries(linearsolver_test_mpi boost_program_options boost_serialization boost_mpi)		
#install(TARGETS flow_simple DESTINATION bin)
install(TARGETS fluid_test DESTINATION bin)
install(TARGETS linearsolver_test DESTINATION bin)
install(TARGETS linearsolver_test_mpi DESTINATION bin)
install(TARGETS flow_blackoil_cpr_tester DESTINATION bin)

	

#include(OpmBashCompletion)
#opm_add_bash_completion(flow_simple)

