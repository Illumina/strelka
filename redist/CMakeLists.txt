#
# Strelka - Small Variant Caller
# Copyright (c) 2009-2018 Illumina, Inc.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#

################################################################################
##
## Configuration file for external tools
##
## author Come Raczy
##
################################################################################

message(STATUS "Building external tools")

include ("${THIS_GLOBALS_CMAKE}")
include ("${THIS_MACROS_CMAKE}")

# convenience macro to set in both current and parent scope:
macro(superset symbol value)
    set(${symbol} "${value}")
    set(${symbol} "${value}" PARENT_SCOPE)
endmacro()


set (THIS_CMAKE_PLATFORM "")
if (WIN32)
    # Sometimes the platform (x64,win32...) is specified as part of the target name
    # this captures the cases where it is specified separately:
    if (CMAKE_GENERATOR_PLATFORM)
        set (THIS_CMAKE_PLATFORM -A "${CMAKE_GENERATOR_PLATFORM}")
    endif ()
endif ()


#
# Find packages required by htslib/samtools:
#
# These packages are needed both for building redist executables and for linking against redist libraries.
# For the latter purpose library information is transmited to PARENT for inclusion in any project binaries
# linking htslib.
#
find_package(ZLIB)
if    (ZLIB_FOUND)
    set (THIS_ADDITIONAL_LIB ${THIS_ADDITIONAL_LIB} z)
    message(STATUS "zlib found")
else  ()
    set(TMP_MSG "zlib not found")
    if (WIN32)
        set (TMP_MSG "${TMP_MSG}. On win32 this can be installed as part of GnuWin32: http://gnuwin32.sourceforge.net/downlinks/zlib.php")
    endif ()
    message(FATAL_ERROR "${TMP_MSG}")
endif ()

superset(THIS_ADDITIONAL_LIB "${THIS_ADDITIONAL_LIB}")


#
# htslib
#
set(HTSLIB_PREFIX "htslib-1.7-6-g6d2bfb7")
superset(HTSLIB_DIR "${CMAKE_CURRENT_BINARY_DIR}/${HTSLIB_PREFIX}")
superset(HTSLIB_LIBRARY "${HTSLIB_DIR}/libhts.a")

## It is occasionally useful to extend debug/asan build options from manta all the way down through htslib,
## the flags below can be uncommented to do so.
##
## Addresss sanitizer build options for htslib/samtools
#set(HTSLIB_CFLAGS '-O0 -g -fsanitize=address -fno-omit-frame-pointer -fno-optimize-sibling-calls')
## Debug build optoins for htslib/samtools
#set(HTSLIB_CFLAGS '-O0 -g')
#set(HTSLIB_CONFIGOPT CFLAGS=${HTSLIB_CFLAGS} LDFLAGS=${HTSLIB_CFLAGS})

add_custom_command(
    OUTPUT ${HTSLIB_DIR}
    COMMAND ${CMAKE_COMMAND} -E remove_directory "${HTSLIB_DIR}"
    COMMAND ${CMAKE_COMMAND} -E tar xjf "${THIS_REDIST_DIR}/${HTSLIB_PREFIX}.tar.bz2"
    COMMENT "Unpacking htslib library")

set (HTSLIB_FINAL_TASK ${HTSLIB_DIR})

if (NOT WIN32)
    add_custom_command(
        OUTPUT ${HTSLIB_LIBRARY}
        COMMAND ./configure CC=${CMAKE_C_COMPILER} ${HTSLIB_CONFIGOPT} --disable-libcurl --disable-bz2 --disable-lzma >htslib.config.log 2>htslib.config.error.log || \(echo "Failed to configure htslib, see ${HTSLIB_DIR}/htslib.config.error.log for details." && exit 1\)
        COMMAND $(MAKE) lib-static bgzip htsfile tabix >htslib.build.log 2>htslib.build.error.log || \(echo "Failed to build htslib, see ${HTSLIB_DIR}/htslib.build.error.log for details." && exit 1\)
        WORKING_DIRECTORY ${HTSLIB_DIR}
        DEPENDS ${HTSLIB_DIR}
        COMMENT "Building htslib library")

    set (HTSLIB_FINAL_TASK ${HTSLIB_LIBRARY})
endif ()

set(THIS_HTSLIB "${THIS_PROJECT_NAME}_htslib")
add_custom_target(${THIS_HTSLIB} DEPENDS "${HTSLIB_FINAL_TASK}")


if (NOT WIN32)
    install(PROGRAMS "${HTSLIB_DIR}/bgzip" DESTINATION "${THIS_LIBEXECDIR}")
    install(PROGRAMS "${HTSLIB_DIR}/htsfile" DESTINATION "${THIS_LIBEXECDIR}")
    install(PROGRAMS "${HTSLIB_DIR}/tabix" DESTINATION "${THIS_LIBEXECDIR}")
endif ()

#
# samtools
#
set(SAMTOOLS_PREFIX "samtools-1.6")
set(SAMTOOLS_DIR "${CMAKE_CURRENT_BINARY_DIR}/${SAMTOOLS_PREFIX}")
set(SAMTOOLS_LIBRARY "${SAMTOOLS_DIR}/libbam.a")
superset(SAMTOOLS_PROG "${SAMTOOLS_DIR}/samtools")

# final directory copy below would ideally be a soft-link, copy is for windows build
add_custom_command(
    OUTPUT ${SAMTOOLS_DIR}
    COMMAND ${CMAKE_COMMAND} -E remove_directory "${SAMTOOLS_DIR}"
    COMMAND ${CMAKE_COMMAND} -E tar xjf "${THIS_REDIST_DIR}/${SAMTOOLS_PREFIX}.tar.bz2"
    DEPENDS ${HTSLIB_FINAL_TASK}
    COMMENT "Unpacking samtools package")

set (SAMTOOLS_FINAL_TASK ${SAMTOOLS_DIR})

if (NOT WIN32)
    add_custom_command(
        OUTPUT ${SAMTOOLS_PROG}
        COMMAND ./configure CC=${CMAKE_C_COMPILER} ${HTSLIB_CONFIGOPT} --without-curses --with-htslib=${HTSLIB_DIR} >samtools.config.log 2>samtools.config.error.log || \(echo "Failed to configure samtools, see ${SAMTOOLS_DIR}/samtools.config.error.log for details." && exit 1\)
        COMMAND $(MAKE) all >samtools.build.log 2>samtools.build.error.log  || \(echo "Failed to build samtools, see ${SAMTOOLS_DIR}/samtools.build.error.log for details." && exit 1\)
        DEPENDS ${HTSLIB_LIBRARY}
        DEPENDS ${SAMTOOLS_DIR}
        WORKING_DIRECTORY ${SAMTOOLS_DIR}
        COMMENT "Building samtools package")

    set (SAMTOOLS_FINAL_TASK ${SAMTOOLS_PROG})
endif ()

set(THIS_SAMTOOLS "${THIS_PROJECT_NAME}_samtools")
add_custom_target(${THIS_SAMTOOLS} DEPENDS "${SAMTOOLS_FINAL_TASK}")

if (NOT WIN32)
    install(PROGRAMS "${SAMTOOLS_PROG}" DESTINATION "${THIS_LIBEXECDIR}")
endif ()

if (NOT WIN32)
    #
    # bgzf_extras
    #
    set(BGZFX_PREFIX "bgzf_extras-v0.1")
    superset(BGZFX_DIR "${CMAKE_CURRENT_BINARY_DIR}/${BGZFX_PREFIX}")
    set(BGCAT_PROG "${BGZFX_DIR}/bgzf_cat")
    set(BGZIP9_PROG "${BGZFX_DIR}/bgzip9")

    add_custom_command(
        OUTPUT ${BGCAT_PROG}
        COMMAND ${CMAKE_COMMAND} -E remove_directory "${BGZFX_DIR}"
        COMMAND ${CMAKE_COMMAND} -E tar xjf "${THIS_REDIST_DIR}/${BGZFX_PREFIX}.tar.bz2"
        COMMAND $(MAKE) -C "${BGZFX_DIR}" all 2>| bgzf_extras.build.log
        COMMENT "Building bgzf_extras")

    set(THIS_BGZFX "${THIS_PROJECT_NAME}_bgzfx")
    add_custom_target(${THIS_BGZFX} DEPENDS "${BGCAT_PROG}")

    install(PROGRAMS "${BGCAT_PROG}" DESTINATION "${THIS_LIBEXECDIR}")
    install(PROGRAMS "${BGZIP9_PROG}" DESTINATION "${THIS_LIBEXECDIR}")
endif ()

#
# rapidjson
#
set(RAPIDJSON_PREFIX "rapidjson-1.1.0")
set(RAPIDJSON_DIR "${CMAKE_CURRENT_BINARY_DIR}/${RAPIDJSON_PREFIX}")
superset(RAPIDJSON_INCLUDE_DIR "${RAPIDJSON_DIR}/include")


add_custom_command(
        OUTPUT ${RAPIDJSON_INCLUDE_DIR}
        COMMAND ${CMAKE_COMMAND} -E remove_directory "${RAPIDJSON_DIR}"
        COMMAND ${CMAKE_COMMAND} -E tar xjf "${THIS_REDIST_DIR}/${RAPIDJSON_PREFIX}.tar.bz2"
        COMMENT "Unpacking rapidjson package")

set (RAPIDJSON_FINAL_TASK ${RAPIDJSON_INCLUDE_DIR})

set(THIS_RAPIDJSON "${THIS_PROJECT_NAME}_rapidjson")
add_custom_target(${THIS_RAPIDJSON} DEPENDS "${RAPIDJSON_FINAL_TASK}")



#
# codemin
#
set(CODEMIN_PREFIX "CodeMin-1.0.5")
superset(CODEMIN_DIR "${CMAKE_CURRENT_BINARY_DIR}/${CODEMIN_PREFIX}/include")
set(CODEMIN_MARK "${CODEMIN_DIR}/minimize_1d.h")

add_custom_command(
    OUTPUT ${CODEMIN_MARK}
    COMMAND ${CMAKE_COMMAND} -E remove_directory "${CODEMIN_DIR}"
    COMMAND ${CMAKE_COMMAND} -E tar xjf "${THIS_REDIST_DIR}/${CODEMIN_PREFIX}.tar.bz2"
    COMMENT "Unpacking CodeMin")

set(THIS_CODEMIN "${THIS_PROJECT_NAME}_codemin")
add_custom_target(${THIS_CODEMIN} DEPENDS "${CODEMIN_MARK}")

#
# pyflow
#

set(PYFLOW_PREFIX "pyflow-1.1.20")
set(PYFLOW_DIR "${CMAKE_CURRENT_BINARY_DIR}/${PYFLOW_PREFIX}")
set(PYFLOW_SCRIPT "${PYFLOW_DIR}/src/pyflow.py")

set (PYFLOW_DEPENDS "${PYFLOW_SCRIPT}")
add_custom_command(
    OUTPUT ${PYFLOW_SCRIPT}
    COMMAND ${CMAKE_COMMAND} -E remove_directory "${PYFLOW_DIR}"
    COMMAND ${CMAKE_COMMAND} -E tar xjf "${THIS_REDIST_DIR}/${PYFLOW_PREFIX}.tar.bz2"
    COMMAND ${CMAKE_COMMAND} -E remove -f "${PYFLOW_DIR}/src/__init__.py"
    COMMENT "Unpacking pyflow")

if (PYTHONINTERP_FOUND)
    set (PYFLOW_DEPENDS "${PYFLOW_SCRIPT}c")
    add_custom_command(
        OUTPUT ${PYFLOW_SCRIPT}c
        COMMAND ${PYTHON_EXECUTABLE} -m compileall -q "${PYFLOW_DIR}/src"
        DEPENDS ${PYFLOW_SCRIPT}
        COMMENT "Building pyflow")
endif()

install_python_lib_dir("${PYFLOW_DIR}/src" "${THIS_PYTHON_LIBDIR}/pyflow")

set (THIS_PYFLOW "${THIS_PROJECT_NAME}_pyflow")
add_custom_target(${THIS_PYFLOW} DEPENDS "${PYFLOW_DEPENDS}")


# tie results back to parent:
#
add_dependencies(${THIS_OPT} ${THIS_SAMTOOLS} ${THIS_RAPIDJSON} ${THIS_CODEMIN} ${THIS_PYFLOW})

if (NOT WIN32)
    add_dependencies(${THIS_OPT} ${THIS_BGZFX})
endif ()

