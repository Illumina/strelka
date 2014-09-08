#
# Starka
# Copyright (c) 2009-2014 Illumina, Inc.
#
# This software is provided under the terms and conditions of the
# Illumina Open Source Software License 1.
#
# You should have received a copy of the Illumina Open Source
# Software License 1 along with this program. If not, see
# <https://github.com/sequencing/licenses/>
#

################################################################################
##
## Configuration file for the opt subfolder (external tools)
##
## author Come Raczy
##
################################################################################

message(STATUS "Creating external tools in subdirectory opt")

include ("${THIS_GLOBALS_CMAKE}")

# convenience macro to set in both current and parent scope:
macro(superset symbol value)
    set(${symbol} "${value}")
    set(${symbol} "${value}" PARENT_SCOPE)
endmacro()

#
# htslib
#
set(HTSLIB_PREFIX "htslib-1.0")
superset(HTSLIB_DIR "${CMAKE_CURRENT_BINARY_DIR}/${HTSLIB_PREFIX}")
superset(HTSLIB_LIBRARY "${HTSLIB_DIR}/libhts.a")

add_custom_command(
    OUTPUT ${HTSLIB_LIBRARY}
    COMMAND ${CMAKE_COMMAND} -E remove_directory "${HTSLIB_DIR}"
    COMMAND ${CMAKE_COMMAND} -E tar xjf "${THIS_REDIST_DIR}/${HTSLIB_PREFIX}.tar.bz2"
    COMMAND $(MAKE) -C "${HTSLIB_DIR}" lib-static tabix bgzip 2>| htslib.build.log
    COMMENT "Building htslib library")
add_custom_target(STARKA_HTSLIB DEPENDS "${HTSLIB_LIBRARY}")

install(PROGRAMS "${HTSLIB_DIR}/tabix" DESTINATION "${THIS_LIBEXECDIR}")
install(PROGRAMS "${HTSLIB_DIR}/bgzip" DESTINATION "${THIS_LIBEXECDIR}")

#
# samtools
#
set(SAMTOOLS_PREFIX "samtools-1.0")
superset(SAMTOOLS_DIR "${CMAKE_CURRENT_BINARY_DIR}/${SAMTOOLS_PREFIX}")
superset(SAMTOOLS_LIBRARY "${SAMTOOLS_DIR}/libbam.a")
superset(SAMTOOLS_PROG "${SAMTOOLS_DIR}/samtools")

add_custom_command(
    OUTPUT ${SAMTOOLS_PROG}
    COMMAND ${CMAKE_COMMAND} -E remove_directory "${SAMTOOLS_DIR}"
    COMMAND ${CMAKE_COMMAND} -E tar xjf "${THIS_REDIST_DIR}/${SAMTOOLS_PREFIX}.tar.bz2"
    COMMAND ${CMAKE_COMMAND} -E chdir "${CMAKE_CURRENT_BINARY_DIR}" ln -sf "${HTSLIB_PREFIX}" htslib
    COMMAND ${CMAKE_COMMAND} -E chdir "${SAMTOOLS_DIR}" $(MAKE) all 2>| samtools.build.log
    DEPENDS "${HTSLIB_LIBRARY}"
    COMMENT "Building samtools package")
add_custom_target(STARKA_SAMTOOLS DEPENDS "${SAMTOOLS_PROG}")

install(PROGRAMS "${SAMTOOLS_PROG}" DESTINATION "${THIS_LIBEXECDIR}")

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
add_custom_target(STARKA_BGZFX DEPENDS "${BGCAT_PROG}")

install(PROGRAMS "${BGCAT_PROG}" DESTINATION "${THIS_LIBEXECDIR}")
install(PROGRAMS "${BGZIP9_PROG}" DESTINATION "${THIS_LIBEXECDIR}")

#
# codemin
#
set(CODEMIN_PREFIX "CodeMin-1.0.2")
superset(CODEMIN_DIR "${CMAKE_CURRENT_BINARY_DIR}/${CODEMIN_PREFIX}/include")
set(CODEMIN_MARK "${CODEMIN_DIR}/minimize_1d.h")

add_custom_command(
    OUTPUT ${CODEMIN_MARK}
    COMMAND ${CMAKE_COMMAND} -E remove_directory "${CODEMIN_DIR}"
    COMMAND ${CMAKE_COMMAND} -E tar xjf "${THIS_REDIST_DIR}/${CODEMIN_PREFIX}.tar.bz2"
    COMMENT "Unpacking CodeMin")
add_custom_target(STARKA_CODEMIN DEPENDS "${CODEMIN_MARK}")

#
# pyflow
#

#make sure unpacking happens regardless of whether ${VCF_TOOLS_DIR} exists
set(PYFLOW_PREFIX "pyflow-1.1.2")
set(PYFLOW_DIR "${CMAKE_CURRENT_BINARY_DIR}/${PYFLOW_PREFIX}")
set(PYFLOW_SCRIPT "${PYFLOW_DIR}/src/pyflow.py")

add_custom_command(
    OUTPUT ${PYFLOW_SCRIPT}
    COMMAND ${CMAKE_COMMAND} -E remove_directory "${PYFLOW_DIR}"
    COMMAND ${CMAKE_COMMAND} -E tar xjf "${THIS_REDIST_DIR}/${PYFLOW_PREFIX}.tar.bz2"
    COMMAND ${CMAKE_COMMAND} -E remove -f "${PYFLOW_DIR}/src/__init__.py"
    COMMAND ${PYTHON_EXECUTABLE} -m compileall -q "${PYFLOW_DIR}/src"
    COMMENT "Extracting pyflow")

add_custom_target(STARKA_PYFLOW DEPENDS "${PYFLOW_SCRIPT}")

set (PYGLOBS "*.py" "*.pyc")
foreach (PYGLOB ${PYGLOBS})
    install (DIRECTORY ${PYFLOW_DIR}/src/ DESTINATION ${THIS_PYTHON_LIBDIR}/pyflow FILES_MATCHING PATTERN ${PYGLOB})
endforeach ()


# tie results back to parent:
#
add_dependencies(STARKA_OPT STARKA_SAMTOOLS STARKA_BGZFX STARKA_CODEMIN STARKA_PYFLOW)


#
# hoedown
#
set(HOEDOWN_PREFIX "hoedown-1.0.2_with_allext")
set(HOEDOWN_DIR "${CMAKE_CURRENT_BINARY_DIR}/${HOEDOWN_PREFIX}")
superset(HOEDOWN_PROG "${HOEDOWN_DIR}/hoedown_allext")

add_custom_command(
    OUTPUT ${HOEDOWN_PROG}
    COMMAND ${CMAKE_COMMAND} -E remove_directory "${HOEDOWN_DIR}"
    COMMAND ${CMAKE_COMMAND} -E tar xjf "${THIS_REDIST_DIR}/${HOEDOWN_PREFIX}.tar.bz2"
    COMMAND $(MAKE) -C "${HOEDOWN_DIR}"
    COMMENT "Building hoedown markdown parser")
add_custom_target(STARKA_MARKDOWN DEPENDS "${HOEDOWN_PROG}")

