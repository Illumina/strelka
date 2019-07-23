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
## Configuration file for boost installation
##
## author Come Raczy
##
################################################################################

# set to TRUE to see more detailed information about the boost find/build procedure:
execute_process(COMMAND uname -p OUTPUT_VARIABLE ARM RESULT_VARIABLE result)
STRING(REGEX REPLACE "(\r?\n)+$" "" ARM "${ARM}")
if ("${ARM}" STREQUAL "aarch64")
    set (aarch64 “ON”)
else()
    set (aarch64 “OFF”)
endif()

set (DEBUG_FINDBOOST FALSE)
if (${DEBUG_FINDBOOST})
    set (Boost_DEBUG "ON")
    set (Boost_DETAILED_FAILURE_MSG "ON")
endif ()

macro (initBoostParams)
    # required boost libraries
    set (THIS_BOOST_VERSION 1.58.0)
    # note we default to alphabetical order here, except where boost libraries depend on
    # each other (timer->chrono)
    set (THIS_BOOST_COMPONENTS date_time filesystem program_options
                               serialization system timer chrono unit_test_framework)

    # the name given to boost.build and the library name are the same for all libraries, except
    # for test, so we need two lists now:
    set (THIS_BOOST_BUILD_COMPONENTS date_time filesystem program_options
                                     serialization system timer chrono test)
    set (Boost_USE_STATIC_LIBS ON)
    if (NOT WIN32)
        # bjam on windows ignores this setting so skip for win32:
        set (Boost_USE_MULTITHREADED OFF)
    endif ()
endmacro()

# simple helper for resetFindBoost
macro(unsetall name)
    unset (${name} CACHE)
    unset (${name})
endmacro()


function(makedir path)
    if(NOT EXISTS "${path}")
        file(MAKE_DIRECTORY "${path}")
    endif()
endfunction()


macro (resetFindBoost)

    set(BOOST_RESET_SYMBOLS FOUND INCLUDE_DIRS INCLUDE_DIR LIBRARIES LIBRARY_DIRS VERSION LIB_VERSION MAJOR_VERSION MINOR_VERSION SUBMINOR_VERSION USE_STATIC_LIBS USE_MULTITHREADED)

    foreach (TAG ${BOOST_RESET_SYMBOLS})
        unsetall (Boost_${TAG})
    endforeach()

    unset (ENV{BOOST_LIBRARYDIR})

    foreach (COMPONENT ${THIS_BOOST_COMPONENTS})
        STRING(TOUPPER ${COMPONENT} UPPERCOMPONENT)
        unsetall (Boost_${UPPERCOMPONENT}_FOUND)
        unsetall (Boost_${UPPERCOMPONENT}_LIBRARY)
        unsetall (Boost_${UPPERCOMPONENT}_LIBRARY_RELEASE)
        unsetall (Boost_${UPPERCOMPONENT}_LIBRARY_DEBUG)
    endforeach ()

    initBoostParams()
endmacro ()


initBoostParams()

if (THIS_FORCE_STATIC_LINK)
    set(Boost_USE_STATIC_LIBS ON)
endif ()

set(BOOST_BOOTSTRAP_INSTALL_DIR ${CMAKE_CURRENT_BINARY_DIR}/bootstrap/boost)
if (EXISTS "${BOOST_BOOTSTRAP_INSTALL_DIR}/boost_install_complete")
    set (BOOST_ROOT "${BOOST_BOOTSTRAP_INSTALL_DIR}")
endif ()

# newer cmake versions (at least cmake 3.5) will issue tons of warning text for
# each component query if no boost version is found at all, for this reason we
# break up the boost version search into two steps:
#
# (1) can you find boost >= min_version at all?
# (2) if so, does that boost version include all of our required components?
#
find_package(Boost ${THIS_BOOST_VERSION})
if (Boost_FOUND)
    find_package(Boost ${THIS_BOOST_VERSION} COMPONENTS ${THIS_BOOST_COMPONENTS} QUIET)
endif ()

# CMAKE_PARALLEL is only used if boost is found, but moving the setting here (outside of the if below) supresses a cmake warning:
if (NOT CMAKE_PARALLEL)
    set (CMAKE_PARALLEL "1")
endif ()

#
# If the right version of boost is not found, it will be built from the distribution
#
if (NOT Boost_FOUND)
    #
    # Provide details on why any user-supplied/auto-detected boost can't be used:
    #
    foreach(COMPONENT ${THIS_BOOST_COMPONENTS})
        STRING(TOUPPER ${COMPONENT} UPPERCOMPONENT)
        if (${Boost_${UPPERCOMPONENT}_FOUND})
            set(FOUND_STATUS "found")
        else()
            set(FOUND_STATUS "")
            if (Boost_USE_STATIC_LIBS)
                set(FOUND_STATUS "STATIC LIBRARY ")
            endif ()
            set(FOUND_STATUS "${FOUND_STATUS}NOT FOUND")
        endif()
        message(STATUS "Boost component: ${COMPONENT}\tstatus: ${FOUND_STATUS}")
    endforeach()

    if (BOOST_ROOT)
        message (STATUS "BOOST_ROOT is set to ${BOOST_ROOT} but boost ${THIS_BOOST_VERSION} or one of its components was not found.")
        message (FATAL_ERROR "Unset BOOST_ROOT or set it to the root location of boost ${THIS_BOOST_VERSION}.")
    endif()

    # Try to find boost in bootstrap installation location
    resetFindBoost()
    message(STATUS "Boost version ${THIS_BOOST_VERSION}+ not found or found without all required components. Boost will be built from source distribution...")

    set(ENV{THIS_BOOST_BUILD_COMPONENTS} "${THIS_BOOST_BUILD_COMPONENTS}")
    set(ENV{THIS_BOOST_VERSION} "${THIS_BOOST_VERSION}")

    set(THIS_BOOTSTRAP_DIR "${THIS_MODULE_DIR}/bootstrap")

    string (REPLACE "." "_" BOOST_FILENAME_VERSION "${THIS_BOOST_VERSION}")
    set (BOOST_SOURCE_PREFIX "boost_${BOOST_FILENAME_VERSION}_subset")
    set (BOOST_INSTALL_DIR "${BOOST_BOOTSTRAP_INSTALL_DIR}")
    set (BOOST_BUILD_DIR "${BOOST_INSTALL_DIR}/build")
    set (BOOST_SRC_DIR "${BOOST_BUILD_DIR}/${BOOST_SOURCE_PREFIX}")
    makedir("${BOOST_BUILD_DIR}")

    if (NOT EXISTS "${BOOST_BUILD_DIR}/boost_unpack_complete")
        message(STATUS "Unpacking boost library source")
        execute_process(
            COMMAND ${CMAKE_COMMAND} -E tar xjf "${THIS_REDIST_DIR}/${BOOST_SOURCE_PREFIX}.tar.bz2"
            WORKING_DIRECTORY "${BOOST_BUILD_DIR}"
            RESULT_VARIABLE TMP_RESULT)

        if (TMP_RESULT)
            message (FATAL_ERROR "Failed to unpack boost library ${THIS_BOOST_VERSION}")
        endif ()
        execute_process(
            COMMAND ${CMAKE_COMMAND} -E touch "${BOOST_BUILD_DIR}/boost_unpack_complete")
    endif ()

    set (BOOST_BOOTSTRAP sh "bootstrap.sh")
    if (WIN32)
        set (BOOST_BOOTSTRAP "bootstrap.bat")
    endif ()

    # boost compile works in windows, but we aren't going to link anyway so we're
    # skipping to save time:
    #
    if (NOT WIN32)

        set (BJAM_OPTIONS "")

        if (${aarch64} STREQUAL "ON")
             set (BJAM_OPTIONS ${BJAM_OPTIONS} "architecture=arm")     
        endif ()

        set (UCONFIG "${BOOST_SRC_DIR}/tools/build/src/user-config.jam")
        if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
            file (WRITE "${UCONFIG}" "using gcc : : \"${CMAKE_CXX_COMPILER}\" ;\n")
        elseif (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
            file (WRITE "${UCONFIG}" "using clang : : \"${CMAKE_CXX_COMPILER}\" ;\n")
            set (BJAM_OPTIONS ${BJAM_OPTIONS} "define=_GLIBCXX_USE_CXX11_ABI=0")
            set (BJAM_OPTIONS ${BJAM_OPTIONS} "toolset=clang")
        endif ()

        message(STATUS "Configuring boost library")
        execute_process(
            COMMAND ${BOOST_BOOTSTRAP}
            WORKING_DIRECTORY "${BOOST_SRC_DIR}"
            RESULT_VARIABLE TMP_RESULT
            OUTPUT_QUIET)

        if (TMP_RESULT)
            message (FATAL_ERROR "Failed to configure boost library ${THIS_BOOST_VERSION}")
        endif ()

        foreach (BOOST_LIBRARY ${THIS_BOOST_BUILD_COMPONENTS})
            list (APPEND BOOST_BJAM_LIBRARY_SELECTION "--with-${BOOST_LIBRARY}")
        endforeach ()

        # Include full path for bjam so that we don't depend on cwd in build user's PATH
        set (BOOST_BJAM "${BOOST_SRC_DIR}/bjam")

        set (BJAM_OPTIONS ${BJAM_OPTIONS} "link=static")
        if (WIN32)
            if (MSVC)
                math (EXPR VS_VERSION "(${MSVC_VERSION}/100) - 6")
                set(TOOLSET "toolset=msvc-${VS_VERSION}.0")
            endif ()

            # windows bjam ignores single/static so don't even try:
            set (BJAM_OPTIONS ${BJAM_OPTIONS} "${TOOLSET}")
        else ()
            set (BJAM_OPTIONS ${BJAM_OPTIONS} "threading=single")
        endif ()

        set (BOOST_BUILD_CMD ${BOOST_BJAM} --prefix=${BOOST_INSTALL_DIR} ${BOOST_BJAM_LIBRARY_SELECTION} -j${CMAKE_PARALLEL} --libdir=${BOOST_INSTALL_DIR}/lib ${BJAM_OPTIONS} install)
        set (BOOST_BUILD_ERROR_LOG "${BOOST_INSTALL_DIR}/boost.build.error.txt")
        message(STATUS "Building boost library")
        execute_process(
            COMMAND ${BOOST_BUILD_CMD}
            WORKING_DIRECTORY "${BOOST_SRC_DIR}"
            RESULT_VARIABLE TMP_RESULT
            OUTPUT_FILE ${BOOST_BUILD_ERROR_LOG}
            ERROR_FILE ${BOOST_BUILD_ERROR_LOG})

        if (TMP_RESULT)
            string (REPLACE ";" " " BOOST_PRETTY_BUILD_CMD "${BOOST_BUILD_CMD}")
            message (STATUS "Boost build command: 'cd ${BOOST_SRC_DIR} && ${BOOST_PRETTY_BUILD_CMD}'")
            message (FATAL_ERROR "Failed to build boost library ${THIS_BOOST_VERSION}. See build log: ${BOOST_BUILD_ERROR_LOG}")
        endif ()

        message (STATUS "Successfuly built boost ${THIS_BOOST_VERSION}")
        execute_process(
            COMMAND ${CMAKE_COMMAND} -E touch "${BOOST_INSTALL_DIR}/boost_install_complete"
            WORKING_DIRECTORY "${BOOST_SRC_DIR}")
    else ()
        # for the time being on windows the goal is only to enable dev and library level compilation, not linking, so we don't need the libraries
        # do a quick copy instead to get the headers in place only:
        message(STATUS "WIN32 - skipping to header only boost installation for non-linked development")
        makedir("${BOOST_INSTALL_DIR}/include")
        makedir("${BOOST_INSTALL_DIR}/lib")
        file(COPY "${BOOST_SRC_DIR}/boost"
            DESTINATION "${BOOST_INSTALL_DIR}/include")
    endif ()

    set (BOOST_ROOT "${BOOST_BOOTSTRAP_INSTALL_DIR}")
    find_package(Boost ${THIS_BOOST_VERSION} COMPONENTS ${THIS_BOOST_COMPONENTS})
endif ()

foreach(COMPONENT ${THIS_BOOST_COMPONENTS})
    STRING(TOUPPER ${COMPONENT} UPPERCOMPONENT)
    set(HAVE_LIBBOOST_${UPPERCOMPONENT} ${Boost_${UPPERCOMPONENT}_FOUND})
endforeach()
