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

# Detect and run python unit tests
#
# All files in this directory following pattern "test_*.py" will be treated as a unittest module.

# infer list of tests modules from files in this directory
file(GLOB PYTHON_TEST_FILES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} "${CMAKE_CURRENT_SOURCE_DIR}/test_*.py")
set(PYTHON_TEST_MODULES "")
foreach(PYTHON_TEST_FILE ${PYTHON_TEST_FILES})
    get_filename_component(PYTHON_TEST_MODULE_NAME ${PYTHON_TEST_FILE} NAME_WE)
    list(APPEND PYTHON_TEST_MODULES ${PYTHON_TEST_MODULE_NAME})
endforeach()

# clean staging area for configured py files:
set(PYSTAGE_TESTS_DIR "${PYSTAGE_DIR}/tests")

# configure test files
file(RELATIVE_PATH THIS_RELATIVE_PYTHON_LIBDIR "${PYSTAGE_TESTS_DIR}" "${PYSTAGE_DIR}")
foreach(PYTHON_TEST_FILE ${PYTHON_TEST_FILES})
    configure_file("${CMAKE_CURRENT_SOURCE_DIR}/${PYTHON_TEST_FILE}" "${PYSTAGE_TESTS_DIR}/${PYTHON_TEST_FILE}" @ONLY)
endforeach()

if (PYTHONINTERP_FOUND)
    set(TEST_TARGET_NAME "${THIS_PROJECT_NAME}_unit_test_python_lib")

    add_test(NAME ${TEST_TARGET_NAME}
             COMMAND ${PYTHON_EXECUTABLE} -m unittest -q ${PYTHON_TEST_MODULES}
             WORKING_DIRECTORY ${PYSTAGE_TESTS_DIR})
endif ()
