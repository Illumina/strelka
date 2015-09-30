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
## file CMakeLists.txt
##
################################################################################

include ("${THIS_GLOBALS_CMAKE}")

set (THIS_STRELKADEMODIR "${THIS_DEMODIR}/strelka")
install(FILES "README.txt" DESTINATION ${THIS_STRELKADEMODIR})
install(DIRECTORY "expectedResults" DESTINATION ${THIS_STRELKADEMODIR})
install(DIRECTORY "expectedResultsStarlingMito" DESTINATION ${THIS_STRELKADEMODIR})

install(FILES "runStrelkaWorkflowDemo.bash"  DESTINATION ${THIS_BINDIR} PERMISSIONS ${THIS_EXECUTABLE_PERMISSIONS})
install(FILES "runStarlingMitoWorkflowDemo.bash"  DESTINATION ${THIS_BINDIR} PERMISSIONS ${THIS_EXECUTABLE_PERMISSIONS})
install(FILES "runStarlingWorkflowDemo.bash"  DESTINATION ${THIS_BINDIR} PERMISSIONS ${THIS_EXECUTABLE_PERMISSIONS})
install(FILES "runDenovoVariantWorkflowDemo.bash"  DESTINATION ${THIS_BINDIR} PERMISSIONS ${THIS_EXECUTABLE_PERMISSIONS})

add_subdirectory (data)

