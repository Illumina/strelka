#
# Manta
# Copyright (c) 2013-2014 Illumina, Inc.
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
## Configuration file for the markdown subfolder
##
## author Roman Petrovski
##
################################################################################

##
## build the documentation when available
##

macro(process_user_guide name)
    set (UG "${name}.md")
    set (UL "${CMAKE_CURRENT_BINARY_DIR}/${UP}.html")

    configure_files (${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR} ${UG})

    add_custom_command(
        OUTPUT ${UL}
        COMMAND ${HOEDOWN_PROG} ${CMAKE_CURRENT_BINARY_DIR}/${SGUG} > ${UL}
        DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/${UG}
        DEPENDS STARKA_MARKDOWN
        COMMENT "Converting ${name} to html")

    add_custom_target(STARKA_USERGUIDE ALL DEPENDS ${UL} ${UL})

    install(FILES ${UL} DESTINATION ${THIS_DOCDIR}/html)
    install(FILES ${UL} DESTINATION ${THIS_DOCDIR}/html)
endmacro()

process_user_guide("starlingUserGuide")
process_user_guide("starkaUserGuide")

install(FILES "userGuide.css" DESTINATION ${THIS_DOCDIR}/html)
install(FILES "tableOfContents.js" DESTINATION ${THIS_DOCDIR}/html)

