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
################################################################################

file (GLOB CONFIG_LIST [a-zA-Z0-9]*.json)

foreach(CONFIG_FILE ${CONFIG_LIST})
    install(FILES ${CONFIG_FILE} DESTINATION ${THIS_CONFIGDIR})
endforeach()
