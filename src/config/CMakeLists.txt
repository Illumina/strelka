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

# starling VQSR config files:
file (GLOB JSON_CONFIG_LIST [a-zA-Z0-9]*.json)

foreach(CONFIG_FILE ${JSON_CONFIG_LIST})
    install(FILES ${CONFIG_FILE} DESTINATION ${THIS_CONFIGDIR})
endforeach()

# strelka user config file templates:
file (GLOB STRELKA_CONFIG_LIST strelka_[a-zA-Z0-9]*.ini)

foreach(CONFIG_FILE ${STRELKA_CONFIG_LIST})
    install(FILES ${CONFIG_FILE} DESTINATION ${THIS_CONFIGDIR})
endforeach()
