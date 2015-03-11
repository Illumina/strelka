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

install(FILES ${JSON_CONFIG_LIST} DESTINATION ${THIS_CONFIGDIR})
