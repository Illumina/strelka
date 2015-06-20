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

set(INSTALL_TO_DIR "${STRELKA_VQSR_LIBDIR}")

include("${THIS_GLOBALS_CMAKE}") # required for *_PERMISSIONS values

install(DIRECTORY "./" DESTINATION "${INSTALL_TO_DIR}" FILE_PERMISSIONS ${THIS_LIBRARY_PERMISSIONS} FILES_MATCHING PATTERN "*.py")

