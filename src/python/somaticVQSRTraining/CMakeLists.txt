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

set(STRELKA_VQSR_DIR "${THIS_DATADIR}/somaticVQSRTraining")
set(STRELKA_VQSR_LIBDIR "${STRELKA_VQSR_DIR}/lib")

install(FILES "README.md" DESTINATION ${STRELKA_VQSR_DIR})

add_subdirectory (bin)
add_subdirectory (lib)

