
# Overview

This directory contains very small datesets which can be used to
demonstrate/verify workflow configuration and run steps using a
demonstration script. To run each demonstration, execute the
corresponding script in the strelka workflow installation:

## Germline

    bash ${STRELKA_INSTALL_PATH}/bin/runStrelkaGermlineWorkflowDemo.bash


## Somatic

    bash ${STRELKA_INSTALL_PATH}/bin/runStrelkaSomaticWorkflowDemo.bash


This script will also check that the final vcfs produced match the
expected result.

This demonstation compares sequences from two normal (non-cancer)
samples -- using NA12892 as the "example normal" and NA12891 as the
"example tumor". For each sample the included bam files contain
sequences for only a small segment of chr20. The region has no
intended significance apart from providing several simulated somatic
snv and indel calls.

Note that this example uses default configurations, except that  the
`--exome` option has been selected, to indicate that this is not WGS
data and therefore standard filters based on locus depth relative to
average chromosome depth are not applied.
