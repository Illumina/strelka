
# Overview

This directory contains small datesets which can be used to
demonstrate/verify workflow configuration and run steps using a
demonstration script.

## Available demos

### Germline

The germline demo runs a joint germline analysis on two samples, and verifies
that the analysis completes without errors. The two samples used for the
demo cover small segments from NA12892 and NA12891. The demo can be run as
follows:

    bash ${STRELKA_INSTALL_PATH}/bin/runStrelkaGermlineWorkflowDemo.bash


### Somatic

The somatic demo runs a tumor-normal subtraction. It verifies that the analysis
completes without errors and that the demo results match expectation. The normal
and tumor samples used for the demo cover small segments from NA12892 and NA12891
respectively. The demo can be run as follows:

    bash ${STRELKA_INSTALL_PATH}/bin/runStrelkaSomaticWorkflowDemo.bash


## Additional Details

All demo read alignments are made to the abstract chromosome 'demo20' which
covers chr20:862001-867000 from the hg19 assembly.

Note that these examples use default configurations, except:
1. The `--exome` option has been selected to indicate that this is not WGS
data. With this option, standard filters based on locus depth relative to
average chromosome depth are not applied.
2. The `--disableSequenceErrorEstimation` option has been selected for the
germline demo. The demo input is too small to dynamically estimate indel
errors. Strelka can automatically detect and disable this feature based on
input data quantity, but adding this flag prevents confusing warnings
from being printed in the demo run.

