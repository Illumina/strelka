
This directory contains very small datesets which can be used to
demonstrate/verify strelka's configuration and run steps using a
demonstration script. To run the demonstration, run the demo
script in the strelka workflow installation: 

Strelka:

"bash $INSTALL_DIR/bin/runStrelkaWorkflowDemo.bash"

This script will also check that the final vcfs produced match the
expected result.

This demonstation compares seqeunces from two normal (non-cancer)
samples -- using NA12892 as the "example normal" and NA12891 as the
"example tumor". For each sample the included bam files contain
sequences for only a small segment of chr20. The region has no
intended significance apart from providing several simulated somatic
snv and indel calls.

Note that this example uses the default BWA configuration file, except
that 'isSkipDepthFilters' has to be set to 1. The depth filter
normally will remove all calls which occur at a depth greater than 3x
the chromosomal mean. It is meant to remove calls near centromeres and
telomeres when the whole genome is sequenced. This filter should be
turned off in any where regions of the genome have been selected or
enriched for, such as exome sequencing or artificial subsets of data
such as this demo.

Starling:

"bash $INSTALL_DIR/bin/runStarlingWorkflowDemo.bash"

Verifies germline variant calling.

"bash $INSTALL_DIR/bin/runStarlingMitoWorkflowDemo.bash"

Verifies the continuous VF variant calling in mitochondria.

