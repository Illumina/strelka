#!/usr/bin/env bash
#
# build latest starka and run a small test on single-sample and tumor-normal
#
# As currently setup, this only runs on sd-qmaster
#

set -o nounset
set -o pipefail
#set -o xtrace

# starka test constants (git checkout not finished yet)
#
starka_git_url=ussd-git.illumina.com:OptimusPrime/starka
starka_version=master #c82a82731e784700df1345dd64cc88360ad83a83
workspace_dir=$(pwd)/workspace


# starling test constants
#
starling_expected_time=53
starling_results_dir=$workspace_dir/starling_results
starling_output_name=Mother_S1_chr15.84M-86M.raw.genome.vcf
starling_expected_dir=/bioinfoSD/csaunders/proj/starka/test/expected_results/starling
starling_result=$starling_results_dir/$starling_output_name
test_data_dir=/bioinfoSD/csaunders/proj/starka/test/data
sample1_bam=$test_data_dir/Mother_S1.bam
test_reference=/illumina/scratch/iGenomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa


# strelka test constants
#
strelka_expected_time=79
strelka_results_dir=$workspace_dir/strelka_results
strelka_expected_dir=/bioinfoSD/csaunders/proj/starka/test/expected_results/strelka
sample2_bam=$test_data_dir/Father_S1.bam


#
# utility functions:
#
rel2abs() {
    ( cd $1 && pwd -P )
}


error() {
    echo ERROR: $@ 1>&2
    exit 1
}


check_time() {
    label=$1
    elapsed=$2
    expect=$3
    printf "%s elapsed_time: %i expected time: %i\n" $label $elapsed $expect 1>&2

    # how much do we allow timings to vary relative to expect
    tol=20

    max=$((($expect*(100+$tol))/100))
    min=$((($expect*(100-$tol))/100))

    if [ $elapsed -gt $max ]; then
        error $(printf "%s elapsed_time %i exceeded maximum time: %i\n" $label $elapsed $max)
    fi
    if [ $elapsed -lt $min ]; then
        error $(printf "%s elapsed_time %i less than minimum time: %i. Could be good... but suspiciious!!!\n" $label $elapsed $min)
    fi
}


check_exit() {
    label=$1
    exit_val=$2
    if [ 0 -ne $exit_val ]; then
        error "$label test failed. Exit val: $exit_val"
    fi
}


filter_variable_metadata() {
    awk '!/^##(fileDate|source_version|startTime|reference|cmdline)/' |\
    awk '$2!~/(START_TIME|CMDLINE|PROGRAM_VERSION)/'
}


diff_dirs() {
    label=$1
    result_dir=$2
    expect_dir=$3

    echo 1>&2
    echo "**** Starting comparison to expected results." 1>&2

    for f in $(ls $expect_dir); do
        if ! [ -f $result_dir/$f ]; then
            cat<<END 1>&2

ERROR: Found difference between current and expected $label results.
       Expected file: $expect_dir/$f
       Missing results file: $result_dir/$f

END
            exit 1
        fi

        diff <(filter_variable_metadata < $expect_dir/$f) <(filter_variable_metadata < $result_dir/$f)
        if [ $? -ne 0 ]; then
            cat<<END 1>&2

ERROR: Found difference between current and expected $label results.
       Expected file: $expect_dir/$f
       Results file: $result_dir/$f

END
            exit 1
        fi
    done

    echo 1>&2
    echo "**** No differences between expected and computed results." 1>&2
    echo "**** Test successfully completed" 1>&2
    echo 1>&2
}


script_dir=$(rel2abs $(dirname $0))
get_depth=$script_dir/util/getBamAvgChromDepth.pl



if ! [ -f $test_reference ]; then
    error "Can't find reference: $test_reference"
fi

#
# 1. put all the mess in one place (can make this a temp directory later:
#

mkdir -p $workspace_dir
cd $workspace_dir


#
# X. create starka tarball and compile (in progress):
#

starka_tarball_key=starka-master
starka_tarball_dir=$starka_tarball_key
starka_tarball_name=$starka_tarball_dir.tar.gz

make_starka_tarball() {
    cd $workspace_dir
    if [ -d $starka_tarball_dir ]; then rm -rf $starka_tarball_dir; fi
    git clone --recursive $starka_git_url

    (
        cd starka
        git checkout $starka_version
        cd scratch
        ./make_release_tarball.bash $starka_tarball_key
        mv $starka_tarball_name $workspace_dir
    )
}

make_starka_tarball
tar -xzf $starka_tarball_name
cd $starka_tarball_dir
make -j4
cd $workspace_dir

starka_bin=$workspace_dir/$starka_tarball_dir/bin


if ! [ -f $starka_bin/starling2 ]; then
    error "Can't find starling executable" 
fi
if ! [ -f $starka_bin/strelka2 ]; then
    error "Can't find strelka executable" 
fi

#
# 3. test starling
#
big_start=84000000
big_end=86000000
small_start=84990000
small_end=85010000

mkdir -p $starling_results_dir

depth_file=$workspace_dir/starling.bam.depth

printf "Starting Starling test\n" 1>&2
start_time=$(date +%s)

starling_command() {
    begin=$1
    end=$2
    echo $starka_bin/starling2 \
--gvcf-min-gqx 30 --gvcf-max-snv-strand-bias 10 --gvcf-max-indel-ref-repeat 8 -min-qscore 17 \
-min-vexp 0.25 -max-window-mismatch 2 20 -max-indel-size 50 -genome-size 2861343702 \
-clobber -min-single-align-score 20 -min-paired-align-score 20 -bsnp-ssd-no-mismatch 0.35 -bsnp-ssd-one-mismatch 0.6 \
--gvcf-file $starling_result \
--chrom-depth-file $depth_file \
-bam-file $sample1_bam \
-samtools-reference $test_reference \
-realigned-read-file $workspace_dir/starling.realigned.bam \
-bam-seq-name chr15 -report-range-begin $begin -report-range-end $end
}

starling_command_big_interval() {
    starling_command $big_start $big_end
}

starling_command_small_interval() {
    starling_command $small_start $small_end
}

$get_depth $sample1_bam >| $depth_file

$(starling_command_big_interval) >| $workspace_dir/starling.log 2>| $workspace_dir/starling.log2

#-bam-seq-name chr15 -report-range-begin 85046000 -report-range-end 85047000 >| log 2>| log2
check_exit Starling $?

end_time=$(date +%s)
starling_elapsed_time=$(($end_time - $start_time))

check_time Starling $starling_elapsed_time $starling_expected_time

diff_dirs Starling $starling_results_dir $starling_expected_dir


#
# 4. test strelka
#

printf "Starting Strelka test\n" 1>&2

mkdir -p $strelka_results_dir

start_time=$(date +%s)

strelka_command() {
    begin=$1
    end=$2
    echo $starka_bin/strelka2 \
-clobber -filter-unanchored -min-paired-align-score 20 \
-min-single-align-score 10 -min-qscore 0 \
-max-window-mismatch 3 20 -print-used-allele-counts -max-indel-size 50 -indel-nonsite-match-prob 0.5 \
--min-contig-open-end-support 35 --somatic-snv-rate 0.000001 --shared-site-error-rate 0.0000005 \
--shared-site-error-strand-bias-fraction 0.5 --somatic-indel-rate 0.000001 --shared-indel-error-rate 0.000001 \
--tier2-min-single-align-score 5 --tier2-min-paired-align-score 5 --tier2-single-align-score-rescue-mode \
--tier2-mismatch-density-filter-count 10 --tier2-no-filter-unanchored \
--tier2-indel-nonsite-match-prob 0.25 --tier2-include-singleton --tier2-include-anomalous \
--somatic-snv-file $strelka_results_dir/somatic.snvs.unfiltered.vcf \
--somatic-indel-file $strelka_results_dir/somatic.indels.unfiltered.vcf \
--variant-window-flank-file 50 $strelka_results_dir/somatic.indels.unfiltered.vcf.window \
--max-input-depth 10000 -genome-size 2861343702 \
-bam-file $sample1_bam \
--tumor-bam-file $sample2_bam \
-samtools-reference $test_reference \
-realigned-read-file $workspace_dir/strelka.realigned.normal.bam \
--tumor-realigned-read-file $workspace_dir/strelka.realigned.normal.bam \
-bam-seq-name chr15 -report-range-begin $begin -report-range-end $end
}

strelka_command_big_interval() {
    strelka_command $big_start $big_end
}

strelka_command_small_interval() {
    strelka_command $small_start $small_end
}


$(strelka_command_big_interval) >| $workspace_dir/strelka.log 2>| $workspace_dir/strelka.log2

check_exit Strelka $?

end_time=$(date +%s)
strelka_elapsed_time=$(($end_time - $start_time))

check_time Strelka $strelka_elapsed_time $strelka_expected_time

diff_dirs Strelka $strelka_results_dir $strelka_expected_dir


# add valgrind tests:
valgrind_prefix() {
    echo valgrind --error-exitcode=1 --tool=memcheck 
}

printf "Starting Starling valgrind test\n" 1>&2
$(valgrind_prefix) $(starling_command_small_interval) >| $workspace_dir/starling.valgrind.log 2>| $workspace_dir/starling.valgrind.log2
check_exit "Starling valgrind" $?

printf "Starting Strelka valgrind test\n" 1>&2
$(valgrind_prefix) $(strelka_command_small_interval) >| $workspace_dir/strelka.valgrind.log 2>| $workspace_dir/strelka.valgrind.log2
check_exit "Strelka valgrind" $?
