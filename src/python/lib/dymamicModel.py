#!/illumina/development/pgtools/pgt-virtualenv/bin/python
########################################
# Tool for estimating indel error paramters from a bam file 
#
# Author: Morten Kallberg
# Date: March 3rd, 2014
#
########################################

parser = OptionParser(usage='%prog ')
parser.add_option('-b', '--inputBam', dest="inputBam", help='A comma delimited list of input bams')
parser.add_option('-s', '--startFromFile', dest="startFromFile", help='Json file that already has sample stats.')
parser.add_option('-p', '--inputSample', dest="inputSample", help='Sample used, in example NA12878', default='NA12878')
parser.add_option('-d', '--', dest="distributionAsTruth", help='Estimate truth site from frequency rather than PG', default=True)
parser.add_option('-o', '--output', dest="output", help='Output path')
parser.add_option('-r', '--reference', dest="reference",help='',default='Path to reference genomes')
parser.add_option('-v', '--verbose', dest="verbose",help='Output debugging information',default=1)

## the following options cause the stored value to be False; the option name/description is inverted from the variable they cause to be set
parser.add_option('--useStrand', dest="ignoreStrand",help='Use strand of alignment to determine homopolymer base',default=False,action='store_false')
parser.add_option('--coverPerAlignment', dest="naiveProject",help='Use CIGAR to determine overlapped reference bases',default=False,action='store_false')
parser.add_option('--requireSpanning', dest="countUnspanned",help='Consider a homopolymer covered by a read only if the read aligns beyond the hpol region on both sides',default=True,action='store_false')
parser.add_option('--countReads', dest="countEveryBase",help='Count reads overlapping a homopolymer, not bases (i.e. count the read once)',default=True,action='store_false')


opt, args = parser.parse_args()
# setup logger
logging.basicConfig(format='%(asctime)s %(levelname)-8s %(message)s', level=logging.INFO)
logging.info("This is the indel error fitting script, we are in business")


