mkdir -p result
echo
echo "================================================================================"
printf " FeatureCounts Corner Case Tests\n http://subread.sourceforge.net/\n"
echo "================================================================================"
echo



SH_CMD=bash
$SH_CMD data/compare.sh data/corner-INDEL.sam data/corner-INDEL.ora data/test-minimum.GTF "-p" "indel reads"
$SH_CMD data/compare.sh data/corner-JUNC.sam data/corner-JUNC.ora data/test-minimum.GTF "-p" "junction reads"
$SH_CMD data/compare.sh data/corner-ONEEND.sam data/corner-ONEEND.ora data/test-minimum.GTF "-p" "paired-end reads (fragment counting)"
$SH_CMD data/compare.sh data/corner-ONEEND.sam data/corner-ONEEND-BOTH.ora data/test-minimum.GTF "-p -B " "paired-end reads (fragment counting, both ends mapped)"
$SH_CMD data/compare.sh data/test-minimum.sam data/test-minimum-O.ora data/test-minimum.GTF "-p -O " "multi-overlapping reads"
$SH_CMD data/compare.sh data/test-minimum.sam data/test-minimum-FL.ora data/test-minimum.GTF "-p -f " "feature-level summarization" FL
$SH_CMD data/compare.sh data/test-minimum.sam data/test-minimum.ora data/test-minimum.GTF "-p " "gene-level summarization"
$SH_CMD data/compare.sh data/corner-NH.sam data/corner-NH.ora data/test-minimum.GTF "-p" "multi-mapping reads"
$SH_CMD data/compare.sh data/corner-NH.sam data/corner-NH-PM.ora data/test-minimum.GTF "-p --primary " "multi-mapping reads (primary only)"

# by default it is in GTF.
$SH_CMD data/compare.sh data/test-minimum.sam data/test-minimum.ora data/test-minimum.GTF "-p " "GTF format annotations"
$SH_CMD data/compare.sh data/test-minimum.sam data/test-minimum.ora data/test-minimum.SAF "-p -F SAF " "SAF format annotations"

# by default it is in SAM.
$SH_CMD data/compare.sh data/test-minimum.sam data/test-minimum.ora data/test-minimum.GTF "-p " "SAM format input"
$SH_CMD data/compare.sh data/test-minimum.bam data/test-minimum.ora data/test-minimum.GTF "-p -b " "BAM format input" 

# by default it is non-strand specific.
$SH_CMD data/compare.sh data/test-minimum.sam data/test-minimum.ora data/test-minimum.GTF "-p -s 0 " "unstranded read summarization"
$SH_CMD data/compare.sh data/test-minimum.sam data/test-minimum-STR.ora data/test-minimum.GTF "-p -s 1 " "stranded read summarization"
$SH_CMD data/compare.sh data/test-minimum.sam data/test-minimum-UNSTR.ora data/test-minimum.GTF "-p -s 2 " "reversely stranded read summarization"

$SH_CMD data/compare.sh data/corner-JUNC.sam data/corner-JUNC-ONLY.ora data/test-minimum.GTF "--countSplitAlignmentsOnly -O -f " "Junction reads only" FL

echo
