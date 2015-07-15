mkdir result

../../bin/subread-buildindex -o ../small1 ../chr901.fa

../../bin/subjunc -i ../small1 -o result/junctions.sam -r data/junction-reads-A.fq -R data/junction-reads-B.fq

../../bin/subjunc -i ../small1 -o result/junctionsNfusions.bam --BAMoutput -r data/junction-reads-A.fq -R data/junction-reads-B.fq --allJunctions
