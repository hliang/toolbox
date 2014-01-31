#!/bin/bash
set -o verbose

# read files
read1="${HOME}/projects/ps/raw/PS001_S1_L001_R1_001.fastq.gz"
read2="${HOME}/projects/ps/raw/PS001_S1_L001_R2_001.fastq.gz"
# ref seq file
ref="${HOME}/projects/ps/B728a.fasta"

# make link files poiting to original files
ln -s $read1 ./raw.r1.fq
ln -s $read2 ./raw.r2.fq
ln -s $ref   ./ref.fa

# exit 0;

# trim 3'-end low quality, keep only first 150bp.
# do not trim 5'-end, because we want to estimate the length of insert.
seqtk trimfq -b 0 -e 150 raw.r1.fq > qc.r1.fq
seqtk trimfq -b 0 -e 150 raw.r2.fq > qc.r2.fq

# mapping using bwa-mem 
#         -t 4 : 4 threads
# index
~/tools/bwa-0.7.5a/bwa index ref.fa 2>bwa.log
if [ ! $? ]; then
    echo "ERROR: there is an issue when indexing reference"
    exit 
fi
# align
~/tools/bwa-0.7.5a/bwa mem -t 4 ref.fa qc.r1.fq qc.r2.fq > all.bwa-pe.sam 2>>bwa.log
if [ ! $? ]; then
    echo "ERROR: there is an issue when mapping reads to ref"
    exit 
fi

# get tlen
cut -f 1-9 all.bwa-pe.sam | awk '$7=="=" && $9!=0 {print $9}' > tlen


echo "INFO: Generating Rscript..."
#######################################################
echo "########### R script for plotting tlen ##################" > insert.R
echo "# read table" >> insert.R
echo "tlen = read.table(\"./tlen\")" >> insert.R
echo "tlen = tlen\$V1" >> insert.R
echo "tlen[tlen>1000] = 1000" >> insert.R
echo "tlen = tlen[tlen > 0 & tlen <= 1000]" >> insert.R
echo "#" >> insert.R
echo "# sample" >> insert.R
echo "tlen = sample(tlen, 10000)" >> insert.R
echo "" >> insert.R
echo "# save plots in pdf" >> insert.R
echo "pdf(\"insertsize.pdf\")" >> insert.R
echo "" >> insert.R
echo "# save the graphics parameters currently used" >> insert.R
echo "def.par <- par(no.readonly = TRUE)" >> insert.R
echo " " >> insert.R
echo "layout(mat = matrix(c(2,1),2,1, byrow=TRUE), height = c(2,1.3))" >> insert.R
echo "par(mar=c(5.1, 4.1, 2.1, 2.1))" >> insert.R
echo " " >> insert.R
echo "leg1 <- paste(\"mean = \", round(mean(tlen), digits = 1))" >> insert.R
echo "leg2 <- paste(\"s.d. = \", round(sd(tlen),digits = 1))" >> insert.R
echo "leg3 <- paste(\"median = \", round(median(tlen), digits = 1))" >> insert.R
echo " " >> insert.R
echo "# make plots" >> insert.R
echo "boxplot(tlen, horizontal = TRUE, ylim=c(0,1000), col=\"gray\", xlab=\"insert size (bp)\")" >> insert.R
echo " " >> insert.R
echo "# lines(density(tlen, kernel=\"gaussian\"), col=\"black\")" >> insert.R
echo "h = hist(tlen, breaks=40, prob=F, xlim=c(0,1000), main=\"\", xlab=\"\", xaxt=\"n\", col=\"gray\", border=F)" >> insert.R
echo "axis(side = 1, at = seq(0, 1000, by=100), col=1)" >> insert.R
echo "legend(x = \"topright\", c(leg1,leg2,leg3), bty = \"n\", text.col=\"black\")" >> insert.R
echo "" >> insert.R
echo "# mimicking bioanalyzer profile." >> insert.R
echo "points(h\$mids+120, h\$counts * (h\$mids+120) * max(h\$counts) / max(h\$counts * (h\$mids+120)), type=\"l\", col=\"blue\")" >> insert.R
echo "legend(x = \"right\", \"mass = (insert + adapter) * counts\n(mimicking bioAnalyzer)\", lty=1, lwd=2, col=\"blue\", bty = \"n\", text.col=\"blue\")" >> insert.R
echo " " >> insert.R
echo "# reset the graphics display to default" >> insert.R
echo "par(def.par)" >> insert.R
echo "" >> insert.R
echo "dev.off()" >> insert.R
##############################################
echo "INFO: Done"

cat insert.R
Rscript --verbose --vanilla insert.R

if [ ! $? ]; then
    echo "ERROR: there is an issue when running R scirpt"
    exit 
fi

exit
