# Charles Plessy, RIKEN
# https://creativecommons.org/publicdomain/zero/1.0/

# curl --silent https://gist.githubusercontent.com/charles-plessy/9dbc8bc98fb773bf71b6/raw | tee getAndParseGencode.bash | bash

# Comments starting with ##!! signal Linux/OSX portability traps

# Setup
# =====

# Default values are for human release 23.

RELEASE=${RELEASE-23}
GENOME=${GENOME-hg38}
BASEURL=${BASEURL-ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_$RELEASE}
GENCODE=${GENCODE-gencode.v$RELEASE.annotation}
ANNOT=${ANNOT-$GENCODE.bed}
JUNCTION_WINDOW=${JUNCTION_WINDOW-10}
PROMOTER_WINDOW=${PROMOTER_WINDOW-500}
SHA_ANNOT=${SHA_ANNOT-496f6f6b25740df218b7886d510c9523209e63cd}

# For mouse: download the script and then run:
# RELEASE=M8 GENOME=mm10 BASEURL=ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_$RELEASE GENCODE=gencode.v$RELEASE.annotation ANNOT=$GENCODE.bed JUNCTION_WINDOW=10 PROMOTER_WINDOW=500 SHA_ANNOT=2d90a3aa234f76ab6908e0bdfc5ed477db8aea5b bash getAndParseGencode.bash
# RELEASE=M1 GENOME=mm9  BASEURL=ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_$RELEASE GENCODE=gencode.v$RELEASE.annotation ANNOT=$GENCODE.bed JUNCTION_WINDOW=10 PROMOTER_WINDOW=500 SHA_ANNOT=7ba9a4ca84d7928271a7f6870550300c65bbe456 bash getAndParseGencode.bash

# Download
# ========

getAndCheck () {  ##!! note that function foo {} is a bashism
  curl --silent --continue-at - --remote-name $1  ##!! no wget on OSX
  echo "$2  $(basename $1)" | shasum --check  ##!! no md5sum on OSX
}
getAndCheck $BASEURL/gencode.v$RELEASE.annotation.gtf.gz $SHA_ANNOT

# Parser functions
# ================

# Filter with Perl expressions.

filterWithPerl () {
  zcat < $GENCODE.gtf.gz |  ##!! zcat foo.gz fails on OSX
    grep $'\t'"$1"$'\t' |  ##!! no grep -P on OSX
    perl -F'\t' -anE "
      (\$chr, \$source, \$feat, \$start, \$end, \$score, \$strand, \$frame, \$meta) = @F ;
      $2" |
    sort --key 1,1 --key 2,2n --key 3,3n --key 4,4 --key 6 --unique  ##!! no "version sort" on OSX
}

# Convert to TSS coordinates

bed2start () { perl -nE ' chomp ;
  ($chr, $start, $end, $name, $score, $strand) = (split "\t") ;
  $strand eq "+" and say join "\t", $chr, $start,  $start +1, $name, $score, $strand ;
  $strand eq "-" and say join "\t", $chr, $end -1, $end     , $name, $score, $strand ; '
}

# Feature names
# =============

# Genes

filterWithPerl gene '
  say join "\t", $chr, $start - 1, $end, $meta =~ /gene_name "(.*?)"/, 0, $strand
  ' | tee $GENCODE.genes.bed | bed2start > $GENCODE.genes-tss.bed

# Transcripts

filterWithPerl transcript '
  say join "\t", $chr, $start - 1, $end, join (",", $meta =~ /gene_name "(.*?)"/, $meta =~ /transcript_id "(.*?)"/), 0, $strand
  ' | tee $GENCODE.transcripts.bed | bed2start > $GENCODE.transcripts-tss.bed

# Annotation file
# ===============

# Gene structure

filterWithPerl gene '
  say join "\t", $chr, $start - 1, $end, "gene", 0, $strand;
  say join "\t", $chr, $start - 1, $end, join ("_", $meta =~ /gene_type "(.*?)"/, $meta =~ /gene_name "(.*?)"/), 0, $strand
  ' > $ANNOT

filterWithPerl exon '
  say join "\t", $chr, $start - 1, $end, "exon", 0, $strand' >> $ANNOT

filterWithPerl exon '
  say join "\t", $chr, $start - 1, $start, "boundary", 0, $strand;
  say join "\t", $chr, $end - 1,   $end,   "boundary", 0, $strand;
  ' | bedtools slop -b $JUNCTION_WINDOW -g /usr/share/bedtools/genomes/*.$GENOME.* >> $ANNOT

filterWithPerl transcript '
  next unless (/transcript_type "protein_coding"/ ||
               /transcript_type "processed_transcript"/ ||
               /transcript_type "lincRNA"/ ||
               /transcript_type "antisense"/ ||
               /transcript_type "processed_pseudogene"/ ||
               /transcript_type "unprocessed_pseudogene"/);
  $strand eq "+" and say join "\t", $chr, $start - 1, $start, "promoter", 0, $strand;
  $strand eq "-" and say join "\t", $chr, $end - 1,   $end,   "promoter", 0, $strand;
  ' |  bedtools slop -b $PROMOTER_WINDOW -g /usr/share/bedtools/genomes/*.$GENOME.* >> $ANNOT

TMPFILE=$(mktemp ./tmpXXXXXXXX)
sort --key 1,1 --key 2,2n --key 3,3n --key 4,4 --key 6 --unique $ANNOT > $TMPFILE
mv $TMPFILE $ANNOT
unset TMPFILE

# Checksum
# ========

echo '5fc86bf51a0b9f6a8ca67400965de6ca5c75f3a4  gencode.v23.annotation.bed'                 | shasum --check
echo '8f3477467284512afb35aa1311b194c4cf1985cf  gencode.v23.annotation.genes.bed'           | shasum --check
echo '61c9e2ef1bf8525fa5423d51ce32f499672b1843  gencode.v23.annotation.genes-tss.bed'       | shasum --check
echo '5a2823fae91ab07c8af5d53ce35f94ef2ec52470  gencode.v23.annotation.transcripts.bed'     | shasum --check
echo '83ca1c473c009365a737e3ba92bdaeab70d01634  gencode.v23.annotation.transcripts-tss.bed' | shasum --check

# for mouse:

# e779c2c7de0a5454570f21a38d22f0c2aacccfac  gencode.vM8.annotation.bed
# d6e36324d25f9604e9cf4fd3f554e1df90647b2e  gencode.vM8.annotation.genes.bed
# fe488a149fffcc30b80a26607fb558a0573a15c2  gencode.vM8.annotation.genes-tss.bed
# 46983f68778fdbb0cdc1938b11fd662292368b72  gencode.vM8.annotation.transcripts.bed
# 47fa967f8a9adc661b01ac2905e0ec6217921f6f  gencode.vM8.annotation.transcripts-tss.bed

# 8bea34df6385e361abf899cb3ba8952dc46c9171  gencode.vM1.annotation.bed
# 739a29e5092e96e1812a43698650c4333e33c160  gencode.vM1.annotation.genes.bed
# 554404527c4744d02f69ddf2548aba97a9b0e105  gencode.vM1.annotation.genes-tss.bed
# 1f865671ef7bd5d1f26a963e76ca3cdd4162ccaf  gencode.vM1.annotation.transcripts.bed
# 220aeee2bf94c9d04e8582f66c0441d45ba877b4  gencode.vM1.annotation.transcripts-tss.bed


