# Command used to create 180215_M00528_0330_000000000-B4GPD.multiplex.txt

# Indexes taken from /sequencedata/MiSeq/180215_M00528_0330_000000000-B4GPD/SampleSheet.csv

printf "samplename\tgroup\tbarcode\tindex\n" > 180215_M00528_0330_000000000-B4GPD.multiplex.txt
for index in TAAGGCGA CGTACTAG AGGCAGAA TCCTGAGC GGACTCCT TAGGCATG CTCTCTAC CGAGGCTG AAGAGGCA GTAGAGGA TAGCGCTC ACTGAGCG CCTAAGAC CGATCAGT TGCAGCTA
do
  for barcode in $(sed 1d multiplex.txt | cut -f3)
  do
    printf "${barcode}_${index}\t${index}\t${barcode}\t${index}\n" >> 180215_M00528_0330_000000000-B4GPD.multiplex.txt
  done
done
