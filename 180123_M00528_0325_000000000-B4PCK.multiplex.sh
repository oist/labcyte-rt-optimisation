# Command used to create 180123_M00528_0325_000000000-B4PCK.multiplex.txt

# Indexes taken from /sequencedata/MiSeq/180123_M00528_0325_000000000-B4PCK/SampleSheet.csv

printf "samplename\tgroup\tbarcode\tindex\n" > 180123_M00528_0325_000000000-B4PCK.multiplex.txt
for index in GTAGAGGA GCTCATGA ATCTCAGG ACTCGCTA GGAGCTAC CTCTCTAC CGGAGCCT TACGCTGC ATGCGCAG TAGCGCTC ACTGAGCG CCTAAGAC CGATCAGT TGCAGCTA TCGACGTC
do
  for barcode in $(sed 1d multiplex.txt | cut -f3)
  do
    printf "${barcode}_${index}\t${index}\t${barcode}\t${index}\n" >> 180123_M00528_0325_000000000-B4PCK.multiplex.txt
  done
done
