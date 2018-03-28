# Command used to create 180326_M00528_0346_000000000-B4GJR.multiplex.txt

# Indexes taken from 180326_M00528_0346_000000000-B4GJR.SampleSheet.csv

printf "samplename\tgroup\tbarcode\tindex\n" > 180326_M00528_0346_000000000-B4GJR.multiplex.txt
for index in TAAGGCGA CGTACTAG AGGCAGAA TCCTGAGC GGACTCCT
do
  for barcode in $(sed 1d multiplex.txt | cut -f3)
  do
    printf "${barcode}_${index}\t${index}\t${barcode}\t${index}\n" >> 180326_M00528_0346_000000000-B4GJR.multiplex.txt
  done
done

for index in TAGGCATG CTCTCTAC CGAGGCTG AAGAGGCA GTAGAGGA
do
  for barcode in $(sed 1d multiplex.txt | cut -f3)
  do
    printf "${barcode}_${index}\t${index}\t${barcode}\t${index}\n" >> 180326_M00528_0346_000000000-B4GJR.multiplex.txt
  done
done
