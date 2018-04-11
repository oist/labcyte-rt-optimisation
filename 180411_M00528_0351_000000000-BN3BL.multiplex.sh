# Command used to create 180411_M00528_0351_000000000-BN3BL.multiplex.txt

# Indexes taken from 180411_M00528_0351_000000000-BN3BL.SampleSheet.csv

printf "samplename\tgroup\tbarcode\tindex\n" > 180411_M00528_0351_000000000-BN3BL.multiplex.txt
for index in $(grep library_ 180411_M00528_0351_000000000-BN3BL.SampleSheet.csv | cut -f6 -d,)
do
  for barcode in $(sed 1d multiplex_model.txt | cut -f3)
  do
    printf "${barcode}_${index}\t${index}\t${barcode}\t${index}\n" >> 180411_M00528_0351_000000000-BN3BL.multiplex.txt
  done
done
