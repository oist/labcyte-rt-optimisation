#!/bin/sh

printf "samplename\tgroup\tbarcode\tindex\n" > 180517_M00528_0364_000000000-BRGK6.multiplex.txt

tomultiplex() {
  while read index barcode
  do
    printf "${barcode}_${index}\t$index\t$barcode\t$index\n"
  done >> 180517_M00528_0364_000000000-BRGK6.multiplex.txt
}

sed 1d plate6a.txt | cut -f 12,16 | tomultiplex
sed 1d plate6b.txt | cut -f 12,16 | tomultiplex
sed 1d plate6c.txt | cut -f 12,16 | tomultiplex
sed 1d plate6d.txt | cut -f 12,16 | tomultiplex
