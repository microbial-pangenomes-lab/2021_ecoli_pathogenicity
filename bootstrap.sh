#!/bin/bash

mkdir -p out
mkdir -p out/logs

echo -e "ID\tPath" > out/unitigs_input.tsv
cat /dev/null > out/mash_input.txt
cat /dev/null > out/panaroo_input.txt
cat /dev/null > out/annotate_input.txt
for strain in $(tail -n+2 data/GWAS_20201217_MEO.csv | awk -F';' '{print $3}');
do
  fasta="data/fastas/"$strain".fasta";
  gff="data/gffs/"$strain".gff";
  if [ -f "$fasta" ] && [ -f "$gff" ];
  then
    echo $strain;
    echo "data/fastas/"$strain".fasta" >> out/mash_input.txt;
    echo "data/gffs/"$strain".gff" >> out/panaroo_input.txt;
    echo -e $strain"\tdata/fastas/"$strain".fasta" >> out/unitigs_input.tsv;
    echo -e "data/fastas/"$strain".fasta\tdata/gffs/"$strain".gff\tdraft" >> out/annotate_input.txt;
  fi;
done

# fix issue with contig IDs (for mapping back)
mkdir -p data/fixed_fastas
for i in $(ls data/fastas/*.fasta);
do
  python workflow/scripts/fix_fastas.py $i data/gffs/$(basename $i .fasta).gff > data/fixed_fastas/$(basename $i);
  echo $i;
done

# prepare snpeff databases
snpEff build -c config/snpeff.config -dataDir ../data/references/ -gff3 ED1a
snpEff build -c config/snpeff.config -dataDir ../data/references/ -gff3 IAI1
snpEff build -c config/snpeff.config -dataDir ../data/references/ -gff3 IAI39
snpEff build -c config/snpeff.config -dataDir ../data/references/ -gff3 K-12
snpEff build -c config/snpeff.config -dataDir ../data/references/ -gff3 UMN026
snpEff build -c config/snpeff.config -dataDir ../data/references/ -gff3 UTI89
snpEff build -c config/snpeff.config -dataDir ../data/references/ -gff3 CFT073
snpEff build -c config/snpeff.config -dataDir ../data/references/ -gff3 536
