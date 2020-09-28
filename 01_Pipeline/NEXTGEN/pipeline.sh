#!/bin/bash

# trimmomatic execution with $1 = MINLEN / $2 = ID / $3 = basefile
function trim() {
  mkdir -p trimmed/$1
  java -jar /bio-bin/Trimmomatic-0.39/trimmomatic-0.39.jar PE \
    -threads 4 \
    -basein $3 \
    -baseout trimmed/$1/$2_trimmed_Q28_L$1.fastq.gz \
    ILLUMINACLIP:/bio-bin/Trimmomatic-0.39/adapters/AllAdapters.fas:2:30:10 \
    LEADING:3 \
    TRAILING:3 \
    SLIDINGWINDOW:4:28 \
    MINLEN:$1

  cat trimmed/$1/$2_trimmed_Q28_L$1_1*.fastq.gz > trimmed/$1/$2_trimmed_Q28_L$1_forward.fastq.gz
  cat trimmed/$1/$2_trimmed_Q28_L$1_2*.fastq.gz > trimmed/$1/$2_trimmed_Q28_L$1_reverse.fastq.gz
  cat trimmed/$1/$2_trimmed_Q28_L$1_forward.fastq.gz trimmed/$1/$2_trimmed_Q28_L$1_reverse.fastq.gz > trimmed/$1/$2_trimmed_Q28_L$1_total.fastq.gz
}

# mitobim execution with $1 = MINLEN / $2 = ID
function mitobim() {
  mkdir -p mitobim/$1/$2
  mkdir -p mitobim/$1/logs
  pushd mitobim/$1/$2
  /bio-bin/MITObim/MITObim.pl \
    -start 1 \
    -end 30 \
    -sample $2 \
    -ref lacerta \
    -readpool ../../../trimmed/$1/$2_trimmed_Q28_L$1_total.fastq.gz \
    --clean \
    --quick ../../../../Daten/references/boehme_AM176577.1_lv.fas > ../logs/$2.log
  popd
}

# bowtie2 execution with
function bowtie() {
  mkdir -p bowtie/$1/$2
  pushd bowtie/$1/$2
  bowtie2 \
    -x ../../db/db_AM176577_1 \
    -1 ../../../trimmed/$1/$2_trimmed_Q28_L$1_1P.fastq.gz \
    -2 ../../../trimmed/$1/$2_trimmed_Q28_L$1_2P.fastq.gz \
    -S $2.sam \
    --no-unal

  samtools view -Sb $2.sam > $2.bam
  samtools sort $2.bam > $2.sorted.bam
  samtools index $2.sorted.bam
  popd
}

# print number of N nucleotids contained or '-' for failed input
function print_count_n_nucleotids() {
  if [ -e "$1" ]
  then
    printf $( grep -o 'N' $1 | wc -l )
  else
    printf '-'
  fi
}

# write alignment: $1 = mitobim result / $2 = id / $3 = extention
function alignment() {
  if [ -e "$1" ]
  then
    echo ">$2_$3" >> alignment.fas
    tail -n +2 $1 >> alignment.fas
  fi
}

# build bowtie index of reference sequence
mkdir -p bowtie/db
pushd bowtie/db
  bowtie2-build -f ../../../Daten/references/boehme_AM176577.1_lv.fas db_AM176577_1
popd

# run with several trim ranges
trim_count=(0 12 25 39)
for i in "${trim_count[@]}"; do
  for file in $(ls ../Daten/Sequencing_07.06.19/190607_M02279_0447_000000000-CDLF5_NS_Lacerta/FastQ/S*R1*fastq.gz); do
    sample_r1=$(basename -- ${file})
    id=$( echo $sample_r1 | cut -d'_' -f 1 )
    trim $i $id $file
    mitobim $i $id
    bowtie $i $id $file
  done;
done;

# show stats
echo 'id, R1, R2, 39, 25, 12, 0'
for file in $(ls ../Daten/Sequencing_07.06.19/190607_M02279_0447_000000000-CDLF5_NS_Lacerta/FastQ/S*R1*fastq.gz); do
  sample_r1=$(basename -- ${file})
  id=$( echo $sample_r1 | cut -d'_' -f 1 )
  assembly39=$( grep "Final assembly result will be written to file:" mitobim/39/logs/$id.log | grep -Eo "[^ ]+$" )
  assembly25=$( grep "Final assembly result will be written to file:" mitobim/25/logs/$id.log | grep -Eo "[^ ]+$" )
  assembly12=$( grep "Final assembly result will be written to file:" mitobim/12/logs/$id.log | grep -Eo "[^ ]+$" )
  assembly0=$( grep "Final assembly result will be written to file:" mitobim/0/logs/$id.log | grep -Eo "[^ ]+$" )
  printf $id
  printf ', '
  print_count_n_nucleotids $assembly39
  printf ', '
  print_count_n_nucleotids $assembly25
  printf ', '
  print_count_n_nucleotids $assembly12
  printf ', '
  print_count_n_nucleotids $assembly0
  echo
done;

# create alignments with reference sequence and all trim count sequences
rm -rf alignment/*
for file in $(ls ../Daten/Sequencing_07.06.19/190607_M02279_0447_000000000-CDLF5_NS_Lacerta/FastQ/S*R1*fastq.gz); do
  sample_r1=$(basename -- ${file})
  id=$( echo $sample_r1 | cut -d'_' -f 1 )
  assembly39=$( grep "Final assembly result will be written to file:" mitobim/39/logs/$id.log | grep -Eo "[^ ]+$" )
  assembly25=$( grep "Final assembly result will be written to file:" mitobim/25/logs/$id.log | grep -Eo "[^ ]+$" )
  assembly12=$( grep "Final assembly result will be written to file:" mitobim/12/logs/$id.log | grep -Eo "[^ ]+$" )
  assembly0=$( grep "Final assembly result will be written to file:" mitobim/0/logs/$id.log | grep -Eo "[^ ]+$" )

  mkdir -p alignment/$id
  pushd alignment/$id
    cat ../../../Daten/references/boehme_AM176577.1_lv.fas >> alignment.fas
    alignment $assembly0 $id "trimmed_MINLEN:0"
    alignment $assembly12 $id "trimmed_MINLEN:12"
    alignment $assembly25 $id "trimmed_MINLEN:25"
    alignment $assembly39 $id "trimmed_MINLEN:39"
    mafft --globalpair alignment.fas > aligment.maffted.fas
  popd
done;

