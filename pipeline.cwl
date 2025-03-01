cwlVersion: v1.2
class: CommandLineTool
baseCommand: sh
requirements:
  - class: InlineJavascriptRequirement
  - class: InitialWorkDirRequirement
    listing:
      - entryname: R1.fastq.gz
        path: input/R1.fastq.gz
      - entryname: R2.fastq.gz
        path: input/R2.fastq.gz
      - entryname: alu_lib.fa
        path: lalu_lib.fa
      - entryname: parseBAMfield_fix.py
        path: parseBAMfield_fix.py
      - entryname: pre_deseq.py
        path: pre_deseq.py
  - class: ResourceRequirement
    ramMin: 16GB
    coresMin: 16
inputs:
  ID:
    type: string
    inputBinding:
      prefix: '-c'
  alu_lib:
    type: string
    inputBinding:
      prefix: '-l'
  fastq_dir:
    type: Directory
    inputBinding:
      prefix: '-i'
  out_dir:
    type: Directory
    inputBinding:
      prefix: '-o'
outputs:
  - id: count3abv_halp
    type: File
    outputBinding:
      glob: count3abv_halp_${ID}.txt
stdout: log.txt
stderr: err.txt
arguments:
  - valueFrom: |
      printf '#!/bin/bash\n'
      printf 'set -euxo pipefail\n'
      printf 'mkdir -p ${TMPDIR:-/tmp}/out\n'
      printf 'pear -f ${fastq_dir}/${ID}_R1.fastq.gz -r ${fastq_dir}/${ID}_R2.fastq.gz -j 16 -o ${TMPDIR:-/tmp}/out/ova1_${ID}\n'
      printf 'bwa mem -t 16 -M -C ${TMPDIR:-/tmp}/lib/${alu_lib}.fa ${TMPDIR:-/tmp}/out/${ID}.assembled.fastq > ${TMPDIR:-/tmp}/out/${ID}.aln.sam\n'
      printf 'sambamba view -S ${TMPDIR:-/tmp}/out/${ID}.aln.sam -f bam -t 16 > ${TMPDIR:-/tmp}/out/${ID}.aln.bam\n'
      printf 'sambamba sort ${TMPDIR:-/tmp}/out/${ID}.aln.bam -t 16 -o ${TMPDIR:-/tmp}/out/${ID}.aln.sort.bam\n'
      printf 'python3 ${TMPDIR:-/tmp}/pre_deseq.py <(cut -f 5 ${TMPDIR:-/tmp}/out/sortuniq_315bp_10var_${ID}.txt) > ${TMPDIR:-/tmp}/out/halp_${ID}.txt\n'
      printf 'awk -F\\t \'BEGIN {OFS = FS} {if($2>=3) print $0}\' ${TMPDIR:-/tmp}/out/halp_${ID}.txt > ${TMPDIR:-/tmp}/out/count3abv_halp_${ID}.txt\n'
      printf 'exit 0\n'
    position: 1
  - valueFrom: |
      printf '#!/bin/bash\n'
      printf 'set -euxo pipefail\n'
      printf 'python3 ${TMPDIR:-/tmp}/parseBAMfield_fix.py ${TMPDIR}
