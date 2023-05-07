#
# Splice Site Score Calculation
#

```bash
conda create -n maxentscan
conda activate maxentscan
conda install -c bioconda maxentscan
```
```bash
wget https://ftp.ebi.ac.uk/pub/databases/RNAcentral/current_release/genome_coordinates/gff3/homo_sapiens.GRCh38.gff3.gz
wget https://ngdc.cncb.ac.cn/lncbook/files/LncBookv2.0_GENCODEv34_GRCh38.gtf.gz
```

#
# lncDC-train.py
#

When ViennaRNA is not properly installed, You can install it by pypl (pip install ViennaRNA). 
```bash
# $ export PATH=$HOME/anaCogentPy37/bin:$PATH 
# $ conda env create -f lncDC.yml
$ source activate lncDC.env
$ pip install ViennaRNA
```
```bash
(master)[uchiyamat@edo: ~/TEST/WORKSPACE/LncDC]$ python3 bin/lncDC-train.py -m /export/home/uchiyamat/TEST/WORKSPACE/reference/gencode.v43.pc_transcripts.fa.gz -c /export/home/uchiyamat/TEST/WORKSPACE/reference/gencode.v43.cds.fa.gz -l /export/home/uchiyamat/TEST/WORKSPACE/reference/gencode.v43.lncRNA_transcripts.fa.gz -o ./test_v01_230425

Process Start.
Checking if the training files exist ...
File /export/home/uchiyamat/TEST/WORKSPACE/reference/gencode.v43.pc_transcripts.fa.gz exist.
File /export/home/uchiyamat/TEST/WORKSPACE/reference/gencode.v43.cds.fa.gz exist.
File /export/home/uchiyamat/TEST/WORKSPACE/reference/gencode.v43.lncRNA_transcripts.fa.gz exist.

Checking if the training files are in fasta format ...
/export/home/uchiyamat/TEST/WORKSPACE/reference/gencode.v43.pc_transcripts.fa.gz format checking: PASS
/export/home/uchiyamat/TEST/WORKSPACE/reference/gencode.v43.cds.fa.gz format checking: PASS
/export/home/uchiyamat/TEST/WORKSPACE/reference/gencode.v43.lncRNA_transcripts.fa.gz format checking: PASS

Initializing dataframe ...
Total Number of transcripts loaded: 169299
Calculating transcript lengths ...

Removing Non-valid transcripts (sequence that have non-ATGCatgc letters & sequence length less than 200 nt) ...
Number of valid transcripts for training: 168442

Extracting SIF and PF features ...
Writing the hexamer table to file: /NGSWORK/TEST/TEST_UCHIYAMA/WORKSPACE/LncDC/./test_v01_230425_hexamer_table.csv
Imputing missing values ...
Writing the imputer to file: /NGSWORK/TEST/TEST_UCHIYAMA/WORKSPACE/LncDC/./test_v01_230425_imputer_SIF_PF.pkl
Standardizing ...
Writing the scaler to file: /NGSWORK/TEST/TEST_UCHIYAMA/WORKSPACE/LncDC/./test_v01_230425_scaler_SIF_PF.pkl
Balancing training data ...
Model fitting ...
Writing the model to file: /NGSWORK/TEST/TEST_UCHIYAMA/WORKSPACE/LncDC/./test_v01_230425_xgb_model_SIF_PF.pkl
Done!
```

#
# Pipeline Test
#

・CogentAPのrefacteringの問題について <br>
　- demuliplxでバーコード分配するため、数多くのfastqが作成される。<br>
　　→ 計算時間はどうか？ まとめてやっているのか？
　- 現在、遺伝子単位と転写単位の算出を１つのツールにまとめることを検討している。 <br>
　- 現状、faetureCount（Gene単位）とrsem（転写産物単位）を使用。　<br>
    - featureCountで、Gene単位では正確だが転写単位ではmultimapによるambiguousな状態に陥りやすいので不正確。ただ、中間ファイルでexon/intron/intergenic由来のリードが分かる。 <br>
    - rsemは、最大最尤法で数理モデルによる算出法なので、近似値を求めているが精度は高い。ただ、exon/intron/intergenic由来のリードがわからない。<br>
　- 中間ファイルを使うことの利点 <br>
　  - 1000細胞分、1000個のfastqファイルに1つずつバーコード配列をつけていく、UMIのduplicateの問題もつぎはぎで付けることができる。<br>
    - 現在、featureCountの中間ファイルから、geneとtranscriptそれぞれでどの領域にassignされたリードなのかを表すstatsが出力される。 <br>
    　その時、geneとtranscriptのそれぞれの、read mapping status(assign/ambiguous/unassign)から、exon/intron/intergenicに張り付いているのかを条件式で調べている。<br>
    　→ ただ、RseQCツールでも代用できるかもしれない。<br>
　- salomnを用いてmultithreadで検討したが、cpuを爆弾に消費し、メモリー3TBを要するレベルに達する見込み。つまり、salomnはBulk用でsingle-cellには不向き。<br>


## Download FASTQ data
・Web : [sra-explorer](https://sra-explorer.info/) <br>
・Search keyword: 「single-cell RNA-Seq data」<br>
```bash
(base) [uchiyamat@edo: ~/TEST/tools]$ wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-centos_linux64.tar.gz
(base) [uchiyamat@edo: ~/TEST/tools]$ tar -zxvf sratoolkit.current-centos_linux64.tar.gz
(base) [uchiyamat@edo: ~/TEST/WORKSPACE]$ mkdir -p FASTQ
```
```bash
# First time 
# Single cell RNA-seq data for iPSC differentiation into Syndetome
# Single-cell RNA-Seq libraries were prepared per the Single Cell 3′ v3.1 Reagent Kits
# User Guide (10x Genomics, Pleasanton, California) using the 10x Genomics Chromium
# Controller. Barcoded sequencing libraries were quantified by quantitative PCR using the Collibri
# Library Quantification Kit (Thermo Fisher Scientific, Waltham, MA). Libraries were sequenced on a NovaSeq 6000 (Illumina, San Diego, CA) as per the Single Cell 3′ v3.1 Reagent Kits User
# Guide, with a sequencing depth of ~40,000 reads/cell. 
#(base) [uchiyamat@edo: ~/TEST/WORKSPACE/FASTQ]$ /export/home/uchiyamat/TEST/tools/sratoolkit.3.0.1-centos_linux64/bin/prefetch --option-file SRR_Acc_List.txt
# /export/home/uchiyamat/TEST/tools/sratoolkit.3.0.1-centos_linux64/bin/fastq-dump --split-files --gzip --aligned -Q 64 SRR24086035
# /export/home/uchiyamat/TEST/tools/sratoolkit.3.0.1-centos_linux64/bin/fastq-dump --split-files --gzip --aligned -Q 64 SRR24086036
# /export/home/uchiyamat/TEST/tools/sratoolkit.3.0.1-centos_linux64/bin/fastq-dump --split-files --gzip --aligned -Q 64 SRR24086031
# /export/home/uchiyamat/TEST/tools/sratoolkit.3.0.1-centos_linux64/bin/fastq-dump --split-files --gzip --aligned -Q 64 SRR24086032

# Each Download file was split into three FASTQ files as unexpected. So, I remove all FASTQ files, then I redownloaded other file.
```
```bash
# Twice time
# random selection data ---
# GSM7134459: PAT2V1, baseline; Homo sapiens; RNA-Seq
(base) [uchiyamat@edo: ~/TEST/WORKSPACE/FASTQ]$ /export/home/uchiyamat/TEST/tools/sratoolkit.3.0.1-centos_linux64/bin/fastq-dump --split-files --gzip SRR24035598
```
```bash
# Singletome data ---
# liver set 1 (GSE115469 (MacParland et al. 2018))
# GSM3178782: Patient 1 Total Liver Homogenate; Homo sapiens; RNA-Seq
# Dissecting the human liver cellular landscape by single cell RNA-seq reveals novel intrahepatic monocyte/ macrophage populations
#(base) [uchiyamat@edo: ~/TEST/WORKSPACE/FASTQ]$ /export/home/uchiyamat/TEST/tools/sratoolkit.3.0.1-centos_linux64/bin/fasterq-dump --split-files　SRR7276474 -e 32 -p
# -> Although SRA reported PE read, Download file was PE read. I removed it.
```
```bash
# Singletome data ---
# liver set 2 (GSE136103 (Ramachandran et al. 2019)) 
[uchiyamat@edo: ~/TEST/WORKSPACE]$ /export/home/uchiyamat/TEST/tools/sratoolkit.3.0.1-centos_linux64/bin/fastq-dump --split-files --gzip SRR10009422 &
[1] 111167
```


```bash
# Submitted by: Harbin Medical University Cancer Hospital
# Study: FN1 regulates aspartate metabolism in breast cancer
# Illumina NovaSeq 6000

[uchiyamat@edo: ~/TEST/WORKSPACE/FASTQ]$ /export/home/uchiyamat/TEST/tools/sratoolkit.3.0.1-centos_linux64/bin/fastq-dump --split-files --gzip SRR24332180 &
[1] 89341
[uchiyamat@edo: ~/TEST/WORKSPACE/FASTQ]$ /export/home/uchiyamat/TEST/tools/sratoolkit.3.0.1-centos_linux64/bin/fastq-dump --split-files --gzip SRR24332179 &
[2] 89415
[uchiyamat@edo: ~/TEST/WORKSPACE/FASTQ]$ /export/home/uchiyamat/TEST/tools/sratoolkit.3.0.1-centos_linux64/bin/fastq-dump --split-files --gzip SRR24332178 &
[3] 89438
```

## Make reference genome
・Reference：[Comparative Analysis of common alignment tools for single cell RNAsequencing](https://www.biorxiv.org/content/10.1101/2021.02.15.430948v1.full.pdf) <br>
Cell Ranger includes STAR, which was designed for bulk RNA-seq data. STAR performs a classical alignment approach by utilizing a maximal mappable seed search, thereby all possible positions of the reads can be determined. <br>
↓↓ <br>
Do I have to address issue of overlapping read between lncRNA and mRNA exons ? <br>

### Sense strand overlap.
Cell Ranger discards reads mapped to overlapping exons on the same strand(red dotted box). 

To avoid miscounting reads to protein-coding genes by the inclusion of additional lncRNAs in the genome, lncRNA genes were discarded if they overlap in sense with protein-coding exons (red x), as it is more difficult to exclude the protein-coding potential of these lncRNAs. 

### Antisense strand overlap.
Cell Ranger prioritizes alignments of sense over antisense reads. If spurious antisense reads are generated from transcripts of protein-coding genes, these could be incorrectly interpreted to indicate expression of an antisense overlapping lncRNA gene. 

To overcome this potential problem, we trimmed the overlapping region (red x) and an additional 100nt of lncRNA exons that were overlapping with protein-coding exons in the antisense direction. We retained the trimmed lncRNA exons if their length was at least 200nt(marked with white check).

Gene and transcripts coordinates were updated accordingly. 

## Prepare for making reference genome
```bash
# install tools
(base) [uchiyamat@edo: ~/TEST/tools]$ wget https://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/bedToGenePred
(base) [uchiyamat@edo: ~/TEST/tools]$ chmod +x bedToGenePred
(base) [uchiyamat@edo: ~/TEST/tools]$ wget https://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/genePredToGtf
(base) [uchiyamat@edo: ~/TEST/tools]$ chmod +x genePredToGtf
```
```bash
# install RNAcentral dataset
(base) [uchiyamat@gifu: ~/TEST/WORKSPACE/reference]$ pwd
/export/home/uchiyamat/TEST/WORKSPACE/reference
(base) [uchiyamat@gifu: ~/TEST/WORKSPACE/reference]$ wget https://ftp.ebi.ac.uk/pub/databases/RNAcentral/current_release/genome_coordinates/bed/homo_sapiens.GRCh38.bed.gz
```
```bash
# Edit RNAcentral dataset
# bedtogtf
(base) [uchiyamat@gifu: ~/TEST/WORKSPACE/reference]$ bedToGenePred homo_sapiens.GRCh38.bed.gz /dev/stdout | genePredToGtf file /dev/stdin rnacentral.GRCh38.gtf
```
```bash
# Modification of gtf file
# (1) Get moleculer type from RNA central bed file
file = 'A:/TEST/TEST_UCHIYAMA/WORKSPACE/reference/homo_sapiens.GRCh38.bed.gz'
type = {}
with gzip.open(file, 'rt') as bed:
    for line in bed.readlines():
        line = line.rstrip('\n').split()
        chrom = line[0]
        start = str(int(line[1]) + 1)
        end =line[2]
        id = line[3].split('_')[0]
        transcripts_id = chrom + ':' + start + ':' + end + ':' + id
        type[transcripts_id] = line[-2]


# (2) Read pre-gtf file
# chr1    /dev/stdin      transcript      10245   10273   .       -       .       gene_id "URS000035F234_9606"; transcript_id "URS000035F234_9606";
# chr1    /dev/stdin      exon    10245   10263   .       -       .       gene_id "URS000035F234_9606"; transcript_id "URS000035F234_9606"; exon_number "2"; exon_id "URS000035F234_9606.2";
#  gene_type "lncRNA"; transcript_type "lncRNA";
file = 'A:/TEST/TEST_UCHIYAMA/WORKSPACE/reference/rnacentral.GRCh38.gtf'
transcript_line = {}
exon_line = {} 
exon_info = {}
with open(file, 'r') as gtf:
    for line in gtf.readlines():
        line = line.replace('/dev/stdin', 'RNAcentral').rstrip('\n').split('\t')
        if '#' in line[0]:
            pass
        else:
            chrom = line[0]
            start = line[3]
            end =line[4]
            strand = line[6]
            if 'transcript' in line[2]:
                    flag = 1
                    for ele in line[-1].split('; '):
                        if 'transcript_id' in ele:
                            transcript_id = re.search(r'transcript_id \"(.*)\"', ele).group(1)
                            tx_id = transcript_id.split('_')[0]
                            id = chrom + ':' + start + ':' + end + ':' + tx_id
                    
                    biotype = type[id]
                    gene_type = 'gene_type \"' + biotype + '\";'
                    transcript_type = 'transcript_type \"' + biotype + '\";'
                    gtf_id = transcript_id + ':' + strand
                    transcript_line[gtf_id] = '\t'.join(line) +  gene_type + ' ' + transcript_type 
                
            elif 'exon' in line[2]:
                if flag == 1:
                    tmp = '\t'.join(line) + gene_type + ' ' + transcript_type 
                    exon = gtf_id + ':' + chrom + ':' + start + ':' + end
                    exon_line[exon] = tmp
                    if exon_info.get(gtf_id):
                        exon_info[gtf_id] += ',' + exon
                    else:
                        exon_info[gtf_id] = exon
            
            elif 'codon' in line[2]:
                flag = 0

# (3) Output modified gtf file 
outfile = 'A:/TEST/TEST_UCHIYAMA/WORKSPACE/reference/rnacentral.GRCh38.modify.gtf'
with open(outfile, 'w') as wgtf:
    for gtf_id, line in transcript_line.items():
        transcript_id, strand = gtf_id.split(':')
        flag = 0
        prev = 0
        for exon in exon_info[gtf_id].split(','):
            _, _, _, start, end = exon.split(':')
            if int(start) <= prev:
                flag = 1
            prev = int(end)
    
        if flag == 0:
            gene_line = line.replace('transcript', 'gene', 1)
            wgtf.write(f'{gene_line}\n')
            wgtf.write(f'{line}\n')
            if strand == '+':
                exons = exon_info[gtf_id].split(',')
            elif strand == '-':
                exons = reversed(exon_info[gtf_id].split(','))
            for exon in exons:
                line = exon_line[exon]
                wgtf.write(f'{line}\n')
```

### Use Singletome method 
```bash
# get mRNA recode from GENCODE.v43
infile = 'A:/TEST/TEST_UCHIYAMA/WORKSPACE/reference/gencode.v43.annotation.gtf.gz'
outfile = 'A:/TEST/TEST_UCHIYAMA/WORKSPACE/reference/gencode.v43.mRNAs.gtf'
flag = 0
with open(outfile, 'w') as wgtf:
    with gzip.open(infile, 'rt') as gtf:
        for line in gtf.readlines():
            tmp = line
            line = line.rstrip('\n').split('\t')
            if '#' in line[0]:
                wgtf.write(f'{tmp}')
            else:
                if line[2] == 'gene':
                    flag = 0
                    for ele in line[-1].split('; '):
                        if 'gene_type' in ele:
                            gene_type = re.search(r'gene_type \"(.*)\"', ele).group(1)
                    if gene_type == 'protein_coding':
                        wgtf.write(f'{tmp}')
                        flag = 1
                if line[2] == 'transcript':
                    if flag == 1:
                        for ele in line[-1].split('; '):
                            if 'transcript_id' in ele:
                                transcript_id = re.search(r'transcript_id \"(.*)\"', ele).group(1)
                            if 'gene_type' in ele:
                                gene_type = re.search(r'gene_type \"(.*)\"', ele).group(1)
                            if 'transcript_type' in ele:
                                transcript_type = re.search(r'transcript_type \"(.*)\"', ele).group(1)
                        if (gene_type == 'protein_coding') & (transcript_type == 'protein_coding'):
                            wgtf.write(f'{tmp}')
                            flag = 2
                else:
                    if flag == 2:
                        wgtf.write(f'{tmp}')
```

### GTF file for IGV
```bash
# GENCODE
[uchiyamat@edo: ~/TEST/WORKSPACE/reference]$ (grep ^"#" ./gencode.v43.annotation.gtf; grep -v ^"#" ./gencode.v43.annotation.gtf | sort -t $'\t' -k1,1V -k4,4n -k5,5n) | bgzip  > ./gencode.v43.annotation.gtf.gz
[uchiyamat@edo: ~/TEST/WORKSPACE/reference]$ tabix -p gff ./gencode.v43.annotation.gtf.gz
```
```bash
# Singletome
[uchiyamat@edo: ~/TEST/WORKSPACE/pipeline/00_dataset]$ (grep ^"#" ./cleanedForSense_and_Antisense.gtf; grep -v ^"#" ./cleanedForSense_and_Antisense.gtf | sort -t $'\t' -k1,1V -k4,4n -k5,5n) | bgzip  > ./cleanedForSense_and_Antisense.gtf.gz
[uchiyamat@edo: ~/TEST/WORKSPACE/pipeline/00_dataset]$ tabix -p gff ./cleanedForSense_and_Antisense.gtf.gz

[uchiyamat@edo: ~/TEST/WORKSPACE/pipeline/00_dataset]$ (grep ^"#" ./cleaned_Same_Strand_Exons.gtf; grep -v ^"#" ./cleaned_Same_Strand_Exons.gtf | sort -t $'\t' -k1,1V -k4,4n -k5,5n) | bgzip  > ./cleaned_Same_Strand_Exons.gtf.gz
[uchiyamat@edo: ~/TEST/WORKSPACE/pipeline/00_dataset]$ tabix -p gff ./cleaned_Same_Strand_Exons.gtf.gz

[uchiyamat@edo: ~/TEST/WORKSPACE/pipeline/00_dataset]$ (grep ^"#" ./same_Strand_Deleted_Genes.gtf; grep -v ^"#" ./same_Strand_Deleted_Genes.gtf | sort -t $'\t' -k1,1V -k4,4n -k5,5n) | bgzip  > ./same_Strand_Deleted_Genes.gtf.gz
[uchiyamat@edo: ~/TEST/WORKSPACE/pipeline/00_dataset]$ tabix -p gff ./same_Strand_Deleted_Genes.gtf.gz
```


## Draft for Pipeline Construction
```bash
(base) [uchiyamat@edo: ~/TEST/WORKSPACE]$ mkdir -p pipeline
```

## (1) Read trimming using fastp
```bash
# cheak
(base) [uchiyamat@edo: ~/TEST/tools]$ which fastp
/NGSWORK/NGS/GridEngine/bin/fastp

# version
fastp: an ultra-fast all-in-one FASTQ preprocessor
version 0.20.1
```
#### Case_01
```bash
# SRR24035598; PAT2V1, baseline; Homo sapiens
(base) [uchiyamat@edo: ~/TEST/WORKSPACE/pipeline]$ mkdir -p 01_fastp
(base) [uchiyamat@edo: ~/TEST/WORKSPACE/pipeline]$ cd 01_fastp
(base) [uchiyamat@edo: ~/TEST/WORKSPACE/pipeline/01_fastp]$ fastp -i /export/home/uchiyamat/TEST/WORKSPACE/FASTQ/SRR24035598_1.fastq.gz -I /export/home/uchiyamat/TEST/WORKSPACE/FASTQ/SRR24035598_2.fastq.gz -o ./SRR24035598_R1.trimed.fastq.gz -O ./SRR24035598_R2.trimed.fastq.gz -h ./SRR24035598_report.html -j ./SRR24035598_report.json -q 30
```
#### Case_02
```bash
# Study: FN1 regulates aspartate metabolism in breast cancer ---
# SRR24332178
[uchiyamat@edo: ~/TEST/WORKSPACE/pipeline/01_fastp]$ fastp -i /export/home/uchiyamat/TEST/WORKSPACE/FASTQ/SRR24332178_1.fastq.gz -I /export/home/uchiyamat/TEST/WORKSPACE/FASTQ/SRR24332178_2.fastq.gz -o ./SRR24332178_R1.trimed.fastq.gz -O ./SRR24332178_R2.trimed.fastq.gz -h ./SRR24332178_report.html -j ./SRR24332178_report.json -q 30 &
[2] 113071

# SRR24332179
[uchiyamat@edo: ~/TEST/WORKSPACE/pipeline/01_fastp]$ fastp -i /export/home/uchiyamat/TEST/WORKSPACE/FASTQ/SRR24332179_1.fastq.gz -I /export/home/uchiyamat/TEST/WORKSPACE/FASTQ/SRR24332179_2.fastq.gz -o ./SRR24332179_R1.trimed.fastq.gz -O ./SRR24332179_R2.trimed.fastq.gz -h ./SRR24332179_report.html -j ./SRR24332179_report.json -q 30 &
[3] 113096

# SRR24332180
[uchiyamat@edo: ~/TEST/WORKSPACE/pipeline/01_fastp]$ fastp -i /export/home/uchiyamat/TEST/WORKSPACE/FASTQ/SRR24332180_1.fastq.gz -I /export/home/uchiyamat/TEST/WORKSPACE/FASTQ/SRR24332180_2.fastq.gz -o ./SRR24332180_R1.trimed.fastq.gz -O ./SRR24332180_R2.trimed.fastq.gz -h ./SRR24332180_report.html -j ./SRR24332180_report.json -q 30 &
[4] 113125
```
#### Case_03
```bash
# DEV1835_37A
[uchiyamat@edo: ~/TEST/WORKSPACE/pipeline/01_fastp]$ fastp -i /export/home/uchiyamat/TEST/WORKSPACE/FASTQ/HFHHHDSXX_DEV1835_37A0302_H2_L004_R1.fastq.gz -I /export/home/uchiyamat/TEST/WORKSPACE/FASTQ/HFHHHDSXX_DEV1835_37A0302_H2_L004_R2.fastq.gz -o ./DEV1835_37A_R1.trimed.fastq.gz -O ./DEV1835_37A_R2.trimed.fastq.gz -h ./DEV1835_37A_report.html -j ./DEV1835_37A_report.json -q 30 &
[1] 117242
```

## (2) mapping and assemble transcripts
```bash
(base) [uchiyamat@edo: ~/TEST/WORKSPACE/pipeline]$ mkdir -p 01_mapping/route01
(base) [uchiyamat@edo: ~/TEST/WORKSPACE/pipeline]$ mkdir -p 01_mapping/route02
```

【route01】<br>
### (1) mapping with STAR
```bash
# cheak version
[uchiyamat@edo: ~/TEST/tools/STAR-2.7.10b]$ STAR --version
STAR_2.5.2b
```
```bash
# Build-up index
[uchiyamat@edo: ~/TEST/WORKSPACE/pipeline/01_mapping/route01]$ mkdir -p 00_index
[uchiyamat@edo: ~/TEST/WORKSPACE/pipeline/01_mapping/route01/00_index]$ /NGSWORK/NGS/GridEngine/TOOLS/CogentAP/CogentAP-v1.5.1/CogentAP/lib/../CogentAP_tools/bin/STAR \
> --runMode genomeGenerate --runThreadN 32 \
> --genomeDir 00_index \
> --genomeFastaFiles /export/home/uchiyamat/TEST/WORKSPACE/reference/GRCh38.primary_assembly.genome.fa \
> --limitGenomeGenerateRAM 55000000000 \
> --sjdbGTFfile /export/home/uchiyamat/TEST/WORKSPACE/reference/gencode.v43.annotation.gtf
```
#### Case_01
```bash
5.2.3 Compatibility with Cufflinks/Cuffdiff.
For unstranded RNA-seq data, Cufflinks/Cuffdiff require spliced alignments with XS strand attribute,
which STAR will generate with --outSAMstrandField intronMotif option. As required, the XS strand attribute will be generated for all alignments that contain splice junctions.
The spliced alignments that have undefined strand (i.e. containing only non-canonical unannotated junctions) will be suppressed.

If you have stranded RNA-seq data, you do not need to use any specific STAR options. Instead,
you need to run Cufflinks with the library option --library-type options. For example, cufflinks 
--library-type fr-firststrand should be used for the “standard” dUTP protocol, including Illumina’s stranded Tru-Seq.
This option has to be used only for Cufflinks runs and not for STAR runs.
In addition, it is recommended to remove the non-canonical junctions for Cufflinks runs using
--outFilterIntronMotifs RemoveNoncanonical.
```
・Reverses strandedness tags (XS) in BAM. Every XS:A:+ tag is changed to XS:A:- and vice versa. <br>
・XS tags are used by the Cuffdiff and Cufflinks tools when assembling transcripts.　<br>
```bash
# Basic parameter --- 
#[uchiyamat@edo: ~/TEST/WORKSPACE/pipeline/01_mapping/route01]$ mkdir -p 01_mapping
#[uchiyamat@edo: ~/TEST/WORKSPACE/pipeline/01_mapping/route01]$ cd 01_mapping/
#[uchiyamat@edo: ~/TEST/WORKSPACE/pipeline/01_mapping/route01/01_mapping]$ STAR --genomeDir /export/home/uchiyamat/TEST/WORKSPACE/pipeline/01_mapping/route01/00_index \
#--runThreadN 24 --sjdbGTFfile /export/home/uchiyamat/TEST/WORKSPACE/reference/gencode.v43.annotation.gtf \
#--outSAMtype BAM SortedByCoordinate \  # similar to samtools sort command
#--outFileNamePrefix SRR24035598. \
#--outSAMstrandField intronMotif \
#--readFilesCommand zcat \
#--readFilesIn /export/home/uchiyamat/TEST/WORKSPACE/pipeline/01_fastp/SRR24035598_R1.trimed.fastq.gz /export/home/uchiyamat/TEST/WORKSPACE/pipeline/01_fastp/SRR24035598_R2.trimed.fastq.gz
```

#### Case_02
```bash
# SRR24332178 ---
[uchiyamat@edo: ~/TEST/WORKSPACE/pipeline/01_mapping/route01/01_mapping]$ /NGSWORK/NGS/GridEngine/TOOLS/CogentAP/CogentAP-v1.5.1/CogentAP/lib/../CogentAP_tools/bin/STAR \
> --genomeDir /export/home/uchiyamat/TEST/WORKSPACE/pipeline/01_mapping/route01/00_index \
> --outFileNamePrefix ./SRR24332178. \
> --runThreadN 24 \
> --outSAMtype BAM Unsorted \
> --genomeLoad LoadAndKeep \
> --outReadsUnmapped Fastx \
> --outSAMstrandField intronMotif \
> --chimSegmentMin 12 \
> --chimJunctionOverhangMin 8 \
> --chimOutJunctionFormat 1 \
> --alignSJDBoverhangMin 10 \
> --alignMatesGapMax 100000 \
> --alignIntronMax 100000 \
> --alignSJstitchMismatchNmax 5 -1 5 5 \
> --outSAMattrRGline ID:GRPundef \
> --chimMultimapScoreRange 3 \
> --chimMultimapNmax 20 \
> --chimNonchimScoreDropMin 10 \
> --peOverlapNbasesMin 12 \
> --peOverlapMMp 0.1 \
> --alignInsertionFlush Right \
> --alignSplicedMateMapLminOverLmate 0 \
> --alignSplicedMateMapLmin 30 \
> --readFilesCommand zcat \
> --readFilesIn /export/home/uchiyamat/TEST/WORKSPACE/pipeline/01_fastp/SRR24332178_R1.trimed.fastq.gz /export/home/uchiyamat/TEST/WORKSPACE/pipeline/01_fastp/SRR24332178_R2.trimed.fastq.gz &
[1] 117666

# SRR24332179 ---
[uchiyamat@edo: ~/TEST/WORKSPACE/pipeline/01_mapping/route01/01_mapping]$ /NGSWORK/NGS/GridEngine/TOOLS/CogentAP/CogentAP-v1.5.1/CogentAP/lib/../CogentAP_tools/bin/STAR \
> --genomeDir /export/home/uchiyamat/TEST/WORKSPACE/pipeline/01_mapping/route01/00_index \
> --outFileNamePrefix ./SRR24332179. \
> --runThreadN 24 \
> --outSAMtype BAM Unsorted \
> --genomeLoad LoadAndKeep \
> --outReadsUnmapped Fastx \
> --outSAMstrandField intronMotif \
> --chimSegmentMin 12 \
> --chimJunctionOverhangMin 8 \
> --chimOutJunctionFormat 1 \
> --alignSJDBoverhangMin 10 \
> --alignMatesGapMax 100000 \
> --alignIntronMax 100000 \
> --alignSJstitchMismatchNmax 5 -1 5 5 \
> --outSAMattrRGline ID:GRPundef \
> --chimMultimapScoreRange 3 \
> --chimMultimapNmax 20 \
> --chimNonchimScoreDropMin 10 \
> --peOverlapNbasesMin 12 \
> --peOverlapMMp 0.1 \
> --alignInsertionFlush Right \
> --alignSplicedMateMapLminOverLmate 0 \
> --alignSplicedMateMapLmin 30 \
> --readFilesCommand zcat \
> --readFilesIn /export/home/uchiyamat/TEST/WORKSPACE/pipeline/01_fastp/SRR24332179_R1.trimed.fastq.gz /export/home/uchiyamat/TEST/WORKSPACE/pipeline/01_fastp/SRR24332179_R2.trimed.fastq.gz &
[2] 118338

# SRR24332180 ---
[uchiyamat@edo: ~/TEST/WORKSPACE/pipeline/01_mapping/route01/01_mapping]$ /NGSWORK/NGS/GridEngine/TOOLS/CogentAP/CogentAP-v1.5.1/CogentAP/lib/../CogentAP_tools/bin/STAR \
> --genomeDir /export/home/uchiyamat/TEST/WORKSPACE/pipeline/01_mapping/route01/00_index \
> --outFileNamePrefix ./SRR24332180. \
> --runThreadN 24 \
> --outSAMtype BAM Unsorted \
> --genomeLoad LoadAndKeep \
> --outReadsUnmapped Fastx \
> --outSAMstrandField intronMotif \
> --chimSegmentMin 12 \
> --chimJunctionOverhangMin 8 \
> --chimOutJunctionFormat 1 \
> --alignSJDBoverhangMin 10 \
> --alignMatesGapMax 100000 \
> --alignIntronMax 100000 \
> --alignSJstitchMismatchNmax 5 -1 5 5 \
> --outSAMattrRGline ID:GRPundef \
> --chimMultimapScoreRange 3 \
> --chimMultimapNmax 20 \
> --chimNonchimScoreDropMin 10 \
> --peOverlapNbasesMin 12 \
> --peOverlapMMp 0.1 \
> --alignInsertionFlush Right \
> --alignSplicedMateMapLminOverLmate 0 \
> --alignSplicedMateMapLmin 30 \
> --readFilesCommand zcat \
> --readFilesIn /export/home/uchiyamat/TEST/WORKSPACE/pipeline/01_fastp/SRR24332180_R1.trimed.fastq.gz /export/home/uchiyamat/TEST/WORKSPACE/pipeline/01_fastp/SRR24332180_R2.trimed.fastq.gz &
[3] 118463
```
```bash
[uchiyamat@edo: ~/TEST/WORKSPACE/pipeline/01_mapping/route02/01_mapping]$ samtools sort -@ 16 -O bam -o ./SRR24332178.sort.bam ./SRR24332178.Aligned.out.bam &
[1] 119809
[uchiyamat@edo: ~/TEST/WORKSPACE/pipeline/01_mapping/route02/01_mapping]$ samtools sort -@ 16 -O bam -o ./SRR24332179.sort.bam ./SRR24332179.Aligned.out.bam &
[2] 119879
[uchiyamat@edo: ~/TEST/WORKSPACE/pipeline/01_mapping/route02/01_mapping]$ samtools sort -@ 16 -O bam -o ./SRR24332180.sort.bam ./SRR24332180.Aligned.out.bam
[3] 119988

[uchiyamat@edo: ~/TEST/WORKSPACE/pipeline/01_mapping/route01/01_mapping]$ samtools index ./SRR24332178.sort.bam &
[1] 120434
[uchiyamat@edo: ~/TEST/WORKSPACE/pipeline/01_mapping/route01/01_mapping]$ samtools index ./SRR24332179.sort.bam &
[2] 120445
[uchiyamat@edo: ~/TEST/WORKSPACE/pipeline/01_mapping/route01/01_mapping]$ samtools index ./SRR24332180.sort.bam &
[3] 120453
```

## (2) Assemble transcripts with cufflinks

```bash
# cufflinks
[uchiyamat@edo: ~/TEST/WORKSPACE/pipeline/01_mapping/route01]$ mkdir -p 02_cufflinks
[uchiyamat@edo: ~/TEST/WORKSPACE/pipeline/01_mapping/route01]$ cd 02_cufflinks/
```
#### Case_01
```bash
# assemble transcripts
[uchiyamat@edo: ~/TEST/WORKSPACE/pipeline/01_mapping/route01/02_cufflinks]$ cufflinks -p 16 \
-g /export/home/uchiyamat/TEST/WORKSPACE/reference/gencode.v43.annotation.gtf \
-o ./SRR24035598 \ 
/export/home/uchiyamat/TEST/WORKSPACE/pipeline/01_mapping/route02/01_mapping/SRR24035598.sort.bam
# -> 「BAM record error: found spliced alignment without XS attribute」?
```
```bash
# cuffmerge
# -> No conducted
```

#### Case_02
```bash
# assemble transcripts
[uchiyamat@edo: ~/TEST/WORKSPACE/pipeline/01_mapping/route01/02_cufflinks]$ cufflinks -p 16 -g /export/home/uchiyamat/TEST/WORKSPACE/reference/gencode.v43.annotation.gtf -o ./SRR24332178 /export/home/uchiyamat/TEST/WORKSPACE/pipeline/01_mapping/route02/01_mapping/SRR24332178.sort.bam &

[uchiyamat@edo: ~/TEST/WORKSPACE/pipeline/01_mapping/route01/02_cufflinks]$ cufflinks -p 16 -g /export/home/uchiyamat/TEST/WORKSPACE/reference/gencode.v43.annotation.gtf -o ./SRR24332179 /export/home/uchiyamat/TEST/WORKSPACE/pipeline/01_mapping/route02/01_mapping/SRR24332179.sort.bam &

[uchiyamat@edo: ~/TEST/WORKSPACE/pipeline/01_mapping/route01/02_cufflinks]$ cufflinks -p 16 -g /export/home/uchiyamat/TEST/WORKSPACE/reference/gencode.v43.annotation.gtf -o ./SRR24332180 /export/home/uchiyamat/TEST/WORKSPACE/pipeline/01_mapping/route02/01_mapping/SRR24332180.sort.bam &
```
```bash
# cuffmerge
[uchiyamat@edo: ~/TEST/WORKSPACE/pipeline/01_mapping/route01/02_cufflinks]$ ls /export/home/uchiyamat/TEST/WORKSPACE/pipeline/01_mapping/route01/02_cufflinks/SRR24332*/transcripts.gtf > assembly_gtf_list.txt
cuffmerge
[uchiyamat@edo: ~/TEST/WORKSPACE/pipeline/01_mapping/route01/02_cufflinks]$ python2.7 /NGSWORK/NGS/GridEngine/bin/cuffmerge -p 16 -o ./cuffmerge -g /export/home/uchiyamat/TEST/WORKSPACE/reference/gencode.v43.annotation.gtf assembly_gtf_list.txt
```

【route02】<br>
### (1) mapping with HISAT2 
```bash
# cheak version
[uchiyamat@edo: ~/TEST/WORKSPACE/pipeline/01_mapping/route02/01_mapping]$ hisat2
No index, query, or output file specified!
HISAT2 version 2.2.1 by Daehwan Kim (infphilo@gmail.com, www.ccb.jhu.edu/people/infphilo)
```
```bash
# Build-up index
[uchiyamat@edo: ~/TEST/WORKSPACE/pipeline/01_mapping/route02]$ mkdir -p 00_index/
[uchiyamat@edo: ~/TEST/WORKSPACE/pipeline/01_mapping/route02]$ cd 00_index/
[uchiyamat@edo: ~/TEST/WORKSPACE/pipeline/01_mapping/route02/00_index]$ hisat2-build /export/home/uchiyamat/TEST/WORKSPACE/reference/GRCh38.primary_assembly.genome.fa GRCh38
```
```bash
# Set up workdir 
[uchiyamat@edo: ~/TEST/WORKSPACE/pipeline/01_mapping/route02]$ mkdir -p ./01_mapping
[uchiyamat@edo: ~/TEST/WORKSPACE/pipeline/01_mapping/route02]$ cd ./01_mapping
```
#### Case_01
```bash
# SRR24035598 ---
[uchiyamat@edo: ~/TEST/WORKSPACE/pipeline/01_mapping/route02/01_mapping]$ hisat2 --no-mixed --dta \
-p 32 \
-x /export/home/uchiyamat/TEST/WORKSPACE/pipeline/01_mapping/route02/00_index/GRCh38 \
-1 /export/home/uchiyamat/TEST/WORKSPACE/pipeline/01_fastp/SRR24035598_R1.trimed.fastq.gz \
-2 /export/home/uchiyamat/TEST/WORKSPACE/pipeline/01_fastp/SRR24035598_R2.trimed.fastq.gz \
-S ./SRR24035598.sam

271920078 reads; of these:
  271920078 (100.00%) were paired; of these:
    263793380 (97.01%) aligned concordantly 0 times
    6591049 (2.42%) aligned concordantly exactly 1 time
    1535649 (0.56%) aligned concordantly >1 times
    ----
    263793380 pairs aligned concordantly 0 times; of these:
      41046 (0.02%) aligned discordantly 1 time
3.00% overall alignment rate
```
```bash
[uchiyamat@edo: ~/TEST/WORKSPACE/pipeline/01_mapping/route02/01_mapping]$ samtools sort -@ 16 -O bam -o ./SRR24035598.sort.bam ./SRR24035598.sam
[uchiyamat@edo: ~/TEST/WORKSPACE/pipeline/01_mapping/route02/01_mapping]$ samtools index ./SRR24035598.sort.bam
[uchiyamat@edo: ~/TEST/WORKSPACE/pipeline/01_mapping/route02/01_mapping]$ rm -rf SRR24035598.sam
```

#### Case_02
```bash
# SRR24332178 ---
[uchiyamat@edo: ~/TEST/WORKSPACE/pipeline/01_mapping/route02/01_mapping]$ hisat2 --no-mixed --dta \
> -p 32 \
> -x /export/home/uchiyamat/TEST/WORKSPACE/pipeline/01_mapping/route02/00_index/GRCh38 \
> -1 /export/home/uchiyamat/TEST/WORKSPACE/pipeline/01_fastp/SRR24332178_R1.trimed.fastq.gz \
> -2 /export/home/uchiyamat/TEST/WORKSPACE/pipeline/01_fastp/SRR24332178_R2.trimed.fastq.gz \
> -S ./SRR24332178.sam &
[1] 115610

# SRR24332179 ---
[uchiyamat@edo: ~/TEST/WORKSPACE/pipeline/01_mapping/route02/01_mapping]$ hisat2 --no-mixed --dta \
> -p 32 \
> -x /export/home/uchiyamat/TEST/WORKSPACE/pipeline/01_mapping/route02/00_index/GRCh38 \
> -1 /export/home/uchiyamat/TEST/WORKSPACE/pipeline/01_fastp/SRR24332179_R1.trimed.fastq.gz \
> -2 /export/home/uchiyamat/TEST/WORKSPACE/pipeline/01_fastp/SRR24332179_R2.trimed.fastq.gz \
> -S ./SRR24332179.sam &
[2] 116037

# SRR24332180 ---
[uchiyamat@edo: ~/TEST/WORKSPACE/pipeline/01_mapping/route02/01_mapping]$ hisat2 --no-mixed --dta \
> -p 32 \
> -x /export/home/uchiyamat/TEST/WORKSPACE/pipeline/01_mapping/route02/00_index/GRCh38 \
> -1 /export/home/uchiyamat/TEST/WORKSPACE/pipeline/01_fastp/SRR24332180_R1.trimed.fastq.gz \
> -2 /export/home/uchiyamat/TEST/WORKSPACE/pipeline/01_fastp/SRR24332180_R2.trimed.fastq.gz \
> -S ./SRR24332180.sam &
[3] 116095

23102370 reads; of these:
  23102370 (100.00%) were paired; of these:
    1051656 (4.55%) aligned concordantly 0 times
    21347669 (92.40%) aligned concordantly exactly 1 time
    703045 (3.04%) aligned concordantly >1 times
    ----
    1051656 pairs aligned concordantly 0 times; of these:
      84832 (8.07%) aligned discordantly 1 time
95.82% overall alignment rate

23230287 reads; of these:
  23230287 (100.00%) were paired; of these:
    906140 (3.90%) aligned concordantly 0 times
    21635204 (93.13%) aligned concordantly exactly 1 time
    688943 (2.97%) aligned concordantly >1 times
    ----
    906140 pairs aligned concordantly 0 times; of these:
      110643 (12.21%) aligned discordantly 1 time
96.58% overall alignment rate

26042094 reads; of these:
  26042094 (100.00%) were paired; of these:
    861026 (3.31%) aligned concordantly 0 times
    24427013 (93.80%) aligned concordantly exactly 1 time
    754055 (2.90%) aligned concordantly >1 times
    ----
    861026 pairs aligned concordantly 0 times; of these:
      75250 (8.74%) aligned discordantly 1 time
96.98% overall alignment rate
```
```bash
[uchiyamat@edo: ~/TEST/WORKSPACE/pipeline/01_mapping/route02/01_mapping]$ samtools sort -@ 16 -O bam -o ./SRR24332178.sort.bam ./SRR24332178.sam &
[1] 118737
[uchiyamat@edo: ~/TEST/WORKSPACE/pipeline/01_mapping/route02/01_mapping]$ samtools sort -@ 16 -O bam -o ./SRR24332179.sort.bam ./SRR24332179.sam &
[2] 118787
[uchiyamat@edo: ~/TEST/WORKSPACE/pipeline/01_mapping/route02/01_mapping]$ samtools sort -@ 16 -O bam -o ./SRR24332180.sort.bam ./SRR24332180.sam
[3] 118918

[uchiyamat@edo: ~/TEST/WORKSPACE/pipeline/01_mapping/route02/01_mapping]$ samtools index ./SRR24332178.sort.bam &
[1] 119239
[uchiyamat@edo: ~/TEST/WORKSPACE/pipeline/01_mapping/route02/01_mapping]$ samtools index ./SRR24332179.sort.bam &
[2] 119256
[uchiyamat@edo: ~/TEST/WORKSPACE/pipeline/01_mapping/route02/01_mapping]$ samtools index ./SRR24332180.sort.bam &
[3] 119264

[uchiyamat@edo: ~/TEST/WORKSPACE/pipeline/01_mapping/route02/01_mapping]$ rm -rf SRR24332178.sam SRR24332179.sam SRR24332180.sam 
```

### (2) Assemble transcripts with stringtie2
```bash
# Set up stringtie2
(base) [uchiyamat@edo: ~/TEST/tools]$ pwd
/export/home/uchiyamat/TEST/tools

# install stringtie2
(base) [uchiyamat@edo: ~/TEST/tools]$ git clone https://github.com/skovaka/stringtie2.git
Cloning into 'stringtie2'...
remote: Enumerating objects: 1669, done.
remote: Total 1669 (delta 0), reused 0 (delta 0), pack-reused 1669
Receiving objects: 100% (1669/1669), 5.75 MiB | 0 bytes/s, done.
Resolving deltas: 100% (896/896), done.
Checking connectivity... done.
Checking out files: 100% (529/529), done.
(base) (master)[uchiyamat@edo: ~/TEST/tools]$ cd stringtie2
(base) (master)[uchiyamat@edo: ~/TEST/tools/stringtie2]$ make release
(base) [uchiyamat@edo: ~/TEST/tools]$ stringtie
StringTie v2.0 usage:
 stringtie <input.bam ..> [-G <guide_gff>] [-l <label>] [-o <out_gtf>] [-p <cpus>]
  [-v] [-a <min_anchor_len>] [-m <min_tlen>] [-j <min_anchor_cov>] [-f <min_iso>]
  [-C <coverage_file_name>] [-c <min_bundle_cov>] [-g <bdist>] [-u] [-L]
  [-e] [-x <seqid,..>] [-A <gene_abund.out>] [-h] {-B | -b <dir_path>}
Assemble RNA-Seq alignments into potential transcripts.
```
```bash
[uchiyamat@edo: ~/TEST/WORKSPACE/pipeline/01_mapping/route02]$ mkdir -p 02_stringtie
[uchiyamat@edo: ~/TEST/WORKSPACE/pipeline/01_mapping/route02]$ cd 02_stringtie
```

#### Case_01
```bash
# Assemble transcripts
[uchiyamat@edo: ~/TEST/WORKSPACE/pipeline/01_mapping/route02/02_stringtie]$ stringtie -p 16 -G /export/home/uchiyamat/TEST/WORKSPACE/reference/gencode.v43.annotation.gtf -o ./SRR24035598.mormal.gtf /export/home/uchiyamat/TEST/WORKSPACE/pipeline/01_mapping/route02/01_mapping/SRR24035598.sort.bam

# Transcript merge, filtering (FPKM > 1 because cufflinks calculate the abundance of transcripts as FPKM) and reference annotation to include in the merging
[uchiyamat@edo: ~/TEST/WORKSPACE/pipeline/01_mapping/route02/02_stringtie]$ stringtie --merge -G /export/home/uchiyamat/TEST/WORKSPACE/reference/gencode.v43.annotation.gtf -F 1 SRR24035598.mormal.gtf -o SRR24035598.merge.gtf
```

#### Case_02
```bash
# Assemble transcripts
[uchiyamat@edo: ~/TEST/WORKSPACE/pipeline/01_mapping/route02/02_stringtie]$ stringtie -p 16 -G /export/home/uchiyamat/TEST/WORKSPACE/reference/gencode.v43.annotation.gtf -o ./SRR24332178.mormal.gtf /export/home/uchiyamat/TEST/WORKSPACE/pipeline/01_mapping/route02/01_mapping/SRR24332178.sort.bam &
[1] 120313
[uchiyamat@edo: ~/TEST/WORKSPACE/pipeline/01_mapping/route02/02_stringtie]$ stringtie -p 16 -G /export/home/uchiyamat/TEST/WORKSPACE/reference/gencode.v43.annotation.gtf -o ./SRR24332179.mormal.gtf /export/home/uchiyamat/TEST/WORKSPACE/pipeline/01_mapping/route02/01_mapping/SRR24332179.sort.bam &
[2] 120322
[uchiyamat@edo: ~/TEST/WORKSPACE/pipeline/01_mapping/route02/02_stringtie]$ stringtie -p 16 -G /export/home/uchiyamat/TEST/WORKSPACE/reference/gencode.v43.annotation.gtf -o ./SRR24332180.mormal.gtf /export/home/uchiyamat/TEST/WORKSPACE/pipeline/01_mapping/route02/01_mapping/SRR24332180.sort.bam &
[3] 120339

# Transcript merge, filtering (FPKM > 1 because cufflinks calculate the abundance of transcripts as FPKM) and reference annotation to include in the merging
[uchiyamat@edo: ~/TEST/WORKSPACE/pipeline/01_mapping/route02/02_stringtie]$ stringtie --merge -G /export/home/uchiyamat/TEST/WORKSPACE/reference/gencode.v43.annotation.gtf -F 1 SRR24332178.mormal.gtf SRR24332179.mormal.gtf SRR24332180.mormal.gtf -o SRR24332178_180.merge.gtf
```

## (3) Merge the gtf files from each route 
The original version of this program was distributed as part of the Cufflinks suite, under the name "CuffCompare". <br>
Another important difference is that the input transcripts are by default no longer discarded when they are found to be "intron redundant", i.e. contained within other, longer isoforms. CuffCompare had the -G option to prevent collapsing of such intron redundant isoforms into their longer "containers", but GffCompare has made this the default mode of operation (hence the -G option is no longer needed and is simply ignored when given). <br>

```bash
# set up gffcompare 
[uchiyamat@edo: ~/TEST/tools]$ wget http://ccb.jhu.edu/software/stringtie/dl/gffcompare-0.10.1.Linux_x86_64.tar.gz
[uchiyamat@edo: ~/TEST/tools]$ tar -zxvf gffcompare-0.10.1.Linux_x86_64.tar.gz
[uchiyamat@edo: ~/TEST/tools]$ gffcompare --version
gffcompare v0.10.1
```
```bash
[uchiyamat@edo: ~/TEST/WORKSPACE/pipeline]$ mkdir -p 02_gffcmp
[uchiyamat@edo: ~/TEST/WORKSPACE/pipeline]$ cd 02_gffcmp
```
```bash
# cheak
[uchiyamat@edo: ~/TEST/WORKSPACE/pipeline/02_gffcmp]$ wc -l /export/home/uchiyamat/TEST/WORKSPACE/pipeline/01_mapping/route01/02_cufflinks/cuffmerge/merged.gtf
2619623 /export/home/uchiyamat/TEST/WORKSPACE/pipeline/01_mapping/route01/02_cufflinks/cuffmerge/merged.gtf
[uchiyamat@edo: ~/TEST/WORKSPACE/pipeline/02_gffcmp]$ wc -l /export/home/uchiyamat/TEST/WORKSPACE/pipeline/01_mapping/route02/02_stringtie/SRR24332178_180.merge.gtf
1967177 /export/home/uchiyamat/TEST/WORKSPACE/pipeline/01_mapping/route02/02_stringtie/SRR24332178_180.merge.gtf
```

#### test_00 (cuffcompare vs gffcompare)
```bash
[uchiyamat@edo: ~/TEST/WORKSPACE/pipeline/02_gffcmp]$ mkdir -p test_00
[uchiyamat@edo: ~/TEST/WORKSPACE/pipeline/02_gffcmp]$ cd test_00/
```
Cufflinks and StringTie were reported to reconstruct many transcript assemblies with a single exon, and most of them are false positives. 
We removed those transcripts assemblies that only have a single exon to avoid having a large number of false positives in the candidates.
#### -> -M discard (ignore) single-exon transfrags and reference transcripts

#### cuffcompare
```bash
[uchiyamat@edo: ~/TEST/WORKSPACE/pipeline/02_gffcmp/test_00]$ cuffcompare \
> -M \
> -s /export/home/uchiyamat/TEST/WORKSPACE/reference/GRCh38.primary_assembly.genome.fa \
> -r /export/home/uchiyamat/TEST/WORKSPACE/pipeline/01_mapping/route01/02_cufflinks/cuffmerge/merged.gtf \
> -o ./cuffcompare 
> /export/home/uchiyamat/TEST/WORKSPACE/pipeline/01_mapping/route02/02_stringtie/SRR24332178_180.merge.gtf

Warning: Could not connect to update server to verify current version. Please check at the Cufflinks website (http://cufflinks.cbcb.umd.edu).
```

#### gffcompare 
```bash
[uchiyamat@edo: ~/TEST/WORKSPACE/pipeline/02_gffcmp/test_00]$ gffcompare \
> -M \
> -s /export/home/uchiyamat/TEST/WORKSPACE/reference/GRCh38.primary_assembly.genome.fa \
> -r /export/home/uchiyamat/TEST/WORKSPACE/pipeline/01_mapping/route01/02_cufflinks/cuffmerge/merged.gtf \
> -o ./gffcompare \
> /export/home/uchiyamat/TEST/WORKSPACE/pipeline/01_mapping/route02/02_stringtie/SRR24332178_180.merge.gtf

No fasta index found for /export/home/uchiyamat/TEST/WORKSPACE/reference/GRCh38.primary_assembly.genome.fa. Rebuilding, please wait..
Fasta index rebuilt.
```
```bash
[uchiyamat@edo: ~/TEST/WORKSPACE/pipeline/02_gffcmp]$ wc -l *gtf
  1634028 cuffcompare.combined.gtf
  1899931 gffcompare.annotated.gtf
```

#### test01 (Gencode vs Stringtie2)
```bash
[uchiyamat@edo: ~/TEST/WORKSPACE/pipeline/02_gffcmp]$ mkdir -p test_01
[uchiyamat@edo: ~/TEST/WORKSPACE/pipeline/02_gffcmp]$ cd test_01/

[uchiyamat@edo: ~/TEST/WORKSPACE/pipeline/02_gffcmp/test_01]$ gffcompare \
> -M \
> -s /export/home/uchiyamat/TEST/WORKSPACE/reference/GRCh38.primary_assembly.genome.fa \
> -r /export/home/uchiyamat/TEST/WORKSPACE/reference/gencode.v43.annotation.gtf \
> -o ./gffcompare \
> /export/home/uchiyamat/TEST/WORKSPACE/pipeline/01_mapping/route02/02_stringtie/SRR24332178_180.merge.gtf

[uchiyamat@edo: ~/TEST/WORKSPACE/pipeline/02_gffcmp/test_01]$ cuffcompare \
> -M \
> -s /export/home/uchiyamat/TEST/WORKSPACE/reference/GRCh38.primary_assembly.genome.fa \
> -r /export/home/uchiyamat/TEST/WORKSPACE/reference/gencode.v43.annotation.gtf \
> -o ./cuffcompare \
> /export/home/uchiyamat/TEST/WORKSPACE/pipeline/01_mapping/route02/02_stringtie/SRR24332178_180.merge.gtf

[uchiyamat@edo: ~/TEST/WORKSPACE/pipeline/02_gffcmp/test_01]$ wc -l *gtf
  1634028 cuffcompare.combined.gtf
  1899931 gffcompare.annotated.gtf
```

#### test02 (Gencode vs Cufflinks)
```bash
[uchiyamat@edo: ~/TEST/WORKSPACE/pipeline/02_gffcmp]$ mkdir -p test_02
[uchiyamat@edo: ~/TEST/WORKSPACE/pipeline/02_gffcmp]$ cd test_02/

[uchiyamat@edo: ~/TEST/WORKSPACE/pipeline/02_gffcmp/test_02]$ gffcompare \
> -M \
> -s /export/home/uchiyamat/TEST/WORKSPACE/reference/GRCh38.primary_assembly.genome.fa \
> -r /export/home/uchiyamat/TEST/WORKSPACE/reference/gencode.v43.annotation.gtf \
> -o ./gffcompare \
> /export/home/uchiyamat/TEST/WORKSPACE/pipeline/01_mapping/route01/02_cufflinks/cuffmerge/merged.gtf

[uchiyamat@edo: ~/TEST/WORKSPACE/pipeline/02_gffcmp/test_02]$ cuffcompare \
> -M \
> -s /export/home/uchiyamat/TEST/WORKSPACE/reference/GRCh38.primary_assembly.genome.fa \
> -r /export/home/uchiyamat/TEST/WORKSPACE/reference/gencode.v43.annotation.gtf \
> -o ./cuffcompare \
> /export/home/uchiyamat/TEST/WORKSPACE/pipeline/01_mapping/route01/02_cufflinks/cuffmerge/merged.gtf

[uchiyamat@edo: ~/TEST/WORKSPACE/pipeline/02_gffcmp/test_02]$ wc -l *gtf
  2096592 cuffcompare.combined.gtf
  2686035 gffcompare.annotated.gtf
```
#### -> The expressed RNA transcript assemblies will be compared with genome annotations using gffcompare, and each results will be merged.
```bash
# route01 ---
[uchiyamat@edo: ~/TEST/WORKSPACE/pipeline/02_gffcmp]$ gffcompare \
> -M \
> -s /export/home/uchiyamat/TEST/WORKSPACE/reference/GRCh38.primary_assembly.genome.fa \
> -r /export/home/uchiyamat/TEST/WORKSPACE/reference/gencode.v43.annotation.gtf \
> -o ./route01 \
> /export/home/uchiyamat/TEST/WORKSPACE/pipeline/01_mapping/route01/02_cufflinks/cuffmerge/merged.gtf

# route02 ---
[uchiyamat@edo: ~/TEST/WORKSPACE/pipeline/02_gffcmp]$ gffcompare \
> -M \
> -s /export/home/uchiyamat/TEST/WORKSPACE/reference/GRCh38.primary_assembly.genome.fa \
> -r /export/home/uchiyamat/TEST/WORKSPACE/reference/gencode.v43.annotation.gtf \
> -o ./route02 \
> /export/home/uchiyamat/TEST/WORKSPACE/pipeline/01_mapping/route02/02_stringtie/SRR24332178_180.merge.gtf

# cheak
[uchiyamat@edo: ~/TEST/WORKSPACE/pipeline/02_gffcmp]$ wc -l *gtf
  2686035 route01.annotated.gtf
  1899931 route02.annotated.gtf
```

## (4) filtering and merge
Only the transcript assemblies longer than 200 nt with Cuffcompare generated class codes of ‘U’, ‘I’, ‘O’, and ‘X’ were selected as lncRNA candidates. ‘U’ stands for unknown intergenic transcripts, ‘I’ indicates that the transcripts fall entirely within a reference intron, ‘O’ stands for generic exonic overlap with reference transcripts, and ‘X’ stands for exonic overlap with reference transcripts on the opposite strand (Supplementary Fig S10). 
The reason we selected transcripts with these codes is that in the GENCODE database most of the lncRNAs are located within the intergenic and protein-coding intronic regions, and the rest of the lncRNAs are overlapped with protein-coding exons on either the same or opposite strand. 
We also used the CD-HIT program to filter out transcript assemblies that exist in the NONCODE Human database and the transcripts from normal control with 80% identity43,66. We conducted filtration with the NONCODE database because it contains lncRNAs obtained from literature mining but may not be verified by GENCODE or NCBI RefSeq yet.



## Quantification of gene and transcripts using salmon
This is bacause both featureCount and rsem is implement in the current CogentAP pipline. featureCount and rsem is used for quantification of genes and transcripts each other. TPM values vary from tool to tool, so I want to use single tool to caluculate TPM for genes and transcripts sitimulately. <br>

Comparative evaluation of full-length isoform quantification from RNA-Seq　<br>
・Reference：[Comparative evaluation of full-length isoform quantification from RNA-Seq](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-021-04198-1) <br>
・salmon：https://combine-lab.github.io/salmon/ <br>

```bash
# Installing salmon via conda is possible, but the version is old.
(base) [uchiyamat@edo: ~/TEST/tools]$ conda install -c bioconda salmon
# Salmon v0.8.1

# install the avalible newest version binary mode (salmon v1.9.0)
(base) [uchiyamat@edo: ~/TEST/tools]$ pwd
/export/home/uchiyamat/TEST/tools
(base) [uchiyamat@edo: ~/TEST/tools]$ wget https://github.com/COMBINE-lab/salmon/releases/download/v1.9.0/salmon-1.9.0_linux_x86_64.tar.gz
(base) [uchiyamat@edo: ~/TEST/tools]$ tar -zxvf salmon-1.9.0_linux_x86_64.tar.gz
(base) [uchiyamat@edo: ~/TEST/tools]$ rm -rf salmon-1.9.0_linux_x86_64.tar.gz
(base) [uchiyamat@edo: ~/TEST/tools]$ salmon -h
salmon v1.9.0

Usage:  salmon -h|--help or
        salmon -v|--version or
        salmon -c|--cite or
        salmon [--no-version-check] <COMMAND> [-h | options]

Commands:
     index      : create a salmon index
     quant      : quantify a sample
     alevin     : single cell analysis
     swim       : perform super-secret operation
     quantmerge : merge multiple quantifications into a single file
```
```bash
  -g [ --geneMap ] arg                  File containing a mapping of
                                        transcripts to genes.  If this file is
                                        provided salmon will output both
                                        quant.sf and quant.genes.sf files,
                                        where the latter contains aggregated
                                        gene-level abundance estimates.  The
                                        transcript to gene mapping should be
                                        provided as either a GTF file, or a in
                                        a simple tab-delimited format where
                                        each line contains the name of a
                                        transcript and the gene to which it
                                        belongs separated by a tab.  The
                                        extension of the file is used to
                                        determine how the file should be
                                        parsed.  Files ending in '.gtf', '.gff'
                                        or '.gff3' are assumed to be in GTF
                                        format; files with any other extension
                                        are assumed to be in the simple format.
                                        In GTF / GFF format, the
                                        "transcript_id" is assumed to contain
                                        the transcript identifier and the
                                        "gene_id" is assumed to contain the
                                        corresponding gene identifier.
```



