import os
import re
import gzip
import pandas as pd
from pybedtools import BedTool

ref_genome = '/Users/tomoyauchiyama/lncRNA/reference/GRCh38.primary_assembly.genome.fa.gz'
gtf_file = '/Users/tomoyauchiyama/lncRNA/reference/LncBookv2.0_GENCODEv34_GRCh38.gtf.gz'

def get_exon_position(gtf_file):
    exon_info = {}
    with gzip.open(gtf_file, 'rt') as gtf:
        for line in gtf.readlines():
            line = line.rstrip('\n').split('\t')
            if '#' in line[0]:
                pass
            else:
                if line[2] == 'transcript':
                    flag = 0
                    for ele in line[-1].split('; '):
                        if 'transcript_id' in ele:
                            transcript_id = re.search(r'transcript_id \"(.*)\"', ele).group(1)
                        if 'gene_type' in ele:
                            gene_type = re.search(r'gene_type \"(.*)\"', ele).group(1)
                        if 'transcript_type' in ele:
                            transcript_type = re.search(r'transcript_type \"(.*)\"', ele).group(1)
                    if (gene_type == 'lncRNA') & (transcript_type == 'lncRNA'):
                        flag = 1
                elif line[2] == 'exon':
                    if flag == 1:
                        chrom = line[0]
                        strand = line[6]
                        start = line[3]
                        end = line[4]
                        exon = chrom + ':' + start + '-' + end
                        tx_id = transcript_id + ',' + strand
                        if exon_info.get(tx_id):
                            exon_info[tx_id] += ',' + exon
                        else:
                            exon_info[tx_id] = exon

    return exon_info

def get_splice_site(exon_info):
    keys = [
        'transcript_id', 'chr', 'exon1', 'exon2',
        'five_ss_start', 'five_ss_end', 'three_ss_start', 'three_ss_end'
    ]
    d = {key: [] for key in keys}

    for tx_id, exon in exon_info.items():
        transcript_id, strand = tx_id.split(',')
        ex_list = exon.split(',')
        exons = sorted(ex_list, key=lambda s: int(re.search(r'.*:(\d+)-\d+', s).group(1)))
        if len(exons) > 1:
            if strand == '+':
                exon1 = exons[0]
                exon2 = exons[1]
                exon1_info = re.search(r'(.*):\d+-(\d+)', exon1)
                chrom = exon1_info.group(1)
                pos = int(exon1_info.group(2))
                five_ss_start = str(pos - 3)
                five_ss_end = str(pos + 6)
                #five_ss = chrom + ':' + five_ss_start + '-' + five_ss_end
                #five_ss_seq = dna[chrom][int(five_ss_start)-1:int(five_ss_end)]

                exon2_info = re.search(r'(.*):(\d+)-\d+', exon2)
                chrom = exon2_info.group(1)
                pos = int(exon2_info.group(2))
                three_ss_start = str((pos - 1) - 20)
                three_ss_end = str((pos - 1) + 3)
                #three_ss = chrom + ':' + three_ss_start + '-' + three_ss_end
                #three_ss_seq = dna[chrom][int(three_ss_start)-1:int(three_ss_end)]

            elif strand == '-':
                exon1 = exons[-1]
                exon2 = exons[-2]
                exon1_info = re.search(r'(.*):(\d+)-\d+', exon1)
                chrom = exon1_info.group(1)
                pos = int(exon1_info.group(2))
                five_ss_start = str((pos - 1) - 6)
                five_ss_end = str((pos - 1) + 3)
                #five_ss = chrom + ':' + five_ss_start + '-' + five_ss_end
                #five_ss_seq = dna[chrom][int(five_ss_start)-1:int(five_ss_end)]
                #five_ss_seq = five_ss_seq[::-1]
                #five_ss_seq = five_ss_seq.translate(nucl.maketrans('ATGC', 'TACG'))

                exon2_info = re.search(r'(.*):\d+-(\d+)', exon2)
                chrom = exon2_info.group(1)
                pos = int(exon2_info.group(2))
                three_ss_start = str(pos - 3)
                three_ss_end = str(pos + 20)
                #three_ss = chrom + ':' + three_ss_start + '-' + three_ss_end
                #three_ss_seq = dna[chrom][int(three_ss_start)-1:int(three_ss_end)]
                #three_ss_seq = three_ss_seq[::-1]
                #three_ss_seq = three_ss_seq.translate(nucl.maketrans('ATGC', 'TACG'))

        append_data = {

            
            'transcript_id': transcript_id, 'chr': chrom, 'exon1': exon1, 'exon2': exon2,
            'five_ss_start': five_ss_start, 'five_ss_end': five_ss_end,
            'three_ss_start': three_ss_start, 'three_ss_end': three_ss_end
        }
        for key in keys:
            d[key].append(append_data[key])

    df = pd.DataFrame(d)
    return df


def run_maxentscan(df, ref_genome, outfa):
    # To work option s,  you need to set strand infomation on 6th column.
    b = BedTool.from_dataframe(df[['chr', 'five_ss_start', 'five_ss_end', 'transcript_id', 'strand']])
    #outfa = '/Users/tomoyauchiyama/lncRNA/reference/LncBookv2.0_lncRNA_5SS.fa'
    b.sequence(fi=ref_genome, s=True, fo=outfa)
    # $ perl score5.pl gencode_v43_lncRNA_5SS.fa > /Users/tomoyauchiyama/lncRNA/maxentpy/gencode.v43.lncRNA/gencode_v43_lncRNA_score5.txt
    # $ perl score3.pl gencode_v43_lncRNA_3SS.fa > /Users/tomoyauchiyama/lncRNA/maxentpy/gencode.v43.lncRNA/gencode_v43_lncRNA_score3.txt


    file = '/Users/tomoyauchiyama/lncRNA/maxentpy/gencode.v43.lncRNA/gencode_v43_lncRNA_score5.txt'
    df5 = pd.read_csv(file, sep='\t', header=None)
    df_split = df5[0].str.split(' ', expand=True)
    df5 = pd.concat([df_split, df5[[1, 2]]], axis=1)
    df5.columns = ['pos', 'ID' , '5ss_seq', '5ss_score']

    file = '/Users/tomoyauchiyama/lncRNA/maxentpy/gencode.v43.lncRNA/gencode_v43_lncRNA_score3.txt'
    df3 = pd.read_csv(file, sep='\t', header=None)
    df_split = df3[0].str.split(' ', expand=True)
    df3 = pd.concat([df_split, df3[[1, 2]]], axis=1)
    df3.columns = ['pos', 'ID' , '3ss_seq', '3ss_score']

    df_merge = pd.merge(df5, df3, how='left', on='ID')

    df_ID = df_merge['ID'].str.split(',', expand=True)
    df_ID.columns = ['transcript_id', 'strand']
    df_ss_score = pd.concat([df_ID, df_merge.drop(columns='ID')], axis=1)

    return df_ss_score

