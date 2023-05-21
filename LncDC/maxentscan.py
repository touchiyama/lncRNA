import os
import re
import sys
import gzip
import datetime
import argparse
import subprocess
import pandas as pd
from pybedtools import BedTool
from logging import Logger, getLogger, Formatter, StreamHandler, FileHandler

DATE = datetime.datetime.now().strftime("%Y_%m_%d_%H_%M_%S")
LOGFILE = os.path.join(
    os.path.abspath('.'),
    os.path.basename(sys.argv[0]).replace('.py', '') + '_' + DATE + '.log'
)

# ref_genome = '/Users/tomoyauchiyama/lncRNA/reference/GRCh38.primary_assembly.genome.fa.gz'
# gtf_file = '/Users/tomoyauchiyama/lncRNA/reference/LncBookv2.0_GENCODEv34_GRCh38.gtf.gz'

def command(cmd, outfile):
    fo = open(outfile, 'w')
    try:
        subprocess.run(
            cmd, check=True, stdout=fo,
            stderr=subprocess.PIPE, universal_newlines=True
        )
    except subprocess.CalledProcessError as e:
        logger.error(f'Failed: {e}')
        logger.error('Exit..')
        sys.exit()

def get_exon_position(gtf_file):
    logger.info('')
    logger.info('Start getting exon position from gtf file...')

    if gtf_file.endswith((".gz", ".Z", ".z")):
        fd = gzip.open(gtf_file, 'rt')
    else:
        fd = open(gtf_file, 'r')

    exon_info = {}
    with fd as gtf:
        for line in gtf.readlines():
            line = line.rstrip('\n').split('\t')
            if '#' in line[0]:
                pass
            else:
                if line[2] == 'transcript':
                    for ele in line[-1].split('; '):
                        if 'gene_id' in ele:
                            gene_id = re.search(r'gene_id \"(.*)\"', ele).group(1)
                        if 'transcript_id' in ele:
                            transcript_id = re.search(r'transcript_id \"(.*)\"', ele).group(1)
                elif line[2] == 'exon':
                    chrom = line[0]
                    strand = line[6]
                    start = line[3]
                    end = line[4]
                    exon = chrom + ':' + start + '-' + end
                    tx_id = gene_id + '//' + transcript_id + '//' + strand
                    if exon_info.get(tx_id):
                        exon_info[tx_id] += '//' + exon
                    else:
                        exon_info[tx_id] = exon

    logger.info('Finish getting exon position from gtf file...')

    return exon_info

def get_splice_site(exon_info):
    logger.info('')
    logger.info('Start getting splice site position...')

    keys = [
        'gene_id', 'transcript_id', 'chr', 'strand', 'exon1', 'exon2',
        'five_ss_start', 'five_ss_end', 'three_ss_start', 'three_ss_end'
    ]
    d = {key: [] for key in keys}

    for tx_id, exon in exon_info.items():
        gene_id, transcript_id, strand = tx_id.split('//')
        ex_list = exon.split('//')
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

                exon2_info = re.search(r'(.*):(\d+)-\d+', exon2)
                chrom = exon2_info.group(1)
                pos = int(exon2_info.group(2))
                three_ss_start = str((pos - 1) - 20)
                three_ss_end = str((pos - 1) + 3)

            elif strand == '-':
                exon1 = exons[-1]
                exon2 = exons[-2]

                exon1_info = re.search(r'(.*):(\d+)-\d+', exon1)
                chrom = exon1_info.group(1)
                pos = int(exon1_info.group(2))
                five_ss_start = str((pos - 1) - 6)
                five_ss_end = str((pos - 1) + 3)

                exon2_info = re.search(r'(.*):\d+-(\d+)', exon2)
                chrom = exon2_info.group(1)
                pos = int(exon2_info.group(2))
                three_ss_start = str(pos - 3)
                three_ss_end = str(pos + 20)

        append_data = {
            'gene_id': gene_id, 'transcript_id': transcript_id, 'chr': chrom,
            'strand': strand, 'exon1': exon1, 'exon2': exon2,
            'five_ss_start': five_ss_start, 'five_ss_end': five_ss_end,
            'three_ss_start': three_ss_start, 'three_ss_end': three_ss_end
        }
        for key in keys:
            d[key].append(append_data[key])
    df = pd.DataFrame(d)

    logger.info('Finish getting splice site position ...')

    return df

def prep_maxentscan(df, ref_genome, outdir, prefix):
    """_summary_

    Args:
        df (_type_): _description_
        ref_genome (_type_): _description_
        outfa (_type_): _description_
        # outfa = '/Users/tomoyauchiyama/lncRNA/reference/LncBookv2.0_lncRNA_5SS.fa'
        # To work option s, you need to set strand infomation on 6th column.
    """

    logger.info('')
    logger.info('Start preparing for maxentscan...')

    five_ss_outfa = os.path.join(outdir, prefix + '_score5.fa')
    b = BedTool.from_dataframe(df[['chr', 'five_ss_start', 'five_ss_end', 'transcript_id', 'gene_id', 'strand']])
    b.sequence(fi=ref_genome, s=True, name=True, fo=five_ss_outfa)

    three_ss_outfa = os.path.join(outdir, prefix + '_score3.fa')
    b = BedTool.from_dataframe(df[['chr', 'three_ss_start', 'three_ss_end', 'transcript_id', 'gene_id', 'strand']])
    b.sequence(fi=ref_genome, s=True, name=True, fo=three_ss_outfa)

    logger.info('Finish preparing for maxentscan...')

    return five_ss_outfa, three_ss_outfa

def run_maxentscan(five_ss_outfa, three_ss_outfa, maxentscan_dir):
    logger.info('')
    logger.info('Start running maxentscan...')

    five_dirname = os.path.dirname(five_ss_outfa)
    five_ss_name = os.path.splitext(os.path.basename(five_ss_outfa))[0]
    five_ss_file = os.path.join(five_dirname, five_ss_name + '.txt')
    five_score = os.path.join(maxentscan_dir, 'score5.pl')

    three_dirname = os.path.dirname(three_ss_outfa)
    three_ss_name = os.path.splitext(os.path.basename(three_ss_outfa))[0]
    three_ss_file = os.path.join(three_dirname, three_ss_name + '.txt')
    three_score = os.path.join(maxentscan_dir, 'score3.pl')

    # shutil.copy(src, dst, *, follow_symlinks=True)
    # os.symlink(src, dst)

    #me2x5 = os.path.join(maxentscan_dir, 'me2x5')
    #logger.info(f'cp {me2x5} {five_dirname}')

    #splicemodels = os.path.join(maxentscan_dir, 'splicemodels')
    #logger.info(f'cp -r {splicemodels} {five_dirname}')

    logger.info(f'perl {five_score} {five_ss_outfa} > {five_ss_file}')
    cmd = ['/usr/local/bin/perl', five_score, five_ss_outfa, '>', five_ss_file]
    command(cmd, five_ss_file)

    logger.info(f'perl {three_score} {three_ss_outfa} > {three_ss_file}')
    cmd = ['/usr/local/bin/perl', three_score, three_ss_outfa, '>', three_ss_file]
    command(cmd, three_ss_file)

    logger.info('Finish running maxentscan...')

def parameters(__desc__):
    """入力されたコマンドライン引数をパースする関数
    Args:
        __desc__ (str): usege文

    Returns:
        args (argparse.Namespace): コマンドライン引数パース結果
    """
    parser = argparse.ArgumentParser(
        description = __desc__,
        #formatter_class=argparse.ArgumentDefaultsHelpFormatter #show default
    )
    parser.add_argument(
        '-g',
        dest='gtf_file',
        help='input genome annotation file (GTF format)',
        type=str,
        required=True
    )
    parser.add_argument(
        '-s',
        dest='seq_file',
        help='input genome sequence file (FASTA format) ',
        type=str,
        required=True
    )
    parser.add_argument(
        '-o',
        dest='out_dir',
        help='input sample name based on the input data header',
        type=str,
        default=os.path.abspath('.')
    )
    parser.add_argument(
        '-p',
        dest='prefix',
        help='input output path (absolution path)',
        type=str,
        default=os.path.basename(sys.argv[0]).replace('.py', '')
    )
    parser.add_argument(
        '-t',
        dest='maxentscan_dir',
        help='input maxentscan_path (absolute path)',
        type=str,
        default='/Users/tomoyauchiyama/lncRNA/maxentscan'
    )
    parser.add_argument(
        '-log',
        dest='loglevel',
        help='choose log level (default: INFO)',
        choices=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'],
        default='INFO'
    )
    args = parser.parse_args()

    return args

def show_param(args: argparse.Namespace):
    """入力として与えられたコマンドライン引数を表示

    Args:
        args (argparse.Namespace): コマンドライン引数パース結果
    """
    logger.info('')
    logger.info('#--------------------------------------------------------------')
    logger.info(f"program: {__file__}")
    cla = vars(args)
    for key in sorted(cla.keys()):
        logger.info(f"{key}: {cla[key]}")
    logger.info("#--------------------------------------------------------------")

def main(args):
    logger.info('')
    logger.info('Output log...')

    # setUp parameter ---
    show_param(args)

    gtf_file = args.gtf_file
    ref_genome = args.seq_file
    outdir = args.out_dir
    prefix = args.prefix
    maxentscan_dir = args.maxentscan_dir

    # run process ---
    exon_info = get_exon_position(gtf_file)
    df = get_splice_site(exon_info)
    five_ss_outfa, three_ss_outfa = prep_maxentscan(df, ref_genome, outdir, prefix)
    run_maxentscan(five_ss_outfa, three_ss_outfa, maxentscan_dir)

    logger.info('Done!')
    logger.info('')

    # move file ---
    os.system(f'mv {LOGFILE} {outdir}')

if __name__ == '__main__':
    __version__ = '1.0'
    __desciption__ = 'Some useful program commands are:'

    parm = parameters(__desciption__)

    logger = getLogger(__name__)
    logger.setLevel(parm.loglevel)

    FORMAT = '%(levelname)s:[%(asctime)s] %(message)s'
    #dt_fmt = '%Y-%m-%d %H:%M:%S'
    formatter = Formatter(FORMAT)

    stream_handler = StreamHandler()
    stream_handler.setLevel(parm.loglevel)
    stream_handler.setFormatter(formatter)

    file_handler = FileHandler(filename=LOGFILE, mode='w', encoding='utf-8')
    file_handler.setLevel(parm.loglevel)
    file_handler.setFormatter(formatter)

    logger.addHandler(stream_handler)
    logger.addHandler(file_handler)

    main(parm)

#
# END
#