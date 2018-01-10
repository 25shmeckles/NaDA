import argparse, subprocess

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Vcf SNP analyzer")
    group = parser.add_mutually_exclusive_group()
    group.add_argument("-b", "--BB200", action="store_true", help="Only files with BB200 will be checked")
    group.add_argument("-a", "--All_B", action="store_true", help="All backbones will be checked")
    parser.add_argument("path", help="path of files")
    parser.add_argument("filename", help="output name")
    parser.add_argument("save_path", help="path to save file")
    parser.add_argument("sequence_size", help="select how many bases sequence analysis should be")
    args = parser.parse_args()
    
    if args.BB200:
        BB = '--BB200'
    elif args.All_B:
        BB = '--All_B'
    else:
        BB = ''

    python3 = '/hpc/local/CentOS7/common/lang/python/3.6.1/bin/python'
    nada_root = '/hpc/cog_bioinf/kloosterman/users/dspaanderman/NaDA'
    
    #Cyclomics pipeline
    cmd0 = f'module load python/3.6.1 && . {nada_root}/env/bin/activate'
    cmd1 = f'{python3} {nada_root}/Script/vcf_snp_variant_analyser.py {BB} {args.path} {args.filename} {args.safe_path} {args.sequence_size}'
    cmd2 = 'deactivate'
    
    echo_cmd = f'echo "{cmd0} && {cmd1} && {cmd2}"'
    qsub_cmd = f'qsub -l h_rt=24:00:00,h_vmem=64G -pe threaded 8 -cwd -M dspaanderman@umcutrecht.nl -m beas -N NaDA'
    
    subprocess.run( f'{echo_cmd} | {qsub_cmd}' , shell=True)