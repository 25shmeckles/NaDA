import argparse, subprocess

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Fastq qualityscore analyzer")
    parser.add_argument("path", help="path of files")
    parser.add_argument("filename", help="output name")
    parser.add_argument("save_path", help="path to save file")
    parser.add_argument("sequence_size", help="select how many bases sequence analysis should be")
    parser.add_argument("overlap", help="select how many bases overlap between sequence shpuld be")
    args = parser.parse_args()
    
    #python3 = '/hpc/local/CentOS7/common/lang/python/3.6.1/bin/python'
    nada_root = '/hpc/cog_bioinf/kloosterman/users/dspaanderman/NaDA'
    
    #Cyclomics pipeline
    cmd0 = f'. {nada_root}/env/bin/activate'
    cmd1 = f'python {nada_root}/Scripts/fastq_qualityscore_analyser.py {args.path} {args.filename} {args.save_path} {args.sequence_size} {args.overlap}'
    cmd2 = 'deactivate'
    
    echo_cmd = f'echo "{cmd0} && {cmd1} && {cmd2}"'
    qsub_cmd = f'qsub -l h_rt=24:00:00,h_vmem=32G -pe threaded 8 -cwd -M dspaanderman@umcutrecht.nl -m beas -N fastq'
    
    subprocess.run( f'{echo_cmd} | {qsub_cmd}' , shell=True)