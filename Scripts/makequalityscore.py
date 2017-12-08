import argparse, douwelib as dl, os, matplotlib, matplotlib.pyplot as plt, glob
from matplotlib.backends.backend_pdf import PdfPages

parser = argparse.ArgumentParser(description="print bases with their qualityscores")
group = parser.add_mutually_exclusive_group()
group.add_argument("-c", "--checkscore", action="store_true", help="check base scores")
parser.add_argument("path", help="path to folder with done_fastq files")
parser.add_argument("size", type=int, help="Size of bases you want to check")
parser.add_argument("overlap", type=int, help="Overlap between bases !Note: can crash if put to high!")
parser.add_argument("startbase", type=int, help="At which bases you want to start, shifts read, starting from 0")
args = parser.parse_args()

if args.checkscore:
    data = dl.parse_fasta_file_error(args.path)
    id_ = list(data.keys())[0]
    error_rate = dl.convert_qualityscore (data[id_]['score'])
    nucleotide_string = (data[id_]['sequence'])
    print(dl.dict_scoremean_bases(nucleotide_string[args.startbase:len(nucleotide_string)], error_rate, args.size, args.overlap))

else:
    print('no command found')