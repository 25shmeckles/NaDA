import argparse, douwelib as dl, os, matplotlib, matplotlib.pyplot as plt, glob
from matplotlib.backends.backend_pdf import PdfPages

parser = argparse.ArgumentParser(description="Plot production of done_fastq files")
group = parser.add_mutually_exclusive_group()
group.add_argument("-pdf1", "--save_to_pdf1", action="store_true", help="Plot and Save files as pdf, 1 = directory")
group.add_argument("-pdf2", "--save_to_pdf2", action="store_true", help="Plot and Save files as pdf, 2 = file")
group.add_argument("-png1", "--save_to_png1", action="store_true", help="Plot and save files as png, 1 = directory")
group.add_argument("-png2", "--save_to_png2", action="store_true", help="Plot and save files as png, 2 = file")
parser.add_argument("path", help="path to folder with done_fastq files")
parser.add_argument("save_file_name", help="Name of the save file")
parser.add_argument("save_file_location", help="Save location path")
args = parser.parse_args()

if args.save_to_pdf1:
    pp = PdfPages('{}\{}.pdf'.format(args.save_file_location ,args.save_file_name))
    os.chdir(args.path)
    for file in glob.iglob('*.done_fastq'):
        fig = dl.plot_done_fastq_files(file)
        pp.savefig(fig)
    pp.close()
    
elif args.save_to_pdf2:
    pp = PdfPages('{}\{}.pdf'.format(args.save_file_location ,args.save_file_name))
    fig = dl.plot_done_fastq_files(args.path)
    pp.savefig(fig)
    pp.close()
    
elif args.save_to_png1:
    os.chdir(args.path)
    points = 1
    for file in glob.iglob('*.done_fastq'):
        fig = dl.plot_done_fastq_files(file)
        plt.savefig('{}\{}_{}.png'.format(args.save_file_location, args.save_file_name, points), bbox_inches='tight')
        points += 1
    plt.close()
    
elif args.save_to_png2:
    fig = dl.plot_done_fastq_files(args.path)
    plt.savefig('{}\{}.png'.format(args.save_file_location, args.save_file_name), bbox_inches='tight')
    plt.close()
    
else:
    print('no command found')