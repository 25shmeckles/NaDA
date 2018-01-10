#!/hpc/local/CentOS7/common/lang/python/3.6.1/bin/python3
#$ -l h_rt=48:00:00,h_vmem=64G
#$ -pe threaded 8
#$ -S /hpc/local/CentOS7/common/lang/python/3.6.1/bin/python3
#$ -cwd
#$ -M d.j.spaanderman@umcutrecht.nl
#$ -m beas

import argparse, glob, os, pandas as pd, statistics, itertools, bokeh.palettes as bp, numpy as np
from bokeh.plotting import figure, output_file, save, ColumnDataSource
from bokeh.layouts import gridplot
from bokeh.models import NumeralTickFormatter, HoverTool, GlyphRenderer, Range1d, LinearColorMapper, BasicTicker, PrintfTickFormatter, ColorBar, ColumnDataSource
from collections import defaultdict, Counter
from math import pi

def data_vcf_file(file_name, backbone_name):
    '''strips data from vcf files
    
    '''
    data = {}
    with open(file_name, 'r') as f:
        loglist = f.readlines()
        found = False
        for lines in loglist:
            if ('{}'.format(backbone_name)) in lines:
                found = True
                
        f.seek(0)
        if found is True:
            id_ = False
            backbone = []
            variance = []
            for line in f:
                if line.startswith('##'):
                    continue
                elif line.startswith('#'):
                        if not id_: #if id_ == False
                            id_ = line.strip()
                        else:
                            print('WARNING: An ID was found without corresponding sequence.', id_)
                            id_ = line.strip()
                elif line.startswith('B'):
                    if './.' in line:
                        continue
                    else:
                        backbone.append(line.strip())
                elif id_:
                    if './.' in line:
                        continue
                    else:
                        variance.append(line.strip())
                        data[id_] = {'variance':variance, 
                                'backbone':backbone}
        else:
            print('{} has other backbone than input'.format(file_name))
    return data

def vcf_whole_sequence_strip(file_name):
    '''get all the sequence data from vcf file
    input is file name.
    
    '''
    data = []
    with open(file_name, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            else:
                data.append(line.strip())
    return data

def filter_score(s):
    ''' filter s if not len >= 10
    
    '''
    if len(s) >= 10:
        return True
    return False

def mutated_reads_vcf_only(variance_or_backbone, data_all, size):
    '''mutations in vcf file with 3 sequences before and
    3 sequences after mutation. Overlap can occur if mutations
    are found beside each other. if location in sequence has
    multiple SNPs like A -> G/T. It will check if score for both
    these SNPs is high enough and will return as A -> G. So
    same location sequence can occur twice in return if A -> G
    and A -> T have both high enough scores
    
    Update 10-1: size can be determined and so size of sequence
    changes
    
    Args:
        either variance or backbone data can be entered and
        data_all which is the strip of all the sequencing data
        in the same vcf file
    
    '''
    
    
    if size == 3:
        lenght_before = 1
        lenght_after = 2
        
    elif size == 5:
        lenght_before = 2
        lenght_after = 3
    
    elif size == 7:
        lenght_before = 3
        lenght_after = 4
        
    else:
        lenght_before = 2
        lenght_after = 3

    score_ = []
    score = []
    for items in variance_or_backbone:
        breakpoint = '\t'
        score_.append(items.split(breakpoint, 9)[9])

    for i in filter(filter_score, score_):
        score.append(i)

    score2 = []
    for s in score:
        score2.append(s.split(':')[1].split(','))

    score3 = []
    for s2 in score:
        score3.append(s2.split(':')[2]+':'+s2.split(':')[3]+':'+s2.split(':')[4])
    
    points = 0
    highmutated = []
    extended = []
    for item in score2:
        if len(item) > 2:
            n1 = int(item[0])
            n2 = int(item[1])
            n3 = int(item[2])
            if n2/(n1+n2+n3) > 0.25:
                for i, items in enumerate(data_all):
                    mutated = ':'+','.join(item[0:3])+':'+score3[points]
                    if mutated in items:
                        r_ = items.split(breakpoint, 5)
                        r = r_[4][0]
                        removed = r_[0]+'\t'+r_[1]+'\t'+r_[2]+'\t'+r_[3]+'\t'+r+'\t'+r_[5]
                        highmutated.append(removed)
                        extended.append(data_all[i-lenght_before:i+lenght_after])
                        continue
            if n3/(n1+n2+n3) > 0.25:
                for i, items in enumerate(data_all):
                    mutated = ':'+','.join(item[0:3])+':'+score3[points]
                    if mutated in items:
                        r_ = items.split(breakpoint, 5)
                        r = r_[4][0]
                        removed = r_[0]+'\t'+r_[1]+'\t'+r_[2]+'\t'+r_[3]+'\t'+r+'\t'+r_[5]
                        highmutated.append(removed)
                        extended.append(data_all[i-lenght_before:i+lenght_after])
        else:
            n1 = int(item[0])
            n2 = int(item[1])
            if n2/(n1+n2) > 0.25:
                for i, items in enumerate(data_all):
                    mutated = ':'+','.join(item[0:2])+':'+score3[points]
                    if mutated in items:
                        highmutated.append(items)
                        extended.append(data_all[i-lenght_before:i+lenght_after])
        points += 1
    points = False
    return extended, highmutated

def list_with_all_combinations(letters, size):
    '''list of all combination letters
    
    '''
    all_combinations = []
    for i in itertools.product('ACTG', repeat = size):
        size_ = '{}'*size
        all_combinations.append(size_.format(*i))
    
    return all_combinations
       
def dictionary_sequence_counter(variance_or_backbone_data, size):
    '''Put in is sequence derived from vcf files. These can be
    from either variance or backbone. Than counts how many times
    all_combinations can be found with overlap in all these vcf
    files and put it into a dictionary. If item is not found in 
    sequence, it will not appear in the dictionary
    {AAAAA: 5}
    
    '''
    x = '-----'.join(variance_or_backbone_data)
    dictionary_sequence = {}
    all_combinations = list_with_all_combinations('ACTG', int(size))

    for item in all_combinations:
        if item in x:
            dictionary_sequence[item] = occurrences(x, item)
        else:
            continue
            
    return dictionary_sequence

def occurrences(string, sub):
    '''counts times something occurs in string
    with sub being what you are looking for.
    This is with overlap!
    
    '''
    count = start = 0
    while True:
        start = string.find(sub, start) + 1
        if start > 0:
            count+=1
        else:
            return count

def vcf_heatmap_snps(data_surrounding, data_variance, size):
    '''data generation for heatmap plot
    creates dictionary with key is surrounding sequence
    and types of variances:
    {'AAAAA':['G','T']}
    where G, T is variance for third A
    
    Args:
        data_surrounding = vcf_all_strip[1]
        data_variance = vcf_all_strip[3]
    
    '''
    sequence =[]
    mutation = []
    points = 0
    for items in data_surrounding:
        for i in items:
            sequence_ = ([j.split('\t')[3] for j in i])
            sequence_2 = [''.join(sequence_)]
            for s in sequence_2:
                sequence.append(s)

    for mut in data_variance:
        mutation_ = ([m.split('\t')[4] for m in mut])
        if not mutation_:
            continue
        else:
            for m in mutation_:
                mutation.append(m) 
    
    dict_variance = {}
    points = 0
    for se in sequence:
        if len(se) == size:
            if se in dict_variance:
                dict_variance[se].append(mutation[points])
            else:
                dict_variance[se] = [mutation[points]]
            points += 1
        else:
            continue
        
    return dict_variance

def pd_df_heatmap_sequence(data_dict, variance_or_backbone_data, size, v_or_b):
    '''From dictionary makes panda's dataframe.
    Firstly makes dictionary with bases to 0
    values and inputs data_dict into that dictionary
    while counting it.
    Next a df is build from this dictionary with index
    is keys from data_dict
    
    Update 10-1, uses dictionary of either insert or 
    backbone to calculate percentage of sequences mutated
    
    Example
    {'AAAA':[A, G, G]} -> {'A':1,'T':0,'C':0,'G':2} ->
    Dataframe with sequence followed by 1 0 0 2
    
    '''
    x = list(data_dict.keys())
    d = defaultdict(list)
    for k in 'ATCG':
        d[k]= 0
    
    data_dictionary = []
    points = 0
    for i in data_dict.values():
        c = Counter(i)
        z = {**d, **c}
        data_dictionary.append(z)
        
    df = pd.DataFrame(data=data_dictionary, index=x)
    df.columns.name = 'Bases'
    df.index.name = 'Sequences'
    
    #update 9-1
    p = pd.DataFrame.from_dict(dictionary_sequence_counter(variance_or_backbone_data, size), orient='index')
        
    try: 
        df['A'] = (df['A']/p[0])*100
        df['C'] = (df['C']/p[0])*100
        df['G'] = (df['G']/p[0])*100
        df['T'] = (df['T']/p[0])*100
        
    except KeyError:
        print("KeyError occurred dataframe for {} hasn't changed to percentages".format(v_or_b))
    
    #print('{} from {}'.format(df, v_or_b))
    return df

def pd_df_heatmap_variance(data):
    '''From variance data of vcf files makes
    pandas dataframe. this dataframe has
    index is always 1. with column has all
    possible variances with occurence of these
    types of variances
      
    '''
    x = list(dict(highmutated_back_variance(data)).keys())
    y = list(dict(highmutated_back_variance(data)).values())

    data_dict = dict(zip(x, y))
    df = pd.DataFrame(data=data_dict, index = [1])
    df.columns.name = 'mutations'
    df.index.name = 'index_'
    
    return df

def highmutated_back_variance(variance_or_backbone_highmutated):
    '''enter highmutated sequence from vcf file into here to count
    which snp occurs often.
    
    Args:
        variance_or_backbone_highmutated:
        [['... \t A \t C'][][]['.. \t C \t T ...']
        
    Returns:
        {A -> C : 1, C -> T: 1}
    
    '''
    mut = False
    mut = []
    for item in variance_or_backbone_highmutated:
        for i in item:
            mut.append((i.split('\t')[3]+'->'+i.split('\t')[4]))

    mutated_counter = Counter(mut)
    return mutated_counter

#plots
def plot_vcf_snps(data, output_name, save_path, v_or_b):
    '''plot single nucleotide polymorphisms
    data should be either mutated sequences from
    backbone or variance
    
    '''
    x = list(dict(highmutated_back_variance(data)).keys())
    y_ = list(dict(highmutated_back_variance(data)).values())
    y = []
    for item in y_:
        y.append(item/sum(y_))

    p = figure(x_range=x, plot_width=400, plot_height=400, title='Distribution of mutations in {}'.format(v_or_b), background_fill_color="#E8DDCB")
    p.vbar(x, width=0.5, bottom=0, top=y, fill_color="#036564", line_color="#033649")

    p.yaxis[0].formatter = NumeralTickFormatter(format="0.0%")
    p.xaxis.major_label_orientation = pi/4
    p.legend.location = "center_right"
    p.legend.background_fill_color = "darkgrey"
    p.xaxis.axis_label = 'Single Nucleotide Polymorphism'
    p.yaxis.axis_label = 'amount of mutations'
    
    output_file("{}/{}_{}_SNP_plot.html".format(save_path, output_name, v_or_b))
    save(p)

def heatmap_vcf_files_snps_with_sequence(df_, output_name, save_path, v_or_b):
    '''Heatmap of single nucleotide polymorphisms 
    from vcf files with surrounding sequence (4 bases
    extra)
    
    '''
    if df_.empty == True:
        print('DataFrame of {} is empty and no plot will be made'.format(v_or_b))
        
    else:
        df = pd.DataFrame(df_.stack(), columns=['scores']).reset_index()

        bases = list(df_.columns)
        sequences = list(df_.index)

        colors = ["#75968f", "#a5bab7", "#c9d9d3", "#e2e2e2", "#dfccce", "#ddb7b1", "#cc7878", "#933b41",
                  "#550b1d"]
        mapper = LinearColorMapper(palette=colors, low=df.scores.min(), high=df.scores.max())

        source = ColumnDataSource(df)

        TOOLS = "hover,reset,xpan,xwheel_zoom"

        p = figure(title='Variant occurence in {}'.format(v_or_b), x_range=sequences,
                   y_range=list(reversed(bases)), x_axis_location='above', plot_width=900, plot_height=400,
                   tools=TOOLS, toolbar_location='below')
        p.grid.grid_line_color = None
        p.axis.axis_line_color = None
        p.axis.major_tick_line_color = None
        p.axis.major_label_text_font_size = "5pt"
        p.axis.major_label_standoff = 0
        p.xaxis.major_label_orientation = pi / 3


        p.rect(y="Bases", x="Sequences", width=1, height=1,
               source=source,
               fill_color={'field': 'scores', 'transform': mapper},
               line_color=None)

        color_bar = ColorBar(color_mapper=mapper, major_label_text_font_size="5pt",
                             ticker=BasicTicker(desired_num_ticks=len(colors)),
                             label_standoff=6, border_line_color=None, location=(0, 0))
        p.add_layout(color_bar, 'right')

        p.select_one(HoverTool).tooltips = [
             ('mutation', '@Sequences -> @Bases'),
             ('occurence', '@scores'),
        ]

        output_file("{}/{}_{}_heatmap_sequences.html".format(save_path, output_name, v_or_b))
        save(p)

def heatmap_vcf_files_snps(df_, output_name, save_path, v_or_b):
    '''Heatmap of single nucleotide polymorphisms 
    from vcf files
    
    '''
    df = pd.DataFrame(df_.stack(), columns=['scores']).reset_index()

    mutations = list(df_.columns)
    index = list(df_.index)

    colors = bp.all_palettes['Viridis'][11]
    mapper = LinearColorMapper(palette=colors, low=0, high=df['scores'].sum()/4)

    source = ColumnDataSource(df)

    TOOLS = "hover,reset"

    p = figure(title='Variant occurence in {}'.format(v_or_b), y_range=list(reversed(mutations)), 
               plot_width=300, x_axis_location='above', plot_height=600, tools=TOOLS, toolbar_location='below')
    p.grid.grid_line_color = None
    p.axis.axis_line_color = None
    p.axis.major_tick_line_color = None
    p.axis.major_label_text_font_size = "5pt"
    p.axis.major_label_standoff = 0
    p.xaxis.ticker = []


    p.rect(y="mutations", x="index_", width=1, height=1,
           source=source,
           fill_color={'field': 'scores', 'transform': mapper},
           line_color=None)

    color_bar = ColorBar(color_mapper=mapper, major_label_text_font_size="5pt",
                         ticker=BasicTicker(desired_num_ticks=1),
                         label_standoff=5, border_line_color=None, height=50, location=(10, -250))
    p.add_layout(color_bar, 'right')

    p.select_one(HoverTool).tooltips = [
         ('mutation', '@mutations'),
         ('occurence', '@scores'),
    ]
    
    output_file("{}/{}_{}_heatmap_SNPs.html".format(save_path, output_name, v_or_b))
    save(p)

#Actual script for running files from dir and subdir
def main(input_folder, output_name, save_path, backbone_name, size):
    #recursive=True
    variance_data = []
    variance_sequence = []
    highmutated_v = []
    backbone_data = []
    backbone_sequence = []
    highmutated_b = []
    points = 1
    points_list = np.arange(500, 500000000, 500)
    
    print('getting inputs')
    for subdir, dirs, files in os.walk(input_folder):
        print('{} files are being processed'.format(len(files)))
        for filename in files:
            if points in points_list:
                print('file {} of {}'.format(points, len(files)))
            points += 1
            if filename.split('.')[-1] == 'vcf':    
                file = os.path.join(input_folder, subdir, filename)
                data = data_vcf_file(file, backbone_name)
                data_all = vcf_whole_sequence_strip(file)
                if bool(data) == True :
                    id_ = list(data.keys())[0]
                    variance = data[id_]['variance']
                    backbone = data[id_]['backbone']
                    
                    #for insert the Data
                    variance_list = []
                    for item in variance:
                        variance_list.append(item.split('\t',4)[3])
                    variance_sequence.append(''.join(variance_list))

                    variance_data.append(mutated_reads_vcf_only(variance, data_all, size)[0])
                    highmutated_v.append(mutated_reads_vcf_only(variance, data_all, size)[1])
                    
                    #for backbone the Data
                    backbone_list = []
                    for item in backbone:
                        backbone_list.append(item.split('\t',4)[3])
                    backbone_sequence.append(''.join(backbone_list))
                    
                    backbone_data.append(mutated_reads_vcf_only(backbone, data_all, size)[0])
                    highmutated_b.append(mutated_reads_vcf_only(backbone, data_all, size)[1])
            
    #plot SNPs
    print('plotting SNPs insert')
    plot_vcf_snps(highmutated_v, output_name, save_path, 'insert')
    print('plotting SNPs backbone')
    plot_vcf_snps(highmutated_b, output_name, save_path, 'backbone')
    
    #plot Heatmap_sequence
    print('plotting Heatmap insert')
    data_dict = vcf_heatmap_snps(variance_data, highmutated_v, size)
    df_ = pd_df_heatmap_sequence(data_dict, variance_sequence, size, 'insert')
    heatmap_vcf_files_snps_with_sequence(df_, output_name, save_path, 'insert')
    
    print('plotting Heatmap backbone')
    data_dict = vcf_heatmap_snps(backbone_data, highmutated_b, size)
    df_ = pd_df_heatmap_sequence(data_dict, backbone_sequence, size, 'backbone')
    heatmap_vcf_files_snps_with_sequence(df_, output_name, save_path, 'backbone')
        
    #plot SNPs_heatmap
    print('plotting SNPs_heatmap insert')
    df_ = pd_df_heatmap_variance(highmutated_v)
    heatmap_vcf_files_snps(df_, output_name, save_path, 'insert')
    
    print('plotting SNPs_heatmap insert')
    df_ = pd_df_heatmap_variance(highmutated_b)
    heatmap_vcf_files_snps(df_, output_name, save_path, 'backbone')

    

if __name__ == '__main__':

    #Argument script
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
        backbone_name = 'BB200'
    
    elif args.All_B:
        backbone_name = 'BB' 

    print('selected backbone: {}'.format(backbone_name))
    
    input_folder = args.path
    output_name = args.filename
    save_path = args.save_path
    size = float(args.sequence_size)

    print('started')
    main(input_folder, output_name, save_path, backbone_name, size)
    print('completed')