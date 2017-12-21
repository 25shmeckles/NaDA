##-IMPORT-##
import matplotlib, numpy as np, matplotlib.mlab as mlab, matplotlib.pyplot as plt, glob, os, statistics, itertools, pandas as pd
from matplotlib.backends.backend_pdf import PdfPages
from collections import defaultdict, Counter
import matplotlib.patches as mpatches, bokeh.palettes as bp
from bokeh.plotting import figure, output_file, show, ColumnDataSource
from bokeh.layouts import gridplot
from math import pi
from bokeh.models import NumeralTickFormatter, HoverTool, GlyphRenderer, Range1d, LinearColorMapper, BasicTicker, PrintfTickFormatter, ColorBar, ColumnDataSource
from sklearn.cluster import KMeans

##-FUNCTIONS-##
def pd_df_heatmap_variance(data):
    '''From highmutated data of vcf files makes
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

def pd_df_heatmap_sequence(data_dict):
    '''From dictionary makes panda's dataframe.
    Firstly makes dictionary with bases to 0
    values and inputs data_dict into that dictionary
    while counting it.
    Next a df is build from this dictionary with index
    is keys from data_dict
    
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
    
    return df

def vcf_heatmap_snps(data_surrounding, data_variance):
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
        if se in dict_variance:
            dict_variance[se].append(mutation[points])
        else:
            dict_variance[se] = [mutation[points]]
        points += 1
        
    return dict_variance

def vcf_all_strip(path, txt_yes_no, txt_name, lenght):
    '''strips vcf files of it's information
    
    Args:
        path where vcf files are for analyzing.
        if you want a txt file of sequences that have
        missing reads.
        what the txt name should be.
        
    Return:
        backbone_data = list of mutated sequence with
        25% or more reads with mutated base and 3 sequences
        before and after = [0]
        variance_data = same for variance data = [1]
        highmutated_b = only the mutated sequence with 25%
        or more reads with the mutated base without other
        sequences = [2]
        highmutated_v = same for variance data = [3]
        
    Notes:
        some sequences (only in backbone) have reads that should
        be there but can't be set on not or mutated. for example
        score is 4, 1:6 here 4 normal, 1 mutated but should be 6 
        reads in total. These sequences can be formatted in a txt
        file if set to true, with id_
    
    '''
    os.chdir(path)
    backbone_data = []
    variance_data = []
    highmutated_v = []
    highmutated_b = []

    for file in glob.iglob('*.vcf'):
        data = data_vcf_file(file)
        data_all = vcf_whole_sequence_strip(file)
        if bool(data) == True:
            id_ = list(data.keys())[0]
            variance = data[id_]['variance']
            backbone = data[id_]['backbone']

            variance_data.append(mutated_reads_vcf_only(variance, data_all, lenght)[0])
            backbone_data.append(mutated_reads_vcf_only(backbone, data_all, lenght)[0]) 
            highmutated_v.append(mutated_reads_vcf_only(variance, data_all, lenght)[1])
            highmutated_b.append(mutated_reads_vcf_only(backbone, data_all, lenght)[1])

            if txt_yes_no == 'yes':
                append_txt_file(txt_name, data, id_)
    return backbone_data, variance_data, highmutated_b, highmutated_v

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

def filter_empty(s):
    '''filter s if empty
    
    '''
    if bool(s) == True:
        return True
    return False

def filter_score(s):
    ''' filter s if not len >= 10
    
    '''
    if len(s) >= 10:
        return True
    return False

def append_txt_file(file_name, data, id_):
    '''append in file_name.txt with input data.
    if file_name doesn't exist, it will make the file
    
    '''
    data_ = np.array(vcf_reads_dissapear(data))

    with open('{}.txt'.format(file_name), 'a') as file:
        file.write('{}'.format(id_))
        file.write('\n')
        file.write('{}'.format(data_))
        file.write('\n')
        file.write('\n')
        file.close()
    
def mutated_reads_vcf_only(variance_or_backbone, data_all, lenght):
    '''mutations in vcf file with 3 sequences before and
    3 sequences after mutation. Overlap can occur if mutations
    are found beside each other. if location in sequence has
    multiple SNPs like A -> G/T. It will check if score for both
    these SNPs is high enough and will return as A -> G. So
    same location sequence can occur twice in return if A -> G
    and A -> T have both high enough scores
    lenght is the amount of bases next to the variance
    
    Args:
        either variance or backbone data can be entered and
        data_all which is the strip of all the sequencing data
        in the same vcf file
    
    '''
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
                        extended.append(data_all[i-lenght:i+(lenght+1)])
                        continue
            if n3/(n1+n2+n3) > 0.25:
                for i, items in enumerate(data_all):
                    mutated = ':'+','.join(item[0:3])+':'+score3[points]
                    if mutated in items:
                        r_ = items.split(breakpoint, 5)
                        r = r_[4][0]
                        removed = r_[0]+'\t'+r_[1]+'\t'+r_[2]+'\t'+r_[3]+'\t'+r+'\t'+r_[5]
                        highmutated.append(removed)
                        extended.append(data_all[i-lenght:i+(lenght+1)])
        else:
            n1 = int(item[0])
            n2 = int(item[1])
            if n2/(n1+n2) > 0.25:
                for i, items in enumerate(data_all):
                    mutated = ':'+','.join(item[0:2])+':'+score3[points]
                    if mutated in items:
                        highmutated.append(items)
                        extended.append(data_all[i-lenght:i+(lenght+1)])
        points += 1
    points = False
    return extended, highmutated

    
def vcf_reads_dissapear(data):
    '''In vcf files, some reads dissapear
    this function picks and strips those sequences 
    from the vcf file.
    
    Args:
        data = from data_vcf_file/dictionary of {id_ :
        {variance: data, backbone: data}}
    
    '''
    not_enough_reads_ = []
    not_enough_reads = []
    n_e_r2 = []
    n_e_r3 = []
    n_e_r4 = []
    n_e_r5 = []
    for v in data.values():
        for b in v.values():
            for items in b:
                breakpoint = '\t'
                not_enough_reads_.append(items.split(breakpoint, 9)[9])
        
    for i in filter(filter_score, not_enough_reads_):
        not_enough_reads.append(i)

    for s in not_enough_reads:
        n_e_r2.append(s.split(':')[1].split(','))

    for s2 in not_enough_reads:
        n_e_r3.append(s2.split(':')[2]+':'+s2.split(':')[3]+':'+s2.split(':')[4])

    to_few_reads_raw = []
    points = 0   
    for i in n_e_r2:
        number1 = int(i[0])
        number2 = int(i[1])
        if (number1+number2) == int(n_e_r3[points]):
            continue 
        elif (number1+number2 ) != (n_e_r3[points]):
            for items in b:
                join_ = ','.join(i[0:2])
                join_2 = ':'+join_+':'+n_e_r3[points]
                if join_2 in items:
                    to_few_reads_raw.append(items)
        points += 1
    points = False
    to_few_reads = list(set(to_few_reads_raw))
    return to_few_reads

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

def data_vcf_file(file_name):
    '''strips data from vcf files
    
    returns:
        data = {id_ of file: {'variance': variance_data,
        'backbone': backbone_data}
    
    '''
    data = {}
    with open(file_name, 'r') as f:
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
    return data

def dict_values_to_other_dict(dict_x, dict_y):
    '''values of dict_y to keys of dict_x
    if dict_x v for k is None change it to v
    of dict_y otherwise append v of dict_y
    
    Example:
        dict_x = {A: None, B: 0.1)
        dict_y = {A: 2, B: 0.3)
        dict_x = {A: 2, B: [0.1, 0.3]
    
    '''
    for k, v in dict_y.items():
        if dict_x[k] is None:
            dict_x[k] = [v]
        else:
            dict_x[k].append(v)
            
            
def list_with_all_combinations(letters, size):
    '''list of all combination letters
    
    '''
    all_combinations = []
    for i in itertools.product('ACTG', repeat = size):
        size_ = '{}'*size
        all_combinations.append(size_.format(*i))
    return all_combinations
        
def zero_to_nan(values):
    """Replace every 0 with 'nan' and return a copy."""
    return [float('nan') if x==0 else x for x in values]

def get_list_hml_scores_fastq(path, size, overlap, letters):
    '''strip done_fastq files from path -> get meanscore of sequence
    with size and puts into high, medium or low
    
    Args:
        get's path with done_fastq files.
        size is the lenght of the sequences which will be chunked, 
        for instance 4 -> ['ACTG']
        overlap gives how many repetative nucleotides occur, for
        instance 1 -> ['ACTG']['GAAA']['ACTA']
        letters gives what you want to select for 'ACTG' most
        likely
    
    Returns:
        x = bases 
        y_h = amount of high scores
        y_m = amount of medium scores
        y_l = amount of low scores
        
    Notes:
        sequences categorized into high, medium or low
        example:
        ['AAAA']['ACCC']['AAAA']['ACCC'] with score [1][2][3][2]
        y_h = [1, 0]
        y_m = [0, 2]
        y_l = [1, 0]
    
    '''
    os.chdir(path)
    all_y_h = [] 
    all_y_m = []
    all_y_l = []
    
    #gets high, medium and low score in list
    for file in glob.iglob('*.done_fastq'):
        X = get_list_fastq(file, size, overlap, letters)[0]
        y_h = (get_list_fastq(file, size, overlap, letters)[1])
        y_m = (get_list_fastq(file, size, overlap, letters)[2])
        y_l = (get_list_fastq(file, size, overlap, letters)[3])
        all_y_h.append(y_h)
        all_y_m.append(y_m)
        all_y_l.append(y_l)
       
    #does sums of list -> [[0,1,2,0][1,1,1,1]]=[1,2,3,1]
    nc_all_y_h = [[i[n] for i in all_y_h] for n in range(len(all_y_h[0]))]
    c_all_y_h = list(sum(x) for x in nc_all_y_h)
    nc_all_y_m = [[i[n] for i in all_y_m] for n in range(len(all_y_m[0]))]
    c_all_y_m = list(sum(x) for x in nc_all_y_m)
    nc_all_y_l = [[i[n] for i in all_y_l] for n in range(len(all_y_l[0]))]
    c_all_y_l = list(sum(x) for x in nc_all_y_l)
    
    #get dict with many high error rates and thus low quality
    all_combinations = list_with_all_combinations(letters, size)

    all_high_scores = dict(zip(all_combinations, c_all_y_h))
    d = dict((k, v) for k, v in all_high_scores.items() if v >= 60)
    
    #changes 0 to none
    y_nan1 = zero_to_nan(c_all_y_h)
    y_nan2 = zero_to_nan(c_all_y_m)
    y_nan3 = zero_to_nan(c_all_y_l)
    
    return X, y_nan1, y_nan2, y_nan3, d

def get_dict_mean_of_seq_fastq_additional_directories(path, size, overlap, letters):
    '''Same as get_dict_mean_of_seq_fastq with other path step.
    
    Usefull if want to run over multiple directories.
    
    '''
    dict_x = dict.fromkeys(list_with_all_combinations(letters, size))
    
    for subdir, dirs, files in os.walk(path):
        for filename in files:
            if filename.split('.')[-1] == 'done_fastq':
                file = os.path.join(path, subdir, filename)
                data = parse_fasta_file_error(file)
                id_ = list(data.keys())[0]
                error_rate = (convert_qualityscore(data[id_]['score']))
                nucleotide_string = ((data[id_]['sequence']))
                n_split = split_overlap(nucleotide_string, size, overlap)
                listed_nucleotides = []
                for l in n_split:
                    if len(l) >= size:
                        listed_nucleotides.append(l)
                    else:
                        #print('string without {} bases were skipped'.format(size))
                        continue
                error_rate_split = mean_calc_list(list(split_overlap(error_rate, size, overlap)), size)
                dict_z = dict(zip(listed_nucleotides, error_rate_split))
                dict_values_to_other_dict(dict_x, dict_z)

    for k, v in dict_x.items():
        if v == None:
            continue
        else:
            dict_x[k] = statistics.mean(v)

    return dict_x

def get_dict_mean_of_seq_fastq(path, size, overlap, letters):
    '''Dictionary with meanscore of quality score for each sequence
    
    Args:
        get's path with done_fastq files.
        size is the lenght of the sequences which will be chunked, 
        for instance 4 -> ['ACTG']
        overlap gives how many repetative nucleotides occur, for
        instance 1 -> ['ACTG']['GAAA']['ACTA']
        letters gives what you want to select for 'ACTG' most
        likely
        
    Returns:
        {'AAAA':1,'AAAC:1,...,'GGGG':2}
    
    '''
    os.chdir(path)
    dict_x = dict.fromkeys(list_with_all_combinations(letters, size))
    
    for file in glob.iglob('*.done_fastq'):
        data = parse_fasta_file_error(file)
        id_ = list(data.keys())[0]
        error_rate = (convert_qualityscore(data[id_]['score']))
        nucleotide_string = ((data[id_]['sequence']))
        n_split = split_overlap(nucleotide_string, size, overlap)
        listed_nucleotides = []
        for l in n_split:
            if len(l) >= size:
                listed_nucleotides.append(l)
            else:
                #print('string without {} bases were skipped'.format(size))
                continue
        error_rate_split = mean_calc_list(list(split_overlap(error_rate, size, overlap)), size)
        dict_z = dict(zip(listed_nucleotides, error_rate_split))
        dict_values_to_other_dict(dict_x, dict_z)

    for k, v in dict_x.items():
        if v == None:
            continue
        else:
            dict_x[k] = statistics.mean(v)

    return dict_x
    
def get_list_fastq(file, size, overlap, letters):
    '''get high, medium and low score array from fastq file
    
    Args:
        file you want to check for scores per base lenght
        
    Returns:
        X = array with all possible base combinations
        Y = array with same score from those bases for high,
        medium and low quality score
    
    '''
    data = parse_fasta_file_error(file)
    id_ = list(data.keys())[0]
    error_rate = (convert_qualityscore(data[id_]['score']))
    nucleotide_string = ((data[id_]['sequence']))
    scoremean_dict = dict_scoremean_bases(nucleotide_string, error_rate, size, overlap)
    plot_data_h = dict_of_allbases_vs_High_Medium_Low(size, 'High', scoremean_dict, letters)
    plot_data_m = dict_of_allbases_vs_High_Medium_Low(size, 'Medium', scoremean_dict, letters)
    plot_data_l = dict_of_allbases_vs_High_Medium_Low(size, 'Low', scoremean_dict, letters)
    X = list(plot_data_h.keys())
    Y_h_none = list(plot_data_h.values())
    Y_h = [0 if v is None else v for v in Y_h_none]
    Y_m_none = list(plot_data_m.values())
    Y_m = [0 if v is None else v for v in Y_m_none]
    Y_l_none = list(plot_data_l.values())
    Y_l = [0 if v is None else v for v in Y_l_none]
    return X, Y_h, Y_m, Y_l #x, y

def dict_of_allbases_vs_High_Medium_Low(size, qualityscore, data, letters):
    '''Makes dictonairy with all possible base combination
    and number of these base combinations found in either the
    High, Medium or Low quality score group. how long the combination
    is can be decided aswell. 
    
    Args:
        size = size of bases
        qualityscore = High, Medium, Low
    
    Return:
        dict with key = base combination and value = 
        amount of base combination in determined qualityscore
        
    Note:
        This also produces a list of all possible base combinations
        
    '''
    try:
        global qualityscore_
        qualityscore_ = data[qualityscore]
    except KeyError:
        print('')
          
    c_qualityscore = Counter(qualityscore_)
    dict_x = dict.fromkeys(list_with_all_combinations(letters, size))
    z = {**dict_x, **c_qualityscore}
    return z
    
def dict_scoremean_bases(nucleotide_string, error_rate, size, overlap):
    '''Makes dictonairy with meanscore against lenght 
    bases.
    
    Args:
        Nucleotide_string = is the string of all the
        nucleotides
        Error_rate = list of all the error scores
        
    Return:
        dict with key = score and value = 'amount of bases'
        
    Note:
        Error_rate and nucleotide string should be 
        same len
    
    '''
    listed_nucleotides_4 = list(split_overlap(nucleotide_string, size, overlap))
    listed_nucleotides = []
    for l in listed_nucleotides_4:
        if len(l) >= size:
            listed_nucleotides.append(l)
        else:
            #print('string without {} bases were skipped'.format(size))
            continue
            
    statistic_mean = mean_calc_list(list(split_overlap(error_rate, size, overlap)), size)
    listed_scores = high_medium_low_scores(statistic_mean, size)
    results_dict = default_dict_bases_scores(listed_scores, listed_nucleotides)
    return results_dict


def default_dict_bases_scores(keys_list, b_list):
    '''Puts list of keys and bases into dict
    without losing any information because of
    non unique keys
    
    Example:
        Keys_list = ['High','High','Low',
        'Medium']
        b_list = ['AGTC','AGTT','CCTC',
        'AAAA']
        
        {'High' : ['AGTC','AGTT'], Medium :
        ['AAAA'], 'Low' : ['CCTC']
    
    '''
    temp = defaultdict(list)
    for keys, b in zip(keys_list, b_list):
        temp[keys].append(b)
    done = dict(temp)
    return done

def high_medium_low_scores(listed_scores, size):
    '''sets list of scores into high, medium
    or low
    size matters because how bigger the size, how
    closer the scores are to each other and lower
    overall (more in medium)
    
    Args:
        listed_scores = list of scores either means
        or individual scores
        size = lenght of chunked sequence
        
    Return: sets scores to high, medium , low
    
    Note:
    high >= 0.40
    0.20 < medium < 0.40
    low <= 0.20
    
    '''
    group_score = []
    for s in listed_scores:
        if s >= (0.40-0.02*size):
            group_score.append('High')
        elif s <= (0.15+0.01*size):
            group_score.append('Low')
        else:
            group_score.append('Medium')
    return group_score
    
def mean_calc_list(error_rate, lenght):
    '''Calculates the mean of multiple lists
    of lenght scores
    
    Args:
        error_rate = list of four scores in list
        of more scores
        
    Returns:
        meanscores of these list in list into a 
        list
        
    Example:
        [[1,2,3,4][5,6,7,8][9,10,11,12][13,14]]=
        [2.5, 6.5, 10.5]
        
        last two are skipped because list >4
    
    '''
    mean_of_all = []
    points = 0 
    for list in error_rate:
        if len(list) >= lenght:
            mean_of_all.append(statistics.mean(error_rate[points]))
            points += 1
        else:
            continue
            #print('string without {} scores were skipped'.format(lenght))
    return mean_of_all

def split_overlap(iterable,size,overlap):
    '''(list,int,int) => [[...],[...],...]
    Split an iterable into chunks of a specific size and overlap.
    Works also on strings! 
    Examples:
        split_overlap(iterable=list(range(10)),size=3,overlap=2)
        >>> [[0, 1, 2, 3], [2, 3, 4, 5], [4, 5, 6, 7], [6, 7, 8, 9]]
        split_overlap(iterable=range(10),size=3,overlap=2)
        >>> [range(0, 3), range(1, 4), range(2, 5), range(3, 6), range(4, 7), range(5, 8), range(6, 9), range(7, 10)]
    '''
    if size < 1 or overlap < 0:
        raise ValueError('"size" must be an integer with >= 1 while "overlap" must be >= 0')
    result = []
    while True:
        if len(iterable) <= size:
            result.append(iterable)
            return result
        else:
            result.append(iterable[:size])
            iterable = iterable[size-overlap:]
            
def chunks(l, n):
    '''Yield successive n-sized chunks from l.
    
    Args:
        l = list which you want to chunk
        n = how larger chunk become
        
    Example:
        [ACGTACGTACGTACGT] ->
        [[ACGT][ACGT][ACGT][ACGT]]
        
    Note: split_overlap is better version
        
    '''
    for i in range(0, len(l), n):
        yield l[i:i + n]

def plot_done_fastq_files(sequence_file):
    '''plots done_fastq files
    
    Args:
        sequence_file = done_fastq files
        
    Return:
        plot with error rate against individual nucleotides
        
    Note:
        combination of multiple functions.
        get's data with key = id_ and value = sequence and bases.
        this is done using the parse_fasta_file_error function.
        Next these items are stripped into lists (id_, error_rate
        and nucleotides_list).
        X and Y are arrays of these last two lists are made for
        individual nucleotide vs error rate distribution.
        
        Plot is made as function plot_error_rate, but writen
        out in full because otherwise error message that ax can
        not be defined.
        
    '''
    data = parse_fasta_file_error(sequence_file)
                                                             
    id_ = list(data.keys())[0]                               
    error_rate = convert_qualityscore (data[id_]['score']) 

    nucleotides_list = [c for c in (data[id_]['sequence'])]  
    X = np.array(nucleotides_list)                           
    Y = np.array(error_rate)                                 
    
    plt.figure(dpi=150)
    ax=plt.subplot()
    ax.scatter(range(len(Y)),Y,color='b',alpha=0.2, s=3)
    ax.set_ylim(0, 1)
    ax.set_xlim(0, len(X)-1)
    ax.set_xticks(range(len(Y)))
    ax.set_xticklabels(X)
    ax.set_title('error rate distribution')
    ax.set_ylabel('error rate')
    ax.set_xlabel('nucleotide')
    ax.text(len(X)-150, 0.90, statistics.mean(Y))
    
def parse_fasta_file_stripscores(path):
    '''Stripscores from every done_fastq in path
    
    Args:
        path = directory location for done_fastq files
        
    Returns:
        data = list base quality scores for every done_fastq
        file in directory
        
    Example:
        [%$^#*, *#&@^] in which two stripped scores from
        done_fastq files are presented.
    
    '''
    os.chdir(path)
    data = []
    for file in glob.iglob('*.done_fastq'):
        #print(file)
        with open(file, 'r') as f:
            id_ = False
            for line in f:
                #print sequence line
                if line.strip() in ['+','\n']:
                    continue
                if not line[0] in '@ATCG':
                    score_ = line.strip()
                    data.append(score_)
                    score_ = False
    return data

def convert_qualityscore(raw_score):
    '''covert symboles in values from quality score
    
    Args:
        raw_score = quality score in symboles
        
    Return:
        quality score in values

    Example:
        [@&#*] = [3,4,5,6]
    
    '''
    error_rate = []
    for symbol in raw_score:
        try:
            error_rate.append(errordict[symbol])
        except KeyError:
            print("not in errordict")
    return error_rate

def parse_fasta_file(file_name):
    '''strip fasta files from id_ and bases
    
    Args:
        file_name = most types of files such as txt, fasta,
        fastq
    
    Return:
        Dictonairy with key = id_ and value = bases 
    
    Note:
        Doesn't strip any quality scores, just bases.
        Also skips any not base letters (anything not ACGTN)
        
    '''
    data = {}
    with open(file_name, 'r') as f:
        for line in f:
            #print(line)
            skip = False
            if not line.startswith('\n'):
                if line.startswith('>'):
                    id_ = line.strip()[1:]
                else:
                    #skip sequences with unexpected chars
                    for letter in line.strip().upper():
                        if letter not in 'ACGTN':
                            print(f'WARNING: The sequence ">{id_}" contained unexpected chars and it was skipped')
                            skip = True
                            break
                    if not skip:
                        try:
                            data[id_] += line.strip()
                        except KeyError:
                            data[id_] = line.strip()
    return data

def parse_fasta_file_error(sequence_file):
    '''strips sequence file of id, bases and scores
    
    Args:
        sequence_file = most types of files such as txt,
        fasta, fastq
        
    Return:
        dictonairy with key = id_ and value another
        dictonairy with key1 = 'score', key2 = 
        'sequence' value1 = actual score and
        value2 = actual bases
        
    Example:
        {'sequence_run_1' : {'score':'#$%*',
        'sequence':'AGCT'}}
    
    '''
    data = {}
    with open(sequence_file, 'r') as f:
        id_ = False
        sequence = False
        for line in f:
            #print(line)
            if line.strip() in ['+','\n']:
                continue
            elif line.startswith('@'):
                #print('id found')
                if not id_: #if id_ == False
                    id_ = line.strip()
                else:
                    print('WARNING: An ID was found without corresponding sequence.', id_)
                    id_ = line.strip()
            elif id_ and not sequence:
                sequence = line.strip()
            elif id_ and sequence:
                    score = line.strip()
                    data[id_] = {'sequence':sequence,
                                   'score':score}
                    sequence = False 
                    id_ = False
    return data

##-FIGURES AND PLOTS-##
def heatmap_vcf_files_snps(df_):
    '''Heatmap of single nucleotide polymorphisms 
    from vcf files
    
    '''
    df = pd.DataFrame(df_.stack(), columns=['scores']).reset_index()

    mutations = list(df_.columns)
    index = list(df_.index)

    #colors = ["#75968f", "#a5bab7", "#c9d9d3", "#e2e2e2", "#dfccce", "#ddb7b1", "#cc7878", "#933b41", "#550b1d"]
    colors = bp.all_palettes['Viridis'][11]
    mapper = LinearColorMapper(palette=colors, low=0, high=df['scores'].sum()/4)

    source = ColumnDataSource(df)

    TOOLS = "hover,reset"

    p = figure(title='Variant occurence in sequence', y_range=list(reversed(mutations)), plot_width=300, 
               x_axis_location='above', plot_height=600, tools=TOOLS, toolbar_location='below')
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

    show(p)
    
def heatmap_vcf_files_snps_with_sequence(df_):
    '''Heatmap of single nucleotide polymorphisms 
    from vcf files with surrounding sequence (4 bases
    extra)
    
    '''
    df = pd.DataFrame(df_.stack(), columns=['scores']).reset_index()

    bases = list(df_.columns)
    sequences = list(df_.index)

    colors = ["#75968f", "#a5bab7", "#c9d9d3", "#e2e2e2", "#dfccce", "#ddb7b1", "#cc7878", "#933b41", "#550b1d"]
    mapper = LinearColorMapper(palette=colors, low=df.scores.min(), high=df.scores.max())

    source = ColumnDataSource(df)

    TOOLS = "hover,reset,xpan,xwheel_zoom"

    p = figure(title='Variant occurence in sequence', x_range=sequences,
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

    show(p)
    
def plot_vcf_snps(data):
    '''plot single nucleotide polymorphisms
    data should be either mutated sequences from
    backbone or variance
    
    '''
    x = list(dict(highmutated_back_variance(data)).keys())
    y_ = list(dict(highmutated_back_variance(data)).values())
    y = []
    for item in y_:
        y.append(item/sum(y_))

    p = figure(x_range=x, plot_width=400, plot_height=400, title='Distribution of mutations in sequence', background_fill_color="#E8DDCB")
    p.vbar(x, width=0.5, bottom=0, top=y, fill_color="#036564", line_color="#033649")

    p.yaxis[0].formatter = NumeralTickFormatter(format="0.0%")
    p.xaxis.major_label_orientation = pi/4
    p.legend.location = "center_right"
    p.legend.background_fill_color = "darkgrey"
    p.xaxis.axis_label = 'Single Nucleotide Polymorphism'
    p.yaxis.axis_label = 'amount of mutations'

    show(p)
    
def plot_raw_values_bokeh(title, dict_of_x_and_y):
    '''bokeh plot which plots dictionary's of sequence
    and score like {AAAA:1} ect.
    
    '''
    keys = []
    values = []
    for key, value in dict_of_x_and_y.items():
        keys.append(key)
        values.append(value)

    x_plot = np.column_stack((keys, values))

    ylim = len(values)
    figsize = 500
    hover = HoverTool(tooltips=[("Bases", "@x"), ("Qualityscore", "@y")], mode='vline')

    p = figure(x_range=keys, y_axis_label='QaulityScore', tools=[hover, 'xwheel_zoom', 'xpan'])
    p.xaxis.ticker = []
    p.title.text = title
    p.xgrid.grid_line_color = None
    p.ygrid.grid_line_color = None

    p.circle(x = keys, y= values, size=5, alpha=1, color='red', legend='QualityScore')

    p.legend.location = "top_left"
         
    return p

def plot_bases_hml_in_bokeh(title, x,  y_high, y_medium, y_low):
    '''Plot sequences with high, medium and low scores
    on a bokeh plot.
    
    '''
    bases = x
    #percentage
    sum_scores = [x + y + z for x, y, z in zip(y_high, y_medium, y_low)]
    per_h_of_all = [x / y * 100 for x, y in zip(y_high, sum_scores)]
    per_m_of_all = [x / y * 100 for x, y in zip(y_medium, sum_scores)]
    per_l_of_all = [x / y * 100 for x, y in zip(y_low, sum_scores)]

    hover = HoverTool(tooltips=[("Bases", "@x"), ("Scores", "@y""%")], mode='vline')

    #plotting
    p = figure(x_range=bases, y_axis_label='How many found bases in category%', tools=[hover, 'xwheel_zoom', 'xpan'])
    p.title.text = title
    p.xaxis.ticker = []
    p.y_range = Range1d(0, 100)
    p.xgrid.grid_line_color = None
    p.ygrid.grid_line_color = None

    #data
    p.circle(x=bases, y=per_h_of_all, size=5, alpha=1, muted_alpha=0.4, color='red', muted_color='red', legend='High scores')
    p.circle(x=bases, y=per_m_of_all, size=5, alpha=1, muted_alpha=0.4, color='blue', muted_color='blue', legend='Medium scores')
    p.circle(x=bases, y=per_l_of_all, size=5, alpha=1, muted_alpha=0.4, color='purple', muted_color='purple', legend='Low scores')

    p.legend.location = "top_left"
    p.legend.click_policy="mute"

    # output to static HTML file
    #output_file("scatterplot.html", title="High, Medium, Low quality score per sequence")

    show(p)
    return p 

def plot_bases_high_medium_low(x, y1, y2, y3, size, ylim):
    '''plot certain size bases and plot high,
    medium and low quality score from similair basis
    against theese bases
    
    Args:
        x = dictionary of certain size bases
        y1 = high quality score
        y2 = medium quality score
        y3 = low quality score
    
    '''
    y_nan1 = zero_to_nan(y1)
    y_nan2 = zero_to_nan(y2)
    y_nan3 = zero_to_nan(y3)

    red_circle = mpatches.Patch(color='red', label='High Scores')
    blue_circle = mpatches.Patch(color='blue', label='Medium Scores')
    yellow_circle = mpatches.Patch(color='yellow', label='Low Scores')


    plt.figure(dpi=size)
    plt.scatter(x, y_nan1, marker='o', color='r')
    plt.scatter(x, y_nan2, marker='o', color='b')
    plt.scatter(x, y_nan3, marker='o', color='y')
    plt.ylim(0, ylim)
    plt.xlim(0, len(x)-1)
    plt.title('High, Medium, Low quality score per sequence')
    plt.ylabel('How many found bases in category')
    plt.legend(bbox_to_anchor=(1, 1), loc='upper right',
                handles=[red_circle, blue_circle, yellow_circle])

    plt.grid(False) #Grid 
    plt.show()


def plot_read_lenght(x,y):
    '''x vs y => print scatter plot
    
    Args:
        x = lenght of a read
        y = amount of the reads
        
    Returns:
        Scatter plot with amount of reads available against
        the lenght of these reads.
        
    Example:
        2 reads of lenght 180
        4 reads of lenght 160
        ect.
    
    '''
    plt.style.use("ggplot")
    plt.figure(dpi=100)
    plt.scatter(x,y)
    plt.axis([0, 250, 0, 12])
    plt.title('sequence lenght distribution')
    plt.ylabel('amount of reads')
    plt.xlabel('read lenght')
    plt.show()
    
def plot_error_rate(X, Y):
    '''X vs Y => print scatter plot
        
    Args:
        X = array of all the nucleotides 
        Y = array of all the error rates
    
    Returns:
        scatter plot with all individual nucleotides on x-axis
        and the quality score of that nucleotide (from 0-1) on
        the y-axis.
        
    Note:
        adds text with mean of all error rates combined.
        Grid is currently turned off. 
        
    '''
    plt.figure(dpi=150)
    plt.plot((range(len(Y))), Y)
    plt.ylim(0, 1)
    plt.xlim(0, len(X)-1)
    plt.xticks(range(len(Y)), '')
    plt.ylabel('error rate')
    plt.xlabel('nucleotides')
    plt.title('error rate distribution')
    plt.text(len(X)-150, 0.90, statistics.mean(Y))
             
    plt.grid(False) #Grid 
    plt.show()
    
def violin_plot(results):
    '''Results => violin plot
    
    Args:
        Results = either meanscores of error rates from multiple files
        or array of error rates from multiple files 
        
    Return:
        In first scenario, gives one violin plot with the meanscores on
        the y-axis against the number of files with that meanscores on
        the x-axis. 
        In the second scenario it gives multiple violin plot of the
        error rates on the y-axis against the number of bases with 
        that error rate on the x-axis.
        
    Note:
        amount of violin plots for the error rate can be changed (1-100)
    
    '''
    fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(6, 6))
    
    axes.violinplot(results, points=20, widths=0.3,
                      showmeans=True, showextrema=True, showmedians=True)
    axes.set_title('mean distribution of quality scores')
    
    plt.show()

##-Dictonaries-##
    '''Dictonairy which can be used for changing
    Quality score symboles to values    
    '''
errordict = {'!': 1.0, '"': 0.79433, '#': 0.63096, '$': 0.50119, 
            '%': 0.39811, '&': 0.31623, "'": 0.25119, ')': 0.15849, 
             '(': 0.19953, '*': 0.12589, '+': 0.1, ',': 0.07943, 
             '-': 0.06310, '.': 0.05012, '/': 0.03981, '0': 0.03162, 
             '1': 0.02512, '2': 0.01995, '3': 0.01584, '4': 0.01259, 
             '5': 0.01, '6': 0.00794, '7': 0.00631, '8': 0.00501, 
             '9': 0.00398, ':':0.00316, ';': 0.00251, '<': 0.002, 
             '=': 0.00158, '>': 0.00126, '?': 0.001, '@': 0.00079, 
             'A': 0.00063, 'B': 0.00050, 'C': 0.00040, 'D': 0.00032, 
             'E': 0.00025, 'F': 0.00020, 'G': 0.00016, 'H': 0.00013, 
             'I': 0.00010, 'J': 0.00008, 'K': 0.00006}


