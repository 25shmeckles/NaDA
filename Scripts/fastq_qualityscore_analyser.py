import sys
print(sys.version)

import argparse, os, numpy as np, statistics, itertools, pandas as pd, bokeh, math, time
from bokeh.plotting import figure, output_file, save, show 
from bokeh.layouts import gridplot
from bokeh.models import HoverTool, ColumnDataSource
from collections import defaultdict, Counter
from sklearn.cluster import KMeans

def parse_fasta_file_error(sequence_file):
    '''strips fastq files from id_, sequence
    and score and return dict = data
    
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

def convert_qualityscore(raw_score):
    '''covert symboles in values from quality score
    
    '''
    error_rate = []
    
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

    for symbol in raw_score:
        try:
            error_rate.append(errordict[symbol])
        except KeyError:
            print("{} not in error dictionary".format(symbol))
    return error_rate

def split_overlap(iterable,size,overlap):
    '''Split an iterable into chunks of a specific size and overlap.
    Works also on strings! 
    
    '''
    if size < 1 or overlap < 0:
        raise ValueError('"size" must be an integer with >= 1 while "overlap" must be >= 0')
    result = []
    
    while True:
        if len(iterable) <= size:
            result.append(iterable)
            return result
        else:
            result.append(iterable[:int(size)])
            iterable = iterable[int(size)-int(overlap):]
    
def mean_calc_list(error_rate, lenght):
    '''Calculates the mean of multiple lists
    of lenght scores
    
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

def high_medium_low_scores(listed_scores, size):
    '''sets list of scores into high, medium
    or low
    size matters because how bigger the size, how
    closer the scores are to each other and lower
    overall (more in medium)
    
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

def default_dict_bases_scores(keys_list, b_list):
    '''Puts list of keys and bases into dict
    without losing any information because of
    non unique keys
    
    '''
    temp = defaultdict(list)
    for keys, b in zip(keys_list, b_list):
        temp[keys].append(b)
    done = dict(temp)
    return done

def list_with_all_combinations(letters, size):
    '''list of all combination letters
    
    '''
    all_combinations = []
    for i in itertools.product('ACTG', repeat = int(size)):
        size_ = '{}'*int(size)
        all_combinations.append(size_.format(*i))
    return all_combinations

def dict_scoremean_bases(nucleotide_string, error_rate, size, overlap):
    '''Makes dictonairy with meanscore against lenght 
    bases.
    
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
    clean_results_dict = dict(zip(listed_nucleotides, statistic_mean))
    listed_scores = high_medium_low_scores(statistic_mean, size)
    results_dict = default_dict_bases_scores(listed_scores, listed_nucleotides)
    return results_dict, clean_results_dict

def dict_of_allbases_vs_normal (size, data, letters):
    ''' makes dictionary with normal score with all possible
    sequences and normalizes to 0 from none'''
    dict_x = dict.fromkeys(list_with_all_combinations(letters, size))
    dict_z = {**dict_x, **data}
    list_z = [0 if x is None else x for x in dict_z.values()]
    z = dict(zip(list_with_all_combinations(letters, size), list_z))
    return z
    
def dict_of_allbases_vs_High_Medium_Low(size, qualityscore, data, letters):
    '''Makes dictonairy with all possible base combination
    and number of these base combinations found in either the
    High, Medium or Low quality score group. how long the combination
    is can be decided aswell. 
    
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

def get_list_fastq(file, size, overlap, letters):
    '''get high, medium and low score array from data
    input from fastq files

    '''
    data = parse_fasta_file_error(file)
    id_ = list(data.keys())[0]
    error_rate = (convert_qualityscore(data[id_]['score']))
    nucleotide_string = (data[id_]['sequence'])
    scoremean_clean_dict = dict_scoremean_bases(nucleotide_string, error_rate, size, overlap)[1]
    Y_n = dict_of_allbases_vs_normal(size, scoremean_clean_dict, letters)
    
    scoremean_dict = dict_scoremean_bases(nucleotide_string, error_rate, size, overlap)[0]
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
    return X, Y_h, Y_m, Y_l, Y_n, error_rate 

def zero_to_nan(values):
    """Replace every 0 with 'nan' and return a copy."""
    return [float('nan') if x==0 else x for x in values]

def make_dataframe(data_high, data_medium, data_low):
    data_high['Scores'] = 'High'
    data_medium['Scores'] = 'Medium'
    data_low['Scores'] = 'Low'
    df = pd.DataFrame.from_dict([data_high, data_medium, data_low])
    
    df = df.set_index('Scores')
    df.columns.name = 'Sequence'
    
    return df

def save_meanscores_boxplot(meanscore, filenames, save_path, output_name):
    index_ = list(range(1, 1+len(meanscore)))
    df = pd.DataFrame(data=meanscore, index=index_)
    df.columns = [output_name]
    df.index.name = 'Name'
    df.to_csv('{}/{}_boxplot.csv'.format(save_path, output_name)) 
    
    print('DataFrame Saved as {}_boxplot in directory {}'.format(output_name, save_path))
    
#Plotting functions
def score_plotting(df_normal, df_changed, save_path, output_name, plot):
        
    TOOLS = 'hover,reset,xpan,xwheel_zoom'
    
    #Plot of normal scoress
    sequences = list(df_normal.columns)
    scores = (df_normal.values.tolist())[0]
    
    source1 = ColumnDataSource(data=dict(Sequence=sequences, Score=scores))
    
    p1 = figure(title='A', x_range=sequences, y_range=[0, 1], tools=TOOLS)
    p1.xaxis.ticker = []
    p1.xgrid.grid_line_color = None
    p1.ygrid.grid_line_color = None
    p1.circle(x='Sequence', y='Score', size=7, alpha=1, color='blue', legend='Error rate', source=source1)
    p1.select_one(HoverTool).tooltips = [
        ('Sequence', '@Sequence'),
        ('Score', '@Score'),
    ]
    p1.xaxis.axis_label = 'Sequences'
    p1.yaxis.axis_label = 'Error Score'
    
    #Plot of changed scores 
    sequences = list(df_changed.columns)
    scores_h = (df_changed.values.tolist())[0]
    scores_m = (df_changed.values.tolist())[1]
    scores_l = (df_changed.values.tolist())[2]
    sum_of_all = [a+b for a,b in zip([0 if math.isnan(x) else x for x in scores_h],
                                     [a+b for a,b in zip([0 if math.isnan(x) else x for x in scores_m], 
                                                         [0 if math.isnan(x) else x for x in scores_l])])]
    sum_of_all = [np.nan if x == 0 else x for x in sum_of_all]
    scores_h = [a/b*100 for a,b in zip(scores_h, sum_of_all)]
    scores_m = [a/b*100 for a,b in zip(scores_m, sum_of_all)]
    scores_l = [a/b*100 for a,b in zip(scores_l, sum_of_all)]

    #standard devation
    points=0
    standard_deviation = []
    this_line = True
    while this_line == True:
        try:
            data = list([scores_h[points], scores_m[points], scores_l[points]])
            standard_deviation.append(statistics.stdev(data))
            points += 1
        except IndexError:
            this_line = False
       
    source2 = ColumnDataSource(data=dict(Sequence=sequences, Score_h=scores_h,
                                        Score_m=scores_m, Score_l=scores_l,
                                        Deviation=standard_deviation))
    
    p2 = figure(title='B', x_range=sequences, y_range=[0, 100], tools=TOOLS)
    p2.xaxis.ticker = []
    p2.xgrid.grid_line_color = None
    p2.ygrid.grid_line_color = None
    p2.circle(x='Sequence', y='Score_h', size=7, alpha=1, muted_alpha=0.4,
              color='red', muted_color='red', legend='High scores', source=source2)
    p2.circle(x='Sequence', y='Score_m', size=7, alpha=1, muted_alpha=0.4,
              color='blue', muted_color='blue', legend='Medium scores', source=source2)
    p2.circle(x='Sequence', y='Score_l', size=7, alpha=1, muted_alpha=0.4,
              color='purple', muted_color='purple', legend='Low scores', source=source2)
    p2.select_one(HoverTool).tooltips = [
        ('Sequence', '@Sequence'),
        ('High Score', '@Score_h%'),
        ('Medium Score', '@Score_m%'),
        ('Low Score', '@Score_l%'),
        ('Standard deviation', '@Deviation%')
    ]
    
    p2.xaxis.axis_label = 'Sequences'
    p2.yaxis.axis_label = 'Percentage in Category'
    p2.legend.location = "top_right"
    p2.legend.click_policy= "mute"

    grid = gridplot([p1, p2], ncols=2, plot_width=600, plot_height=500, title='sequence error rate in {}'.format(output_name))
    
    
    if plot == True:
        output_file('{}/{}_score_plotting.html'.format(save_path, output_name))
        save(grid)
    return scores_h, scores_l
            
def cluster_plot(scores_h, scores_l, df_changed , save_path, output_name):
    
    TOOLS = 'hover,reset,xpan,xwheel_zoom'
    
    sequences = list(df_changed.columns)
    x_plot = np.column_stack(([[0 if math.isnan(x) else x for x in scores_l],
                               [0 if math.isnan(x) else x for x in scores_h]]))
    kmeans= KMeans(n_clusters=4)
    kmeans.fit(x_plot)
    y_kmeans = kmeans.predict(x_plot)
    colors = np.array([x for x in ('purple', 'blue', 'red', 'green')])
    source = ColumnDataSource(data=dict(Sequence=sequences, Score_h=scores_h,
                                        Score_l=scores_l, color=colors[y_kmeans].tolist()))

    p = figure(title='Clustering of sequence error rate in {}'.format(output_name),
               x_range=[0,100], y_range=[0,100], tools=TOOLS)
    p.circle('Score_l', 'Score_h', color='color', fill_alpha=0.2, source=source)
    centers = kmeans.cluster_centers_
    p.circle(centers[:, 0], centers[:, 1], color='black', size=10, alpha=0.4);
    p.select_one(HoverTool).tooltips = [
        ('Sequence', '@Sequence'),
        ('High Score', '@Score_h%'),
        ('Low Score', '@Score_l%'),
    ]
    
    output_file('{}/{}_score_clustering.html'.format(save_path, output_name))
    save(p)
    
def main(input_folder, output_name, save_path, size, overlap):
    letters = 'ACTG'
    filenames = []
    all_y_h = [] 
    all_y_m = []
    all_y_l = []
    all_y_n = {}
    meanscore = []
    points = 1
    points_list = np.arange(500, 500000000, 500)
    
    print('getting inputs')
    for subdir, dirs, files in os.walk(input_folder):
        print('{} files are being processed'.format(len(files)))
        for filename in files:
            points += 1
            if points in points_list:
                print('file {} of {}'.format(points, len(files)))
            if filename.split('.')[-1] == 'done_fastq':
                filenames.append(filename)
                file = os.path.join(input_folder, subdir, filename)
                X = get_list_fastq(file, size, overlap, letters)[0]
                y_h = get_list_fastq(file, size, overlap, letters)[1]
                y_m = get_list_fastq(file, size, overlap, letters)[2]
                y_l = get_list_fastq(file, size, overlap, letters)[3]
                all_y_h.append(y_h)
                all_y_m.append(y_m)
                all_y_l.append(y_l)
                
                #normal score dict
                y_n = get_list_fastq(file, size, overlap, letters)[4]
                points = -1
                for k, v in y_n.items():
                    points += 1
                    if k in all_y_n:
                        if v == 0:
                            continue
                        else:
                            all_y_n[k] = (v+list(all_y_n.values())[points])/2
                    else:
                        all_y_n[k] = v
                        
                #meanscore
                error_score = get_list_fastq(file, size, overlap, letters)[5]
                meanscore_ = sum(error_score)/len(error_score)
                meanscore.append(meanscore_)
                 
    #does sums of list -> [[0,1,2,0][1,1,1,1]]=[1,2,3,1]
    nc_all_y_h = [[i[n] for i in all_y_h] for n in range(len(all_y_h[0]))]
    c_all_y_h = list(sum(x) for x in nc_all_y_h)
    nc_all_y_m = [[i[n] for i in all_y_m] for n in range(len(all_y_m[0]))]
    c_all_y_m = list(sum(x) for x in nc_all_y_m)
    nc_all_y_l = [[i[n] for i in all_y_l] for n in range(len(all_y_l[0]))]
    c_all_y_l = list(sum(x) for x in nc_all_y_l)
    
    #get dict of error_score
    all_combinations = list_with_all_combinations(letters, size)

    all_high_scores = dict(zip(all_combinations, c_all_y_h))
    d_h = dict((k, v) for k, v in all_high_scores.items() if v >= 60)
    all_medium_scores = dict(zip(all_combinations, c_all_y_m))
    d_m = dict((k, v) for k, v in all_medium_scores.items() if v >= 60)
    all_low_scores = dict(zip(all_combinations, c_all_y_l))
    d_l = dict((k, v) for k, v in all_low_scores.items() if v >= 60)
    
    df_s  = make_dataframe(all_high_scores, all_medium_scores, all_low_scores)
    
    #save data of meanscore
    save_meanscores_boxplot(meanscore, filenames, save_path, output_name)
    
    #score plot
    print('getting score plot')
    df_n = pd.DataFrame(data=all_y_n, index=['Scores'])
    df_n = df_n.replace(0, np.nan)
    df_s = df_s.replace(0, np.nan)

    score_plotting(df_n, df_s, save_path, output_name, True)
    print('score plot compleet')
    
    #clustering
    print('clustering data')
    s_h = score_plotting(df_n, df_s, save_path, output_name, False)[0]
    s_l = score_plotting(df_n, df_s, save_path, output_name, False)[1]
    cluster_plot(s_h, s_l, df_s, save_path, output_name)
    
    
if __name__ == '__main__':
    
    #Argument script
    parser = argparse.ArgumentParser(description="Fastq qualityscore analyzer")
    parser.add_argument("path", help="path of files")
    parser.add_argument("filename", help="output name")
    parser.add_argument("save_path", help="path to save file")
    parser.add_argument("sequence_size", help="select how many bases sequence analysis should be")
    parser.add_argument("overlap", help="select how many bases overlap between sequence shpuld be")
    args = parser.parse_args()

    input_folder = args.path
    output_name = args.filename
    save_path = args.save_path
    size = float(args.sequence_size)
    overlap = float(args.overlap)
    print(input_folder)
    print(output_name)
    print(save_path)

    print('started')
    start = time.time()
    main(input_folder, output_name, save_path, size, overlap)
    end = time.time()
    print('completed in {} seconds'.format(end-start))