#Important, need csv files from fastq_variant_analyser

import argparse, pandas as pd, time, os, numpy as np
from bokeh.plotting import figure, show, output_file, save

def boxplot(df_, output_name, save_path):
    '''boxplot of meanscore from fastq files.
    groups are devided between nanopore runs
    
    '''
    df_.columns.name = 'group'
    cats = list(df_.columns.values)
    df = pd.DataFrame(df_.stack(), columns=['score']).reset_index()
    df = df.drop('Name', axis=1)
    
    # find the quartiles and IQR for each category
    groups = df.groupby('group')
    q1 = groups.quantile(q=0.25)
    q2 = groups.quantile(q=0.5)
    q3 = groups.quantile(q=0.75)
    iqr = q3 - q1
    upper = q3 + 1.5*iqr
    lower = q1 - 1.5*iqr

    # find the outliers for each category
    def outliers(group):
        cat = group.name
        return group[(group.score > upper.loc[cat]['score']) | (group.score < lower.loc[cat]['score'])]['score']
    out = groups.apply(outliers).dropna()

    # prepare outlier data for plotting, we need coordinates for every outlier.
    if not out.empty:
        outx = []
        outy = []
        for cat in cats:
            # only add outliers if they exist
            if not out.loc[cat].empty:
                for value in out[cat]:
                    outx.append(cat)
                    outy.append(value)

    p = figure(tools="save", background_fill_color="#EFE8E2", title="", x_range=cats)

    # if no outliers, shrink lengths of stems to be no longer than the minimums or maximums
    qmin = groups.quantile(q=0.00)
    qmax = groups.quantile(q=1.00)
    upper.score = [min([x,y]) for (x,y) in zip(list(qmax.loc[:,'score']),upper.score)]
    lower.score = [max([x,y]) for (x,y) in zip(list(qmin.loc[:,'score']),lower.score)]

    # stems
    p.segment(cats, upper.score, cats, q3.score, line_color="black")
    p.segment(cats, lower.score, cats, q1.score, line_color="black")

    # boxes
    p.vbar(cats, 0.7, q2.score, q3.score, fill_color="#E08E79", line_color="black")
    p.vbar(cats, 0.7, q1.score, q2.score, fill_color="#3B8686", line_color="black")

    # whiskers (almost-0 height rects simpler than segments)
    p.rect(cats, lower.score, 0.2, 0.01, line_color="black")
    p.rect(cats, upper.score, 0.2, 0.01, line_color="black")

    # outliers
    if not out.empty:
        p.circle(outx, outy, size=6, color="#F38630", fill_alpha=0.6)

    p.xgrid.grid_line_color = None
    p.ygrid.grid_line_color = "white"
    p.grid.grid_line_width = 2
    p.xaxis.major_label_text_font_size="12pt"

    output_file('{}\{}_boxplot.html'.format(save_path, output_name))
    save(p)
    
def main(input_folder, output_name, save_path):
    '''get df out of csv file and make it into
    boxplot, for every file in directory'''
    points = 0
    df = pd.DataFrame({'A' : []})
    points_list = np.arange(500, 500000000, 500)
    
    print('getting inputs')
    for subdir, dirs, files in os.walk(input_folder):
        print('{} files are being processed'.format(len(files)))
        for filename in files:
            points += 1
            if points in points_list:
                print('file {} of {}'.format(points, len(files)))
            if filename.split('.')[-1] == 'csv':
                file = os.path.join(input_folder, subdir, filename)
                df_ = pd.read_csv(file)
                df_ = df_.set_index('Name')
                if df.empty:
                    df = df_
                else:
                    df = pd.merge(df, df_, left_index=True, right_index=True, how='outer')
                
    #boxplot
    boxplot(df, output_name, save_path)
    
if __name__ == '__main__':

    #Argument script
    parser = argparse.ArgumentParser(description="meanscore analyser")
    parser.add_argument("path", help="path of files")
    parser.add_argument("filename", help="output name")
    parser.add_argument("save_path", help="path to save file")
    args = parser.parse_args()
    
    input_folder = args.path
    output_name = args.filename
    save_path = args.save_path
    print(input_folder)
    print(output_name)
    print(save_path)

    print('started')
    start = time.time()
    main(input_folder, output_name, save_path)
    end = time.time()
    print('completed in {}'.format(end-start))