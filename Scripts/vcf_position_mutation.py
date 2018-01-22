import argparse, time, os, pandas as pd, numpy as np, math
from openpyxl import load_workbook

def make_df_percentages(df, size):
    '''Makes the data in df into percentages
    to visualize the occurence of a mutation in a
    specific position better
    '''
    df_index = list(df.index.values) 
    df_columns_name = list(df.columns.values)[0]
    df_columns = list(df.loc[:,df_columns_name])
    
    points = 0
    data_dict = {}
    for items in df_index:
        if items[-8:] in data_dict:
            data_dict[items[-8:]] += df_columns[points]
        else:
            data_dict[items[-8:]] = df_columns[points]

        points += 1
        
    all_dict = dict(zip(df_index, df_columns))
    data = {}
    for i, items in enumerate(all_dict):
        for j, item in enumerate(data_dict):
            if item in items:
                data[list(all_dict.keys())[i]] = ((list(all_dict.values())[i]/list(data_dict.values())[j])*100) 
    
    df = pd.DataFrame({df_columns_name:list(data.values())},
                     index=data.keys())
    df.columns.name = '{}'.format(df_columns_name)
    
    return df
    
def main(input_folder, output_name, save_path, size, i_or_b):
    '''inputs df of csv files derived from vcf files to identify
    mutations % for that position.
    Possible to identify if a mutation is more certain or if it occurs
    in more variaty troughout the sequence
    
    '''
    points = 0
    df = pd.DataFrame({'A' : []})
    points_list = np.arange(500, 500000000, 500)
    
    if size == 5:
        not_size = 3
    elif size == 3:
        not_size = 5
    else:
        not_size = 3
        
    if i_or_b == 'insert':
        in_file = 'backbone'
    elif i_or_b == 'backbone':
        in_file = 'insert'
    else:
        in_file = 'backbone'
    
    print('getting inputs')
    for subdir, dirs, files in os.walk(input_folder):
        print('{} files are being processed'.format(len(files)))
        for filename in files:
            points += 1
            if points in points_list:
                print('file {} of {}'.format(points, len(files)))
            if filename.split('.')[-1] == 'csv':
                if '{}'.format(not_size) in filename:
                    continue
                elif '{}'.format(in_file) in filename:
                    continue
                else:
                    file = os.path.join(input_folder, subdir, filename)
                    df_ = pd.read_csv(file)
                    df_.columns = ['Name','{}'.format(filename)]
                    df_ = df_.set_index('Name')
                    df_p = make_df_percentages(df_, size)
                    
                    if df.empty:
                        df = df_p
                    else:
                        df = pd.merge(df, df_p, left_index=True, right_index=True, how='outer')
                    
    #print to excel
    book = load_workbook('{}/{}_location.xlsx'.format(save_path, output_name))
    writer = pd.ExcelWriter('{}/{}_location.xlsx'.format(save_path, output_name), engine='openpyxl') 
    writer.book = book
    writer.sheets = dict((ws.title, ws) for ws in book.worksheets)
    
    df.to_excel(writer,'{}_{}'.format(i_or_b, size))
    writer.save()
    
    #writer = pd.ExcelWriter('{}/{}_location.xlsx'.format(save_path, output_name))
    #df.to_excel(writer,'{}_{}'.format(i_or_b, size))
    #writer.save()
    print('DataFrame saved as {}_variance_backbone_location at {}'.format(output_name, save_path))
    
if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description="Vcf SNP analyzer")
    parser.add_argument("path", help="path of files")
    parser.add_argument("filename", help="output name")
    parser.add_argument("save_path", help="path to save file")
    parser.add_argument("sequence_size", help="select how many bases sequence analysis should be")
    parser.add_argument("insert_or_backbone", help="select if you want to analyse backbone or insert")
    args = parser.parse_args()
    
    input_folder = args.path
    output_name = args.filename
    save_path = args.save_path
    size = float(args.sequence_size)
    i_or_b = args.insert_or_backbone
    print(input_folder)
    print(output_name)
    print(save_path)
    print('{} is selected size'.format(size))

    print('started')
    start = time.time()
    main(input_folder, output_name, save_path, size, i_or_b)
    end = time.time()
    print('completed in {}'.format(end-start))