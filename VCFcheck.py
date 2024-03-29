import io
import os
import time
import gzip
import numpy as np
import pandas as pd
#import seaborn as sns
#import matplotlib.pyplot as plt
#import feather
from scipy.stats import chisquare
from scipy.stats import zscore
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

import dash
from dash.dependencies import Input, Output, State
import dash_core_components as dcc
import dash_html_components as html
#import dash_table_experiments as dt
import dash_table
from plotly.offline import plot
import plotly.figure_factory as ff
import plotly.graph_objs as go
#from flask_caching import Cache
#import multiprocessing as mp
import urllib.parse
import dask.dataframe as dd

#######################
## Data Manipulation ##
#######################

def read_header(path):
    
    '''Returns the header of the VCF'''
    
    if path.endswith('.gz'):
        f = gzip.open(path, 'rt')
    else:
        f = open(path, 'r')

    header = [l.rstrip() for l in f if l.startswith('##')]

    f.close()

    contig = {}
    info_fields = []
    genotype_fields = []
    filters_fields = []
    for h in header:
        type_h,desc_h = h.split('=',1)
        if type_h == '##contig':
            ncontig = desc_h.split(',',1)[0].split('=')[1]
            lcontig = desc_h.split(',',1)[1].split('=')[1].split('>')[0]
            contig[ncontig] = lcontig
        elif type_h == '##INFO':
            infofield = desc_h.split(',',1)[0].split('=')[1]
            info_fields.append(infofield)
        elif type_h == '##FORMAT':
            formatfield = desc_h.split(',',1)[0].split('=')[1]
            genotype_fields.append(formatfield)
        elif type_h == '##FILTER':
            filterfield = desc_h.split(',',1)[0].split('=')[1]
            filters_fields.append(filterfield)

    return contig, info_fields, genotype_fields, filters_fields, header

#contig, info_fields, genotype_fields, filters_fields = read_header(vcfFile)

def read_vcf(path):

    '''Returns the VCF without header as DataFrame of pandas'''
    dict_dtypes = {}

    if path.endswith('.gz'):
        f = gzip.open(path, 'rt')
    else:
        f = open(path, 'r')

    fout = open('vcftemp','a') 
    #lines = [l for l in f if not l.startswith('##')]
    for l in f:
        if l.startswith('#C'):
            header = l.rstrip('\n')
        if not l.startswith('##'):
            fout.write(l)
    
    #for c in lines[0].split('\t'):
    for c in header.split('\t'):
        #c = c.rstrip('\n')
        if c == 'POS':
            dict_dtypes[c] = 'int'
        elif c == 'INFO':
            dict_dtypes[c] = 'str'
        else:
            dict_dtypes[c] = 'category'

    f.close()
    fout.close()

    #return pd.read_csv(
    #    io.StringIO(''.join(lines)),
    #    sep='\t',
    #    dtype=dict_dtypes,
    #    low_memory=False
    #    ).rename(columns={'#CHROM': 'CHROM'})
    return dd.read_csv(
        'vcftemp',
        sep='\t',
        header=0,
        dtype=dict_dtypes
        ).compute().rename(columns={'#CHROM': 'CHROM'})


#contig, info_fields, genotype_fields, filters_fields = read_header(vcfFile)
#vcf_table = read_vcf(vcfFile)

def obtain_samples(table):

    '''Returns the list of samples'''

    #samples = pd.Series(table.columns[9:])
    samples = table.columns[9:].tolist()

    return samples

#samples = obtain_samples(vcf_table)

def table_indiv(popfile, samples):

    '''Returns a table of individuals'''

    if popfile:
        indiv_table = pd.read_csv(popfile, sep='\t', names=['Individual','Population'])
    else:
        indiv_table = pd.DataFrame({'Individual':samples, 'Population':np.NaN})
    
    return indiv_table

#indiv_table = table_indiv(pop_indFile,samples)

def missing_ind(table,samples):

    '''Returns one list with the missing by individual'''

    vcf_samples = table.iloc[:,9:9+len(samples)].replace('./.', np.NaN)
    missing_individual = (vcf_samples.isnull().sum()/(vcf_samples.isnull().sum()+vcf_samples.notnull().sum()))*100
    
    del vcf_samples
    
    return missing_individual

#missing_individual = missing_ind(vcf_table,samples)
#indiv_table.index = indiv_table.Individual
#indiv_table = pd.concat([indiv_table, missing_individual], axis=1, sort=False, join='inner')
#indiv_table.rename(columns = {0:"Missing"}, inplace=True)

def missing_snp(table,samples):

    '''Returns one list with the missing by SNP'''

    vcf_samples = pd.DataFrame()
    for s in samples:
        colname = s + '_GT'
        if colname in table.columns:
            vcf_samples[s] = table.loc[:,colname].replace('./.', np.NaN)
        else:
            vcf_samples[s] = table.loc[:,s].replace('./.', np.NaN)

    missing_snp = (vcf_samples.T.isnull().sum()/(vcf_samples.T.isnull().sum()+vcf_samples.T.notnull().sum()))*100

    del vcf_samples
    
    return missing_snp

#vcf_table['Missing'] = missing(vcf_table,samples)

def type_pos(table, info_fields):

    '''Generate columns for INDEL and END'''

    if 'INDEL' in info_fields:
        table['INDEL'] = table.INFO.str.extract('(INDEL)',expand=False).astype('category')
    else:
        table['INDEL'] = pd.np.nan

    if 'END' in info_fields:
        string1 = '(END=.*$)'
        string2 = '(?P<END>[-+]?\d*\.?\d+;?)'
        table['END'] = table.INFO.str.extract(string1, expand=False).str.extract(string2,expand=False).str.replace(";", "").astype(float)
    else:
        table['END'] = pd.np.nan

    table.fillna(value=pd.np.nan, inplace=True)
    #table.fillna(value=pd.np.nan)

def add_info_cols(table, info_fields):

    '''Splits the INFO column into different columns'''

    for i in info_fields:
        string1 = '(' + i + '=.*$)'
        string2 = '(?P<' + i + '>[-+]?\d*\.?\d+;?)'
        if i == 'DP4':
            string2 = '(?P<' + i + '>[-+]?\d*\,?\d*,?\d*\,?\d+;)'
            table[i] = table.INFO.str.extract(string1, expand=False).str.extract(string2,expand=False).str.replace(";", "")
        if (i != 'INDEL') and (i != 'DP4'):
            table[i] = table.INFO.str.extract(string1, expand=False).str.extract(string2,expand=False).str.replace(";", "").astype(float)

    table.fillna(value=pd.np.nan, inplace=True)
    #table.drop('INFO', axis=1, inplace=True)

#add_info_cols(vcf_table, info_fields)

def add_genotype_cols(table, gt_fields, samples):

    '''Splits the FORMAT (Genotype) column into different columns'''

    for s in samples:
        colGT = s + '_GT'
        colPL = s + '_PL'
        colDP = s + '_DP'
        df_genotype = table[s].str.split(':', expand=True)
        if len(df_genotype.columns) == 3:
            df_genotype.fillna(value='-1', inplace=True)
            table[colGT] = df_genotype[0].astype('category')
            #table[colPL] = df_genotype[1]
            table[colDP] = df_genotype[2].astype(int)
            table[colDP].replace(-1, pd.np.nan, inplace=True)
            #table.drop(s, axis=1, inplace=True)
            #table.drop('FORMAT', axis=1, inplace=True)
        #else:
        #    table.drop('FORMAT', axis=1, inplace=True)

    table.fillna(value=pd.np.nan, inplace=True)

#add_genotype_cols(vcf_table, genotype_fields)

def table_indels(table):

    '''Returns a table only with the INDELS of the VCF'''

    indel_table = {}

    if 'INDEL' in table.columns:
        indel_table = table[table.INDEL.notnull()]
    else:
        indel_table = table[table.POS.isnull()]

    return indel_table

#indel_table = table_indels(vcf_table)

def table_rohs(table):

    '''Returns a table only with the Homozigous blocks of the VCF'''

    roh_table = {}

    if 'END' in table.columns:
        roh_table = table[table.END.notnull()]
    else:
        roh_table = table[table.POS.isnull()]

    return roh_table

#roh_table = table_rohs(vcf_table)

def table_snps(table):

    '''Returns a table only with the SNPs of the VCF'''

    snp_table = {}

    if set(['INDEL','END']).issubset(table.columns):
        snp_table = table[(table.END.isnull()) & (table.INDEL.isnull())]
    else:
        snp_table = table

    return snp_table

#snp_table = table_snps(vcf_table)
                      
def biallelic_snps(table):

    '''Returns a table only with the biallelic SNPs of the VCF'''

    snp_table = {}

    if set(['INDEL','END']).issubset(table.columns):
        snp_table = table[(table.END.isnull()) & (table.INDEL.isnull())]
    else:
        snp_table = table

    biallelic_table = snp_table[snp_table.ALT.str.contains(',') == False]

    return biallelic_table

#biallelic_table = biallelic_snps(snp_table)

def multiallelic_snps(table):

    '''Returns one table only with the multiallelic SNPs of the VCF'''

    snp_table = {}

    if set(['INDEL','END']).issubset(table.columns):
        snp_table = table[(table.END.isnull()) & (table.INDEL.isnull())]
    else:
        snp_table = table

    multiallelic_table = snp_table[snp_table.ALT.str.contains(',')]

    return multiallelic_table

#multiallelic_table = multiallelic_snps(snp_table)

def number_roh_positions(table):

    '''Returns the number of positions that are homozigous and equal to the ancestral allele'''

    if len(table) != 0:
        table_roh_positions = pd.numeric(table.END) - pd.numeric(table.POS)
        total_number_ROHs = table_roh_positions.sum()
    else:
        total_number_ROHs = 0

    return total_number_ROHs

#number_positions_roh = number_roh_positions(roh_table)

#########################
# Dashboard Layout / View
#########################

## Set up Dashboard and create layout:
app = dash.Dash()
#app.css.append_css({'external_url': 'https://codepen.io/chriddyp/pen/bWLwgP.css'})
app.css.config.serve_locally = True
app.scripts.config.serve_locally = True

#cache = Cache(app.server, config={
#    'CACHE_TYPE': 'filesystem',
#    'CACHE_DIR': 'cache-directory'
#})

## Layout of the web application:
app.layout = html.Div([
    html.Link(
        rel='stylesheet',
        href='/static/bWLwgP.css'
    ),

    ## Page Header
    html.Div([
        html.H3(
            'VCFcheck: A tool for VCF files diagnostics',
            id='title',
            style={'font-weight': 'bold'}
        )
    ]),

    html.Div([

        ## Upload the VCF file:
        dcc.Upload(
            id='upload-vcf',
            children=html.Div([
                html.A('Select the VCF file')
            ]),
            style={
                'width': '12%',
                'height': '60px',
                'lineHeight': '60px',
                'borderWidth': '1px',
                'borderStyle': 'dashed',
                'borderRadius': '5px',
                'textAlign': 'center',
                'margin': '10px',
                'background-color':'rgb(245,245,245)'
            },
            multiple=False, # Don't allow multiple files to be uploaded
            className='three columns',
        ),

        ## Upload the file with the samples and its populations:
        dcc.Upload(
            id='upload-poplist',
            children=html.Div([
                html.A('Select the Population list')
            ]),
            style={
                'width': '12%',
                'height': '60px',
                'lineHeight': '60px',
                'borderWidth': '1px',
                'borderStyle': 'dashed',#'solid'
                'borderRadius': '5px',
                'textAlign': 'center',
                'margin': '10px',
                'background-color':'rgb(245,245,245)'
            },
            multiple=False, # Don't allow multiple files to be uploaded
            className='three columns',
        ),

        ## Select the type of mutations to show:
        dcc.RadioItems(
            id='radio-typeSNP',
            options=[
                {'label': 'All positions', 'value': 'all'},
                {'label': 'SNPs and INDELs', 'value': 'snpindel'},
                {'label': 'SNPs', 'value': 'snps'},
                {'label': 'Biallelic SNPs', 'value': 'biallelic'},
                {'label': 'Multiallelic SNPs', 'value': 'multiallelic'},
                {'label': 'INDELs', 'value': 'indels'},
                {'label': 'ROHs', 'value': 'rohs'}
            ],
            value='all',
            className='three columns',
        ),

        ## Select the range of Depth by Individual:
        html.Div([
            dcc.RangeSlider(
                id='sliderDP',
                min=0,
                max=0,
                step=1,
                value=[0,0],
                updatemode='drag',
                className='two columns',
                disabled=True
            ),
            html.Div(id='output-container-range-sliderDP', className='three columns')
        ], style={'padding': 20}),

        ## Select the range of Map quality by SNP:
        html.Div([
            dcc.RangeSlider(
                id='sliderMQ',
                min=0,
                max=0,
                step=1,
                value=[0,0],
                updatemode='drag',
                className='two columns',
                disabled=True
            ),
            html.Div(id='output-container-range-sliderMQ', className='three columns')
        ], style={'padding': 20}),

        ## Select the range of Missing by SNP:
        html.Div([
            dcc.Slider(
                id='sliderMiss',
                min=0,
                max=100,
                step=1,
                value=100,
                updatemode='drag',
                className='two columns',
                disabled=True
            ),
            html.Div(id='output-container-range-sliderMiss', className='three columns')
        ], style={'padding': 20}),

    ],id='div1',className='row'),
    
    ## Objects to go between callbacks:
    html.Div(id='output-vcf', style={'display': 'none'}),
    html.Div(id='output-vcf-filt', style={'display': 'none'}),
    html.Div(id='output-samples', style={'display': 'none'}),
    html.Div(id='header', style={'display': 'none'}),
    html.Div(id='nsnpsbi', style={'display': 'none'}),
    html.Div(id='nsnpsmulti', style={'display': 'none'}),
    html.Div(id='nindels', style={'display': 'none'}),
    html.Div(id='nrohs', style={'display': 'none'}),
    html.Div(id='list_avg_miss', style={'display': 'none'}),
    html.Div(id='list_avg_dp_gt', style={'display': 'none'}),
    html.Div(id='list_dict_perc_gt', style={'display': 'none'}),

    ## Warning Pop-ups:
    #dcc.ConfirmDialog(
    #    id='warningVCF',
    #    message='VCF file not uploaded!',
    #),

    ## Name of inputs files:
    html.Div('Input VCF file: -', id='output-name'),
    html.Div('Input Population file: -', id='output-popname'),

    ## Button for start the app:
    html.Div([
        html.Button('Start', id='button-start', n_clicks=0),
    ], style={'padding': 20}),

    ## Horizontal line:
    html.Hr(),

    html.Div([

        ## Datatable of VCF:
        html.Div(id='datatable-upload', className='six columns'),

        ## Select Plot in the Dropdown:
        html.Div([
            html.Div([
                html.Div([
                    html.Div('Select Plot', className='two columns'),
                    html.Div(dcc.Dropdown(id='plot-selector', 
                             ),className='four columns')
                ]),
            ], className='six columns'),
        ], style={'padding': 50}),

        ## Plot:
        html.Div(id='plot-graph', className='six columns'),

        ## Button for download the filtered VCF:
        html.Div([
            html.Button('Download VCF', id='button-download', n_clicks=0),
            html.Div(id='output-download'),
        ], style={'padding': 20}),

        ## Button for download the summary of results:
        html.Div([
            html.Button('Download summary', id='button-summary', n_clicks=0),
            html.Div(id='output-summary'),
        ], style={'padding': 20}),

    ],id='div2',className='row'),

    ## Horizontal line:
    html.Hr(),

    ## Summary of VCF (number of mutations, missing per population and average depth per coverage
    html.Div(id='summary'),

    ## Horizontal line:
    html.Hr(),

])

@app.server.route('/static/<path:path>')
def static_file(path):
    static_folder = os.path.join(os.getcwd(), 'static')
    return send_from_directory(static_folder, path)

#############################################
# Interaction Between Components / Controller
#############################################

def generate_table(df):

    '''Given a dataframe, return a table generated in the web using Dash components'''

    df.replace(pd.np.nan, 'NA', inplace=True)

    return html.Div([
        dash_table.DataTable(
            data=df.to_dict('rows'),
            columns=[{'name': i, 'id': i} for i in df.columns],
            #n_fixed_columns=9,
            #n_fixed_rows=1,
            style_as_list_view=True,
            style_table={
                'overflowX': 'scroll',
                'maxWidth': '800',
                #'overflowY': 'scroll',
                'maxHeight': '600',  
            },
            pagination_settings={
                'current_page': 0,
                'page_size': 50
            },
            #pagination_mode='be',
            #style_cell={
            #    'minWidth': '60px', 'maxWidth': '180px',
            #},
        ),
    ])

def draw_missing_pop(plot,table):

    '''Representation of the distribution of missing data by population'''

    pops = table['Population'].dropna().unique()

    if len(pops) != 0:
        figure = ff.create_distplot([table[table['Population'] == p]['Missing'] for p in pops], pops, show_hist=False)
        return figure

def draw_missing_snp(plot,table):

    '''Representation of the distribution of missing data by SNP'''

    if not table.empty:
        figure = ff.create_distplot([table], ['Missing'], show_hist=False)
        return figure

def draw_raf(plot,table):

    '''Representation of the distribution of the reference allele frequency'''

    if not table.empty:
        figure = ff.create_distplot([table], ['Reference allele frequency'], show_hist=False)
        return figure

def draw_depth_ind(plot,table,samples):

    '''Representation of the distribution of the depth by individual'''

    if (samples[0] + '_DP') in table.columns:    
        figure = ff.create_distplot([table[s + '_DP'].dropna() for s in samples], samples, show_hist=False)
        return figure

def draw_depth_gt(plot,table,samples):

    '''Representation of the distribution of the depth by genotype'''

    if (samples[0] + '_DP') in table.columns:    

        gtposib = ['0/0','0/1','1/1','0/2','1/2','2/2','0/3','1/3','2/3','3/3']
    
        list_gt = {}
    
        for s in samples:
            if s + '_GT' in table.columns:
                for gt in gtposib:
                    gtlistdp = table[table[s + '_GT'] == gt][s + '_DP'].dropna().values.tolist()
                    if len(gtlistdp) != 0:
                        if gt in list_gt:
                            list_gt[gt] += gtlistdp
                        else:
                            list_gt[gt] = gtlistdp
        if list_gt != {}:
            figure = ff.create_distplot([list_gt[key] for key in list_gt if len(list_gt[key]) > 1], [key for key in list_gt.keys() if len(list_gt[key]) > 1], show_hist=False)
            return figure

def draw_depth_vs_gt(plot,df,samples,valuedepth):

    '''Amount of each genotype per sample depth'''
    
    vcf_samples = pd.DataFrame()
    sum00_list = []
    sum01_list = []
    sum11_list = []
    ranges_depth = []

    length_of_range = round(((valuedepth[1] - valuedepth[0] + 1) / 10) + 0.5)
    for i in range(1,11):
        start_of_range = length_of_range * (i-1) + valuedepth[0]
        end_of_range = start_of_range + length_of_range - 1
        if end_of_range > valuedepth[1]:
            end_of_range = valuedepth[1]
        if start_of_range >= valuedepth[1]:
            break
        ranges_depth.append(str(str(start_of_range) + '-' + str(end_of_range)))

        ngt = {}
        sum00 = 0
        sum01 = 0
        sum11 = 0
        for s in samples:

            s_gt = s + '_GT'
            s_dp = s + '_DP'

            if s_gt in df.columns:
                ngt[s] = df.loc[(df[s_dp] >= start_of_range) & (df[s_dp] <= end_of_range)].groupby([s_gt])[s_dp].count()
                if '0/0' in ngt[s]:
                    sum00 += ngt[s]['0/0']
                if '0/1' in ngt[s]:
                    sum01 += ngt[s]['0/1']
                if '1/1' in ngt[s]:
                    sum11 += ngt[s]['1/1']

        sum00_list.append(sum00)
        sum01_list.append(sum01)
        sum11_list.append(sum11)

    sum_total = [sum00_list,sum01_list,sum11_list]

    gt_types = ['0/0','0/1','1/1']

    data = []

    for i in range(0,len(gt_types)):
        trace = go.Bar(
            x=ranges_depth,
            y=sum_total[i],
            name=gt_types[i]
        )
        data.append(trace)

    return data

def draw_pca(plot,table):

    '''Representation of the PCA'''

    pops = table['Population'].unique()
    data = []
    #colors = ['red', 'green', 'blue']
    
    for name in pops:
        indicesToKeep = table['Population'] == name
        data.append(go.Scatter(x = table.loc[indicesToKeep, 'principal component 1'],
                               y = table.loc[indicesToKeep, 'principal component 2'],
                               name=name,
                               mode = 'markers',
                               marker= dict(size = 10,
                                            line = dict(width = 2)
                                            )
                               )
                    )
    return data

def draw_hwe(plot,table,pops):

    '''Representation of the distribution of the Hardy-Weinberg p-Value'''
    if len(pops) != 0:
        figure = ff.create_distplot([table['HWE_pval_' + p].dropna() for p in pops], pops, show_hist=False)
        return figure

def draw_inbreed(plot,table,pops):

    '''Representation of the distribution of the Inbreeding coefficient'''

    if len(pops) != 0:
        figure = ff.create_distplot([table[p].dropna() for p in pops], pops, show_hist=False)
        return figure

## Write the VCF filename:
@app.callback([Output('output-name', 'children'),
              Output('output-popname', 'children')],
              [Input('upload-vcf', 'filename'),
              Input('upload-poplist', 'filename')])
def update_nameoutput(filename, popfile):
    if filename:
        if popfile:
            return 'Input VCF file: ' + filename, 'Input Population file: ' + popfile
        else:
            return 'Input VCF file: ' + filename, 'Input Population file: -'
    else:
        if popfile:
            return 'Input VCF file: -', 'Input Population file: ' + popfile
        else:
            return 'Input VCF file: -', 'Input Population file: -'

## Warning if VCF file is not uploaded:
#@app.callback(Output('warningVCF', 'displayed'),
#              [Input('button-start', 'n_clicks')],
#              [State('upload-vcf', 'filename')])
#def update_nameoutput(nclicks, filename):
#    if nclicks > 0 and not filename:
#        return True
#    return False

## Upload VCF file:
@app.callback([Output('output-samples', 'children'),
              Output('header', 'children'),
              Output('sliderDP', 'min'),
              Output('sliderDP', 'max'),
              Output('sliderDP', 'disabled'),
              Output('sliderDP', 'value'),
              Output('sliderMQ', 'min'),
              Output('sliderMQ', 'max'),
              Output('sliderMQ', 'disabled'),
              Output('sliderMQ', 'value'),
              Output('sliderMiss', 'max'),
              Output('sliderMiss', 'value'),
              Output('sliderMiss', 'disabled'),
              Output('nsnpsbi', 'children'),
              Output('nsnpsmulti', 'children'),
              Output('nindels', 'children'),
              Output('nrohs', 'children'),
              Output('plot-selector','options'),
              Output('list_avg_miss','children'),
              Output('list_dict_perc_gt','children'),
              Output('list_avg_dp_gt','children')],
              [Input('button-start', 'n_clicks')],
              [State('radio-typeSNP', 'value'),
              State('upload-vcf', 'filename'),
              State('upload-poplist', 'filename')])
def update_nameoutput(nclicks, type_snps, filename, popfile):

    if nclicks > 0 and filename:

        starttime = time.time()

        ## Read the input VCF file (mutations and header separately)
        df = read_vcf(filename)
        contig, info_fields, genotype_fields, filters_fields, header = read_header(filename)
        #df_dd = dd.from_pandas(df, npartitions=3)
        endtime = time.time()
        print('Open VCF: ' + str(endtime-starttime))
        print('0:', time.time())
        ## Obtain a list of the VCF samples
        samples = obtain_samples(df)
        if os.path.exists('vcftemp'):
            os.remove('vcftemp')
        #else:
        #    print('File does not exists')

        ## Obtain the possible type of positions from the header (SNPs, INDELs or ROHs)
        type_pos(df, info_fields)

        nSNPsbi = len(df[(df.END.isnull()) & (df.INDEL.isnull()) & (df.ALT.str.contains(',') == False)])
        nSNPsmulti = len(df[(df.END.isnull()) & (df.INDEL.isnull()) & (df.ALT.str.contains(','))])
        nINDELs = len(df[df.INDEL.notnull()])
        nROHs = ((df[df.END.notnull()].END) -(df[df.END.notnull()].POS)).sum()
        
        add_genotype_cols(df, genotype_fields, samples)
        add_info_cols(df, info_fields)

        list_avg_miss = {}
        list_avg_dp_gt = {}
        list_dict_perc_gt = {}

        plots_options = {}

        if (samples[0] + '_DP') in df.columns:
            gtposib = ['0/0','0/1','1/1','0/2','1/2','2/2','0/3','1/3','2/3','3/3']
            list_gt = {}
            for s in samples:
                if s + '_GT' in df.columns:
                    for gt in gtposib:
                        gtlistdp = df[df[s + '_GT'] == gt][s + '_DP'].dropna().values.tolist()
                        if len(gtlistdp) != 0:
                            if gt in list_gt:
                                list_gt[gt] += gtlistdp
                            else:
                                list_gt[gt] = gtlistdp
            for g in list_gt:
                avg_dp_gt = {}
                avg_dp_gt['Genotype'] = g
                avg_dp_gt['Depth'] = round(sum(list_gt[g])/len(list_gt[g]),1)
                for k, v in avg_dp_gt.items():
                    if k in list_avg_dp_gt:
                        list_avg_dp_gt[k].append(v)
                    else:
                        list_avg_dp_gt[k] = [v]


            if popfile is not None:
                plots = ["Missing by Population", 
                         "Missing by SNP", 
                         "Reference allele frequency",
                         "Depth by individual",
                         "Depth by genotype",
                         "Count genotypes per Depth",
                         "Principal Component Analyisis (PCA)",
                         "Hardy-Weinberg test",
                         "Inbreeding coefficient"]
                plots_options = ([{'label': p, 'value': p} for p in plots])

                missing_individual = missing_ind(df, samples)
                indiv_table = table_indiv(popfile,samples)
                indiv_table.index = indiv_table.Individual
                indiv_table = pd.concat([indiv_table, missing_individual], axis=1, sort=False, join='inner')
                indiv_table.rename(columns = {0:"Missing"}, inplace=True)
                tablepop = indiv_table.replace('./.', np.NaN)
                del indiv_table
                pops = tablepop['Population'].unique()
    
                gtposib = ['0/0','0/1','1/1','0/2','1/2','2/2','0/3','1/3','2/3','3/3','./.']

                for p in pops:
                    avg_miss = {}
                    avg_miss['Population'] = p
                    avg_miss['Missing'] = str(round(tablepop[tablepop['Population'] == p]['Missing'].mean(),2)) + '%'
                    for k, v in avg_miss.items():
                        if k in list_avg_miss:
                            list_avg_miss[k].append(v)
                        else:
                            list_avg_miss[k] = [v]

                    dict_perc_gt = {}
                    dict_perc_gt['Population'] = p
                    samples_p = tablepop[tablepop['Population'] == p]['Individual']
                    ngt = {}
                    count_gt = {}
                    total_gt = 0
                    for s in samples_p:
                        ngt[s] = df.groupby([s + '_GT'])[s + '_GT'].count()
                        total_gt += ngt[s].sum()
                        for gt in gtposib:
                            if gt in ngt[s]:
                                if gt in count_gt:
                                    count_gt[gt] += ngt[s][gt]
                                else:
                                    count_gt[gt] = ngt[s][gt]
                    for g in gtposib:
                        if g in count_gt:
                            dict_perc_gt[g] = str(round((count_gt[g]/total_gt)*100,2)) + '%'
                        else:
                            dict_perc_gt[g] = str(0.00) + '%'

                    for k, v in dict_perc_gt.items():
                        if k in list_dict_perc_gt:
                            list_dict_perc_gt[k].append(v)
                        else:
                            list_dict_perc_gt[k] = [v]

            else:
                popfile = '-'
  
                plots = ["Missing by SNP", 
                         "Reference allele frequency",
                         "Depth by individual",
                         "Depth by genotype",
                         "Count genotypes per Depth"]

                plots_options = ([{'label': p, 'value': p} for p in plots])

        else:
            if popfile is not None:
                plots = ["Missing by Population", 
                         "Missing by SNP", 
                         "Reference allele frequency",
                         "Principal Component Analyisis (PCA)",
                         "Hardy-Weinberg test",
                         "Inbreeding coefficient"]
                plots_options = ([{'label': p, 'value': p} for p in plots])

                missing_individual = missing_ind(df, samples)
                indiv_table = table_indiv(popfile,samples)
                indiv_table.index = indiv_table.Individual
                indiv_table = pd.concat([indiv_table, missing_individual], axis=1, sort=False, join='inner')
                indiv_table.rename(columns = {0:"Missing"}, inplace=True)
                tablepop = indiv_table.replace('./.', np.NaN)
                del indiv_table
                pops = tablepop['Population'].unique()
    
                gtposib = ['0/0','0/1','1/1','0/2','1/2','2/2','0/3','1/3','2/3','3/3','./.']

                for p in pops:
                    avg_miss = {}
                    avg_miss['Population'] = p
                    avg_miss['Missing'] = str(round(tablepop[tablepop['Population'] == p]['Missing'].mean(),2)) + '%'
                    for k, v in avg_miss.items():
                        if k in list_avg_miss:
                            list_avg_miss[k].append(v)
                        else:
                            list_avg_miss[k] = [v]

                    dict_perc_gt = {}
                    dict_perc_gt['Population'] = p
                    samples_p = tablepop[tablepop['Population'] == p]['Individual']
                    ngt = {}
                    count_gt = {}
                    total_gt = 0
                    for s in samples_p:
                        ngt[s] = df.groupby([s])[s].count()
                        total_gt += ngt[s].sum()
                        for gt in gtposib:
                            if gt in ngt[s]:
                                if gt in count_gt:
                                    count_gt[gt] += ngt[s][gt]
                                else:
                                    count_gt[gt] = ngt[s][gt]
                    for g in gtposib:
                        if g in count_gt:
                            dict_perc_gt[g] = str(round((count_gt[g]/total_gt)*100,2)) + '%'
                        else:
                            dict_perc_gt[g] = str(0.00) + '%'

                    for k, v in dict_perc_gt.items():
                        if k in list_dict_perc_gt:
                            list_dict_perc_gt[k].append(v)
                        else:
                            list_dict_perc_gt[k] = [v]

            else:
                popfile = '-'
  
                plots = ["Missing by SNP", 
                         "Reference allele frequency"]
                plots_options = ([{'label': p, 'value': p} for p in plots])

        ## Filter by the selected type of position
        df_filt = {}
        if type_snps == 'all':
            df_filt = df
            df_filt.drop('INDEL', axis=1, inplace=True)
        elif type_snps == 'snpindel':
            df_filt = df
            df_filt.drop('INDEL', axis=1, inplace=True)
            df_filt.drop('END', axis=1, inplace=True)
        elif type_snps == 'snps':
            df_filt = table_snps(df)
            df_filt.drop('INDEL', axis=1, inplace=True)
            df_filt.drop('END', axis=1, inplace=True)
        elif type_snps == 'biallelic':
            df_filt = biallelic_snps(df)
            df_filt.drop('INDEL', axis=1, inplace=True)
            df_filt.drop('END', axis=1, inplace=True)
        elif type_snps == 'multiallelic':
            df_filt = multiallelic_snps(df)
            df_filt.drop('INDEL', axis=1, inplace=True)
            df_filt.drop('END', axis=1, inplace=True)
        elif type_snps == 'indels':
            df_filt = table_indels(df)
            df_filt.drop('INDEL', axis=1, inplace=True)
            df_filt.drop('END', axis=1, inplace=True)
        elif type_snps == 'rohs':
            df_filt = table_rohs(df)
            df_filt.drop('INDEL', axis=1, inplace=True)

        del df

        for c in df_filt.columns:
            if df_filt[c].notnull().sum() == 0:
                df_filt.drop(c, axis=1, inplace=True)

        ## Configure the sliders by the VCF values of Sample Depth, Mapping quality and Missing data
        minDP_i = 100
        maxDP_i = 0
        minDP_f = 0
        maxDP_f = 0
        for s in samples:
            colname = s + '_DP'
            if colname in df_filt.columns:
                disableDP = False
                minDP_f = minDP_i
                maxDP_f = maxDP_i
                minDP = df_filt[colname].min()
                maxDP = df_filt[colname].max()
                if minDP < minDP_f:
                    minDP_f = minDP
                if maxDP > maxDP_f:
                    maxDP_f = maxDP
            else:
                disableDP = True
        valueDP = [minDP_f,maxDP_f]

        if 'MQ' in df_filt.columns:
            minMQ = df_filt.MQ.min()
            maxMQ = df_filt.MQ.max()
            disableMQ = False
        else:
            minMQ = 0
            maxMQ = 0
            disableMQ = True
        valueMQ = [minMQ,maxMQ]

        if type_snps == 'all':
            disableDP = True
            disableMQ = True

        df_filt['Missing'] = missing_snp(df_filt, samples)
        maxMiss = df_filt.Missing.max()
        valueMiss = maxMiss
        disableMiss = False

        ## Convert the Pandas DataFrame into Json table to the transport of the table to other callback
        #df_json = df_filt.to_json(date_format='iso', orient='table')
        df_json_list_avg_miss = pd.DataFrame.from_dict(list_avg_miss).to_json(date_format='iso', orient='table')
        df_json_list_dict_perc_gt = pd.DataFrame.from_dict(list_dict_perc_gt).to_json(date_format='iso', orient='table')
        df_json_list_avg_dp_gt = pd.DataFrame.from_dict(list_avg_dp_gt).to_json(date_format='iso', orient='table')
        
        print('1:' + str(time.time()))
        #feather.write_dataframe(df_filt, 'df.feather')
        df_filt.reset_index().to_feather('df.feather')

        return samples, header, minDP_f, maxDP_f, disableDP, valueDP, minMQ, maxMQ, disableMQ, valueMQ, maxMiss, valueMiss, disableMiss, nSNPsbi, nSNPsmulti, nINDELs, nROHs, plots_options, df_json_list_avg_miss, df_json_list_dict_perc_gt, df_json_list_avg_dp_gt

## Display the datatable of the filtered VCF
@app.callback(Output('datatable-upload', 'children'),
              [Input('output-samples', 'children'),
              Input('sliderDP', 'value'),
              Input('sliderMQ', 'value'),
              Input('sliderMiss', 'value')],
              [State('radio-typeSNP', 'value')])
def showing_table(samples,valueDP,valueMQ,valueMiss,type_snps):

    if samples:

        ## Read the Json table in Pandas DataFrame format:
        #df = pd.read_json(jsonified_dataframe, orient='table')
        #df = feather.read_dataframe('df.feather')
        df = pd.read_feather('df.feather', use_threads=True)
        print('2:' + str(time.time()))

        ## Apply the selected ranges of the sliders:
        if type_snps != 'all':
            for s in samples:
                colnameDP = s + '_DP'
                colnameGT = s + '_GT'
                if colnameDP in df.columns:
                    df[colnameGT] = np.where(((df[colnameDP] >= valueDP[0]) & (df[colnameDP] <= valueDP[1])) | (df[colnameDP].notnull() == False), df[colnameGT], './.')

            if 'MQ' in df.columns:
                df = df.loc[(df['MQ'] >= valueMQ[0]) & (df['MQ'] <= valueMQ[1])]
    
        df = df.loc[(df['Missing'] <= valueMiss)]

        ## Convert the filtered DataFrame into Json table:
        #df_json = df.to_json(date_format='iso', orient='table')
            #feather.write_dataframe(df, 'df2.feather')
        df.reset_index().to_feather('df2.feather')
        print('3:' + str(time.time()))

        if (os.path.getsize('df.feather') < 100000000):
            return generate_table(df)#, df_json

## Display the summary statistics of the filtered VCF
@app.callback(Output('summary', 'children'),
              [Input('nsnpsbi', 'children'),
              Input('upload-poplist', 'filename'),
              Input('nsnpsmulti', 'children'),
              Input('nindels', 'children'),
              Input('nrohs', 'children'),
              Input('list_avg_miss', 'children'),
              Input('list_avg_dp_gt', 'children'),
              Input('list_dict_perc_gt', 'children')])
def summary_stats(nsnpsbi,popfile,nsnpsmulti,nindels,nrohs,json_list_avg_miss,json_list_avg_dp_gt,json_list_dict_perc_gt):

    if (nsnpsbi is not None) and (nsnpsmulti is not None) and (nindels is not None) and (nrohs is not None) and json_list_avg_miss and json_list_avg_dp_gt and json_list_dict_perc_gt:
        df_list_avg_miss = pd.read_json(json_list_avg_miss, orient='table')
        df_list_dict_perc_gt = pd.read_json(json_list_dict_perc_gt, orient='table')
        df_list_avg_dp_gt = pd.read_json(json_list_avg_dp_gt, orient='table')

        if df_list_avg_dp_gt.empty == False:
            if popfile and (df_list_dict_perc_gt.empty == False) and (df_list_avg_miss.empty == False):
                dict_perc_gt = df_list_dict_perc_gt.to_dict(orient='records')
                warning_proportion = ""
                warning_missing = ""
                warning_multi = ""
                if nsnpsmulti > nsnpsbi:
                    warning_multi = "Greater number of multiallelic SNPs than biallelic"

                warning_dpgt1 = ""
                warning_dpgt2 = ""
                for p in df_list_avg_dp_gt[(zscore(df_list_avg_dp_gt['Depth']) > 3)]['Genotype']:
                    if warning_dpgt1 == "":
                        warning_dpgt1 = "Greater depth than expected in genotype " + p
                    else:
                        warning_dpgt1 += ", " + p
                if warning_dpgt1 != "":
                    warning_dpgt1 += " (zscore > 3)"
                for p in df_list_avg_dp_gt[(zscore(df_list_avg_dp_gt['Depth']) < -3)]['Genotype']:
                    if warning_dpgt2 == "":
                        warning_dpgt2 = "Lower depth than expected in genotype " + p
                    else:
                        warning_dpgt2 += ", " + p
                if warning_dpgt2 != "":
                    warning_dpgt2 += " (zscore > 3)"

                warning_diff_prop = ""
                df_list_dict_perc_gt['0/0'] = df_list_dict_perc_gt['0/0'].str.replace('%','').astype(float)
                df_list_dict_perc_gt['0/1'] = df_list_dict_perc_gt['0/1'].str.replace('%','').astype(float)
                df_list_dict_perc_gt['1/1'] = df_list_dict_perc_gt['1/1'].str.replace('%','').astype(float)
                df_list_dict_perc_gt['0/2'] = df_list_dict_perc_gt['0/2'].str.replace('%','').astype(float)
                df_list_dict_perc_gt['1/2'] = df_list_dict_perc_gt['1/2'].str.replace('%','').astype(float)
                df_list_dict_perc_gt['2/2'] = df_list_dict_perc_gt['2/2'].str.replace('%','').astype(float)
                df_list_dict_perc_gt['0/3'] = df_list_dict_perc_gt['0/3'].str.replace('%','').astype(float)
                df_list_dict_perc_gt['1/3'] = df_list_dict_perc_gt['1/3'].str.replace('%','').astype(float)
                df_list_dict_perc_gt['2/3'] = df_list_dict_perc_gt['2/3'].str.replace('%','').astype(float)
                df_list_dict_perc_gt['3/3'] = df_list_dict_perc_gt['3/3'].str.replace('%','').astype(float)
                df_list_dict_perc_gt['./.'] = df_list_dict_perc_gt['./.'].str.replace('%','').astype(float)
                for p in df_list_dict_perc_gt[(np.abs(zscore(df_list_dict_perc_gt.iloc[:,1:11])) > 3)]['Population']:
                    if warning_diff_prop == "":
                        warning_diff_prop = "Different frequency in at least one genotype in the population/s " + p
                    else:
                        warning_diff_prop += ", " + p
                if warning_diff_prop != "":
                    warning_diff_prop += " (zscore > 3)"

                for i in range(0,len(dict_perc_gt)):
                    dict_perc_gt[i]['0/0'] = dict_perc_gt[i].pop('_1')
                    dict_perc_gt[i]['0/1'] = dict_perc_gt[i].pop('_2')
                    dict_perc_gt[i]['1/1'] = dict_perc_gt[i].pop('_3')
                    dict_perc_gt[i]['0/2'] = dict_perc_gt[i].pop('_4')
                    dict_perc_gt[i]['1/2'] = dict_perc_gt[i].pop('_5')
                    dict_perc_gt[i]['2/2'] = dict_perc_gt[i].pop('_6')
                    dict_perc_gt[i]['0/3'] = dict_perc_gt[i].pop('_7')
                    dict_perc_gt[i]['1/3'] = dict_perc_gt[i].pop('_8')
                    dict_perc_gt[i]['2/3'] = dict_perc_gt[i].pop('_9')
                    dict_perc_gt[i]['3/3'] = dict_perc_gt[i].pop('_10')
                    dict_perc_gt[i]['./.'] = dict_perc_gt[i].pop('_11')
                    if float(dict_perc_gt[i]['0/0'].replace('%', '')) - (float(dict_perc_gt[i]['1/1'].replace('%', ''))+float(dict_perc_gt[i]['1/2'].replace('%', ''))+float(dict_perc_gt[i]['2/2'].replace('%', ''))+float(dict_perc_gt[i]['1/3'].replace('%', ''))+float(dict_perc_gt[i]['2/3'].replace('%', ''))+float(dict_perc_gt[i]['3/3'].replace('%', ''))) < 0:
                        if warning_proportion == "":
                            warning_proportion += 'Higher proportion of alternative allele in ' + dict_perc_gt[i]['Population']
                        else:
                            warning_proportion += ', ' + dict_perc_gt[i]['Population']
                if warning_proportion != "":
                    warning_proportion += '!!'

                for p in df_list_avg_miss['Population']:
                    if float(df_list_avg_miss.set_index('Population').loc[p,'Missing'].replace('%', '')) > 30:
                        if warning_missing == "":
                            warning_missing += 'Missing in ' + p
                        else:
                            warning_missing += ', ' + p
                if warning_missing != "":
                    warning_missing += ' higher than 30%!!'

                return html.Div([
                            html.Div('- Number of biallelic SNPs = {}'.format(nsnpsbi)),
                            html.Div('- Number of multiallelic SNPs = {}'.format(nsnpsmulti)),
                            html.Div(warning_multi,style={'color': 'red','padding-left':'2%'}),
                            html.Div('- Number of INDELs = {}'.format(nindels)),
                            html.Div('- Number of homozigous positions = {}'.format(nrohs)),
                            html.Div('- Average missing per population:'),
                            html.Div([
                                dash_table.DataTable(
                                    data=df_list_avg_miss.to_dict(orient='records'),
                                    columns=[{'name': i, 'id': i} for i in df_list_avg_miss.columns],
                                    style_as_list_view=True,
                                    style_table={
                                        'maxWidth': '150',
                                        'maxHeight': '600',  
                                    },
                                ),
                            ], style={'marginTop':'1em', 'marginBottom':'1em', 'marginLeft':'2em'}),
                            html.Div(warning_missing,style={'color': 'red','padding-left':'2%'}),
                            html.Div('- Average depth by genotype:'),
                            html.Div([
                                dash_table.DataTable(
                                    data=df_list_avg_dp_gt.to_dict(orient='records'),
                                    columns=[{'name': i, 'id': i} for i in df_list_avg_dp_gt.columns],
                                    style_as_list_view=True,
                                    style_table={
                                        'maxWidth': '150',
                                        'maxHeight': '600',  
                                    },
                                ),
                            ], style={'marginTop':'1em', 'marginBottom':'1em', 'marginLeft':'2em'}),
                            html.Div(warning_dpgt1,style={'color': 'red','padding-left':'2%'}),
                            html.Div(warning_dpgt2,style={'color': 'red','padding-left':'2%'}),
                            html.Div('- Proportion of genotypes per population:'),
                            html.Div([
                                dash_table.DataTable(
                                    data=dict_perc_gt,
                                    columns=[{'name': i, 'id': i} for i in df_list_dict_perc_gt.columns],
                                    style_as_list_view=True,
                                    style_table={
                                        'maxWidth': '800',
                                        'maxHeight': '600',  
                                    },
                                ),
                            ], style={'marginTop':'1em', 'marginBottom':'1em', 'marginLeft':'2em'}),
                            html.Div(warning_proportion,style={'color': 'red','padding-left':'2%'}),
                            html.Div(warning_diff_prop,style={'color': 'red','padding-left':'2%'}),
                ])
            elif popfile and (df_list_dict_perc_gt.empty == False) and (df_list_avg_miss.empty == True):
                dict_perc_gt = df_list_dict_perc_gt.to_dict(orient='records')
                warning_proportion = ""
                warning_multi = ""
                if nsnpsmulti > nsnpsbi:
                    warning_multi = "Greater number of multiallelic SNPs than biallelic"

                warning_dpgt1 = ""
                warning_dpgt2 = ""
                for p in df_list_avg_dp_gt[(zscore(df_list_avg_dp_gt['Depth']) > 3)]['Genotype']:
                    if warning_dpgt1 == "":
                        warning_dpgt1 = "Greater depth than expected in genotype " + p
                    else:
                        warning_dpgt1 += ", " + p
                if warning_dpgt1 != "":
                    warning_dpgt1 += " (zscore > 3)"
                for p in df_list_avg_dp_gt[(zscore(df_list_avg_dp_gt['Depth']) < -3)]['Genotype']:
                    if warning_dpgt2 == "":
                        warning_dpgt2 = "Lower depth than expected in genotype " + p
                    else:
                        warning_dpgt2 += ", " + p
                if warning_dpgt2 != "":
                    warning_dpgt2 += " (zscore > 3)"

                warning_diff_prop = ""
                df_list_dict_perc_gt['0/0'] = df_list_dict_perc_gt['0/0'].str.replace('%','').astype(float)
                df_list_dict_perc_gt['0/1'] = df_list_dict_perc_gt['0/1'].str.replace('%','').astype(float)
                df_list_dict_perc_gt['1/1'] = df_list_dict_perc_gt['1/1'].str.replace('%','').astype(float)
                df_list_dict_perc_gt['0/2'] = df_list_dict_perc_gt['0/2'].str.replace('%','').astype(float)
                df_list_dict_perc_gt['1/2'] = df_list_dict_perc_gt['1/2'].str.replace('%','').astype(float)
                df_list_dict_perc_gt['2/2'] = df_list_dict_perc_gt['2/2'].str.replace('%','').astype(float)
                df_list_dict_perc_gt['0/3'] = df_list_dict_perc_gt['0/3'].str.replace('%','').astype(float)
                df_list_dict_perc_gt['1/3'] = df_list_dict_perc_gt['1/3'].str.replace('%','').astype(float)
                df_list_dict_perc_gt['2/3'] = df_list_dict_perc_gt['2/3'].str.replace('%','').astype(float)
                df_list_dict_perc_gt['3/3'] = df_list_dict_perc_gt['3/3'].str.replace('%','').astype(float)
                df_list_dict_perc_gt['./.'] = df_list_dict_perc_gt['./.'].str.replace('%','').astype(float)
                for p in df_list_dict_perc_gt[(np.abs(zscore(df_list_dict_perc_gt.iloc[:,1:11])) > 3)]['Population']:
                    if warning_diff_prop == "":
                        warning_diff_prop = "Different frequency in at least one genotype in the population/s " + p
                    else:
                        warning_diff_prop += ", " + p
                if warning_diff_prop != "":
                    warning_diff_prop += " (zscore > 3)"

                for i in range(0,len(dict_perc_gt)):
                    dict_perc_gt[i]['0/0'] = dict_perc_gt[i].pop('_1')
                    dict_perc_gt[i]['0/1'] = dict_perc_gt[i].pop('_2')
                    dict_perc_gt[i]['1/1'] = dict_perc_gt[i].pop('_3')
                    dict_perc_gt[i]['0/2'] = dict_perc_gt[i].pop('_4')
                    dict_perc_gt[i]['1/2'] = dict_perc_gt[i].pop('_5')
                    dict_perc_gt[i]['2/2'] = dict_perc_gt[i].pop('_6')
                    dict_perc_gt[i]['0/3'] = dict_perc_gt[i].pop('_7')
                    dict_perc_gt[i]['1/3'] = dict_perc_gt[i].pop('_8')
                    dict_perc_gt[i]['2/3'] = dict_perc_gt[i].pop('_9')
                    dict_perc_gt[i]['3/3'] = dict_perc_gt[i].pop('_10')
                    dict_perc_gt[i]['./.'] = dict_perc_gt[i].pop('_11')
                    if float(dict_perc_gt[i]['0/0'].replace('%', '')) - (float(dict_perc_gt[i]['1/1'].replace('%', ''))+float(dict_perc_gt[i]['1/2'].replace('%', ''))+float(dict_perc_gt[i]['2/2'].replace('%', ''))+float(dict_perc_gt[i]['1/3'].replace('%', ''))+float(dict_perc_gt[i]['2/3'].replace('%', ''))+float(dict_perc_gt[i]['3/3'].replace('%', ''))) < 0:
                        if warning_proportion == "":
                            warning_proportion += 'Higher proportion of alternative allele in ' + dict_perc_gt[i]['Population']
                        else:
                            warning_proportion += ', ' + dict_perc_gt[i]['Population']
                if warning_proportion != "":
                    warning_proportion += '!!'

                return html.Div([
                            html.Div('- Number of biallelic SNPs = {}'.format(nsnpsbi)),
                            html.Div('- Number of multiallelic SNPs = {}'.format(nsnpsmulti)),
                            html.Div(warning_multi,style={'color': 'red','padding-left':'2%'}),
                            html.Div('- Number of INDELs = {}'.format(nindels)),
                            html.Div('- Number of homozigous positions = {}'.format(nrohs)),
                            html.Div('- Average depth by genotype:'),
                            html.Div([
                                dash_table.DataTable(
                                    data=df_list_avg_dp_gt.to_dict(orient='records'),
                                    columns=[{'name': i, 'id': i} for i in df_list_avg_dp_gt.columns],
                                    style_as_list_view=True,
                                    style_table={
                                        'maxWidth': '150',
                                        'maxHeight': '600',  
                                    },
                                ),
                            ], style={'marginTop':'1em', 'marginBottom':'1em', 'marginLeft':'2em'}),
                            html.Div(warning_dpgt1,style={'color': 'red','padding-left':'2%'}),
                            html.Div(warning_dpgt2,style={'color': 'red','padding-left':'2%'}),
                            html.Div('- Proportion of genotypes per population:'),
                            html.Div([
                                dash_table.DataTable(
                                    data=dict_perc_gt,
                                    columns=[{'name': i, 'id': i} for i in df_list_dict_perc_gt.columns],
                                    style_as_list_view=True,
                                    style_table={
                                        'maxWidth': '800',
                                        'maxHeight': '600',  
                                    },
                                ),
                            ], style={'marginTop':'1em', 'marginBottom':'1em', 'marginLeft':'2em'}),
                            html.Div(warning_proportion,style={'color': 'red','padding-left':'2%'}),
                            html.Div(warning_diff_prop,style={'color': 'red','padding-left':'2%'}),
                ])
            elif popfile and (df_list_dict_perc_gt.empty == True) and (df_list_avg_miss.empty == False):
                warning_missing = ""
                warning_multi = ""
                if nsnpsmulti > nsnpsbi:
                    warning_multi = "Greater number of multiallelic SNPs than biallelic"

                for p in df_list_avg_miss['Population']:
                    if float(df_list_avg_miss.set_index('Population').loc[p,'Missing'].replace('%', '')) > 30:
                        if warning_missing == "":
                            warning_missing += 'Missing in ' + p
                        else:
                            warning_missing += ', ' + p
                if warning_missing != "":
                    warning_missing += ' higher than 30%!!'

                warning_dpgt1 = ""
                warning_dpgt2 = ""
                for p in df_list_avg_dp_gt[(zscore(df_list_avg_dp_gt['Depth']) > 3)]['Genotype']:
                    if warning_dpgt1 == "":
                        warning_dpgt1 = "Greater depth than expected in genotype " + p
                    else:
                        warning_dpgt1 += ", " + p
                if warning_dpgt1 != "":
                    warning_dpgt1 += " (zscore > 3)"
                for p in df_list_avg_dp_gt[(zscore(df_list_avg_dp_gt['Depth']) < -3)]['Genotype']:
                    if warning_dpgt2 == "":
                        warning_dpgt2 = "Lower depth than expected in genotype " + p
                    else:
                        warning_dpgt2 += ", " + p
                if warning_dpgt2 != "":
                    warning_dpgt2 += " (zscore > 3)"

                return html.Div([
                            html.Div('- Number of biallelic SNPs = {}'.format(nsnpsbi)),
                            html.Div('- Number of multiallelic SNPs = {}'.format(nsnpsmulti)),
                            html.Div(warning_multi,style={'color': 'red','padding-left':'2%'}),
                            html.Div('- Number of INDELs = {}'.format(nindels)),
                            html.Div('- Number of homozigous positions = {}'.format(nrohs)),
                            html.Div('- Average missing per population:'),
                            html.Div([
                                dash_table.DataTable(
                                    data=df_list_avg_miss.to_dict(orient='records'),
                                    columns=[{'name': i, 'id': i} for i in df_list_avg_miss.columns],
                                    style_as_list_view=True,
                                    style_table={
                                        'maxWidth': '150',
                                        'maxHeight': '600',  
                                    },
                                ),
                            ], style={'marginTop':'1em', 'marginBottom':'1em', 'marginLeft':'2em'}),
                            html.Div(warning_missing,style={'color': 'red','padding-left':'2%'}),
                            html.Div('- Average depth by genotype:'),
                            html.Div([
                                dash_table.DataTable(
                                    data=df_list_avg_dp_gt.to_dict(orient='records'),
                                    columns=[{'name': i, 'id': i} for i in df_list_avg_dp_gt.columns],
                                    style_as_list_view=True,
                                    style_table={
                                        'maxWidth': '150',
                                        'maxHeight': '600',  
                                    },
                                ),
                            ], style={'marginTop':'1em', 'marginBottom':'1em', 'marginLeft':'2em'}),
                            html.Div(warning_dpgt1,style={'color': 'red','padding-left':'2%'}),
                            html.Div(warning_dpgt2,style={'color': 'red','padding-left':'2%'}),
                ])

            else:

                warning_dpgt1 = ""
                warning_dpgt2 = ""
                for p in df_list_avg_dp_gt[(zscore(df_list_avg_dp_gt['Depth']) > 3)]['Genotype']:
                    if warning_dpgt1 == "":
                        warning_dpgt1 = "Greater depth than expected in genotype " + p
                    else:
                        warning_dpgt1 += ", " + p
                if warning_dpgt1 != "":
                    warning_dpgt1 += " (zscore > 3)"
                for p in df_list_avg_dp_gt[(zscore(df_list_avg_dp_gt['Depth']) < -3)]['Genotype']:
                    if warning_dpgt2 == "":
                        warning_dpgt2 = "Lower depth than expected in genotype " + p
                    else:
                        warning_dpgt2 += ", " + p
                if warning_dpgt2 != "":
                    warning_dpgt2 += " (zscore > 3)"

                return html.Div([
                            html.Div('- Number of biallelic SNPs = {}'.format(nsnpsbi)),
                            html.Div('- Number of multiallelic SNPs = {}'.format(nsnpsmulti)),
                            html.Div(warning_multi,style={'color': 'red','padding-left':'2%'}),
                            html.Div('- Number of INDELs = {}'.format(nindels)),
                            html.Div('- Number of homozigous positions = {}'.format(nrohs)),
                            html.Div('- Average depth by genotype:'),
                            html.Div([
                                dash_table.DataTable(
                                    data=df_list_avg_dp_gt.to_dict(orient='records'),
                                    columns=[{'name': i, 'id': i} for i in df_list_avg_dp_gt.columns],
                                    style_as_list_view=True,
                                    style_table={
                                        'maxWidth': '150',
                                        'maxHeight': '600',  
                                    },
                                ),
                            ], style={'marginTop':'1em', 'marginBottom':'1em', 'marginLeft':'2em'}),
                            html.Div(warning_dpgt1,style={'color': 'red','padding-left':'2%'}),
                            html.Div(warning_dpgt2,style={'color': 'red','padding-left':'2%'}),
                ])
    
        else:
            if popfile and (df_list_dict_perc_gt.empty == False) and (df_list_avg_miss.empty == False):
                dict_perc_gt = df_list_dict_perc_gt.to_dict(orient='records')
                warning_proportion = ""
                warning_missing = ""
                warning_multi = ""
                if nsnpsmulti > nsnpsbi:
                    warning_multi = "Greater number of multiallelic SNPs than biallelic"

                warning_diff_prop = ""
                df_list_dict_perc_gt['0/0'] = df_list_dict_perc_gt['0/0'].str.replace('%','').astype(float)
                df_list_dict_perc_gt['0/1'] = df_list_dict_perc_gt['0/1'].str.replace('%','').astype(float)
                df_list_dict_perc_gt['1/1'] = df_list_dict_perc_gt['1/1'].str.replace('%','').astype(float)
                df_list_dict_perc_gt['0/2'] = df_list_dict_perc_gt['0/2'].str.replace('%','').astype(float)
                df_list_dict_perc_gt['1/2'] = df_list_dict_perc_gt['1/2'].str.replace('%','').astype(float)
                df_list_dict_perc_gt['2/2'] = df_list_dict_perc_gt['2/2'].str.replace('%','').astype(float)
                df_list_dict_perc_gt['0/3'] = df_list_dict_perc_gt['0/3'].str.replace('%','').astype(float)
                df_list_dict_perc_gt['1/3'] = df_list_dict_perc_gt['1/3'].str.replace('%','').astype(float)
                df_list_dict_perc_gt['2/3'] = df_list_dict_perc_gt['2/3'].str.replace('%','').astype(float)
                df_list_dict_perc_gt['3/3'] = df_list_dict_perc_gt['3/3'].str.replace('%','').astype(float)
                df_list_dict_perc_gt['./.'] = df_list_dict_perc_gt['./.'].str.replace('%','').astype(float)
                for p in df_list_dict_perc_gt[(np.abs(zscore(df_list_dict_perc_gt.iloc[:,1:11])) > 3)]['Population']:
                    if warning_diff_prop == "":
                        warning_diff_prop = "Different frequency in at least one genotype in the population/s " + p
                    else:
                        warning_diff_prop += ", " + p
                if warning_diff_prop != "":
                    warning_diff_prop += " (zscore > 3)"

                for i in range(0,len(dict_perc_gt)):
                    dict_perc_gt[i]['0/0'] = dict_perc_gt[i].pop('_1')
                    dict_perc_gt[i]['0/1'] = dict_perc_gt[i].pop('_2')
                    dict_perc_gt[i]['1/1'] = dict_perc_gt[i].pop('_3')
                    dict_perc_gt[i]['0/2'] = dict_perc_gt[i].pop('_4')
                    dict_perc_gt[i]['1/2'] = dict_perc_gt[i].pop('_5')
                    dict_perc_gt[i]['2/2'] = dict_perc_gt[i].pop('_6')
                    dict_perc_gt[i]['0/3'] = dict_perc_gt[i].pop('_7')
                    dict_perc_gt[i]['1/3'] = dict_perc_gt[i].pop('_8')
                    dict_perc_gt[i]['2/3'] = dict_perc_gt[i].pop('_9')
                    dict_perc_gt[i]['3/3'] = dict_perc_gt[i].pop('_10')
                    dict_perc_gt[i]['./.'] = dict_perc_gt[i].pop('_11')
                    if float(dict_perc_gt[i]['0/0'].replace('%', '')) - (float(dict_perc_gt[i]['1/1'].replace('%', ''))+float(dict_perc_gt[i]['1/2'].replace('%', ''))+float(dict_perc_gt[i]['2/2'].replace('%', ''))+float(dict_perc_gt[i]['1/3'].replace('%', ''))+float(dict_perc_gt[i]['2/3'].replace('%', ''))+float(dict_perc_gt[i]['3/3'].replace('%', ''))) < 0:
                        if warning_proportion == "":
                            warning_proportion += 'Higher proportion of alternative allele in ' + dict_perc_gt[i]['Population']
                        else:
                            warning_proportion += ', ' + dict_perc_gt[i]['Population']
                if warning_proportion != "":
                    warning_proportion += '!!'

                for p in df_list_avg_miss['Population']:
                    if float(df_list_avg_miss.set_index('Population').loc[p,'Missing'].replace('%', '')) > 30:
                        if warning_missing == "":
                            warning_missing += 'Missing in ' + p
                        else:
                            warning_missing += ', ' + p
                if warning_missing != "":
                    warning_missing += ' higher than 30%!!'

                return html.Div([
                            html.Div('- Number of biallelic SNPs = {}'.format(nsnpsbi)),
                            html.Div('- Number of multiallelic SNPs = {}'.format(nsnpsmulti)),
                            html.Div(warning_multi,style={'color': 'red','padding-left':'2%'}),
                            html.Div('- Number of INDELs = {}'.format(nindels)),
                            html.Div('- Number of homozigous positions = {}'.format(nrohs)),
                            html.Div('- Average missing per population:'),
                            html.Div([
                                dash_table.DataTable(
                                    data=df_list_avg_miss.to_dict(orient='records'),
                                    columns=[{'name': i, 'id': i} for i in df_list_avg_miss.columns],
                                    style_as_list_view=True,
                                    style_table={
                                        'maxWidth': '150',
                                        'maxHeight': '600',  
                                    },
                                ),
                            ], style={'marginTop':'1em', 'marginBottom':'1em', 'marginLeft':'2em'}),
                            html.Div(warning_missing,style={'color': 'red','padding-left':'2%'}),
                            html.Div('- Proportion of genotypes per population:'),
                            html.Div([
                                dash_table.DataTable(
                                    data=dict_perc_gt,
                                    columns=[{'name': i, 'id': i} for i in df_list_dict_perc_gt.columns],
                                    style_as_list_view=True,
                                    style_table={
                                        'maxWidth': '800',
                                        'maxHeight': '600',  
                                    },
                                ),
                            ], style={'marginTop':'1em', 'marginLeft':'2em'}),
                            html.Div(warning_proportion,style={'color': 'red','padding-left':'2%'}),
                            html.Div(warning_diff_prop,style={'color': 'red','padding-left':'2%'}),
                ])
            elif popfile and (df_list_dict_perc_gt.empty == False) and (df_list_avg_miss.empty == True):
                dict_perc_gt = df_list_dict_perc_gt.to_dict(orient='records')
                warning_proportion = ""
                warning_multi = ""
                if nsnpsmulti > nsnpsbi:
                    warning_multi = "Greater number of multiallelic SNPs than biallelic"

                warning_diff_prop = ""
                df_list_dict_perc_gt['0/0'] = df_list_dict_perc_gt['0/0'].str.replace('%','').astype(float)
                df_list_dict_perc_gt['0/1'] = df_list_dict_perc_gt['0/1'].str.replace('%','').astype(float)
                df_list_dict_perc_gt['1/1'] = df_list_dict_perc_gt['1/1'].str.replace('%','').astype(float)
                df_list_dict_perc_gt['0/2'] = df_list_dict_perc_gt['0/2'].str.replace('%','').astype(float)
                df_list_dict_perc_gt['1/2'] = df_list_dict_perc_gt['1/2'].str.replace('%','').astype(float)
                df_list_dict_perc_gt['2/2'] = df_list_dict_perc_gt['2/2'].str.replace('%','').astype(float)
                df_list_dict_perc_gt['0/3'] = df_list_dict_perc_gt['0/3'].str.replace('%','').astype(float)
                df_list_dict_perc_gt['1/3'] = df_list_dict_perc_gt['1/3'].str.replace('%','').astype(float)
                df_list_dict_perc_gt['2/3'] = df_list_dict_perc_gt['2/3'].str.replace('%','').astype(float)
                df_list_dict_perc_gt['3/3'] = df_list_dict_perc_gt['3/3'].str.replace('%','').astype(float)
                df_list_dict_perc_gt['./.'] = df_list_dict_perc_gt['./.'].str.replace('%','').astype(float)
                for p in df_list_dict_perc_gt[(np.abs(zscore(df_list_dict_perc_gt.iloc[:,1:11])) > 3)]['Population']:
                    if warning_diff_prop == "":
                        warning_diff_prop = "Different frequency in at least one genotype in the population/s " + p
                    else:
                        warning_diff_prop += ", " + p
                if warning_diff_prop != "":
                    warning_diff_prop += " (zscore > 3)"

                for i in range(0,len(dict_perc_gt)):
                    dict_perc_gt[i]['0/0'] = dict_perc_gt[i].pop('_1')
                    dict_perc_gt[i]['0/1'] = dict_perc_gt[i].pop('_2')
                    dict_perc_gt[i]['1/1'] = dict_perc_gt[i].pop('_3')
                    dict_perc_gt[i]['0/2'] = dict_perc_gt[i].pop('_4')
                    dict_perc_gt[i]['1/2'] = dict_perc_gt[i].pop('_5')
                    dict_perc_gt[i]['2/2'] = dict_perc_gt[i].pop('_6')
                    dict_perc_gt[i]['0/3'] = dict_perc_gt[i].pop('_7')
                    dict_perc_gt[i]['1/3'] = dict_perc_gt[i].pop('_8')
                    dict_perc_gt[i]['2/3'] = dict_perc_gt[i].pop('_9')
                    dict_perc_gt[i]['3/3'] = dict_perc_gt[i].pop('_10')
                    dict_perc_gt[i]['./.'] = dict_perc_gt[i].pop('_11')
                    if float(dict_perc_gt[i]['0/0'].replace('%', '')) - (float(dict_perc_gt[i]['1/1'].replace('%', ''))+float(dict_perc_gt[i]['1/2'].replace('%', ''))+float(dict_perc_gt[i]['2/2'].replace('%', ''))+float(dict_perc_gt[i]['1/3'].replace('%', ''))+float(dict_perc_gt[i]['2/3'].replace('%', ''))+float(dict_perc_gt[i]['3/3'].replace('%', ''))) < 0:
                        if warning_proportion == "":
                            warning_proportion += 'Higher proportion of alternative allele in ' + dict_perc_gt[i]['Population']
                        else:
                            warning_proportion += ', ' + dict_perc_gt[i]['Population']
                if warning_proportion != "":
                    warning_proportion += '!!'

                return html.Div([
                            html.Div('- Number of biallelic SNPs = {}'.format(nsnpsbi)),
                            html.Div('- Number of multiallelic SNPs = {}'.format(nsnpsmulti)),
                            html.Div(warning_multi,style={'color': 'red','padding-left':'2%'}),
                            html.Div('- Number of INDELs = {}'.format(nindels)),
                            html.Div('- Number of homozigous positions = {}'.format(nrohs)),
                            html.Div('- Proportion of genotypes per population:'),
                            html.Div([
                                dash_table.DataTable(
                                    data=dict_perc_gt,
                                    columns=[{'name': i, 'id': i} for i in df_list_dict_perc_gt.columns],
                                    style_as_list_view=True,
                                    style_table={
                                        'maxWidth': '800',
                                        'maxHeight': '600',  
                                    },
                                ),
                            ], style={'marginTop':'1em', 'marginLeft':'2em'}),
                            html.Div(warning_proportion,style={'color': 'red','padding-left':'2%'}),
                            html.Div(warning_diff_prop,style={'color': 'red','padding-left':'2%'}),
                ])
            elif popfile and (df_list_dict_perc_gt.empty == True) and (df_list_avg_miss.empty == False):
                warning_missing = ""
                warning_multi = ""
                if nsnpsmulti > nsnpsbi:
                    warning_multi = "Greater number of multiallelic SNPs than biallelic"

                for p in df_list_avg_miss['Population']:
                    if float(df_list_avg_miss.set_index('Population').loc[p,'Missing'].replace('%', '')) > 30:
                        if warning_missing == "":
                            warning_missing += 'Missing in ' + p
                        else:
                            warning_missing += ', ' + p
                if warning_missing != "":
                    warning_missing += ' higher than 30%!!'

                return html.Div([
                            html.Div('- Number of biallelic SNPs = {}'.format(nsnpsbi)),
                            html.Div('- Number of multiallelic SNPs = {}'.format(nsnpsmulti)),
                            html.Div(warning_multi,style={'color': 'red','padding-left':'2%'}),
                            html.Div('- Number of INDELs = {}'.format(nindels)),
                            html.Div('- Number of homozigous positions = {}'.format(nrohs)),
                            html.Div('- Average missing per population:'),
                            html.Div([
                                dash_table.DataTable(
                                    data=df_list_avg_miss.to_dict(orient='records'),
                                    columns=[{'name': i, 'id': i} for i in df_list_avg_miss.columns],
                                    style_as_list_view=True,
                                    style_table={
                                        'maxWidth': '150',
                                        'maxHeight': '600',  
                                    },
                                ),
                            ], style={'marginTop':'1em', 'marginBottom':'1em', 'marginLeft':'2em'}),
                            html.Div(warning_missing,style={'color': 'red','padding-left':'2%'}),
                ])
            else:
                warning_multi = ""
                if nsnpsmulti > nsnpsbi:
                    warning_multi = "Greater number of multiallelic SNPs than biallelic"
                return html.Div([
                            html.Div('- Number of biallelic SNPs = {}'.format(nsnpsbi)),
                            html.Div('- Number of multiallelic SNPs = {}'.format(nsnpsmulti)),
                            html.Div(warning_multi,style={'color': 'red','padding-left':'2%'}),
                            html.Div('- Number of INDELs = {}'.format(nindels)),
                            html.Div('- Number of homozigous positions = {}'.format(nrohs)),
                ])

## Download the filtered VCF:
@app.callback(Output('output-download', 'children'),
             [Input('button-download', 'n_clicks')],
             [State('output-samples', 'children'),
             State('header', 'children'),
             State('output-name','children')])
def vcf_download_link(nclicks,samples,header,filename):

    if nclicks > 0:
        
        ## For the output name we use the input name + ".filtered.vcf":
        outputfile = filename.split(' ')[-1]
        f = open(outputfile + '.filtered.vcf', 'w')

        ## Write the header in the output VCF:
        f.write('## Filtered by VCFcheck' + '\n')
        for h in header:
            f.write(h + '\n')

        ## Read the Json table in Pandas DataFrame format:
        #df = pd.read_json(jsonified_dataframe, orient='table')
        #df = feather.read_dataframe('df2.feather')
        df = pd.read_feather('df2.feather', use_threads=True)

        ## Select the standard columns of a VCF to write them in the output VCF:
        columns_str = [col for col in ['CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT'] + samples]
        df = df[columns_str].rename(columns={'CHROM': '#CHROM'})
        csv_string = df.to_csv(sep='\t', index=False, encoding='utf-8',na_rep='NA')
        f.write(csv_string)
        #csv_string = "data:text/csv;charset=utf-8," + urllib.parse.quote(bytes(csv_string))

        return outputfile + '.filtered.vcf' + ' downloaded!'

## Download the filtered VCF:
@app.callback(Output('output-summary', 'children'),
             [Input('button-summary', 'n_clicks')],
             [State('summary','children'),
             State('output-name','children')])
def summary_download_link(nclicks,summary,filename):

    if nclicks > 0:
        
        ## For the output name we use the input name + ".summary.txt":
        outputfile = filename.split(' ')[-1]

        # We create the file with the header (if exists, this replace it)
        f = open(outputfile + '.summary.txt', 'w+')
        f.write('## Summary from VCFcheck' + '\n')
        f.close()

        ## Write the summary in the existing file with the header:
        f = open(outputfile + '.summary.txt', 'a')

        for c in summary['props']['children']:
            if type(c['props']['children']) == list:
                f.write('\n\n') 
                f.write(pd.DataFrame(c['props']['children'][0]['props']['data']).to_csv(sep='\t', index=False, encoding='utf-8',na_rep='NA'))
            else:
                f.write('\n\n')
                f.write(c['props']['children'])

        f.close()

        return outputfile + '.summary.txt' + ' downloaded!'

## Print value of the selected range of DP next to the Range Slider:
@app.callback(Output('output-container-range-sliderDP', 'children'),
              [Input('sliderDP', 'value')])
def update_output1(value):
    range_selected = 'Sample depth (DP) selected between ' + str(value[0]) + ' and ' + str(value[1])
    return range_selected

## Print value of the selected range of MQ next to the Range Slider:
@app.callback(Output('output-container-range-sliderMQ', 'children'),
              [Input('sliderMQ', 'value')])
def update_output2(value):
    range_selected = 'Average mapping quality (MQ) selected between ' + str(value[0]) + ' and ' + str(value[1])
    return range_selected

## Print value of the selected range of Missing next to the Range Slider:
@app.callback(Output('output-container-range-sliderMiss', 'children'),
              [Input('sliderMiss', 'value')])
def update_output3(value):
    range_selected = 'Max. percentage of missing data by SNP: ' + str(value) + '%'
    return range_selected

## Display the selected plot in the dropdown
@app.callback(Output('plot-graph', 'children'),
              [Input('plot-selector', 'value')],
              [State('output-samples', 'children'),
              State('upload-poplist', 'filename'),
              State('sliderDP', 'value')])
def load_distribution_graph(plot, samples, popfile, valueDP):

    if plot is not None:

        ## Read the Json table in Pandas DataFrame format:
        #df = pd.read_json(jsonified_dataframe, orient='table')
        #df = feather.read_dataframe('df2.feather')
        df = pd.read_feather('df2.feather', use_threads=True)

        if plot == "Missing by SNP":

            ## Density plot of the distribution of missing data by SNP
            if (df.Missing.max() - df.Missing.min()) > 0:
                fig = draw_missing_snp(plot, df['Missing'])
                if fig:
                    fig['layout'].update(title = 'Missing by SNP', xaxis=dict(title='Percentage of missing data'))
                    return html.Div([dcc.Graph(figure=fig)])
                else:
                    layout = dict(title = 'Missing by SNP', xaxis=dict(title='Percentage of missing data'))
                    return html.Div([dcc.Graph(figure=dict(layout=layout)), html.Div('It is not possible to display the distribution of missing',style={'color': 'red','padding-left':'3%','padding-bottom':'2%'})])
            else:
                layout = dict(title = 'Missing by SNP', xaxis=dict(title='Percentage of missing data'))
                return html.Div([dcc.Graph(figure=dict(layout=layout)), html.Div('It is not possible to display the distribution of missing because all SNPs have the same missing (' + str(int(df.Missing.max())) + '%)',style={'color': 'red','padding-left':'3%','padding-bottom':'2%'})])

        elif plot == "Missing by Population":

            ## Density plot of the distribution of missing data by population (it is necessary the input file of individuals-populations)
            if popfile is not None:
                if (df.Missing.max() - df.Missing.min()) > 0:
                    missing_individual = missing_ind(df, samples)
                    indiv_table = table_indiv(popfile,samples)
                    indiv_table.index = indiv_table.Individual
                    indiv_table = pd.concat([indiv_table, missing_individual], axis=1, sort=False, join='inner')
                    indiv_table.rename(columns = {0:"Missing"}, inplace=True)
            
                    tablepop = indiv_table.replace('./.', np.NaN)
                    del indiv_table
    
                    if len(tablepop) > 0:
                        fig = draw_missing_pop(plot, tablepop)
                    if fig:
                        fig['layout'].update(title = 'Missing by Population', xaxis=dict(title='Percentage of missing data'))
                        return html.Div([dcc.Graph(figure=fig)])
                    else:
                        layout = dict(title = 'Missing by Population', xaxis=dict(title='Percentage of missing data'))     
                        return html.Div([dcc.Graph(figure=dict(layout=layout)), html.Div('It is not possible to display the distribution of missing',style={'color': 'red','padding-left':'3%','padding-bottom':'2%'})])
                else:
                    layout = dict(title = 'Missing by Population', xaxis=dict(title='Percentage of missing data'))     
                    return html.Div([dcc.Graph(figure=dict(layout=layout)), html.Div('It is not possible to display the distribution of missing because all SNPs have the same missing (' + str(int(df.Missing.max())) + '%)',style={'color': 'red','padding-left':'3%','padding-bottom':'2%'})])
            else:
                layout = dict(title = 'Missing by Population', xaxis=dict(title='Percentage of missing data'))     
                return html.Div([dcc.Graph(figure=dict(layout=layout)), html.Div('It is not possible to display the distribution of missing because the Population file is not uploaded',style={'color': 'red','padding-left':'3%','padding-bottom':'2%'})])

        elif plot == "Reference allele frequency":

            ## Density plot of the distribution of the reference allele frequency

            vcf_samples = pd.DataFrame()

            for s in samples:
                s_gt = s + '_GT'

                if s_gt in df.columns:
                    vcf_samples[s_gt] = df[s_gt].astype('str').replace({ "\.\/\." : np.NaN, "0\/0" : 1, "0\/." : 0.5, ".\/." : 0}, regex=True)
                else:
                    vcf_samples[s] = df[s].astype('str').replace({ "\.\/\." : np.NaN, "0\/0" : 1, "0\/." : 0.5, ".\/." : 0}, regex=True)

            vcf_samples = vcf_samples.T

            reffreq = vcf_samples.sum()/len(samples)

            del vcf_samples

            fig = draw_raf(plot, reffreq)

            del reffreq

            if fig:
                fig['layout'].update(title = 'Reference allele frequency', xaxis=dict(title='Frequency'))
                return html.Div([dcc.Graph(figure=fig)])
            else:
                layout = dict(title = 'Reference allele frequency', xaxis=dict(title='Frequency'))     
                return html.Div([dcc.Graph(figure=dict(layout=layout)), html.Div('It is not possible to display the distribution of the reference allele frequency',style={'color': 'red','padding-left':'3%','padding-bottom':'2%'})])

        elif plot == "Depth by individual":

            ## Density plot of the depth distribution by individual

            fig = draw_depth_ind(plot, df, samples)
            if fig:
                fig['layout'].update(title = 'Depth by individual', xaxis=dict(title='Depth'))
                return html.Div([dcc.Graph(figure=fig)])
            else:
                layout = dict(title = 'Depth by individual', xaxis=dict(title='Depth'))
                return html.Div([dcc.Graph(figure=dict(layout=layout)), html.Div('It is not possible to display the distribution of depth for lack of information',style={'color': 'red','padding-left':'3%','padding-bottom':'2%'})])

        elif plot == "Depth by genotype":

            ## Density plot of the depth distribution by genotype

            fig = draw_depth_gt(plot, df, samples)
            if fig:
                fig['layout'].update(title = 'Depth by genotype', xaxis=dict(title='Depth'))
                return html.Div([dcc.Graph(figure=fig)])
            else:
                layout = dict(title = 'Depth by genotype', xaxis=dict(title='Depth'))
                return html.Div([dcc.Graph(figure=dict(layout=layout)), html.Div('It is not possible to display the distribution of depth for lack of information',style={'color': 'red','padding-left':'3%','padding-bottom':'2%'})])

        elif plot == "Count genotypes per Depth":

            ## Count number of each genotype by depth
            data = draw_depth_vs_gt(plot, df, samples, valueDP)
            layout = go.Layout(barmode='stack', title = 'Count of genotypes per sample depth', xaxis=dict(title='Sample depth',type='category'), yaxis=dict(title='Count of genotypes'))
            if data:
                return html.Div([dcc.Graph(figure=dict(data=[n for n in data], layout=layout))])
            else:
                return html.Div([dcc.Graph(figure=dict(layout=layout)), html.Div('It is not possible to display the distribution of depth for lack of information',style={'color': 'red','padding-left':'3%','padding-bottom':'2%'})])

        elif plot == "Principal Component Analyisis (PCA)":

            ## Representation of the PCA

            if popfile is not None:

                ## First of all, we change the genotypes by a weight (0/0 = 0, 0/1 = 0.5 y 1/1 = 1),
                ## if the genotype is missing (./.) we replace it by the intermediate frequency of the SNP.
                ## We calculated this intermediate frequency as VCFtools, which do it as:
                ##    1. Sum the weights of each individual by SNP and divide by the total number of non-missing SNPs (intermediate frequency using only non-missing SNPs)
                ##    2. Replace the missing genotype by this intermediate frequency
                ##    3. Recalculate the intermediate frequency by the sum of all frequencies (using the previous intermediate frequency) and dividing by the total number of SNPs

                vcf_samples = pd.DataFrame()

                for s in samples:
                    s_gt = s + '_GT'
    
                    if s_gt in df.columns:
                        vcf_samples[s] = df[s_gt].astype('str').replace({ "\.\/\." : np.NaN, "0\/0" : 1, "0\/." : 0.5, ".\/." : 0}, regex=True)
                    else:
                        vcf_samples[s] = df[s].astype('str').replace({ "\.\/\." : np.NaN, "0\/0" : 1, "0\/." : 0.5, ".\/." : 0}, regex=True)

                vcf_samples.dropna(how = 'all',inplace=True)

                freq_nomiss = vcf_samples.T.sum()/vcf_samples.T.notnull().sum()
                nnull = vcf_samples.T.isnull().sum()
                weight_miss = (freq_nomiss + (freq_nomiss*nnull))/len(samples)
    
                vcf_samples.fillna(value = "NA", inplace = True)

                for s in samples:
                    vcf_samples[s] = np.where(vcf_samples[s] == "NA", weight_miss, vcf_samples[s])

                ## Transpose the table to have the individuals as the Row names (index of the table) and Add a column with the population of each individual:
                indiv_table = table_indiv(popfile,samples)
                indiv_table.set_index('Individual', inplace=True)

                vcf_samples_T = pd.concat([vcf_samples.T,indiv_table.Population], axis=1, sort=False, join='inner')
                vcf_samples_T.reset_index(drop=True,inplace=True)

                ## Principal Component Analysis:

                x = vcf_samples_T.iloc[:, 0:-1].astype('float') # Separating out the features (frequencies)
                y = vcf_samples_T.loc[:,['Population']].values.astype('str') # Separating out the target (populations)

                x = StandardScaler().fit_transform(x) # Standardizing the features

                pca = PCA(n_components = 2)
                principalComponents = pca.fit_transform(x)

                principalDf = pd.DataFrame(data = principalComponents, columns = ['principal component 1', 'principal component 2'])
                finalDf = pd.concat([principalDf, vcf_samples_T[['Population']]], axis = 1)

                layout = dict(title = 'PCA decomposition', xaxis=dict(title='PC 1'), yaxis=dict(title='PC 2'))
                data = draw_pca(plot, finalDf)

                if data:
                    return html.Div([dcc.Graph(figure=dict(data=[n for n in data], layout=layout))])
                else:
                    return html.Div([dcc.Graph(figure=dict(layout=layout)), html.Div('It is not possible to perform the PCA',style={'color': 'red','padding-left':'3%','padding-bottom':'2%'})])

        elif plot == "Hardy-Weinberg test":

            ## Plot of the Hardy Weinberg p-Value distribution

            if popfile is not None:

                hwe = pd.DataFrame()

                indiv_table = table_indiv(popfile,samples)
                pops = indiv_table['Population'].unique()

                warning_hw = ""

                for pop in pops:
                    if samples[0] + '_GT' in df.columns:
                        counts_gt_snps = df[indiv_table[indiv_table.Population == pop].Individual + '_GT'].replace('./.', np.NaN).apply(pd.value_counts,axis=1).fillna(0)
                    else:
                        counts_gt_snps = df[indiv_table[indiv_table.Population == pop].Individual].replace('./.', np.NaN).apply(pd.value_counts,axis=1).fillna(0)
                    counts_gt_snps['total'] = counts_gt_snps.sum(axis=1)

                    ## We only use the genotype 0/0, 0/1, 1/1:

                    if {'0/0', '0/1', '1/1'}.issubset(counts_gt_snps.columns):
                        counts_gt_snps['p0'] = ((2 * counts_gt_snps['0/0']) + counts_gt_snps['0/1']) / (2 * counts_gt_snps['total'])
                        counts_gt_snps['p1'] = ((2 * counts_gt_snps['1/1']) + counts_gt_snps['0/1']) / (2 * counts_gt_snps['total'])
                        counts_gt_snps['p0/0exp'] = counts_gt_snps['p0']**2
                        counts_gt_snps['p0/1exp'] = 2*counts_gt_snps['p0']*counts_gt_snps['p1']
                        counts_gt_snps['p1/1exp'] = counts_gt_snps['p1']**2
                        counts_gt_snps['p0/0obs'] = counts_gt_snps['0/0']/counts_gt_snps['total']
                        counts_gt_snps['p0/1obs'] = counts_gt_snps['0/1']/counts_gt_snps['total']
                        counts_gt_snps['p1/1obs'] = counts_gt_snps['1/1']/counts_gt_snps['total']

                        ## Chi-Square of Hardy-Weinberg test:
                        name_col1 = 'HWE_chisq_' + pop
                        hwe[name_col1] = np.where(((counts_gt_snps['p0/0exp'] != 0) & (counts_gt_snps['p0/1exp'] != 0) & (counts_gt_snps['p1/1exp'] != 0) & (counts_gt_snps['p0/0obs'] != 0) & (counts_gt_snps['p0/1obs'] != 0) & (counts_gt_snps['p1/1obs'] != 0)), chisquare([counts_gt_snps['p0/0exp'], counts_gt_snps['p0/1exp'], counts_gt_snps['p1/1exp']], f_exp=[counts_gt_snps['p0/0obs'], counts_gt_snps['p0/1obs'], counts_gt_snps['p1/1obs']])[0], np.NaN)

                        ## p-Value of Hardy-Weinberg test:
                        name_col2 = 'HWE_pval_' + pop
                        hwe[name_col2] = np.where(((counts_gt_snps['p0/0exp'] != 0) & (counts_gt_snps['p0/1exp'] != 0) & (counts_gt_snps['p1/1exp'] != 0) & (counts_gt_snps['p0/0obs'] != 0) & (counts_gt_snps['p0/1obs'] != 0) & (counts_gt_snps['p1/1obs'] != 0)), chisquare([counts_gt_snps['p0/0exp'], counts_gt_snps['p0/1exp'], counts_gt_snps['p1/1exp']], f_exp=[counts_gt_snps['p0/0obs'], counts_gt_snps['p0/1obs'], counts_gt_snps['p1/1obs']])[1], np.NaN)
                        
                        if hwe[name_col2].mean() < 0.05:
                            if warning_hw == "":
                                warning_hw = "Population/s " + pop
                            else:
                                warning_hw += ", " + pop

                        del counts_gt_snps

                if warning_hw != "":
                    warning_hw += " is/are not in HW equilibrium"

                if hwe.dropna().empty == False:
                    fig = draw_hwe(plot, hwe, pops)
                    if fig:
                        fig['layout'].update(title = 'HWE p-value by population', xaxis=dict(title='p-Value'))
                        return html.Div([dcc.Graph(figure=fig), html.Div(warning_hw,style={'color': 'red','padding-left':'3%','padding-bottom':'2%'})])
                    else:
                        layout = dict(title = 'HWE p-value by population', xaxis=dict(title='p-Value'))
                        return html.Div([dcc.Graph(figure=dict(layout=layout)), html.Div('It is not possible to perform the HWE analysis',style={'color': 'red','padding-left':'3%','padding-bottom':'2%'})])
                else:
                    layout = dict(title = 'HWE p-value by population', xaxis=dict(title='p-Value'))
                    return html.Div([dcc.Graph(figure=dict(layout=layout)), html.Div('It is not possible to perform the HWE analysis',style={'color': 'red','padding-left':'3%','padding-bottom':'2%'})])


        elif plot == "Inbreeding coefficient":

            ## Plot of the Inbreeding coefficient distribution (from the Hardy-Weinberg test, F = 1 - (n(0/1)observed / n(0/1)expected) )

            if popfile is not None:

                inbreed = pd.DataFrame()

                indiv_table = table_indiv(popfile,samples)
                pops = indiv_table['Population'].unique()

                for pop in pops:
                    if samples[0] + '_GT' in df.columns:
                        counts_gt_snps = df[indiv_table[indiv_table.Population == pop].Individual + '_GT'].replace('./.', np.NaN).apply(pd.value_counts,axis=1).fillna(0)
                    else:
                        counts_gt_snps = df[indiv_table[indiv_table.Population == pop].Individual].replace('./.', np.NaN).apply(pd.value_counts,axis=1).fillna(0)
                    counts_gt_snps['total'] = counts_gt_snps.sum(axis=1)
                    if {'0/0', '0/1', '1/1'}.issubset(counts_gt_snps.columns):
                        counts_gt_snps['p0'] = ((2 * counts_gt_snps['0/0']) + counts_gt_snps['0/1']) / (2 * counts_gt_snps['total'])
                        counts_gt_snps['p1'] = ((2 * counts_gt_snps['1/1']) + counts_gt_snps['0/1']) / (2 * counts_gt_snps['total'])
                        counts_gt_snps['p0/1exp'] = 2*counts_gt_snps['p0']*counts_gt_snps['p1']

                        inbreed[pop] = np.where((counts_gt_snps['p0/1exp'] != 0), 1-((counts_gt_snps['0/1']/counts_gt_snps['total'])/counts_gt_snps['p0/1exp']), np.NaN)

                        del counts_gt_snps

                layout = dict(title = 'Inbreeding coefficient by population')
                fig = draw_inbreed(plot, inbreed, pops)
                if fig:
                    fig['layout'].update(title = 'Inbreeding coefficient by population', xaxis=dict(title='Inbreeding coefficient'))
                    return html.Div([dcc.Graph(figure=fig)])
                else:
                    layout = dict(title = 'Inbreeding coefficient by population', xaxis=dict(title='Inbreeding coefficient'))
                    return html.Div([dcc.Graph(figure=dict(layout=layout)), html.Div('It is not possible to obtain the Inbreeding coefficient',style={'color': 'red','padding-left':'3%','padding-bottom':'2%'})])

# start Flask server
if __name__ == '__main__':
    app.run_server(
        threaded=True,
        debug=True,
        host='0.0.0.0',
        port=8050
    )

