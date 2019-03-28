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

    lines = [l for l in f if not l.startswith('##')]

    for c in lines[0].split('\t'):
        c = c.rstrip('\n')
        if c == 'POS':
            dict_dtypes[c] = 'int'
        elif c == 'INFO':
            dict_dtypes[c] = 'str'
        else:
            dict_dtypes[c] = 'category'

    return pd.read_csv(
        io.StringIO(''.join(lines)),
        sep='\t',
        dtype=dict_dtypes
        ).rename(columns={'#CHROM': 'CHROM'})


#contig, info_fields, genotype_fields, filters_fields = read_header(vcfFile)
#vcf_table = read_vcf(vcfFile)

def obtain_samples(table):

    '''Returns the list of samples'''

    samples = pd.Series(table.columns[9:])
    
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
    missing_individual = vcf_samples.isnull().sum()/(vcf_samples.isnull().sum()+vcf_samples.notnull().sum())
    
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

    missing_snp = vcf_samples.T.isnull().sum()/(vcf_samples.T.isnull().sum()+vcf_samples.T.notnull().sum())

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
        if table[i].notnull().sum() == 0:
            table.drop(i, axis=1, inplace=True)

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
            df_genotype.fillna(value=0, inplace=True)
            table[colGT] = df_genotype[0].astype('category')
            #table[colPL] = df_genotype[1]
            table[colDP] = df_genotype[2].astype(int)
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

def dropdown_options():

    '''Possible plots to perform'''

    plots = ["Missing by Population", 
             "Missing by SNP", 
             "Reference allele frequency",
             "Depth by individual",
             "Depth by genotype",
             "Principal Component Analyisis (PCA)",
             "Hardy-Weinberg test",
             "Inbreeding coefficient"]

    plots_options = (
        [{'label': p, 'value': p}
         for p in plots]
    )
    return plots_options

## Set up Dashboard and create layout:
app = dash.Dash()
app.css.append_css({
    "external_url": "https://codepen.io/chriddyp/pen/bWLwgP.css"
})
#cache = Cache(app.server, config={
#    'CACHE_TYPE': 'filesystem',
#    'CACHE_DIR': 'cache-directory'
#})

## Layout of the web application:
app.layout = html.Div([

    ## Page Header
    html.Div([
        html.H2(
            'VCFcheck: A tool for VCF files diagnostics',
            id='title'
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
                'width': '15%',
                'height': '60px',
                'lineHeight': '60px',
                'borderWidth': '1px',
                'borderStyle': 'dashed',
                'borderRadius': '5px',
                'textAlign': 'center',
                'margin': '10px'
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
                'width': '15%',
                'height': '60px',
                'lineHeight': '60px',
                'borderWidth': '1px',
                'borderStyle': 'dashed',
                'borderRadius': '5px',
                'textAlign': 'center',
                'margin': '10px'
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
            className='three columns',
        ),

        ## Select the range of Depth by Individual:
        html.Div([
            dcc.RangeSlider(
                id='sliderDP',
                min=0,
                max=100,
                step=1,
                value=[0,100],
                updatemode='drag',
                className='two columns',
            ),
            html.Div(id='output-container-range-sliderDP', className='three columns')
        ], style={'padding': 20}),

        ## Select the range of Map quality by SNP:
        html.Div([
            dcc.RangeSlider(
                id='sliderMQ',
                min=0,
                max=100,
                step=1,
                value=[0,100],
                updatemode='drag',
                className='two columns',
            ),
            html.Div(id='output-container-range-sliderMQ', className='three columns')
        ], style={'padding': 20}),

        ## Select the range of Missing by SNP:
        html.Div([
            dcc.RangeSlider(
                id='sliderMiss',
                min=0,
                max=1,
                step=0.01,
                value=[0,1],
                updatemode='drag',
                className='two columns',
            ),
            html.Div(id='output-container-range-sliderMiss', className='three columns')
        ], style={'padding': 20}),

    ],id='div1',className='row'),
    
    ## Objects to go between callbacks:
    html.Div(id='output-vcf', style={'display': 'none'}),
    html.Div(id='output-vcf-filt', style={'display': 'none'}),
    html.Div(id='output-samples', style={'display': 'none'}),
    html.Div(id='header', style={'display': 'none'}),

    ## Name of inputs files:
    html.Div(id='output-name'),
    html.Div(id='output-popname'),

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
                    html.Div(dcc.Dropdown(id='plot-selector', options=dropdown_options()), className='four columns')
                ]),
            ], className='six columns'),
        ], style={'padding': 50}),

        ## Plot:
        html.Div(id='plot-graph', className='six columns'),

        ## Button for download the filtered VCF:
        html.Div([
            html.Button('Download VCF', id='button-download'),
            html.Div(id='output-download'),
        ], style={'padding': 20}),

    ],id='div2',className='row'),

])

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
        html.Hr(),  # horizontal line
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

    if (samples[0] + '_DP') in table:    
        figure = ff.create_distplot([table[s + '_DP'].dropna() for s in samples], samples, show_hist=False)
        return figure

def draw_depth_gt(plot,table,samples):

    '''Representation of the distribution of the depth by genotype'''

    if (samples[0] + '_DP') in table:    

        gtposib = ['0/0','0/1','1/1','0/2','1/2','2/2','0/3','1/3','2/3','3/3']
    
        list_gt = {}
    
        for s in samples:
            if s + '_GT' in table.columns:
                for gt in gtposib:
                    gtlistdp = table[table[s + '_GT'] == gt][s + '_DP'].dropna().values
                    if len(gtlistdp) != 0:
                        if gt in list_gt:
                            list_gt[gt] += gtlistdp
                        else:
                            list_gt[gt] = gtlistdp
    
        if list_gt != {}:
            figure = ff.create_distplot([list_gt[key] for key in list_gt], [key for key in list_gt.keys()], show_hist=False)
            return figure

def draw_pca(plot,table):

    '''Representation of the PCA'''

    pops = table['Population'].unique()
    data = []
    colors = ['red', 'green', 'blue']
    
    for name, col in zip(pops, colors):
        indicesToKeep = table['Population'] == name
        data.append(go.Scatter(x = table.loc[indicesToKeep, 'principal component 1'],
                               y = table.loc[indicesToKeep, 'principal component 2'],
                               name=name,
                               mode = 'markers',
                               marker= dict(size = 10,
                                            color = col,
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

## Upload VCF file:
@app.callback([Output('output-name', 'children'),
              Output('output-vcf', 'children'),
              Output('output-samples', 'children'),
              Output('header', 'children'),
              Output('sliderDP', 'min'),
              Output('sliderDP', 'max'),
              Output('sliderDP', 'disabled'),
              Output('sliderMQ', 'min'),
              Output('sliderMQ', 'max'),
              Output('sliderMQ', 'disabled'),
              Output('sliderMiss', 'min'),
              Output('sliderMiss', 'max')],
              [Input('radio-typeSNP', 'value')],
              [State('upload-vcf', 'filename')])
def update_nameoutput(type_snps,filename):

    if (filename is not None) and (type_snps is not None):

        starttime = time.time()

        ## Read the input VCF file (mutations and header separately)
        df = read_vcf(filename)
        contig, info_fields, genotype_fields, filters_fields, header = read_header(filename)

        endtime = time.time()
        print('1: ' + str(endtime-starttime))

        ## Obtain a list of the VCF samples
        samples = obtain_samples(df)

        ## Obtain the possible type of positions from the header (SNPs, INDELs or ROHs)
        type_pos(df, info_fields)

        ## Filter by the selected type of position
        df_filt = {}
        if type_snps == 'all':
            df_filt = df
            df_filt.drop('END', axis=1, inplace=True)
            df_filt.drop('INDEL', axis=1, inplace=True)
        elif type_snps == 'snpindel':
            df_filt = df
            if not df_filt.empty:
                add_genotype_cols(df_filt, genotype_fields, samples)
                add_info_cols(df_filt, info_fields)
            df_filt.drop('END', axis=1, inplace=True)
            df_filt.drop('INDEL', axis=1, inplace=True)

        elif type_snps == 'snps':
            df_filt = table_snps(df)
            if not df_filt.empty:
                add_genotype_cols(df_filt, genotype_fields, samples)
                add_info_cols(df_filt, info_fields)
        elif type_snps == 'biallelic':
            df_filt = biallelic_snps(df)
            if not df_filt.empty:
                add_genotype_cols(df_filt, genotype_fields, samples)
                add_info_cols(df_filt, info_fields)
        elif type_snps == 'multiallelic':
            df_filt = multiallelic_snps(df)
            if not df_filt.empty:
                add_genotype_cols(df_filt, genotype_fields, samples)
                add_info_cols(df_filt, info_fields)
        elif type_snps == 'indels':
            df_filt = table_indels(df)
            if not df_filt.empty:
                add_genotype_cols(df_filt, genotype_fields, samples)
                add_info_cols(df_filt, info_fields)
                df_filt.drop('INDEL', axis=1, inplace=True)
        elif type_snps == 'rohs':
            df_filt = table_rohs(df)
            if not df_filt.empty:
                add_genotype_cols(df_filt, genotype_fields, samples)
                add_info_cols(df_filt, info_fields)

        del df

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

        if 'MQ' in df_filt.columns:
            minMQ = df_filt.MQ.min()
            maxMQ = df_filt.MQ.max()
            disableMQ = False
        else:
            minMQ = 0
            maxMQ = 0
            disableMQ = True

        df_filt['Missing'] = missing_snp(df_filt, samples)
        minMiss = df_filt.Missing.min()
        maxMiss = df_filt.Missing.max()

        ## Convert the Pandas DataFrame into Json table to the transport of the table to other callback
        df_json = df_filt.to_json(date_format='iso', orient='table')

        #feather.write_dataframe(df_filt, 'df.feather')

        return 'Input VCF file: ' + filename, df_json, samples, header, minDP_f, maxDP_f, disableDP, minMQ, maxMQ, disableMQ, minMiss, maxMiss

## Upload Poplist name:
@app.callback(Output('output-popname', 'children'),
              [Input('upload-poplist', 'filename')])
def update_poplist(filename):
    if filename is not None:
        return 'Input list of individuals and breeds: ' + filename

## Display the datatable of the filtered VCF
@app.callback([Output('datatable-upload', 'children'),
              Output('output-vcf-filt', 'children')],
              [Input('output-vcf', 'children'),
              Input('output-samples', 'children'),
              Input('sliderDP', 'value'),
              Input('sliderMQ', 'value'),
              Input('sliderMiss', 'value')])
def showing_table(jsonified_dataframe,samples,valueDP,valueMQ,valueMiss):

    if jsonified_dataframe:

        #df = feather.read_dataframe('df.feather')

        ## Read the Json table in Pandas DataFrame format:
        df = pd.read_json(jsonified_dataframe, orient='table')

        ## Apply the selected ranges of the sliders:
        for s in samples:
            colnameDP = s + '_DP'
            colnameGT = s + '_GT'
            if colnameDP in df.columns:
                df[colnameGT] = np.where((df[colnameDP] >= valueDP[0]) & (df[colnameDP] <= valueDP[1]), df[colnameGT], './.')

        if 'MQ' in df.columns:
            df = df.loc[(df['MQ'] >= valueMQ[0]) & (df['MQ'] <= valueMQ[1])]

        df = df.loc[(df['Missing'] >= valueMiss[0]) & (df['Missing'] <= valueMiss[1])]

        ## Convert the filtered DataFrame into Json table:
        df_json = df.to_json(date_format='iso', orient='table')

        return generate_table(df), df_json

## Download the filtered VCF:
@app.callback(Output('output-download', 'children'),
             [Input('button-download', 'n_clicks')],
             [State('output-vcf-filt', 'children'),
             State('output-samples', 'children'),
             State('header', 'children'),
             State('output-name','children')])
def update_download_link(nclicks,jsonified_dataframe,samples,header,filename):

    if jsonified_dataframe and (nclicks > 0):
        
        ## For the output name we use the input name + ".filtered.vcf":
        outputfile = filename.split(' ')[-1]
        f = open(outputfile + '.filtered.vcf', 'w')

        ## Write the header in the output VCF:
        f.write('## Filtered by VCFcheck' + '\n')
        for h in header:
            f.write(h + '\n')

        ## Read the Json table in Pandas DataFrame format:
        df = pd.read_json(jsonified_dataframe, orient='table')

        ## Select the standard columns of a VCF to write them in the output VCF:
        columns_str = [col for col in ['CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT'] + samples]
        df = df[columns_str].rename(columns={'CHROM': '#CHROM'})
        csv_string = df.to_csv(sep='\t', index=False, encoding='utf-8',na_rep='NA')
        f.write(csv_string)
        #csv_string = "data:text/csv;charset=utf-8," + urllib.parse.quote(bytes(csv_string))

        return outputfile + '.filtered.vcf' + ' downloaded!'

## Print value of the selected range of DP next to the Range Slider:
@app.callback(Output('output-container-range-sliderDP', 'children'),
              [Input('sliderDP', 'value')])
def update_output(value):
    range_selected = 'Sample depth (DP) selected between ' + str(value[0]) + ' and ' + str(value[1])
    return range_selected

## Print value of the selected range of MQ next to the Range Slider:
@app.callback(Output('output-container-range-sliderMQ', 'children'),
              [Input('sliderMQ', 'value')])
def update_output(value):
    range_selected = 'Average mapping quality (MQ) selected between ' + str(value[0]) + ' and ' + str(value[1])
    return range_selected

## Print value of the selected range of Missing next to the Range Slider:
@app.callback(Output('output-container-range-sliderMiss', 'children'),
              [Input('sliderMiss', 'value')])
def update_output(value):
    range_selected = 'Percentage selected of missing data between ' + str(value[0]) + ' and ' + str(value[1])
    return range_selected

## Display the selected plot in the dropdown
@app.callback(Output('plot-graph', 'children'),
              [Input('plot-selector', 'value'),
              Input('output-vcf-filt', 'children'),
              Input('sliderDP', 'value'),
              Input('sliderMQ', 'value'),
              Input('sliderMiss', 'value')],
              [State('output-samples', 'children'),
              State('upload-poplist', 'filename')])
def load_distribution_graph(plot, jsonified_dataframe, valueDP, valueMQ, valueMiss, samples, popfile):

    if plot is not None:

        ## Read the Json table in Pandas DataFrame format:
        df = pd.read_json(jsonified_dataframe, orient='table')
        #df = feather.read_dataframe('df.feather')

        if plot == "Missing by SNP":

            ## Density plot of the distribution of missing data by SNP

            fig = draw_missing_snp(plot, df['Missing'])
            layout = dict(title = 'Missing by SNP')     
            return html.Div([dcc.Graph(figure=dict(data=fig, layout=layout))])

        elif plot == "Missing by Population":

            ## Density plot of the distribution of missing data by population (it is necessary the input file of individuals-populations)

            if popfile is not None:
                missing_individual = missing_ind(df, samples)
    
                indiv_table = table_indiv(popfile,samples)
                indiv_table.index = indiv_table.Individual
                indiv_table = pd.concat([indiv_table, missing_individual], axis=1, sort=False, join='inner')
                indiv_table.rename(columns = {0:"Missing"}, inplace=True)
        
                tablepop = indiv_table.replace('./.', np.NaN)
    
                if len(tablepop) > 0:
                    fig = draw_missing_pop(plot, tablepop)

                layout = dict(title = 'Missing by Population')

                return html.Div([dcc.Graph(figure=dict(data=fig, layout=layout))])

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
            layout = dict(title = 'Reference allele frequency')

            del reffreq

            return html.Div([dcc.Graph(figure=dict(data=fig, layout=layout))])

        elif plot == "Depth by individual":

            ## Density plot of the depth distribution by individual

            fig = draw_depth_ind(plot, df, samples)
            layout = dict(title = 'Depth by individual')

            return html.Div([dcc.Graph(figure=dict(data=fig, layout=layout))])

        elif plot == "Depth by genotype":

            ## Density plot of the depth distribution by genotype

            fig = draw_depth_gt(plot, df, samples)
            layout = dict(title = 'Depth by genotype')

            return html.Div([dcc.Graph(figure=dict(data=fig, layout=layout))])

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
                        vcf_samples[s_gt] = df[s_gt].astype('str').replace({ "\.\/\." : np.NaN, "0\/0" : 1, "0\/." : 0.5, ".\/." : 0}, regex=True)
                    else:
                        vcf_samples[s] = df[s].astype('str').replace({ "\.\/\." : np.NaN, "0\/0" : 1, "0\/." : 0.5, ".\/." : 0}, regex=True)

                vcf_samples.dropna(how = 'all',inplace=True)

                freq_nomiss = vcf_samples.T.sum()/vcf_samples.T.notnull().sum()
                nnull = vcf_samples.T.isnull().sum()
                weight_miss = (freq_nomiss + (freq_nomiss*nnull))/len(samples)
    
                vcf_samples.fillna(value = "NA", inplace = True)

                for s in samples:
                    s_gt = s + '_GT'
    
                    if s_gt in df.columns:
                        vcf_samples[s_gt] = np.where(vcf_samples[s_gt] == "NA", weight_miss, vcf_samples[s_gt])
                    else:
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

                return html.Div([dcc.Graph(figure=dict(data=[n for n in data], layout=layout))])

        elif plot == "Hardy-Weinberg test":

            ## Plot of the Hardy Weinberg p-Value distribution

            if popfile is not None:

                hwe = pd.DataFrame()

                indiv_table = table_indiv(popfile,samples)
                pops = indiv_table['Population'].unique()

                for pop in pops:
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

                        del counts_gt_snps

                fig = draw_hwe(plot, hwe, pops)
                layout = dict(title = 'HWE p-value by population')

                return html.Div([dcc.Graph(figure=dict(data=fig, layout=layout))])

        elif plot == "Inbreeding coefficient":

            ## Plot of the Inbreeding coefficient distribution (from the Hardy-Weinberg test, F = 1 - (n(0/1)observed / n(0/1)expected) )

            if popfile is not None:

                inbreed = pd.DataFrame()

                indiv_table = table_indiv(popfile,samples)
                pops = indiv_table['Population'].unique()

                for pop in pops:
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

                return html.Div([dcc.Graph(figure=dict(data=fig, layout=layout))])

# start Flask server
if __name__ == '__main__':
    app.run_server(
        threaded=True,
        debug=True,
        host='0.0.0.0',
        port=8050
    )

