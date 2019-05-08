
# VCFcheck
A tool for VCF files diagnostics
&nbsp;

&nbsp;
  
### Main features

- Input VCF file can be compressed (.gz) and uncompressed
- Tool adapted to work with multi-individual VCF and with VCF of only one individual
- Positions of the VCF can be SNPs, INDELS or ROHs
- Can filter the VCF by sample depth, mapping quality and missing data by SNP
- Can generate distribution plots of the missing data by SNP, the reference allele frequency and the sample depth by individual
- Can generate distribution plot of missing data by population, perform a PCA, test the Hardy-Weinberg equilibrium and generate the distribution plot of the inbreeding coefficient, if the VCF is multi-individual and the list of individuals and its populations is uploaded
- It allows to download the filtered VCF and the plots
&nbsp;

&nbsp;
  
### Requirements and execution

#### Pre-install:

1. Create virtual environment:

```conda create --name vcfcheck-app python=3.7```
  
2. Activate virtual environment:

```source vcfcheck-app/bin/activate```
  or
```conda activate vcfcheck-app```

3. Install dependencies:

```pip install -r requirements.txt```
  
#### Execute the program
  
1. Activate virtual environment:

```source vcfcheck-app/bin/activate```
  or
```conda activate vcfcheck-app```

2. Execute program:

```python VCFcheck.py```
  
3. Copy url on the web browser (http://0.0.0.0:8050/)
  
#### Exit from the virtual environment:
  
```deactivate```
&nbsp;

&nbsp;
  
### Usage
  
As input, the software requires a (compressed or uncompressed) VCF file:

```
##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##contig=<ID=1,length=315321322>
##contig=<ID=2,length=162569375>
##ALT=<ID=X,Description="Represents allele(s) other than observed.">
##INFO=<ID=INDEL,Number=0,Type=Flag,Description="Indicates that the variant is an INDEL.">
##INFO=<ID=IDV,Number=1,Type=Integer,Description="Maximum number of reads supporting an indel">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Raw read depth">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="List of Phred-scaled genotype likelihoods">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Number of high-quality bases">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##INFO=<ID=MQ,Number=1,Type=Integer,Description="Average mapping quality">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">
##FILTER=<ID=LOWQUAL,Description="Set if true: %QUAL<10 || %MAX(DP)<5 || %MAX(DP)>15">
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  IBGU1330
1       4040059 .       T       .       .       PASS    END=4040089     GT      0/0
1       4624048 .       T       .       .       PASS    END=4624079     GT      0/0
1       7900743 .       G       .       .       PASS    END=7900743     GT      0/0
1       14739858        .       A       .       .       PASS    END=14739886    GT      0/0
1       15676648        .       T       .       .       PASS    END=15676649    GT      0/0
1       17520198        .       G       .       .       PASS    END=17520206    GT      0/0
1       19562975        .       A       .       .       PASS    END=19563008    GT      0/0
1       21407709        .       A       T       118     PASS    DP=7;VDB=0.788273;SGB=-0.590765;MQSB=0.5;MQ0F=0;AC=2;AN=2;DP4=0,0,2,3;MQ=53     GT:PL:DP        1/1:145,15,0:5
1       22086814        .       A       G       177     PASS    DP=7;VDB=0.86172;SGB=-0.636426;MQSB=1;MQ0F=0;AC=2;AN=2;DP4=0,0,5,2;MQ=60        GT:PL:DP        1/1:204,21,0:7
```
  
As optional input, it is possible to upload a file with all the samples in a column and another column with the populations to which each sample belongs (separated by tab).

```
SAMPLE1	POP1
SAMPLE2	POP1
SAMPLE3	POP1
SAMPLE4	POP1
SAMPLE5	POP2
SAMPLE6	POP2
SAMPLE7	POP2
SAMPLE8	POP3
SAMPLE9	POP3
SAMPLE10	POP3
SAMPLE11	POP3
```

When the inputs are uploaded, the type of mutations to analyze must be selected before starting the analysis. 

![img1b](https://user-images.githubusercontent.com/30473077/57376107-81077880-719f-11e9-8ab3-7b3675d4c705.png)
&nbsp;

&nbsp;
  
Then, the table from the VCF and a summary will be displayed. In that moment, it is possible to filter the file by sample depth, mapping quality or missing data using the corresponding sliders.

![img2](https://user-images.githubusercontent.com/30473077/57376314-fd9a5700-719f-11e9-91b1-78f81253aaf6.png)
&nbsp;

&nbsp;
  
Once the table is generated, we can perform different plots and analysis by its selection in the dropdown and download the filtered VCF or the summary using the corresponding buttons.

![img3](https://user-images.githubusercontent.com/30473077/57376321-fffcb100-719f-11e9-839b-42f8abd79580.png)
