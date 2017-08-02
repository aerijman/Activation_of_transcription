'''import libraries'''
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import glob


''' Create the dictionary to match wiggle or other files '''
roma = {'chrI':'1','chrII':'2','chrIII':'3','chrIV':'4','chrV':'5','chrVI':'6','chrVII':'7',\
'chrVIII':'8','chrIX':'9','chrX':'10','chrXI':'11','chrXII':'12','chrXIII':'13','chrXIV':'14',\
'chrXV':'15','chrXVI':'16','chrXVII':'17','chr17':'17','chrM':'M'}


''' 
    Create DataFrame of SGD which already contains the TSSs from Park et al. 
'''
def SGD():
    global sgd
    #sgd = pd.read_csv('/home/aerijman/ChIPseq/scripts/sgd_tss.csv')
    #sgd.columns = np.insert(sgd.columns.values.tolist()[1:],0,'locus')
    #sgd.dropna(inplace=True)
    #sgd.set_index('locus', inplace=True)
    tss = pd.read_excel('./data/supp_gkt1366_nar-02868-n-2013-File009.xls', sheetname='TSS')
    tes = pd.read_excel('./data/supp_gkt1366_nar-02868-n-2013-File009.xls', sheetname='PAS')

    # IN order to join them I remove duplicated columns names
    tss.rename(columns={'ORF':'locus'}, inplace=True)
    tss.set_index('locus', inplace=True)
    tss.rename(columns={'coordinate':'tss', 'chr':'chromosome'}, inplace=True)
    tss = tss[['chromosome', 'tss']]
    tss.chromosome = [roma[i] for i in tss.chromosome]

    tes.rename(columns={'ORF':'locus'}, inplace=True)
    tes.set_index('locus', inplace=True)
    tes.drop('chr', axis=1, inplace=True)
    tes.rename(columns={'coordinate':'tes'}, inplace=True)

    # join tss and tes and add features to ease downstream processing
    tss_tes = pd.concat([tss,tes], join='inner', axis=1)
    tss_tes['strand'] = (tss_tes.tes-tss_tes.tss) / abs(tss_tes.tes-tss_tes.tss)
    tss_tes['tss-50'] = tss_tes.tss - 50*tss_tes.strand 
    
    tss_tes['half'] = tss_tes.tss + (abs(tss_tes.tss - tss_tes.tes)/2)*tss_tes.strand
    tss_tes['min_cassette'] = tss_tes[['half','tss-50']].min(axis=1)
    tss_tes['max_cassette'] = tss_tes[['half','tss-50']].max(axis=1)
    tss_tes['min'], tss_tes['max'] = tss_tes[['tss','tes']].min(axis=1), tss_tes[['tss','tes']].max(axis=1)

    sgd = tss_tes.sort_values(['chromosome','min'])
    
    return sgd

''' 
    Call the f(x) in case the output is needed for downstream f(x)s
'''
SGD()


''' 
    wiggle class includes the wig dataframe (wiggle file) 
    as object and the chrom array (with pointers to the 
    chromosome positions within the wig file).
    It can retrieve a dataframe with absolute cumulative values
    of reads for each gene.
'''
class wiggle:
    
    def __init__(self, f):  
        
        ''' 
            Read wiggle file into a dataframe.This f(x) takes a wiggle file as argument
        '''
        def read_wig(archivo):
            fg, chrom, m = 0, [0 for i in xrange(17)], []
            print archivo
            for i in open(archivo):
                fg+=1
                if fg<2:
                    continue
                if i[:3]=='var':
                    chrom[int(roma[i.strip()[19:]])]=fg-2
                    fg-=1
                    continue
                posRead = np.array(i.strip().split("\t")).astype(float)
                m.append(posRead)
                
            return pd.DataFrame(m, columns=['position','reads']), chrom
        
        self.wig, self.chrom = read_wig(f)
        self.sample=f

    
    
        ''' 
            f(x) retrive array [upstreamOfTSS...downstreamOfTSS] for a locus    
            that is user-defined  
        '''        
    
    def chromWig(wig, locus, upstreamOfTSS, downstreamOfTSS):
        [chromo, tss, strand] = sgd.loc[locus,['chromosome','tss','strand']].values.astype(float)
        tes, chromo, tss = tss+(downstreamOfTSS * strand), int(chromo), tss-(upstreamOfTSS * strand)
        pointer_i, pointer_f = wig.chrom[chromo], [wig.chrom[chromo+1] if chromo<16 else len(wig.wig)-1][0]
        [tmin, tmax] = [[tss,tes] if tes>tss else [tes,tss]][0]
        tmp = wig.wig.iloc[pointer_i:pointer_f]
        if tes>tss:
            return tmp[(tmp.position<tmax) & (tmp.position>tmin)]
        else:
            return tmp[(tmp.position<tmax) & (tmp.position>tmin)].iloc[::-1]
        
        
    ''' 
        f(x) merges SGD and TSS and returns a single value (cumulative reads) 
        for each gene (from wiggle.reads) in a dataframe 
        This f(x) takes 4 arguments: A wiggle object, complete sgd or a slice of it 
        (including only selected genes), #bp down and upstream 
    '''
    def wig2log(self, sgd_slice=sgd, upstreamOfTSS=1, downstreamOfTSS=100):
        
        reads = []
        for i in sgd_slice.itertuples():
            tmp = self.chromWig(i.Index, upstreamOfTSS, downstreamOfTSS)
            reads.append([i.Index,tmp.reads.sum()])
        df = pd.DataFrame(reads, columns=['locus','reads'])
        df.set_index('locus', inplace=True)
            
        return df
    
    
    ''' f(x) merges SGD and TSS and returns a cummulative value of reads\
        for ALL genes (from wiggle.reads) in a matrix [upstreamOfTSS...downstreanOfTSS]  
        This f(x) takes 4 arguments: A wiggle object, complete sgd or a slice of it 
        (including only selected genes), #bp down and upstream 
    '''    
    def meta_genomic(self, sgd_slice=sgd, upstreamOfTSS=250, downstreamOfTSS=750):
        
        mnew, n = pd.DataFrame(np.array(np.arange(-upstreamOfTSS,downstreamOfTSS)).T), 0
        mnew.set_index(0, inplace=True)

        for i in sgd.itertuples():
            tmp = self.chromWig(i.Index, upstreamOfTSS, downstreamOfTSS)    
            tmp.columns=['position', i.Index]

            tmp['zero_position'] = (tmp.position - i.tss) * i.strand 
            tmp.set_index('zero_position', inplace=True)        
            mnew = pd.concat([mnew,tmp[i.Index]], axis=1)
    
            n+=1

        return mnew