import sys,os,yaml,warnings,logging,functools,pickle,configparser
from datetime import datetime
import pandas as pd
import scanpy as sc
import numpy as np
#import paste as pst
#import rpy2.robjects as robjects
#from rpy2.robjects import pandas2ri
#readRDS = robjects.r['readRDS']
#from ipython_exit import exit

#print=functools.partial(print, flush=True)
warnings.simplefilter(action='ignore', category=FutureWarning)
logging.disable(level=logging.INFO)
batchKey="library_id"
strPipePath=os.path.dirname(os.path.realpath(__file__))

def msgError(msg):
    print(msg)
    exit()
def MsgInit():
    print("\n\n*****",datetime.now().strftime("%Y-%m-%d %H:%M:%S"),"*****")
    appPath=os.path.dirname(strPipePath)
    if os.path.isdir(os.path.join(appPath,".git")):
        gitConfig = configparser.ConfigParser()
        tmp = gitConfig.read(os.path.join(appPath,".git","config"))
        url = gitConfig['remote "origin"']['url']
    
        gitLog = pd.read_csv(os.path.join(appPath,".git","logs","HEAD"),sep="\t",header=None)
        gitLog = gitLog.iloc[-1,0].split(" ")
        print("###########\n## scRNAsequest: %s"%url)
        print("## Pipeline Path: %s"%appPath)
        print("## Pipeline Date: %s %s"%(datetime.fromtimestamp(int(gitLog[-2])).strftime('%Y-%m-%d %H:%M:%S'),gitLog[-1]))
        print("## git HEAD: %s\n###########\n"%gitLog[1])
    a=getSys()#testing the exists
def getConfig(strConfig):
    print("Initializing ...")
    if not os.path.isfile(strConfig):
        msgError("Missing config: %s"%strConfig)    
    with open(strConfig,"r") as f:
        config = yaml.safe_load(f)
    if config.get("sample_meta") is not None and os.path.isfile(config["sample_meta"]):
        sInfo = pd.read_csv(config["sample_meta"])
    else:
        msgError("Missing sample meta: %s"%config.get("sample_meta"))
    if config.get('output') is None:
        msgError("The output path is required in config")
    if config.get('prj_name') is None:
        msgError("The project name (prj_name) is required in config")
    if config.get('clusterN') is None:
        msgError("The cluster numbers (clusterN) are required in config")
    clusterN = [int(a) for a in config.get('clusterN') if isinstance(a,int) or isinstance(a,float) or a.isdigit()]
    if len(clusterN)==0:
        msgError("The integer list of cluster numbers (clusterN) are required in config")
    config['clusterN'] = clusterN
    return config, sInfo
def readLines(strF):
    with open(strF,"r") as f:
        lines = f.readlines()
    return lines
def getSys():
    strSys=os.path.join(strPipePath,"sys.yml")
    if not os.path.isfile(strSys):
        msgError("missing sys.yml: Please contact admin to setup system config file!")
    with open(strSys,"r") as f:
        config = yaml.safe_load(f)
    return config

def readPkl(strPkl):
    with open(strPkl,'rb') as f:
        obj=pickle.load(f)
    return obj
def writePkl(obj,strPkl):  
    with open(strPkl, 'wb') as fp:
        pickle.dump(obj, fp) 

def fillNA(x,v):#for panda series
    if x.isna().any():
        if x.dtype=='category' and not x.cat.categories.isin([v]).any():
            x=x.cat.add_categories(v)
        x.fillna(v,inplace=True)
    return x
## visium SpaceRanger (sr) functions
sr_path_column="SpaceRanger_path"
def sr_checkMeta(meta,config):
    sName = config.get("sample_name")
    if sName is None:
        msgError("sample_name is required in the config file!")
    if not sName in meta.columns:
        msgError("%s specified in config is not a column in the sample sheet!"%sName)
    if not sr_path_column in meta.columns:
        msgError("%s is required column in the sample sheet!"%sr_path_column)
def sr_extract(strH5ad,sName):
    D = sc.read_h5ad(strH5ad)
    adata = D[D.obs[batchKey]==sName].copy()
    adata.uns['spatial'] = {sName:adata.uns['spatial'][sName]}
    return adata
def sr_merge(strPkl,strH5ad,config):
    if strH5ad is None or strPkl is None:
        return
    if not strH5ad is None and os.path.isfile(strH5ad):
        print("The raw merged file exists: %s\n\tIf a new merge is required, please rename/remove the above file!"%strH5ad)
        return
    with open(strPkl,"rb") as f:
        D_dict = pickle.load(f)
    D_list = list(D_dict.values())
    adata = sc.AnnData.concatenate(*D_list,
      join="outer",
      batch_categories=list(D_dict.keys()),
      batch_key=batchKey,
      uns_merge="unique")
    sr_filter(adata,config)
    adata.write(strH5ad)
def sr_addGP(adata,config):
    groupKey = 'gene_group'
    rmGene=np.full(adata.shape[1],False)
    varKey=[]
    if groupKey in config.keys():
        for k in config[groupKey]:
            if type(config[groupKey][k]['startwith']) is not list:
                msgError("config error with %s, 'startwith' has to be a list"%k)
            gList = np.full(adata.shape[1],False)
            for one in config[groupKey][k]['startwith']:
                if len(one)>1:
                    gList |= adata.var_names.str.startswith(one)
            if gList.sum()==0:
                print("\t\tNo genes found for %s"%k)
                continue
            print("\t\tGene group %s contains %d genes"%(k,gList.sum()))
            adata.var[k]=gList
            varKey += [k]
            if config[groupKey][k]['rm']:
                print("\t\t%s genes will be removed"%k)
                rmGene |= gList
    if len(varKey)>0:
        sc.pp.calculate_qc_metrics(adata,qc_vars=varKey,inplace=True)
        selObs = np.full(adata.shape[0],True)
        for one in varKey:
            adata.obs.drop("total_counts_%s"%one,axis=1,inplace=True)
            adata.obs.drop("log1p_total_counts_%s"%one,axis=1,inplace=True)
            selObs = np.logical_and(selObs,adata.obs["pct_counts_%s"%one]<config[groupKey][one]["cutoff"])
            print("\t Cell, %s cutoff %d: %d cells"%(one,config[groupKey][one]["cutoff"],selObs.sum()))
        adata._inplace_subset_obs(selObs)
    else:
        sc.pp.calculate_qc_metrics(adata,inplace=True)
    
    if rmGene.sum()>0:
        adata._inplace_subset_var(np.invert(rmGene))
        print("\t Total of %d genes are removed",rmGene.sum())
def sr_filter(D,config):
    print("*** filter ***")
    print("\t Total: %d cells with %d genes"%(D.shape))
    sc.pp.filter_genes(D,min_counts=config.get('UMI_count_gene'))
    print("\t Gene, min UMI %d: %d cells with %d genes"%((config.get('UMI_count_gene'),)+D.shape))
    sc.pp.filter_genes(D,min_cells=config.get('min_cells'))
    print("\t Gene, min cells %d: %d cells with %d genes"%((config.get('min_cells'),)+D.shape))
    sc.pp.filter_cells(D,min_counts=config.get('UMI_count_cell'))
    print("\t Cell, min UMI %d: %d cells with %d genes"%((config.get('UMI_count_cell'),)+D.shape))
    sc.pp.filter_cells(D,min_genes=config.get('min_features'))
    print("\t Cell, min gene %d: %d cells with %d genes"%((config.get('min_features'),)+D.shape))
    sc.pp.filter_cells(D,max_counts=config.get('highCount_cutoff'))
    print("\t Cell, max UMI %d: %d cells with %d genes"%((config.get('highCount_cutoff'),)+D.shape))
    sc.pp.filter_cells(D,max_genes=config.get('highGene_cutoff'))
    print("\t Cell, max gene %d: %d cells with %d genes"%((config.get('highGene_cutoff'),)+D.shape))
    sr_addGP(D,config)
    return

def sr_read(meta,config,strPkl=None,strH5ad=None):
    sName = config['sample_name']
    print("Reading samples ...")
    if not strPkl is None and os.path.isfile(strPkl):
        print("The raw file exists: %s\n\tIf a new read is required, please rename/remove the above file!"%strPkl)
        sr_merge(strPkl,strH5ad,config)
        return
    all_slices = {}
    for i in range(meta.shape[0]):
        print("\t",meta[sName][i])
        oneD = sc.read_visium(meta[sr_path_column][i])
        oneD.var_names_make_unique()
        sr_addGP(oneD,config)
        all_slices[meta[sName][i]]=oneD
    if strPkl is None:
        return all_slices
    with open(strPkl,"wb") as f:
        pickle.dump(all_slices,f)
    if not strH5ad is None and os.path.isfile(strH5ad):
        os.remove(strH5ad)
    sr_merge(strPkl,strH5ad,config)
    return

## preprocess
def filterSlice(slices,min_gene=None,max_gene=None,min_count_cell=100,min_cell=None,min_count_gene=15):
    for i in range(len(slices)):
        sc.pp.filter_genes(slices[i],min_counts=min_count_gene,min_cells=min_cell)
        sc.pp.filter_cells(slices[i],min_counts=min_count_cell,min_genes=min_gene,max_genes=max_gene)

def filterD(D,min_gene=None,max_gene=None,min_count_cell=100,min_cell=None,min_count_gene=15):
    sc.pp.filter_genes(D,min_counts=min_count_gene,min_cells=min_cell)
    sc.pp.filter_cells(D,min_counts=min_count_cell,min_genes=min_gene,max_genes=max_gene)
  
def preAlign(all_slices):
    pass

def main():
    if len(sys.argv)<2:
        print("Unknown utility call!")
        return()
    task = sys.argv[1]
    if task=="mergeRDS":
        mergeRDS(sys.argv[2],sys.argv[3])
    else:
        print("Unknown task: %s"%task)

if __name__ == "__main__":
    main()
