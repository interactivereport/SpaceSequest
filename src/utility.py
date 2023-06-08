import os,yaml,warnings,logging,functools,pickle
import pandas as pd
import scanpy as sc
import paste as pst
#from ipython_exit import exit

#print=functools.partial(print, flush=True)
warnings.simplefilter(action='ignore', category=FutureWarning)
logging.disable(level=logging.INFO)
batchKey="library_id"
strPipePath=os.path.dirname(os.path.realpath(__file__))

def msgError(msg):
    print(msg)
    exit()
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
        msgError("The output path is required in config")
    return config, sInfo

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

def sr_merge(strPkl,strH5ad):
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
      batch_key=batchKey)
    adata.write(strH5ad)
    
def sr_read(meta,sName,strPkl=None,strH5ad=None):
    print("Reading samples ...")
    if not strPkl is None and os.path.isfile(strPkl):
        print("The raw file exists: %s\n\tIf a new read is required, please rename/remove the above file!"%strPkl)
        sr_merge(strPkl,strH5ad)
        return
    all_slices = {}
    for i in range(meta.shape[0]):
        print("\t",meta[sName][i])
        oneD = sc.read_visium(meta[sr_path_column][i])
        oneD.var_names_make_unique()
        all_slices[meta[sName][i]]=oneD
    if strPkl is None:
        return all_slices
    with open(strPkl,"wb") as f:
        pickle.dump(all_slices,f)
    if not strH5ad is None and os.path.isfile(strH5ad):
        os.remove(strH5ad)
    sr_merge(strPkl,strH5ad)
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
