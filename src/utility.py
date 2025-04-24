import sys,os,yaml,warnings,logging,functools,pickle,configparser,shutil,h5py,re,math
from datetime import datetime
import pandas as pd
import scanpy as sc
import numpy as np
import anndata as ad
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
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

def msgError(msg=""):
    if not msg.startswith('Questions'):
        print("Error:",end=" ")
    print("Error: ",msg)
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
    if config.get("sample_name") is not None and config["sample_name"] in sInfo.columns:
        sInfo['Sample_Name'] = sInfo['Sample_Name'].astype(str)
        sInfo.index = sInfo['Sample_Name'].tolist() #list(sInfo[config["sample_name"]])
    else:
        msgError("Missing sample name column (%s)"%config.get("sample_name"))
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

def writeHDF(df,strF):
    df.columns = [re.sub(" ","_",_) for _ in df.columns]
    try:
        os.remove(strF)
    except:
        pass
    with h5py.File(strF,'w') as f:
        f.create_dataset('axis_row',data=list(df.index))
        f.create_dataset('axis_col',data=list(df.columns))
        for one in df.columns:
            print(one)
            if isinstance(df[one].dtype, pd.CategoricalDtype):
                g = f.create_group(one)
                g.create_dataset('value', data=list(df[one].cat.categories))
                g.create_dataset('key', data=[1+_ for _ in df[one].cat.codes])
            elif pd.api.types.is_bool_dtype(df[one]):
                f.create_dataset(one, data=list(df[one].astype(int)))
            else:
                f.create_dataset(one, data=list(df[one]))
def strDecode(x):
    if isinstance(x[0],bytes):
        return([_.decode('utf-8') for _ in x])
    else:
        return(list(x))
def readHDF(strF):
    with h5py.File(strF,'r') as f:
        rName = strDecode(f['axis_row'])
        cName = strDecode(f['axis_col'])
        df = {}
        for one in f.keys():
            if one in ['axis_row','axis_col']:
                continue
            if isinstance(f[one], h5py.Dataset):
                df[one] = strDecode(f[one])
                if set(df[one]).issubset({0,1}):
                    df[one] = [bool(i) for i in df[one]]
            elif isinstance(f[one], h5py.Group):
                if all(_ in f[one].keys() for _ in ['key','value']):
                    val = strDecode(f[one]['value'])
                    df[one] = pd.Categorical([val[int(_)-1] for _ in f[one]['key']],val)
            else:
                msgError('readHDF: unknow type '+str(type(f[one])))
    return pd.DataFrame(df,index=rName)[cName]

def writeFeather(df,strF):
    df.reset_index().rename(columns={'index': 'rowNames'}).to_feather(strF)
def readFeather(strF):
    D = pd.read_feather(strF)
    if 'rowNames' in D.columns:
        D = D.set_index('rowNames')
        D.index.name = None
    return D
    

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
    if D.obs[batchKey].dtype == 'category':
        asType = D.obs[batchKey].cat.categories.dtype
    else:
        asType = D.obs[batchKey].dtype
    adata = D[D.obs[batchKey]==asType.type(sName)].copy()
    adata.uns['spatial'] = {sName:adata.uns['spatial'][sName]}
    return adata
def sr_merge(strPkl,strH5ad,config):
    if strH5ad is None or strPkl is None:
        return
    if not strH5ad is None and os.path.isfile(strH5ad):
        print("The raw merged file exists: %s\n\tIf a new merge is required, please rename/remove the above file!"%strH5ad)
        return
    D_dict = readPkl(strPkl)
    D_list = list(D_dict.values())
    #adata = sc.AnnData.concatenate(*D_list,
    #  join="outer",
    # batch_categories=list(D_dict.keys()),
    # batch_key=batchKey,
    #  uns_merge="unique")
    adata = ad.concat(D_list,
        join="outer",
        keys=list(D_dict.keys()),
        label=batchKey,
        merge="unique",
        uns_merge="unique",
        index_unique="_")
    sr_filter(adata,config)
    sr_plotQC(adata,os.path.join(config['output'],"QC","QC.pdf"))
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
            print("\t Spot, %s cutoff %d: %d spots"%(one,config[groupKey][one]["cutoff"],selObs.sum()))
        adata._inplace_subset_obs(selObs)
    else:
        sc.pp.calculate_qc_metrics(adata,inplace=True)
    if rmGene.sum()>0:
        adata._inplace_subset_var(np.invert(rmGene))
        print("\t Total of %d genes are removed",rmGene.sum())
def sr_filter(D,config):
    print("*** filter ***")
    print("\t Total: %d spots with %d genes"%(D.shape))
    sc.pp.filter_genes(D,min_counts=config.get('UMI_count_gene'))
    print("\t Gene, min UMI %d: %d spots with %d genes"%((config.get('UMI_count_gene'),)+D.shape))
    sc.pp.filter_genes(D,min_cells=config.get('min_cells'))
    print("\t Gene, min cells %d: %d spots with %d genes"%((config.get('min_cells'),)+D.shape))
    sc.pp.filter_cells(D,min_counts=config.get('UMI_count_cell'))
    print("\t Cell, min UMI %d: %d spots with %d genes"%((config.get('UMI_count_cell'),)+D.shape))
    sc.pp.filter_cells(D,min_genes=config.get('min_features'))
    print("\t Cell, min gene %d: %d spots with %d genes"%((config.get('min_features'),)+D.shape))
    sc.pp.filter_cells(D,max_counts=config.get('highCount_cutoff'))
    print("\t Cell, max UMI %d: %d spots with %d genes"%((config.get('highCount_cutoff'),)+D.shape))
    sc.pp.filter_cells(D,max_genes=config.get('highGene_cutoff'))
    print("\t Cell, max gene %d: %d spots with %d genes"%((config.get('highGene_cutoff'),)+D.shape))
    sr_addGP(D,config)
    return
def sr_metrics(meta,config):
    df = []
    for id in meta.index:
        strF = os.path.join(meta[sr_path_column][id],"metrics_summary.csv")
        if os.path.isfile(strF):
            df.append(pd.read_csv(strF))
    if len(df)==0:
        print("No metrics information is available!")
        return
    strMetrics = os.path.join(config['output'],"QC","metrics_summary.csv")
    os.makedirs(os.path.dirname(strMetrics),exist_ok=True)
    pd.concat(df,ignore_index=True).to_csv(strMetrics,index=False)
def sr_plotQC(D,strPDF):
    os.makedirs(os.path.dirname(strPDF),exist_ok=True)
    w = max(6.4,D.obs[batchKey].nunique()*2/10)
    plt.rcParams["figure.figsize"] = (w,4.8)
    with PdfPages(strPDF) as pdf:
        for qc in ["total_counts","n_genes_by_counts"]+[_ for _ in D.obs.columns if _.startswith('pct_counts_')]:
            ax = sc.pl.violin(D,keys=qc,groupby=batchKey,rotation=90,show=False)
            #ax.get_legend().remove()
            #sns.displot(data=D.obs,x=qc,hue=batchKey,kind="kde")
            pdf.savefig(bbox_inches="tight")
            plt.close()
    plt.rcParams['figure.figsize'] = plt.rcParamsDefault['figure.figsize']
def sr_read(meta,config,strPkl=None,strH5ad=None):
    sName = config['sample_name']
    print("Reading samples ...")
    if not strPkl is None and os.path.isfile(strPkl):
        print("The raw file exists: %s\n\tIf a new read is required, please rename/remove the above file!"%strPkl)
        sr_merge(strPkl,strH5ad,config)
        return
    all_slices = {}
    for i in meta.index:
        print("\t",meta[sName][i])
        oneD = sc.read_visium(meta[sr_path_column][i],library_id=meta[sName][i])
        oneD.var_names_make_unique()
        if config.get('UMI_dtype') is None:
            config['UMI_dtype'] = 'int32'
        print("\t\tMake UMI as X %s!"%config['UMI_dtype'])
        oneD.X = oneD.X.astype(config['UMI_dtype'])
        if 'spatial' in oneD.obsm.keys():
            oneD.obsm['X_spatial'] = oneD.obsm['spatial']
        #    del oneD.obsm['spatial'] # cannot delete, since the SpaGCN needs it 
        sr_addGP(oneD,config)
        #add additional sample information
        for j in meta.columns:
            if j in [sName,sr_path_column]:
                continue
            oneD.obs[j] = meta[j][i]
        all_slices[meta[sName][i]]=oneD # meta[sName] is always str in getConfig step
    sr_metrics(meta,config)
    if strPkl is None:
        return all_slices
    writePkl(all_slices,strPkl)
    if not strH5ad is None and os.path.isfile(strH5ad):
        os.remove(strH5ad)
    sr_merge(strPkl,strH5ad,config)
    return

## plot annotation
def adjustAlpha(x,adjA):
    if x.max()<=0:
        return adjA
    x = x/x.max() #(alpha-min(alpha))/(max(alpha)-min(alpha))
    x[x==0] = x[x!=0].min()/2
    a = 1-adjA
    b = 20**a
    alpha =(1+a**b)/(1+(a/x)**b)
    return alpha
def plotVisiumOne(adata,sid,col,ax,fig,alpha=1,nMap='viridis',cMap='Set1',dotsize=4):
    keys = {'slide_column':'library_id','img':['images','hires'],'scale':['scalefactors','tissue_hires_scalef'],'coordinates':'X_spatial'}
    img = adata.uns['spatial'][sid]
    for v in keys['img']:
        img = img[v]
    scaler = adata.uns['spatial'][sid]
    for v in keys['scale']:
        scaler = scaler[v]
    subD = adata[adata.obs[keys['slide_column']]==sid,]
    a=ax.imshow(img)
    x = subD.obsm[keys['coordinates']][:,0]*scaler
    y = subD.obsm[keys['coordinates']][:,1]*scaler
    g_column = None
    if 'feature_name' in subD.var.columns:
        g_column = 'feature_name'
    elif 'name_0' in subD.var.columns:
        g_column = 'name_0'
    df = sc.get.obs_df(subD,[col],gene_symbols=g_column)
    if pd.api.types.is_numeric_dtype(df[col]):
        a=ax.scatter(x,y,
            c=df[col].to_numpy(),cmap=nMap,s=dotsize,
            alpha=adjustAlpha(df[col].to_numpy(),alpha))
        a=fig.colorbar(a,ax=ax)
    else:
        try:
            grps = sorted(adata.obs[col].unique().tolist(),key=int)
        except:
            grps = sorted(adata.obs[col].unique().tolist())
        colors = dict(zip(grps,sns.color_palette(cMap,n_colors=len(grps)).as_hex()))#husl
        a=ax.scatter(x,y,
            color=[colors[_] for _ in df[col]],
            s=dotsize,alpha=alpha)
        legend_handles = [Line2D([],[],marker='.',color=colors[i],linestyle='None') for i in colors.keys()]
        a=ax.legend(handles=legend_handles,labels=colors.keys(),ncols=math.ceil(len(grps)/11),
            bbox_to_anchor=(1.05, 1),loc='upper left')
    a=ax.set_title("%s: %s"%(sid,col))
def plotVisiumFun(adata,strPDF=None,obs=None,genes=None,alpha=0.9,subSize=4,ncol=4,nMap='viridis',cMap='Set3',dotsize=2,sIDs=None,sortID=True):
    sel = []
    if obs is not None:
        sel += adata.obs.columns[adata.obs.columns.isin(obs)].tolist()
    if genes is not None:
        sel += adata.var_names[adata.var_names.isin(genes)].tolist()
    if len(sel)==0:
        print("\t Nothing select to plot")
        return
    nSel = len(sel)
    #keys = adata.uns['visium']['keys']
    slides = adata.obs[batchKey].unique() if sIDs is None else sIDs
    nSlides = len(slides)
    if sortID:
        nrowSection = math.ceil(nSel/ncol)
        nrow = nSlides*nrowSection
        fig, axs = plt.subplots(nrow, ncol, figsize=(ncol*subSize,nrow*subSize))
        if isinstance(axs,np.ndarray):
            axs_flat = axs.flatten()
        else:
            axs_flat = np.array(axs)
        for i in range(nSlides):
            for j in range(nSel):
                plotVisiumOne(adata,slides[i],sel[j],
                    axs_flat[i*nrowSection*ncol+j],fig,
                    alpha=alpha,
                    nMap=nMap,
                    cMap=cMap,
                    dotsize=dotsize)
    else:
        nrowSection = math.ceil(nSlides/ncol)
        nrow = nSel*nrowSection
        fig, axs = plt.subplots(nrow, ncol, figsize=(ncol*subSize,nrow*subSize))
        if isinstance(axs,np.ndarray):
            axs_flat = axs.flatten()
        else:
            axs_flat = np.array(axs)
        for j in range(nSel):
            for i in range(nSlides):
                plotVisiumOne(adata,slides[i],sel[j],
                    axs_flat[j*nrowSection*ncol+i],fig,
                    alpha=alpha,
                    nMap=nMap,
                    cMap=cMap,
                    dotsize=dotsize)
    for ax in axs_flat:
        a=ax.axis('off')
    if strPDF is None:
        return fig
    else:
        fig.savefig(strPDF,bbox_inches="tight")
        plt.close()
def plotVisium(adata,strOut=None,obs=None,genes=None,alpha=0.9,subSize=4,ncol=4,nMap='viridis',cMap='Set3',dotsize=2):
    for one in adata.obs[batchKey].unique():
        plotVisiumFun(adata,os.path.join(strOut,one+".pdf"),obs=obs,genes=genes,alpha=alpha,subSize=subSize,ncol=ncol,nMap=nMap,cMap=cMap,dotsize=dotsize,sIDs=[one])
    
## standard logNormal normalization
def logNormal(strH5ad,config):
    logH5ad = os.path.join(config['output'],"logNormal",config['prj_name']+".h5ad")
    if os.path.isfile(logH5ad):
        print("Using previous *logNormal*: %s\n\tIf a new process is wanted, please rename/remove the above file"%logH5ad)
        #shutil.copy(logH5ad,strOut)
        return {'logNormal':logH5ad}
    print("*** logNormal ***")
    os.makedirs(os.path.dirname(logH5ad),exist_ok=True)
    D = sc.read_h5ad(strH5ad)
    D.raw = D.copy()
    sc.pp.highly_variable_genes(D,flavor='seurat_v3',n_top_genes=3000,inplace=True)
    sc.pp.normalize_total(D,target_sum=config['normScale'])
    sc.pp.log1p(D)
    sc.pp.pca(D, n_comps=50, use_highly_variable=True, svd_solver='arpack')
    sc.pp.neighbors(D)
    sc.tl.umap(D)
    D.obsm["X_umap_lognormal"] = D.obsm["X_umap"].copy()
    del D.obsm["X_umap"]
    sc.tl.louvain(D,key_added='lognormal_louvain_clusters')
    # harmony
    sc.external.pp.harmony_integrate(D,key=batchKey,max_iter_harmony=50)#X_pca_harmony
    sc.pp.neighbors(D,use_rep="X_pca_harmony")
    sc.tl.umap(D)
    D.obsm["X_umap_harmony"] = D.obsm["X_umap"].copy()
    del D.obsm["X_umap"]
    sc.tl.louvain(D,key_added='harmony_louvain_clusters')
    D.write(logH5ad)
    return {'logNormal':logH5ad}
    #shutil.copy(logH5ad,strOut)

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
