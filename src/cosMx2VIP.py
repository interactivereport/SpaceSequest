import sys,os,warnings,logging,shutil,random,re,functools,yaml,pickle
import utility as ut
import pandas as pd
import anndata as ad
import matplotlib.pyplot as plt

warnings.simplefilter(action='ignore', category=FutureWarning)
logging.disable(level=logging.INFO)

strPipePath=os.path.dirname(os.path.dirname(os.path.realpath(__file__)))

meta_cols={'exp_col':'exp_path',
    'meta_col':'meta_path',
    'cell_col':'cell_def_path',
    'tx_col':'transcript_loc_path',
    'his_col':'histology_img_path'

}
cellID_col = 'cell_ID'
suf_px={'local':'_local_px','global':'_global_px'}
sample_col = 'library_id'
cosMx_uns = 'cosMx'
cell_uns = 'cell_polygons'
transc_uns = 'tx_loc'
img_uns = 'images'

def init(strDir):
    strMeta=os.path.join(strDir,"sample.csv")
    with open(strMeta,"w") as f:
        f.write("Sample_Name")
        for one in meta_cols.values():
            f.write(",%s"%one)
        f.write("\n")
    config = ut.readLines(os.path.join(strPipePath,"src","cosMx.yml"))
    config = [one.replace("OUTPUT",strDir)
                    .replace("SAMPLE",strMeta)
                    .replace("JOBID","j%d"%random.randint(10,99)) for one in config]
    strConfig = os.path.join(strDir,"config.yml")
    with open(strConfig,"w") as f:
        f.writelines(config)
    print("Config file is created: ",strConfig)
def getConfig(strConfig):
    print("Checking config ...")
    if not os.path.isfile(strConfig):
        ut.msgError("Missing config: %s"%strConfig)    
    with open(strConfig,"r") as f:
        config = yaml.safe_load(f)
    if config.get("sample_meta") is not None and os.path.isfile(config["sample_meta"]):
        sInfo = pd.read_csv(config["sample_meta"])
    else:
        ut.msgError("Error: Missing sample meta: %s"%config.get("sample_meta"))
    for one in meta_cols.values():
        if not one in sInfo.columns:
            ut.msgError('Error: Missing requied column: %s'%one)
    if config.get("sample_name") is not None and config["sample_name"] in sInfo.columns:
        sInfo.index = list(sInfo[config["sample_name"]])
    else:
        ut.msgError("Error: Missing sample name column (%s)"%config.get("sample_name"))
    return config,sInfo
def readData(strF,sID,possibleIDs,cID2index=False,suf_loc=None):
    if not os.path.isfile(strF):
        ut.msgError('\tError: missing file %s'%strF)
    X = pd.read_csv(strF)
    cID_col = [one for one in possibleIDs if one in X.columns]
    if len(cID_col)==0:
        ut.msgError('\tError: cannot find cell id column in %s'%strF)
    print("\tSet %s as %s, and add sample id (%s) in cell ID"%(cID_col[0],cellID_col,sID))
    X.rename(columns={cID_col[0]:cellID_col},inplace=True)
    X[cellID_col] = sID+"_"+X[cellID_col].astype('str')
    print("\tUpdate px coordinates columns")
    if cID2index:
        print("\tSet %s as index and drop the column"%cellID_col)
        X.index=list(X[cellID_col])
        X.drop(cellID_col,axis=1,inplace=True)
    if suf_loc is not None:
        print("\tCheck if coordinates are available")
        cols = X.columns
        loc_cols=[]
        for k in suf_loc:
            for a in ['x','y']:
                cols = [re.sub('.*%s.*%s'%(a,suf_loc[k]),'%s%s'%(a,suf_px[k]),one,flags=re.I) for one in cols]
                loc_cols.append('%s%s'%(a,suf_px[k]))
        X.columns = cols
        if(any(x not in X.columns for x in loc_cols)):
            ut.msgError('\tError: missing coordinates in file %s'%strF)
    return X
def readOneCapture(sID,meta,config):
    print("Reading ",sID," ...")
    suffixLoc = dict(zip(suf_px.keys(),[config['suffix_local_px'],config['suffix_global_px']]))
    # exp
    print(" exp")
    gExp = readData(meta[meta_cols['exp_col']][sID],sID,
        config['cell_id_col'],cID2index=True)
    # meta
    print(" meta")
    cInfo = readData(meta[meta_cols['meta_col']][sID],sID,
        config['cell_id_col'],cID2index=True,suf_loc=suffixLoc)
    for one in meta.columns:
        if one in [config['sample_name']]+list(meta_cols.values()):
            continue
        cInfo[one] = meta[one][sID]
    cID = cInfo.index.intersection(gExp.index)
    # create anndata
    print(" create AnnData")
    D = ad.AnnData(X=gExp.loc[cID,:],obs=cInfo.loc[cID,:])
    D.uns[cosMx_uns] = {sID:{}}
    del gExp
    # add px coordinates
    print(" add px coordinates")
    for k in suf_px:
        if k=='local' and not config['saveLocalX']:
            continue
        D.obsm['X_%s'%k] = cInfo[["%s%s"%(i,suf_px[k]) for i in ['x','y']]].to_numpy('float32')
    # add cell boundaries
    print(" add cell boundaries")
    D.uns[cosMx_uns][sID][cell_uns] = readData(meta[meta_cols['cell_col']][sID],sID,
        config['cell_id_col'],suf_loc=suffixLoc).reset_index(drop=True)
    # add transcript coordinates
    print(" add transcript coordinates")
    D.uns[cosMx_uns][sID][transc_uns] = readData(meta[meta_cols['tx_col']][sID],sID,
        config['cell_id_col'],suf_loc=suffixLoc).reset_index(drop=True)
    if not config['tx_col'] in D.uns[cosMx_uns][sID][transc_uns].columns:
        ut.msgError('\tError: missing transcript column %s in %s'%(config['tx_col'],strF))
    # add histology image
    if not os.path.isfile(meta[meta_cols['his_col']][sID]):
        ut.msgError('\tError: missing histology image file %s'%meta[meta_cols['his_col']][sID])
    D.uns[cosMx_uns][sID][img_uns]=plt.imread(meta[meta_cols['his_col']][sID])
    return D
def saveVIP(strConfig):
    config, meta = getConfig(strConfig)
    adatas = []
    for sID in meta.index:
        adatas.append(readOneCapture(sID,meta,config))
    if len(adatas)==1:
        adata = adatas[0]
        adata.obs[batchKey]=list(meta.index)[0]
    else:   
        adata = ad.AnnData.concatenate(*adatas,
            join="outer",
            batch_categories=list(meta.index),
            batch_key=sample_col,
            index_unique=None,
            uns_merge='unique')
    print("Saving ...")
    adata.uns[cosMx_uns]['keys'] = {'%s_px'%i:{j:'%s%s'%(j,suf_px[i]) for j in ['x','y']} for i in suf_px}
    adata.uns[cosMx_uns]['keys']['cell'] = cell_uns
    adata.uns[cosMx_uns]['keys']['tx_loc'] = transc_uns
    adata.uns[cosMx_uns]['keys']['img'] = img_uns
    adata.uns[cosMx_uns]['keys']['tx_col'] = config['tx_col']
    adata.write(os.path.join(config['output'],config['prj_name']+".h5ad"))
    return()
    
def main():
    strPath = os.path.realpath(sys.argv[1])
    if os.path.isdir(strPath):
        init(strPath)
    elif os.path.isfile(strPath):
        saveVIP(strPath)
    else:
        print("Unknown input: ",strPath)

    print("\n\n=== cosMx2VIP process is completed! ===")
if __name__ == "__main__":
    if len(sys.argv)==1:
        print("=== Welcome to 'cosMx2VIP' from SpaceRequest! ===")
        print("This application will process spatial data from cosMx into h5ad which can be visualized in cellxgene VIP.")
        print("\tPlease provide either a path to a folder or a config file.")
        print("\tAn empty config file will be created if a path to a folder is provided.")
        ut.msgError(ut.getSys().get('powerMsg'))
    else:
        main()
        ut.msgError(ut.getSys().get('powerMsg'))
