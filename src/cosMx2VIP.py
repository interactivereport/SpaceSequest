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
sample_col = 'FOV_id'
fov_img_col = 'fov_img'
cosMx_uns = 'cosMx'
cell_uns = 'cell_polygons'
transc_uns = 'tx_loc'
img_uns = 'images'
fov_prefix = 'FOV'

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
    for checkF in ['expression_file','meta_file','cell_polygon_file','transcript_loc_file','fov_image_file']:
        if config.get(checkF) is None or not os.path.isfile(config[checkF]):
            ut.msgError("Error: Missing %s: %s"%(checkF,config.get(checkF)))
    fov = pd.read_csv(config['fov_image_file'])
    for oneF in fov[fov_img_col]:
        if not os.path.isfile(oneF):
            ut.msgError("Error: Missing FOV image meta: %s"%oneF)
    return config,fov
def readData(strF,fov_col,possibleIDs,cID2index=False,suf_loc=None,rm_fov_col=False):
    X = pd.read_csv(strF)
    if not fov_col in X.columns:
        ut.msgError('\tError: missing fov coloumn (%s) in %s'%(fov_col,strF))
    cID_col = [one for one in possibleIDs if one in X.columns]
    if len(cID_col)==0:
        ut.msgError('\tError: cannot find cell id column in %s'%strF)
    print("\tSet %s as %s, and add fov column (%s) in cell ID"%(cID_col[0],cellID_col,fov_col))
    X.rename(columns={cID_col[0]:cellID_col},inplace=True)
    X[cellID_col] = fov_prefix+X[fov_col].astype('str')+"_"+X[cellID_col].astype('str')
    if rm_fov_col:
        X.drop(fov_col,axis=1,inplace=True)
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
def saveVIP(strConfig):
    config, fov = getConfig(strConfig)
    suffixLoc = dict(zip(suf_px.keys(),[config['suffix_local_px'],config['suffix_global_px']]))
    fov_col = fov.columns[0]
    print(" read expression")
    gExp = readData(config['expression_file'],fov_col,config['cell_id_col'],cID2index=True,rm_fov_col=True)
    print(" read cell meta information")
    cInfo = readData(config['meta_file'],fov_col,config['cell_id_col'],cID2index=True,suf_loc=suffixLoc)
    #cInfo[sample_col] = fov_prefix+cInfo[fov_col]
    cID = cInfo.index.intersection(gExp.index)
    print(" create AnnData")
    D = ad.AnnData(X=gExp.loc[cID,:],obs=cInfo.loc[cID,:])
    D.uns[cosMx_uns] = {}
    for i in D.obs[fov_col].astype('str').unique():
        D.uns[cosMx_uns][fov_prefix+i] = {}
    del gExp
    print(" add px coordinates")
    for k in suf_px:
        if k=='local' and not config['saveLocalX']:
            continue
        D.obsm['X_%s'%k] = cInfo[["%s%s"%(i,suf_px[k]) for i in ['x','y']]].to_numpy('float32')
    print(" add cell boundaries")
    cellBD = readData(config['cell_polygon_file'],fov_col,config['cell_id_col'],suf_loc=suffixLoc)
    cellBD = cellBD[cellBD[fov_col].isin(D.obs[fov_col])]
    for i in cellBD[fov_col].astype('str').unique():
        D.uns[cosMx_uns][fov_prefix+i][cell_uns] = cellBD[cellBD[fov_col].astype('str')==i].drop(fov_col,axis=1).reset_index(drop=True)
    del cellBD
    print(" add transcript coordinates")
    trX = readData(config['transcript_loc_file'],fov_col,config['cell_id_col'],suf_loc=suffixLoc)
    trX = trX[trX[fov_col].isin(D.obs[fov_col])]
    if not config['tx_col'] in trX.columns:
        ut.msgError('\tError: missing transcript column %s in %s'%(config['tx_col'],config['transcript_loc_file']))
    for i in trX[fov_col].astype('str').unique():
        D.uns[cosMx_uns][fov_prefix+i][transc_uns] = trX[trX[fov_col].astype('str')==i].drop(fov_col,axis=1).reset_index(drop=True)
    del trX
    print(" add FOV image")
    fov.index = fov[fov_col].astype('str')
    for i in range(fov.shape[0]):
        if fov[fov_col][i] in D.obs[fov_col].values:
            D.uns[cosMx_uns][fov_prefix+str(fov[fov_col][i])][img_uns]=plt.imread(fov[fov_img_col][i])

    print("Saving ...")
    D.uns[cosMx_uns]['keys'] = {'%s_px'%i:{j:'%s%s'%(j,suf_px[i]) for j in ['x','y']} for i in suf_px}
    D.uns[cosMx_uns]['keys']['cell'] = cell_uns
    D.uns[cosMx_uns]['keys']['tx_loc'] = transc_uns
    D.uns[cosMx_uns]['keys']['img'] = img_uns
    D.uns[cosMx_uns]['keys']['tx_col'] = config['tx_col']
    D.write(os.path.join(config['output'],config['prj_name']+".h5ad"))
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
