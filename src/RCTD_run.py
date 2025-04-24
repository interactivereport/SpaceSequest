import os
import cmdUtility as cu
import anndata as ad
import pandas as pd
import utility as ut

strPipePath=os.path.dirname(os.path.realpath(__file__))
mKey="RCTD"

def run(strConfig,strH5ad_raw):
    config,sInfo = ut.getConfig(strConfig)
    strOut=os.path.join(config['output'],mKey,config['prj_name']+".feather")
    if os.path.isfile(strOut):
        print("\n\nUsing previous *%s* results: %s\n\tIf a new process is wanted, please rename/remove the above file"%(mKey,strOut))
        return {mKey:strOut}
    else:
        cmd = {}
        cmd[mKey] = "Rscript %s/RCTD_run.R %s %s %s"%(strPipePath,strConfig,strH5ad_raw,strOut)
        cu.submit_cmd(cmd,config)
    return {mKey:strOut}

## merge obs
def merge(D,allRes):
    if not mKey in allRes.keys():
        return
    strObs = allRes[mKey]
    print("merging %s: %s"%(mKey,strObs))
    if not os.path.isfile(strObs):
        print("\tSkip: the above file is missing!")
        return
    #D = ad.read_h5ad(strH5ad)#,backed="r+"
    obs = ut.readFeather(strObs)
    #obs = pandas2ri.rpy2py_dataframe(readRDS(strObs))
    selCol=~obs.columns.isin(D.obs.columns)
    if (~selCol).sum()>0:
        print("\tSkip exists:",",".join(obs.columns[~selCol]))
    if selCol.sum()>0:
        D.obs = D.obs.merge(obs.loc[:,selCol],"left",left_index=True,right_index=True)
        modCol = D.obs.columns.isin(obs.columns[selCol])
        D.obs.loc[:,modCol]=D.obs.loc[:,modCol].apply(lambda x: ut.fillNA(x,0))
        ut.plotVisium(D,strOut=os.path.dirname(strObs),obs=D.obs.columns[modCol].tolist())
        #D.write(strH5ad)

