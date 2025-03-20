import sys, os, time, random, warnings, logging, re #, skimage
warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.simplefilter(action='ignore', category=UserWarning)
logging.disable(level=logging.INFO)

import anndata as ad
import cmdUtility as cu
import utility as ut
from scipy.sparse import csc_matrix
from rpy2 import robjects as robj
from rpy2.robjects import pandas2ri

strPipePath=os.path.dirname(os.path.realpath(__file__))
c2l_path='c2l'
mKey="SpaTalk"

def check_config(config,sInfo):
    if config.get("st_scH5ad") is None or config.get("st_annotation_obs") is None:
        print("___ Skip SpaTalk: no st_scH5ad (or st_annotation_obs) is provided! ___")
        return False
    if not os.path.isfile(config["st_scH5ad"]):
        print("___ Skip SpaTalk: st_scH5ad (%s) does not exist! ___"%config["st_scH5ad"])
        return False
    D = ad.read_h5ad(config["tg_scH5ad"],backed="r")
    if not config["st_annotation_obs"] in D.obs.columns:
        print("___ Skip SpaTalk: st_annotation_obs (%s) does not exist in scH5ad! ___"%config["st_annotation_obs"])
        return False
    if not config.get("st_matchColumn") is None:
        if not config["st_matchColumn"] in D.obs.columns:
            print("___ Skip tangram: st_matchColumn (%s) does not exist in scH5ad! ___"%config["st_matchColumn"])
            return False
        if not config["st_matchColumn"] in sInfo.columns:
            print("___ Skip tangram: st_matchColumn (%s) does not exist in sample meta! ___"%config["st_matchColumn"])
            return False
        missName = [str(one) for one in sInfo[config["st_matchColumn"]] if len(str(one))>0 and not D.obs[config["st_matchColumn"]].eq(one).any() and not D.obs[config["st_matchColumn"]].eq(str(one)).any()]
        if len(missName)>0:
            print("___ Skip tangram: missing st_matchColumn values (%s) in scH5ad! ___"%', '.join(missName))
            return False
    #st_species: Mouse or Human
    if not config.get("st_species") in ['Mouse','Human']:
        print("Error: unknown st_species (%s), only support Mouse or Human"%config.get("st_species"))
        return False
    if config.get('st_use_cell2location') and not 'cell2location' in config['methods']:
        print("Error: cell2location is NOT set in methods but st_use_cell2location: True")
        return False
    return True

def run(strConfig,strH5ad):
    config,sInfo = ut.getConfig(strConfig)
    strOut = os.path.join(config['output'],mKey)
    os.makedirs(strOut,exist_ok=True)
    strPkl=os.path.join(strOut,config['prj_name']+".pkl")
    
    if os.path.isfile(strPkl):
        print("\n\nUsing previous *SpaTalk* results: %s\n\tIf a new process is wanted, please rename/remove the above file"%strPkl)
    else:
        strRDS = re.sub("pkl","rds",strPkl)
        if not os.path.isfile(strRDS):
            if not check_config(config,sInfo):
                return
            strC2L=''
            if config.get('st_use_cell2location'):
                strC2L = os.path.join(config['output'],c2l_path,config['prj_name']+'.pkl')
                st = time.time()
                while not os.path.isfile(strC2L) and (time.time()-st)<3600*10: #max 10 hours delay
                    time.sleep(60)
                if os.path.isfile(strC2L):
                    C2L = ut.readPkl(strC2L)
                    C2L.columns = [re.sub('^c2l_','',one) for one in C2L.columns.to_list()]
                else:
                    print("___ Skip SpaTalk: C2L results (%s) is NOT available ___! "%strC2L)
                    return
                strC2L = os.path.join(strOut,config['prj_name']+'_C2L.h5ad')
                A = ad.AnnData(C2L)
                A.X = csc_matrix(A.X)
                A.write(strC2L)
            cmd = {}
            cmd['SpaTalk'] = "Rscript %s/SpaTalk_run.R %s %s %s %s %s"%(strPipePath,strConfig,strH5ad,strOut,strRDS,strC2L)
            cu.submit_cmd(cmd,config)
        if not os.path.isfile(strRDS):
            print("___ Error SpaTalk: missing result file %s! "%strRDS)
            return
        res = robj.r('readRDS("%s")'%strRDS)
        resKey = list(res.names)
        saveObj={}
        for i in range(len(resKey)):
            saveObj[resKey[i]] = robj.pandas2ri.rpy2py_dataframe(res[i])
        ut.writePkl(saveObj,strPkl)
    return {mKey:strPkl}

def merge(D,allRes):
    if not mKey in allRes.keys():
        return
    strF = allRes[mKey]
    print("merging %s: %s"%(mKey,strF))
    if not os.path.isfile(strF):
        print("\tSkip, the above file is missing!")
        return
    #D = ad.read_h5ad(strH5ad)#,backed="r+"
    if mKey in D.uns.keys():
        print("Over-writing exitsting uns key: %s"%mKey)
    D.uns[mKey] = ut.readPkl(strF)
    #print("\tsaving")
    #D.write(strH5ad)
