
import sys,os,warnings,logging,shutil
import cmdUtility as cu
import utility as ut

warnings.simplefilter(action='ignore', category=FutureWarning)
logging.disable(level=logging.INFO)

strPipePath=os.path.dirname(os.path.dirname(os.path.realpath(__file__)))

def saveRaw(meta,sName,strPkl,strH5ad):
    os.makedirs(os.path.dirname(strPkl),exist_ok=True)
    ut.sr_read(meta,sName,strPkl,strH5ad)
def main():
    if len(sys.argv)<2:
        print("Test running...")
        strConfig = os.path.join(strPipePath,"example","visium","config.yml")
    else:
        strConfig=os.path.realpath(sys.argv[1])
    config,sInfo = ut.getConfig(strConfig)
    ut.sr_checkMeta(sInfo,config)
    
    # read in samples
    strSep = os.path.join(config['output'],"raw",config['prj_name']+".pkl")
    strH5ad_raw = os.path.join(config['output'],"raw",config['prj_name']+".h5ad")
    cu.submit_cmd({'rawRead':"python -u %s/src/visium.py saveRaw %s %s %s"%(strPipePath,strConfig,strSep,strH5ad_raw)},config)
    if not os.path.isfile(strSep):
        ut.msgError("Error in reading files!")
    
    # apply SpaGCN
    strSpaGCN = os.path.join(config['output'],"SpaGCN",config['prj_name']+".pkl")
    cu.submit_cmd({'SpaGCN':"python -u %s/src/SpaGCN_run.py %s %s %s"%(strPipePath,strConfig,strSep,strSpaGCN)},config)
    if not os.path.isfile(strSpaGCN):
        ut.msgError("Error in spaGCN!")    
    
    # merge all
    strH5ad = os.path.join(config['output'],config['prj_name']+".h5ad")
    if os.path.isfile(strH5ad):
        print("Final h5ad file exists: ",strH5ad)
    else:
        shutil.copyfile(strH5ad_raw, strH5ad)
    cu.submit_cmd({'mergeSpaGCN':"python -u %s/src/SpaGCN_run.py %s %s"%(strPipePath,strH5ad,strSpaGCN)},config)
    
    print("=== visium process is completed! ===")
if __name__ == "__main__":
    if len(sys.argv)==1:
        print("=== Welcome to 'visium' from SpaceRequest! ===")
        print("\tPlease provide either a path to a folder or a config file.")
        print("\tAn empty config file will be created if a path to a folder is provided.")
        ut.msgError(ut.getSys().get('powerMsg'))
    elif len(sys.argv)>1 and sys.argv[1]=="saveRaw":
        config,sInfo = ut.getConfig(sys.argv[2])
        saveRaw(sInfo,config['sample_name'],sys.argv[3],sys.argv[4])
    else:
        main()
        ut.msgError(ut.getSys().get('powerMsg'))
