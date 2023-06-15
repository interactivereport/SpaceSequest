
import sys,os,warnings,logging,shutil,random,re
import cmdUtility as cu
import utility as ut

warnings.simplefilter(action='ignore', category=FutureWarning)
logging.disable(level=logging.INFO)

strPipePath=os.path.dirname(os.path.dirname(os.path.realpath(__file__)))

def init(strDir):
    strMeta=os.path.join(strDir,"sample.csv")
    with open(strMeta,"w") as f:
        f.write("Sample_Name,%s\n"%ut.sr_path_column)
    config = ut.readLines(os.path.join(strPipePath,"src","visium.yml"))
    config = [one.replace("OUTPUT",strDir)
                    .replace("SAMPLE",strMeta)
                    .replace("JOBID","j%d"%random.randint(10,99)) for one in config]
    strConfig = os.path.join(strDir,"config.yml")
    with open(strConfig,"w") as f:
        f.writelines(config)
    print("Config file is created: ",strConfig)
def saveRaw(meta,sName,strPkl,strH5ad):
    os.makedirs(os.path.dirname(strPkl),exist_ok=True)
    ut.sr_read(meta,sName,strPkl,strH5ad)
def main():
    strPath = os.path.realpath(sys.argv[1])
    if os.path.isdir(strPath):
        init(strPath)
        return()
    elif os.path.isfile(strPath):
        strConfig=strPath
    else:
        print("Unknown input: ",strPath)
        return()
    
    ut.MsgInit()
    config,sInfo = ut.getConfig(strConfig)
    ut.sr_checkMeta(sInfo,config)
    
    # read in samples
    strSep = os.path.join(config['output'],"raw",config['prj_name']+".pkl")
    strH5ad_raw = os.path.join(config['output'],"raw",config['prj_name']+".h5ad")
    cu.submit_cmd({'rawRead':"python -u %s/src/visium.py saveRaw %s %s %s"%(strPipePath,strConfig,strSep,strH5ad_raw)},config)
    if not os.path.isfile(strSep):
        ut.msgError("Error in reading files!")
    
    cmd = {}
    strFinals = {}
    # apply SpaGCN
    strFinals["SpaGCN"] = os.path.join(config['output'],"SpaGCN",config['prj_name']+".pkl")
    if os.path.isfile(strFinals["SpaGCN"]):
        print("\n\nUsing previous SpaGCN results: %s\n\tIf a new process is wanted, please rename/remove the above file"%strFinals["SpaGCN"])
    else:
        cmd['SpaGCN'] = "python -u %s/src/SpaGCN_run.py %s %s %s"%(strPipePath,strConfig,strSep,strFinals["SpaGCN"])
    
    # apply BayesSpace
    strFinals["BayesSpace"] = os.path.join(config['output'],"BayesSpace",config['prj_name']+".rds")
    if os.path.isfile(strFinals["BayesSpace"]):
        print("\n\nUsing previous BayesSpace results: %s\n\tIf a new process is wanted, please rename/remove the above file"%strFinals["BayesSpace"])
    else:
        cmd['BayesSpace'] = "Rscript %s/src/BayesSpace_run.R %s %s"%(strPipePath,strConfig,strFinals["BayesSpace"])
    
    # Run all methods
    cu.submit_cmd(cmd,config)
    for one in strFinals:
        if not os.path.isfile(strFinals[one]):
            ut.msgError("Error in %s!"%one)

    # merge all
    strH5ad = os.path.join(config['output'],config['prj_name']+".h5ad")
    if os.path.isfile(strH5ad):
        print("Final h5ad file exists: ",strH5ad)
    else:
        shutil.copyfile(strH5ad_raw, strH5ad)
    cu.submit_cmd({'mergeSpaGCN':"python -u %s/src/SpaGCN_run.py %s %s"%(strPipePath,strH5ad,strFinals["SpaGCN"])},config)
    cu.submit_cmd({'mergeBayesSpace':"python -u %s/src/utility.py mergeRDS %s %s"%(strPipePath,strH5ad,strFinals["BayesSpace"])},config)
    
    print("\n\n=== visium process is completed! ===")
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
