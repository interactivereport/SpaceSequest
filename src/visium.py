import sys,os,warnings,logging,shutil,random,re,functools
import cmdUtility as cu
import utility as ut
import SpaGCN_run as spa
import BayesSpace_run as bay
import Cell2location_run as c2l


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
def saveRaw(meta,config,strPkl,strH5ad):
    os.makedirs(os.path.dirname(strPkl),exist_ok=True)
    ut.sr_read(meta,config,strPkl,strH5ad)
def main():
    strPath = os.path.realpath(sys.argv[1])
    if os.path.isdir(strPath):
        init(strPath)
        return
    elif os.path.isfile(strPath):
        strConfig=strPath
    else:
        print("Unknown input: ",strPath)
        return
    
    ut.MsgInit()
    config,sInfo = ut.getConfig(strConfig)
    ut.sr_checkMeta(sInfo,config)
    
    strH5ad = os.path.join(config['output'],config['prj_name']+".h5ad")
    if os.path.isfile(strH5ad):
        print("Final h5ad file exists: ",strH5ad)
        return
    # read in samples
    strSep = os.path.join(config['output'],"raw",config['prj_name']+".pkl")
    strH5ad_raw = os.path.join(config['output'],"raw",config['prj_name']+".h5ad")
    cu.submit_cmd({'rawRead':"python -u %s/src/visium.py saveRaw %s %s %s"%(strPipePath,strConfig,strSep,strH5ad_raw)},config)
    if not os.path.isfile(strSep):
        ut.msgError("Error in reading files!")
    
    methods = []
    # logNormal
    methods.append(functools.partial(ut.logNormal,strH5ad_raw,strH5ad,config['normScale']))
    
    # apply SpaGCN
    methods.append(functools.partial(spa.run,strConfig,strSep,config))
    
    # apply BayesSpace
    methods.append(functools.partial(bay.run,strConfig,config))

    # apply cell2location
    methods.append(functools.partial(c2l.run,strConfig,strH5ad_raw))
    
    # Run all methods
    print("\n\t===== Running methods =====")
    strFinals=cu.submit_funs(methods,len(methods) if config['parallel'] else 1)

    # merge all
    print("\n\t===== Merging methods =====")
    spa.merge(strH5ad,strFinals)
    bay.merge(strH5ad,strFinals)
    c2l.merge(strH5ad,strFinals)
    print("\n\n=== visium process is completed! ===")
if __name__ == "__main__":
    if len(sys.argv)==1:
        print("=== Welcome to 'visium' from SpaceRequest! ===")
        print("\tPlease provide either a path to a folder or a config file.")
        print("\tAn empty config file will be created if a path to a folder is provided.")
        ut.msgError(ut.getSys().get('powerMsg'))
    elif len(sys.argv)>1 and sys.argv[1]=="saveRaw":
        config,sInfo = ut.getConfig(sys.argv[2])
        saveRaw(sInfo,config,sys.argv[3],sys.argv[4])
    else:
        main()
        ut.msgError(ut.getSys().get('powerMsg'))
