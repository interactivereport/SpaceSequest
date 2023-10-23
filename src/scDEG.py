import os, time, sys, re, yaml, sqlite3, json, math
import utility as ut
import cmdUtility as cU
import pandas as pd

strPipePath=os.path.dirname(os.path.realpath(__file__))

def MsgHelp():
    print("\nspDEG /path/to/a/folder === or === spDEG /path/to/a/config/file\n")
    ut.msgError("An empty config file will be generated automatically when a folder is provided")
  
def initProject(strInput):
    strInput=os.path.realpath(strInput)
    print("*****\nCreating an empty config file and a empty comparison definition file in\n\t",strInput)
    print("\t Please Fill in the required information in the config (DEG_config.yml) and comparison (DEG_info.csv) files")
    # save empty DE csv table
    strDEG = os.path.join(strInput,"DEG_info.csv")
    with open(strDEG,"w") as f:
        f.write("comparisonName,sample,cluster,group,alt,ref,covars[+ separated],method[default NEBULA],model[default HL]\n")
    DEGconfig = [re.sub("initDEG",os.path.join(strInput,"DEG_info.csv"),
        re.sub("initJob","j%d"%np.random.choice(range(100),1)[0],
        re.sub('strOutput',strInput,one))) for one in ut.readLines(os.path.join(strPipePath,'deg.config'))]
    strConfig = os.path.join(strInput,"DEG_config.yml")
    with open(strConfig,"w") as f:
        f.writelines(DEGconfig)
    ut.msgError()

def getConfig(strConfig):
    if not os.path.isfile(strConfig):
        ut.msgError("Missing config: %s"%strConfig)    
    with open(strConfig,"r") as f:
        config = yaml.safe_load(f)
    dInfo = checkConfig(config)
    return config, dInfo

def checkConfig(config):
    if config.get('UMI') is None or not os.path.isfile(config.get('UMI')):
        ut.msgError("Skip spDEG: UMI file required to exist")
    if config.get('meta') is None or not os.path.isfile(config.get('meta')):
        ut.msgError("Skip spDEG: meta file required to exist")

    if config['DEG_desp'] is None or not os.path.isfile(config['DEG_desp']):
        ut.msgError("Skip spDEG: Missing DEG description file!")
    D = pd.read_csv(config['DEG_desp'],header=0)
    if D.shape[0]==0:
        ut.msgError("Skip spDEG: Empty DEG description file!")
    return D

def runDEG(strConfig):
    config, dInfo = getConfig(strConfig)
    prefix = os.path.join(config['output'],config['DBname'])
    if os.path.isfile("%s.db"%prefix) and not config["newProcess"]:
        ut.msgError("Skip spDEG: db file exists: %s.db"%prefix)

    cmd = "Rscript %s/scRNAseq_DE.R %s"%(strPipePath,strConfig)
    msg = cU.run_cmd(cmd).stdout.decode("utf-8")
    #msg="spDEG task creation completed"
    if "spDEG task creation completed" in msg:
        with open("%s_spDEG.cmd.json"%prefix,"r") as f:
            spDEGtask = json.load(f)
        if not config.get('memory') is None:
            memG=int(re.sub("G$","",config.get('memory')))
        else:
            memG = math.ceil(os.path.getsize(config.get('UMI'))*50/1e9)
        cU.submit_cmd(spDEGtask,config,math.ceil(memG/16),memG)
        formatDEG(prefix)

def formatDEG(prefix):
    print("=== Formating spDEG results to create the db file ===")
    with open("%s_spDEG.cmd.json"%prefix,"r") as f:
        DEGcmds = json.load(f)
    DEGpaths = list(set([one.split(";")[0].replace("cd ","") for k,one in DEGcmds.items()]))
    csv = []
    for onePath in DEGpaths:
        for f in os.listdir(onePath):
            strCSV = os.path.join(onePath,f)
            if not os.path.isfile(strCSV) or not f.endswith("csv"):
                continue
            print("\tprocessing: ",f)
            tab = pd.read_csv(strCSV).iloc[:,0:4]
            tab.columns = ["gene","log2fc","pval","qval"]
            tags = f[:-4].split("__")
            tab["contrast"] = tags[0]
            tab["tags"] = ";".join(tags[1:]+[os.path.basename(onePath)])
            csv += [tab]
    if len(csv)==0:
        return
    data = pd.concat(csv,ignore_index=True)
    data = data.dropna()
    D = data[["log2fc","pval","qval"]]
    D.index = pd.MultiIndex.from_frame(data[["gene","contrast","tags"]])
    conn = sqlite3.connect('%s.db'%prefix)
    D.to_sql("DEG",conn,if_exists="replace")
    conn.close()
def main():
    if len(sys.argv)<2:
        MsgHelp()
    strPath = sys.argv[1]
    if os.path.isdir(strPath):
        initProject(strPath)
    elif os.path.isfile(strPath):
        runDEG(os.path.realpath(strPath))
    else:
        print("The config file is required, and %s doesn't exist!"%strPath)

if __name__ == "__main__":
    start_time = time.time()
    main()
    print("--- total time passed %s seconds ---" % (time.time() - start_time))
