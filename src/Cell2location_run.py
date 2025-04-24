import sys,os,re,warnings,logging,random#,scvi
from numba.core.errors import NumbaDeprecationWarning, NumbaPendingDeprecationWarning
warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.simplefilter(action='ignore', category=UserWarning)
warnings.simplefilter(action='ignore', category=NumbaDeprecationWarning)
logging.disable(level=logging.INFO)
import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np

import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import rcParams
rcParams['pdf.fonttype'] = 42 # enables correct plotting of text

import cmdUtility as cu
import utility as ut
#import scvi

strPipePath=os.path.dirname(os.path.realpath(__file__))
mKey = "c2l"
## took from cell2location without plt.show() to save the plot outside
def filter_genes(adata, cell_count_cutoff=15, cell_percentage_cutoff2=0.05, nonz_mean_cutoff=1.12):
    r"""Plot the gene filter given a set of cutoffs and return resulting list of genes.

    Parameters
    ----------
    adata :
        anndata object with single cell / nucleus data.
    cell_count_cutoff :
        All genes detected in less than cell_count_cutoff cells will be excluded.
    cell_percentage_cutoff2 :
        All genes detected in at least this percentage of cells will be included.
    nonz_mean_cutoff :
        genes detected in the number of cells between the above mentioned cutoffs are selected
        only when their average expression in non-zero cells is above this cutoff.

    Returns
    -------
    a list of selected var_names
    """

    adata.var["n_cells"] = np.array((adata.X > 0).sum(0)).flatten()
    adata.var["nonz_mean"] = np.array(adata.X.sum(0)).flatten() / adata.var["n_cells"]

    cell_count_cutoff = np.log10(cell_count_cutoff)
    cell_count_cutoff2 = np.log10(adata.shape[0] * cell_percentage_cutoff2)
    nonz_mean_cutoff = np.log10(nonz_mean_cutoff)

    gene_selection = (np.array(np.log10(adata.var["n_cells"]) > cell_count_cutoff2)) | (
        np.array(np.log10(adata.var["n_cells"]) > cell_count_cutoff)
        & np.array(np.log10(adata.var["nonz_mean"]) > nonz_mean_cutoff)
    )
    gene_selection = adata.var_names[gene_selection]
    adata_shape = adata[:, gene_selection].shape

    fig, ax = plt.subplots()
    ax.hist2d(
        np.log10(adata.var["nonz_mean"]),
        np.log10(adata.var["n_cells"]),
        bins=100,
        norm=mpl.colors.LogNorm(),
        range=[[0, 0.5], [1, 4.5]],
    )
    ax.axvspan(0, nonz_mean_cutoff, ymin=0.0, ymax=(cell_count_cutoff2 - 1) / 3.5, color="darkorange", alpha=0.3)
    ax.axvspan(
        nonz_mean_cutoff,
        np.max(np.log10(adata.var["nonz_mean"])),
        ymin=0.0,
        ymax=(cell_count_cutoff - 1) / 3.5,
        color="darkorange",
        alpha=0.3,
    )
    plt.vlines(nonz_mean_cutoff, cell_count_cutoff, cell_count_cutoff2, color="darkorange")
    plt.hlines(cell_count_cutoff, nonz_mean_cutoff, 1, color="darkorange")
    plt.hlines(cell_count_cutoff2, 0, nonz_mean_cutoff, color="darkorange")
    plt.xlabel("Mean non-zero expression level of gene (log)")
    plt.ylabel("Number of cells expressing gene (log)")
    plt.title(f"Gene filter: {adata_shape[0]} cells x {adata_shape[1]} genes")
    return gene_selection

def buildModel(config,strOne,selSample=None,use_gpu=False):
    from cell2location.models import RegressionModel
    print("Building model: use_gpu (%s)"%use_gpu)
    if selSample is None:# if all h5ad cells were used to build a model, all slices will use the same model
        strModel = os.path.join(os.path.dirname(strOne),"model","full_sc_model.h5ad")
    else:
        strModel = os.path.join(os.path.dirname(strOne),"model",re.sub('pkl$','h5ad',os.path.basename(strOne)))
    os.makedirs(os.path.dirname(strModel),exist_ok=True)
    if os.path.isfile(strModel):
        print("\tLoading exist model file: ",strModel)
        adata_ref=ad.read_h5ad(strModel)
    else:    
        D = ad.read_h5ad(config["scH5ad"],backed="r")
        if not config.get("matchColumn") is None and not selSample is None:
            D=D[D.obs[config["matchColumn"]]==selSample] #np.logical_and(D.obs[config["matchColumn"]]==selSample,pd.Series(random.choices([True]+[False]*49,k=D.shape[0]),index=D.obs.index))
        if D.shape[0] <10:
            print("Too few reference cells (%d) selected!"%D.shape[0])
            return None
        if not D.raw is None and not D.raw.var is None:
            gInfo = D.raw.var.copy()
        else:
            gInfo = D.var.copy()
        if 'feature_name' in gInfo.columns:
            print("\tUsing feature_name as gene names")
            gInfo.index = gInfo.feature_name.tolist()
        else:
            print("\tUsing var_names as gene names")
        print("\tUsing %d cells to build the model"%D.shape[0])
        if not D.raw is None and not D.raw.X is None:
            print("\traw is available and used to build the cell2location model")
            adata_ref = ad.AnnData(D.raw.X.value.copy(),obs=D.obs.copy(),var=gInfo,dtype='int32')
        else:
            print("\traw is NOT available. X is assumed to contain UMI and used to build the cell2location model")
            if 'value' in dir(D.X):
                adata_ref = ad.AnnData(D.X.value.copy(),obs=D.obs.copy(),var=gInfo,dtype='int32')
            else:
                adata_ref = ad.AnnData(D.X.copy(),obs=D.obs.copy(),var=gInfo,dtype='int32')
        del D
        sc.pp.filter_cells(adata_ref, min_genes=1)
        sc.pp.filter_genes(adata_ref, min_cells=1)
        with PdfPages(re.sub("h5ad$","pdf",strModel)) as pdf:
            selected = filter_genes(adata_ref, cell_count_cutoff=config['ref_cell_count_cutoff'],
                                    cell_percentage_cutoff2=config['ref_cell_percentage_cutoff2'],
                                    nonz_mean_cutoff=config['ref_nonz_mean_cutoff'])       
            pdf.savefig(bbox_inches="tight")
            plt.close()
            adata_ref = adata_ref[:, selected].copy()
            adata_ref.obs[config['annotation_obs']] =adata_ref.obs[config['annotation_obs']].astype(str)
            RegressionModel.setup_anndata(adata=adata_ref,batch_key=config['batch'],labels_key=config['annotation_obs'])
            mod = RegressionModel(adata_ref)
            mod.view_anndata_setup()
            print("\ttraining")
            mod.train(max_epochs=300,batch_size=15000,train_size=1,lr=0.002,use_gpu=use_gpu)  # , use_gpu=True     ,batch_size=2500
            # plot ELBO loss history during training, removing first 20 epochs from the plot
            fig = plt.figure()
            mod.plot_history(20)
            pdf.savefig(bbox_inches="tight")
            plt.close()
            fig = plt.figure()
            adata_ref = mod.export_posterior(adata_ref, sample_kwargs={'num_samples': 1000, 'batch_size': 15000, 'use_gpu': use_gpu}) #'batch_size': 2500
            mod.plot_QC()
            pdf.savefig(bbox_inches="tight")
            plt.close()
        # Save model
        adata_ref.write(strModel)
        #ut.writePkl(mod,adata_ref],re.sub("h5ad$","pkl",strModel))

    if 'means_per_cluster_mu_fg' in adata_ref.varm.keys():
        inf_aver = adata_ref.varm['means_per_cluster_mu_fg'][[f'means_per_cluster_mu_fg_{i}'
                                                              for i in adata_ref.uns['mod']['factor_names']]].copy()
    else:
        inf_aver = adata_ref.var[[f'means_per_cluster_mu_fg_{i}'
                                  for i in adata_ref.uns['mod']['factor_names']]].copy()
    inf_aver.columns = [mKey+"_"+_ for _ in adata_ref.uns['mod']['factor_names']]
    return inf_aver

def applyModel(adata,inf_aver,config,useGPU,strOne):
    import cell2location
    print("Applying model: use_gpu (%s)"%useGPU)
    if os.path.isfile(strOne):
        print("\tUsing previous results:",strOne)
        print("\t\tPlease remove/rename the above file if a new run is wanted!")
    else:
        intersect = adata.var_names.intersection(inf_aver.index) #np.intersect1d(adata.var_names, inf_aver.index)
        if len(intersect) < 10:
            ut.msgError("in cell2location, less than 10 genes overlap with reference!")
        adata = adata[:, intersect].copy()
        inf_aver = inf_aver.loc[intersect, :].copy()
        
        ## prepare anndata for cell2location model
        cell2location.models.Cell2location.setup_anndata(adata=adata, batch_key=ut.batchKey)
        # create model
        mod = cell2location.models.Cell2location(
            adata, cell_state_df=inf_aver,
            # the expected average cell abundance: tissue-dependent
            # hyper-prior which can be estimated from paired histology:
            N_cells_per_location=config['N_cells_per_location'],
            # hyperparameter controlling normalisation of
            # within-experiment variation in RNA detection (using default here):
            detection_alpha=config['detection_alpha']
        )
        # train model
        mod.train(max_epochs=30000,
                  # train using full data (batch_size=None)
                  batch_size=None,
                  # use all data points in training because
                  # we need to estimate cell abundance at all locations
                  train_size=1,
                  use_gpu=useGPU,
                  progress_bar_refresh_rate=0.2)
        adata = mod.export_posterior(adata, sample_kwargs={'num_samples': 1000, 'batch_size': mod.adata.n_obs, 'use_gpu': useGPU})
        adata.obs[adata.uns['mod']['factor_names']] = adata.obsm['q05_cell_abundance_w_sf']
        # Compute expected expression per cell type
        expected_dict = mod.module.model.compute_expected_per_cell_type(mod.samples["post_sample_q05"], mod.adata_manager)
        # Add to anndata layers
        for i, n in enumerate(mod.factor_names_):
            adata.layers[n] = expected_dict['mu'][i]
        adata.write(re.sub("pkl$","h5ad",strOne))
        ut.writePkl(adata.obs.copy(),strOne)
        plotCell2Location(mod,adata,strOne)
    
def plotCell2Location(mod,adata,strOne):
    print("\tPlotting the final cell2location")
    with PdfPages(re.sub("pkl$","pdf",strOne)) as pdf:
        mod.plot_history(1000)
        plt.legend(labels=['full data training'])
        pdf.savefig(bbox_inches="tight")
        plt.close()
        celltypes = adata.uns['mod']['factor_names']
        with mpl.rc_context({'axes.facecolor':  'black',
                     'figure.figsize': [4.5, 5]}):
            fig=sc.pl.spatial(adata, cmap='magma',
                          color = celltypes,
                          ncols=4, size=1.3,
                          img_key='lowres',
                          # limit color scale at 99.2% quantile of cell abundance
                          vmin=0, vmax='p99.2',
                          return_fig=True,show=False)
            pdf.savefig(bbox_inches="tight")
            plt.close()

def oneRun(strOut,strConfig,strH5ad,sName):
    print("*** cell2location: %s ***"%sName)
    strOne = os.path.join(strOut,sName+".pkl")
    if os.path.isfile(strOne):
        print("\tUsing existing result: %s\n\tIf a new run for this sample is wanted, please rename/remove the above file!"%strOne)
        return()
    config,sInfo = ut.getConfig(strConfig)
    sInfo.index = sInfo[config['sample_name']]
    useGPU = True
    if config['parallel']=="slurm":
        sysConfig = ut.getSys()
        useGPU = False if sysConfig.get("use_gpu") is None else sysConfig.get("use_gpu")

    inf_aver=buildModel(config,strOne,
        None if config["matchColumn"] is None else sInfo.loc[sName,config["matchColumn"]],
        useGPU)
    adata_vis = ut.sr_extract(strH5ad,sName)
    adata_vis = applyModel(adata_vis,inf_aver,config,useGPU,strOne)

def mergeInd(strFinal,sNamesList):
    strOut=os.path.dirname(strFinal)
    sNames=sNamesList.split(',')
    obs = []
    for sName in sNames:
        strOne = os.path.join(strOut,sName+".pkl")
        if not os.path.isfile(strOne):
            print("\tSkip %s: missing final result: %s"%(sName,strOne))
        oneObs = ut.readPkl(strOne)
        #D1 = ad.read_h5ad(strOne,backed="r")
        obs.append(oneObs[[i for i in oneObs.columns if i.startswith(mKey+"_")]].copy())
        del oneObs
    obs = pd.concat(obs)
    obs.fillna(0,inplace=True)
    ut.writePkl(obs,strFinal)

def check(config):
    if config.get("scH5ad") is None or config.get("annotation_obs") is None:
        print("___ Skip cell2location: no scH5ad (or annotation_obs) is provided! ___")
        return False
    if config.get("batch") is None:
        print("___ Skip cell2location: 'batch' is required for reference data! ___")
        return False
    if not os.path.isfile(config["scH5ad"]):
        print("___ Skip cell2location: scH5ad (%s) does not exist! ___"%config["scH5ad"])
        return False
    D = ad.read_h5ad(config["scH5ad"],backed="r")
    if not config["annotation_obs"] in D.obs.columns:
        print("___ Skip cell2location: annotation_obs (%s) does not exist in scH5ad! ___"%config["annotation_obs"])
        return False
    if not config.get("matchColumn") is None and not config["matchColumn"] in D.obs.columns:
        print("___ Skip cell2location: matchColumn (%s) does not exist in scH5ad! ___"%config["matchColumn"])
        return False
    return True

def run(strConfig,strH5ad):
    config,sInfo = ut.getConfig(strConfig)
    strOut = os.path.join(config['output'],mKey)
    strFinal = os.path.join(strOut,config['prj_name']+".pkl")
    if os.path.isfile(strFinal):
        print("\n\nUsing previous *Cell2location* results: %s\n\tIf a new process is wanted, please rename/remove the above file"%strFinal)
        return {mKey:strFinal}
    if not check(config):
        return
    os.makedirs(os.path.dirname(strOut),exist_ok=True)
    obsAnno = config.get("matchColumn")
    if config['parallel']=="slurm":
        sysConfig = ut.getSys()
        config['gpu'] = False if sysConfig.get("use_gpu") is None else sysConfig.get("use_gpu")
    strCMD={}
    for one in sInfo[config['sample_name']]:
        strCMD['C2L_%s'%one] = "python -u %s/Cell2location_run.py oneRun %s %s %s %s"%(strPipePath,strOut,strConfig,strH5ad,one)
    #print("\n".join(strCMD.values()))
    cu.submit_cmd(strCMD,config,condaEnv="condaEnv_C2L")
    
    strCMD['C2L_Extract'] = "python -u %s/Cell2location_run.py Extract %s %s"%(strPipePath,strFinal,','.join([str(_) for _ in sInfo[config['sample_name']]]))
    config['gpu'] = False
    cu.submit_cmd(strCMD,config,condaEnv="condaEnv_C2L")
    
    return {mKey:strFinal}

def merge(D,allRes):
    if not mKey in allRes.keys():
        return
    strF = allRes[mKey]
    print("merging %s: %s"%(mKey,strF))
    if not os.path.isfile(strF):
        print("\tSkip, the above file is missing!")
        return
    obs = ut.readPkl(strF)
    #D = ad.read_h5ad(strH5ad)#,backed="r+"
    selCol=~obs.columns.isin(D.obs.columns)
    if (~selCol).sum()>0:
        print("\tSkip exists:",",".join(obs.columns[~selCol]))
    if selCol.sum()>0:
        D.obs = D.obs.merge(obs.loc[:,selCol],"left",left_index=True,right_index=True)
        modCol = D.obs.columns.isin(obs.columns[selCol])
        D.obs.loc[:,modCol]=D.obs.loc[:,modCol].apply(lambda x: ut.fillNA(x,0))
        ut.plotVisium(D,strOut=os.path.dirname(strF),obs=D.obs.columns[modCol].tolist())
        #D.write(strH5ad)

def main():
    task=sys.argv[1]
    if task=="oneRun":
        if not len(sys.argv)==6:
            ut.msgError("The number of input for oneRun in cell2location is incorrect!")
        a = sys.argv[2:].copy()
        oneRun(*a)
    elif task=="Extract":
        if not len(sys.argv)==4:
            ut.msgError("The number of input for Extract in cell2location is incorrect!")
        a = sys.argv[2:].copy()
        mergeInd(*a)
    else:
        ut.msgError("unknown task for cell2location: %s!"%task)

if __name__ == "__main__":
    main()
