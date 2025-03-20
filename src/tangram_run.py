import sys, os, time, random, warnings, logging, re #, skimage
warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.simplefilter(action='ignore', category=UserWarning)
logging.disable(level=logging.INFO)

import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

import cmdUtility as cu
import utility as ut

mKey="tangram"
cellN='tg_cell_count'
tg_pred = 'tangram_ct_pred'
tg_count = 'tangram_ct_count'
tg_cID = 'tangram_cid_count'
tg_suffix = "_processed.pkl"
tmp_dir= 'tmp'
capture_id = 'capture_id'
strPipePath=os.path.dirname(os.path.realpath(__file__))

def savePDF(plot,pdf):
    pdf.savefig(bbox_inches="tight")
    plt.close()
def getRef(strRef,cluster_label,pair_sample=None):
    print("\tget sc/sn data")
    D = ad.read_h5ad(strRef,backed="r")
    if pair_sample is not None:
        print('\t\tsubset paired sample:',pair_sample)
        k = list(pair_sample.keys())[0]
        D = D[D.obs[k]==pair_sample[k]]
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
        ad_ref = ad.AnnData(D.raw.X.value.copy(),obs=D.obs.copy(),var=gInfo,dtype='int32')
    else:
        print("\traw is NOT available. X is assumed to contain UMI and used to build the cell2location model")
        if 'value' in dir(D.X):
            ad_ref = ad.AnnData(D.X.value.copy(),obs=D.obs.copy(),var=gInfo,dtype='int32')
        else:
            ad_ref = ad.AnnData(D.X.copy(),obs=D.obs.copy(),var=gInfo,dtype='int32')
    del D
    cluster_freq = ad_ref.obs[cluster_label].value_counts()
    cluster_keep = cluster_freq.index[cluster_freq>=2]
    if len(cluster_freq)!=len(cluster_keep):
        cluster_drop = cluster_freq.index[cluster_freq<2]
        print("\tFiltering out the following label due to lack of cells\n\t\t",', '.join(cluster_drop))
        ad_ref = ad_ref[ad_ref.obs[cluster_label].isin(cluster_keep)]
    print("\t\tcell number:",ad_ref.shape[1])
    return ad_ref
def getQuery(strQuery,sID,pdf,rmG=None,strImg=None):
    print("\tget visium data")
    ad_qur = ut.sr_extract(strQuery,sID)
    if not rmG is None:
        gDel = np.full(ad_qur.shape[1],False)
        for g in rmG:
            gDel |= ad_qur.var_names.str.startswith(g)
        if gDel.sum()>0:
            print("\t\tFiltering",gDel.sum(),"genes starts with",', '.join(rmG))
            ad_qur._inplace_subset_var(np.invert(gDel))
    return ad_qur
def qcQueryImg(adata,pdf,strImg=None,image_sf=1,seg_method=None):
    import squidpy as sq
    print("\tProcess visium image")
    sID = list(adata.uns['spatial'].keys())[0]
    img = sq.im.ImageContainer(adata.uns['spatial'][sID]['images']['hires'], layer="image")
    fig, axs = plt.subplots(1, 1, figsize=(6, 6))
    savePDF(axs.imshow(img["image"][:,: , 0, :3],interpolation="none"),pdf)
    sf = adata.uns['spatial'][sID]['scalefactors']['tissue_hires_scalef']
    savePDF(sc.pl.spatial(adata, alpha=0.7, frameon=False, show=False,
        scale_factor=sf),
        pdf)
    if strImg is not None and os.path.isfile(strImg):
        print("\t\treplacing with high resolution image: ",strImg)
        img = np.array(Image.open(strImg))
        adata.uns['spatial'][sID]['images']['hires']=img
        savePDF(sc.pl.spatial(adata, alpha=0.7, frameon=False, show=False,scale_factor=1),pdf)
        adata.uns['spatial'][sID]['scalefactors']['tissue_hires_scalef']=image_sf
        sf = adata.uns['spatial'][sID]['scalefactors']['tissue_hires_scalef']
        img = sq.im.ImageContainer(img, layer="image")
    sq.im.process(img=img, layer="image", method="smooth")
    savePDF(img.show("image", channelwise=True),pdf)
    print("\t\tsegmenting image")
    if seg_method is None or seg_method=='watershed':
        sq.im.segment(img=img, layer="image_smooth", method="watershed", geq=False)
        seg_name = 'segmented_watershed'
    elif seg_method=='cellpose':
        sq.im.segment(img=img,layer="image",channel=None,method=cellpose_he,channel_cellpose=1)
        seg_name = 'segmented_custom'
    else:
        ut.msgError("Unknown segment method: ",seg_method)
    print("\t\tsegmenting image completed with Number: %d"%len(np.unique(img[seg_name])))
    #seg_binary=np.array(img[seg_name])
    # the following produce an error message: cannot reindex or align along dimension 'channels' because of conflicting dimension sizes: {1, 3}
    #img.add_img(np.multiply((seg_binary>0), 1), layer="segmented_binary")
    adata.obsm['spatial']=adata.obsm['spatial']*sf
    print('\t\tcalculate segmentation features')
    features_kwargs = {"segmentation": {
        "label_layer": seg_name,
        "props": ["label", "centroid"],
        "channels": [0,1,2]
    }}
    sq.im.calculate_image_features(adata,
        img,
        layer="image",
        key_added="image_features",
        features_kwargs=features_kwargs,
        features="segmentation",
        mask_circle=True)
    adata.obs[cellN] = adata.obsm["image_features"]["segmentation_label"]
    savePDF(sc.pl.spatial(adata, color=cellN,frameon=False,scale_factor=1,show=False),pdf)
def cellpose_he(img, min_size=15, flow_threshold=1, channel_cellpose=0):
    model = models.Cellpose(model_type='nuclei', gpu=True)
    res, _, _, _ = model.eval(
        img,
        channels=[channel_cellpose, 0],
     #   diameter=None,
        diameter=5, #can customize the diameter if the estimated one does not fit the image
        min_size=min_size,
        invert=True,
        flow_threshold=flow_threshold,
    )
    return res
def normRef(ad_ref,ad_qur,cluster_label):
    print("\t\tNormalize data")
    sc.pp.normalize_total(ad_qur)
    sc.pp.normalize_total(ad_ref)
    ad_ref.layers['log_norm']=sc.pp.log1p(ad_ref.X,copy=True)
    sc.tl.rank_genes_groups(ad_ref, groupby=cluster_label,layer='log_norm',use_raw=False)
    markers_df_log = pd.DataFrame(ad_ref.uns["rank_genes_groups"]["names"]).iloc[0:100, :]
    markers_log = list(np.unique(markers_df_log.melt().value.values))
    return markers_log
def applyTangram(ad_ref,ad_qur,cluster_label,pdf,parallel='cpu'):
    import tangram as tg
    print("\tApply Tangram")
    markers = normRef(ad_ref,ad_qur,cluster_label)
    tg.pp_adatas(ad_ref,ad_qur,genes=markers)
    print("\t\tmap cell to space")
    ad_map = tg.map_cells_to_space(
        ad_ref,
        ad_qur,
        mode="constrained",
        target_count=ad_qur.obs[cellN].sum(),
        density_prior=np.array(ad_qur.obs[cellN])/ad_qur.obs[cellN].sum(),
        num_epochs=1000,
        device=parallel)
    #return ad_qur, ad_ref, ad_map
    print("\t\tproject cell annotations")
    tg.project_cell_annotations(ad_map,ad_qur,annotation=cluster_label)
    annotation_list = list(ad_ref.obs[cluster_label].unique())
    savePDF(tg.plot_cell_annotation_sc(ad_qur,annotation_list,perc=0.008,scale_factor=1),pdf)
    savePDF(tg.plot_training_scores(ad_map, bins=20, alpha=.5),pdf)
    ad_ge = tg.project_genes(adata_map=ad_map, adata_sc=ad_ref)
    df_all_genes = tg.compare_spatial_geneexp(ad_ge, ad_qur, ad_ref)
    savePDF(tg.plot_auc(df_all_genes),pdf)
    tg.create_segment_cell_df(ad_qur)
    print("\t\tcount cell annoations:",cluster_label)
    tg.count_cell_annotations(
        ad_map,
        ad_ref,
        ad_qur,
        annotation=cluster_label)
    adata_segment = tg.deconvolve_cell_annotations(ad_qur)
    fig, ax = plt.subplots(1, 1, figsize=(8, 8))
    sc.pl.spatial(
        adata_segment,
        color="cluster",
        size=0.8,
        show=False,
        frameon=False,
        alpha_img=0.2,
        legend_fontsize=20,
        ax=ax,
        scale_factor=1
    )
    savePDF(None,pdf)
    ct_count_tmp = ad_qur.obsm[tg_count].copy()# keep it in temp, it will be replaced by the following cell id 
    ad_ref.obs['cell_id']=ad_ref.obs.index
    print("\t\tcount cell annoations: cell_id")
    tg.count_cell_annotations(
        ad_map,
        ad_ref,
        ad_qur,
        annotation="cell_id",
    )
    ad_qur.obsm[tg_cID] = ad_qur.obsm[tg_count].copy()
    ad_qur.obsm[tg_count]=ct_count_tmp
    return ad_qur, ad_ref, adata_segment
def mapRef(ctMap,segLoc,ref_label):
    print("\t\tmaping scRef loc")
    segLoc['spot'] = segLoc.apply(lambda x:'_'.join(x['centroids'].split("_")[:-1]),axis=1)
    ref_loc = pd.DataFrame({'y':np.nan,'x':np.nan,'spot':""},index=ctMap.columns[4:])
    ctMap['spot'] = list(ctMap.index)
    ct_num = ctMap[ref_loc.index].sum()
    for cID in [one for one in ct_num.index if ct_num[one]>0]:
        spot = list(ctMap.index[ctMap[cID]>0])
        ix = np.logical_and(segLoc.cluster==ref_label[cID],segLoc.spot.isin(spot))
        if ix.sum()>0:
            ref_loc.loc[cID,['y','x','spot']]=segLoc[ix][['y','x','spot']].iloc[random.choice(range(ix.sum())),]
        else:
            ref_loc.loc[cID,['y','x','spot']]=ctMap.loc[ctMap[cID]>0,['y','x','spot']].iloc[0,]
    return ref_loc
def returnOne(ad_qur,ad_ref,adata_segment,cluster_label,prefix=None):
    print("\tsaving")
    adata_segment.uns['tangram_cell_segmentation'].centroids = adata_segment.uns['tangram_cell_segmentation'].centroids.astype(str)
    sc_loc = mapRef(ad_qur.obsm[tg_cID].copy(),adata_segment.obs.copy(),ad_ref.obs[cluster_label])
    qur_cluster=adata_segment.obs
    qur_cluster.index = list(adata_segment.obs.centroids)
    if prefix is not None and os.path.exists(os.path.dirname(prefix)):
        adata_segment.write(prefix+"_segment.h5ad")
        ut.writePkl([ad_qur.obsm[tg_pred],ad_qur.obsm[tg_count],qur_cluster,sc_loc],prefix+tg_suffix)
    else:
        return ad_qur.obsm['tg_pred'], ad_qur.obsm['tg_count'], qur_cluster, sc_loc
def createSegment(tg_seg):
    D=ad.AnnData(pd.DataFrame({"FakeG%d"%i:0 for i in range(2)},index=tg_seg.index))
    D.obs=tg_seg
    for one in tg_seg[capture_id].unique():
        D.obsm["X_%s"%one] = tg_seg.apply(lambda x: x[['y','x']] if x[capture_id]==one else pd.Series({'y':0,'x':0}),axis=1).to_numpy()
    return(D)
def createRefLoc(strRef,sc_loc,strF):
    D = ad.read_h5ad(strRef)
    D.obs['tg_capture_id']=""
    D.obs['tg_spot']=""
    for one in sc_loc.keys():
        D.obs.loc[sc_loc[one].index,'tg_capture_id'] = one
        D.obs.loc[sc_loc[one].index,'tg_spot'] = sc_loc[one]['spot']
        D.obsm["X_%s"%one] = pd.DataFrame(index=D.obs_names).merge(sc_loc[one][['y','x']],'left',left_index=True, right_index=True).to_numpy()
    D.write(strF)

#parallel: cpu or cuda:0
def processOne(strRef,cluster_label,strQuery,sID,strOut,pair_sample=None,rmG=None,strImg=None,image_sf=1,seg_method=None,parallel='cpu'):
    print("---",sID,'---')
    with PdfPages(os.path.join(strOut,sID+".pdf")) as pdf:
        ad_ref = getRef(strRef,cluster_label,pair_sample)
        ad_qur = getQuery(strQuery,sID,pdf,rmG,strImg)
        #ad_ref._inplace_subset_var(ad_ref.var_names.isin(ad_qur.var_names))
        #ad_qur._inplace_subset_var(ad_qur.var_names.isin(ad_ref.var_names))
        ad_ref = ad_ref[:,ad_ref.var_names.isin(ad_qur.var_names)]
        ad_qur = ad_qur[:,ad_qur.var_names.isin(ad_ref.var_names)]
        qcQueryImg(ad_qur,pdf,strImg,image_sf,seg_method)
        ad_qur, ad_ref, adata_segment = applyTangram(ad_ref,ad_qur,cluster_label,pdf,parallel)
    return returnOne(ad_qur, ad_ref, adata_segment,cluster_label,os.path.join(strOut,sID))
def oneRun(strConfig,strH5ad,sID):
    config,sInfo = ut.getConfig(strConfig)
    strOut = os.path.join(config['output'],mKey,tmp_dir)
    if os.path.isfile(os.path.join(strOut,sID+tg_suffix)):
        print("\t Using previous %s results: %s"%(sID,os.path.join(strOut,sID+tg_suffix)))
        return
    pair_sample = None
    if config.get('tg_matchColumn') is not None and sInfo.loc[sID,config['tg_matchColumn']] is not None and len(str(sInfo.loc[sID,config['tg_matchColumn']]))>0:
        pair_sample={config["matchColumn"]:sInfo.loc[sID,config["matchColumn"]]}
    strImg = None
    if config.get('tg_image_column') is not None and len(sInfo.loc[sID,config['tg_image_column']])>0:
        strImg = sInfo.loc[sID,config['tg_image_column']]
    if config['parallel']=="slurm":
        sysConfig = ut.getSys()
        config['gpu'] = False if sysConfig.get("use_gpu") is None else sysConfig.get("use_gpu")
    processOne(config['tg_scH5ad'],config['tg_annotation_obs'],strH5ad,sID,strOut,
        pair_sample=pair_sample,rmG=config['tg_rmGeneStart'],strImg=strImg,seg_method=config['tg_segment_method'],
        parallel='cuda:0' if config['gpu'] else 'cpu')
def mergeInd(strConfig,strPkl):
    print("--- merge individual ---")
    config,sInfo = ut.getConfig(strConfig)
    strOut = os.path.join(config['output'],mKey,tmp_dir)
    tg_ct_pred=[]
    tg_ct_count=[]
    tg_segment=[]
    sc_loc={}
    for one in sInfo[config['sample_name']]:
        print("\textracting",one)
        strOne = os.path.join(strOut,one+tg_suffix)
        if not os.path.isfile(strOne):
            print("\t\tSkip:",one,"missing",strOne)
        one_pred, one_count, one_segment, one_loc = ut.readPkl(strOne)
        tg_ct_pred.append(one_pred)
        tg_ct_count.append(one_count)
        one_segment[capture_id] = one
        tg_segment.append(one_segment)
        sc_loc[one]=one_loc
    print("\tmerging")
    tg_ct_pred = pd.concat(tg_ct_pred)
    tg_ct_count = pd.concat(tg_ct_count)
    print("\tsaving")
    createSegment(pd.concat(tg_segment)).write(re.sub(".pkl$","_segment.h5ad",strPkl))
    createRefLoc(config['tg_scH5ad'],sc_loc,re.sub(".pkl$","_sc_location.h5ad",strPkl))
    ut.writePkl([tg_ct_pred,tg_ct_count],strPkl)
def check(config,sInfo):
    if config.get("tg_scH5ad") is None or config.get("tg_annotation_obs") is None:
        print("___ Skip tangram: no tg_scH5ad (or tg_annotation_obs) is provided! ___")
        return False
    if not os.path.isfile(config["tg_scH5ad"]):
        print("___ Skip tangram: tg_scH5ad (%s) does not exist! ___"%config["tg_scH5ad"])
        return False
    D = ad.read_h5ad(config["tg_scH5ad"],backed="r")
    if not config["tg_annotation_obs"] in D.obs.columns:
        print("___ Skip tangram: tg_annotation_obs (%s) does not exist in scH5ad! ___"%config["tg_annotation_obs"])
        return False
    if not config.get("tg_matchColumn") is None:
        if not config["tg_matchColumn"] in D.obs.columns:
            print("___ Skip tangram: tg_matchColumn (%s) does not exist in scH5ad! ___"%config["tg_matchColumn"])
            return False
        if not config["tg_matchColumn"] in sInfo.columns:
            print("___ Skip tangram: tg_matchColumn (%s) does not exist in sample meta! ___"%config["tg_matchColumn"])
            return False
        missName = [str(one) for one in sInfo[config["tg_matchColumn"]] if len(str(one))>0 and not D.obs[config["tg_matchColumn"]].eq(one).any() and not D.obs[config["tg_matchColumn"]].eq(str(one)).any()]
        if len(missName)>0:
            print("___ Skip tangram: missing tg_matchColumn values (%s) in scH5ad! ___"%', '.join(missName))
            return False
    if not config.get('tg_image_column') is None:
        missImg = [one for one in sInfo[config["tg_image_column"]] if len(one)>0 and not os.path.isfile(one)]
        if len(missImg)>0:
            print("___ Skip tangram: missing tg_image_column files: ___\n%s"%'\n'.join(missImg))
            return False
    return True
    
def run(strConfig,strH5ad):
    config,sInfo = ut.getConfig(strConfig)
    if not check(config,sInfo):
        return
    strOut=os.path.join(config['output'],mKey,config['prj_name']+".pkl")
    if os.path.isfile(strOut):
        print("\n\nUsing previous *Tangram* results: %s\n\tIf a new process is wanted, please rename/remove the above file"%strOut)
        return {mKey:strOut}
    os.makedirs(os.path.join(os.path.dirname(strOut),tmp_dir),exist_ok=True)

    if config['parallel']=="slurm":
        sysConfig = ut.getSys()
        config['gpu'] = False if sysConfig.get("use_gpu") is None else sysConfig.get("use_gpu")
    cmd = {}
    for one in sInfo[config['sample_name']]:
        cmd['%s_%s'%(mKey,one)] = "python -u %s/tangram_run.py oneRun %s %s %s"%(strPipePath,strConfig,strH5ad,one)
    cu.submit_cmd(cmd,config,condaEnv="condaEnv_C2L")

    cmd= {'%s_Extract'%mKey:"python -u %s/tangram_run.py Extract %s %s"%(strPipePath,strConfig,strOut)}
    config['gpu'] = False
    cu.submit_cmd(cmd,config)#,condaEnv="condaEnv_C2L"
    
    return {mKey:strOut}

def merge(D,allRes):
    if not mKey in allRes.keys():
        return
    strF = allRes[mKey]
    print("merging %s: %s"%(mKey,strF))
    if not os.path.isfile(strF):
        print("\tSkip, the above file is missing!")
        return
    tg_ct_pred,tg_ct_count = ut.readPkl(strF)
    #D = ad.read_h5ad(strH5ad)#,backed="r+"
    if not tg_pred in D.obsm.keys():
        print("\tmerging tg_pred into obsm")
        D.obsm[tg_pred] = pd.DataFrame(index=D.obs_names).merge(tg_ct_pred,'left',left_index=True, right_index=True)#.to_numpy()
    if not tg_count in D.obsm.keys():
        print("\tmerging tg_count into obsm")
        D.obsm[tg_count] = pd.DataFrame(index=D.obs_names).merge(tg_ct_count.select_dtypes(include=['number']),'left',left_index=True, right_index=True)#.to_numpy()
    tg_ct_pred.columns = ['tg_'+one for one in tg_ct_pred.columns]
    if tg_ct_pred.columns.isin(D.obs.columns).sum()==0:
        print("\tmerging labels into obs")
        D.obs = D.obs.merge(tg_ct_pred,'left',left_index=True, right_index=True)
        D.obs[tg_ct_pred.columns] = D.obs[tg_ct_pred.columns].fillna(0)
        ut.plotVisium(D,strOut=os.path.dirname(strF),selObs=tg_ct_pred.columns.tolist())
    #print("\tsaving")
    #D.write(strH5ad)

def main():
    task = sys.argv[1]
    if task == 'oneRun':
        oneRun(sys.argv[2],sys.argv[3],sys.argv[4])
    elif task == 'Extract':
        mergeInd(sys.argv[2],sys.argv[3])
    else:
        ut.msgError("unknown task for tangram: %s!"%task)

if __name__ == "__main__":
    main()
