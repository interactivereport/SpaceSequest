#The implementation from https://github.com/interactivereport/Mouse_Cuprizone_Spatial/blob/main/SpaGCN_multisample.py
import os,sys,warnings,logging,yaml,pickle,random,torch,math,functools,io,colorsys,re
import numpy as np
import scanpy as sc
import paste as pst
import utility as ut
import cmdUtility as cu
import seaborn as sns
from copy import deepcopy
import SpaGCN as spg
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.backends.backend_pdf import PdfPages

print=functools.partial(print, flush=True)
warnings.simplefilter(action='ignore', category=FutureWarning)
logging.disable(level=logging.INFO)
predKey = "SpaGCN_pred"
#sys.stdout = sys.__stdout__
strPipePath=os.path.dirname(os.path.realpath(__file__))
mKey="SpaGCN"

def preprocess(slices,config):
    print("*** filter ***")
    # SpaGCN dosen't return normalized expression value but only return the clustering information (one column in obs)
    for one in slices.values():
        sc.pp.filter_genes(one,min_counts=config.get('UMI_count_gene'))
        sc.pp.filter_genes(one,min_cells=config.get('min_cells'))
        sc.pp.filter_cells(one,min_counts=config.get('UMI_count_cell'))
        sc.pp.filter_cells(one,min_genes=config.get('min_features'))
        sc.pp.filter_cells(one,max_counts=config.get('highCount_cutoff'))
        sc.pp.filter_cells(one,max_genes=config.get('highGene_cutoff'))
    return

def plotAlign(slices,strWK):
    print("*** plot align ***")
    sKey = list(slices.keys())
    center_color='orange'
    palette = dict(zip(sKey,sns.color_palette(None, len(sKey))))
    plt.figure(figsize=(7,7))
    for i in slices.keys():
        slice_i = slices[i].copy()
        # reversing x-y coordinate
        slice_i.obsm['spatial'] = slice_i.obsm['spatial'][:, [1, 0]]
        pst.plot_slice(slice_i,palette[i],s=10)
    handle_list = [mpatches.Patch(color=palette[i], label=i) for i in sKey]
    plt.legend(handles=handle_list)
    plt.gca().invert_yaxis()
    plt.axis('off')
    plt.savefig(os.path.join(strWK,"align.pdf"),bbox_inches="tight")
    plt.close()

def generate_rgb_colors(n):
    colors = []
    for i in range(n):
        hue = i / n
        rgb = colorsys.hsv_to_rgb(hue, 1, 1)
        colors.append(rgb)
    return colors

def plotSpaGCN(adata_all,strWK):
    print("*** plot SpaGCN ***")
    SpaGCN_list={}
    batchN=adata_all.obs.batch.nunique()
    nc = math.ceil(math.sqrt(batchN))
    nr = round(math.sqrt(batchN))
    with PdfPages(os.path.join(strWK,"SpaGCN.pdf")) as pdf:
        for oneCluster in [a for a in adata_all.obs.columns if a.startswith(predKey)]:
            cm = dict(zip(adata_all.obs[oneCluster].unique().tolist(), generate_rgb_colors(adata_all.obs[oneCluster].nunique())))
            fig, axs = plt.subplots(nr, nc, figsize=(nc*5+2,nr*5))
            axs = axs.reshape(-1)
            n=0
            legends=[]
            for oneKey in adata_all.obs.batch.unique():
                RJdata = adata_all[adata_all.obs.batch==oneKey,:].copy()
                ax  = sc.pl.scatter(RJdata,
                                    alpha=1,
                                    x="pixel_x_align",
                                    y="pixel_y_align",
                                    color=oneCluster,
                                    palette=cm,
                                    show=False,size=100000/RJdata.shape[0],
                                    ax = axs[n])  
                ax.set_aspect('equal', 'box')
                ax.axes.invert_yaxis()
                ax.axes.set_xlim([-1000, 1000])
                ax.axes.set_ylim([-1000, 1000])
                ax.set(title=oneKey)
                legends.append(ax.get_legend_handles_labels())
                ax.get_legend().remove()
                n +=1
            while n<len(axs):
                axs[n].axis('off')
                n +=1
            handles, labels = [sum(one, []) for one in zip(*legends)]
            ix=np.unique(labels,return_index=True)[1]
            fig.legend([handles[i] for i in ix], [labels[i] for i in ix],markerscale=2)
            pdf.savefig(bbox_inches="tight")
            #plt.savefig(os.path.join(strWK,"SpaGCN.pdf"))#,bbox_inches="tight"
            plt.close()

def alignSlices(D_dict,strWK):
    print("*** slice align ***")
    strAlign = os.path.join(strWK,"align.pkl")
    if os.path.isfile(strAlign):
        print("Using the previous align result %s\n\tIf a new alignment is needed, please remove/rename the above file!"%strAlign)
        center,new_slices_all=ut.readPkl(strAlign)
    else:
        slices = list(D_dict.values()) # this is not a copy, changes in slices will change D_dict
        initial_slice = slices[0].copy()
        lmbda = len(slices)*[1/len(slices)]
        center_slice,pis = pst.center_align(initial_slice,slices,lmbda)
        #W = center_slice.uns['paste_W']
        #H = center_slice.uns['paste_H']
        center, new_slices_all = pst.stack_slices_center(center_slice, slices, pis)
        ut.writePkl([center, new_slices_all],strAlign)
    dKey=list(D_dict.keys())
    for i in range(len(dKey)):
        D_dict[dKey[i]] = new_slices_all[i]
    plotAlign(D_dict,strWK)     
    return center

def find_resolution(n_clusters,adata_list,adj_list,l_list):
    obtained_clusters = -1
    iteration = 0
    resolutions = [0., 1.]
    oneContinue = 0
    while obtained_clusters != n_clusters and iteration < 50:
        current_res = sum(resolutions)/2
        clf=spg.multiSpaGCN()
        # setting seeds! Need to explicitly set it for each source of randomness. 
        r_seed=t_seed=n_seed=random.randint(1000,9999)
        random.seed(r_seed)
        torch.manual_seed(t_seed)
        np.random.seed(n_seed)
        clf.train(adata_list,adj_list,l_list,init_spa=True,init="louvain",res=current_res, tol=5e-3, lr=0.05, max_epochs=200)
        y_pred, prob = clf.predict()
        labels = y_pred.astype('str')
        obtained_clusters = len(np.unique(labels))
        if obtained_clusters < n_clusters:
            resolutions[0] = current_res
        else:
            resolutions[1] = current_res        
        iteration = iteration + 1
        print("iter %d: res=%.5f found %d clusters\n\n"%(iteration,current_res,obtained_clusters))
        # following due to sometimes it stacked at reporting 1 cluster, (not sure the reason) restart seems help
        if obtained_clusters==1:
            oneContinue +=1
        else:
            oneContinue=0
        if oneContinue>3:
            return None,None
    return clf,current_res

def find_resolution_multispagcn(n_clusters,adata_list,adj_list,l_list,strWK):
    adata_all = None
    res_all =[]
    for clusterN in n_clusters:
        print("*** cluster prediction for spg_clusterN=%d ***"%clusterN)
        strF = os.path.join(strWK,"res_clusterN%d.pkl"%clusterN)
        if os.path.isfile(strF):
            print("Using the previous clustering result %s\n\tIf a new clustering is needed, please remove/rename the above file!"%strF)
            clf,current_res=ut.readPkl(strF)
        else:
            nTry = 0
            while nTry<5:
                nTry += 1
                clf, current_res=find_resolution(clusterN,adata_list,adj_list,l_list)
                if not clf is None:
                    break
                print("\tStuck on the only one cluster: try again at %d times"%nTry)
            if clf is None:
                ut.msgError("SpaGCN failed with clustering %d times!"%nTry)
            ut.writePkl([clf, current_res],strF)
        if adata_all is None:
            adata_all = clf.adata_all
        y_pred, prob = clf.predict()
        adata_all.obs["%s%d"%(predKey,clusterN)] = ["C%d"%i for i in y_pred]
        res_all.append(current_res)
    return adata_all, res_all

def formatOutput(adata_all,strOut):
    plotSpaGCN(adata_all,os.path.dirname(strOut))
    ut.writePkl(adata_all.obs.copy(),strOut)
    return

def process(D_dict,p):
    print("*** SpaGCN spg_p=%.2f ***"%p)
    l_dict = {}#[None]*len(D_dict)
    adj_dict = {}#[None]*len(DlD_dictist)
    text_trap = io.StringIO()
    for i in D_dict.keys():
        print("\t",i)
        adata = D_dict[i]
        sc.pp.normalize_per_cell(adata, min_counts=0)
        sc.pp.log1p(adata)
        sc.pp.scale(adata, max_value=10)
        adata.obs['batch'] = i
        adata.obs['pixel_x_align'] = adata.obsm['spatial'][:,0].tolist()
        adata.obs['pixel_y_align'] = adata.obsm['spatial'][:,1].tolist()
        x = adata.obs['pixel_x_align']
        y = adata.obs['pixel_y_align']
        # Running SpaGCN to obtain the spatial domains
        adj=spg.calculate_adj_matrix(x=x,y=y,histology=False)
        adj_dict[i] = adj        
        sys.stdout = text_trap
        l=spg.search_l(p, adj, start=0.01, end=1000, tol=0.01, max_run=100)
        sys.stdout = sys.__stdout__
        l_dict[i] = l
    return adj_dict,l_dict

def SpaGCN(strConfig,strRaw,strOut):
    with open(strConfig,"r") as f:
        config = yaml.safe_load(f)
    with open(strRaw,"rb") as f:
        D_dict = pickle.load(f)
    os.makedirs(os.path.dirname(strOut),exist_ok=True)
    preprocess(D_dict,config)
    center = alignSlices(D_dict,os.path.dirname(strOut))
    adj_dict,l_dict=process(D_dict,config['spg_p'])
    adata_all,res_all = find_resolution_multispagcn(config['clusterN'],
                                          list(D_dict.values()),
                                          list(adj_dict.values()),
                                          list(l_dict.values()),
                                          os.path.dirname(strOut))
    return formatOutput(adata_all,strOut)

def run(strConfig,strRaw,config):
    strOut=os.path.join(config['output'],mKey,config['prj_name']+".pkl")
    if os.path.isfile(strOut):
        print("\n\nUsing previous *SpaGCN* results: %s\n\tIf a new process is wanted, please rename/remove the above file"%strOut)
        return {mKey:strOut}
    else:
        cmd = {}
        cmd[mKey] = "python -u %s/SpaGCN_run.py %s %s %s"%(strPipePath,strConfig,strRaw,strOut)
        cu.submit_cmd(cmd,config)
    return {mKey:strOut}

def merge(strH5ad,allRes):
    if not mKey in allRes.keys():
        return
    strPkl=allRes[mKey]
    print("merging %s: %s"%(mKey,strPkl))
    if not os.path.isfile(strPkl):
        print("\tSkip: the above file is missing!")
        return
    obs = ut.readPkl(strPkl)
    obs.index = [re.sub(obs["dataset_batch"][i]+"$",obs["batch"][i],obs.index[i]) for i in range(obs.shape[0])]
    D = sc.read_h5ad(strH5ad)#,backed="r+"
    selCol=~obs.columns.isin(D.obs.columns)
    if (~selCol).sum()>0:
        print("\tSkip exists:",",".join(obs.columns[~selCol]))
    if selCol.sum()>0:
        D.obs = D.obs.merge(obs.loc[:,selCol],"left",left_index=True,right_index=True)
        modCol = D.obs.columns.isin(obs.columns[selCol])
        D.obs.loc[:,modCol]=D.obs.loc[:,modCol].apply(lambda x: ut.fillNA(x,'Missing'))
        D.write(strH5ad)

def main():
    if len(sys.argv)==3:
        mergeSpaGCN(sys.argv[1],sys.argv[2])
        return()
    if len(sys.argv)<4:
        ut.msgError("3 arguments are required in SpaGCN.")
    SpaGCN(sys.argv[1],sys.argv[2],sys.argv[3])

if __name__ == "__main__":
    main()
