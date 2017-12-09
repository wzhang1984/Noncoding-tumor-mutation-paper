import numpy as np
import networkx as nx
import os
import mkl
import pandas as pd
import scipy as sp
from scipy.linalg import inv
from sklearn.decomposition import SparsePCA,NMF,PCA
from sklearn.cluster import SpectralClustering,KMeans,AgglomerativeClustering
from multiprocessing import Pool
from scipy.stats import fisher_exact,ttest_ind,mannwhitneyu,rankdata
from sklearn import metrics


def renorm(network):    
    degree = np.sum(network,axis=0)
    return network*1.0/degree  


def run_diffusion(network,rst_prob,mutation_profile,converge_rate,max_iter=50,normalize_mutations=True):
    
    P = renorm(network).T
    if normalize_mutations:
        mutation_profile = renorm(mutation_profile.T).T
    Q = mutation_profile.copy()
    
    for i in range(max_iter):
        Q_new = rst_prob * mutation_profile + (1-rst_prob) * np.dot(Q,P)
        
        delta_Q = Q - Q_new
        delta = np.sqrt(np.sum(np.square(delta_Q)))
        
        print i,'iteration: delta is',delta
        
        Q = Q_new
        if delta < converge_rate:
            break
    
    return Q


def create_ppr_matrix(network,rst_prob,output_dir):
    
    ## Create the Personalized PageRank (PPR) matrix using Scipy
    # Create "walk" matrix (normalized adjacency matrix)
    print "* Creating PPR  matrix..."
    W = renorm(network).T
    
    ## Create PPR matrix using Python
    n = network.shape[0]
    PPR = (1.-rst_prob)*inv(sp.eye(n)-rst_prob*W)
    
    os.system( 'mkdir -p ' + output_dir )
    pprfile = "{}/ppr_{:g}.npy".format(output_dir, rst_prob)
    np.save(pprfile, PPR)
    
    return PPR


def run_diffusion_PPR(PPR,mutation_profile,normalize_mutations=False):

    if normalize_mutations:
        mutation_profile = renorm(mutation_profile.T).T

    Q = np.dot(mutation_profile,PPR)

    return Q
    

def load_network(file_name,output_dir,gene_list=set()):
    
    # Load graph
    print "* Loading PPI..."
    with open(file_name) as file_handle:
        if gene_list:
            gene_pairs = [(g1, g2) for g1, g2 in [line.split()[:2] for line in file_handle.read().splitlines()] 
                          if (g1 in gene_list and g2 in gene_list)]
        else:
            gene_pairs = [(g1, g2) for g1, g2 in [line.split()[:2] for line in file_handle.read().splitlines()]]
    
    gene_set = set([g for gene_pair in gene_pairs for g in gene_pair])
    
    G = nx.Graph()
    G.add_nodes_from( gene_set ) # in case any nodes have degree zero
    G.add_edges_from( gene_pairs )
    
    print "\t- Edges:", len(G.edges())
    print "\t- Nodes:", len(G.nodes())
    
    # Remove self-loops and restrict to largest connected component
    print "* Removing self-loops, multi-edges, and restricting to",
    print "largest connected component..."
    selfLoops = [(u, v) for u, v in G.edges() if u == v]
    G.remove_edges_from( selfLoops )
    G = G.subgraph( sorted(nx.connected_components( G ), key=lambda cc: len(cc),
                           reverse=True)[0] )
    nodes = sorted(G.nodes())
    print "\t- Largest CC Edges:", len( G.edges() )
    print "\t- Largest CC Nodes:", len( G.nodes() )
    
    # Set up output directory
    print "* Saving updated node list to file..."
    os.system( 'mkdir -p ' + output_dir )
    
    # Index mapping for genes
    index_map = [ "{}\t{}".format(i, nodes[i]) for i in range(len(nodes)) ]
    with open("{}/index_genes".format(output_dir), 'w') as outfile:
        outfile.write( "\n".join(index_map) )
    gene2index = { nodes[i]: i for i in range(len(nodes)) }

    network = nx.to_numpy_matrix( G , nodelist=nodes, dtype=np.float64 )
    network = np.asarray(network)
    
    return network, gene2index


def load_mutation(file_name,output_dir,gene2index):
    
    with open(file_name) as file_handle:
        pat_gene_pairs = [(p, g) for p, g in [line.split()[:2] for line in file_handle.read().splitlines()] if g in gene2index]
    
    pats = sorted(set([p for p,g in pat_gene_pairs]))
    geneset=set([g for p,g in pat_gene_pairs])
    print "\t- Genes in adjacency matrix:", len(geneset)
    
    # Set up output directory
    print "* Saving patient list to file..."
    os.system( 'mkdir -p ' + output_dir )
    
    # Index mapping for genes
    index_map = [ "{}\t{}".format(i, pats[i]) for i in range(len(pats)) ]
    with open("{}/index_patients".format(output_dir), 'w') as outfile:
        outfile.write( "\n".join(index_map) )
    pat2index = { pats[i]: i for i in range(len(pats)) }
    
    mutation_profile = np.zeros((len(pats), len(gene2index)))
    mutation_profile[zip(*[(pat2index[p],gene2index[g]) for p,g in pat_gene_pairs])] = 1.
    return mutation_profile, pat2index


def load_mutation_from_df(df,output_dir,gene2index):
    
    # df should be a patient by gene matrix
    pats = df.index
    geneset = set(df.columns)&set(gene2index.keys())
    print "\t- Genes in adjacency matrix:", len(geneset)
    
    # Set up output directory
    print "* Saving patient list to file..."
    os.system( 'mkdir -p ' + output_dir )
    
    # Index mapping for genes
    index_map = [ "{}\t{}".format(i, pats[i]) for i in range(len(pats)) ]
    with open("{}/index_patients".format(output_dir), 'w') as outfile:
        outfile.write( "\n".join(index_map) )
    pat2index = { pats[i]: i for i in range(len(pats)) }
    
    mutation_profile = pd.DataFrame(index=df.index,columns=sorted(gene2index.keys()))
    mutation_profile.update(df)
    mutation_profile = mutation_profile.fillna(0).as_matrix()
    return mutation_profile, pat2index

################################################################################################################################


def run_pca(propagated_profile):
    
    propagated_profile = propagated_profile.subtract(propagated_profile.mean())
    
    pca = PCA()
    pca.fit(propagated_profile)
    
    propagated_profile_pca = pca.transform(propagated_profile)
    propagated_profile_pca = pd.DataFrame(data=propagated_profile_pca,index=propagated_profile.index)
    PCs=['PC{}'.format(i+1) for i in propagated_profile_pca.columns]
    propagated_profile_pca.columns = PCs
    
    pca_components = pca.components_
    pca_components = pd.DataFrame(data=pca_components,columns=propagated_profile.columns)
    pca_components.index = PCs
    
    explained_variance_ratio = pca.explained_variance_ratio_
    
    return propagated_profile_pca, pca_components, explained_variance_ratio


def run_SpectralClustering(args):
    [propagated_profile_pca, n_clusters] = args[:2]
    cluster = SpectralClustering(affinity='nearest_neighbors', n_clusters=n_clusters, n_init=1000, gamma=0.5, 
                                 n_neighbors=170, assign_labels='discretize')
    cluster.fit(propagated_profile_pca)
#    print "Calinski-Harabasz Score with n_clusters=", n_clusters,"score:", metrics.calinski_harabaz_score(propagated_profile_pca, cluster.labels_) 
    return cluster.labels_


def run_KMeans(args):
    [propagated_profile_pca, n_clusters] = args[:2]
    cluster = KMeans(n_clusters=n_clusters, n_init=1000)
    cluster.fit(propagated_profile_pca)
#    print "Calinski-Harabasz Score with n_clusters=", n_clusters,"score:", metrics.calinski_harabaz_score(propagated_profile_pca, cluster.labels_) 
    return cluster.labels_


def run_AgglomerativeClustering(args):
    [propagated_profile_pca, n_clusters] = args[:2]
    cluster = AgglomerativeClustering(n_clusters=n_clusters, affinity='cosine', linkage='average')
    cluster.fit(propagated_profile_pca)
#    print "Calinski-Harabasz Score with n_clusters=", n_clusters,"score:", metrics.calinski_harabaz_score(propagated_profile_pca, cluster.labels_) 
    return cluster.labels_


def run_clustering_mp(propagated_profile_pca, maxK, func):
    
    n_processes = maxK-1
    pool = Pool(processes=n_processes)
    
    args = zip([propagated_profile_pca]*(maxK-1),range(2,maxK+1))
    labels = pool.map(func, args)
    
    labels = pd.DataFrame(data=np.array(labels).T,index=propagated_profile_pca.index,
                          columns=['K{}'.format(i) for i in range(2,maxK+1)])
    
    pool.close()
    pool.join()
    
    return labels


################################################################################################################################


def run_coxph(pat2surv_fn, labels, output_dir):
    
    os.system( 'mkdir -p ' + output_dir )
    
    pat2surv = pd.read_table(pat2surv_fn,index_col=0)
    pat2surv = pd.concat([pat2surv,labels],join='inner',axis=1)
    pat2surv.to_csv('{}/pat2surv2labels.txt'.format(output_dir),sep='\t')
    
    os.system("Rscript label2coxph.R {}".format(output_dir))
    return 0

################################################################################################################################


def run_fisher(gene,K):
    if sum(gene)==0:
        p_values = [1]*len(K.unique())
    else:
        p_values = []
        for k in np.sort(K.unique()):
            tab = pd.crosstab(gene > 0, K == k)
            oddsratio, p_value = fisher_exact(tab)
            if oddsratio < 1:
                p_value = 1
            p_values.append(p_value)
    return p_values


def run_ttest(gene,K):
    p_values = []
    for k in np.sort(K.unique()):
        cat1 = gene.loc[K == k]
        cat2 = gene.loc[K != k]
        t, p_value = ttest_ind(cat1, cat2)
        if t < 0:
            p_value = 1
        p_values.append(p_value)
    return p_values


def apply_fisher(args):
    [df, K] = args[:2]
    return df.apply(run_fisher, args=(K,))


def apply_ttest(args):
    [df, K] = args[:2]
    return df.apply(run_ttest, args=(K,))


def apply_test_mp(df, n_processes, func, K):
    
    n_processes = n_processes
    pool = Pool(processes=n_processes)

    df_split = np.array_split(df, n_processes, axis=1)
    args = zip(df_split,[K]*n_processes)
    df_test = pool.map(func, args)
    df_test = pd.concat(df_test)
    df_test = pd.DataFrame(list(df_test),index=df_test.index)
    df_test.sort_index(inplace =True)
    
    pool.close()
    pool.join()
    
    return df_test


def summarize_subtype_signatures(args):
    [output_dir, K, k, M_prop_ttest, M_fisher, pat2mut, ttest_fdr_cut] = args[:7]
    header="gene\tsubtype\tttest_p\tttest_fdr\tfisher_p\tfisher_fdr\tnPats\tperc_Muts\tcds\tnc\tboth\tup_down\n"
    ttest_p = M_prop_ttest.loc[:,k] / 2.
    ttest_fdr = ttest_p * len(ttest_p) / rankdata(ttest_p)
    ttest_fdr = ttest_fdr[ttest_fdr < ttest_fdr_cut]
    fisher_p = M_fisher.loc[ttest_fdr.index,k] / 2.
    fisher_fdr = fisher_p * len(fisher_p) / rankdata(fisher_p)
    signatures = ttest_fdr.index.tolist()
    subtype2signatures = {k:signatures}
    line_out=header
    pats_k=set(K[K==k].index)
    for gene in signatures:
        line_out+='{}\t{}\t{}\t{}\t{}\t{}\t{}\t'.format(gene,k,ttest_p[gene],ttest_fdr[gene],fisher_p[gene],fisher_fdr[gene],
                                                        len(pats_k))
        pats_cds=set(pat2mut.loc[(pat2mut.loc[:,1] == gene) & (pat2mut.loc[:,2] == 'cds'),0].unique()) & pats_k
        pats_nc=set(pat2mut.loc[(pat2mut.loc[:,1] == gene) & (pat2mut.loc[:,2] == 'nc'),0].unique()) & pats_k
        up_down=','.join(pat2mut.loc[(pat2mut.loc[:,1] == gene),3].unique())
        line_out+='{}\t{}\t{}\t{}\t{}'.format(len(pats_nc|pats_cds)/float(len(pats_k)),len(pats_cds-pats_nc),
                                              len(pats_nc-pats_cds),len(pats_nc&pats_cds),up_down)
        line_out+='\n'
    open('{}/subtype_{}.noa'.format(output_dir,k),'wb').write(line_out)
    return subtype2signatures


def summarize_subtype_signatures_mp(output_dir, K, M_prop_ttest, M_fisher, pat2mut, ttest_fdr_cut):
    
    n_subtypes = len(K.unique())
    n_processes = n_subtypes
    pool = Pool(processes=n_processes)

    k = range(n_subtypes)
    args = zip([output_dir]*n_subtypes, [K]*n_subtypes, k, [M_prop_ttest]*n_subtypes, [M_fisher]*n_subtypes, 
               [pat2mut]*n_subtypes, [ttest_fdr_cut]*n_subtypes)
    tmp = pool.map(summarize_subtype_signatures, args)

    subtype2signatures={}
    for i in tmp:
        subtype2signatures[i.keys()[0]]=i.values()[0]
    del tmp
    
    pool.close()
    pool.join()
    
    return subtype2signatures


def summarize_subnetworks(network_fn, subtype2signatures, output_dir):
    network = pd.read_table(network_fn)
    for subtype in subtype2signatures:
        subtype2network = network.loc[network.loc[:,'Gene1'].isin(subtype2signatures[subtype]) & 
                                      network.loc[:,'Gene2'].isin(subtype2signatures[subtype]),:]
        subtype2network.to_csv('{}/subtype_{}.txt'.format(output_dir, subtype),sep='\t',index=False)
    return 0


def subnetwork_wrapper(M, M_prop, K, test_n_processes, pat2mut_fn, network_fn, output_dir, ttest_fdr_cut):
    
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        
    M.sort_index(inplace =True)
    M_prop.sort_index(inplace =True)
    K.sort_index(inplace =True)

    M_prop_ttest = apply_test_mp(M_prop, test_n_processes, apply_ttest, K)
    print 'Finish ttest'
    M_fisher = apply_test_mp(M, test_n_processes, apply_fisher, K)
    print 'Finish Fisher exact test'

    pat2mut = pd.read_table(pat2mut_fn,header=None)
    print 'Summarizing subtype signatures'
    subtype2signatures = summarize_subtype_signatures_mp(output_dir, K, M_prop_ttest, M_fisher, pat2mut, ttest_fdr_cut)
    
    print 'Summarizing subnetworks'
    summarize_subnetworks(network_fn, subtype2signatures, output_dir)
    
    return 0

