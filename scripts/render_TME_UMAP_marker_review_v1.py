"""Visual review of a frozen candidate tree on audited existing representations.

UMAP is display only. Do not concatenate independently fitted source PCAs.
All cell-level axes/expression stay server-side; export static figures only.
"""
import argparse,collections,datetime,gc,hashlib,html,json,os,subprocess,sys,time,warnings
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import Normalize
from TME_visualization_sampling_v1 import balanced_indices

def load(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,v):
 with Path(p).open('x') as f:json.dump(v,f,ensure_ascii=False,indent=2,allow_nan=False);f.write('\n')
def checkpoint(root,stage,**v):
 p=root/'checkpoints';p.mkdir(exist_ok=True);save(p/(str(time.time_ns())+'.json'),dict(stage=stage,observed_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),**v))

MARKERS=['CD3D','CD3E','TRAC','CD4','IL7R','CD8A','CD8B','FOXP3','IL2RA','NKG7','KLRD1','GNLY','CD79A','MS4A1','MZB1','JCHAIN','LST1','LYZ','C1QA','C1QB','FCN1','S100A8','CLEC9A','XCR1','CD1C','FCER1A','CSF3R','FCGR3B','EPCAM','KRT8','KRT19','COL1A1','DCN','PECAM1','VWF','PROX1','CCL21','RGS5','MYH11','TPSAB1','CPA3','SOX10','PLP1','OLIG1','OLIG2','GFAP','AQP4','PMEL','MLANA','MITF','TYR','DCT','TYRP1','MKI67','GZMB','PDCD1','IFIT1']
MAJOR_COLORS={'T':'#377eb8','NK':'#1b9e77','T_NK':'#66c2a5','B':'#984ea3','Plasma':'#c994c7','Myeloid':'#ff7f00','Mast':'#e7298a','Epithelial':'#e6ab02','Fibroblast':'#4daf4a','Endothelial':'#e41a1c','Mural':'#a6761d','Glial':'#8c564b','Astrocyte':'#bcbd22','Neuron':'#17becf','Erythroid':'#9e0142','Platelet':'#f781bf','unresolved':'#bdbdbd'}

def fit_umap(x,cfg,seed):
 import umap
 assert x.ndim==2 and np.isfinite(x).all() and len(x)>=4
 model=umap.UMAP(n_neighbors=min(cfg['n_neighbors'],len(x)-1),min_dist=cfg['min_dist'],n_components=2,metric='euclidean',random_state=seed,n_jobs=1,n_epochs=cfg['n_epochs'],low_memory=True)
 return model.fit_transform(np.asarray(x,np.float32)).astype(np.float32)

def scatter(ax,xy,labels,title,colors=None,legend=True,point_size=2):
 labels=np.asarray(labels,str);levels=sorted(set(labels),key=lambda k:(k=='unresolved',k));cmap=plt.get_cmap('tab20');palette={k:(colors[k] if colors and k in colors else cmap(j%20)) for j,k in enumerate(levels)}
 order=['unresolved']+[k for k in levels if k!='unresolved'] if 'unresolved' in levels else levels
 for k in order:
  m=labels==k;ax.scatter(xy[m,0],xy[m,1],s=point_size,c=[palette[k]],label=k,linewidths=0,alpha=.72,rasterized=True)
 ax.set(title=title,xlabel='UMAP 1',ylabel='UMAP 2');ax.set_xticks([]);ax.set_yticks([]);ax.spines[['top','right']].set_visible(False)
 if legend:ax.legend(loc='center left',bbox_to_anchor=(1.01,.5),fontsize=7,markerscale=3,frameon=False)

def main(root):
 started=time.monotonic();cfg=load(root/'config/spec_v1.json');P=Path(cfg['producer']);D=Path(cfg['delivery']);A=Path(cfg['auditor']);J=Path(cfg['joined']);Q=Path(cfg['root_producer'])
 for r in [P,A,D,Q]:assert load(r/'provenance/supervisor_terminal_v1.json')['actual_child_exit_code']==0
 t=subprocess.run([sys.executable,'-B',str(root/'tests/test_TME_visualization_sampling_v1.py')],capture_output=True,text=True,timeout=60);(root/'provenance/sampling_test_stderr.txt').write_text(t.stderr);assert t.returncode==0 and 'Ran 3 tests' in t.stderr and 'skipped' not in t.stderr
 resource=load(load(Q/'config/spec_v1.json')['resource_policy_source']);df=subprocess.run(['df','-P',str(root)],capture_output=True,text=True,timeout=15);quota=subprocess.run(['quota','-g','-v','-w',resource['quota_group']],capture_output=True,text=True,timeout=15);assert df.returncode==quota.returncode==0;fs=df.stdout.splitlines()[-1].split()[0];qr=[r.split() for r in quota.stdout.splitlines() if r.split() and r.split()[0]==fs];assert len(qr)==1;v=os.statvfs(root);available=min((int(qr[0][3])-int(qr[0][1].rstrip('*')))*1024,v.f_bavail*v.f_frsize);assert available>resource['workspace_reserve_bytes']+cfg['output_cap_bytes'];save(root/'provenance/storage_preflight_v1.json',dict(minimum_physical_hard_quota_free_bytes=available,output_cap_bytes=cfg['output_cap_bytes'],existing_reserve_bytes=resource['workspace_reserve_bytes'],counts_not_copied=True))
 assert load(A/'provenance/terminal_v1.json')['producer_terminal_SHA256']==sha(P/'provenance/terminal_v1.json');pt=load(P/'provenance/terminal_v1.json');assert sha(P/'outputs/cluster_hierarchy_named_catalog_v1.json')==pt['catalog_SHA256'];pc=load(P/'config/spec_v1.json');nodes=load(P/'outputs/cluster_hierarchy_named_catalog_v1.json')['nodes'];leaves=[r for r in nodes if r['final_leaf']];lookup={(r['source'],r['source_node_code']):r for r in nodes};summary=load(P/'outputs/cluster_hierarchy_summary_v1.json');baseline=load(Path(pc['working_delivery'])/'outputs/accepted_TME_scoped_working_annotation_v1.json');w=load(baseline['producer_manifest']['path']);wb=Path(w['index']['path']).parent;wi={r['source']:r for r in load(w['index']['path'])}
 jm=load(J/'outputs/joined_native_atlas_delivery_v1.json');proj=load(jm['qualified_projection']);pi={(r['source'],r['model_id']):r for r in proj['records']};qi={r['source']:r for r in load(Q/'outputs/TME_qualified_root_space_index_v1.json')};sys.path.insert(0,str(J/'scripts'))
 from native_archive_collection_v1 import NativeArchiveCollection
 from native_panel_aggregate_v1 import gather_csr_rows
 native=NativeArchiveCollection(jm['archive_manifest']);O=root/'outputs';F=O/'figures';F.mkdir();figures=[];sample=[];source_cache={};coverage=[];by_type=collections.defaultdict(list);plt.rcParams.update({'font.size':9,'axes.titlesize':11,'figure.dpi':100,'savefig.dpi':140,'pdf.fonttype':42})
 def emit(fig,name,title,pdf=None):
  fig.savefig(F/(name+'.png'),dpi=cfg['dpi'],bbox_inches='tight',facecolor='white')
  if pdf:pdf.savefig(fig,dpi=cfg['dpi'],bbox_inches='tight')
  figures.append(dict(name=name+'.png',title=title));plt.close(fig)
 # Align the entire current fit scope before taking a source-balanced display.
 for entry in load(P/'outputs/cluster_hierarchy_index_v1.json'):
  source=entry['source'];mp=P/'outputs'/entry['metadata'];assert sha(mp)==entry['metadata_SHA256'];m=load(mp);assert sha(m['arrays']['path'])==m['arrays']['SHA256'] and sha(m['statistics']['path'])==m['statistics']['SHA256'];wr=wi[source];assert sha(wb/wr['arrays'])==wr['arrays_SHA256'];wm=load(wb/wr['metadata'])
  with np.load(m['arrays']['path'],allow_pickle=False) as f:pool=f['source_scope_position'];nativepos=f['source_native_position'];rows=f['native_row'];leaf=f['final_leaf_node']
  with np.load(wb/wr['arrays'],allow_pickle=False) as f:total=f['native_total_RNA'][pool];np.testing.assert_array_equal(rows,f['native_row'][pool]);oldmajor=np.asarray(wm['candidate_dictionaries']['consensus_major'],object)[f['consensus_major_code'][pool]];n_qualified=len(f['native_row'])
  record=pi[source,'global_overview__shared']['base_projection_record'];ap=Path(proj['base_projection'])/'outputs'/record['arrays'];assert sha(ap)==record['arrays_sha256']
  with np.load(ap,allow_pickle=False) as f:
   domain=f['source_native_position'];idx=np.searchsorted(domain,nativepos);hit=idx<len(domain);hit[hit]&=domain[idx[hit]]==nativepos[hit];raw=np.full((len(pool),50),np.nan,np.float32);corrected=raw.copy();supported=np.zeros(len(pool),bool);raw[hit]=f['raw_PCA'][idx[hit]];corrected[hit]=f['candidate_Harmony'][idx[hit]];supported[hit]=f['projection_numerical_supported'][idx[hit]];np.testing.assert_array_equal(f['native_row'][idx[hit]],rows[hit])
  eligible=hit&supported&np.isfinite(raw).all(axis=1)&np.isfinite(corrected).all(axis=1);seed=cfg['seed']+int(hashlib.sha256(source.encode()).hexdigest()[:6],16);pick=balanced_indices(leaf,eligible,cfg['global_per_source'],cfg['minimum_per_leaf'],seed);major=np.asarray([lookup[source,int(k)]['cluster_RNA_proposal']['major'] for k in leaf]);fine=np.asarray([lookup[source,int(k)]['cluster_RNA_proposal']['celltype'] for k in leaf]);status=np.asarray([lookup[source,int(k)]['cluster_RNA_proposal']['status'] for k in leaf]);labels=np.asarray([lookup[source,int(k)]['display_label'] for k in leaf]);genevals=np.full((len(pick),len(MARKERS)),np.nan,np.float32)
  with np.load(m['statistics']['path'],allow_pickle=False) as f:stats={k:f[k] for k in ['n','log_sum','detected']};cols=f['native_columns']
  by={g:j for j,g in enumerate(m['native_measured_genes'])};measured=[(j,int(cols[by[g]])) for j,g in enumerate(MARKERS) if g in by];csr=native.store(source)[native.legacy_config['counts_key'] if native.is_legacy(source) else native.primary.partitions[source]['counts_key']]
  for lo in range(0,len(pick),1024):
   inds=np.arange(lo,min(lo+1024,len(pick)));x=gather_csr_rows(csr,rows[pick[inds]],max_span=8192);rna=np.asarray(x.sum(axis=1)).ravel();np.testing.assert_array_equal(rna,total[pick[inds]]);v=x[:,[c for _,c in measured]].toarray();genevals[np.ix_(inds,[j for j,_ in measured])]=np.log1p(v*(10000/rna)[:,None])
  for typ in sorted(set(fine)):
   group_nodes=[r['source_node_code'] for r in leaves if r['source']==source and r['cluster_RNA_proposal']['celltype']==typ];n=int(stats['n'][group_nodes].sum());means=np.full(len(MARKERS),np.nan);pct=means.copy()
   for j,g in enumerate(MARKERS):
    if g in by:means[j]=stats['log_sum'][group_nodes,by[g]].sum()/n;pct[j]=stats['detected'][group_nodes,by[g]].sum()/n
   by_type[typ].append(dict(n=n,mean=means,pct=pct))
  sample.append(dict(source=source,cancer=cfg['source_cancer_broad'][source],raw=raw[pick],corrected=corrected[pick],major=major[pick],fine=fine[pick],status=status[pick],leaf=leaf[pick],labels=labels[pick],baseline_major=oldmajor[pick],genes=genevals,pool_index=pick));coverage.append(dict(source=source,n_current_qualified=n_qualified,n_fit_records=len(pool),n_unfitted_retained=n_qualified-len(pool),n_projection_eligible=int(eligible.sum()),n_display_sample=len(pick),n_leaf_groups=len(set(leaf)),n_leaf_groups_represented=len(set(leaf[pick])),n_leaves_without_eligible_global_projection=len(set(leaf)-set(leaf[eligible]))))
  if source in cfg['native_sources']:
   rootmeta=Q/'outputs'/qi[source]['metadata'];assert sha(rootmeta)==qi[source]['metadata_SHA256'];qm=load(rootmeta);assert sha(qm['model']['path'])==qm['model']['SHA256']
   with np.load(qm['model']['path'],allow_pickle=False) as f:np.testing.assert_array_equal(f['source_scope_position'],pool);nativepick=balanced_indices(leaf,np.ones(len(pool),bool),cfg['native_per_source'],cfg['minimum_per_leaf'],seed+1);embedding=f['fit_embedding'][nativepick]
   source_cache[source]=dict(embedding=embedding,major=major[nativepick],fine=fine[nativepick],leaf=leaf[nativepick],labels=labels[nativepick])
  checkpoint(root,'display_input_and_marker_source_complete',source=source,n_display_sample=len(pick));native._fields.clear();del raw,corrected,stats;gc.collect()
 arrays={k:np.concatenate([r[k] for r in sample]) for k in ['raw','corrected','major','fine','status','leaf','labels','baseline_major','genes']};arrays['source']=np.concatenate([np.repeat(r['source'],len(r['major'])) for r in sample]);arrays['cancer']=np.concatenate([np.repeat(r['cancer'],len(r['major'])) for r in sample]);assert sum(r['n_fit_records'] for r in coverage)==pt['totals']['n_pool_records'];assert sum(r['n_current_qualified'] for r in coverage)==pc['frozen_preflight']['n_current_qualified_records']
 save(O/'UMAP_display_scope_v1.json',dict(full_qualified_records=sum(r['n_current_qualified'] for r in coverage),full_fit_records=sum(r['n_fit_records'] for r in coverage),n_display_records=len(arrays['major']),n_leaf_groups=len(leaves),n_leaf_groups_represented=sum(r['n_leaf_groups_represented'] for r in coverage),source_coverage=coverage,display_is_source_balanced_and_not_cell_abundance=True,embedding_is_existing_transductive_reference_not_new_identity_fit=True,local_cluster_labels_are_overlay_not_joint_UMAP_clusters=True))
 # Full-profile marker dotplots are available before the slower UMAP layout.
 types=sorted(by_type,key=lambda t:(t=='unresolved',t));means=[];pcts=[]
 for typ in types:
  means.append(np.nanmean(np.stack([r['mean'] for r in by_type[typ]]),axis=0));pcts.append(np.nanmean(np.stack([r['pct'] for r in by_type[typ]]),axis=0))
 means=np.asarray(means);pcts=np.asarray(pcts);z=(means-np.nanmean(means,axis=0))/np.maximum(np.nanstd(means,axis=0),.05)
 with PdfPages(O/'TME_UMAP_marker_review.pdf') as pdf:
  fig,ax=plt.subplots(figsize=(23,10));xx,yy=np.meshgrid(np.arange(len(MARKERS)),np.arange(len(types)));ok=np.isfinite(z)&np.isfinite(pcts);sc=ax.scatter(xx[ok],yy[ok],s=110*pcts[ok],c=z[ok],cmap='RdBu_r',vmin=-2,vmax=2,linewidths=.2,edgecolors='#555555');ax.scatter(xx[~ok],yy[~ok],marker='x',s=10,c='#bbbbbb');ax.set_xticks(range(len(MARKERS)),MARKERS,rotation=90,fontsize=8);ax.set_yticks(range(len(types)),types);ax.invert_yaxis();ax.set_title('Full-fit native marker evidence | equal source weight within each proposed type');ax.set_xlabel('Color: gene-wise standardized mean log1p(CP10K); dot area: detection fraction; x: not measurable');fig.colorbar(sc,ax=ax,shrink=.45,label='Gene-wise z score');fig.tight_layout();emit(fig,'01_full_marker_dotplot','Full native marker evidence; classification markers are not independent validation',pdf)
  totals=summary['final_leaf_record_proposal_counts'];fig,axs=plt.subplots(1,2,figsize=(16,8));names=sorted(totals,key=totals.get);axs[0].barh(names,[totals[n] for n in names],color=['#bbbbbb' if n=='unresolved' else '#377eb8' for n in names]);axs[0].set(xlabel='All fitted records (not display sample)',title='Current RNA proposals: full denominator');sizes=np.array([r['n_records'] for r in leaves]);axs[1].hist(np.log10(sizes),bins=35,color='#4c78a8');axs[1].set(xlabel='log10(records per terminal cluster)',ylabel='Clusters',title=f'{len(leaves):,} terminal clusters; not shared biological types');fig.tight_layout();emit(fig,'02_annotation_completeness','Full-scope annotation coverage and terminal-cluster sizes',pdf)
  checkpoint(root,'global_Harmony_UMAP_started',n_display_records=len(arrays['major']));xy=fit_umap(arrays['corrected'],cfg,cfg['seed']);arrays['UMAP_Harmony']=xy
  fig,axs=plt.subplots(1,2,figsize=(16,6));scatter(axs[0],xy,arrays['major'],'Current cluster-major RNA proposals',MAJOR_COLORS);scatter(axs[1],xy,arrays['baseline_major'],'Baseline v5 majors on identical points',MAJOR_COLORS);fig.suptitle(f'Existing global source-Harmony reference; {len(xy):,} source-balanced display points',fontsize=12);fig.tight_layout();emit(fig,'03_UMAP_major_current_vs_v5','Same UMAP coordinates, current candidate versus baseline majors',pdf)
  fig,axs=plt.subplots(1,2,figsize=(17,7));scatter(axs[0],xy,arrays['fine'],'Current fine/coarse RNA proposals');scatter(axs[1],xy,arrays['cancer'],'Cancer context');fig.tight_layout();emit(fig,'04_UMAP_celltypes_and_cancer','Candidate fine identities and cancer context',pdf)
  fig,axs=plt.subplots(1,2,figsize=(15,6));scatter(axs[0],xy,arrays['source'],'Source context',legend=False);flags=np.where(arrays['major']=='unresolved','unresolved','provisional_named');scatter(axs[1],xy,flags,'Unresolved group locations',{'unresolved':'#d73027','provisional_named':'#d9d9d9'});fig.tight_layout();emit(fig,'05_UMAP_source_and_unresolved','Source context and unclassified groups',pdf)
  feature_genes=['CD3D','CD8A','FOXP3','LST1','EPCAM','COL1A1','PECAM1','SOX10','PLP1','MLANA','PMEL','OLIG2'];fig,axs=plt.subplots(3,4,figsize=(16,11))
  for ax,g in zip(axs.flat,feature_genes):
   v=arrays['genes'][:,MARKERS.index(g)];ok=np.isfinite(v);ax.scatter(xy[:,0],xy[:,1],s=1,c='#dddddd',linewidths=0,rasterized=True);q=np.argsort(np.nan_to_num(v,nan=-1));positive=v[ok&(v>0)];vmax=float(np.quantile(positive,.99)) if len(positive) else 1.;sc=ax.scatter(xy[q[ok[q]],0],xy[q[ok[q]],1],s=1.3,c=v[q[ok[q]]],cmap='Blues',vmin=0,vmax=max(vmax,.1),linewidths=0,rasterized=True);ax.set(title=f'{g} | measured {int(ok.sum()):,}');ax.set_xticks([]);ax.set_yticks([]);fig.colorbar(sc,ax=ax,shrink=.6)
  fig.suptitle('Native single-gene log1p(CP10K); missing assays remain grey',fontsize=12);fig.tight_layout();emit(fig,'06_UMAP_marker_features','Marker expression on the same frozen UMAP, including melanoma/glial boundary',pdf)
  checkpoint(root,'global_raw_PCA_UMAP_started');rawxy=fit_umap(arrays['raw'],cfg,cfg['seed']);fig,axs=plt.subplots(1,2,figsize=(16,6));scatter(axs[0],rawxy,arrays['major'],'Existing global raw-PCA UMAP',MAJOR_COLORS);scatter(axs[1],xy,arrays['major'],'Existing global candidate-Harmony UMAP',MAJOR_COLORS);fig.suptitle('Separate 2D layouts; visual comparison does not select an integration winner');fig.tight_layout();emit(fig,'07_UMAP_raw_vs_corrected','Existing raw and source-corrected representations; no new integration selection',pdf)
  sources=sorted(set(arrays['source']));fig,axs=plt.subplots(8,6,figsize=(18,22));limits=(np.quantile(xy[:,0],[.001,.999]),np.quantile(xy[:,1],[.001,.999]))
  for ax,source in zip(axs.flat,sources):
   m=arrays['source']==source;ax.scatter(xy[:,0],xy[:,1],s=.2,c='#dddddd',rasterized=True);ax.scatter(xy[m,0],xy[m,1],s=1,c='#377eb8',rasterized=True);ax.set_title(source.replace('CXG_','CXG:')[:31],fontsize=8);ax.set_xlim(limits[0]);ax.set_ylim(limits[1]);ax.set_xticks([]);ax.set_yticks([])
  for ax in list(axs.flat)[len(sources):]:ax.set_visible(False)
  fig.suptitle('Each source highlighted on the identical global UMAP');fig.tight_layout();emit(fig,'08_all_source_facets','All source distributions on common coordinates',pdf)
  for route,majors in [('T_NK',{'T','NK','T_NK'}),('Myeloid',{'Myeloid'}),('Epithelial',{'Epithelial'}),('Stromal',{'Fibroblast','Endothelial','Mural'})]:
   take=np.flatnonzero(np.isin(arrays['major'],list(majors)))
   if len(take)<30:continue
   if len(take)>cfg['lineage_cap']:take=np.sort(np.random.default_rng(cfg['seed']).choice(take,cfg['lineage_cap'],replace=False))
   checkpoint(root,'lineage_display_UMAP_started',route=route,n_records=len(take));u=fit_umap(arrays['corrected'][take],cfg,cfg['seed']);fig,axs=plt.subplots(1,2,figsize=(16,6));scatter(axs[0],u,arrays['fine'][take],route+' subset: candidate identities');scatter(axs[1],u,arrays['source'][take],route+' subset: source',legend=False);fig.suptitle('Subset UMAP of existing common representation; no new clustering or annotation');fig.tight_layout();emit(fig,'09_'+route+'_detail',route+' detail with source context',pdf)
  # Real source-local root PCA coordinates, each fitted in its own source.
  for source,c in source_cache.items():
   checkpoint(root,'source_native_UMAP_started',source=source,n_records=len(c['major']));u=fit_umap(c['embedding'],cfg,cfg['seed']);fig,axs=plt.subplots(1,2,figsize=(15,6));scatter(axs[0],u,c['major'],'Major RNA proposal',MAJOR_COLORS);leafcolors=np.asarray([plt.get_cmap('turbo')((int(k)*.61803398875)%1) for k in c['leaf']]);axs[1].scatter(u[:,0],u[:,1],s=3,c=leafcolors,linewidths=0,rasterized=True);axs[1].set(title='Terminal source-local clusters (colors reused)',xlabel='UMAP 1',ylabel='UMAP 2');axs[1].set_xticks([]);axs[1].set_yticks([]);fig.suptitle(source+' | '+cfg['source_cancer_detail'][source]+' | actual source-root PCA basis',fontsize=11);fig.tight_layout();emit(fig,'10_native_'+source,source+' native source-root space and final leaves',pdf)
 # Cell-level arrays are private server artifacts, excluded from export.
 np.savez_compressed(root/'private_display_arrays_v1.npz',**{k:v for k,v in arrays.items() if k not in ['labels']})
 assert sum(p.stat().st_size for p in root.rglob('*') if p.is_file())<=cfg['output_cap_bytes']
 counts=summary['final_leaf_record_proposal_counts'];unresolved=counts.get('unresolved',0);text=f"""<!doctype html><html lang='zh-CN'><meta charset='utf-8'><title>TME UMAP 与 marker 审阅</title><style>body{{font:16px system-ui;max-width:1500px;margin:32px auto;color:#17212b}}img{{width:100%;height:auto}}section{{margin:32px 0;border-top:1px solid #ddd;padding-top:16px}}table{{border-collapse:collapse}}td,th{{padding:8px 16px;border:1px solid #ddd}}</style><h1>TME 分群与注释可视化审阅</h1><p>完整来源合格范围 {sum(r['n_current_qualified'] for r in coverage):,} 条；参与本轮分群 {pt['totals']['n_pool_records']:,} 条；末端群 {len(leaves):,} 个。未决 {unresolved:,} 条（{unresolved/pt['totals']['n_pool_records']:.1%}）。本轮绘图 {len(xy):,} 个来源平衡代表点，覆盖的群见统计；抽样图不代表细胞丰度。</p><p>目前尚未完成所有大类定类。当前名称是 RNA 候选，Glial 在黑色素瘤来源的命名边界需结合 SOX10/PLP1 与 PMEL/MLANA/OLIG2 等基因重审；不由谱系名称判断良恶性。</p><p>UMAP 底图复用已审全局参考的原始 PCA/来源校正候选；颜色叠加本轮来源内群。不同来源独立 PCA 未直接拼接。来源局部图各自使用实际来源根群 PCA，坐标不可跨图比较。UMAP 不是聚类质量、批次消除或身份正确性的独立证明。</p><p>定类依据同细胞多基因模块和群级支持，不由展示 marker 单独定类。展示 marker 来自本群相对同来源其他根群或同父兄弟群的原生表达比较，数字是群编号。本轮排序是描述性表达证据，不是患者重复感知 DEG。dotplot 大小为检出比例、颜色为基因内标准化平均 log1p(CP10K)，缺测不补零。</p><p><a href='TME_UMAP_marker_review.pdf'>下载完整 PDF</a></p>"""
 text+=''.join('<section><h2>'+html.escape(f['title'])+'</h2><img loading="lazy" src="figures/'+f['name']+'"></section>' for f in figures);(O/'UMAP_marker_review.html').write_text(text+'</html>')
 exported=[O/'TME_UMAP_marker_review.pdf',O/'UMAP_marker_review.html',O/'UMAP_display_scope_v1.json',*sorted(F.glob('*.png'))];save(O/'visualization_manifest_v1.json',dict(files=[dict(path=str(p.relative_to(O)),bytes=p.stat().st_size,SHA256=sha(p)) for p in exported],n_figures=len(figures),cell_level_arrays_not_exported=True));save(root/'provenance/terminal_v1.json',dict(status='FROZEN_TME_UMAP_MARKER_VISUAL_REVIEW_RENDERED_NOT_NEW_ANNOTATION',completed_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),n_figures=len(figures),n_display_records=len(xy),n_leaf_groups_represented=sum(r['n_leaf_groups_represented'] for r in coverage),manifest_SHA256=sha(O/'visualization_manifest_v1.json'),no_counts_identity_or_CM_changes=True,wall_seconds=time.monotonic()-started))

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True);root=Path(p.parse_args().output)
 try:main(root)
 except BaseException:
  import traceback
  save(root/'provenance/failure_v1.json',dict(traceback=traceback.format_exc()));raise
