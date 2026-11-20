import pathlib,sys,json,hashlib,datetime
import h5py,numpy as np
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
root,out=map(pathlib.Path,sys.argv[1:]);out.mkdir(parents=True,exist_ok=True)
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
 return h.hexdigest()
rows=[];fig,axes=plt.subplots(3,4,figsize=(16,12))
for i,s in enumerate(['S01','S02','S03','S04','S05','S06']):
 p=root/s/'native_PRIVATE.h5ad'
 with h5py.File(p) as f:
  g=next(iter(f['uns/spatial'].values()));sc={k:float(v[()]) for k,v in g['scalefactors'].items()};im=g['images/hires'][:];obs=np.column_stack([f['obs/x_centroid'][:],f['obs/y_centroid'][:]]);xy=f['obsm/spatial'][:];u=f['obsm/xenium_um_spatial'][:]
  assert np.isfinite(obs).all() and np.isfinite(xy).all()
  design=np.column_stack([obs,np.ones(len(obs))]);coef=np.linalg.lstsq(design,xy,rcond=None)[0];res=np.linalg.norm((design@coef-xy)*sc['tissue_hires_scalef'],axis=1)
  points=xy*sc['tissue_hires_scalef'];inside=(points[:,0]>=0)&(points[:,0]<im.shape[1])&(points[:,1]>=0)&(points[:,1]<im.shape[0])
  area=f['obs/cell_area'][:];narea=f['obs/nucleus_area'][:];nc=f['obs/nucleus_count'][:]
  row=dict(section=s,cells=len(obs),image_shape=list(im.shape),scales=sc,input_sha256=sha(p),image_array_sha256=hashlib.sha256(im.tobytes()).hexdigest(),in_image=int(inside.sum()),out_image=int((~inside).sum()),affine_coefficients=coef.tolist(),affine_rms_hires_pixels=float(np.sqrt(np.mean(res**2))),affine_max_hires_pixels=float(res.max()),um_label_times_pixel_size_vs_native_max=float(abs(u*sc['microns_per_pixel']-obs).max()),cell_area_nonpositive=int((area<=0).sum()),cell_area_nonfinite=int((~np.isfinite(area)).sum()),nuclear_area_missing=int((~np.isfinite(narea)).sum()),zero_nuclei=int((nc==0).sum()),nuclei_gt1=int((nc>1).sum()),nuclear_area_gt_cell=int((narea>area).sum()),cell_area_quantiles=np.nanquantile(area,[0,.01,.5,.99,1]).tolist(),target_markers_used=False)
  rows.append(row)
  for j in [0,1]:
   ax=axes.flat[(i//2)*4+(i%2)*2+j];ax.imshow(im);ax.set_title(s+(' | image' if j==0 else ' | all-cell sample'));ax.axis('off')
   if j:
    # Deterministic evenly spaced row sample; no expression/identity selection.
    idx=np.linspace(0,len(points)-1,min(3000,len(points)),dtype=int);ax.scatter(points[idx,0],points[idx,1],s=.2,c='cyan',alpha=.45);ax.set_xlim(0,im.shape[1]);ax.set_ylim(im.shape[0],0)
 print(s,'complete',flush=True)
fig.suptitle('Embedded images and supplied image-coordinate overlay | not segmentation acceptance');fig.tight_layout();fig.savefig(out/'embedded_image_overview_PRIVATE.jpg',dpi=140);plt.close(fig)
report=dict(status='EMBEDDED_IMAGES_AND_COORDINATE_DIAGNOSTIC_ONLY',time_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),rows=rows,scope='All six pre-existing arrays; no cell boundaries or pathology labels; no target-driven ROI selection',coordinate_policy='Overlay stored obsm/spatial with stored hires scale. Affine fit is a consistency diagnostic, not independently validated ground truth.',biological_tests_run=False)
(out/'complete.json').write_text(json.dumps(report,indent=2)+'\n')
