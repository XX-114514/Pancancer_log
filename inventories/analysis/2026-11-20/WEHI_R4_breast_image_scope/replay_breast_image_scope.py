import pathlib,sys,json,hashlib,datetime
import h5py,numpy as np
native,results,out=map(pathlib.Path,sys.argv[1:]);r=json.loads((out/'complete.json').read_text());rows=[]
for q in r['rows']:
 s=q['section'];p=native/s/'native_PRIVATE.h5ad'
 receipt,=native.parent.glob('outputs*/spatial_controls_v1/'+s+'/complete.json');old=json.loads(receipt.read_text());assert q['input_sha256']==old['native_source_file_sha256']
 with h5py.File(p) as f:
  g=next(iter(f['uns/spatial'].values()));scale=float(g['scalefactors/tissue_hires_scalef'][()]);height,width,_=g['images/hires'].shape;xy=f['obsm/spatial'][:]*scale;raw=np.column_stack([f['obs/x_centroid'][:],f['obs/y_centroid'][:]]);nc=f['obs/nucleus_count'][:];na=f['obs/nucleus_area'][:];area=f['obs/cell_area'][:]
  outside=(xy[:,0]<0)|(xy[:,0]>=width)|(xy[:,1]<0)|(xy[:,1]>=height)
  assert int(outside.sum())==q['out_image'];assert int((nc==0).sum())==q['zero_nuclei'];assert np.array_equal(~np.isfinite(na),nc==0)
  assert int((nc>1).sum())==q['nuclei_gt1'] and int((na>area).sum())==q['nuclear_area_gt_cell']
  # Replay producer transform using explicit coordinate equations.
  c=np.asarray(q['affine_coefficients']);pred=np.column_stack([raw[:,0]*c[0,k]+raw[:,1]*c[1,k]+c[2,k] for k in [0,1]])*scale
  assert abs(np.sqrt(np.mean(np.sum((pred-xy)**2,axis=1)))-q['affine_rms_hires_pixels'])<1e-10
  for state,folder in [('broad_mural','imaging_breast_core_cells_v1'),('pc_low_MYH11','imaging_breast_mural_states_v1/pc_low_MYH11'),('MYH11_contractile','imaging_breast_mural_states_v1/MYH11_contractile')]:
   z=np.load(results/folder/(s+'_core_PRIVATE.npz'));assert z['xy'].shape==raw.shape and np.max(abs(z['xy']-raw))<1e-10
   for label,mask in [('EC',z['ec']),('mural',z['mural']),('ESM1_EC',z['ec']&(z['ESM1']>=1)),('RGS5_mural',z['mural']&(z['RGS5']>=1))]:
    rows.append(dict(section=s,state=state,population=label,coordinate_max_error_um=float(np.max(abs(z['xy']-raw))),cells=int(mask.sum()),out_image=int((mask&outside).sum()),zero_nuclei=int((mask&(nc==0)).sum()),nuclei_gt1=int((mask&(nc>1)).sum())))
report=dict(status='PASS_EXISTING_SOURCE_HASHES_NATIVE_COORDINATES_AND_ALL_GROUP_COUNTS',time_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),complete_sha256=hashlib.sha256((out/'complete.json').read_bytes()).hexdigest(),rows=rows,scope='All six source hashes match prior accepted receipts; all native coordinates agree with previously analyzed arrays within 1e-10 um; source rows and group counts replayed. No independent histology ground truth.',new_biological_tests=False)
(out/'independent_scope_replay.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(dict(status=report['status'],groups=len(rows),outside_core=sum(v['out_image'] for v in rows))))
