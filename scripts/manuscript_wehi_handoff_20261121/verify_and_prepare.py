#!/usr/bin/env python3
"""Verify and prepare the frozen manuscript bundle on its receiving server."""
import argparse,hashlib,json,os,shutil,stat,subprocess,sys,zipfile
from pathlib import Path,PurePosixPath
from datetime import datetime,timezone

def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()

def prepare(root):
    root=root.resolve();meta=json.loads((root/'transfer_manifest.json').read_text());archive=root/meta['archive']
    assert archive.stat().st_size==meta['archive_bytes'],'archive size mismatch'
    assert sha(archive)==meta['archive_sha256'],'archive SHA256 mismatch'
    target=root/'manuscript'
    with zipfile.ZipFile(archive) as z:
        assert z.testzip() is None,'archive CRC failure'
        assert len(z.infolist())==meta['archive_entries'],'entry count mismatch'
        for entry in z.infolist():
            p=PurePosixPath(entry.filename)
            assert not p.is_absolute() and '..' not in p.parts,'unsafe ZIP member'
            assert not stat.S_ISLNK(entry.external_attr>>16),'symlink not allowed'
        if not target.exists():
            temporary=root/'manuscript_extract_incoming'
            temporary.mkdir(exist_ok=False)
            for entry in z.infolist():
                q=temporary.joinpath(*PurePosixPath(entry.filename).parts)
                assert os.path.commonpath([str(temporary),str(q.resolve())])==str(temporary),'escaped path'
                if entry.is_dir():q.mkdir(parents=True,exist_ok=True);continue
                q.parent.mkdir(parents=True,exist_ok=True)
                with z.open(entry) as src,q.open('xb') as dst:shutil.copyfileobj(src,dst)
            temporary.rename(target)
    frozen=json.loads((target/'file_manifest.json').read_text())
    assert len(frozen['files'])==meta['manifest_files'],'file manifest count mismatch'
    for row in frozen['files']:
        p=(target/row['path']).resolve()
        assert os.path.commonpath([str(target),str(p)])==str(target),'unsafe manifest path'
        assert p.stat().st_size==row['bytes'] and sha(p)==row['sha256'],'file mismatch: '+row['path']
    process=subprocess.run([sys.executable,str(target/'plotting/code/verify_plot_inputs.py'),'--package-root',str(target/'plotting')],check=True,capture_output=True,text=True)
    (root/'claude_outputs').mkdir(exist_ok=True)
    receipt=root/'WEHI_RECEIPT.json'
    data={'status':'WEHI_RECEIVED_EXTRACTED_AND_VERIFIED','verified_UTC':datetime.now(timezone.utc).isoformat(),'archive_sha256':meta['archive_sha256'],'archive_bytes':archive.stat().st_size,'archive_entries':meta['archive_entries'],'verified_manifest_files':len(frozen['files']),'plot_input_verifier':process.stdout.strip(),'destination':str(root),'large_expression_matrices_included':False,'plotting_executed':False}
    if receipt.exists():
        prior=json.loads(receipt.read_text());assert prior['archive_sha256']==data['archive_sha256'] and prior['destination']==data['destination'],'existing receipt belongs to a different transfer'
        data=prior
    else:receipt.write_text(json.dumps(data,indent=2)+'\n')
    print(json.dumps(data,indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--handoff-root',type=Path,required=True);prepare(p.parse_args().handoff_root)
