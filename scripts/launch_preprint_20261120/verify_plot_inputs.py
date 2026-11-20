#!/usr/bin/env python3
"""Check the delivered input manifest before rendering a private plot package."""
import argparse,hashlib,json
from pathlib import Path

def verify(root):
    root=root.resolve()
    manifest=json.loads((root/'input_manifest.json').read_text())
    for row in manifest['files']:
        path=(root/row['path']).resolve()
        if not path.is_relative_to(root):
            raise ValueError('Manifest path escapes the package')
        if path.stat().st_size != row['bytes']:
            raise ValueError('Size mismatch: '+row['path'])
        h=hashlib.sha256()
        with path.open('rb') as handle:
            for block in iter(lambda:handle.read(1048576),b''):h.update(block)
        if h.hexdigest()!=row['sha256']:
            raise ValueError('SHA256 mismatch: '+row['path'])
    print('INPUT_MANIFEST_VERIFIED',len(manifest['files']))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--package-root',type=Path,required=True)
    verify(p.parse_args().package_root)
