"""Fetch the exact public viewer inventory directly on the execution server."""
import argparse
import json
import urllib.request
from pathlib import Path
from types import SimpleNamespace

from prepare_CM02_panel_bridge_reference_scope_v1 import budget
from recover_SCP2510_public_processed_files_v1 import main as recover


def main(a):
    budget(a.root,512*1024**2)
    metadata=a.root/'server_public_source_metadata_v1';assert not metadata.exists();metadata.mkdir(mode=0o700)
    urls=dict(page='https://singlecell.broadinstitute.org/single_cell/study/SCP2510/human-lung-cancer-harbors-spatially-organized-stem-immunity-hubs-that-associate-with-response-to-immunotherapy',
              inventory='https://singlecell.broadinstitute.org/single_cell/api/v1/site/studies/SCP2510',
              explore='https://singlecell.broadinstitute.org/single_cell/api/v1/studies/SCP2510/explore')
    paths={}
    for name,url in urls.items():
        request=urllib.request.Request(url,headers={'Accept':'application/json' if name!='page' else 'text/html'})
        with urllib.request.urlopen(request,timeout=45) as r:data=r.read(8*1024**2+1)
        assert len(data)<=8*1024**2
        path=metadata/(name+('.html' if name=='page' else '.json'));path.write_bytes(data);paths[name]=path
    recover(SimpleNamespace(**paths,output=a.root/'SCP2510_author_processed_server_direct_v1'))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);main(p.parse_args())
