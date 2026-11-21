"""Render only the three requested figures from pinned branch source blocks.

This adapter selects complete AST statements, avoiding unrelated analyses.
Inputs are the external frozen portable plotting package; no fitting occurs.
"""
import argparse
import ast
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
PINNED = {
    'render_R4_portable.py': '52a1b610f7aa9025dd016fe706b70f97c2cf97a5e5fea741085597f38e0a84b1',
    'render_story_panels.py': 'a6c10058a79381dbf5129e0a234e6afa69488cd63655a5da93036c7355d7a515',
}

def subset(filename, function, ranges):
    path = HERE / filename
    if hashlib.sha256(path.read_bytes()).hexdigest() != PINNED[filename]:
        raise ValueError('Renderer changed: review target AST ranges before running')
    tree = ast.parse(path.read_text())
    for node in tree.body:
        if isinstance(node, ast.FunctionDef) and node.name == function:
            node.body = [s for s in node.body if any(a <= s.lineno <= b for a,b in ranges)]
            assert node.body
    tree.body = [n for n in tree.body if not isinstance(n, ast.If)]
    ns = {'__file__': str(path), '__name__': 'target_renderer'}
    exec(compile(tree, str(path), 'exec'), ns)
    return ns

def main(inputs, output):
    output.mkdir(parents=True, exist_ok=False)
    r = subset('render_R4_portable.py', 'render', [(109,128),(136,136),(170,206),(370,374)])
    r['render'](inputs/'r4', inputs/'r4/R4_preprint_20261120_v1', output/'r4')
    s = subset('render_story_panels.py', 'run', [(25,25),(27,27),(57,66)])
    s['run'](inputs, output/'story')
    import numpy as np
    import pandas as pd
    with np.load(inputs/'r4/R4_preprint_20261120_v1/atlas/atlas_embedding_PRIVATE.npz') as z:
        nodes=z['nodes'].tolist()
        assert len(nodes)==44 and len(z['xy'])==1974309
    pd.DataFrame({'number':range(1,45),'identity':nodes}).to_csv(output/'identity_key_1-44.tsv', sep='\t',index=False)
    with np.load(inputs/'immune/immune_embedding.npz') as z:
        majors=z['majors'].tolist()
    excluded={'Unknown','Ambiguous','Unresolved','T_NK_boundary'}
    pd.DataFrame({'identity':majors,'direct_label':[x not in excluded for x in majors]}).to_csv(output/'immune_label_key.tsv',sep='\t',index=False)
    names={'r4/Figure2_unified_atlas':'Figure1_resolved_atlas','r4/FigureS10_full_identity_UMAP':'FigureS06_identity_key','story/FigureS07_immune_atlas':'FigureS07_immune_atlas'}
    for old,new in names.items():
        for ext in ['png','pdf','svg']:
            (output/(new+'.'+ext)).write_bytes((output/(old+'.'+ext)).read_bytes())
    sources=[inputs/'r4/R4_preprint_20261120_v1/atlas/atlas_embedding_PRIVATE.npz',inputs/'immune/immune_embedding.npz']
    receipt={'coordinates_refitted':False,'identity_count':44,'major_count':8,'input_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in sources},'source_branch':'codex/0906-umap-labels','files':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in output.glob('*') if p.is_file()}}
    (output/'redraw_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--inputs',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();main(a.inputs,a.output)
