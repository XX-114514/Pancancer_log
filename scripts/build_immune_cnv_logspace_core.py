#!/usr/bin/env python3
"""Version a minimal original-core correction; preserve all inference rules."""
import argparse,ast,copy,hashlib,json
from datetime import datetime,timezone
from pathlib import Path


class CoreInput(ast.NodeTransformer):
    def __init__(self):self.changed=0
    def visit_Call(self,node):
        self.generic_visit(node)
        if ast.unparse(node.func)=='cnv.tl.infercnv':
            layers=[k for k in node.keywords if k.arg=='layer']
            if len(layers)!=1 or not isinstance(layers[0].value,ast.Constant) or layers[0].value.value!='counts':
                raise ValueError('Unexpected original infercnv input layer')
            layers[0].value=ast.parse('str(_cfg_value(cfg,"infercnv_expression_layer","counts"))',mode='eval').body
            if any(k.arg=='n_jobs' for k in node.keywords):raise ValueError('Original jobs already configured')
            node.keywords.append(ast.keyword(arg='n_jobs',value=ast.parse('int(_cfg_value(cfg,"infercnv_n_jobs",2))',mode='eval').body))
            self.changed+=1
        return node


def build(a):
    original=a.original.read_bytes();tree=ast.parse(original.decode())
    patch=CoreInput();tree=patch.visit(tree);ast.fix_missing_locations(tree)
    if patch.changed!=1:raise ValueError('Expected exactly one original infercnv entrypoint')
    restored=copy.deepcopy(tree)
    for node in ast.walk(restored):
        if isinstance(node,ast.Call) and ast.unparse(node.func)=='cnv.tl.infercnv':
            for keyword in node.keywords:
                if keyword.arg=='layer':keyword.value=ast.Constant(value='counts')
            node.keywords=[k for k in node.keywords if k.arg!='n_jobs']
    if ast.dump(restored,include_attributes=False)!=ast.dump(ast.parse(original.decode()),include_attributes=False):
        raise ValueError('Unexpected scientific-rule AST change')
    text=ast.unparse(tree)+'\n';compile(text,str(a.output),'exec')
    with a.output.open('x') as f:f.write(text)
    report=dict(observed_utc=datetime.now(timezone.utc).isoformat(),original_core_sha256=hashlib.sha256(original).hexdigest(),
        corrected_core_sha256=hashlib.sha256(text.encode()).hexdigest(),infercnv_call_sites_changed=patch.changed,
        changes=['explicit_sample_log_expression_layer_parameter','bounded_infercnv_process_workers_default2'],
        unchanged=['candidate_reference_masks','genomic_window','reference_quantile','CNV_clustering_and_small_cluster_merge','cluster_malignant_rules','CopyKAT_disabled','all_diagnostic_plots'],
        original_file_untouched=True,restored_AST_exactly_original=True)
    with a.report.open('x') as f:json.dump(report,f,indent=2)
    print(json.dumps(report))


if __name__=='__main__':
    p=argparse.ArgumentParser()
    for n in ['original','output','report']:p.add_argument('--'+n,type=Path,required=True)
    build(p.parse_args())
