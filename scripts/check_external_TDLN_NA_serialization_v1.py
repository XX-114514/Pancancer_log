"""Exercise the actual R serializer call with positive/negative/missing fixtures."""
from pathlib import Path
from copy import deepcopy
import argparse
import json
import subprocess
from check_external_TDLN_fixed_acceptance_v1 import artificial_summary
from accept_external_TDLN_fixed_spatial_v1 import validate_summary, reject_individual_payload, STATE_MARKERS, ENDS

R_CHECK = r'''
args<-commandArgs(trailingOnly=TRUE)
input<-jsonlite::fromJSON(file("stdin"),simplifyVector=FALSE)
result<-input
result$summary$phenotype$state_records[[6]]$cell_rows<-NA_integer_
for(j in c(9L,10L,11L,12L,15L)) {
  for(k in c("median_effect","bootstrap_lower","bootstrap_upper","primary_p","BH","BY"))
    result$summary$spatial$primary_rows[[j]][[k]]<-NA_real_
}
calls<-list()
walk<-function(e) {
 if(is.call(e)&&identical(e[[1]],quote(jsonlite::toJSON)))calls[[length(calls)+1L]]<<-e
 if(is.call(e)||is.expression(e))for(z in as.list(e))walk(z)
}
for(e in parse(args[1]))walk(e)
stopifnot(length(calls)==1L)
cat(eval(calls[[1]]),"\n")
'''

def check(Rscript, directory):
    directory=Path(directory)
    config=dict(physical_panel=sorted({m for k,v in STATE_MARKERS.items() if k!=ENDS[2] for m in v}),
                spatial=dict(minimum_cases=5,primary_grid=4,permutations=999))
    fixture={'kind':'SOURCE_ONLY_EXTERNAL_FIXED_WORKFLOW_COMPLETE','summary':artificial_summary()}
    results={}
    for version in ('v1','v2'):
        command=[Rscript,'--vanilla','-e',R_CHECK,str(directory/f'run_external_TDLN_fixed_spatial_native_RAM_{version}.R')]
        run=subprocess.run(command,input=json.dumps(fixture).encode(),stdout=subprocess.PIPE,stderr=subprocess.PIPE,timeout=30)
        assert run.returncode==0, 'ARTIFICIAL_R_SERIALIZER_FAILED'
        value=json.loads(run.stdout)
        if version=='v1':
            assert value['summary']['phenotype']['state_records'][5]['cell_rows']=='NA'
            try:reject_individual_payload(value)
            except AssertionError:pass
            else:raise AssertionError('LEGACY_MISSING_VALUE_FAILURE_NOT_REPRODUCED')
            results['legacy_failure_reproduced']=True
        else:
            reject_individual_payload(value)
            accepted=validate_summary(value['summary'],config)
            assert accepted['Python_aggregate_audit']['estimable_primary_slots']==10
            assert value['summary']['phenotype']['state_records'][5]['cell_rows'] is None
            for j in (8,9,10,11,14):
                assert value['summary']['spatial']['primary_rows'][j]['primary_p'] is None
            bad=deepcopy(value);bad['summary']['patient_ids']=['ARTIFICIAL_CANARY']
            try:reject_individual_payload(bad)
            except AssertionError:pass
            else:raise AssertionError('INDIVIDUAL_GUARD_WEAKENED')
            results.update(actual_v2_R_call_passed_full_15_slot_consumer=True,
                           positive_negative_and_structural_missing_preserved=True,
                           individual_canary_still_rejected=True)
    print(json.dumps(dict(status='PASS_ACTUAL_R_NA_SERIALIZATION_AND_FULL_AGGREGATE_CONSUMER',
                         checks=results,no_real_native_or_individual_data_read=True)))

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--Rscript',required=True)
    parser.add_argument('--code-directory',default=str(Path(__file__).parent));args=parser.parse_args()
    check(args.Rscript,args.code_directory)
