"""Membership keys, append-only display numbers and version crosswalks.

This module does not infer cell identities or qualify scientific evidence. Callers
must bind accepted memberships and reviewed evidence before creating a release.
Hashes stay in dedicated columns; display labels contain one gene and a number.
"""
from collections import Counter, defaultdict
import hashlib
import re


def membership_key(source, member_ids):
    """An order-independent key of source-namespaced exact archive membership."""
    if not isinstance(source,str) or not source:
        raise ValueError('Stable source namespace required')
    ids=list(member_ids)
    if not ids or any(not isinstance(v,str) or not v for v in ids) or len(ids)!=len(set(ids)):
        raise ValueError('Nonempty unique native archive identities required')
    h=hashlib.sha256()
    for value in ['atlas_exact_membership_v1',source,*sorted(ids)]:
        raw=value.encode('utf-8');h.update(len(raw).to_bytes(8,'big'));h.update(raw)
    return h.hexdigest()


def assign_display_numbers(current_keys, previous=None):
    """Keep existing numbers; allocate new numbers after every historical one."""
    previous={} if previous is None else dict(previous)
    keys=list(current_keys)
    if len(keys)!=len(set(keys)) or any(not re.fullmatch(r'[0-9a-f]{64}',k) for k in keys):
        raise ValueError('Current membership keys must be unique SHA256 values')
    if any(not re.fullmatch(r'[0-9a-f]{64}',k) for k in previous):
        raise ValueError('Historical membership key invalid')
    if any(type(v) is not int or v<1 for v in previous.values()) or len(previous.values())!=len(set(previous.values())):
        raise ValueError('Historical display numbers must be unique positive integers')
    registry=dict(previous);next_number=max(registry.values(),default=0)+1
    for key in sorted(set(keys)-set(registry)):
        registry[key]=next_number;next_number+=1
    return {'current':{key:registry[key] for key in keys},'registry':registry,
            'inactive_historical_keys':sorted(set(previous)-set(keys))}


def select_display_marker(ordered_records, own_membership_key, allow_parent=False):
    """Evidence precedence without treating descriptive markers as DEG.

    Within each qualified evidence tier, keep the caller's bound ranking. A
    record without explicit accepted arithmetic/evidence cannot supply a marker.
    """
    priority={'own_subset_valid_DEG':0,'own_subset_expression_specificity':1,'parent_group_borrowed_explicit':2}
    candidates=[]
    for index,record in enumerate(ordered_records):
        kind=record.get('marker_evidence_type')
        if kind not in priority or record.get('review_status')!='accepted':
            continue
        gene=record.get('gene')
        if not isinstance(gene,str) or not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9.-]*',gene):
            raise ValueError('Display marker must be one explicit gene symbol')
        if kind.startswith('own_subset_') and record.get('membership_key')!=own_membership_key:
            raise ValueError('Own marker belongs to another exact subset')
        if kind=='parent_group_borrowed_explicit':
            if not allow_parent:continue
            if not record.get('parent_membership_key') or not record.get('borrow_reason'):
                raise ValueError('Borrowed parent marker requires explicit origin and reason')
        candidates.append((priority[kind],index,record))
    if not candidates:
        return {'marker':None,'marker_origin':'unavailable','marker_evidence_type':'no_accepted_own_marker','parent_marker_borrowed':False}
    _,_,record=min(candidates,key=lambda v:(v[0],v[1]))
    return {'marker':record['gene'],'marker_origin':record.get('membership_key') if record['marker_evidence_type'].startswith('own_subset_') else record['parent_membership_key'],
            'marker_evidence_type':record['marker_evidence_type'],'parent_marker_borrowed':record['marker_evidence_type']=='parent_group_borrowed_explicit',
            'evidence_record_sha256':record.get('evidence_record_sha256'),
            'borrow_reason':record.get('borrow_reason') if record['marker_evidence_type']=='parent_group_borrowed_explicit' else None}


def display_label(celltype, marker, number):
    if not isinstance(celltype,str) or not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_+.-]*',celltype):
        raise ValueError('Celltype display token must be explicit and pre-reviewed')
    if type(number) is not int or number<1:raise ValueError('Positive display number required')
    if marker is not None and (not isinstance(marker,str) or not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9.-]*',marker)):
        raise ValueError('Display marker must contain exactly one gene token')
    return f'{celltype}_{marker or "marker_unresolved"}_{number}'


def version_crosswalk(source, native_ids, old_labels, new_labels, new_status):
    """Count the complete unchanged archive axis; split and merge are separate.

    Both flags may hold for a many-to-many mapping. Unresolved is an explicit
    reviewed status, never guessed from an underscore or a display number.
    """
    ids,old,new,status=map(list,[native_ids,old_labels,new_labels,new_status])
    if not source or not len(ids) or any(len(a)!=len(ids) for a in [old,new,status]):
        raise ValueError('Crosswalk full source archive axes differ')
    if len(ids)!=len(set(ids)) or any(not isinstance(i,str) or not i for i in ids):
        raise ValueError('Crosswalk native IDs must be unique within source')
    if any(not isinstance(v,str) or not v for a in [old,new,status] for v in a):
        raise ValueError('Missing label/status must be an explicit nonempty category')
    destinations=defaultdict(set);origins=defaultdict(set)
    for before,after in zip(old,new):destinations[before].add(after);origins[after].add(before)
    counts=Counter(zip(old,new,status));records=[]
    for (before,after,state),count in sorted(counts.items()):
        split=len(destinations[before])>1;merge=len(origins[after])>1
        relation='many_to_many' if split and merge else 'split' if split else 'merge' if merge else 'one_to_one'
        records.append({'source':source,'old_label':before,'new_label':after,'new_identity_status':state,'n_archive_records':count,
                        'old_label_splits':split,'new_label_merges':merge,'relation':relation,
                        'mapping_is_identity_truth':False})
    if sum(r['n_archive_records'] for r in records)!=len(ids):raise AssertionError('Crosswalk dropped archive rows')
    return records
