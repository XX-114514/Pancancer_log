"""Fixed author-name semantics; mixed T/NK is not a pure negative reference."""
import re

KINDS = ['EXPLICIT_GAMMA_DELTA_AUTHOR', 'DELTA_SUBTYPE_TEXT_ONLY',
         'GAMMA_MARKER_TEXT_ONLY', 'EXPLICIT_ALPHA_BETA_AUTHOR',
         'EXPLICIT_NK_AUTHOR', 'OTHER_T_AUTHOR', 'UNRESOLVED_OR_OTHER',
         'CONFLICTING_GAMMA_ALPHA_BETA_LABEL', 'MIXED_T_NK_AUTHOR']


def label_kind(value):
    s = str(value).casefold()
    gamma = bool(re.search(r'γ\s*δ|gamma[\s_\-/]*delta|(?:^|[^a-z0-9])(?:gd[\s_-]*t|t[\s_-]*gd|gdt)(?:$|[^a-z0-9])', s))
    ab = bool(re.search(r'α\s*β|alpha[\s_\-/]*beta|(?:^|[^a-z0-9])(?:abt|tab)(?:$|[^a-z0-9])', s))
    if gamma and ab:
        return 7
    if gamma:
        return 0
    if re.search(r'v\s*δ\s*[12]|v[\s_-]*(?:delta|d)[\s_-]*[12](?:$|[^0-9])', s):
        return 1
    if re.search(r'(?:^|[^a-z0-9])tr(?:dc|gc[12]|dv[123]|gv9)(?:$|[^a-z0-9])', s):
        return 2
    if ab:
        return 3
    nk = bool(re.search(r'(?:^|[^a-z0-9])nk(?:$|[^a-z0-9])|natural[\s_-]*killer', s))
    t = bool(re.search(r'(?:^|[^a-z0-9])t(?:$|[^a-z0-9])|(?:^|[^a-z0-9])(?:cd[48]|treg)(?:$|[^a-z0-9])', s))
    if nk and t:
        return 8
    if nk:
        return 4
    if t:
        return 5
    return 6


def fixture():
    examples = {'T_gdT_TRDV1':0, 'γδ T':0, 'gamma-delta T cells':0,
                'Vdelta1':1, 'TRDC+ T':2, 'alpha-beta T':3, 'NK':4,
                'CD8 T':5, 'T cells':5, 'unknown':6,
                'gdT/alpha-beta mixture':7, 'T_NK_boundary':8, 'T/NK cells':8}
    assert all(label_kind(k)==v for k,v in examples.items())
    return dict(status='PASS_MIXED_T_NK_NOT_PURE_REFERENCE_SEMANTICS', examples=len(examples))
