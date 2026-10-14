"""Fetch public reference gene sets only; no patient genes sent to services."""
import concurrent.futures
import hashlib,json,os,re,urllib.request
from pathlib import Path
from datetime import datetime,timezone


def get(url):
    return urllib.request.urlopen(url,timeout=25).read()


def main():
    out=Path(os.environ['PANCANCER_PUBLIC_HALLMARK_OUTPUT'])
    assert not out.exists();out.mkdir(parents=True)
    page='https://www.gsea-msigdb.org/gsea/msigdb/human/genesets.jsp?collection=H'
    html=get(page).decode()
    names=sorted(set(re.findall(r'>(HALLMARK_[A-Z0-9_]+)</a>',html)))
    assert len(names)==50, 'Official collection must contain exactly 50 sets'
    (out/'public_collection_page.html').write_text(html)
    def one(name):
        url='https://www.gsea-msigdb.org/gsea/msigdb/human/download_geneset.jsp?geneSetName='+name+'&fileType=json'
        b=get(url);obj=json.loads(b);assert set(obj)=={name}
        assert len(obj[name]['geneSymbols'])>0 and obj[name]['collection']=='H'
        (out/(name+'.json')).write_bytes(b)
        return dict(name=name,url=url,sha256=hashlib.sha256(b).hexdigest(),genes=len(obj[name]['geneSymbols']))
    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
        records=list(pool.map(one,names))
    (out/'reference_provenance.json').write_text(json.dumps(dict(status='PUBLIC_50_HALLMARK_REFERENCES_FETCHED',
      observed_UTC=datetime.now(timezone.utc).isoformat(),collection_page=page,records=records,
      patient_or_learned_genes_sent=False),indent=2)+'\n')
    print('PUBLIC_50_HALLMARK_REFERENCES_FETCHED_NO_PATIENT_INPUTS')


if __name__=='__main__':main()
