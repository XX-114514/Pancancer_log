"""Synthetic integrity, budget and biological-scope safeguards."""
import io,sys,unittest,zipfile
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from Danenberg_bounded_member_transport_v1 import RangeFile,IntegrityError,BudgetError,RangeStatusError
from qualify_Danenberg_native_cells_v1 import qc_mask,MARKERS,author_alias
import pandas as pd

URL='https://zenodo.org/api/records/7324285/files/test.zip/content'

class Response:
    status=206
    def __init__(self,body,headers):self.body=body;self.headers=headers
    def getheader(self,k):return self.headers.get(k)
    def read(self,n):return self.body[:n]
    def close(self):pass

class Connection:
    def __init__(self,data,change=None):self.data=data;self.change=change;self.calls=0
    def request(self,method,path,headers):self.range=headers['Range']
    def getresponse(self):
        self.calls+=1;a,b=map(int,self.range[6:].split('-'))
        h={'Content-Range':f'bytes {a}-{b}/{len(self.data)}','Content-Length':str(b-a+1),'ETag':'same'}
        r=Response(self.data[a:b+1],h)
        if self.change:self.change(r,self.calls)
        return r
    def close(self):pass

def reader(data=b'0123456789abcdef',change=None):
    c=Connection(data,change);r=RangeFile(URL,len(data),7897,connection_factory=lambda:c)
    r.BLOCK=4;return r,c

class Tests(unittest.TestCase):
    def test_cache_and_seek(self):
        r,c=reader();self.assertEqual(r.read(6),b'012345');self.assertEqual(c.calls,2)
        r.seek(1);self.assertEqual(r.read(4),b'1234');self.assertEqual(c.calls,2)
        r.seek(-2,2);self.assertEqual(r.read(),b'ef')
    def test_status_200_rejected(self):
        r,c=reader(change=lambda x,i:setattr(x,'status',200))
        with self.assertRaises(RangeStatusError):r.read(1)
        self.assertEqual(c.calls,1)
    def test_identity_change_rejected(self):
        r,c=reader(change=lambda x,i:x.headers.update(ETag=str(i)))
        with self.assertRaises(IntegrityError):r.read(8)
    def test_lengths_and_range_rejected(self):
        for change in [lambda x,i:setattr(x,'body',b''),lambda x,i:setattr(x,'body',b'12345'),
                       lambda x,i:x.headers.update({'Content-Range':'invalid'})]:
            r,c=reader(change=change)
            with self.assertRaises(IntegrityError):r.read(1)
            self.assertEqual(c.calls,1)
    def test_budget_and_source_scope(self):
        r,c=reader();r.MAX_REQUESTED=3
        with self.assertRaises(BudgetError):r.read(1)
        self.assertEqual(c.calls,0)
        with self.assertRaises(IntegrityError):RangeFile(URL.replace('zenodo.org','example.org'),8,7897)
    def test_CRC_rejection(self):
        b=io.BytesIO()
        with zipfile.ZipFile(b,'w',zipfile.ZIP_STORED) as z:z.writestr('x',b'hello')
        raw=bytearray(b.getvalue());raw[31]^=1
        with zipfile.ZipFile(io.BytesIO(raw)) as z:
            with self.assertRaises(zipfile.BadZipFile):z.read('x')
    def test_QC_sample_not_malignancy_and_negative_preserved(self):
        rows=[]
        for tumour,artifact,epithelial in [('1','0','0'),('1','1','0'),('0','0','1'),('1','unknown','0')]:
            rows.append(dict(is_tumour=tumour,is_hotAggregate=artifact,is_epithelial=epithelial,
                is_normal='0',is_dcis='0',**{m:'-1' for m in MARKERS},
                Location_Center_X='0',Location_Center_Y='1',AreaShape_Area='2'))
        good,v,reasons=qc_mask(pd.DataFrame(rows));self.assertEqual(good.tolist(),[True,False,False,False])
        self.assertEqual(v.FOXP3.iloc[0],-1);self.assertEqual(sum(reasons.values()),3)
    def test_author_labels_stay_mixed(self):
        self.assertEqual(author_alias('Treg & Tex'),'T_{Reg} & T_{Ex}')
        self.assertEqual(author_alias('CD4+ T cells & APCs'),'CD4^{+} T cells & APCs')

if __name__=='__main__':unittest.main()
