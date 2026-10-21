from pathlib import Path
import sys,unittest,io,tarfile,hashlib,urllib.request
from unittest.mock import patch
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import qualify_CM02_three_public_sources_RAM_v1 as m

class Receiver(unittest.TestCase):
    def archive(self,names):
        stream=io.BytesIO()
        with tarfile.open(fileobj=stream,mode='w:') as t:
            for name,body in names:
                item=tarfile.TarInfo(name);item.size=len(body);t.addfile(item,io.BytesIO(body))
        stream.seek(0);return stream
    def test_exact_set_hash_and_duplicate_rejection(self):
        bodies={n:b'\x1f\x8bsynthetic'+n.encode() for n in m.SOURCES}
        expected={n:hashlib.sha256(b).hexdigest() for n,b in bodies.items()}
        with patch.object(m,'SOURCES',expected):
            cache,proof=m.receive_exact_sources(self.archive(bodies.items()));self.assertEqual(len(proof),3)
            with self.assertRaises(AssertionError):m.receive_exact_sources(self.archive(list(bodies.items())[:2]))
            with self.assertRaises(AssertionError):m.receive_exact_sources(self.archive([list(bodies.items())[0]]*3))
            broken=dict(bodies);broken[next(iter(broken))]+=b'x'
            with self.assertRaises(AssertionError):m.receive_exact_sources(self.archive(broken.items()))
    def test_request_cache_restored_on_error(self):
        old=urllib.request.urlopen
        def callback():
            self.assertEqual(urllib.request.urlopen('public').read(),b'synthetic')
            urllib.request.urlopen('undeclared')
        with self.assertRaises(AssertionError):m.with_cached_sources({'public':b'synthetic'},callback)
        self.assertIs(urllib.request.urlopen,old)
if __name__=='__main__':unittest.main()
