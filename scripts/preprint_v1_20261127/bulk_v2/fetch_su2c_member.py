"""Fetch one member of the public SU2C-MARK Zenodo archive via HTTP range requests; verify CRC32 and SHA256 against the original receipt."""
import io, zipfile, urllib.request, hashlib, zlib, sys
URL = 'https://zenodo.org/records/11179623/files/Source%20Data%203.zip?download=1'
MEMBER = 'Source Data/RNA/SU2C-MARK_Harmonized_rnaseqc_tpm_v1.gct'
SHA = '892f8687ed42eb017894aa9c99d4d4a04f5af52037f88709ee7913df36166160'
OUT = sys.argv[1]
class RangeFile(io.RawIOBase):
    def __init__(self, url):
        self.url = url; self.pos = 0
        r = urllib.request.urlopen(urllib.request.Request(url, method='HEAD'), timeout=60); self.size = int(r.headers['Content-Length'])
        self.cache = {}
    def seekable(self): return True
    def readable(self): return True
    def tell(self): return self.pos
    def seek(self, off, whence=0):
        self.pos = off if whence == 0 else self.pos + off if whence == 1 else self.size + off; return self.pos
    def readinto(self, b):
        n = min(len(b), self.size - self.pos)
        if n <= 0: return 0
        req = urllib.request.Request(self.url, headers={'Range': f'bytes={self.pos}-{self.pos + n - 1}', 'User-Agent': 'generic-research-client'})
        data = urllib.request.urlopen(req, timeout=300).read()
        b[:len(data)] = data; self.pos += len(data); return len(data)
f = io.BufferedReader(RangeFile(URL), buffer_size=8 * 1024 * 1024)
z = zipfile.ZipFile(f)
info = z.getinfo(MEMBER); print('member', info.file_size, 'compressed', info.compress_size, 'crc', info.CRC, flush=True)
data = z.read(MEMBER)
assert zlib.crc32(data) & 0xffffffff == info.CRC
h = hashlib.sha256(data).hexdigest(); print('sha256', h, 'match_receipt', h == SHA)
assert h == SHA
open(OUT, 'xb').write(data); print('written', len(data))
