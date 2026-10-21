"""Bounded read-only file adapter for metadata qualification; cache stays RAM."""
import io


class RangeReader(io.RawIOBase):
    def __init__(self,size,fetch,block=65536,budget=16*1024**2):
        super().__init__();self.size=int(size);self.fetch=fetch;self.block=block;self.budget=budget
        self.position=0;self.cache={};self.transferred=0

    def readable(self):return True
    def seekable(self):return True
    def writable(self):return False
    def tell(self):return self.position

    def seek(self,offset,whence=0):
        point=offset if whence==0 else self.position+offset if whence==1 else self.size+offset if whence==2 else -1
        if point<0:raise ValueError('INVALID_SEEK')
        self.position=int(point);return self.position

    def read(self,size=-1):
        if size<0:size=self.size-self.position
        if size>1024**2:raise ValueError('SINGLE_METADATA_READ_OVER_1MIB')
        end=min(self.size,self.position+size);values=[]
        while self.position<end:
            block_start=(self.position//self.block)*self.block
            if block_start not in self.cache:
                length=min(self.block,self.size-block_start)
                if self.transferred+length>self.budget:raise MemoryError('METADATA_RANGE_BUDGET_EXCEEDED')
                value=self.fetch(block_start,length)
                if len(value)!=length:raise ValueError('INCOMPLETE_METADATA_RANGE')
                self.cache[block_start]=value;self.transferred+=length
            value=self.cache[block_start];lo=self.position-block_start;take=min(len(value)-lo,end-self.position)
            values.append(value[lo:lo+take]);self.position+=take
        return b''.join(values)

    def readinto(self,b):
        value=self.read(len(b));b[:len(value)]=value;return len(value)


def synthetic():
    data=bytes(range(251))*500;calls=[]
    def fetch(start,length):calls.append((start,length));return data[start:start+length]
    f=RangeReader(len(data),fetch,block=1024,budget=len(data))
    assert f.read(1030)==data[:1030];f.seek(15);assert f.read(50)==data[15:65] and len(calls)==2
    f.seek(-17,2);a=bytearray(30);assert f.readinto(a)==17 and a[:17]==data[-17:]
    f.seek(len(data)+10);assert f.read(1)==b''
    try:f.read(2*1024**2)
    except ValueError:pass
    else:raise AssertionError('READ_CAP_NOT_ENFORCED')
    limited=RangeReader(len(data),fetch,block=1024,budget=1024)
    try:limited.read(1030)
    except MemoryError:pass
    else:raise AssertionError('BUDGET_CAP_NOT_ENFORCED')


if __name__=='__main__':synthetic();print('PASS_RANGE_CACHE_SEEK_READINTO_AND_CAPS')
