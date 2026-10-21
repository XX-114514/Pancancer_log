"""Length-framed RAM input; completion does not require an SSH EOF signal."""


def read_exact(stream,size):
    assert isinstance(size,int) and 0<=size<=20*1024**2
    chunks=[];remaining=size
    while remaining:
        block=stream.read(min(65536,remaining))
        if not block:raise EOFError('incomplete_declared_frame')
        chunks.append(block);remaining-=len(block)
    return b''.join(chunks)


def read_frame(stream,maximum):
    length=int.from_bytes(read_exact(stream,8),'big')
    assert 0<length<=maximum
    return read_exact(stream,length)
