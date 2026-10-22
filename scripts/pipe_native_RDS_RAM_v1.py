"""Stream a verified compressed RAM buffer to an R binary stdin connection.

No native or expanded payload files are created. The caller owns the deadline
and input integrity checks. Compressed file names and individual data are never
part of the returned technical receipt.
"""
import bz2
import gzip
import lzma
import subprocess


def pipe_native_RDS(buffer, command, *, maximum_expanded_bytes, timeout, **kwargs):
    buffer.seek(0)
    magic = buffer.read(6)
    buffer.seek(0)
    if magic.startswith(b'\x1f\x8b'):
        codec, stream = 'gzip', gzip.GzipFile(fileobj=buffer, mode='rb')
    elif magic.startswith(b'BZh'):
        codec, stream = 'bzip2', bz2.BZ2File(buffer, mode='rb')
    elif magic.startswith(b'\xfd7zXZ\x00'):
        codec, stream = 'xz', lzma.LZMAFile(buffer, mode='rb')
    else:
        codec, stream = 'none', buffer
    child = subprocess.Popen(command, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                             stderr=subprocess.DEVNULL, **kwargs)
    expanded = 0
    try:
        while True:
            chunk = stream.read(1024**2)
            if not chunk:
                break
            expanded += len(chunk)
            if expanded > maximum_expanded_bytes:
                raise ValueError('EXPANDED_STREAM_RESOURCE_CAP')
            child.stdin.write(chunk)
        child.stdin.close()
        child.stdin = None
        stdout, _ = child.communicate(timeout=timeout)
        if len(stdout) > 1024**2:
            raise ValueError('AGGREGATE_OUTPUT_RESOURCE_CAP')
        return subprocess.CompletedProcess(command, child.returncode, stdout), {
            'codec': codec, 'expanded_stream_bytes': expanded,
            'no_native_or_expanded_payload_files': True}
    finally:
        if child.poll() is None:
            child.kill()
            child.wait(timeout=15)
        if stream is not buffer:
            stream.close()
