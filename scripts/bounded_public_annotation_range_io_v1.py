"""Fixed-size verified HTTP ranges; compressed annotation bytes stay in RAM."""
import hashlib
import time
import urllib.request


def receive(url, object_bytes, wanted_bytes=None, progress=None):
    wanted = object_bytes if wanted_bytes is None else wanted_bytes
    assert 0 < wanted <= object_bytes <= 32 * 1024**2
    parts = []
    receipts = []
    validators = None
    for start in range(0, wanted, 1024**2):
        stop = min(start + 1024**2, wanted) - 1
        for attempt in range(3):
            try:
                request = urllib.request.Request(url, headers={'Range': 'bytes=%d-%d' % (start, stop), 'Accept-Encoding': 'identity'})
                with urllib.request.urlopen(request, timeout=30) as response:
                    assert response.status == 206
                    assert response.headers.get('Content-Range') == 'bytes %d-%d/%d' % (start, stop, object_bytes)
                    part = response.read(stop - start + 2)
                    assert len(part) == stop - start + 1
                    tags = (response.headers.get('ETag'), response.headers.get('Last-Modified'))
                    if validators is None:
                        validators = tags
                    else:
                        for a, b in zip(validators, tags):
                            assert a is None or b is None or a == b
                    receipt = dict(start=start, stop=stop, bytes=len(part), sha256=hashlib.sha256(part).hexdigest(), attempt=attempt + 1)
                break
            except Exception as error:
                if progress:
                    progress(dict(stage='RANGE_TRANSPORT_RETRY', completed_bytes=sum(len(p) for p in parts), current_start=start, attempt=attempt + 1, error_class=type(error).__name__))
                if attempt == 2:
                    raise RuntimeError('PUBLIC_ANNOTATION_RANGE_ATTEMPTS_EXHAUSTED') from None
                time.sleep(.5)
        parts.append(part)
        receipts.append(receipt)
        if progress:
            progress(dict(stage='VERIFIED_COMPRESSED_RANGES_IN_RAM', completed_bytes=sum(len(p) for p in parts), requested_bytes=wanted, completed_ranges=len(parts)))
    result = b''.join(parts)
    assert len(result) == wanted
    return result, dict(ETag=validators[0], Last_Modified=validators[1], verified_ranges=receipts)
