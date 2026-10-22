"""Exercise the changed input-size boundary without native bytes or network."""
from pathlib import Path
from types import SimpleNamespace


def check_guard(code, size):
    ns = {}
    exec(code, ns)
    calls, allocated = [], []
    class VirtualOutput:
        def seek(self, position):
            allocated.append(position)
        def write(self, value):
            assert value == b'\0'
        def close(self):
            pass
    class ReachedRangeStage(Exception):
        pass
    class NoNetwork:
        def open(self, request, timeout):
            calls.append(request.get_header('Range'))
            raise ReachedRangeStage()
    ns['io'] = SimpleNamespace(BytesIO=VirtualOutput)
    try:
        ns['fetch_ranges']('https://zenodo.org/records/0/files/SYNTHETIC',
            expected_bytes=size, expected_MD5='0'*32, chunk_bytes=1024**2,
            workers=1, request_timeout=1, maximum_range_seconds=2, opener_factory=NoNetwork)
    except (AssertionError, ReachedRangeStage) as error:
        return type(error).__name__, calls, allocated
    raise AssertionError('SYNTHETIC_NO_NETWORK_PROBE_DID_NOT_STOP')


def main():
    root = Path(__file__).parent
    old = (root/'fetch_public_RDS_ranges_RAM_v1.py').read_text()
    new = (root/'fetch_public_RDS_ranges_RAM_v2.py').read_text()
    # Only the bounded byte/range-count guards change; HTTP/hash logic is equal.
    old_code = '\n'.join(old.splitlines()[1:])
    new_code = '\n'.join(new.splitlines()[1:])
    expected = old_code.replace('0 < expected_bytes <= 1024**3', '0 < expected_bytes <= 2*1024**3')
    expected = expected.replace('len(ranges) <= 1024', 'len(ranges) <= 2048')
    assert expected == new_code
    size = 1024**3+12345
    result, calls, allocated = check_guard(old, size)
    assert result == 'AssertionError' and not calls and not allocated
    result, calls, allocated = check_guard(new, size)
    assert result == 'ReachedRangeStage' and calls and allocated == [size-1]
    for size in [0, -1, 2*1024**3+1]:
        result, calls, allocated = check_guard(new, size)
        assert result == 'AssertionError' and not calls and not allocated
    result, calls, allocated = check_guard(new, 2*1024**3)
    assert result == 'ReachedRangeStage' and calls and allocated == [2*1024**3-1]
    print('{"status":"PASS_2GiB_GUARD_BOUNDARIES_NO_NETWORK_NO_NATIVE_ALLOCATION"}')


if __name__ == '__main__':
    main()
