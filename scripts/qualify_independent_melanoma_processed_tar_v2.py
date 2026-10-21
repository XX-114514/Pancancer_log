"""Same RAM archive/schema rules; demote HTTP-timeout peers, use 256KiB ranges."""
from pathlib import Path
import hashlib
import importlib.util
import json
import sys


def load(path, expected, name):
    assert hashlib.sha256(path.read_bytes()).hexdigest()==expected
    spec=importlib.util.spec_from_file_location(name,path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module


def main():
    folder=Path(__file__).parent
    core=load(folder/'qualify_independent_melanoma_processed_tar_v1.py',
        '3c5efec16e03040e588b9ef36fda4fe72ba9bde713fc8e83b0d20130a371804f','fixed_tar_schema_core')
    original_load=core.load
    def selected_load(path,expected,name):
        if path.name!='melanoma_IMC_verified_dns_transport_v2.py':
            return original_load(path,expected,name)
        assert expected=='12ffe4ac9c1c37efefbe42e8153928473da3cee8997a2086118d46ba81d7c25c'
        current=load(path.with_name('melanoma_IMC_verified_dns_transport_v3.py'),
            '3ab47b580b95a5517277737516ff372841b4b132bd1b3dbba98a03147c5221b6','same_host_HTTP_failover')
        ordinary=current.RangeFile
        class SmallRangeFile(ordinary):
            def __init__(self,*args,**kwargs):
                super().__init__(*args,**kwargs)
                self.max_chunk=256*1024
            def close(self):
                super().close()
                destination=Path(sys.argv[sys.argv.index('--output')+1])
                destination.joinpath('private_current_DNS_HTTP_transport_diagnostics.json').write_text(json.dumps(dict(
                    diagnostics=self.tls_diagnostics(),successful_bytes=self.transferred,
                    requested_bytes_including_failed=self.requested,no_addresses_or_payloads_written=True,
                    wrapper_SHA256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()),indent=2)+'\n')
        current.RangeFile=SmallRangeFile
        return current
    core.load=selected_load
    core.main()


if __name__=='__main__':main()
