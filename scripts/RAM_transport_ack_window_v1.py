"""Bound queued SSH input using numeric receive acknowledgements only.

No payload is decoded or saved. This is a per-invocation transport window,
not a global proxy/SSH configuration change. Acknowledgements are technical
byte counts, and a launch receipt is distinct from scientific completion.
"""
import json


class Window:
    def __init__(self,window=262144):
        assert isinstance(window,int) and window>0
        self.window=window;self.acknowledged=0;self.launch=None

    def allowance(self,sent,total):
        assert isinstance(sent,int) and 0<=self.acknowledged<=sent<=total
        return max(0,min(total-sent,self.window-(sent-self.acknowledged)))

    def receive(self,line,sent,total):
        assert len(line)<=65536
        value=json.loads(line)
        if set(value)=={'transport_ack_received_bytes'}:
            count=value['transport_ack_received_bytes']
            assert type(count) is int and self.acknowledged<=count<=sent<=total
            self.acknowledged=count;return 'ACK'
        assert value.get('status','').startswith('STARTED_') and sent==total
        assert all(k in value for k in ['pid','startticks','started_UTC']) and self.launch is None
        self.launch=value;return 'SCIENTIFIC_LAUNCH_NOT_COMPLETION'
