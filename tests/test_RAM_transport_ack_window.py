from pathlib import Path
import sys,json,unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from RAM_transport_ack_window_v1 import Window


class AckWindow(unittest.TestCase):
    def test_upload_cannot_queue_full_payload_before_receiving_ack(self):
        flow=Window(128);total=1500;sent=0
        for unused in range(4):sent+=min(32,flow.allowance(sent,total))
        self.assertEqual(sent,128);self.assertEqual(flow.allowance(sent,total),0)
        flow.receive(json.dumps({'transport_ack_received_bytes':64}),sent,total)
        self.assertEqual(flow.allowance(sent,total),64)

    def test_wrong_future_regressing_and_noninteger_ack_rejected(self):
        for count in [129,1.2,True]:
            with self.assertRaises(AssertionError):Window(128).receive(json.dumps({'transport_ack_received_bytes':count}),128,1000)
        flow=Window(128);flow.receive(json.dumps({'transport_ack_received_bytes':100}),128,1000)
        with self.assertRaises(AssertionError):flow.receive(json.dumps({'transport_ack_received_bytes':90}),128,1000)

    def test_complete_partial_ack_sequence_and_distinct_launch(self):
        flow=Window(128);total=1000;sent=0
        while sent<total:
            sent+=min(32,flow.allowance(sent,total))
            if sent%64==0:flow.receive(json.dumps({'transport_ack_received_bytes':sent}),sent,total)
            self.assertLessEqual(sent-flow.acknowledged,128)
        receipt=dict(status='STARTED_SYNTHETIC_JOB',pid=123,startticks='456',started_UTC='synthetic')
        self.assertEqual(flow.receive(json.dumps(receipt),sent,total),'SCIENTIFIC_LAUNCH_NOT_COMPLETION')
        self.assertEqual(flow.launch,receipt)
        with self.assertRaises(AssertionError):flow.receive(json.dumps(receipt),sent,total)


if __name__=='__main__':unittest.main()
