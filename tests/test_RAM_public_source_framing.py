import importlib.util
import io
from pathlib import Path
import unittest

p=Path(__file__).resolve().parents[1]/'scripts/RAM_public_source_framing_v1.py'
s=importlib.util.spec_from_file_location('RAM_frames',p);m=importlib.util.module_from_spec(s);s.loader.exec_module(m)


class Frames(unittest.TestCase):
    def test_declared_completion_needs_no_EOF_and_next_frame_retained(self):
        a=b'synthetic_A';b=b'synthetic_B'
        stream=io.BytesIO(len(a).to_bytes(8,'big')+a+len(b).to_bytes(8,'big')+b+b'unread')
        self.assertEqual(m.read_frame(stream,100),a)
        self.assertEqual(m.read_frame(stream,100),b)
        self.assertEqual(stream.read(),b'unread')
    def test_partial_reads_and_truncation_or_oversize(self):
        class Partial(io.BytesIO):
            def read(self,size=-1):return super().read(min(size,3))
        self.assertEqual(m.read_exact(Partial(b'0123456789'),10),b'0123456789')
        with self.assertRaises(EOFError):m.read_exact(io.BytesIO(b'a'),2)
        with self.assertRaises(AssertionError):m.read_frame(io.BytesIO((101).to_bytes(8,'big')),100)


if __name__=='__main__':unittest.main()
