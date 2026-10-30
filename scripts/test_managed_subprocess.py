"""A timed-out shell must not leave a running child behind."""
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import time
import unittest
from managed_subprocess import run_group


class ProcessTreeTests(unittest.TestCase):
    def test_timeout_stops_term_ignoring_child(self):
        with tempfile.TemporaryDirectory() as name:
            root=Path(name)
            (root/'child.py').write_text('import os,signal,time\nfrom pathlib import Path\n'
                'signal.signal(signal.SIGTERM,signal.SIG_IGN)\nPath("child.pid").write_text(str(os.getpid()))\ntime.sleep(60)\n')
            with (root/'log').open('w') as log:
                with self.assertRaises(subprocess.TimeoutExpired):
                    run_group(['/bin/sh','-c', '"'+sys.executable+'" child.py & wait'],cwd=root,
                              env=os.environ.copy(),log=log,timeout=.5)
            pid=int((root/'child.pid').read_text())
            for _ in range(30):
                stat=Path(f'/proc/{pid}/stat')
                if not stat.exists() or stat.read_text().split()[2]=='Z':
                    break
                time.sleep(.01)
            else:
                self.fail('Timed-out descendant remains running')


if __name__=='__main__':
    unittest.main()
