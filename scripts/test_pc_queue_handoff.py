"""Packing handoff respects live locks and does not silently retry failures."""
import fcntl
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import complete_pc_screen as screen


class QueueHandoffTests(unittest.TestCase):
    def test_live_queue_is_never_duplicated(self):
        with tempfile.TemporaryDirectory() as temp,patch.object(screen,'ROOT',Path(temp)),patch.object(screen.subprocess,'Popen') as launch:
            base=Path(temp)/'runs/pc-dry-packing-v4';base.mkdir(parents=True)
            with (base/'queue.lock').open('a') as lock:
                fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
                self.assertEqual(screen.packing_ready_or_launch(),'running')
            launch.assert_not_called()

    def test_unreviewed_failure_stays_held(self):
        with tempfile.TemporaryDirectory() as temp,patch.object(screen,'ROOT',Path(temp)),patch.object(screen.subprocess,'Popen') as launch:
            case=Path(temp)/'runs/pc-dry-packing-v4/case';case.mkdir(parents=True)
            (case/'status.json').write_text(json.dumps({'state':'blocked_preserved'}))
            self.assertEqual(screen.packing_ready_or_launch(),'idle_completed_or_held')
            launch.assert_not_called()


if __name__=='__main__':unittest.main()
