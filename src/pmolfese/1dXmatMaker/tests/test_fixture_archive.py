# Developed at the NIH
#
# Part of AFNI. See the repository LICENSE.txt for licensing terms and
# disclaimers (bundled with 1dXmatMaker as licenses/AFNI.txt).
# Third-party and adapted portions retain their existing attribution and license
# terms; see THIRD_PARTY.md.

"""Archive integrity, reproducibility, update tracking and relocated runner checks."""
import io
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tarfile
import tempfile
import unittest
from concurrent.futures import ThreadPoolExecutor

from fixture_archive import ARCHIVE, MANIFEST, ROOT, pack, sha256, unpack


class FixtureArchiveTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.source = self.root / 'source'
        (self.source / 'afni/basis').mkdir(parents=True)
        (self.source / 'afni/basis/matrix.1D').write_bytes(b'0 1\n2 3\n')
        (self.source / 'LICENSE').write_bytes(b'Preserve this notice.\n')
        self.package = self.root / 'package'
        self.manifest = pack(self.source, self.package, '2026-11-07', 'Initial inventory')

    def save_manifest(self):
        (self.package / MANIFEST).write_text(json.dumps(self.manifest))

    def test_round_trip_and_deterministic_repack(self):
        destination = unpack(self.package)
        for entry in self.manifest['files']:
            self.assertEqual((self.source / entry['path']).read_bytes(),
                             (destination / entry['path']).read_bytes())
        second = self.root / 'second'
        repacked = pack(destination, second, '2026-11-08', 'Repacked unchanged', self.package / MANIFEST)
        self.assertEqual((self.package / ARCHIVE).read_bytes(), (second / ARCHIVE).read_bytes())
        self.assertTrue(all(entry['last_updated'] == '2026-11-07' for entry in repacked['files']))
        self.assertEqual(repacked['updates'][-1]['changed'], 0)
        self.assertEqual(len(repacked['updates']), 2)
        self.assertEqual(unpack(self.package), destination)

    def test_update_dates_and_counts(self):
        (self.source / 'afni/basis/matrix.1D').write_bytes(b'0 2\n4 6\n')
        (self.source / 'LICENSE').unlink()
        (self.source / 'new.txt').write_bytes(b'new')
        updated = pack(self.source, self.root / 'updated', '2026-11-08', 'Reviewed update', self.package / MANIFEST)
        self.assertEqual([entry['last_updated'] for entry in updated['files']], ['2026-11-08'] * 2)
        self.assertEqual([updated['updates'][-1][key] for key in ['added', 'changed', 'removed']], [1, 1, 1])

    def test_archive_corruption(self):
        with (self.package / ARCHIVE).open('ab') as stream:
            stream.write(b'corruption')
        with self.assertRaisesRegex(ValueError, 'archive checksum mismatch'):
            unpack(self.package)
        self.assertFalse((self.package / 'target').exists())

    def test_cache_corruption_missing_extra_and_symlink(self):
        for kind in ['modified', 'missing', 'extra', 'symlink']:
            with self.subTest(kind=kind):
                destination = unpack(self.package, self.root / kind)
                target = destination / 'LICENSE'
                if kind == 'modified':
                    target.write_bytes(b'altered')
                elif kind == 'missing':
                    target.unlink()
                elif kind == 'extra':
                    (destination / 'extra').touch()
                else:
                    target.unlink()
                    target.symlink_to(self.source / 'LICENSE')
                with self.assertRaises(ValueError):
                    unpack(self.package, destination)

    def test_manifest_path_traversal(self):
        self.manifest['files'][0]['path'] = '../outside'
        self.save_manifest()
        with self.assertRaisesRegex(ValueError, 'invalid or duplicate fixture path'):
            unpack(self.package)

    def test_member_hash_and_type_checked_even_with_valid_archive_hash(self):
        for kind in ['bytes', 'symlink', 'duplicate', 'unlisted', 'missing']:
            with self.subTest(kind=kind):
                archive_path = self.package / ARCHIVE
                with tarfile.open(archive_path, 'w:gz') as archive:
                    for entry in self.manifest['files']:
                        if kind == 'missing':
                            continue
                        data = (self.source / entry['path']).read_bytes()
                        member = tarfile.TarInfo(entry['path'])
                        member.size = len(data)
                        if kind == 'bytes':
                            data = b'x' * len(data)
                        elif kind == 'symlink':
                            member.type = tarfile.SYMTYPE
                            member.linkname = '../outside'
                        elif kind == 'unlisted':
                            member.name = '../outside'
                        archive.addfile(member, io.BytesIO(data))
                        if kind == 'duplicate':
                            archive.addfile(member, io.BytesIO(data))
                data = archive_path.read_bytes()
                self.manifest['archive'].update(sha256=sha256(data), size_bytes=len(data))
                self.save_manifest()
                destination = self.root / kind
                with self.assertRaises(ValueError):
                    unpack(self.package, destination)
                self.assertFalse(destination.exists())
                self.assertFalse((self.root / 'outside').exists())

    def test_concurrent_unpack(self):
        with ThreadPoolExecutor(max_workers=2) as pool:
            paths = list(pool.map(lambda _: unpack(self.package), range(2)))
        self.assertEqual(paths[0], paths[1])
        self.assertEqual(unpack(self.package), paths[0])

    def test_relocated_runner_help_and_full_corpus_audit(self):
        package = self.root / 'relocated'
        shutil.copytree(ROOT / 'tests', package / 'tests', ignore=shutil.ignore_patterns('__pycache__'))
        for name in [ARCHIVE, MANIFEST]:
            shutil.copyfile(ROOT / name, package / name)
        runner = package / 'tests/run.py'
        for args in [['--help'], ['reference', '--help']]:
            subprocess.run([sys.executable, str(runner), *args], cwd=self.root,
                           check=True, capture_output=True, text=True)
        self.assertFalse((package / 'target').exists())
        result = subprocess.run([sys.executable, str(runner), 'reference', 'audit'],
                                cwd=self.root, check=True, capture_output=True, text=True)
        self.assertIn('PASS:', result.stdout)
        self.assertTrue((package / 'target/fixtures').is_dir())


if __name__ == '__main__':
    unittest.main()
