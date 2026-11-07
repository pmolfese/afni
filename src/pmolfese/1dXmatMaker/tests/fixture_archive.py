#!/usr/bin/env python3
# Developed at the NIH
#
# Part of AFNI. See the repository LICENSE.txt for licensing terms and
# disclaimers (bundled with 1dXmatMaker as licenses/AFNI.txt).
# Third-party and adapted portions retain their existing attribution and license
# terms; see THIRD_PARTY.md.

"""Pack, verify and unpack frozen fixtures using only the Python standard library.

Packing writes a new archive/manifest pair for review; it never captures reference
values or replaces the distributed pair. Extraction checks every member before
publishing a complete directory. The JSON inventory also verifies cached copies.
"""
import argparse
from datetime import date
import gzip
import hashlib
import io
import json
from pathlib import Path, PurePosixPath
import tarfile
import tempfile

ROOT = Path(__file__).resolve().parents[1]
ARCHIVE = 'fixtures.tgz'
MANIFEST = 'fixtures.manifest.json'
GROUPS = {
    'afni/baselines': 'Polynomial baselines, sampled signals and nuisance regressors.',
    'afni/basis': 'Response waveforms, expressions and function coverage.',
    'afni/compatibility': 'Shared AFNI compatibility cases and executable/capture provenance.',
    'afni/events': 'Event timing, event columns and event identities.',
    'afni/semantics': 'Modulation, censoring, contrasts, selectors and automatic polort.',
    'glmsingle': 'Expected sampled HRFs, upstream reference functions and their license.',
}


def sha256(data):
    return hashlib.sha256(data).hexdigest()


def read_manifest(package):
    """Load the inventory, rejecting ambiguous or unsafe relative member names."""
    manifest = json.loads((package / MANIFEST).read_text())
    if manifest['schema_version'] != 1 or manifest['archive']['name'] != ARCHIVE:
        raise ValueError('unsupported fixture manifest')
    digest = manifest['archive']['sha256']
    if len(digest) != 64 or any(c not in '0123456789abcdef' for c in digest):
        raise ValueError('invalid archive SHA-256')
    names = set()
    for entry in manifest['files']:
        name = entry['path']
        path = PurePosixPath(name)
        if (not name or path.is_absolute() or '..' in path.parts or
                path.as_posix() != name or name == '.' or '\\' in name or name in names):
            raise ValueError(f'invalid or duplicate fixture path: {name}')
        names.add(name)
    if len(names) != manifest['file_count']:
        raise ValueError('fixture count differs from manifest')
    if sum(entry['size_bytes'] for entry in manifest['files']) != manifest['uncompressed_bytes']:
        raise ValueError('fixture sizes differ from manifest')
    return manifest


def fixture_directory(package=ROOT):
    """Return the cache location without creating files, including for --help."""
    return package / 'target/fixtures' / read_manifest(package)['archive']['sha256']


def verify_tree(directory, manifest):
    """Require the exact regular-file inventory and byte hashes; reject symlinks."""
    entries = {entry['path']: entry for entry in manifest['files']}
    found = set()
    if directory.is_symlink():
        raise ValueError(f'fixture directory is a symlink: {directory}')
    for path in directory.rglob('*'):
        if path.is_symlink():
            raise ValueError(f'fixture is a symlink: {path}')
        if path.is_dir():
            continue
        name = path.relative_to(directory).as_posix()
        if name not in entries or not path.is_file():
            raise ValueError(f'unexpected fixture: {path}')
        data = path.read_bytes()
        entry = entries[name]
        if len(data) != entry['size_bytes'] or sha256(data) != entry['sha256']:
            raise ValueError(f'fixture checksum mismatch: {path}')
        found.add(name)
    if found != entries.keys():
        raise ValueError(f'missing fixtures in {directory}: {sorted(entries.keys() - found)}')


def unpack(package=ROOT, destination=None):
    """Verify the archive and unpack once, or verify all bytes in an existing cache.

    Only listed regular files are accepted; no archive-supplied paths or links
    can escape the destination. A temporary sibling avoids partial caches and
    allows concurrent test runners to publish the same content safely.
    """
    manifest = read_manifest(package)
    data = (package / ARCHIVE).read_bytes()
    if len(data) != manifest['archive']['size_bytes'] or sha256(data) != manifest['archive']['sha256']:
        raise ValueError('fixture archive checksum mismatch')
    if destination is None:
        destination = fixture_directory(package)
    if destination.exists() or destination.is_symlink():
        verify_tree(destination, manifest)
        return destination
    destination.parent.mkdir(parents=True, exist_ok=True)
    entries = {entry['path']: entry for entry in manifest['files']}
    with tempfile.TemporaryDirectory(prefix='.fixtures-', dir=destination.parent) as temporary:
        tree = Path(temporary) / 'contents'
        tree.mkdir()
        seen = set()
        with tarfile.open(fileobj=io.BytesIO(data), mode='r:gz') as archive:
            for member in archive:
                if member.name not in entries or member.name in seen or not member.isfile():
                    raise ValueError(f'unexpected archive member: {member.name}')
                entry = entries[member.name]
                if member.size != entry['size_bytes']:
                    raise ValueError(f'fixture size mismatch: {member.name}')
                content = archive.extractfile(member).read()
                if sha256(content) != entry['sha256']:
                    raise ValueError(f'fixture checksum mismatch: {member.name}')
                path = tree / member.name
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_bytes(content)
                seen.add(member.name)
        if seen != entries.keys():
            raise ValueError('archive is missing manifest entries')
        try:
            tree.rename(destination)
        except OSError:
            # Another runner may have published this same complete cache.
            if not destination.exists():
                raise
            verify_tree(destination, manifest)
    return destination


def pack(source, output, updated, note, previous=None, source_revision=None):
    """Create a deterministic tgz plus a dated, per-file SHA-256 JSON inventory.

    Dates record when bytes were first inventoried or changed, not inferred
    capture dates. Existing scientific provenance stays inside the archive.
    Supply the previous manifest to retain unchanged-file dates and update notes.
    Output must be a new directory, so the distributed pair is never overwritten.
    """
    date.fromisoformat(updated)
    source = source.resolve(strict=True)
    output = output.resolve()
    if output == source or source in output.parents:
        raise ValueError('archive output must be outside the fixture input tree')
    old = json.loads(previous.read_text()) if previous else {}
    old_entries = {entry['path']: entry for entry in old.get('files', [])}
    entries = []
    paths = sorted(source.rglob('*'))
    for path in paths:
        if path.is_symlink() or not (path.is_dir() or path.is_file()):
            raise ValueError(f'not a regular fixture file/directory: {path}')
    files = [path for path in paths if path.is_file()]
    if not files:
        raise ValueError('fixture input is empty')
    output.mkdir(parents=True, exist_ok=False)
    with (output / ARCHIVE).open('xb') as raw:
        with gzip.GzipFile(filename='', fileobj=raw, mode='wb', mtime=0, compresslevel=9) as compressed:
            with tarfile.open(fileobj=compressed, mode='w', format=tarfile.USTAR_FORMAT) as archive:
                for path in files:
                    name = path.relative_to(source).as_posix()
                    data = path.read_bytes()
                    digest = sha256(data)
                    prior = old_entries.get(name, {})
                    last_updated = prior.get('last_updated', updated) if prior.get('sha256') == digest else updated
                    entries.append(dict(path=name, size_bytes=len(data), sha256=digest, last_updated=last_updated))
                    info = tarfile.TarInfo(name)
                    info.size = len(data)
                    info.mode = 0o644
                    info.mtime = 0
                    archive.addfile(info, io.BytesIO(data))
    archive_bytes = (output / ARCHIVE).read_bytes()
    names = {entry['path'] for entry in entries}
    update = dict(date=updated, note=note,
                  added=sum(entry['path'] not in old_entries for entry in entries),
                  changed=sum(entry['path'] in old_entries and entry['sha256'] != old_entries[entry['path']]['sha256'] for entry in entries),
                  removed=len(old_entries.keys() - names))
    if source_revision:
        update['source_revision'] = source_revision
    groups = []
    for prefix, description in GROUPS.items():
        selected = [entry for entry in entries if entry['path'].startswith(prefix + '/')]
        groups.append(dict(path=prefix, description=description, file_count=len(selected),
                           uncompressed_bytes=sum(entry['size_bytes'] for entry in selected)))
    manifest = dict(schema_version=1, description='Frozen independent scientific regression fixtures for 1dXmatMaker.',
                    last_updated=updated,
                    date_semantics='Dates record first inventory or subsequent byte changes; original capture provenance is preserved inside the archive.',
                    archive=dict(name=ARCHIVE, sha256=sha256(archive_bytes), size_bytes=len(archive_bytes)),
                    file_count=len(entries), uncompressed_bytes=sum(entry['size_bytes'] for entry in entries),
                    groups=groups, updates=old.get('updates', []) + [update], files=entries)
    (output / MANIFEST).write_text(json.dumps(manifest, indent=2) + '\n')
    # Validate the serialized inventory against the source before handing it off.
    verify_tree(source, read_manifest(output))
    return manifest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest='command', required=True)
    extract = commands.add_parser('unpack', help='verify and unpack the distributed fixtures')
    extract.add_argument('--output-dir', type=Path, help='default: content-addressed cache under target/fixtures')
    create = commands.add_parser('pack', help='create a NEW archive/manifest pair for review')
    create.add_argument('--input', type=Path, required=True)
    create.add_argument('--output-dir', type=Path, required=True)
    create.add_argument('--date', required=True, help='update date, YYYY-MM-DD')
    create.add_argument('--note', required=True, help='describe the reviewed change and its reason')
    create.add_argument('--source-revision', help='optional revision identifying the source fixtures')
    create.add_argument('--previous-manifest', type=Path, default=ROOT / MANIFEST)
    args = parser.parse_args()
    if args.command == 'unpack':
        print(unpack(destination=args.output_dir))
    else:
        manifest = pack(args.input, args.output_dir, args.date, args.note,
                        args.previous_manifest, args.source_revision)
        print(f"Packed {manifest['file_count']} files into {args.output_dir}")


if __name__ == '__main__':
    main()
