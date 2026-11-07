#!/usr/bin/env python3
# Developed at the NIH
#
# Part of AFNI. See the repository LICENSE.txt for licensing terms and
# disclaimers (bundled with 1dXmatMaker as licenses/AFNI.txt).
# Third-party and adapted portions retain their existing attribution and license
# terms; see THIRD_PARTY.md.

"""Create reproducible native binary and offline source archives, without builds.

The caller supplies the tested executable and its Rust target triple. This tool
never cross-compiles or certifies a platform. Sorted paths, fixed timestamps,
and normalized ownership make archives repeatable for identical input bytes.
Only the explicit source directories below are packaged: no target trees,
working datasets, credentials, or local build logs can enter the source archive.
"""
import argparse,gzip,hashlib,io,json,os,tarfile,tomllib
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]


def archive(output,files,created):
    """Write one new gzip/tar with stable metadata; executable bits are explicit."""
    with output.open('xb') as raw:
        # Register only a file we successfully created, and remember its inode
        # so cleanup cannot remove a replacement created by another process.
        created.append((output,os.fstat(raw.fileno()).st_ino))
        with gzip.GzipFile(filename='',mode='wb',fileobj=raw,mtime=0) as compressed:
            with tarfile.open(fileobj=compressed,mode='w') as tar:
                for name,(path,mode) in sorted(files.items()):
                    data=path.read_bytes();info=tarfile.TarInfo(name)
                    info.size=len(data);info.mode=mode;info.mtime=0
                    tar.addfile(info,io.BytesIO(data))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--binary',type=Path,required=True);p.add_argument('--target',required=True)
    p.add_argument('--output-dir',type=Path,required=True);a=p.parse_args()
    binary=a.binary.resolve(strict=True);out=a.output_dir.resolve();out.mkdir(parents=True,exist_ok=True)
    if not a.target or any(c not in 'abcdefghijklmnopqrstuvwxyz0123456789_-' for c in a.target):p.error('invalid target triple')
    version=tomllib.loads((ROOT/'Cargo.toml').read_text())['package']['version'];stem=f'1dXmatMaker-{version}'
    native=out/f'{stem}-{a.target}.tar.gz';source=out/f'{stem}-source.tar.gz';checksums=out/'SHA256SUMS.json'
    # Preflight all named outputs before creating any; never replace an earlier release.
    for path in [native,source,checksums]:
        if path.exists() or path.is_symlink():raise FileExistsError(path)
    source_files={}
    roots=['src','tests','examples','fixtures','data','licenses']
    paths=[ROOT/name for name in ['Cargo.toml','Cargo.lock','rust-toolchain.toml','Makefile','CMakeLists.txt','README.md','THIRD_PARTY.md']]
    for name in roots:
        paths.extend(p for p in (ROOT/name).rglob('*') if p.is_file() and '__pycache__' not in p.parts and p.suffix!='.pyc')
    for path in paths:
        if path.is_symlink():raise ValueError(f'symlinks are not source release inputs: {path}')
        source_files[f'{stem}/'+path.relative_to(ROOT).as_posix()]=(path,0o644)
    native_files={f'{stem}/bin/1dXmatMaker':(binary,0o755)}
    for path in paths:
        relative=path.relative_to(ROOT)
        # Binary users already have the waveform bank inside the executable.
        # Ship its license/provenance, but keep scientific fixtures source-only.
        if relative.parts[0]=='licenses' or (relative.parts[0]=='data' and path.name in ['LICENSE','provenance.json']) or (relative.parts[0]=='examples' and path.suffix in ['.tsv','.1D']) or (path.parent==ROOT and path.name in ['README.md','THIRD_PARTY.md']):
            native_files[f'{stem}/share/1dXmatMaker/{relative.as_posix()}']=(path,0o644)
    native_files[f'{stem}/share/1dXmatMaker/licenses/GLMsingle.txt']=(ROOT/'data/glmsingle/LICENSE',0o644)
    created=[]
    try:
        for output,files in [(native,native_files),(source,source_files)]:
            archive(output,files,created)
        hashes={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p,_ in created}
        with checksums.open('x') as f:
            created.append((checksums,os.fstat(f.fileno()).st_ino))
            json.dump(hashes,f,indent=2);f.write('\n')
    except BaseException:
        for path,inode in created:
            if path.exists() and path.lstat().st_ino==inode:path.unlink()
        raise
    print(json.dumps(hashes,indent=2))
if __name__=='__main__':main()
