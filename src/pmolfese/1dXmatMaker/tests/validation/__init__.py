# Developed at the NIH
#
# Part of AFNI. See the repository LICENSE.txt for licensing terms and
# disclaimers (bundled with 1dXmatMaker as licenses/AFNI.txt).
# Third-party and adapted portions retain their existing attribution and license
# terms; see THIRD_PARTY.md.

"""Scientific test support, loaded only by tests/run.py.

ROOT is the standalone package directory, including when this tree is extracted
from a source release. Keeping it here avoids fragile cwd-dependent fixture paths.
No test runs or expected-file captures happen merely by importing the package.
"""
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
