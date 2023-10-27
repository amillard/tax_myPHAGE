#!/usr/bin/env python3
# Created: Fri Oct 27 13:17:42 2023
# Last changed: Time-stamp: <Last changed 2023-10-27 14:45:52 by Thomas Sicheritz, thomas>

import sys
import re
from pathlib import Path

sys.path.insert(0, Path(__file__).parent.parent.parent.as_posix())
from taxmyphage.__main__ import main

if __name__ == "__main__":
    sys.argv[0] = re.sub(r"(-script\.pyw|\.exe)?$", "", sys.argv[0])
    sys.exit(main())
