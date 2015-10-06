#!/usr/bin/env python

import json
import sys

data = json.load(sys.stdin)
json.dump(data,sys.stdout)
