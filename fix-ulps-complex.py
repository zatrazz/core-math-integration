#!/usr/bin/env python3

# usage: fix-ulps.py atan2 float libm-test-ulps

import sys
import tempfile
import shutil

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

def isfunc(func, func2rem):
  f = func.split('_')
  if len(f) >= 2:
    return func2rem == f[0] and f[1] in ['downward', 'towardzero', 'upward']
  else:
    return func2rem == f[0]


class temporary_copy(object):
  def __init__(self, original_path):
    self.original_path = original_path

  def __enter__(self):
    self.temp = tempfile.NamedTemporaryFile()
    return self.temp

  def __exit__(self, exc_type, exc_val, exc_tb):
    self.temp.flush()
    shutil.copy2 (self.temp.name, self.original_path)


if __name__ == "__main__":
  func2rem = sys.argv[1]
  ftype2rem = sys.argv[2]
  fpart2rem = sys.argv[3]

  for arg in sys.argv[4:]:
    with open(arg, "rb") as fr, temporary_copy(arg) as fw:
      func = ''
      ftype = ''
      fpart = ''
      for line in fr:
        dline = line.decode()
        if not dline.startswith("#") and not dline in ['\n', '\r\n']:
          if dline.startswith("Function: "):
            tofilter = False
            fields = dline.split()
            if len(fields) == 5 and fpart2rem == fields[1]:
              func = fields[4].replace('\"','').replace(':', '')
              tofilter = True
          else:
            ftype = dline.split()[0].replace(':', '')
            if isfunc(func, func2rem) and ftype2rem == ftype and tofilter:
              continue
        #eprint(line)
        fw.write(line)
