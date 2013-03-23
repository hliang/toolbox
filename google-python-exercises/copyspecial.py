#!/usr/bin/python
# Copyright 2010 Google Inc.
# Licensed under the Apache License, Version 2.0
# http://www.apache.org/licenses/LICENSE-2.0

# Google's Python Class
# http://code.google.com/edu/languages/google-python-class/

import sys
import re
import os
import shutil
import commands

"""Copy Special exercise
"""

# +++your code here+++
# Write functions and modify main() to call them

def get_special_paths(dirs):
  '''input a list of dirr, output a list of special files
  special files have pattern '__\w+__' in file names
  '''
  special_paths = [ ]
  for dir in dirs:
    fnames = os.listdir(dir)
    for fname in fnames:
      if re.search(r'__\w+__', fname):
        special_paths.append( os.path.abspath(os.path.join(dir, fname)) )
  return special_paths

def copyto(path, dir):
  '''copy a LIST of files to destination dir
  source: path
  target: dir
  '''
  if not os.path.exists(dir):
    os.mkdir(dir)
  for f in path:
    # dst_path = os.path.join(dir, os.path.basename(f))  # not necessary
    # shutil.copy(f, dst_path )  # not necessary, because dst_path can be /dir/ or /dir/file
    shutil.copy(f, dir )

def zipto(paths, zippath):
  '''create a zipfile 'zippath' containing the files in the list 'paths'
  '''
  cmd = 'zip -j ' + zippath 
  for path in paths:
    cmd = cmd + ' ' + path
  print "## Command I'm going to do:\n##", cmd

  (status, output) = commands.getstatusoutput(cmd)
  if status:
    print status
    print output
    # or write to stderr:
    # sys.stderr.write(output)
    sys.exit(1)
  else:
    print '## Completed without error'

def main():
  # This basic command line argument parsing code is provided.
  # Add code to call your functions below.

  # Make a list of command line arguments, omitting the [0] element
  # which is the script itself.
  args = sys.argv[1:]
  if not args:
    print "usage: [--todir dir][--tozip zipfile] dir [dir ...]";
    sys.exit(1)

  # todir and tozip are either set from command line
  # or left as the empty string.
  # The args array is left just containing the dirs.
  todir = ''
  if args[0] == '--todir':
    todir = args[1]
    del args[0:2]

  tozip = ''
  if args[0] == '--tozip':
    tozip = args[1]
    del args[0:2]

  if len(args) == 0:
    print "error: must specify one or more dirs"
    sys.exit(1)

  # +++your code here+++
  # Call your functions
  to_check_dirs = args
  special_paths = get_special_paths(to_check_dirs)

  if todir != '':
    copyto(special_paths, todir)
  elif tozip != '':
    zipto(special_paths, tozip)
  else:
    print "\n".join(special_paths)
  
if __name__ == "__main__":
  main()
