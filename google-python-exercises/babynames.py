#!/usr/bin/python
# Copyright 2010 Google Inc.
# Licensed under the Apache License, Version 2.0
# http://www.apache.org/licenses/LICENSE-2.0

# Google's Python Class
# http://code.google.com/edu/languages/google-python-class/

import sys
import re

"""Baby Names exercise

Define the extract_names() function below and change main()
to call it.

For writing regex, it's nice to include a copy of the target
text for inspiration.

Here's what the html looks like in the baby.html files:
...
<h3 align="center">Popularity in 1990</h3>
....
<tr align="right"><td>1</td><td>Michael</td><td>Jessica</td>
<tr align="right"><td>2</td><td>Christopher</td><td>Ashley</td>
<tr align="right"><td>3</td><td>Matthew</td><td>Brittany</td>
...

Suggested milestones for incremental development:
 -Extract the year and print it
 -Extract the names and rank numbers and just print them
 -Get the names data into a dict and print it
 -Build the [year, 'name rank', ... ] list and print it
 -Fix main() to use the extract_names list
"""

def extract_names(filename):
  """
  Given a file name for baby.html, returns a list starting with the year string
  followed by the name-rank strings in alphabetical order.
  ['2006', 'Aaliyah 91', Aaron 57', 'Abagail 895', ' ...]
  """
  # +++your code here+++
  year=None
  name_ranks={}
  f = open(filename, 'r')
#  year = re.findall(r'Popularity in (\d+)', f.read())
#  name_ranks = re.findall(r'<tr align="right"><td>(\d+)</td><td>(\w+)</td><td>(\w+)</td>', f.read())
#  print 'year', year
#  print 'names', name_ranks
  for line in f:
    # search for year
    yearinline = re.search(r'Popularity in (\d+)' , line)
    if yearinline:
      year = yearinline.group(1)
    # search for name and rank, store them in dict
    this_name_rank = re.findall(r'<tr align="right"><td>(\d+)</td><td>(\w+)</td><td>(\w+)</td>' , line)
    if this_name_rank:
      for e in this_name_rank: # e.g. tuples [('4', 'Joshua', 'Madison')]
        name_ranks[e[1]]=e[0]
        name_ranks[e[2]]=e[0]

  # sort alphabetically
  name_ranks_sorted=[year]
  for name in sorted(name_ranks.keys()):
    name_ranks_sorted.append(name+' '+name_ranks[name])

  return name_ranks_sorted # return list ['2006', 'Aaliyah 91', Aaron 57', 'Abagail 895', ' ...] 


def main():
  # This command-line parsing code is provided.
  # Make a list of command line arguments, omitting the [0] element
  # which is the script itself.
  args = sys.argv[1:]

  if not args:
    print 'usage: [--summaryfile] file [file ...]'
    sys.exit(1)

  # Notice the summary flag and remove it from args if it is present.
  summary = False
  if args[0] == '--summaryfile':
    summary = True
    del args[0]

  # +++your code here+++
  # For each filename, get the names, then either print the text output
  # or write it to a summary file

  for thisfile in args:
    if summary:
      fout=open(thisfile+'.sum', 'w')  # open file for writing
      for e in extract_names(thisfile):
        print >> fout, e
    else:
      for e in extract_names(thisfile):
        print e
  
if __name__ == '__main__':
  main()
