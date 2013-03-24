#!/usr/bin/python
# Copyright 2010 Google Inc.
# Licensed under the Apache License, Version 2.0
# http://www.apache.org/licenses/LICENSE-2.0

# Google's Python Class
# http://code.google.com/edu/languages/google-python-class/

import os
import re
import sys
import urllib

"""Logpuzzle exercise
Given an apache logfile, find the puzzle urls and download the images.

Here's what a puzzle url looks like:
10.254.254.28 - - [06/Aug/2007:00:13:48 -0700] "GET /~foo/puzzle-bar-aaab.jpg HTTP/1.0" 302 528 "-" "Mozilla/5.0 (Windows; U; Windows NT 5.1; en-US; rv:1.8.1.6) Gecko/20070725 Firefox/2.0.0.6"
"""

def sec_wd (s):
  """returns the second word if the string ends in the pattern "-wordchars-wordchars.jpg" """
  patn = re.search(r'-\w+-(\w+)\.jpg', s)
  if patn:
    return patn.group(1)
  else:
    return s

def read_urls(filename):
  """Returns a list of the puzzle urls from the given log file,
  extracting the hostname from the filename itself.
  Screens out duplicate urls and returns the urls sorted into
  increasing order."""
  # +++your code here+++
  us_index = filename.find('_')
  domain = filename[us_index+1:]

  urls = []
  logfile = open(filename, 'r')
  for line in logfile:
    path = re.search(r'GET (\S+\.jpg) HTTP', line, re.IGNORECASE) 
    if path and 'puzzle' in path.group(1): # path found
      newurl = 'http://' + domain + path.group(1)
      if not newurl in urls:  # no duplicates
        urls.append(newurl )  # add to url list
  logfile.close()

  # urls = sorted(urls)  # sorted into alphabetical order
  urls = sorted(urls, key=sec_wd)
  return urls

def download_images(img_urls, dest_dir):
  """Given the urls already in the correct order, downloads
  each image into the given directory.
  Gives the images local filenames img0, img1, and so on.
  Creates an index.html in the directory
  with an img tag to show each local image file.
  Creates the directory if necessary.
  """
  # +++your code here+++
  # create dir if it doesn't exist
  if not os.path.exists(dest_dir):
    os.mkdir(dest_dir)

  # download and store img in dest_dir
  dest_img_paths = [ ]
  for i in range(len(img_urls)):
    dest_img_path = os.path.join(dest_dir, 'img'+str(i)+'.jpg')
    dest_img_paths.append('img'+str(i)+'.jpg')
    print 'Retrieving ', img_urls[i], '... saving as:', dest_img_path
    urllib.urlretrieve(img_urls[i], dest_img_path)

  # create index.html to display images
  index_html_path = os.path.join(dest_dir, 'index.html')
  index_html = open(index_html_path, 'w')
  print >> index_html, r'''
<verbatim>
<html>
<body>'''

  img_src = ''
  for i in range(len(dest_img_paths)):
    img_src = img_src + r'<img src="' + dest_img_paths[i] + r'">'
  print >> index_html, img_src
  
  print >> index_html, r'''
</body>
</html>'''

  index_html.close()

def main():
  args = sys.argv[1:]

  if not args:
    print 'usage: [--todir dir] logfile '
    sys.exit(1)

  todir = ''
  if args[0] == '--todir':
    todir = args[1]
    del args[0:2]

  img_urls = read_urls(args[0])

  if todir:
    download_images(img_urls, todir)
  else:
    print '\n'.join(img_urls)

if __name__ == '__main__':
  main()
