from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import os
import sys

if len(sys.argv) != 2:
  print('Error: usage is {} <input file>'.format(sys.argv[0]))
  sys.exit(1)

if not os.path.isfile(sys.argv[1]):
  print('Error: specified input file is not a file')
  sys.exit(1)

text = open(sys.argv[1], 'r').read().splitlines()

if len(text) != 257 or text[0] != 'FD ACL,FD AHD,TR ACL,TR AHD':
  print('Error: specified file not formatted like output of `automata.cpp`')
  sys.exit(1)

data = []
for line in text[1:]:
  row = line.split(',')
  if len(row) != 4:
    print('Error: specified file not formatted like output of `automata.cpp`')
    sys.exit(1)
  try:
    data.append(tuple(map(float, row)))
  except ValueError:
    print('Error: rows of specified file do not contain floating points')
    sys.exit(1)

num_fd_more_ordered = 0
num_tr_more_ordered = 0
num_inconclusive = 0

for row in data:
  if row[0] < row[2] and row[1] <= row[3]:
    num_fd_more_ordered += 1
  elif row[2] > row[0] and row[3] >= row[1]:
    num_tr_more_ordered += 1
  else:
    num_inconclusive += 1

print('Number of more ordered  deterministic  automata:', num_fd_more_ordered)
print('Number of more ordered time-reversible automata:', num_tr_more_ordered)
assert num_inconclusive == 256 - num_fd_more_ordered - num_tr_more_ordered
print('Number of inconclusive results:', num_inconclusive)