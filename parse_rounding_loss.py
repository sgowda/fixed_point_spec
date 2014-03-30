#!/usr/bin/python
"""
parser for rounding loss log generated by matlab fixed-point library rounding functions
"""
import numpy as np
import re
import argparse
from collections import defaultdict

parser = argparse.ArgumentParser(description="process rounding loss log")
parser.add_argument('filename', default='/nas/users/sgowda/rounding_loss')
args = parser.parse_args()
filename = args.filename

rounding_loss = defaultdict(list)
file = open(filename, 'r')

pattern = '^(\w+): originally (\d+) bits, now (\d+) bits, loss = (.*?)$' 
regex = re.compile(pattern)
for line in file:
    m = regex.match(line)

    block         = m.group(1)
    orig_bitwidth = float(m.group(2))
    new_bitwidth  = float(m.group(3))
    loss          = float(m.group(4))

    op = (block, orig_bitwidth, new_bitwidth)
    rounding_loss[op].append(loss)

# separate rounding loss into 25bit keys and 35bit keys
keys_25bit = filter(lambda x: x[2] == 25, rounding_loss.keys())
keys_35bit = filter(lambda x: x[2] == 35, rounding_loss.keys())

def extract_loss(keys):
    loss = []
    for key in keys:
        loss += rounding_loss[key]
    return loss

loss_25bit = extract_loss(keys_25bit)
loss_35bit = extract_loss(keys_35bit)
