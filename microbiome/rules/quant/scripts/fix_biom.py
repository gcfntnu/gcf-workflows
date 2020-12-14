import sys
import re

fn = sys.argv[-1]

with open(fn, encoding='utf-8') as fh:
    txt = fh.read()

txt = txt.replace('[null]', 'null')
txt = txt.replace('"]', '"')       
txt = txt.replace('["', '"')       

sys.stdout.write(txt)

