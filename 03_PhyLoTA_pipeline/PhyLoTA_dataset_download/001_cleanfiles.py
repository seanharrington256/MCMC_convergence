#!/usr/bin/python

import os
replacement = """<html><pre>"""
for dname, dirs, files in os.walk("./DATASETS"):
    for fname in files:
        fpath = os.path.join(dname, fname)
        with open(fpath) as f:
            s = f.read()
        s = s.replace("<html><pre>", "")
        with open(fpath, "w") as f:
            f.write(s)
            
import os
replacement = """</pre></html>"""
for dname, dirs, files in os.walk("./DATASETS"):
    for fname in files:
        fpath = os.path.join(dname, fname)
        with open(fpath) as f:
            s = f.read()
        s = s.replace("</pre></html>", "")
        with open(fpath, "w") as f:
            f.write(s)

import os
replacement = """-"""
for dname, dirs, files in os.walk("./DATASETS"):
    for fname in files:
        fpath = os.path.join(dname, fname)
        with open(fpath) as f:
            s = f.read()
        s = s.replace("-", "_")
        with open(fpath, "w") as f:
            f.write(s)

import os
replacement = """'"""
for dname, dirs, files in os.walk("./DATASETS"):
    for fname in files:
        fpath = os.path.join(dname, fname)
        with open(fpath) as f:
            s = f.read()
        s = s.replace("'", "_")
        with open(fpath, "w") as f:
            f.write(s)

import os
replacement = """ """
for dname, dirs, files in os.walk("./DATASETS"):
    for fname in files:
        fpath = os.path.join(dname, fname)
        with open(fpath) as f:
            s = f.read()
        s = s.replace(" ", "_")
        with open(fpath, "w") as f:
            f.write(s)


import os
replacement = """|"""
for dname, dirs, files in os.walk("./DATASETS"):
    for fname in files:
        fpath = os.path.join(dname, fname)
        with open(fpath) as f:
            s = f.read()
        s = s.replace("|", "_")
        with open(fpath, "w") as f:
            f.write(s)