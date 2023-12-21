# -*- coding: utf-8 -*-
"""
Created on Thu Dec 21 16:28:54 2023

@author: treyd
"""

import os

os.system("./raxmlHPC-PTHREADS-SSE3 -p 1234 -f a -x 1234 -s builtSeq.fasta -n pythreads.raxml -m GTRGAMMA -# 1000")
os.makedir("RAxML_results/")
os.system("mv *.raxml RAxML_results/")