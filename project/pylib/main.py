# *-* coding: utf-8 *-*

"""
Created on mar 02 mar 2021 13:55:46 UTC

@author: vekemans

"""

import sys,io

old_stdout = sys.stdout
sys.stdout = _buffer = io.StringIO()

print('Hellow World !')
a = 123
print(a)

sys.stdout = old_stdout

whatWasPrinted = _buffer.getvalue()
print(whatWasPrinted)
_buffer.close()
