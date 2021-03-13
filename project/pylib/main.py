# *-* coding: utf-8 *-*

"""
Created on mar 02 mar 2021 13:55:46 UTC

@author: vekemans

"""

from utils import *
from websites import *

import time
import socket
import multiprocessing
 
NUM_WORKERS = 4
 
start_time = time.time()
 
with multiprocessing.Pool(processes=NUM_WORKERS) as pool:
    results = pool.map_async(check_website, WEBSITE_LIST)
    results.wait()

end_time = time.time()        
 
print("Time for MultiProcessingSquirrel: %ssecs" % (end_time - start_time))
 
# WARNING:root:Timeout expired for website http://really-cool-available-domain.com
# WARNING:root:Timeout expired for website http://another-really-interesting-domain.co
# WARNING:root:Website http://bing.com returned status_code=405
# Time for MultiProcessingSquirrel: 2.8224599361419678secs

# import sys,io
# 
# print(len(sys.argv))
# for i in sys.argv:
# 	print(i)
# 
# print(__file__ == sys.argv[0])
# print(__name__)
# 
# quit()
# 
# old_stdout = sys.stdout
# sys.stdout = _buffer = io.StringIO()
# 
# print('Hellow World !')
# a = 123
# print(a)
# 
# sys.stdout = old_stdout
# 
# whatWasPrinted = _buffer.getvalue()
# print(whatWasPrinted)
# _buffer.close()
