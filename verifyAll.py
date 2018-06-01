#!/usr/bin/env python
## Run on all .txt files the doctest, convert to doxygen and run doxygen
import doctest
import os

for fileName in os.listdir("."):
  if fileName.endswith(".doctest"):
    result = doctest.testfile(fileName)
    print "%-30s %-30s"%(fileName,result,)
    txtFile = open(fileName,'r')
    doxyFile= open(fileName[:-7]+"doxy",'w')
    for line in txtFile:
      if "#doctest: +SKIP" in line:  line=line.replace('#doctest: +SKIP','')
      elif "doctest"       in line:  continue
      doxyFile.write(line)
    txtFile.close()
    doxyFile.close()

os.unlink("doxygenOutput.txt")
os.system("doxygen")
os.system("rm *.doxy")
os.unlink('doctest.png')