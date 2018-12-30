#!/usr/bin/env python3
## Run on all .txt files the doctest, convert to doxygen and run doxygen
import doctest, os, shutil, hashlib
import matplotlib.pyplot as plt
from PIL import Image

if __name__ == "__main__":
  noFailure = True
  for fileName in os.listdir("."):
    if fileName.endswith(".doctest"):
      result = doctest.testfile(fileName)  #, optionflags=doctest.REPORT_ONLY_FIRST_FAILURE)
      print(("%-30s %-30s"%(fileName,result,)))
      if result.failed>0:
        noFailure = False
      txtFile = open(fileName,'r')
      doxyFile= open(fileName[:-7]+"doxy",'w')
      for line in txtFile:
        if "#doctest: +SKIP" in line:  line=line.replace('#doctest: +SKIP','')
        if ", doctest=True)" in line:  line=line.replace(', doctest=True)',')')
        elif "doctest"       in line:  continue
        doxyFile.write(line)
      txtFile.close()
      doxyFile.close()

  if os.path.exists('doxygenOutput.txt'):
    os.unlink("doxygenOutput.txt")
  os.system("doxygen")
  os.system("rm *.doxy")
  if os.path.exists('doctest.png'):
    os.unlink('doctest.png')

  #commit to github
  if noFailure:
    message = input("Enter github message? [empty: no commit] ")
    if message != "":
      gitignoreLocal = ".directory\n.gitignore\ndoxygenOutput.txt\n*.pyc\ndocs/\n"
      gitignoreRemote= ".directory\n.gitignore\ndoxygenOutput.txt\n*.pyc\n"
      with open(".gitignore", "w") as fout:
        fout.write(gitignoreRemote)
      os.system("git add -A")
      os.system('git commit -m "'+message+'"')
      os.system('git push -u origin master')
      with open(".gitignore", "w") as fout:
        fout.write(gitignoreLocal)


def doctestImage(fileName):
  orgFile = "doctestImages/"+fileName+".png"
  if os.path.exists(orgFile):
    md5New = hashlib.md5(open("doctest.png","rb").read()).hexdigest()
    md5Org = hashlib.md5(open(orgFile,"rb").read()).hexdigest()
    if md5New==md5Org:  #THE SAME
      print("doctest 1")
    else:  #not the same
      plt.subplots_adjust(left=0.0,right=1.0,wspace=0.0)
      plt.subplot(121)
      plt.imshow(Image.open("doctest.png"))
      plt.title("the same?")
      plt.subplot(122)
      plt.imshow(Image.open(orgFile))
      plt.show()
      answer = input("")
      if answer=="n":
        print("doctest: image not the same")
      else:
        shutil.copy("doctest.png", orgFile)
        print("doctest 1")
      
  else:  #Failue, old does not exist
    plt.imshow(Image.open("doctest.png"))
    plt.title("correct?")
    plt.show()
    answer = input("")
    if answer=="n":
      print("doctest: image not correct")
    else:
      shutil.copy("doctest.png", orgFile)
      print("doctest 1")
  return