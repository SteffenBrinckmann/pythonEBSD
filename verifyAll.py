#!/usr/bin/env python3
## Run on all .txt files the doctest, convert to doxygen and run doxygen
import doctest, os, shutil, hashlib, sys
from distutils import dir_util
import matplotlib.pyplot as plt
from PIL import Image

if __name__ == "__main__":
  noFailure = True
  for fileName in os.listdir("."):
    if fileName.endswith(".doctest"):
      if len(sys.argv)==1 or sys.argv[1]=="skipDoctest":
        result = doctest.testfile(fileName)  #, optionflags=doctest.REPORT_ONLY_FIRST_FAILURE)
        print(("%-30s %-30s"%(fileName,result,)))
        if result.failed>0:
          noFailure = False
      txtFile = open(fileName,'r')
      doxyFile= open(fileName[:-7]+"doxy",'w')
      for line in txtFile:
        if "#doctest: +SKIP" in line:  line=line.replace('#doctest: +SKIP','')
        if ", doctest=True)" in line:  line=line.replace(', doctest=True)',')')
        if "doctestImage("   in line:  
          line=line.replace('>>> doctestImage("','')
          line=line.replace('")\n','')
          line="\\endverbatim\n"+\
               "\\image html "+line+".png width=40%\n"+\
               "\\verbatim"
        elif "doctest"       in line:  continue
        doxyFile.write(line)
      txtFile.close()
      doxyFile.close()

  if os.path.exists('doxygenOutput.txt'):
    os.unlink("doxygenOutput.txt")
  dir_util.copy_tree("HTMLInputStatic","docs/HTMLInputStatic")
  dir_util.copy_tree("HTMLInputDynamic","docs/HTMLInputDynamic")
  os.system("doxygen")
  os.system("rm *.doxy")
  if os.path.exists('doctest.png'):
    os.unlink('doctest.png')

  #commit to github
  if noFailure:
    message = input("Enter github message? [empty: no commit] ")
    if message != "":
      gitignoreMaster = ".directory\n.gitignore\ndoxygenOutput.txt\n*.pyc\n\n.vscode/\ndocs/\nHTMLInputDynamic_TEMP/\nBackupEBSD/\n"
      gitignoreGHPages= ".directory\n.gitignore\ndoxygenOutput.txt\n*.pyc\n*.py\n*.doctest\nDoxyfile\narial.ttf\n\n.vscode/\nHTMLInputDynamic/\nHTMLInputStatic/\nExamples/\nBackupEBSD/\n"

      #Master branch
      with open(".gitignore", "w") as fout:
        fout.write(gitignoreMaster)
      os.system("mv HTMLInputDynamic HTMLInputDynamic_TEMP")
      os.system("mkdir HTMLInputDynamic")
      with open("HTMLInputDynamic/empty.txt", "w") as fout:
        fout.write("")
      os.system("git rm -r --cached .")
      os.system("git add .")
      os.system('git commit -m "'+message+'"')
      os.system('git push -u origin master')

      #GH-pages branch
      os.system('git checkout -b gh-pages')
      with open(".gitignore", "w") as fout:
        fout.write(gitignoreGHPages)
      os.system("mv HTMLInputDynamic_TEMP HTMLInputDynamic")
      os.system("git rm -r --cached .")
      os.system("git add .")
      os.system('git commit -m "gh-pages"')
      os.system('git push origin gh-pages')

      #Back to master
      with open(".gitignore", "w") as fout:
        fout.write(gitignoreMaster)
      os.system("git checkout master")
      os.system("git branch -D gh-pages")


def doctestImage(fileName):
  orgFile = "HTMLInputDynamic/"+fileName+".png"
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
