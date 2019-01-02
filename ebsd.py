# -*- coding: utf-8 -*-
##
# @file
# @brief Class to allow for read EBSD data
#
import math
from PIL import Image, ImageDraw, ImageFont, ImageChops
import time, os, sys, io
import numpy as np
import scipy.ndimage as ndi
from scipy.stats import gaussian_kde
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib import mlab, colors
from ebsd_Orientation import Orientation
from ebsd_Symmetry import Symmetry
from ebsd_Quaternion import Quaternion


class EBSD:
  """Class to allow for read EBSD data.

  Uses quaternions and symmetries but not orientations
  """

  ##
  # @name INPUT METHODS
  # @{
  def __init__(self, fileName, doctest=False):
    """
    read input file <br>
    initialize things<br>
    .ang or .osc file format
    - Header: ASCII information starting with #

    Args:
       fileName: file name in the present directory
    """
    # initialize
    startTime = time.time()
    fontName = 'arial.ttf'
    self.fontFile = ""
    for path in sys.path:
      fontFile = os.path.join(path, fontName)
      if os.path.isfile(fontFile):
        # print "Found font file:", fontFile
        self.fontFile = fontFile
    self.scanUnit = "um"
    self.sym = []
    self.doctest = doctest

    # read input file header and parse it
    self.fileName = fileName
    if (self.fileName[-3:] == "ang"):
      self.loadANG()
    elif (self.fileName[-3:] == "osc"):
      self.loadOSC()
    elif (self.fileName[-3:] == "txt"):
      self.loadTXT()
    elif (self.fileName[-3:] == "crc"):
      self.loadCRC()
    elif (self.fileName[:4] == 'void'):
      print("Void mode", self.fileName[4:])
      self.loadVoid(self.fileName[4:])
    else:
      print("This file-extension is not implemented yet")
      sys.exit(2)

    #basic tests
    if self.stepSizeY is None: self.stepSizeY=self.stepSizeX
    if self.stepSizeX<0.00001 and self.stepSizeY>0.00001:
      self.stepSizeX = self.stepSizeY
    print ("   Read file with step size:",self.stepSizeX, self.stepSizeY)
    print ("   Optimal image pixel size:",int(self.width/self.stepSizeX))
    print ("   Number of points:",len(self.x))
    if np.max(self.x)>10000.0 or np.max(self.y)>10000.0 :
      print ("Error in reading in ebsd.py (possibly latest version of .osc)")
      return

    # convert into quaternions and only use that
    eulers = np.vstack((self.phi1, self.PHI, self.phi2))
    self.quaternions = Quaternion.fromEulers(eulers)
    del self.phi1; del self.PHI; del self.phi2

    # for plotting: determine image and imageSize once, use multiple times
    self.image = None
    self.mask = self.CI > -1  # all are visible initially
    self.vMask = np.ones_like(self.x, dtype=np.bool)
    self.periodicLen = np.where( (self.x[1:]-self.x[:-1])<0 )[0][0]+1
    if not self.doctest:
      print("   Duration init: ", int(np.round(time.time()-startTime)), "sec")
    return

  def loadANG(self, fileName=None):
    """
    Load .ang file: filename saved in self. No need to use it
    """
    if fileName != None: self.fileName = fileName
    print("Load .ang file: ", self.fileName)
    keys = ['MaterialName', 'LatticeConstants', 'WorkingDistance', 'SEMVoltage', 'GRID:', "Symmetry"]
    fileHandle = open(self.fileName, 'r')
    keyValues = [''] * len(keys)  # actual values
    for line in fileHandle:
      if line[0:10] == "# OPERATOR":
        break
      for key in keys:
        searchTerm = "# "+key
        if searchTerm == line[0:len(searchTerm)]:
          index = keys.index(key)
          value = line.rstrip().split()[2:]
          if len(value)==1:
            value = value[0]
            try:
              value = float(value)
            except:
              pass
          keyValues[index] = value
          break
    meta = dict(list(zip(keys,keyValues)))
    if meta['Symmetry'] == 43:
      self.sym.append( Symmetry('cubic'))
    else:
      print("ERROR: no symmetry found")
      return
    # read data: print "Reading file, this can take a bit..."
    data = np.loadtxt(fileHandle)
    self.phi1      = data[:,0].astype(np.float)
    self.PHI       = data[:,1].astype(np.float)
    self.phi2      = data[:,2].astype(np.float)
    self.x         = data[:,3].astype(np.float)
    self.y         = data[:,4].astype(np.float)
    self.IQ        = data[:,5].astype(np.float)
    self.CI        = data[:,6].astype(np.float)
    self.phaseID   = data[:,7].astype(np.uint8)
    self.SEMsignal = data[:,8].astype(np.uint8)
    self.fit       = data[:,9].astype(np.float)
    self.width     = max(self.x)
    self.height    = max(self.y)
    self.ratio     = self.width/self.height
    self.stepSizeX = self.x[1] - self.x[0]
    self.stepSizeY = None
    fileHandle.close()
    del data
    return


  def loadTXT(self, fileName=None, update=False):
    """
    read txt file and possibly update data. Warning, this resets the mask to the one of the file

    Update makes more sense if you have original data and update it with some partial information,
    since the partial information is incomplete (stepSize, width and height are in many cases wrong).
    Update will keep the old data and overwrite the new. THIS IS SLOWER THAN CREATING NEW

    Args:
      fileName: fileName to load (partition data from OIM)
      update: update data or read new (read-new: default)
    """
    print("TODO: Symmetry has to be read and used")
    startTime = time.time()
    print("Load .txt file:",fileName)
    if fileName is None:
      fileName = self.fileName
    fileHandle = open(fileName,'r')
    foundKeys = {}
    for line in fileHandle:
      if line[0] != "#": break
      parts = line.split()
      if len(parts)<2: continue
      if parts[1]=="Header:":
        print("   Header: ",parts[2])
        continue
      foundKeys[parts[3]] = int( parts[2].split(":")[0].split("-")[0] )
    print("   Found data:",foundKeys)
    if "Grain" in foundKeys:  # open new array if data exists
      self.grainID= -np.ones_like(self.phaseID)

    # read data
    data = np.loadtxt(fileName)
    print("   Reading file of size ",data.shape,"  this can take a bit...")
    if update:
      self.mask[:] = False
      print("Warning: this is too slow")
      """
      be intelligent where you seearch, check if old and new data monotonically increases
      then search in sections of equal y
      or subdivide into half, of half of half
      """
      for i in range(data.shape[0]):
        x    = data[i,foundKeys["x,"]   - 1].astype(np.float32)
        y    = data[i,foundKeys["x,"]   - 0].astype(np.float32)
        # identify index: nice and much much slowes
        idx = np.argmax(np.logical_and(np.abs(self.x-x)<self.stepSizeX/10.0, \
                                       np.abs(self.y-y)<self.stepSizeX/10.0))     #very save error of 10th of STEPSIZE
        # update
        self.mask[idx] = True
        self.phi1[idx] = data[i,foundKeys["phi1,"] - 1].astype(np.float16)
        self.PHI[idx]  = data[i,foundKeys["phi1,"] - 0].astype(np.float16)
        self.phi2[idx] = data[i,foundKeys["phi1,"] + 1].astype(np.float16)
        if "IQ" in foundKeys:
          self.IQ[idx] = data[i,foundKeys["IQ"] - 1].astype(np.float16)
        if "CI" in foundKeys:
          self.CI[idx] = data[i,foundKeys["CI"] - 1].astype(np.float16)
        if "Fit" in foundKeys:
          self.fit[idx]= data[i,foundKeys["Fit"] - 1].astype(np.float16)
        if "Phase" in foundKeys:
          self.phaseID[idx]   = data[i,foundKeys["Phase"] - 1].astype(np.float16)
        if "sem" in foundKeys:
          self.SEMsignal[idx] = data[i,foundKeys["sem"] - 1].astype(np.float16)
        if "Grain" in foundKeys:
          self.grainID[idx] = data[i,foundKeys["Grain"] - 1].astype(np.float16)
        # stepSizeX, width, height etc do not change
    else: #read new
      self.mask = True
      self.phi1 = data[:,foundKeys["phi1,"] - 1].astype(np.float16)
      self.PHI  = data[:,foundKeys["phi1,"] - 0].astype(np.float16)
      self.phi2 = data[:,foundKeys["phi1,"] + 1].astype(np.float16)
      self.x    = data[:,foundKeys["x,"]   - 1].astype(np.float32)
      self.y    = data[:,foundKeys["x,"]   - 0].astype(np.float32)
      if "IQ" in foundKeys:
        self.IQ = data[:,foundKeys["IQ"] - 1].astype(np.float16)
      if "CI" in foundKeys:
        self.CI = data[:,foundKeys["CI"] - 1].astype(np.float16)
      if "Fit" in foundKeys:
        self.fit= data[:,foundKeys["Fit"] - 1].astype(np.float16)
      if "Phase" in foundKeys:
        self.phaseID   = data[:,foundKeys["Phase"] - 1].astype(np.float16)
      if "sem" in foundKeys:
        self.SEMsignal = data[:,foundKeys["sem"] - 1].astype(np.float16)
      self.mask = np.ones_like(self.x, dtype=np.bool)
      self.width     = max(self.x)
      self.height    = max(self.y)
      self.ratio     = self.width/self.height
      delta = self.x[1:] - self.x[:-1]
      self.stepSizeX = np.min(delta[delta>0]) #estimate since not ordered
      self.stepSizeY = None
    fileHandle.close()
    print("Duration loadTXT: ",int(np.round(time.time()-startTime)),"sec")
    return


  def writeANG(self, fileName):
    """
    write body of ang file

    Args:
       fileName: file name
    """
    startTime = time.time()
    fileOut = open(fileName, 'w')
    fileOut.write("# MaterialName void\n")
    fileOut.write("# Formula \n")
    fileOut.write("# Symmetry 43\n")  #adopt for HCP (fcc and bcc the same)
    fileOut.write("# LatticeConstants 1.0 1.0 1.0 90.0 90.0 90.0\n")
    fileOut.write("# NumberFamilies 4\n")
    fileOut.write("# khlFamilies 1 1 1 1 0.0\n")  #adopt for HCP
    fileOut.write("# khlFamilies 2 0 0 1 0.0\n")
    fileOut.write("# khlFamilies 2 2 0 1 0.0\n")
    fileOut.write("# khlFamilies 3 1 1 1 0.0\n")
    fileOut.write("#\n# GRID: HexGrid\n#\n")   #TODO adopt
    for i in range(len(self.x)):
      phi1, PHI, phi2 = tuple(self.quaternions[i].asEulers())
      fileOut.write(" %8.5f %8.5f %8.5f %12.5f %12.5f %8.3f %6.3f %2d %6d %7.3f\n" % \
	      (phi1,PHI,phi2,self.x[i],self.y[i],self.IQ[i],self.CI[i],self.phaseID[i],self.SEMsignal[i],self.fit[i]) )
    fileOut.close()
    print("Duration writeANG: ",int(np.round(time.time()-startTime)),"sec")
    return



  def loadOSC(self, fileName=None):
    """
    Load .osc file; filename saved in self. No need to use it.
    Copied from mtex and translated into python
    Warning: SEMsignal not parsed correctly
    """
    print("TODO: Symmetry has to be read and used")
    if fileName!=None: self.fileName = fileName
    print("Load .osc file: ",self.fileName)
    def find_subsequence(seq, subseq):
      target = np.dot(subseq, subseq)
      candidates = np.where(np.correlate(seq, subseq, mode='valid') == target)[0]
      # some of the candidates entries may be false positives, double check
      check = candidates[:, np.newaxis] + np.arange(len(subseq))
      mask = np.all((np.take(seq, check) == subseq), axis=-1)
      return candidates[mask]

    f = open(self.fileName,"r")
    header = np.fromfile(f, dtype=np.uint32, count=8)
    n = header[6]   #number of data points

    # find start position by using startByte-pattern
    bufferLength = int(math.pow(2,20))
    startBytes = np.array( [ int(i,16)  for i in ['B9', '0B', 'EF', 'FF', '02', '00', '00', '00'] ],dtype=np.uint8 )
    startPos = 0
    f.seek(startPos)
    startData = np.fromfile(f, dtype=np.uint8, count=bufferLength)
    # startPos += [x for x in xrange(len(startData)-len(startBytes)) if (startData[x:x+len(startBytes)] == startBytes).all() ][0]
    startPos += find_subsequence( startData, startBytes)[0]
    f.seek(startPos+8)

    # there are different osc file versions, one does have some count of data, the other proceeds with xStep and yStep (!=1)
    dn = np.double( np.fromfile(f,dtype=np.uint32, count=1))
    if round(((dn/4-2)/10)/n) != 1:
      f.seek(startPos+8)
    self.stepSizeX = np.double(np.fromfile(f, dtype=np.float32, count=1))
    self.stepSizeY = np.double(np.fromfile(f, dtype=np.float32, count=1))

    data = np.reshape( np.double(np.fromfile(f, count=n*10, dtype=np.float32)) , (n,10) )
    self.phi1    = data[:,0].astype(np.float16)
    self.PHI     = data[:,1].astype(np.float16)
    self.phi2    = data[:,2].astype(np.float16)
    self.x       = data[:,3].astype(np.float32)
    self.y       = data[:,4].astype(np.float32)
    self.IQ      = data[:,5].astype(np.float16)
    self.CI      = data[:,6].astype(np.float16)
    self.phaseID = data[:,7].astype(np.float16)
    self.SEMsignal= data[:,8].astype(np.float16) #SEMSignal
    self.fit     = data[:,9].astype(np.float16) #Fit
    self.width   = max(self.x)
    self.height  = max(self.y)
    self.ratio   = self.width/self.height
    f.close()
    del data
    return



  def loadCRC(self, fileName=None):
    """
    Load .crc file; filename saved in self. No need to use it.
    Copied from mtex and translated into python
    """
    import struct
    if fileName!=None: self.fileName = fileName
    cprFileName = self.fileName[:-4]+".cpr"
    print("Load .crc file: ",self.fileName,cprFileName)
    if not os.path.exists(cprFileName):
      print("CPR file does not exist")
    cprFile = open(cprFileName,"r")
    cprData = {}
    for line in cprFile:
      line = line.strip()
      if line[0] == '[':
        title = line[1:-1].lower()
        cprData[title] = {}
        continue
      key, value = line.split('=')[0],line.split('=')[1]
      try:
        cprData[title][key.lower()] = float(value)
      except:
        cprData[title][key.lower()] = value.lower()
    cprFile.close()
    # print "META DATA",cprData
    self.stepSizeX = np.double(cprData['job']['griddistx'])
    self.stepSizeY = np.double(cprData['job']['griddisty'])
    xcells         = int(cprData['job']['xcells'])
    ycells         = int(cprData['job']['ycells'])
    numDataPoints  = xcells * ycells
    self.width     = xcells * self.stepSizeX
    self.height    = ycells * self.stepSizeY
    self.ratio     = self.width/self.height
    if cprData['phase1']['lauegroup'] == 11:
      self.sym.append( Symmetry()        )  #phase 0: default = not identified
      self.sym.append( Symmetry('cubic') )  #phase 1: cubic
    else:
      print("ERROR: no symmetry found")
      return

    # verify that data in correct order
    allColumnNames = [
      'X',                  # 1    4 bytes
      'Y',                  # 2       "
      'phi1',               # 3       "
      'Phi',                # 4       "
      'phi2',               # 5       "
      'MAD',                # 6       "
      'BC',                 # 7    1 byte
      'BS',                 # 8       "
      'Unknown',            # 9       "
      'Bands',              # 10      "
      'Error',              # 11      "
      'ReliabilityIndex']    # 12      "
    allDataType = np.ones((12,),dtype=np.int)
    allDataType[:6]=4; allDataType[-1]=4
    columnNames, columnType = ['Phase'], [1]
    for k in range(int(cprData['fields']['count'])):
      order = int(cprData['fields']['field'+str(k+1)])-1
      if order <= 12:
        columnNames.append( allColumnNames[order] )
        columnType.append(  allDataType[order]    )
      else:
        columnNames.append( 'Unknown'+str(order) )
        columnType.append( 4 )
    if columnNames == ['Phase', 'phi1', 'Phi', 'phi2', 'MAD', 'BC', 'BS', 'Bands', 'Error', 'ReliabilityIndex']:
      print("  CRC-Data in correct order")
    else:
      print("  WARNING! CRC-Data not in correct order! WARNING")
      print("    should be ['Phase', 'phi1', 'Phi', 'phi2', 'MAD', 'BC', 'BS', 'Bands', 'Error', 'ReliabilityIndex']")
      print("    is       ",columnNames)
      print("    if data missing at end, no problem")
    # print columnType

    # coordinates
    x_ = np.arange(xcells)*self.stepSizeX
    y_ = np.arange(ycells)*self.stepSizeY
    self.x, self.y = np.meshgrid(x_,y_)
    self.x, self.y = self.x.flatten(), self.y.flatten()

    # read data from crcFile
    crcFile = open(self.fileName,'rb')
    self.phaseID = np.zeros((numDataPoints),dtype=np.uint8)
    self.BC,self.BS,self.Bands,self.Error=np.zeros_like(self.phaseID),np.zeros_like(self.phaseID),np.zeros_like(self.phaseID),np.zeros_like(self.phaseID)
    self.phi1    = np.zeros((numDataPoints),dtype=np.float)
    self.PHI,self.phi2,self.CI,self.RI  =np.zeros_like(self.phi1),np.zeros_like(self.phi1),np.zeros_like(self.phi1),np.zeros_like(self.phi1)
    self.IQ, self.SEMsignal, self.fit   =np.zeros_like(self.phi1),np.zeros_like(self.phi1),np.zeros_like(self.phi1)
    for i in range(numDataPoints):
      self.phaseID[i] = struct.unpack('B', crcFile.read(1))[0]
      self.phi1[i]    = struct.unpack('f', crcFile.read(4))[0]
      self.PHI[i]     = struct.unpack('f', crcFile.read(4))[0]
      self.phi2[i]    = struct.unpack('f', crcFile.read(4))[0]
      self.CI[i]     = struct.unpack('f', crcFile.read(4))[0]
      self.BC[i]      = struct.unpack('B', crcFile.read(1))[0]
      self.BS[i]      = struct.unpack('B', crcFile.read(1))[0]
      self.Bands[i]   = struct.unpack('B', crcFile.read(1))[0]
      self.Error[i]   = struct.unpack('B', crcFile.read(1))[0]
      if 'ReliabilityIndex' in columnNames:
        self.RI[i]      = struct.unpack('f', crcFile.read(4))[0]
    crcFile.close()
    if not len(self.sym) == np.max(self.phaseID)-np.min(self.phaseID)+1:
      print("ERRRO in reading CRC: symmetries do not match",len(self.sym), np.max(self.phaseID)-np.min(self.phaseID)+1)
    return


  def loadVoid(self, rotation):
    """
    rotation angles in degree
    """
    numPerAxis, distrib = 6, 0
    if "|" in rotation:
      rotation = [float(i) for i in rotation.split('|')]
      if len(rotation)==3:
        phi1, PHI, phi2 = np.radians(rotation)
      elif len(rotation)==4:
        phi1, PHI, phi2 = np.radians(rotation[:3])
        distrib         = rotation[-1]
      elif len(rotation)==5:
        phi1, PHI, phi2 = np.radians(rotation[:3])
        distrib, numPerAxis  = rotation[-2:]
      else:
        print("ERROR")
        return
      print("   Euler angles:",np.round(phi1,2), np.round(PHI,2), np.round(phi2,2),\
         "| distribution:",distrib,"| numberPerAxis:",numPerAxis)
    else:
      phi1, PHI, phi2 = 0,0,0
    if distrib < 0.001: distrib=0.001
    self.sym.append( Symmetry('cubic') )
    self.stepSizeX = 1.
    numDataPoints = int(numPerAxis**2)
    x_ = np.arange(numPerAxis)*self.stepSizeX
    self.x, self.y = np.meshgrid(x_,x_)
    self.x, self.y = self.x.flatten(), self.y.flatten()
    self.phaseID = np.ones((numDataPoints),dtype=np.uint8)
    self.phi1    = np.zeros((numDataPoints),dtype=np.float)+phi1 #+ np.random.normal(loc=0,scale=distrib,size=numDataPoints)
    self.PHI     = np.zeros((numDataPoints),dtype=np.float)+PHI #+ np.random.normal(loc=0,scale=distrib,size=numDataPoints)
    self.phi2    = np.zeros((numDataPoints),dtype=np.float)+phi2 #+ np.random.normal(loc=0,scale=distrib,size=numDataPoints)
    self.CI      = np.ones((numDataPoints),dtype=np.float)
    self.stepSizeY = self.stepSizeX
    self.width     = np.max(self.x)
    self.height    = np.max(self.y)
    self.ratio   = self.width/self.height
    return


  # @}
  ##
  # @name Mask and Path routines
  # masked areas are plotted in black. Hence initially no point is part of the mask, i.e. all points are false
  # @{

  def maskCI(self, CI):
    """
    masked all points off, which have a CI less than: good points=False, bad points=True

    Args:
       CI: critical CI
    """
    self.mask = self.CI > CI
    return


  def maskReset(self):
    """
    reset mask
    """
    self.mask = self.CI > -1
    return

  def removePointsOutsideMask(self):
    """
    set all data-points outside of mask to invalid such that after export to OIM, it will be read there as non-existing points
    """
    self.CI[~self.mask]   = -1.0
    self.fit[~self.mask]  = 180.0
    self.phi1[~self.mask] = 12.5625
    self.PHI[~self.mask]  = 12.5625
    self.phi2[~self.mask] = 12.5625
    return



  def setVMask(self,every=1):
    """
    mask every kth point off, for fast plotting, does not influence results in any way

    Args:
       every: use only every k-th point. Improves plotting speed. every=1 resets.
    """
    self.vMask[:] = False
    self.vMask[::every] = True
    return


  def cropVMask(self, xmin=None, ymin=None, xmax=None, ymax=None):
    """
    crop visible area

    Args:
       xmin: minimum x-coordinate
       ymin: minimum y-coordinate
       xmax: maximum x-coordinate
       ymax: maximum y-coordinate
    """
    if not xmin: xmin=0
    if not xmax: xmax=np.max(self.x)
    if not ymin: ymin=0
    if not ymax: ymax=np.max(self.y)
    # print xmin,xmax,ymin, ymax
    self.vMask = np.logical_and(self.vMask,  self.x>=xmin)
    self.vMask = np.logical_and(self.vMask,  self.x<=xmax)
    self.vMask = np.logical_and(self.vMask,  self.y<=ymax)
    self.vMask = np.logical_and(self.vMask,  self.y>=ymin)
    return


  def neighbors(self, idx=None,layers=1):
    """
    identify neighboring indexes

    Args:
       idx: index to find [if None: calculate all]
       layers: number of neighboring layers

    Returns:
     array of neighbors; invalid points have a value=-10
    """
    l = self.periodicLen
    if layers==1:
      if idx is None:
        original = np.outer(np.arange(len(self.x),dtype=np.int), np.ones([6,],dtype=np.int))
        neighbors = original.copy()
        neighbors[:,0] += -l
        neighbors[:,1] += -l+1
        neighbors[:,2] += -1
        neighbors[:,3] += +1
        neighbors[:,4] += +l-1
        neighbors[:,5] += +l
      else:
        neighbors = np.array([-l,-l+1,  -1,+1,  +l-1,+l])+idx
    elif layers==2:
      neighbors = np.array([-l-1,-l,-l+1,-l+2,  -2,-1,+1,+2,  +l-2,+l-1,+l,+l+1])+idx
    else:
      print("number of layers not implemented")
      return None
    neighbors[neighbors<0]            = -10
    neighbors[neighbors>=len(self.x)] = -10
    mask = (self.x[neighbors]-self.x[original])>(self.stepSizeX*1.1)  #only check for distance in x; b/c check in y already done by previous lines
    neighbors[mask] = -10
    return neighbors


  def calcKAM(self, layers=1):
    """
    calculate Kerner Average Misorientation in DEGREES (because user focused)

    Args:
       layers: number of neighboring layers used for KAM (more: slower)
    """
    startTime    = time.time()
    sym          = self.sym[0]
    neighbors    = self.neighbors()
    fzThreshold  = math.sqrt(2.0)-1.0
    qConj        = self.quaternions.conjugated()
    angles       = np.empty_like(neighbors, dtype=np.float)
    neighborSymQ = sym.symmetryQuats()
    symQ         = neighborSymQ[0]
    for iNeighbor in range(6):
      neighborQ    = self.quaternions[neighbors[:,iNeighbor]]
      misQ         = (self.quaternions.conjugated() * neighborQ).copy()
      foundAngle   = np.zeros( (len(self.x)), dtype=np.bool )
      for nSQ in neighborSymQ:
        theQ = symQ.conjugated()*misQ*nSQ
        for k in xrange(2): #try both conjugated versions
          theQ.conjugate()  #verified before
          theQ_Rod = abs(theQ.asRodrigues())
          inFZ = np.logical_and( \
            np.logical_and(fzThreshold>=theQ_Rod[0],fzThreshold>=theQ_Rod[1]) , \
            np.logical_and(fzThreshold>=theQ_Rod[2],1.0>=np.sum(theQ_Rod, axis=0)) )
          #angle = theQ.asAngleAxis()[0]  #much slower: requires additional class; slight differences to faster version
          angle= 2.0*np.arctan(  np.linalg.norm(theQ_Rod, axis=0)  )
          foundAngle[inFZ] = True
          angles[inFZ,iNeighbor] = angle[inFZ]
          mask = self.CI[ neighbors[:,iNeighbor] ]==-1.0
          angles[mask,iNeighbor] = np.nan
        if np.all(foundAngle): break   #stop looking for alternatives if filled already all
    self.kam = np.degrees( np.nanmean(angles, axis=1) )
    self.kam[ self.CI==-1.0 ] = np.nan
    print("Duration KAM evaluation: ",int(np.round(time.time()-startTime)),"sec")
    return


  # @}
  ##
  # @name PLOT METHODS
  #@{
  def plot(self, vector, widthPixel=None, vmax="", vmin="", interpolationType="nearest", cmap=None, show=True, cbar=True):
    """
    given a class-vector, plot the vector as an image<br>
    the x and y are given by the class-vector x and y

    Args:
       vector: vector to be plotted as a 2D image
       widthPixel: rescale to horizontal size of the image [default: optimal pixel width]
       vmax: rescale z-scale to maximal value
       vmin: rescale z-scale to minimal value
       interpolationType: interpolation type [default: "nearest" next-neighbor]
    """
    startTime = time.time()
    if widthPixel is None:
      widthPixel = int(self.width/self.stepSizeX)

    # create a special cmap palette with blacK as value for bad-numbers
    if cmap is None:
      cmap = cm.Spectral
      cmap.set_bad('k', 1.0)
    # create a new grid with the given resolution, and interpolate
    xMax = np.max(self.x[self.vMask])
    xMin = np.min(self.x[self.vMask])
    yMax = np.max(self.y[self.vMask])
    yMin = np.min(self.y[self.vMask])
    self.ratio       = (xMax-xMin)/ (yMax-yMin)
    heightPixel = int(widthPixel/self.ratio)
    xAxis = np.linspace(xMin, xMax, widthPixel)
    yAxis = np.linspace(yMin, yMax, heightPixel)
    x, y = np.meshgrid(xAxis, yAxis)
    points = np.vstack(  (self.x[self.vMask], self.y[self.vMask]) ).T
    z    = griddata( points, vector[self.vMask],     (x,y), interpolationType)
    mask = griddata( points, ~self.mask[self.vMask], (x,y), interpolationType)
    # plot if/if-not the maximum and minimum are given
    if vmax!="" and vmin!="":
      plt.imshow( np.ma.masked_where(mask, z.astype(np.float32)) , extent=[xMin,xMax, yMax,yMin], cmap=cmap, vmax=vmax, vmin=vmin, origin='upper')
    else:
      plt.imshow( np.ma.masked_where(mask, z.astype(np.float32)) , extent=[xMin,xMax, yMax,yMin], cmap=cmap, origin='upper')
    if cbar:
      plt.colorbar()
    if self.doctest: plt.savefig('doctest.png'); plt.close(); return
    print("   Plot with x and y axis in [um]")
    print("Duration plot: ",int(np.round(time.time()-startTime)),"sec")
    if show:
      plt.show()
    z *= 255/np.max(z)
    self.image  = Image.fromarray( z.astype(np.float32) )
    return


  def plotRGB(self, rgb, widthPixel=256, interpolationType="nearest", fileName=None):
    """
    given a RGB vector (same size as the other class vectors)
    plot the vector as an image<br>
    the x and y are given by the class-vector x and y
    USED INTERNALLY

    Args:
       rgb: matrix [3, classVectorSize] to be plotted as a 2D image
       widthPixel: horizontal size of the image [default: 256 pixel]
       interpolationType: interpolation type [default: "nearest"]
       fileName: save to file instead of showing
    """
    # create a new grid with the given resolution
    xMax = np.max(self.x[self.vMask])
    xMin = np.min(self.x[self.vMask])
    yMax = np.max(self.y[self.vMask])
    yMin = np.min(self.y[self.vMask])
    self.ratio       = (xMax-xMin)/ (yMax-yMin)
    heightPixel = int(widthPixel/self.ratio)
    xAxis = np.linspace(xMin, xMax,  widthPixel)
    yAxis = np.linspace(yMin, yMax,  heightPixel)
    x, y = np.meshgrid(xAxis, yAxis)
    # filter out using the mask: assign 0 to the mask on the rgb values
    #  ensure that the left hand right side of = have the same mask
    rgb[0,:][ ~self.mask ] = np.zeros( (len(self.x)) )[ ~self.mask ]
    rgb[1,:][ ~self.mask ] = np.zeros( (len(self.x)) )[ ~self.mask ]
    rgb[2,:][ ~self.mask ] = np.zeros( (len(self.x)) )[ ~self.mask ]
    # interpolate the rbg onto the red,blue,green
    points = np.vstack(  (self.x[self.vMask], self.y[self.vMask]) ).T
    red   = np.uint8(griddata( points, rgb[0,self.vMask], (x,y), interpolationType)*255)
    green = np.uint8(griddata( points, rgb[1,self.vMask], (x,y), interpolationType)*255)
    blue  = np.uint8(griddata( points, rgb[2,self.vMask], (x,y), interpolationType)*255)
    # put them all in one array and then reshape it and transpose by changing the order to 0->2->1 (determined by try and error)
    allColors = np.concatenate( (red, green, blue), axis=1)
    imageArray = np.transpose( allColors.reshape(heightPixel, 3, widthPixel), (0,2,1) )
    # finally plot
    self.image  = Image.fromarray( imageArray )
    plt.imshow( self.image, extent=[xMin,xMax, yMax, yMin], origin='upper')
    return


  def plotIPF(self, direction="ND", widthPixel=None, fileName=None, interpolationType="nearest"):
    """
    plot Inverse Pole Figure (IPF)

    Args:
       direction: default.."ND", "RD", "TD"
       widthPixel: horizontal size of the image [default: optimal size based on data]
       interpolationType: interpolation type [default: "nearest"]
       fileName: save to file instead of showing
    """
    startTime = time.time()
    if   direction=="RD": axis = [1,0,0]
    elif direction=="TD": axis = [0,1,0]
    elif direction=="ND": axis = [0,0,1]
    else:              #if first argument specifies widthPixel
      widthPixel = direction
      axis = [0,0,1]
    if widthPixel is None:
      widthPixel = int(self.width/self.stepSizeX)

    flags = np.zeros( (len(self.x)), dtype=np.bool)
    rgbs  = np.zeros( (3,len(self.x)), dtype=np.float)
    for sym in self.sym:
      if sym.__repr__() == "None": continue
      equivQuaternions = sym.equivalentQuaternions( self.quaternions )
      for equivQuaternion in equivQuaternions:
        pole          = equivQuaternion.conjugated()*axis
        flags_, rgbs_ = sym.inSST( pole[:,~flags], color=True, proper=False)
        if len(rgbs_.shape)==2:
          rgbs[:,~flags] = rgbs_
          flags[~flags] = flags_
    self.plotRGB( rgbs, widthPixel, interpolationType, fileName)
    if self.doctest: plt.savefig('doctest.png'); plt.close(); return
    print("Duration plotIPF: ",int(np.round(time.time()-startTime)),"sec")
    if fileName == None:
      plt.show()
    else:
      plt.savefig(fileName, dpi=150, bbox_inches='tight')
      plt.close()
    return


  def addSymbol(self, x, y, fileName=None, scale=1., colorCube='black'):
    """
    TODO: use version in ebsd_Orientation
    Add symbol of crystal orientation (symmetry and rotation) to IPF at given location

    Args:
       x: x-coordinate
       y: y-coordinate
       fileName: export to file
       scale: scale of symbol
       colorCube: color of symbol
    """
    def plotLine(ax, start,delta,color='k',lw=1):
      ax.plot( [start[0]]+[start[0]+delta[0]],
               [start[1]]+[start[1]+delta[1]],
               color=color,lw=lw)
      return
    def trim(im):
      bg = Image.new(im.mode, im.size, im.getpixel((0,0)))
      diff = ImageChops.difference(im, bg)
      diff = ImageChops.add(diff, diff, 2.0, -100)
      bbox = diff.getbbox()
      if bbox:
        return im.crop(bbox)

    fig = plt.figure()
    ax  = fig.add_subplot (111)
    xMax = np.max(self.x[self.vMask])
    xMin = np.min(self.x[self.vMask])
    yMax = np.max(self.y[self.vMask])
    yMin = np.min(self.y[self.vMask])
    ax.imshow( self.image, extent=[xMin,xMax, yMax, yMin], origin='upper')

    iClose      = np.argmin((self.x-x)**2 + (self.y-y)**2)
    iQuaternion = self.quaternions[iClose]
    if not self.doctest: print("Euler angles at point:",iQuaternion.asEulers(degrees=True, round=1))
    loc         = np.array([x,y,0])
    for sym in self.sym:
      if sym.__repr__() == None: continue
      for line in sym.unitCell():
        start = iQuaternion*(np.array(line[:3],dtype=np.float)*scale)
        end   = iQuaternion*(np.array(line[3:],dtype=np.float)*scale)
        # use OIM coordinate system: up-left: new vector (-y, x, z)
        # use imshow with upper origin: second coordinate negative -> (-y, -x, z)
        start = np.array([-start[1], -start[0],  start[2]])
        end   = np.array([-end[1],   -end[0],    end[2]])
        # once the orientation of crystal is correct: add location
        if start[2]<0 and end[2]<0:
          plotLine(ax, start+loc, end-start,color=colorCube,lw=0.2)
        elif start[2]>0 and end[2]>0:
          plotLine(ax, start+loc, end-start,color=colorCube,lw=2)
        else:
          delta = end-start
          k     = -start[2]/delta[2]
          mid   = start+k*delta
          if start[2]>0:
            plotLine(ax, start+loc, mid-start,color=colorCube,lw=2)
            plotLine(ax, mid+loc,   end-mid,color=colorCube,lw=0.2)
          else:
            plotLine(ax, start+loc, mid-start,color=colorCube,lw=0.2)
            plotLine(ax, mid+loc,   end-mid,color=colorCube,lw=2)
    ax.set_xticks([]); ax.set_yticks([])
    ax.axis('off')
    fig.subplots_adjust(left=0.0, right=1.0, top=1.0, bottom=0.0)
    # fig.tight_layout()
    # plt.show()
    buf = io.BytesIO()
    plt.savefig(buf, format='png')
    buf.seek(0)
    self.image = trim(Image.open(buf)).convert("RGB")
    buf.close()
    plt.close()
    plt.imshow( self.image, extent=[xMin,xMax, yMax, yMin], origin='upper')
    if self.doctest: plt.savefig('doctest.png'); plt.close(); return
    if fileName == None:
      plt.show()
    else:
      plt.savefig(fileName, dpi=150, bbox_inches='tight')
      plt.close()
    return



  def addScaleBar(self, fileName=None, site="BL", barLength=None, scale = -1, alpha=0.5):
    """
    Add scale-bar to image

    Args:
       fileName: if given, save to file
       site: where to put the scale bar<br> bottom-left "BL" (default)<br>
             bottom-right "BR"<br> top-left "TL"<br> top-right "TR"
       barLength: length of scale bar. It is calculate if not given
       scale: of font and rectangle. Default: widthInPixel / 16, which is for a 1024x786 image = 64
       alpha: transparency of scale bar background
    """
    widthPixel, heightPixel = self.image.size
    if barLength is None:
      digits = int(math.log10(round(self.width/4.)))
      barLength = round(   max(self.width,self.height)   /6., -digits)
    barPixel = int(widthPixel * barLength/self.width)
    image = self.image.copy()
    draw = ImageDraw.Draw(image, 'RGBA')
    if scale < 0:
      if widthPixel>heightPixel:  scale = widthPixel  / 32
      else:                       scale = heightPixel / 16
    font = ImageFont.truetype(self.fontFile,int(scale/5*3)  )
    # identify top-left corner of scale bar section
    if   site=="BL":  offsetX = 0;                                   offsetY = heightPixel-scale
    elif site=="BR":  offsetX = widthPixel-barPixel-scale/5; offsetY = heightPixel-scale
    elif site=="TL":  offsetX = 0;                                   offsetY = 0
    elif site=="TR":  offsetX = widthPixel-barPixel-scale/5; offsetY = 0
    else:             offsetX = 0;                                   offsetY = heightPixel-scale
    textString = str(barLength)+" "+'\u03BC'+"m"
    textWidth, textHeight = draw.textsize( textString, font=font)
    draw.rectangle((offsetX,        offsetY,         offsetX+barPixel+scale/5,  offsetY+scale    ), (255, 255, 255, int(alpha*255)))  #white background
    draw.rectangle((offsetX+scale/10, offsetY+scale*7/10, offsetX+barPixel+scale/10, offsetY+scale*9/10), 'black')    #black bar
    draw.text( (offsetX+(barPixel+scale/5-textWidth)/2,offsetY), textString, 'black', font=font)
    xMax, xMin = np.max(self.x[self.vMask]), np.min(self.x[self.vMask])
    yMax, yMin = np.max(self.y[self.vMask]), np.min(self.y[self.vMask])
    plt.imshow( image, origin='upper')
    plt.xticks([])   ; plt.yticks([])
    plt.axis('off')
    if self.doctest: fileName='doctest.png'
    if fileName == None:
      plt.show()
    else:
      plt.savefig(fileName, dpi=150, bbox_inches='tight')
      plt.close()
    return


  def plotPF(self, axis=[1,0,0], points=False, fileName=None, color='#1f77b4', alpha=1.0, show=True, density=256, size=2, proj2D='up-left', vmin=0.0, vmax=1.0):
    """
    plot pole figure

    Projection onto 2D: cooradinate systems are given as xDirection-yDirection (z follows)
    - down-right: [default in text books, mTex] RD = x = down; TD = y = right; ND = z = outOfPlane
    - up-left: [default in OIM and here] RD = x = up; TD = y = left; ND = z = outOfPlane

    Args:
      axis:    axis to plot: default: axis=1,0,0
      points:  plot individual points [default], or plot distribution
      fileName: if given, save to file
      color:   plot color
      alpha:   alpha transparency
      show:    show figure [default], False for subsequent plotting
      density: how many points to plot on the distribution
      size:    points: point size; distribution: amount of smoothing: higher more smoothing
      proj2D:  orientation of 2D projection: [down-right, up-left, None]
      vmin:    minimum value plotted, used as cut-off for transparency
      vmax:    max. used in color coding, allows to focus on minor texture
    """
    startTime = time.time()
    maxColor = tuple(np.array(colors.hex2color(color))*0.5)
    for sym in self.sym:
      if sym.__repr__() == None: continue
      oHelp = Orientation(Eulers=np.array([0.,0.,0.]), symmetry=sym.__repr__())
      axis = np.array(axis, dtype=np.float)
      axis /= np.linalg.norm(axis)
      mask = np.logical_and(self.mask, self.vMask)
      x, y = None, None
      for q in oHelp.symmetry.equivalentQuaternions(oHelp.quaternion):
        conjAxis  = q*axis
        direction = self.quaternions*conjAxis
        direction = direction[:,mask]               #filter mask
        direction = direction[:, direction[2,:]>0]  #filter upward dome
        direction[0,:] /= direction[2,:]+1.
        direction[1,:] /= direction[2,:]+1.
        if x is None:
          x,y = direction[0,:], direction[1,:]
        else:
          x,y = np.hstack((x,direction[0,:])), np.hstack((y,direction[1,:]))
    if points:
      if proj2D=='down-right':
        plt.plot(-x, y,'.', color=maxColor, markersize=size)  #markersize=0.05
      elif proj2D=='up-left':
        plt.plot(-y, x,'.', color=maxColor, markersize=size)  #markersize=0.05
      else:
        return
      plt.plot( np.cos(np.linspace(0,2*np.pi,100)), np.sin(np.linspace(0,2*np.pi,100)), 'k-')
      plt.plot( [-1,1],[0,0],'k--')
      plt.plot( [0,0],[-1,1],'k--')
    else:
      cmap = colors.LinearSegmentedColormap.from_list('my', [(1,1,1),maxColor])
      center = (density - 1)/2
      imgDim = density+2*size
      img = np.zeros((imgDim,imgDim))
      x,y = np.nan_to_num(x), np.nan_to_num(y)
      if proj2D=='down-right':   zippedList = list(zip(-x,y))
      elif proj2D=='up-left':    zippedList = list(zip(-y,x))
      else:                      return
      for x_, y_ in zippedList:
        ix = int((x_ - -1.) * center) + size
        iy = int((y_ - -1.) * center) + size
        if 0 <= ix < imgDim and 0 <= iy < imgDim:
            img[iy][ix] += 1
      img = ndi.gaussian_filter(img, (size,size))  # gaussian convolution
      img /= np.max(img)                               # normalize
      img[img<vmin] = np.nan                           #filter out low values to make transparent
      plt.imshow(img, cmap=cmap, alpha=alpha, vmin=0.0, vmax=vmax,origin='lower')
      plt.plot( center*np.cos(np.linspace(0,2*np.pi,100))+center+size,
                center*np.sin(np.linspace(0,2*np.pi,100))+center+size, 'k-', lw=2)
      plt.plot( [center+size,center+size], [size,imgDim-size], 'k--', lw=1)
      plt.plot( [size,imgDim-size], [center+size,center+size], 'k--', lw=1)
      # plt.colorbar()
    plt.xlim([-1,1]) ; plt.ylim([-1,1])
    plt.xticks([])   ; plt.yticks([])
    plt.axis('equal'); plt.axis('off')
    if self.doctest: plt.savefig('doctest.png'); plt.close(); return
    print("Duration plotPF: ",int(np.round(time.time()-startTime)),"sec")
    if fileName == None:
      plt.show()
    else:
      plt.savefig(fileName, dpi=150, bbox_inches='tight')
      plt.clf();  plt.cla()
    return
  # @}

