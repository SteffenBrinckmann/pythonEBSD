# -*- coding: UTF-8 no BOM -*-
##
# @file
# @brief Material symmetry class; this class only stores lattice: cubic, hexagonal
#
# Part of DAMASK (http:\\damask.mpie.de); Martin Diehl, Philip Eisenlohr, Franz Roters
# Other part: Steffen Brinckmann
#
import numpy as np
import matplotlib.pyplot as plt
import math
from ebsd_Quaternion import Quaternion

##
#  Material symmetry class; this class only stores lattice: cubic, hexagonal
class Symmetry:

  lattices = [None,'orthorhombic','tetragonal','hexagonal','cubic',]

  # @name CONVENTIONAL ROUTINES
  #@{

  def __init__(self, symmetry = None):
    if isinstance(symmetry, str) and symmetry.lower() in Symmetry.lattices:
      self.lattice = symmetry.lower()
    else:
      self.lattice = None


  def __copy__(self):
    return self.__class__(self.lattice)
  copy = __copy__


  def __repr__(self):
    return '%s' % (self.lattice)


  def __eq__(self, other):
    return self.lattice == other.lattice


  def __neq__(self, other):
    return not self.__eq__(other)


  def __cmp__(self,other):
    return cmp(Symmetry.lattices.index(self.lattice),Symmetry.lattices.index(other.lattice))


  #@}
  ##
  # @name MATERIAL SPECIFIC ROUTINES
  #@{

  def symmetryQuats(self, who = []):
    '''
    List of symmetry operations as quaternions.
    '''
    if self.lattice == 'cubic':
      symQuats =  [
                    [ 1.0,              0.0,              0.0,              0.0              ],
                    [ 0.0,              1.0,              0.0,              0.0              ],
                    [ 0.0,              0.0,              1.0,              0.0              ],
                    [ 0.0,              0.0,              0.0,              1.0              ],
                    [ 0.0,              0.0,              0.5*math.sqrt(2), 0.5*math.sqrt(2) ],
                    [ 0.0,              0.0,              0.5*math.sqrt(2),-0.5*math.sqrt(2) ],
                    [ 0.0,              0.5*math.sqrt(2), 0.0,              0.5*math.sqrt(2) ],
                    [ 0.0,              0.5*math.sqrt(2), 0.0,             -0.5*math.sqrt(2) ],
                    [ 0.0,              0.5*math.sqrt(2),-0.5*math.sqrt(2), 0.0              ],
                    [ 0.0,             -0.5*math.sqrt(2),-0.5*math.sqrt(2), 0.0              ],
                    [ 0.5,              0.5,              0.5,              0.5              ],
                    [-0.5,              0.5,              0.5,              0.5              ],
                    [-0.5,              0.5,              0.5,             -0.5              ],
                    [-0.5,              0.5,             -0.5,              0.5              ],
                    [-0.5,             -0.5,              0.5,              0.5              ],
                    [-0.5,             -0.5,              0.5,             -0.5              ],
                    [-0.5,             -0.5,             -0.5,              0.5              ],
                    [-0.5,              0.5,             -0.5,             -0.5              ],
                    [-0.5*math.sqrt(2), 0.0,              0.0,              0.5*math.sqrt(2) ],
                    [ 0.5*math.sqrt(2), 0.0,              0.0,              0.5*math.sqrt(2) ],
                    [-0.5*math.sqrt(2), 0.0,              0.5*math.sqrt(2), 0.0              ],
                    [-0.5*math.sqrt(2), 0.0,             -0.5*math.sqrt(2), 0.0              ],
                    [-0.5*math.sqrt(2), 0.5*math.sqrt(2), 0.0,              0.0              ],
                    [-0.5*math.sqrt(2),-0.5*math.sqrt(2), 0.0,              0.0              ],
                  ]
    elif self.lattice == 'hexagonal':
      symQuats =  [
                    [ 1.0,0.0,0.0,0.0 ],
                    [-0.5*math.sqrt(3), 0.0, 0.0,-0.5 ],
                    [ 0.5, 0.0, 0.0, 0.5*math.sqrt(3) ],
                    [ 0.0,0.0,0.0,1.0 ],
                    [-0.5, 0.0, 0.0, 0.5*math.sqrt(3) ],
                    [-0.5*math.sqrt(3), 0.0, 0.0, 0.5 ],
                    [ 0.0,1.0,0.0,0.0 ],
                    [ 0.0,-0.5*math.sqrt(3), 0.5, 0.0 ],
                    [ 0.0, 0.5,-0.5*math.sqrt(3), 0.0 ],
                    [ 0.0,0.0,1.0,0.0 ],
                    [ 0.0,-0.5,-0.5*math.sqrt(3), 0.0 ],
                    [ 0.0, 0.5*math.sqrt(3), 0.5, 0.0 ],
                  ]
    elif self.lattice == 'tetragonal':
      symQuats =  [
                    [ 1.0,0.0,0.0,0.0 ],
                    [ 0.0,1.0,0.0,0.0 ],
                    [ 0.0,0.0,1.0,0.0 ],
                    [ 0.0,0.0,0.0,1.0 ],
                    [ 0.0, 0.5*math.sqrt(2), 0.5*math.sqrt(2), 0.0 ],
                    [ 0.0,-0.5*math.sqrt(2), 0.5*math.sqrt(2), 0.0 ],
                    [ 0.5*math.sqrt(2), 0.0, 0.0, 0.5*math.sqrt(2) ],
                    [-0.5*math.sqrt(2), 0.0, 0.0, 0.5*math.sqrt(2) ],
                  ]
    elif self.lattice == 'orthorhombic':
      symQuats =  [
                    [ 1.0,0.0,0.0,0.0 ],
                    [ 0.0,1.0,0.0,0.0 ],
                    [ 0.0,0.0,1.0,0.0 ],
                    [ 0.0,0.0,0.0,1.0 ],
                  ]
    else:
      symQuats =  [
                    [ 1.0,0.0,0.0,0.0 ],
                  ]
    return list(map(Quaternion,
               np.array(symQuats)[np.atleast_1d(np.array(who)) if who != [] else range(len(symQuats))]))


  def unitCell(self):
    """Return unit cell edges
    """
    if self.lattice == 'cubic':
      unitCell =  [[0,0,0,0,0,1], #z
                   [0,0,1,0,1,1], #y
                   [0,1,1,0,1,0], #z
                   [0,1,0,0,0,0], #y
                   [1,0,0,1,0,1], #z
                   [1,0,1,1,1,1], #y
                   [1,1,1,1,1,0], #z
                   [1,1,0,1,0,0], #y
                   [0,0,0,1,0,0], #x
                   [0,0,1,1,0,1], #x
                   [0,1,1,1,1,1], #x
                   [0,1,0,1,1,0]] #x
      unitCell = np.array(unitCell)-np.array([0.5,0.5,0.5,0.5,0.5,0.5])
    elif self.lattice == 'hexagonal':
      thh      = np.sqrt(3)/2
      unitCell =  [[-1.0, 0.0, 0.0,   -0.5,-thh, 0.0],
                   [-0.5,-thh, 0.0,    0.5,-thh, 0.0],
                   [ 0.5,-thh, 0.0,    1.0, 0.0, 0.0],
                   [-1.0, 0.0, 0.0,   -0.5, thh, 0.0],
                   [-0.5, thh, 0.0,    0.5, thh, 0.0],
                   [ 0.5, thh, 0.0,    1.0, 0.0, 0.0],
                   [-1.0, 0.0, 2.0,   -0.5,-thh, 2.0],
                   [-0.5,-thh, 2.0,    0.5,-thh, 2.0],
                   [ 0.5,-thh, 2.0,    1.0, 0.0, 2.0],
                   [-1.0, 0.0, 2.0,   -0.5, thh, 2.0],
                   [-0.5, thh, 2.0,    0.5, thh, 2.0],
                   [ 0.5, thh, 2.0,    1.0, 0.0, 2.0],
                   [-1.0, 0.0, 0.0,   -1.0, 0.0, 2.0],
                   [ 1.0, 0.0, 0.0,    1.0, 0.0, 2.0],
                   [ 0.5, thh, 0.0,    0.5, thh, 2.0],
                   [ 0.5,-thh, 0.0,    0.5,-thh, 2.0],
                   [-0.5,-thh, 0.0,   -0.5,-thh, 2.0],
                   [-0.5, thh, 0.0,   -0.5, thh, 2.0],
                   ]
      unitCell = list(np.array(unitCell)/2.)
    else:
      print("Unit cell not implemented")
      unitCell = [[None]]
    return unitCell


  def equivalentQuaternions(self, quaternion, who = []):
    '''
    List of symmetrically equivalent quaternions based on own symmetry.
    '''
    return [quaternion*q for q in self.symmetryQuats(who)]


  def inFZ(self,R):
    '''
    Check whether given Rodrigues vector falls into fundamental zone of own symmetry.
    '''
    if isinstance(R, Quaternion): R = R.asRodrigues()      # translate accidentially passed quaternion
    R = abs(R)                                             # fundamental zone in Rodrigues space is point symmetric around origin
    if self.lattice == 'cubic':
      return     math.sqrt(2.0)-1.0 >= R[0] \
             and math.sqrt(2.0)-1.0 >= R[1] \
             and math.sqrt(2.0)-1.0 >= R[2] \
             and 1.0 >= R[0] + R[1] + R[2]
    elif self.lattice == 'hexagonal':
      return     1.0 >= R[0] and 1.0 >= R[1] and 1.0 >= R[2] \
             and 2.0 >= math.sqrt(3)*R[0] + R[1] \
             and 2.0 >= math.sqrt(3)*R[1] + R[0] \
             and 2.0 >= math.sqrt(3) + R[2]
    elif self.lattice == 'tetragonal':
      return     1.0 >= R[0] and 1.0 >= R[1] \
             and math.sqrt(2.0) >= R[0] + R[1] \
             and math.sqrt(2.0) >= R[2] + 1.0
    elif self.lattice == 'orthorhombic':
      return     1.0 >= R[0] and 1.0 >= R[1] and 1.0 >= R[2]
    else:
      return True


  def inDisorientationSST(self,R):
    '''
    Check whether given Rodrigues vector (of misorientation) falls into standard stereographic triangle of own symmetry.

    Determination of disorientations follow the work of A. Heinz and P. Neumann:
    Representation of Orientation and Disorientation Data for Cubic, Hexagonal, Tetragonal and Orthorhombic Crystals
    Acta Cryst. (1991). A47, 780-789
    '''
    if isinstance(R, Quaternion): R = R.asRodrigues()       # translate accidentially passed quaternion
    epsilon = 0.0
    if self.lattice == 'cubic':
      return R[0] >= R[1]+epsilon                and R[1] >= R[2]+epsilon    and R[2] >= epsilon
    elif self.lattice == 'hexagonal':
      return R[0] >= math.sqrt(3)*(R[1]-epsilon) and R[1] >= epsilon         and R[2] >= epsilon
    elif self.lattice == 'tetragonal':
      return R[0] >= R[1]-epsilon                and R[1] >= epsilon         and R[2] >= epsilon
    elif self.lattice == 'orthorhombic':
      return R[0] >= epsilon                     and R[1] >= epsilon         and R[2] >= epsilon
    else:
      return True


  def inSST(self, vector, proper = False,  color = False):
    '''
    Check whether given vector falls into standard stereographic triangle of own symmetry.

    Args:
       vector: hkl vector tested/converted to rbg; considers only this vector, no crystal orientation
       proper: considers only vectors with z >= 0, hence uses two neighboring SSTs to determine if in SST. i.e. allows more positive results, rgb-value does not depend on this
       color: if true, return also color

    Returns:
       whether in SST
       (if requested) inverse pole figure color if in SST in RGB values
    '''
#     basis = {'cubic' :        np.linalg.inv(np.array([[0.,0.,1.],                                 # direction of red
#                                                       [1.,0.,1.]/np.sqrt(2.),                     # direction of green
#                                                       [1.,1.,1.]/np.sqrt(3.)]).transpose()),      # direction of blue
#              'hexagonal' :    np.linalg.inv(np.array([[0.,0.,1.],                                 # direction of red
#                                                       [1.,0.,0.],                                 # direction of green
#                                                       [np.sqrt(3.),1.,0.]/np.sqrt(4.)]).transpose()),      # direction of blue
#              'tetragonal' :   np.linalg.inv(np.array([[0.,0.,1.],                                 # direction of red
#                                                       [1.,0.,0.],                                 # direction of green
#                                                       [1.,1.,0.]/np.sqrt(2.)]).transpose()),      # direction of blue
#              'orthorhombic' : np.linalg.inv(np.array([[0.,0.,1.],                                 # direction of red
#                                                       [1.,0.,0.],                                 # direction of green
#                                                       [0.,1.,0.]]).transpose()),                  # direction of blue
#             }
    if self.lattice == 'cubic':
      basis = {'improper':np.array([ [-1.            ,  0.            ,  1. ],
                                   [ np.sqrt(2.)   , -np.sqrt(2.)   ,  0. ],
                                   [ 0.            ,  np.sqrt(3.)   ,  0. ] ]),
             'proper':np.array([ [ 0.            , -1.            ,  1. ],
                                   [-np.sqrt(2.)   , np.sqrt(2.)    ,  0. ],
                                   [ np.sqrt(3.)   ,  0.            ,  0. ] ]),
              }
    elif self.lattice == 'hexagonal':
      basis = {'improper':np.array([ [ 0.            ,  0.            ,  1. ],
                                   [ 1.            , -np.sqrt(3.),     0. ],
                                   [ 0.            ,  2.            ,  0. ] ]),
             'proper':np.array([ [ 0.            ,  0.            ,  1. ],
                                   [-1.            ,  np.sqrt(3.)   ,  0. ],
                                   [ np.sqrt(3)    , -1.            ,  0. ] ]),
              }
    elif self.lattice == 'tetragonal':
      basis = {'improper':np.array([ [ 0.            ,  0.            ,  1. ],
                                   [ 1.            , -1.            ,  0. ],
                                   [ 0.            ,  np.sqrt(2.),     0. ] ]),
             'proper':np.array([ [ 0.            ,  0.            ,  1. ],
                                   [-1.            ,  1.            ,  0. ],
                                   [ np.sqrt(2.)   ,  0.            ,  0. ] ]),
              }
    elif self.lattice == 'orthorhombic':
      basis = {'improper':np.array([ [ 0., 0., 1.],
                                   [ 1., 0., 0.],
                                   [ 0., 1., 0.] ]),
             'proper':np.array([ [ 0., 0., 1.],
                                   [-1., 0., 0.],
                                   [ 0., 1., 0.] ]),
              }
    else:
      basis = {'improper':np.zeros((3,3),dtype=float),
             'proper':np.zeros((3,3),dtype=float),
              }
    #theComponents = color components
    #inSST = vector of true/false for each data-point
    if np.all(basis == 0.0):  #default in basis not set
      theComponents = -np.ones(3,'d')
      inSST = False   #np.all(theComponents >= 0.0), this statement is always false
    else:
      v = np.array(vector,dtype = float)
      if proper:                                                                                  # check both improper ...
        theComponents = np.dot(basis['improper'],v)
        inSST = np.all(theComponents >= 0.0, axis=0)
        if not np.any(inSST):                                                                               # ... and proper SST
          theComponents = np.dot(basis['proper'],v)
      else:
        v[2] = abs(v[2])                                                                            # z component projects identical for positive and negative values
        theComponents = np.dot(basis['improper'],v)
      inSST = np.all(theComponents >= 0.0, axis=0)
    if color:                                                                                       # have to return color array
      if np.any(inSST):
        theComponentsNorm = theComponents/np.linalg.norm(theComponents, axis=0)
        rgb = np.power(np.abs(theComponentsNorm),0.5)                             # smoothen color ramps
        if rgb.ndim>1: rgb[:,~inSST] = 0.0
        rgb[rgb>1.0] = 1.0                                                        # limit to maximum intensity
        if rgb.ndim>1: rgb[:,inSST] /= np.max(rgb[:,inSST], axis=0)                                                                             # normalize to (HS)V = 1
        else         : rgb          /= np.max(rgb)
      else:
        rgb = np.zeros(3,'d')
      return (inSST,rgb)
    else:
      return inSST


  #@}
  ##
  # @name Math routines to convert between HKL, XY, RGB
  #@{

  def xyToHKL(self, inPlane):
      """
      convert xy (coordinates in projection plane (-1<x<1,-1<y<1)) -> hkl

      required for plotting the standardTriangle

      Theory:
      - unscaled vector (x; y; 1) has to be scaled by l such that the scaled vector - unit_vector_Z has the unit length
      - scaled: ( xl; yl; l)
      - |(xl; yl; l-1)| =!=1
      - -> l(x^2+y^2+1)=2 #2 is a random number
      - hkl = scaled - unit_vector_Z = (xl; yl; l-1)

      Args:
         inPlane: vector coordinates in the plane of the stereograhic projection plane

      Returns:
         hkl vector
      """
      if self.lattice != 'cubic':
        print("ERROR: only implemented for cubic lattice")
        return
      l = 2.0/ ( inPlane[0]*inPlane[0] + inPlane[1]*inPlane[1] + 1.)
      if inPlane.ndim==1:
        print(l)
        hkl =np.zeros((3), dtype=np.float)
        hkl[:2]=l*inPlane
        hkl[2]=l-1.0
      else:
        hkl = np.zeros((3,inPlane.shape[1]), dtype=np.float)
        hkl[:2,:]=l*inPlane
        hkl[2,:]=l-1.0
      return hkl


  def standardTriangle(self, fileName=None, show=True, stepSize = 0.01):
    """
    Plot standard triangle with background, discrete points in color, save to file, add text

    Corner points of triangle
    - x=0.41421 =sqrt(2)-1 -> [110] corner at (0.41421,0)
    - [111] corner at (0.366025403784,0.366025403784)

    Args:
       fileName: save to file
       show: True [default] shows figure, else not
       stepSize: plotting accuracy: lower value=better quality
    """
    from scipy.interpolate import interp1d
    if self.lattice != 'cubic':
      print("ERROR: only implemented for cubic lattice")
      return
    #create border
    border = [[0,0]]
    xTemp = []
    yTemp = []
    for a in range(16):
      border.append([math.cos(a/180.0*math.pi)*math.sqrt(2.0)-1,math.sin(a/180.0*math.pi)*math.sqrt(2) ])
      xTemp.append(math.cos(a/180.0*math.pi)*math.sqrt(2.0)-1)
      yTemp.append(math.sin(a/180.0*math.pi)*math.sqrt(2))
    border.append([0,0])
    border = np.array(border)
    func = interp1d(yTemp,xTemp)  #function yTemp=func(xTemp)
    #create colored background
    xy, rgb = [], []
    for iY in np.arange(0,0.366025403784,stepSize):
      for iX in np.arange(iY,0.41421,stepSize):
        if func(iY)<=iX: continue
        xy.append( [iX,iY])
    xy = np.array(xy)
    hkl=self.xyToHKL(xy.T)
    flags, rgb = self.inSST(hkl, color=True, proper=False)
    colors = []
    for i in range(rgb.shape[1]):
      values = rgb[:,i]*255
      string = '#%02x%02x%02x' % (values[0], values[1], values[2])
      colors.append(string)
    #plotting: background, border, labels
    plt.scatter(xy[:,0], xy[:,1], c=colors, s=15000.*stepSize, linewidths=0)#, alpha=0.05)
    plt.plot(border[:,0],border[:,1],'-k',linewidth=3)#, alpha=0.5)
    plt.rcParams['font.size'] = 18.
    plt.text(0,0,'[100]',horizontalalignment='right')  #, zorder=40
    plt.text(0.42,0,'[110]')
    plt.text(0.37,0.37,'[111]')
    plt.axis('off')
    plt.axis('equal')
    if fileName:
      plt.savefig(fileName)
    if show:
      plt.show()
    return

  #@}
