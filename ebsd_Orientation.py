# -*- coding: utf-8 -*-
##
# @file
# @brief Orientation class: combination of material symmetry and specific rotation.
# Copyright:
#   Original version was part of DAMASK <damask.mpie.de> (Martin Diehl, Philip Eisenlohr, Franz Roters)
#
import numpy as np
import random
from ebsd_Symmetry import Symmetry
from ebsd_Quaternion import Quaternion


class Orientation:
  """Orientation class: combination of material symmetry and specific rotation
  """
  __slots__ = ['quaternion', 'symmetry', 'plot2D', 'eps', 'doctest']

  # @name CONVENTIONAL ROUTINES
  #@{

  def __init__(self,
               quaternion = Quaternion.fromIdentity(),
               Rodrigues  = None,
               angleAxis  = None,
               matrix     = None,
               Eulers     = None,
               random     = False,                                                                  # put any integer to have a fixed seed or True for real random
               symmetry   = None,
              ):
    if random:                                                                                      # produce random orientation
      if isinstance(random, bool ):
        self.quaternion = Quaternion.fromRandom()
      else:
        self.quaternion = Quaternion.fromRandom(randomSeed=random)
    elif isinstance(Eulers, np.ndarray) and Eulers.shape == (3,):                                   # based on given Euler angles
      self.quaternion = Quaternion.fromEulers(Eulers,'bunge')
    elif isinstance(matrix, np.ndarray):                                                            # based on given rotation matrix
      self.quaternion = Quaternion.fromMatrix(matrix)
    elif isinstance(angleAxis, np.ndarray) and angleAxis.shape == (4,):                             # based on given angle and rotation axis
      self.quaternion = Quaternion.fromAngleAxis(angleAxis[0],angleAxis[1:4])
    elif isinstance(Rodrigues, np.ndarray) and Rodrigues.shape == (3,):                             # based on given Rodrigues vector
      self.quaternion = Quaternion.fromRodrigues(Rodrigues)
    elif isinstance(quaternion, Quaternion):                                                        # based on given quaternion
      self.quaternion = quaternion.homomorphed()
    elif isinstance(quaternion, np.ndarray) and quaternion.shape == (4,):                           # based on given quaternion
      self.quaternion = Quaternion(quaternion).homomorphed()
    self.symmetry = Symmetry(symmetry)
    self.plot2D = "down-right"
    self.eps = 1e-6
    self.doctest = False
    return



  def __copy__(self):
    return self.__class__(quaternion=self.quaternion,symmetry=self.symmetry.lattice)
  copy = __copy__


  def __repr__(self):
    return 'Symmetry: %s\n' % (self.symmetry) + \
           'Quaternion: %s\n' % (self.quaternion) + \
           'Matrix:\n%s\n' % ( '\n'.join(['\t'.join(map(str,self.asMatrix()[i,:])) for i in range(3)]) ) + \
           'Bunge Eulers / deg: %s' % ('\t'.join(map(str,self.asEulers('bunge',degrees=True))) )


  def asQuaternion(self):
    return self.quaternion.asList()


  def asEulers(self,
               notation = 'bunge',
               degrees = False,
               standardRange = False):
    return self.quaternion.asEulers(notation, degrees, standardRange)
  eulers = property(asEulers)


  def asRodrigues(self):
    return self.quaternion.asRodrigues()
  rodrigues = property(asRodrigues)


  def asAngleAxis(self,
                  degrees = False):
    return self.quaternion.asAngleAxis(degrees)
  angleAxis = property(asAngleAxis)


  def asMatrix(self):
    return self.quaternion.asMatrix()
  matrix = property(asMatrix)


  def inFZ(self):
    """Check whether given Rodrigues vector falls into fundamental zone of own symmetry.
    """
    return self.symmetry.inFZ(self.quaternion.asRodrigues())
  infz = property(inFZ)


  def equivalentQuaternions(self,
                            who = []):
    return self.symmetry.equivalentQuaternions(self.quaternion,who)


  def equivalentOrientations(self,
                             who = []):
    return [Orientation(quaternion = q, symmetry = self.symmetry.lattice) for q in self.equivalentQuaternions(who)]


  def reduced(self):
    '''
    Transform orientation to fall into fundamental zone according to symmetry
    '''
    for me in self.symmetry.equivalentQuaternions(self.quaternion):
      if self.symmetry.inFZ(me.asRodrigues()): break
    return Orientation(quaternion=me,symmetry=self.symmetry.lattice)


  #@}
  ##
  # @name MATERIAL SPECIFIC ROUTINES
  #@{
  def disorientation(self, other, SST = True):
    """Disorientation between myself and given other orientation.

    Rotation axis falls into SST if SST == True.
      (Currently requires same symmetry for both orientations.
      Look into A. Heinz and P. Neumann 1991 for cases with differing sym.)

     Args:
      other: other orientation
      SST: True (rotation axis falls into SST); False

    Returns:
      disorientation quaternion, idx of equivalent orientation1, idx of equivalent orientation2, num. of conjugated
    """
    if self.symmetry != other.symmetry: raise TypeError('disorientation between different symmetry classes not supported yet.')
    misQ = self.quaternion.conjugated()*other.quaternion
    mySymQs    =  self.symmetry.symmetryQuats() if SST else self.symmetry.symmetryQuats()[:1]       # take all or only first sym operation
    otherSymQs = other.symmetry.symmetryQuats()
    for i,sA in enumerate(mySymQs):  #if not in SST: only one sA
      for j,sB in enumerate(otherSymQs): #changes always
        theQ = sA.conjugated()*misQ*sB
        for k in range(2):
          theQ.conjugate()
          breaker = self.symmetry.inFZ(theQ) and (not SST or other.symmetry.inDisorientationSST(theQ))
          if breaker: break
        if breaker: break
      if breaker: break
    return (Orientation(quaternion = theQ,symmetry = self.symmetry.lattice),
            i,j,k == 1)                                                                             # disorientation, own sym, other sym, self-->other: True, self<--other: False


  def inversePole(self, axis, proper = False, SST = True):
    """axis rotated according to orientation (using crystal symmetry to ensure location falls into SST)

    Args:
      axis: vector in crystal orientation, e.g. [100]
      proper: considers only vectors with z >= 0, hence uses two neighboring SSTs to determine if in SST. i.e. allows more positive results, rgb-value does not depend on this
      SST: iterate through all equivalent and find the one in the SST

    Returns:
      vector of axis
    """
    if SST:                                                                                         # pole requested to be within SST
      for i,q in enumerate(self.symmetry.equivalentQuaternions(self.quaternion)):                   # test all symmetric equivalent quaternions
        pole = q.conjugated()*axis                                                                  # align crystal direction to axis
        if self.symmetry.inSST(pole,proper): break                                                # found SST version
    else:
      pole = self.quaternion.conjugated()*axis                                                      # align crystal direction to axis
    return (pole,i if SST else 0)



  def IPFcolor(self,axis, proper=False):
    """color of inverse pole figure for given axis

    Args:
       axis: axis of pole figure (ND=001)
       proper: considers only vectors with z >= 0, hence uses two neighboring SSTs to determine if in SST. i.e. allows more positive results, rgb-value does not depend on this

    Returns:
       vector of color (rgb)
    """
    color = np.zeros(3,'d')
    for q in self.symmetry.equivalentQuaternions(self.quaternion):
      pole = q.conjugated()*axis                                                                    # align crystal direction to axis
      inSST,color = self.symmetry.inSST(pole,color=True,proper=proper)
      if inSST: break
    return color


  @classmethod
  def average(cls,
              orientations,
              multiplicity = []):
    """Return the average orientation

    ref: F. Landis Markley, Yang Cheng, John Lucas Crassidis, and Yaakov Oshman,
      Averaging Quaternions,
      Journal of Guidance, Control, and Dynamics, Vol. 30, No. 4 (2007), pp. 1193-1197.
      doi: 10.2514/1.28949

    Usage:
      * a = Orientation(Eulers=np.radians([10, 10, 0]), symmetry='hexagonal')
      * b = Orientation(Eulers=np.radians([20, 0, 0]),  symmetry='hexagonal')
      * avg = Orientation.average([a,b])

    Args:
      cls: class method (void)
      orientations: list of orientations
      multiplicity: --

    Returns:
      average orientation (not rotation) in radians
    """
    if not all(isinstance(item, Orientation) for item in orientations):
      raise TypeError("Only instances of Orientation can be averaged.")
    N = len(orientations)
    if multiplicity == [] or not multiplicity:
      multiplicity = np.ones(N,dtype='i')
    reference = orientations[0]                                                                     # take first as reference
    for i,(o,n) in enumerate(zip(orientations,multiplicity)):
      closest = o.equivalentOrientations(reference.disorientation(o,SST = False)[2])[0]             # select sym orientation with lowest misorientation
      M = closest.quaternion.asM() * n if i == 0 else M + closest.quaternion.asM() * n              # add (multiples) of this orientation to average
    eig, vec = np.linalg.eig(M/N)
    return Orientation(quaternion = Quaternion(quatArray = np.real(vec.T[eig.argmax()])),
                       symmetry = reference.symmetry.lattice)

  #@}
  ##
  # @name PLOTTING, PRINTING
  #@{
  def project(self, x,y,z):
    """

    down-right: y, -x
    up-left   : -y, x
    right-up   : x,y
    left-down: -x,-y
    3D        : x,y,z
    """
    if type(x) == list:
      x,y,z = np.array(x),np.array(y),np.array(z)
      if   self.plot2D=='down-right':  return  y,-x
      elif self.plot2D=='up-left'   :  return -y, x
      elif self.plot2D=='right-up'  :  return  x, y
      elif self.plot2D=='left-down' :  return -x,-y
      elif self.plot2D=='3D':          return  x, y, z
      else:   print("Error: plot2D not well defined: plotLine")
    else:
      if   self.plot2D=='down-right':  return  y,-x
      elif self.plot2D=='up-left'   :  return -y, x
      elif self.plot2D=='right-up'  :  return  x, y
      elif self.plot2D=='left-down' :  return -x,-y
      elif self.plot2D=='3D':          return  x, y, z
      else:   print("Error: plot2D not well defined: plotLine")
    return


  def plotLine(self, ax, start,delta,color='k',lw=1, ls='solid', markerSize=None):
    """
    Plot one line using given projection

    Args:
       ax: axis to plot into
       start: start coordinate
       delta: delta cooradinate (end-start)
       color: color
       lw: line width
       ls: line style "solid",'dashed'
       markerSize: size of marker (only used for non-lines: delta>0)
    """
    if np.linalg.norm(delta)<self.eps:
      marker='o'
      if markerSize is None:
        markerSize=7
    else:
      marker=None; markerSize=0
    if self.plot2D=='3D':
      ax.plot( [start[0]]+[start[0]+delta[0]],
               [start[1]]+[start[1]+delta[1]],
               [start[2]]+[start[2]+delta[2]],
               color=color,lw=lw, marker=marker, ls=ls, markersize=markerSize)
    else:
      x,y = self.project([start[0]]+[start[0]+delta[0]],
                    [start[1]]+[start[1]+delta[1]],
                    [start[2]]+[start[2]+delta[2]])
      ax.plot(x,y,  color=color,lw=lw, marker=marker, ls=ls, markersize=markerSize)
    return


  def plotUnit(self, ax, xlabel,ylabel,zlabel,  x=0,y=0,z=0,   s=1):  # unit axis
    """
    Coordinate systems: see plotLine

    Args:
       ax: axis to used for plotting
       xlabel: x-label
       ylabel: y-label
       zlabel: z-label
       x: x-coordinate of origin
       y: y-coordinate of origin
       z: z-coordinate of origin
       s: scale
    """
    self.plotLine(ax, [x,y,z], [s,0,0], 'k',lw=3)
    self.plotLine(ax, [x,y,z], [0,s,0], 'k',lw=3)
    self.plotLine(ax, [x,y,z], [0,0,s], 'k',lw=3)
    if self.plot2D=='3D':
      ax.text(x+s,   y+0.1, z+0.1,xlabel )
      ax.text(x+0.1, y+s,   z+0.1,ylabel )
      ax.text(x+0.1, y+0.1, z+s  ,zlabel )
    else:
      ax.text(*(self.project(x+s,   y  , z+0.1)+(xlabel,)) )
      ax.text(*(self.project(x+0.1, y+s, z+0.1)+(ylabel,{"ha":"right"})) )
      ax.text(*(self.project(x+0.1, y  , z+s  )+(zlabel,)) )
    return


  def plot(self, poles=None, unitCell=True, cos=True, annotate=False, plot2D=None, scale=2, fileName=None):
    """Plot rotated unit-cell in 3D, and possibly the pole-figure and specific poles

    Projection onto 2D: cooradinate systems are given as xDirection-yDirection (z follows)
    - down-right: [default in text books] RD = x = down; TD = y = right; ND = z = outOfPlane
    - up-left: [default in OIM] RD = x = up; TD = y = left; ND = z = outOfPlane

    Args:
       poles: if given (e.g. [1,0,0]), plot pole-figure and the corresponding poles
       unitCell: plot unit cell
       cos: plot coordinate system
       annotate: annotate poles in pole figure (requires poles given)
       plot2D: do a normal projection onto 2D plane: [down-right, up-left, None]
       scale: scale of pole-figure dome over crystal
       fileName: fileName for image output (if given, image not shown)
    """
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    if plot2D is not None:
      self.plot2D = plot2D
    if self.plot2D == "3D":
      fig = plt.figure()
      ax = fig.gca(projection='3d')
      #ax.view_init(90,0)
    else:
      fig, ax = plt.subplots()

    if unitCell:
      for line in self.symmetry.unitCell():
        start, end = np.array(line[:3],dtype=np.float), np.array(line[3:],dtype=np.float)
        start = self.quaternion*start
        end   = self.quaternion*end
        if start[2]<0 and end[2]<0:
          self.plotLine(ax, start, end-start,color="b",lw=0.2)
        elif start[2]>0 and end[2]>0:
          self.plotLine(ax, start, end-start,color="b",lw=2)
        else:
          delta = end-start
          k     = -start[2]/delta[2]
          mid   = start+k*delta
          if start[2]>0:
            self.plotLine(ax, start, mid-start,color="b",lw=2)
            self.plotLine(ax, mid,   end-mid,color="b",lw=0.2)
          else:
            self.plotLine(ax, start, mid-start,color="b",lw=0.2)
            self.plotLine(ax, mid,   end-mid,color="b",lw=2)

    if cos: self.plotUnit(ax, "RD [100]","TD [010]","ND [001]")

    if poles is not None:
      #plot sphere
      if self.plot2D=="3D":
        u = np.linspace(0,2*np.pi,50)
        v = np.linspace(0,np.pi/2,50)
        x = scale*np.outer(np.cos(u)      , np.sin(v))
        y = scale*np.outer(np.sin(u)      , np.sin(v))
        z = scale*np.outer(np.ones_like(u), np.cos(v))
        ax.plot_surface(x, y, z, color='gray',alpha=0.7, cstride=10, rstride=10, lw=0)
      else:
        ax.plot( scale*np.cos(np.linspace(0.,2.*np.pi,100)), scale*np.sin(np.linspace(0.,2.*np.pi,100)), 'k--' )
      # plot poles
      oHelp = Orientation(Eulers=np.array([0.,0.,0.]), symmetry=self.symmetry.__repr__())
      poles = np.array(poles, dtype=np.float)
      poles /= np.linalg.norm(poles)
      for idx,q in enumerate(oHelp.symmetry.equivalentQuaternions(oHelp.quaternion)):
        conjAxis      = q*poles  #e.g. [100]
        direction = self.quaternion*conjAxis
        if direction[2]<-self.eps: continue                                                        #prevent rounding errors
        fromBase = direction+np.array([0,0,1])
        #self.plotLine(ax, [0,0,0], direction*scale, color='c', lw=1) #in plane lines: not needed
        if self.plot2D=="3D":
          self.plotLine(ax, -scale*np.array([0,0,1]), (fromBase)*scale, 'c')  #lines from bottom base to points
        xy  = fromBase/fromBase[2]*scale
        xy[2]=0.0
        self.plotLine(ax, xy, [0.,0.,0.], color='c')  #plot point
        if annotate:
          x_,y_ = self.project(xy[0],xy[1],xy[2])
          label_ = str(np.array(conjAxis,dtype=np.int))[1:-1]
          label_ = label_.replace(" ","")
          ax.text(x_+0.05,y_+0.05, label_)

    #finalize plot
    ax.axis('equal'); ax.axis('off')
    ax.set_xlim([-scale*1.1,scale*1.1])
    ax.set_ylim([-scale*1.1,scale*1.1])
    if self.plot2D=="3D":
      ax.set_zlabel(''); ax.set_zticks([])
    if self.doctest: plt.savefig('doctest.png'); plt.close(); return
    if fileName is not None:
      plt.savefig(fileName, dpi=150, bbox_inches='tight')
    else:
      plt.show()
    return


  def toScreen(self, equivalent=True):
    """
    print Euler angles and HKL /UVW

    Args:
      equivalent: print also equivalent orientations
    """
    print("Euler angles:",np.round(self.quaternion.asEulers(degrees=True),1))
    rotM = self.quaternion.asMatrix()
    print("HKL",np.array( rotM[2,:]/np.min(rotM[2,:]) ,dtype=np.int))
    print("UVW",-np.array( rotM[0,:]/np.min(rotM[0,:]) ,dtype=np.int))
    if equivalent:
      print("Equivalent orientations - Euler angles:")
      for q in self.symmetry.equivalentQuaternions(self.quaternion):
        angles = q.asEulers(degrees=True)
        angles[angles<0] += 360.
        print("   [%5.1f  %5.1f  %5.1f]"%tuple(angles))
    return



  #@}
  ##
  # @name MISC
  #@{
  def related(self,
              relationModel,
              direction,
              targetSymmetry = None):
    """Related

    Models:
      * KS from S. Morito et al./Journal of Alloys and Compounds 5775 (2013) S587-S592 DOES THIS PAPER EXISTS?
      * GT from Y. He et al./Journal of Applied Crystallography (2006). 39, 72-81
      * GT' from Y. He et al./Journal of Applied Crystallography (2006). 39, 72-81
      * NW from H. Kitahara et al./Materials Characterization 54 (2005) 378-386
      * Pitsch from Y. He et al./Acta Materialia 53 (2005) 1179-1190
      * Bain from Y. He et al./Journal of Applied Crystallography (2006). 39, 72-81

    Args:
      relationModel: --
      direction: --
      targetSymmetry: --

    Returns:
      vector
    """
    if relationModel not in ['KS','GT','GTdash','NW','Pitsch','Bain']:  return None
    if int(direction) == 0:  return None
    variant  = int(abs(direction))-1
    (me,other)  = (0,1) if direction > 0 else (1,0)
    planes = {'KS': \
                    np.array([[[  1,  1,  1],[  0,  1,  1]],
                              [[  1,  1,  1],[  0,  1,  1]],
                              [[  1,  1,  1],[  0,  1,  1]],
                              [[  1,  1,  1],[  0,  1,  1]],
                              [[  1,  1,  1],[  0,  1,  1]],
                              [[  1,  1,  1],[  0,  1,  1]],
                              [[  1, -1,  1],[  0,  1,  1]],
                              [[  1, -1,  1],[  0,  1,  1]],
                              [[  1, -1,  1],[  0,  1,  1]],
                              [[  1, -1,  1],[  0,  1,  1]],
                              [[  1, -1,  1],[  0,  1,  1]],
                              [[  1, -1,  1],[  0,  1,  1]],
                              [[ -1,  1,  1],[  0,  1,  1]],
                              [[ -1,  1,  1],[  0,  1,  1]],
                              [[ -1,  1,  1],[  0,  1,  1]],
                              [[ -1,  1,  1],[  0,  1,  1]],
                              [[ -1,  1,  1],[  0,  1,  1]],
                              [[ -1,  1,  1],[  0,  1,  1]],
                              [[  1,  1, -1],[  0,  1,  1]],
                              [[  1,  1, -1],[  0,  1,  1]],
                              [[  1,  1, -1],[  0,  1,  1]],
                              [[  1,  1, -1],[  0,  1,  1]],
                              [[  1,  1, -1],[  0,  1,  1]],
                              [[  1,  1, -1],[  0,  1,  1]]]),
              'GT': \
                    np.array([[[  1,  1,  1],[  1,  0,  1]],
                              [[  1,  1,  1],[  1,  1,  0]],
                              [[  1,  1,  1],[  0,  1,  1]],
                              [[ -1, -1,  1],[ -1,  0,  1]],
                              [[ -1, -1,  1],[ -1, -1,  0]],
                              [[ -1, -1,  1],[  0, -1,  1]],
                              [[ -1,  1,  1],[ -1,  0,  1]],
                              [[ -1,  1,  1],[ -1,  1,  0]],
                              [[ -1,  1,  1],[  0,  1,  1]],
                              [[  1, -1,  1],[  1,  0,  1]],
                              [[  1, -1,  1],[  1, -1,  0]],
                              [[  1, -1,  1],[  0, -1,  1]],
                              [[  1,  1,  1],[  1,  1,  0]],
                              [[  1,  1,  1],[  0,  1,  1]],
                              [[  1,  1,  1],[  1,  0,  1]],
                              [[ -1, -1,  1],[ -1, -1,  0]],
                              [[ -1, -1,  1],[  0, -1,  1]],
                              [[ -1, -1,  1],[ -1,  0,  1]],
                              [[ -1,  1,  1],[ -1,  1,  0]],
                              [[ -1,  1,  1],[  0,  1,  1]],
                              [[ -1,  1,  1],[ -1,  0,  1]],
                              [[  1, -1,  1],[  1, -1,  0]],
                              [[  1, -1,  1],[  0, -1,  1]],
                              [[  1, -1,  1],[  1,  0,  1]]]),
              'GTdash': \
                    np.array([[[  7, 17, 17],[ 12,  5, 17]],
                              [[ 17,  7, 17],[ 17, 12,  5]],
                              [[ 17, 17,  7],[  5, 17, 12]],
                              [[ -7,-17, 17],[-12, -5, 17]],
                              [[-17, -7, 17],[-17,-12,  5]],
                              [[-17,-17,  7],[ -5,-17, 12]],
                              [[  7,-17,-17],[ 12, -5,-17]],
                              [[ 17, -7,-17],[ 17,-12, -5]],
                              [[ 17,-17, -7],[  5,-17,-12]],
                              [[ -7, 17,-17],[-12,  5,-17]],
                              [[-17,  7,-17],[-17, 12, -5]],
                              [[-17, 17, -7],[ -5, 17,-12]],
                              [[  7, 17, 17],[ 12, 17,  5]],
                              [[ 17,  7, 17],[  5, 12, 17]],
                              [[ 17, 17,  7],[ 17,  5, 12]],
                              [[ -7,-17, 17],[-12,-17,  5]],
                              [[-17, -7, 17],[ -5,-12, 17]],
                              [[-17,-17,  7],[-17, -5, 12]],
                              [[  7,-17,-17],[ 12,-17, -5]],
                              [[ 17, -7,-17],[ 5, -12,-17]],
                              [[ 17,-17,  7],[ 17, -5,-12]],
                              [[ -7, 17,-17],[-12, 17, -5]],
                              [[-17,  7,-17],[ -5, 12,-17]],
                              [[-17, 17, -7],[-17,  5,-12]]]),
              'NW': \
                    np.array([[[  1,  1,  1],[  0,  1,  1]],
                              [[  1,  1,  1],[  0,  1,  1]],
                              [[  1,  1,  1],[  0,  1,  1]],
                              [[ -1,  1,  1],[  0,  1,  1]],
                              [[ -1,  1,  1],[  0,  1,  1]],
                              [[ -1,  1,  1],[  0,  1,  1]],
                              [[  1, -1,  1],[  0,  1,  1]],
                              [[  1, -1,  1],[  0,  1,  1]],
                              [[  1, -1,  1],[  0,  1,  1]],
                              [[ -1, -1,  1],[  0,  1,  1]],
                              [[ -1, -1,  1],[  0,  1,  1]],
                              [[ -1, -1,  1],[  0,  1,  1]]]),
              'Pitsch': \
                    np.array([[[  0,  1,  0],[ -1,  0,  1]],
                              [[  0,  0,  1],[  1, -1,  0]],
                              [[  1,  0,  0],[  0,  1, -1]],
                              [[  1,  0,  0],[  0, -1, -1]],
                              [[  0,  1,  0],[ -1,  0, -1]],
                              [[  0,  0,  1],[ -1, -1,  0]],
                              [[  0,  1,  0],[ -1,  0, -1]],
                              [[  0,  0,  1],[ -1, -1,  0]],
                              [[  1,  0,  0],[  0, -1, -1]],
                              [[  1,  0,  0],[  0, -1,  1]],
                              [[  0,  1,  0],[  1,  0, -1]],
                              [[  0,  0,  1],[ -1,  1,  0]]]),
              'Bain': \
                    np.array([[[  1,  0,  0],[  1,  0,  0]],
                              [[  0,  1,  0],[  0,  1,  0]],
                              [[  0,  0,  1],[  0,  0,  1]]]),
              }
    normals = {'KS': \
                    np.array([[[ -1,  0,  1],[ -1, -1,  1]],
                              [[ -1,  0,  1],[ -1,  1, -1]],
                              [[  0,  1, -1],[ -1, -1,  1]],
                              [[  0,  1, -1],[ -1,  1, -1]],
                              [[  1, -1,  0],[ -1, -1,  1]],
                              [[  1, -1,  0],[ -1,  1, -1]],
                              [[  1,  0, -1],[ -1, -1,  1]],
                              [[  1,  0, -1],[ -1,  1, -1]],
                              [[ -1, -1,  0],[ -1, -1,  1]],
                              [[ -1, -1,  0],[ -1,  1, -1]],
                              [[  0,  1,  1],[ -1, -1,  1]],
                              [[  0,  1,  1],[ -1,  1, -1]],
                              [[  0, -1,  1],[ -1, -1,  1]],
                              [[  0, -1,  1],[ -1,  1, -1]],
                              [[ -1,  0, -1],[ -1, -1,  1]],
                              [[ -1,  0, -1],[ -1,  1, -1]],
                              [[  1,  1,  0],[ -1, -1,  1]],
                              [[  1,  1,  0],[ -1,  1, -1]],
                              [[ -1,  1,  0],[ -1, -1,  1]],
                              [[ -1,  1,  0],[ -1,  1, -1]],
                              [[  0, -1, -1],[ -1, -1,  1]],
                              [[  0, -1, -1],[ -1,  1, -1]],
                              [[  1,  0,  1],[ -1, -1,  1]],
                              [[  1,  0,  1],[ -1,  1, -1]]]),
              'GT': \
                    np.array([[[ -5,-12, 17],[-17, -7, 17]],
                              [[ 17, -5,-12],[ 17,-17, -7]],
                              [[-12, 17, -5],[ -7, 17,-17]],
                              [[  5, 12, 17],[ 17,  7, 17]],
                              [[-17,  5,-12],[-17, 17, -7]],
                              [[ 12,-17, -5],[  7,-17,-17]],
                              [[ -5, 12,-17],[-17,  7,-17]],
                              [[ 17,  5, 12],[ 17, 17,  7]],
                              [[-12,-17,  5],[ -7,-17, 17]],
                              [[  5,-12,-17],[ 17, -7,-17]],
                              [[-17, -5, 12],[-17,-17,  7]],
                              [[ 12, 17,  5],[  7, 17, 17]],
                              [[ -5, 17,-12],[-17, 17, -7]],
                              [[-12, -5, 17],[ -7,-17, 17]],
                              [[ 17,-12, -5],[ 17, -7,-17]],
                              [[  5,-17,-12],[ 17,-17, -7]],
                              [[ 12,  5, 17],[  7, 17, 17]],
                              [[-17, 12, -5],[-17,  7,-17]],
                              [[ -5,-17, 12],[-17,-17,  7]],
                              [[-12,  5,-17],[ -7, 17,-17]],
                              [[ 17, 12,  5],[ 17,  7, 17]],
                              [[  5, 17, 12],[ 17, 17,  7]],
                              [[ 12, -5,-17],[  7,-17,-17]],
                              [[-17,-12,  5],[-17,  7, 17]]]),
              'GTdash': \
                    np.array([[[  0,  1, -1],[  1,  1, -1]],
                              [[ -1,  0,  1],[ -1,  1,  1]],
                              [[  1, -1,  0],[  1, -1,  1]],
                              [[  0, -1, -1],[ -1, -1, -1]],
                              [[  1,  0,  1],[  1, -1,  1]],
                              [[  1, -1,  0],[  1, -1, -1]],
                              [[  0,  1, -1],[ -1,  1, -1]],
                              [[  1,  0,  1],[  1,  1,  1]],
                              [[ -1, -1,  0],[ -1, -1,  1]],
                              [[  0, -1, -1],[  1, -1, -1]],
                              [[ -1,  0,  1],[ -1, -1,  1]],
                              [[ -1, -1,  0],[ -1, -1, -1]],
                              [[  0, -1,  1],[  1, -1,  1]],
                              [[  1,  0, -1],[  1,  1, -1]],
                              [[ -1,  1,  0],[ -1,  1,  1]],
                              [[  0,  1,  1],[ -1,  1,  1]],
                              [[ -1,  0, -1],[ -1, -1, -1]],
                              [[ -1,  1,  0],[ -1,  1, -1]],
                              [[  0, -1,  1],[ -1, -1,  1]],
                              [[ -1,  0, -1],[ -1,  1, -1]],
                              [[  1,  1,  0],[  1,  1,  1]],
                              [[  0,  1,  1],[  1,  1,  1]],
                              [[  1,  0, -1],[  1, -1, -1]],
                              [[  1,  1,  0],[  1,  1, -1]]]),
              'NW': \
                    np.array([[[  2, -1, -1],[  0, -1,  1]],
                              [[ -1,  2, -1],[  0, -1,  1]],
                              [[ -1, -1,  2],[  0, -1,  1]],
                              [[ -2, -1, -1],[  0, -1,  1]],
                              [[  1,  2, -1],[  0, -1,  1]],
                              [[  1, -1,  2],[  0, -1,  1]],
                              [[  2,  1, -1],[  0, -1,  1]],
                              [[ -1, -2, -1],[  0, -1,  1]],
                              [[ -1,  1,  2],[  0, -1,  1]],
                              [[ -1,  2,  1],[  0, -1,  1]],
                              [[ -1,  2,  1],[  0, -1,  1]],
                              [[ -1, -1, -2],[  0, -1,  1]]]),
              'Pitsch': \
                    np.array([[[  1,  0,  1],[  1, -1,  1]],
                              [[  1,  1,  0],[  1,  1, -1]],
                              [[  0,  1,  1],[ -1,  1,  1]],
                              [[  0,  1, -1],[ -1,  1, -1]],
                              [[ -1,  0,  1],[ -1, -1,  1]],
                              [[  1, -1,  0],[  1, -1, -1]],
                              [[  1,  0, -1],[  1, -1, -1]],
                              [[ -1,  1,  0],[ -1,  1, -1]],
                              [[  0, -1,  1],[ -1, -1,  1]],
                              [[  0,  1,  1],[ -1,  1,  1]],
                              [[  1,  0,  1],[  1, -1,  1]],
                              [[  1,  1,  0],[  1,  1, -1]]]),
              'Bain': \
                    np.array([[[  0,  1,  0],[  0,  1,  1]],
                              [[  0,  0,  1],[  1,  0,  1]],
                              [[  1,  0,  0],[  1,  1,  0]]]),
              }
    myPlane   = [float(i) for i in planes[relationModel][variant,me]]                               # map(float, planes[...]) does not work in python 3
    myPlane  /= np.linalg.norm(myPlane)
    myNormal  = [float(i) for i in normals[relationModel][variant,me]]                              # map(float, planes[...]) does not work in python 3
    myNormal /= np.linalg.norm(myNormal)
    myMatrix  = np.array([myPlane,myNormal,np.cross(myPlane,myNormal)])
    otherPlane   = [float(i) for i in planes[relationModel][variant,other]]                         # map(float, planes[...]) does not work in python 3
    otherPlane  /= np.linalg.norm(otherPlane)
    otherNormal  = [float(i) for i in normals[relationModel][variant,other]]                        # map(float, planes[...]) does not work in python 3
    otherNormal /= np.linalg.norm(otherNormal)
    otherMatrix  = np.array([otherPlane,otherNormal,np.cross(otherPlane,otherNormal)])
    rot=np.dot(otherMatrix.T,myMatrix)
    return Orientation(matrix=np.dot(rot,self.asMatrix()))                                      # no symmetry information ??

  #@}
