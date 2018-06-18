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

  __slots__ = ['quaternion','symmetry']

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
    elif isinstance(matrix, np.ndarray) :                                                           # based on given rotation matrix
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
    self.doctest  = False


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
    return map(lambda q: Orientation(quaternion = q, symmetry = self.symmetry.lattice),
               self.equivalentQuaternions(who))


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
      disorientation angle
    """

    if self.symmetry != other.symmetry: raise TypeError('disorientation between different symmetry classes not supported yet.')
    misQ = self.quaternion.conjugated()*other.quaternion
    mySymQs    =  self.symmetry.symmetryQuats() if SST else self.symmetry.symmetryQuats()[:1]       # take all or only first sym operation
    otherSymQs = other.symmetry.symmetryQuats()
    for i,sA in enumerate(mySymQs):
      for j,sB in enumerate(otherSymQs):
        theQ = sA.conjugated()*misQ*sB
        for k in xrange(2):
          theQ.conjugate()
          breaker = self.symmetry.inFZ(theQ) \
                    and (not SST or other.symmetry.inDisorientationSST(theQ))
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
  def plot(self, axis=None, scale=2, proj2D=None, show=True):
    """Plot rotated unit-cell in 3D, and possibly the pole-figure and specific poles

    Projection onto 2D: cooradinate systems are given as xDirection-yDirection (z follows)
    - down-right: [default in text books] RD = x = down; TD = y = right; ND = z = outOfPlane
    - up-left: [default in OIM] RD = x = up; TD = y = left; ND = z = outOfPlane

    Args:
       axis: if given (e.g. [1,0,0]), plot pole-figure and the corresponding poles
       scale: scale of pole-figure dome over crystal
       proj2D: do a normal projection onto 2D plane: [down-right, up-left, None]
    """
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    eps = 1e-6
    if proj2D is not None:
      fig, ax = plt.subplots()
    else:
      fig = plt.figure()
      ax = fig.gca(projection='3d')
      #ax.view_init(90,0)

    def plotLine(start,delta,color='k',lw=1):
      if np.linalg.norm(delta)<eps:
        marker='o';  markerSize=7
      else:
        marker=None; markerSize=0
      if proj2D=='down-right':
        ax.plot( [start[1]]+[start[1]+delta[1]],
                 [-start[0]]+[-start[0]-delta[0]],
                 color=color,lw=lw, marker=marker, markersize=markerSize)
      elif proj2D=='up-left':
        ax.plot( [-start[1]]+[-start[1]-delta[1]],
                 [start[0]]+[start[0]+delta[0]],
                 color=color,lw=lw, marker=marker, markersize=markerSize)
      else:
        ax.plot( [start[0]]+[start[0]+delta[0]],
                 [start[1]]+[start[1]+delta[1]],
                 [start[2]]+[start[2]+delta[2]],
                 color=color,lw=lw, marker=marker, markersize=markerSize)

    for line in self.symmetry.unitCell():
      start, end = np.array(line[:3],dtype=np.float), np.array(line[3:],dtype=np.float)
      start = self.quaternion*start
      end   = self.quaternion*end
      if start[2]<0 and end[2]<0:
        plotLine(start, end-start,color="b",lw=0.2)
      elif start[2]>0 and end[2]>0:
        plotLine(start, end-start,color="b",lw=2)
      else:
        delta = end-start
        k     = -start[2]/delta[2]
        mid   = start+k*delta
        if start[2]>0:
          plotLine(start, mid-start,color="b",lw=2)
          plotLine(mid,   end-mid,color="b",lw=0.2)
        else:
          plotLine(start, mid-start,color="b",lw=0.2)
          plotLine(mid,   end-mid,color="b",lw=2)
    plotLine([0,0,0], [1,0,0], 'k',lw=3)
    plotLine([0,0,0], [0,1,0], 'k',lw=3)
    plotLine([0,0,0], [0,0,1], 'k',lw=3)
    if proj2D=='down-right':
      ax.text(0.1,-1,"RD [100]")
      ax.text(1,-0.1,"TD [010]")
      ax.text(0.1,-0.1,"ND [001]")
    elif proj2D=='up-left':
      ax.text(-0.1,1,"RD [100]")
      ax.text(-1,0.1,"TD [010]")
      ax.text(-0.1,0.1,"ND [001]")
    else:
      ax.text(1,0.1,0.1,"RD [100]")
      ax.text(0.1,1,0.1,"TD [010]")
      ax.text(0.1,0.1,1,"ND [001]")

    if axis is not None:
      #plot sphere
      u = np.linspace(0,2*np.pi,50)
      v = np.linspace(0,np.pi/2,50)
      x = scale*np.outer(np.cos(u)      , np.sin(v))
      y = scale*np.outer(np.sin(u)      , np.sin(v))
      z = scale*np.outer(np.ones_like(u), np.cos(v))
      if proj2D is None:
        ax.plot_surface(x, y, z, color='gray',alpha=0.7, cstride=10, rstride=10, lw=0)
      else:
        ax.plot( scale*np.cos(np.linspace(0.,2.*np.pi,100)), scale*np.sin(np.linspace(0.,2.*np.pi,100)), 'k--' )
      # plot poles
      oHelp = Orientation(Eulers=np.array([0.,0.,0.]), symmetry=self.symmetry.__repr__())
      axis = np.array(axis, dtype=np.float)
      axis /= np.linalg.norm(axis)
      for q in oHelp.symmetry.equivalentQuaternions(oHelp.quaternion):
        conjAxis      = q*axis
        direction = self.quaternion*conjAxis
        if direction[2]<-eps: continue                                                        #prevent rounding errors
        fromBase = direction+np.array([0,0,1])
        plotLine([0,0,0], direction*scale, color='b', lw=1)
        plotLine(-scale*np.array([0,0,1]), (fromBase)*scale, 'skyblue')
        xy  = fromBase/fromBase[2]*scale
        xy[2]=0.0
        plotLine( xy, [0.,0.,0.], color='skyblue')
    ax.axis('equal'), ax.axis('off')
    ax.set_xlabel(''); ax.set_xticks([])
    ax.set_ylabel(''); ax.set_yticks([])
    if not proj2D:   ax.set_zlabel(''); ax.set_zticks([])
    if self.doctest: plt.savefig('doctest.png'); plt.close(); return
    if show is True:
      plt.show()
    return


  def toScreen(self, equivalent=True):
    """
    print Euler angles and HKL /UVW

    Args:
      equivalent: print also equivalent orientations
    """
    print "Euler angles:",np.round(self.quaternion.asEulers(degrees=True),1)
    rotM = self.quaternion.asMatrix()
    print "HKL",np.array( rotM[2,:]/np.min(rotM[2,:]) ,dtype=np.int)
    print "UVW",-np.array( rotM[0,:]/np.min(rotM[0,:]) ,dtype=np.int)
    if equivalent:
      print "Equivalent orientations - Euler angles:"
      for q in self.symmetry.equivalentQuaternions(self.quaternion):
        angles = q.asEulers(degrees=True)
        angles[angles<0] += 360.
        print "   [%5.1f  %5.1f  %5.1f]"%tuple(angles)
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
