class depricatedMethods:


  #@}
  ##
  # @name GET - SAVE METHODS
  #@{

  def getAngles(self, x,y):
    """
    get Euler angles at (vector of) coodinates

    Args:
       x: (vector of) x-coordinates
       y: (vector of) y-coordinates

    Returns:
       euler angles
    """
    startTime = time.time()
    from matplotlib.tri import Triangulation, LinearTriInterpolator, CubicTriInterpolator
    triang = Triangulation(self.x, self.y)
    triangInterpolator_phi1 = LinearTriInterpolator(triang, self.phi1)
    phi1 = triangInterpolator_phi1.__call__( x, y).filled()
    triangInterpolator_PHI = LinearTriInterpolator(triang, self.PHI)
    PHI = triangInterpolator_PHI.__call__( x, y).filled()
    triangInterpolator_phi2 = LinearTriInterpolator(triang, self.phi2)
    phi2 = triangInterpolator_phi2.__call__( x, y).filled()
    print "Duration getAngles: ",np.round(time.time()-startTime,2),"sec"
    return [phi1, PHI, phi2]

  #@}
  ##
  # @name Math routines to convert between HKL, XY, RGB
  #@{

  @classmethod
  def rotationMatrix(self, phi1, PHI, phi2):
    """
    calculate rotation matrix by giving Bunge's Euler Angles

    Args:
       phi1: phi' rotation angle
       PHI: PHI rotation angle
       phi2: phi'' rotation angle

    Returns:
       rotation matrix [3x3]
    """
    #check length
    if not ( (len(phi1)==len(PHI)) and (len(phi1)==len(phi2))  ):
      print "ERROR: length of input does not match"
      return []
    R = np.zeros(( 3,3,len(phi1) ))
    R[0,0,:] =  np.cos(phi1)*np.cos(phi2) - np.sin(phi1)*np.sin(phi2)*np.cos(PHI)
    R[1,0,:] = -np.cos(phi1)*np.sin(phi2) - np.sin(phi1)*np.cos(phi2)*np.cos(PHI)
    R[2,0,:] =  np.sin(phi1)*np.sin(PHI)
    R[0,1,:] =  np.sin(phi1)*np.cos(phi2) + np.cos(phi1)*np.sin(phi2)*np.cos(PHI)
    R[1,1,:] = -np.sin(phi1)*np.sin(phi2) + np.cos(phi1)*np.cos(phi2)*np.cos(PHI)
    R[2,1,:] = -np.cos(phi1)*np.sin(PHI)
    R[0,2,:] =  np.sin(phi2)*np.sin(PHI)
    R[1,2,:] =  np.cos(phi2)*np.sin(PHI)
    R[2,2,:] =               np.cos(PHI)
    return R


  @classmethod
  def HKLtoRGB(self, hkl):
    """
    DEPRICATED, calculate RBG values from HKL vector; DOES WORK FOR VECTOR OF HKLs

    this method makes sure that all coefficients are positiv and that L>H>K

    Args:
       hkl (possibly vector of) HKLs

    Returns:
       corresponding vector of RGB values
    """
    print "DEPRICATED VERSION THAT USES DIFFERENT COLORING SCHEME. ONLY USE FOR TESTING"
    xy = self.HKLtoXY(hkl)
    shapeMatrix = hkl.shape
    angle = np.arctan2( xy[1,:], xy[0,:] )
    length= -np.cos(angle)+np.sqrt( np.cos(angle)*np.cos(angle)+1.0)
    #determine the color coding
    colors = np.zeros( shapeMatrix )
    colors[0,:] = 1.0-np.linalg.norm(xy, axis=0) / length
    colors[2,:] = angle/(np.pi/4.)  * ( 1.0-colors[0,:] )
    colors[1,:] = (1.0-colors[2,:]) * ( 1.0-colors[0,:] )
    colors = colors / np.max(colors, axis=0)
    return colors


  @classmethod
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
      l = 2.0/ ( inPlane[0]*inPlane[0] + inPlane[1]*inPlane[1] + 1.)
      if inPlane.ndim==1:
        hkl =np.zeros((3), dtype=np.float)
        hkl[:2]=l*inPlane
        hkl[2]=l-1.0
      else:
        hkl = np.zeros((3,inPlane.shape[1]), dtype=np.float)
        hkl[:2,:]=l*inPlane
        hkl[2,:]=l-1.0
      return hkl


  @classmethod
  def XYtoRGB(self,xy):
    """
    DEPRICATED, do not use, calculate RBG values from XY vector; DOES WORK FOR VECTOR OF XYs<br>

    Args:
       xy (possibly vector of) XYs

    Returns:
       corresponding vector of RGB values
    """
    print "DEPRICATED VERSION THAT USES DIFFERENT COLORING SCHEME. ONLY USE FOR TESTING"
    shapeMatrix = (3,xy.shape[1])
    angle = np.arctan2( xy[1,:], xy[0,:] )
    length= -np.cos(angle)+np.sqrt( np.cos(angle)*np.cos(angle)+1.0)
    #determine the color coding
    colors = np.zeros( shapeMatrix )
    colors[0,:] = 1.0-np.linalg.norm(xy, axis=0) / length
    colors[2,:] = angle/(np.pi/4.)  * ( 1.0-colors[0,:] )
    colors[1,:] = (1.0-colors[2,:]) * ( 1.0-colors[0,:] )
    colors = colors / np.max(colors, axis=0)
    return colors


  @classmethod
  def HKLtoXY(self,hkl):
    """
    calculate XY values from HKL vector; DOES WORK FOR VECTOR OF HKLs

    this method makes sure that all coefficients are positiv and that L>H>K

    Args:
       hkl (possibly vector of) HKLs

    Returns:
       corresponding vector of XY values in the unit-sphere
    """
    shapeMatrix = hkl.shape
    #prepare hkl: absolute values for all hkls
    #   and correctly sorted: l is the largest value, h the second largest
    hkl = np.abs(hkl)
    hkl = np.sort(hkl, axis=0)
    hkl[:2,:] = np.sort( hkl[:2,:], axis=0 )[::-1,:]
    #add the unit vector to hkl
    fromBase = np.zeros( shapeMatrix )
    fromBase[2] = np.ones( shapeMatrix[1] )
    fromBase = fromBase + hkl
    #in the drawing plane: length, angle
    xy = (fromBase[:,:]/fromBase[2,:])[:2,:]
    return xy



  def plotPointDensity(self, show=True):
    """
    Calculate the point density using the Gaussian kernel density

    Args:
      show: show figure after calculation; it is plotted anyhow
    """
    # Calculate the point density
    startTime = time.time()
    xy = np.vstack([self.x[~self.mask],self.y[~self.mask]])
    density = gaussian_kde(xy)(xy)
    # Sort the points by density, so that the densest points are plotted last
    idx = density.argsort()
    x, y, density = xy[0,:][idx], xy[1,:][idx], density[idx]
    # plot
    scpl = plt.scatter(x, y, c=density, s=1, edgecolor='', cmap=plt.cm.get_cmap('jet'))
    plt.xlim([np.min(self.x),  np.max(self.x)])
    plt.ylim([np.max(self.y),  np.min(self.y)])
    plt.colorbar(scpl)
    plt.xticks([])
    plt.yticks([])
    if show:
      plt.show()
    print "Duration plotPointDensity: ",np.round(time.time()-startTime,2),"sec"
    return



  #@}
  ##
  # @name Neighborhood functions, clean data: ALL DEPRICATED, use OIM better
  #@{

  def findNeighbors(self, i):
    """
    DEPRICATED. USE mtex system.
    find geometric neighbors to the point index i

    Args:
       i: index to identify neighbors

    Returns:
       list of neighbors
    """
    # fast method
    print "    DEPRICATED. USE mtex system. "
    possibleNeighbors = [ i-self.neighborLength-1, i-self.neighborLength, i-1, i+1, i+self.neighborLength, i+self.neighborLength+1]
    realNeighbors = []
    myX = self.x[i]
    myY = self.y[i]
    dist2 = self.stepSizeX*self.stepSizeX * 1.2  # ensure next neighbors
    for j in possibleNeighbors:
      if j<0:
        continue
      if j>=len(self.x):
        continue
      if (self.x[j]-myX)*(self.x[j]-myX) + (self.y[j]-myY)*(self.y[j]-myY) < dist2:
        realNeighbors.append(j)
    return realNeighbors


  def cleanPixelsCI(self, thresholdCI=0.1):
    """
    clean data: I suggest to use the OIM software instead
    - if CI is below threshold
    - take the values from the neighboring point, which has the highest CI
    - can be used iteratively

    Args:
       thresholdCI theshold CI [default=0.1]

    Returns:
       number of pixels that are changed
    """
    startTime = time.time()
    repairpoint = 0
    self.repaired = np.zeros( self.x.shape).astype(np.int)
    tempPhi1 = np.array( self.phi1 )
    tempPHI  = np.array( self.PHI )
    tempPhi2 = np.array( self.phi2 )
    tempCI   = np.array( self.CI )
    for i in range(len(self.x)):  #for all pixel
      neighbors = self.findNeighbors(i)
      if self.CI[i] < thresholdCI:
        bestNeighbor = np.argmax( self.CI[neighbors])
        if self.CI[bestNeighbor] < thresholdCI:
          continue
        tempPhi1[i] = self.phi1[neighbors][bestNeighbor]
        tempPHI[i]  = self.PHI[neighbors][bestNeighbor]
        tempPhi2[i] = self.phi2[neighbors][bestNeighbor]
        tempCI[i]   = self.CI[neighbors][bestNeighbor]
        self.repaired[i] = 1
        repairpoint += 1
    self.phi1 = tempPhi1
    self.PHI  = tempPHI
    self.phi2 = tempPhi2
    self.CI   = tempCI
    print "Duration cleanPixelsCI: ",np.round(time.time()-startTime,2),"sec"
    return repairpoint


  def cleanPixels(self, deltaRad=5.0, minDistantNN=6):
    """
    clean data: I sugget to use the OIM software instead
    - if neighbors are too far rotated from the center
    - take the values from the neighboring point, which has the highest CI
    - can be used iteratively
    - takes 3times as long as the cleanPixelsCI

    Args:
       deltaRad: theshold angle [default=5], that would define the grain interiour
       minDistantNN: number of neighbors that must be distant [default=6: one single difference in center is eliminated]

    Returns:
       number of pixels that are changed
    """
    startTime = time.time()
    deltaHKL = np.sin(np.pi* deltaRad/180.)
    repairpoint = 0
    self.repaired = np.zeros( self.x.shape).astype(np.int)
    tempPhi1 = np.array( self.phi1 )
    tempPHI  = np.array( self.PHI )
    tempPhi2 = np.array( self.phi2 )
    tempCI   = np.array( self.CI )
    for i in range(len(self.x)):  #for all pixel
      iHKL = np.zeros( (3,1) )
      iHKL[0,0] =  np.sin(self.phi2[i])*np.sin(self.PHI[i])
      iHKL[1,0] =  np.cos(self.phi2[i])*np.sin(self.PHI[i])
      iHKL[2,0] =               np.cos(self.PHI[i])
      iHKL = np.abs(iHKL)
      iHKL = np.sort(iHKL, axis=0)
      iHKL[:2,:] = np.sort( iHKL[:2,:], axis=0 )[::-1,:]
      # find neighbors and calculate its jHKL, and the difference
      neighbors = self.findNeighbors(i)
      jHKL = np.zeros( (3, len(neighbors) ) )
      jHKL[0,:] =  np.sin(self.phi2[neighbors])*np.sin(self.PHI[neighbors])
      jHKL[1,:] =  np.cos(self.phi2[neighbors])*np.sin(self.PHI[neighbors])
      jHKL[2,:] =               np.cos(self.PHI[neighbors])
      jHKL = np.abs(jHKL)
      jHKL = np.sort(jHKL, axis=0)
      jHKL[:2,:] = np.sort( jHKL[:2,:], axis=0 )[::-1,:]
      error = np.abs(jHKL - iHKL)
      #threshold (>deltaHKL) using the maximal HKL, and find number of neighbor which have higher than allowed distance
      distantNN = np.sum(np.sum( error>deltaHKL ,axis=0) > 1 )
      if distantNN>=minDistantNN:
        repairpoint += 1
        bestNeighbor = np.argmax( self.CI[neighbors])
        tempPhi1[i] = self.phi1[neighbors][bestNeighbor]
        tempPHI[i]  = self.PHI[neighbors][bestNeighbor]
        tempPhi2[i] = self.phi2[neighbors][bestNeighbor]
        tempCI[i]   = self.CI[neighbors][bestNeighbor]
        self.repaired[i] = 1
    self.phi1 = tempPhi1
    self.PHI  = tempPHI
    self.phi2 = tempPhi2
    self.CI   = tempCI
    print "Duration cleanPixels: ",np.round(time.time()-startTime,2),"sec"
    return repairpoint


  def __init__():


    if self.neighborLength<0:
      #find periodic length for determining neighbors during cleaning
      myX = self.x[1].astype(np.float16)
      myY = self.y[1].astype(np.float16)
      dist2 = (self.stepSizeX*self.stepSizeX * 1.2).astype(np.float16)  #ensure next neighbors
      oldSettings = np.seterr(over="ignore")
      for j in range(3, len(self.x)):
        if (self.x[j]-myX)*(self.x[j]-myX) + (self.y[j]-myY)*(self.y[j]-myY) < dist2:
          self.neighborLength = j-1
          break
      np.seterr(**oldSettings)
    print "   Periodic length in x-direction:",self.neighborLength
    if self.stepSizeY==None:
      self.stepSizeY=self.y[self.neighborLength+1] - self.y[0]
      print '   Step Size Y adopted:',self.stepSizeY

    #set up other things
    if self.stepSizeX != self.stepSizeY:  #probably hexagonal grid
      self.unitCell = [
        [-self.stepSizeX/2. , -1.*self.stepSizeY/3],
        [-self.stepSizeX/2. ,     self.stepSizeY/3],
        [0                  ,  2.*self.stepSizeY/3],
        [ self.stepSizeX/2. ,     self.stepSizeY/3],
        [ self.stepSizeX/2. , -1.*self.stepSizeY/3],
        [0                  , -2.*self.stepSizeY/3]  ]
    else:
      self.unitCell = [                 #cubic grid
        [ self.stepSizeX/2. , -self.stepSizeY/2],
        [ self.stepSizeX/2. ,  self.stepSizeY/2],
        [-self.stepSizeX/2. ,  self.stepSizeY/2],
        [-self.stepSizeX/2. , -self.stepSizeY/2] ]


