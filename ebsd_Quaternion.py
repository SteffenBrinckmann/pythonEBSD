# -*- coding: utf-8 -*-
##
# @file
# @brief Quaternion = mathematical rotation description
#
# Part of DAMASK (http:\\damask.mpie.de); Martin Diehl, Philip Eisenlohr, Franz Roters
# vector versions by Steffen Brinckmann
import numpy as np
import math,random,os

##
#  Quaternion: mathematical rotation description
############################## <HR>
# All methods and naming conventions based of
# http://www.euclideanspace.com/maths/algebra/realNormedAlgebra/quaternions <br>
# code derived from http://pyeuclid.googlecode.com/svn/trunk/euclid.py <br>
# suggested reading: http://web.mit.edu/2.998/www/QuaternionReport1.pdf <br>
#
# w is the real part, (x, y, z) are the imaginary parts
#
# Representation of rotation is in ACTIVE form!<br>
# (derived directly or through angleAxis, Euler angles, or active matrix)<br>
# vector "a" (defined in coordinate system "A") is actively rotated to new coordinates "b"<br>
# b = Q * a<br>
# b = np.dot(Q.asMatrix(),a)
#
class Quaternion:

    # @name CONVENTIONAL ROUTINES
    #@{

    def __init__(self,
                 quatArray = [1.0,0.0,0.0,0.0]):
      self.w, self.x, self.y, self.z = quatArray
      self.homomorph()


    def __iter__(self):
      return iter([self.w,self.x,self.y,self.z])


    def __copy__(self):
      Q = Quaternion([self.w,self.x,self.y,self.z])
      return Q
    copy = __copy__


    def __repr__(self):
      if "float" in type(self.w).__name__:
        return 'Q(r=%+.6f, i=<%+.6f, %+.6f, %+.6f>)' % (self.w, self.x, self.y, self.z)
      else:
	      return 'Quaternion array of length '+str(len(self.w))+str(type(self.w))


    def __ne__(self,other):
      return not self.__eq__(self,other)


    def __cmp__(self,other):
      return cmp(self.Rodrigues(),other.Rodrigues())


    def __eq__(self,other):
      return (abs(self.w-other.w) < 1e-8 and \
              abs(self.x-other.x) < 1e-8 and \
              abs(self.y-other.y) < 1e-8 and \
              abs(self.z-other.z) < 1e-8) \
              or \
             (abs(-self.w-other.w) < 1e-8 and \
              abs(-self.x-other.x) < 1e-8 and \
              abs(-self.y-other.y) < 1e-8 and \
              abs(-self.z-other.z) < 1e-8)

    def __getitem__(self,item):
      return Quaternion( [self.w[item],
                          self.x[item],
                          self.y[item],
                          self.z[item]]  )


    #@}
    ##
    # @name MATHEMATICAL SPECIFIC ROUTINES
    #@{

    def __pow__(self, exponent):
      omega = math.acos(self.w)
      vRescale = math.sin(exponent*omega)/math.sin(omega)
      Q = Quaternion()
      Q.w = math.cos(exponent*omega)
      Q.x = self.x * vRescale
      Q.y = self.y * vRescale
      Q.z = self.z * vRescale
      return Q


    def __ipow__(self, exponent):
      omega    = math.acos(self.w)
      vRescale = math.sin(exponent*omega)/math.sin(omega)
      self.w  = np.cos(exponent*omega)
      self.x *= vRescale
      self.y *= vRescale
      self.z *= vRescale
      return self


    def __mul__(self, other):
      try:                                                          # quaternion
          Aw = self.w
          Ax = self.x
          Ay = self.y
          Az = self.z
          Bw = other.w
          Bx = other.x
          By = other.y
          Bz = other.z
          Q = Quaternion()
          Q.w = - Ax * Bx - Ay * By - Az * Bz + Aw * Bw
          Q.x = + Ax * Bw + Ay * Bz - Az * By + Aw * Bx
          Q.y = - Ax * Bz + Ay * Bw + Az * Bx + Aw * By
          Q.z = + Ax * By - Ay * Bx + Az * Bw + Aw * Bz
          return Q
      except: pass
      try:                                                         # vector (perform active rotation, i.e. q*v*q.conjugated)
          w = self.w
          x = self.x
          y = self.y
          z = self.z
          Vx = other[0]
          Vy = other[1]
          Vz = other[2]
          return np.array([\
             w * w * Vx + 2 * y * w * Vz - 2 * z * w * Vy + \
             x * x * Vx + 2 * y * x * Vy + 2 * z * x * Vz - \
             z * z * Vx - y * y * Vx,
             2 * x * y * Vx + y * y * Vy + 2 * z * y * Vz + \
             2 * w * z * Vx - z * z * Vy + w * w * Vy - \
             2 * x * w * Vz - x * x * Vy,
             2 * x * z * Vx + 2 * y * z * Vy + \
             z * z * Vz - 2 * w * y * Vx - y * y * Vz + \
             2 * w * x * Vy - x * x * Vz + w * w * Vz ])
      except: pass
      try:                                                        # scalar
          Q = self.copy()
          Q.w *= other
          Q.x *= other
          Q.y *= other
          Q.z *= other
          return Q
      except:
          return self.copy()


    def __imul__(self, other):
      try:                                                        # Quaternion
          Ax = self.x
          Ay = self.y
          Az = self.z
          Aw = self.w
          Bx = other.x
          By = other.y
          Bz = other.z
          Bw = other.w
          self.x =  Ax * Bw + Ay * Bz - Az * By + Aw * Bx
          self.y = -Ax * Bz + Ay * Bw + Az * Bx + Aw * By
          self.z =  Ax * By - Ay * Bx + Az * Bw + Aw * Bz
          self.w = -Ax * Bx - Ay * By - Az * Bz + Aw * Bw
      except: pass
      return self


    def __div__(self, other):
      if isinstance(other, (int,float)):
        w = self.w / other
        x = self.x / other
        y = self.y / other
        z = self.z / other
        return self.__class__([w,x,y,z])
      else:
          return NotImplemented


    def __idiv__(self, other):
      if isinstance(other, (int,float)):
          self.w /= other
          self.x /= other
          self.y /= other
          self.z /= other
      return self


    def __add__(self, other):
      if isinstance(other, Quaternion):
        w = self.w + other.w
        x = self.x + other.x
        y = self.y + other.y
        z = self.z + other.z
        return self.__class__([w,x,y,z])
      else:
          return NotImplemented


    def __iadd__(self, other):
      if isinstance(other, Quaternion):
          self.w += other.w
          self.x += other.x
          self.y += other.y
          self.z += other.z
      return self


    def __sub__(self, other):
      if isinstance(other, Quaternion):
          Q = self.copy()
          Q.w -= other.w
          Q.x -= other.x
          Q.y -= other.y
          Q.z -= other.z
          return Q
      else:
          return self.copy()


    def __isub__(self, other):
      if isinstance(other, Quaternion):
          self.w -= other.w
          self.x -= other.x
          self.y -= other.y
          self.z -= other.z
      return self


    def __neg__(self):
      self.w = -self.w
      self.x = -self.x
      self.y = -self.y
      self.z = -self.z
      return self


    def __abs__(self):
      return np.sqrt(self.w ** 2 + \
                     self.x ** 2 + \
                     self.y ** 2 + \
                     self.z ** 2)
    magnitude = __abs__


    def magnitude_squared(self):
      return self.w ** 2 + \
             self.x ** 2 + \
             self.y ** 2 + \
             self.z ** 2


    def identity(self):
      self.w = 1.
      self.x = 0.
      self.y = 0.
      self.z = 0.
      return self


    def normalize(self):
      d = self.magnitude()
      if "float" in type(self.w).__name__:
        if d > 0.0:
            self /= d
      else:
        mask = d > 0.0
        self.w[mask] /= d[mask]
        self.x[mask] /= d[mask]
        self.y[mask] /= d[mask]
        self.z[mask] /= d[mask]
      return self


    def conjugate(self):
      self.x = -self.x
      self.y = -self.y
      self.z = -self.z
      return self


    def inverse(self):
      d = self.magnitude()
      if d > 0.0:
        self.conjugate()
        self /= d
      return self


    def homomorph(self):
      if "float" in type(self.w).__name__:
        if self.w < 0.0:
          self.w = -self.w
          self.x = -self.x
          self.y = -self.y
          self.z = -self.z
      else:
        mask = self.w<0.0
        self.w[mask] = -self.w[mask]
        self.x[mask] = -self.x[mask]
        self.y[mask] = -self.y[mask]
        self.z[mask] = -self.z[mask]
      return self


    #@}
    ##
    # @name GET ROTATION IN DIFFERENT FORMS
    #@{

    def normalized(self):
      return self.copy().normalize()


    def conjugated(self):
      return self.copy().conjugate()


    def inversed(self):
      return self.copy().inverse()


    def homomorphed(self):
      return self.copy().homomorph()


    def asList(self):
      return [i for i in self]


    def asM(self):                                                # to find Averaging Quaternions (see F. Landis Markley et al.)
      return np.outer([i for i in self],[i for i in self])


    def asMatrix(self):
      return np.array([[1.0-2.0*(self.y*self.y+self.z*self.z),     2.0*(self.x*self.y-self.z*self.w),     2.0*(self.x*self.z+self.y*self.w)],
                       [    2.0*(self.x*self.y+self.z*self.w), 1.0-2.0*(self.x*self.x+self.z*self.z),     2.0*(self.y*self.z-self.x*self.w)],
                       [    2.0*(self.x*self.z-self.y*self.w),     2.0*(self.x*self.w+self.y*self.z), 1.0-2.0*(self.x*self.x+self.y*self.y)]])


    def asAngleAxis(self, degrees = False):
      """
      Returns: tuple(angle,axis)
      """
      if "float" in type(self.w).__name__:
        if self.w > 1:
            self.normalize()
        s = math.sqrt(1. - self.w**2)
        x = 2*self.w**2 - 1.
        y = 2*self.w * s
        angle = math.atan2(y,x)
        if angle < 0.0:
          angle *= -1.
          s     *= -1.
        return (np.degrees(angle) if degrees else angle,
                np.array([1.0, 0.0, 0.0] if np.abs(angle) < 1e-6 else [self.x / s, self.y / s, self.z / s]))
      else:
        toNorm= np.abs(self.w)>1.0
        w,x,y,z = self.w.copy()[~toNorm],self.x.copy()[~toNorm],self.y.copy()[~toNorm],self.z.copy()[~toNorm]
        self.normalized()
        self.w[~toNorm], self.x[~toNorm], self.y[~toNorm], self.z[~toNorm] = w,x,y,z
        s = np.sqrt(1. - self.w**2)
        x = 2*self.w**2 - 1.
        y = 2*self.w * s
        angle = np.arctan2(y,x)
        mask = angle<0.0
        s[angle<0.0] *= -1.
        angle[angle<0.0] *= -1.
        axis = np.vstack((self.x,self.y,self.z))
        mask = np.abs(angle)<1e-6
        axis[0,mask]  = 1.
        axis[1:,mask] = 0.
        axis[:,~mask]/= s[~mask]
        return (np.degrees(angle) if degrees else angle, axis)


    def asRodrigues(self):
      if "float" in type(self.w).__name__:
        return np.inf*np.ones(3) if self.w == 0.0 else np.array([self.x, self.y, self.z])/self.w
      else:
        result = np.vstack((self.x,self.y,self.z))
        mask   = np.abs(self.w)<1e-6
        result[:,mask]   = np.inf
        result[:,~mask] /= self.w[~mask]
        return result


    def asEulers(self,
                 notation = 'bunge',
                 degrees = False,
                 standardRange = False,
                 round=8):
      '''
      conversion of ACTIVE rotation to Euler angles taken from.

      Melcher, A.; Unser, A.; Reichhardt, M.; Nestler, B.; Pötschke, M.; Selzer, M.
      Conversion of EBSD data by a quaternion based algorithm to be used for grain structure simulations
      Technische Mechanik 30 (2010) pp 401--413
      '''
      if notation.lower() != 'bunge' and notation.lower() != 'zxz':
	      angles = [0.0,0.0,0.0]
	      return np.degrees(angles) if degrees else angles
      if "float" in type(self.w).__name__:  #scalar version
        angles = [0.0,0.0,0.0]
        if   abs(self.x) < 1e-4 and abs(self.y) < 1e-4:
          x = self.w**2 - self.z**2
          y = 2.*self.w*self.z
          angles[0] = math.atan2(y,x)
        elif abs(self.w) < 1e-4 and abs(self.z) < 1e-4:
          x = self.x**2 - self.y**2
          y = 2.*self.x*self.y
          angles[0] = math.atan2(y,x)
          angles[1] = math.pi
        else:
          chi = math.sqrt((self.w**2 + self.z**2)*(self.x**2 + self.y**2))
          x = (self.w * self.x - self.y * self.z)/2./chi
          y = (self.w * self.y + self.x * self.z)/2./chi
          angles[0] = math.atan2(y,x)
          x = self.w**2 + self.z**2 - (self.x**2 + self.y**2)
          y = 2.*chi
          angles[1] = math.atan2(y,x)
          x = (self.w * self.x + self.y * self.z)/2./chi
          y = (self.z * self.x - self.y * self.w)/2./chi
          angles[2] = math.atan2(y,x)
        if standardRange:
          angles[0] %= 2*math.pi
          if angles[1] < 0.0:
            angles[1] += math.pi
            angles[2] *= -1.0
          angles[2] %= 2*math.pi
        return np.round(np.degrees(angles),round) if degrees else angles
      else:  						#vector version
        angles = np.zeros( (3,len(self.w) ))
        mask1 = np.logical_and( np.abs(self.x)<1e-4, np.abs(self.y)<1e-4 )
        x = self.w[mask1]**2 - self.z[mask1]**2
        y = 2.*self.w[mask1]*self.z[mask1]
        angles[0,mask1] = np.arctan2(y,x)
        mask2 = np.logical_and( np.abs(self.w) < 1e-4, abs(self.z) < 1e-4)
        mask2 = np.logical_and( ~mask1, mask2)
        x = self.x[mask2]**2 - self.y[mask2]**2
        y = 2.*self.x[mask2]*self.y[mask2]
        angles[0,mask2] = np.arctan2(y,x)
        angles[1,mask2] = math.pi
        mask3 = np.logical_and(~mask1, ~mask2)
        chi = np.sqrt((self.w[mask3]**2 + self.z[mask3]**2)*(self.x[mask3]**2 + self.y[mask3]**2))
        x = (self.w[mask3] * self.x[mask3] - self.y[mask3] * self.z[mask3])/2./chi
        y = (self.w[mask3] * self.y[mask3] + self.x[mask3] * self.z[mask3])/2./chi
        angles[0,mask3] = np.arctan2(y,x)
        x = self.w[mask3]**2 + self.z[mask3]**2 - (self.x[mask3]**2 + self.y[mask3]**2)
        y = 2.*chi
        angles[1,mask3] = np.arctan2(y,x)
        x = (self.w[mask3] * self.x[mask3] + self.y[mask3] * self.z[mask3])/2./chi
        y = (self.z[mask3] * self.x[mask3] - self.y[mask3] * self.w[mask3])/2./chi
        angles[2,mask3] = np.arctan2(y,x)
        if standardRange:
          angles[0,:] %= 2*math.pi
          mask_ = angles[1,:]<0.0
          angles[1,mask_] += math.pi
          angles[2,mask_] *= -1.0
          angles[2,:] %= 2*math.pi
        return np.round(np.degrees(angles),round) if degrees else angles





    #@}
    ##
    # @name STATIC CLASS METHODS
    #@{

    @classmethod
    def fromIdentity(cls):
      return cls()


    @classmethod
    def fromRandom(cls,randomSeed = None):
      if randomSeed == None:
        randomSeed = int(os.urandom(4).encode('hex'), 16)
      random.seed(randomSeed)
      r1 = random.random()
      r2 = random.random()
      r3 = random.random()
      w = math.cos(2.0*math.pi*r1)*math.sqrt(r3)
      x = math.sin(2.0*math.pi*r2)*math.sqrt(1.0-r3)
      y = math.cos(2.0*math.pi*r2)*math.sqrt(1.0-r3)
      z = math.sin(2.0*math.pi*r1)*math.sqrt(r3)
      return cls([w,x,y,z])


    @classmethod
    def fromRodrigues(cls, rodrigues):
      if not isinstance(rodrigues, np.ndarray): rodrigues = np.array(rodrigues)
      halfangle = math.atan(np.linalg.norm(rodrigues))
      c = math.cos(halfangle)
      w = c
      x,y,z = c*rodrigues
      return cls([w,x,y,z])


    @classmethod
    def fromAngleAxis(cls, angle, axis):
      if not isinstance(axis, np.ndarray): axis = np.array(axis,dtype='d')
      axis = axis.astype(float)/np.linalg.norm(axis)
      s = math.sin(angle / 2.0)
      w = math.cos(angle / 2.0)
      x = axis[0] * s
      y = axis[1] * s
      z = axis[2] * s
      return cls([w,x,y,z])


    @classmethod
    def fromEulers(cls, eulers, type = 'Bunge'):
      eulers *= 0.5                       # reduce to half angles
      if eulers.ndim==1:
        c1 = math.cos(eulers[0])
        s1 = math.sin(eulers[0])
        c2 = math.cos(eulers[1])
        s2 = math.sin(eulers[1])
        c3 = math.cos(eulers[2])
        s3 = math.sin(eulers[2])
      else:
        c1 = np.cos(eulers[0,:])
        s1 = np.sin(eulers[0,:])
        c2 = np.cos(eulers[1,:])
        s2 = np.sin(eulers[1,:])
        c3 = np.cos(eulers[2,:])
        s3 = np.sin(eulers[2,:])
      if type.lower() == 'bunge' or type.lower() == 'zxz':
        w =   c1 * c2 * c3 - s1 * c2 * s3
        x =   c1 * s2 * c3 + s1 * s2 * s3
        y = - c1 * s2 * s3 + s1 * s2 * c3
        z =   c1 * c2 * s3 + s1 * c2 * c3
      else:
#          print 'unknown Euler convention'
        w = c1 * c2 * c3 - s1 * s2 * s3
        x = s1 * s2 * c3 + c1 * c2 * s3
        y = s1 * c2 * c3 + c1 * s2 * s3
        z = c1 * s2 * c3 - s1 * c2 * s3
      return cls([w,x,y,z])


    @classmethod
    def fromMatrix(cls, m):
      """
      Modified Method to calculate Quaternion from Orientation Matrix

      Source: http://www.euclideanspace.com/maths/geometry/rotations/conversions/matrixToQuaternion/
      """
      if m.shape != (3,3) and np.prod(m.shape) == 9:
        m = m.reshape(3,3)
      tr = np.trace(m)
      if tr > 1e-8:
        s = math.sqrt(tr + 1.0)*2.0
        return cls(
          [ s*0.25,
            (m[2,1] - m[1,2])/s,
            (m[0,2] - m[2,0])/s,
            (m[1,0] - m[0,1])/s,
          ])
      elif m[0,0] > m[1,1] and m[0,0] > m[2,2]:
        t = m[0,0] - m[1,1] - m[2,2] + 1.0
        s = 2.0*math.sqrt(t)
        return cls(
          [ (m[2,1] - m[1,2])/s,
            s*0.25,
            (m[0,1] + m[1,0])/s,
            (m[2,0] + m[0,2])/s,
          ])
      elif m[1,1] > m[2,2]:
        t = -m[0,0] + m[1,1] - m[2,2] + 1.0
        s = 2.0*math.sqrt(t)
        return cls(
          [ (m[0,2] - m[2,0])/s,
            (m[0,1] + m[1,0])/s,
            s*0.25,
            (m[1,2] + m[2,1])/s,
          ])
      else:
        t = -m[0,0] - m[1,1] + m[2,2] + 1.0
        s = 2.0*math.sqrt(t)
        return cls(
          [ (m[1,0] - m[0,1])/s,
            (m[2,0] + m[0,2])/s,
            (m[1,2] + m[2,1])/s,
            s*0.25,
          ])


    @classmethod
    def new_interpolate(cls, q1, q2, t):
        """
        Interpolation

        see http://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/20070017872_2007014421.pdf for (another?) way to interpolate quaternions
        """
        assert isinstance(q1, Quaternion) and isinstance(q2, Quaternion)
        Q = cls()
        costheta = q1.w * q2.w + q1.x * q2.x + q1.y * q2.y + q1.z * q2.z
        if costheta < 0.:
            costheta = -costheta
            q1 = q1.conjugated()
        elif costheta > 1:
            costheta = 1
        theta = math.acos(costheta)
        if abs(theta) < 0.01:
            Q.w = q2.w
            Q.x = q2.x
            Q.y = q2.y
            Q.z = q2.z
            return Q
        sintheta = math.sqrt(1.0 - costheta * costheta)
        if abs(sintheta) < 0.01:
            Q.w = (q1.w + q2.w) * 0.5
            Q.x = (q1.x + q2.x) * 0.5
            Q.y = (q1.y + q2.y) * 0.5
            Q.z = (q1.z + q2.z) * 0.5
            return Q
        ratio1 = math.sin((1 - t) * theta) / sintheta
        ratio2 = math.sin(t * theta) / sintheta
        Q.w = q1.w * ratio1 + q2.w * ratio2
        Q.x = q1.x * ratio1 + q2.x * ratio2
        Q.y = q1.y * ratio1 + q2.y * ratio2
        Q.z = q1.z * ratio1 + q2.z * ratio2
        return Q

