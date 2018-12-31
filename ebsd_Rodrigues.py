# -*- coding: utf-8 -*-
##
# @file
# @brief Rodrigues representation of rotation
#
# Part of DAMASK (http:\\damask.mpie.de); Martin Diehl, Philip Eisenlohr, Franz Roters
#
import numpy as np
from ebsd_Quaternion import Quaternion

##
#  Rodrigues representation of rotation
############################## <HR>
#
class Rodrigues:

    def __init__(self, vector = np.zeros(3)):
      self.vector = vector


    def asQuaternion(self):
      norm = np.linalg.norm(self.vector)
      halfAngle = np.arctan(norm)
      return Quaternion(np.cos(halfAngle),np.sin(halfAngle)*self.vector/norm)


    def asAngleAxis(self):
      norm = np.linalg.norm(self.vector)
      halfAngle = np.arctan(norm)
      return (2.0*halfAngle,self.vector/norm)
