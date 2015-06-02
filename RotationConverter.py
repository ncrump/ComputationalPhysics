# Created on Fri Oct 25 13:25:30 2013
"""
For Arm Performance Testing
-------------------------
Converts from Rx,Ry,Rz angles to rotation matrix using X-Y-Z Fixed Angle convention
Converts back from X-Y-Z Fixed Angle rotation matrix to Rx,Ry,Rz angles
"""

import numpy as np
from math import sin,cos,atan2,pi


# --------------------------------------------------------------
# R = rotation matrix entered as 3x3 array
# --------------------------------------------------------------
def MatrixToAngles(R):

    ry = atan2(-R[2][0],(R[0][0]**2 + R[1][0]**2)**0.5)
    rz = atan2(R[1][0]/cos(ry),R[0][0]/cos(ry))
    rx = atan2(R[2][1]/cos(ry),R[2][2]/cos(ry))
    
    # return X-Y-Z fixed angles rotated about x,then y,then z in radians
    return rx,ry,rz
# --------------------------------------------------------------


# --------------------------------------------------------------
# rx = angle about x (degrees)
# ry = angle about y (degrees)
# rz = angle about z (degrees)
# --------------------------------------------------------------
def AnglesToMatrix(rx,ry,rz):

    rx = rx*(pi/180)
    ry = ry*(pi/180)
    rz = rz*(pi/180)
    
    Rx = np.matrix(([1,0,0],[0,cos(rx),-sin(rx)],[0,sin(rx),cos(rx)]))
    Ry = np.matrix(([cos(ry),0,sin(ry)],[0,1,0],[-sin(ry),0,cos(ry)]))
    Rz = np.matrix(([cos(rz),-sin(rz),0],[sin(rz),cos(rz),0],[0,0,1]))

    R = np.round(Rz*Ry*Rx,12)
        
    # return compound rotation matrix in X-Y-Z Fixed Angle convention
    return R
# --------------------------------------------------------------


# test angles
R = AnglesToMatrix(30,60,90)
print R, '\n'

rx, ry, rz = MatrixToAngles(R)
print rx*(180/pi), ry*(180/pi), rz*(180/pi)