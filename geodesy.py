# -*- coding: utf-8 -*-
""""
# Coordinate Systems and Geodesy I - Fall Term, 2018
# Practice - 1
#       @creator_of_ell_class: Mustafa Serkan Isik
#       @rest_of_the_script  : Deniz Gökduman
#
Acronyms=======================================================================
#
# AT            ----> Average Terrestial
# IT            ----> Instantenous Terrestial
# G             ----> Geodetic
#
# LA            ----> Local Astronomic
# LG            ----> Local Geodetic
#
# xyz           ----> Cartesian   Coordinates
# cur           ----> Curvelinear Coordinates
#
# DOV           ----> Deflection on vertices(E, n)
#
# Upper-case    ---->  Geoid-based coordinate systems (AT/IT/LA/...)
# lower-case    ---->  Ellipsoid-based coordinate systems (G/LG/...)
# 
# R1, R2, R3    ----> Rotation matrices for x, y, and z axes
# P1, P2, P3    ----> Reflection matrices for x, y, and z axes
#
Examples=======================================================================
#
#   ATxyz       ----> Cartesian Coordinates on Average Terr. Coor Sys.
# |   X  |
# |   Y  |      ----> ATxyz Example
# |   Z  |AT
#
#   ATcur       ----> Curvelinear Coordinates on Average Terr. Coor Sys.
# | LAT  |
# | LON  |      ----> ATcur Example
# |  H   |AT
#
# ----------------------------------------------------------------------------- 
#
#   ITxyz       ----> Cartesian Coordinates on Intantenous Terr. Coor Sys.
# |   X* |
# |   Y* |      ----> ITxyz Example
# |   Z* |IT

#   ITcur       ----> Curvelinear Coordinates on Intantenous Terr. Coor Sys.
# | LAT* |
# | LON* |      ----> ITcur Example
# |  H*  |IT
#
# -----------------------------------------------------------------------------
#
#   Gxyz        ----> Cartesian Coordinates on Geodetic Coor Sys.
# |   x  |
# |   y  |      ----> Gxyz Example
# |   z  |G
#
#   Gcur        ----> Curvelinear Coordinates on Geodetic Coor Sys.
# | lat  |
# | lon  |      ----> Gcur Example
# |  h   |G
#
# ----------------------------------------------------------------------------- 
#
#   LAxyz       ----> Cartesian Coordinates on Local Astronomic Coor Sys.
# |   X_ |
# |   Y_ |      ----> LAxyz Example
# |   Z_ |LA
#
#   LAcur       ----> Curvelinear Coordinates on Local Astronomic Coor Sys.
# | LAT_ |
# | LON_ |      ----> LAcur Example
# |  H_  |LA
#
# -----------------------------------------------------------------------------
#   LGxyz       ----> Cartesian Coordinates on Local Geodetic Coor Sys.
# |   x_ |
# |   y_ |      ----> LGxyz Example
# |   z_ |LG
# 
#   LGcur       ----> Curvelinear Coordinates on Local Geodetic Coor Sys.
# | lat_ |
# | lon_ |      ----> LGxyz Example
# |  h_  |LG
# 
"""
# -----------------------------------------------------------------------------
from numpy import (pi, eye, cos, sin, tan, arccos as acos, arctan2 as atan, matrix)
from difflib import SequenceMatcher as seq
# -----------------------------------------------------------------------------
P1 = matrix([[-1,  0,  0],
             [ 0,  1,  0],
             [ 0,  0,  1]])

P2 = matrix([[ 1,  0,  0],
             [ 0, -1,  0],
             [ 0,  0,  1]])

P3 = matrix([[ 1,  0,  0],
             [ 0,  1,  0],
             [ 0,  0, -1]])
# -----------------------------------------------------------------------------
R1 = lambda ex : matrix([[       1,             0,               0],
                         [       0,   cos(ex.rad),     sin(ex.rad)],
                         [       0,  -sin(ex.rad),     cos(ex.rad)]])
    
R2 = lambda ey : matrix([[ cos(ey.rad),         0,    -sin(ey.rad)],
                         [           0,         1,               0],
                         [ sin(ey.rad),         0,     cos(ey.rad)]])
    
R3 = lambda ez : matrix([[ cos(ez.rad),   sin(ez.rad),           0],
                         [-sin(ez.rad),   cos(ez.rad),           0],
                         [           0,             0,           1]])

R  = lambda f, ex, ey, ez : f * eye(3) * R3(ez) * R2(ey) * R1(ex)
# -----------------------------------------------------------------------------
class Ellipsoid:
    #
    # @author: Mustafa Serkan Isik
    # 
	def __init__(self,a,b):
		self.a    = a
		self.b    = b
		self.f    = (a - b) / a # flattening
		self.e1   = ((a**2 - b**2) / a**2) ** .5 # first eccentricity
		self.e2   = ((a**2 - b**2) / b**2) ** .5 # second eccentricity
	def radiusOfCurvature(self,lat):
		# Meridian radius of curvature (M)
		M         = self.a * (1 - self.e1**2) / (1 - self.e1**2 * sin(lat)**2)**(3/2) 
		# Prime vertical radius of curvature (N)
		N         = self.a / (1 - self.e1**2 * sin(lat)**2)**(1/2)
		# Gaussian mean radius of curvature
		rMean     = (M * N) ** .5
		# Raidus of curvature for parallel of latitude
		r_lat     = N * cos(lat) 
		return M, N, rMean, r_lat
# -----------------------------------------------------------------------------
def ellipsoidInstance(ellipsoidName):
    #
    # @author: Mustafa Serkan Isik
    #
	axes = {
            'GRS80'  : [6378137.000, 6356752.3141],
		    'WGS84'  : [6378137.000, 6356752.3141],
		    'Hayford': [6378388.000, 6356911.9460]
			}
		
	a, b = axes[ellipsoidName]
	return Ellipsoid(a,b)
# -----------------------------------------------------------------------------
class angle:
    # A class dedicated to handle all those angle
    def __init__(self, angle, base = "deg"):
        self.base    = base
        self.angle   = angle
        self.unit    = {   
                         "deg"   : 180,
                         "rad"   : pi,
                         "gon"   : 200
                       }
        self.dms2deg(self.angle) if base == "dms" else 1 # to bypass "dms" base usage error
        self.deg     = self.unit["deg"] * self.angle / self.unit[self.base]
        self.rad     = self.unit["rad"] * self.angle / self.unit[self.base]
        self.gon     = self.unit["gon"] * self.angle / self.unit[self.base]
        self.dms     = self.deg2dms(self.deg)
        del self.base, self.angle
    def deg2dms(self, deg):
        d            = int(deg)
        m            = int((deg - d)*60)
        s            = (deg - d - m/60)*3600
        return "{:3d}d{:02d}m{:03.5f}s".format(d, m, s) # edit to get more digits 
    def dms2deg(self, dms):
        d, m, s      = list(map(float, dms.split(' ')))
        self.deg     = d + m / 60 + s / 3600
        self.base, self.angle = "deg", self.deg
    def __add__(self, other):
        return angle(self.deg + other.deg)
    def __sub__(self, other):
        return angle(self.deg - other.deg)
    def __mul__(self, cons):
        return angle(self.deg * cons)
    def __rmul__(self, cons):
        return angle(cons * self.deg)
    def __truediv__(self, cons):
        return angle(self.deg / cons)
    def __floordiv__(self, cons):
        return angle(self.deg // cons)
    def __pow__(self, cons):
        return angle(self.deg ** cons)
    def __pos__(self):
        return angle(self.deg)
    def __neg__(self):
        return angle(-self.deg)
# -----------------------------------------------------------------------------
def ITxyz2ATxyz(pole, ITxyz):
    """
     |   X* |                               |   X  |
     |   Y* |      ----------+--------->    |   Y  |
     |   Z* |IT              |              |   Z  |AT
                             |
                             |
                         |  xP  |
                         |  yP  |AT
    """
    # given====================================================================
    # ITxyz     ---->   Car. Coor. (X*, Y*, Z*)[3x1] of 
    #                   Station Point on Inst. Terr. Coor. Sys. 
    # pole      ---->   Car. Coor. (xP, yP)[2x1] of 
    #                   Pole Point on the Given Date 
    #
    # asked====================================================================
    # ATxyz     ---->   Car.Coor. (X, Y, Z)[3x1] of 
    #                   Station Point on Avg. Terr. Coor. Sys. 
    #
    return R2(-pole[0,0])*R1(pole[1,0])*ITxyz
"""
# ITxyz2ATxyz Control
ITxyz = matrix([[3987696.2562],
                [2778939.2664],
                [4117699.6076]])
pole  = matrix([[angle("0 0 0.18767", "dms")],
                [angle("0 0 0.11004", "dms")]])
    
ATxyz = ITxyz2ATxyz(pole, ITxyz)
"""
def ATxyz2ITxyz(pole, ATxyz):
    """
     |   X  |                             |   X* |
     |   Y  |      ----------+--------->  |   Y* |
     |   Z  |AT              |            |   Z* |IT
                             |
                             |
                         |  xP  |
                         |  yP  |AT
    """
    # given====================================================================
    # ATxyz     ---->   Car. Coor. (X , Y , Z )[3x1] of 
    #                   Station Point on Avg. Terr. Coor. Sys. 
    # pole      ---->   Car. Coor. (xP, yP)[2x1] of 
    #                   Pole Point on the Given Date 
    #
    # asked====================================================================
    # ITxyz     ---->   Car.Coor. (X*, Y*, Z*)[3x1] of 
    #                   Station Point on Inst. Terr. Coor. Sys. 
    #
    return (R2(pole[0,0]) * R1(pole[1,0])).I * ATxyz
"""
# ATxyz2ATxyz Control
ATxyz = matrix([[3987700.00268801],
                [2778937.06965226],
                [4117697.46192639]])

pole  = matrix([[angle("0 0 0.18767", "dms")],
                [angle("0 0 0.11004", "dms")]])
    
ITxyz = ATxyz2ITxyz(pole, ATxyz)
"""
def ITcur2ATcur(pole, ITcur, A):
    """
     | LAT* |   , Astr.     ----------+--------->  | LAT  |   , Astr.
     | LON* |IT   Azimuth        |                 | LON  |AT   Azimuth 
                  on IT CS.      |                              on AT CS.
                  (A*)           |                              (A)
                             |  xP  |
                             |  yP  |AT
    """
    # given====================================================================
    # ITcur     ---->   Cur. Coor. (LAT*, LON*)[2x1] of 
    #                   Station Point on Inst. Terr. Coor. Sys. 
    # pole      ---->   Car. Coor. (xP, yP)[2x1] of 
    #                   Pole Point on the Given Date 
    #
    # asked====================================================================
    # ATcur     ---->   Cur. Coor. (LAT , LON)[2x1] of 
    #                   Station Point on Avg. Terr. Coor. Sys. 
    #
    dA      = pole.transpose()*(matrix([[sin(ITcur[1,0].rad)], [cos(ITcur[1,0].rad)]])/-cos(ITcur[0,0].rad))
    
    dCur    = pole.transpose()*matrix([[-cos(ITcur[1,0].rad), -sin(ITcur[1,0].rad)*tan(ITcur[0,0].rad)],
                                       [ sin(ITcur[1,0].rad), -cos(ITcur[1,0].rad)*tan(ITcur[0,0].rad)]])
    return ITcur+dCur.transpose(), A+dA[0,0]
"""      
# ITcur2ATcur Control
ITcur = matrix([[angle("40 27 21.7251", "dms")],
                [angle("34 52 32.3467", "dms")]])
A     = angle("74 32 56.7409", "dms")

pole  = matrix([[angle("0 0 0.18676", "dms")],
                [angle("0 0 0.11004", "dms")]])
    
ATcur = ITcur2ATcur(pole, ITcur, A)
"""
def Gcur2Gxyz(Gcur, ell = 'GRS80'):
    """
    @author: Mustafa Serkan Isik
     | lat  |                             |   x  |
     | lon  |      ----------+--------->  |   y  |
     |  h   |G               |            |   z  |G
                             |
                             |
                   <DATUM PARAMETERS>
    """
    # given====================================================================
    # Gcur      ---->   Cur.Coor. (lat, lon, h)[3x1] of
    #           |||     Station Point on Geodetic Coor. Sys.
    #           |||
    #           ||----lat   ----> Ellipsoidal latitude on Geodetic
    #           ||                  Coor. Sys.
    #           ||
    #           |-----lon   ----> Ellipsoidal longtitude on Geodetic
    #           |                  Coor. Sys.
    #           |
    #           ----- h     ----> Ellipsoidal height on Geodetic
    #                             Coor. Sys.
    #
    # asked====================================================================
    # Gxyz      ---->   Car.Coor. (xG, yG, zG)[3x1] of
    #                   Station Point on Geodetic Coor. Sys. 
    #
    ell_    = ellipsoidInstance(ell)  # create an ellipsoid instance
    N       = ell_.radiusOfCurvature(Gcur[0,0].rad)[1] # Radius of curvature in the prime vertical N
    return matrix([[   (N+Gcur[2,0])*cos(Gcur[0,0].rad)*cos(Gcur[1,0].rad)],
                   [   (N+Gcur[2,0])*cos(Gcur[0,0].rad)*sin(Gcur[1,0].rad)],
                   [   ((1-ell_.e1**2)*N+Gcur[2,0])*sin(Gcur[0,0].rad)    ]])
"""
# Gcur2Gxyz Control    
Gcur = matrix([[angle("40 27 35.8741", "dms")],
               [angle("34 52 23.6714", "dms")],
               [623.845]])
Gxyz = Gcur2Gxyz(Gcur, ell = 'Hayford')
"""
def Gxyz2GcurI(Gxyz, iter_ = 999,ell = 'GRS80'):
    """
    @author: Mustafa Serkan Isik
     |   x  |                            +--------> | lato |
     |   y  |      ----------+-----------+(iter.)   | lono |
     |   z  |G               |           |          |  ho  |G
                             |           +--------<
                             |                          |
                   <DATUM PARAMETERS>                   v
                                                    | lat  |
                                                    | lon  |
                                                    |  h   |
    """
    # given====================================================================
    # Gxyz      ---->   Car. Coor. (x, y, z)[3x1] of 
    #                   Station Point on Geodetic Coor. Sys.
    # ell       ---->   Reference ellipsoid
    #
    # asked====================================================================
    # Gcur      ---->   Cur.Coor. (lat, lon, h)[3x1] of
    #           |||     Station Point on Geodetic Coor. Sys.
    #           |||
    #           ||----lat   ----> Ellipsoidal latitude on Geodetic
    #           ||                Coor. Sys.
    #           ||
    #           |-----lon   ----> Ellipsoidal longtitude on Geodetic
    #           |                 Coor. Sys.
    #           |
    #           ----- h     ----> Ellipsoidal height on Geodetic
    #                             Coor. Sys.
    #
    ell_    = ellipsoidInstance(ell)  # create an ellipsoid instance
    lon     = angle(atan(Gxyz[1,0], Gxyz[0,0]), "rad")
    p       = (Gxyz[0,0]**2+Gxyz[1,0]**2)**.5
    N       = ell_.a
    h       = (Gxyz[0,0]**2+Gxyz[1,0]**2+Gxyz[2,0]**2)**.5-(ell_.a * ell_.b)**.5
    lat     = angle(atan(Gxyz[2,0], p*(1-ell_.e1**2*N/(N  + h ))), "rad")
    i = 0
    while True:
        print("===========================[{:02d}]".format(i))
        N   = ell_.a/(1-ell_.e1**2*sin(lat.rad)**2)**.5
        h   = (p/cos(lat.rad))-N
        if (abs(atan(Gxyz[2,0], p*(1-ell_.e1**2*N/(N+h)))-lat.rad) < 1e-10  and abs((p/cos(atan(Gxyz[2,0], p*(1-ell_.e1**2*N/(N+h)))))-N-h) < 1e-10) or iter_ == i:
            lat = angle(atan(Gxyz[2,0], p * (1 - ell_.e1**2 * N / (N + h))), "rad")
            h   = (p / cos(lat.rad)) - N
            print("lat[{:02d}]\t:  {}".format(i, lat.dms))
            print("lon[{:02d}]\t:  {}".format(i, lon.dms))
            print("h  [{:02d}]\t:  {}".format(i, h))
            print("===========================[##]")
            return matrix([[lat],
                           [lon],
                           [h]])
        print(N, h, lat)
        lat = angle(atan(Gxyz[2,0], p*(1-ell_.e1**2*N/(N+h))), "rad")
        print("lat[{:02d}]\t:  {}".format(i, lat.dms))
        print("lon[{:02d}]\t:  {}".format(i, lon.dms))
        print("h  [{:02d}]\t:  {}".format(i, h))
        i += 1
"""
# Gxyz2GcurI Control        
Gxyz = matrix([[3987577.612],
               [2779005.035],
               [4117452.563]])
a = Gxyz2GcurI(Gxyz, ell = 'Hayford')
"""
def Gxyz2GcurII(Gxyz, ell = 'GRS80'):
    """
    @author: Mustafa Serkan Isik
     |   x  |                  (direct)   | lat  |  
     |   y  |      ----------+--------->  | lon  |
     |   z  |G               |            |  h   |G
                             |
                             |
                   <DATUM PARAMETERS>
    """
    # given====================================================================
    # Gxyz      ---->   Car. Coor. (x, y, z)[3x1] of 
    #                   Station Point on Geodetic Coor. Sys.
    # ell       ---->   Reference ellipsoid
    #
    # asked====================================================================
    # Gcur      ---->   Cur.Coor. (lat, lon, h)[3x1] of
    #           |||     Station Point on Geodetic Coor. Sys.
    #           |||
    #           ||----lat   ----> Ellipsoidal latitude on Geodetic
    #           ||                Coor. Sys.
    #           ||
    #           |-----lon   ----> Ellipsoidal longtitude on Geodetic
    #           |                 Coor. Sys.
    #           |
    #           ----- h     ----> Ellipsoidal height on Geodetic
    #                             Coor. Sys.
    #
    ell_    = ellipsoidInstance(ell)  # create an ellipsoid instance
    lon     = atan(Gxyz[1,0],Gxyz[0,0])
    p       = (Gxyz[0,0]**2+Gxyz[1,0]**2)**.5
    B       = atan(ell_.a*Gxyz[2,0], ell_.b*p)
    lat     = atan(Gxyz[2,0]+ell_.e2**2*ell_.b*sin(B)**3, p-ell_.e1**2*ell_.a*cos(B)**3)
    h       = p/cos(lat)-ell_.radiusOfCurvature(lat)[1] 
    return matrix([[angle(lat, "rad")],
                   [angle(lon, "rad")],
                   [h]])
"""
# Gxyz2GcurII Control        
Gxyz = matrix([[3987577.612],
               [2779005.035],
               [4117452.563]])
a = Gxyz2GcurII(Gxyz, ell = 'Hayford')
"""
def Gxyz2ATxyz(Gxyz, C, fo, E):
    """
    @author: Mustafa Serkan Isik
     |   x  |                             |   X  |  
     |   y  |      ----------+--------->  |   Y  |
     |   z  |G               |            |   Z  |AT
                             |
                             |
               Shifting,   Scale,   Rotation
               vector      factor   angles
               (C)         (fo)     (ex, ey, ez)
    """
    # given====================================================================
    # C         ---->   Shifting matrix (xo, yo, zo)[3x1]
    # fo        ---->   Scale factor
    # Gxyz      ---->   Coor.(x, y, z)[3x1] on the Geodetic Coor. Sys.
    # E         ---->   Axes rotation angles(eX, eY, eZ)[3x1]
    #
    # asked====================================================================
    # Axyz      ----> Coor. (X, Y, Z)[3x1] on the Avg. Terr. Coor. Sys.
    return C + R(1+fo,E[0,0].rad, E[1,0].rad, E[2,0].rad)*Gxyz
"""
# Gxyz2ATxyz Control  
Gxyz = matrix([[3987577.612],
               [2779005.035],
               [4117452.563]])
C    = matrix([[135.763],
               [ 81.345],
               [113.928]])
fo   = .00000175
E    = matrix([[angle("0 0 -2.3", "dms")],
               [angle("0 0 4.8", "dms")],
               [angle("0 0 5.6", "dms")]])
a = Gxyz2ATxyz(Gxyz, C, fo, E)
"""
def ATxyz2Gxyz(ATxyz, C, fo, E):
    """
     |   X  |                             |   x  |  
     |   Y  |      ----------+--------->  |   y  |
     |   Z  |AT              |            |   z  |G
                             |
                             |
               Shifting,   Scale,   Rotation
               vector      factor   angles
               (C)         (fo)     (ex, ey, ez)
    """
    # given====================================================================
    # C         ---->   Shifting matrix (xo, yo, zo)[3x1]
    # fo        ---->   Scale factor
    # AT        ---->   Coor.(X, Y, Z)[3x1] on Avg. Terr. Coor. Sys.
    # E         ---->   Rotation angles(eX, eY, eZ)[3x1]
    #
    # asked====================================================================
    # G         ----> Coor. (x, y, z)[3x1] on Geodetic Coor. Sys. 
    
    return R(1+fo,E[0,0].rad, E[1,0].rad, E[2,0].rad).I*(ATxyz-C)
"""
# ATxyz2Gxyz Control
ATxyz= matrix([[3987699.98198479],
               [2778937.06721651],
               [4117697.47897314]])
C    = matrix([[135.763],
               [ 81.345],
               [113.928]])
fo   = .00000175
E    = matrix([[angle("0 0 -2.3", "dms")],
               [angle("0 0 4.8", "dms")],
               [angle("0 0 5.6", "dms")]])
a = ATxyz2Gxyz(ATxyz, C, fo, E)  
"""
def DOV(ATcur, A, Gcur, a):
    """
     |  LAT |      ----------+---------- Astr. 
     |  LON |                |           Azimuth
     |   H  |AT              |           (A)
                             |
                             |
                             +-----------------> Deflection of   ,  Geoid
                             |                   Vertices           ondulation
                             |                   (E, nlon, na)      (N)
                             |
     |  lat |      ----------+---------- Ell. 
     |  lon |                            Azimuth   
     |   h  |G                           (a)
    """
    # given====================================================================
    # ATcur     ---->   Cur.Coor. (lat_a, lon_a, h_ort)[3x1] of
    #           |||     Station Point on Avg. Terr. Coor. Sys.
    #           |||
    #           ||----lat_a ----> Astronomic latitude on Avg. Terr.
    #           ||                Coor. Sys.
    #           ||
    #           |-----lon_a ----> Astronomic latitude on Avg. Terr.
    #           |                 Coor. Sys.
    #           |
    #           ----- h_ort ----> Orthometric height on Avg. Terr.
    #                             Coor. Sys.
    # A         ---->   Astronomic azimuth angle from point P to K
    # Gcur      ---->   Cur.Coor. (lat, lon, h)[3x1] of
    #           |||     Station Point on Geodetic Coor. Sys.
    #           |||
    #           ||----lat   ----> Ellipsoidal latitude on Geodetic
    #           ||                  Coor. Sys.
    #           ||
    #           |-----lon   ----> Ellipsoidal longtitude on Geodetic
    #           |                  Coor. Sys.
    #           |
    #           ----- h     ----> Ellipsoidal height on Geodetic
    #                             Coor. Sys.
    # a         ---->   Ellipsoidal azimuth angle from point P to K
    #
    # asked====================================================================
    # N         ---->   Geoid Ondulation
    # E         ---->   North-south component
    # nlon, na  ---->   East-west components
    #
    dCur = ATcur - Gcur
    E, nlon, na, N = dCur[0,0], dCur[1,0]*cos(Gcur[0,0].rad), (A-a)/tan(Gcur[0,0].rad), -dCur[2,0]
    return E, nlon, na, N
"""
# DOV Control
ATcur       = matrix([[angle("40 27 21.6348" ,"dms")],
                      [angle("34 52 32.1786" ,"dms")],
                      [587.218]])

Gcur        = matrix([[angle("40 27 35.8741" ,"dms")],
                      [angle("34 52 23.6714" ,"dms")],
                      [623.845]])

A, a        = angle("74 32 56.4819" ,"dms"), angle("74 32 50.9614" ,"dms")
b = DOV(ATcur,  A, Gcur, a)
"""
def LAxyz(Apk, Bpk, Spk, LAT, k = .13, ell = 'GRS80'):
    """
    Astronomic      ----------+     Ell. parameters   
    Azimuth                   |     the measurements    
    (Apk),                    |     are done
                              |     (ell)
                              |     |
                              |     |               | LAT_ |
    Zenith          ----------+--+--+----------->   | LON_ |
    Angle                     |  |  |               |  H_  |LA_K(Obs. Point)
    w/o corr.                 |  |  |
    (Bpk),                    |  |  |
                              |  |  Refraction
    Slope          -----------+  |  Constant                
    Distance                     |  (k)
    (Spk)                        |
                                 Astronomic Latitude of
   +-------------------+         sta. point(P) on 
  *| values are        |         Avg. Terr. Coor. Sys.
   | measured from     |         (LAT)
   | sta. point(P) on  |         
   | local astronomic  |     
  *| coor. sys.        |*
   +-------------------+
    """
    # given====================================================================
    # LAT       ---->   Astronomic latitude of 
    #                   station point(P) on Local Ast. Coor. Sys.
    # Spk       ---->   Slope distance between point P and K
    # Apk       ---->   Astronomic azimuth angle from point P to K
    # Bpk       ---->   Zenith angle (w/o any corrections)
    # k         ---->   Refraction constant
    # ell       ---->   Reference ellipsoid
    #
    # asked====================================================================
    # LAxyz      ---->  Car. Coor. (xG, yG, zG)[3x1] of 
    #                   Obs. Point(K) on Local Ast. Coor. Sys.
    #
    ell_    = ellipsoidInstance(ell)
    M, N    = ell_.radiusOfCurvature(LAT.rad)[:2]
    Rpk     = M*N/(N*cos(Apk.rad)**2+M*sin(Apk.rad)**2)
    dB      = angle(Spk*k/(2*Rpk), "rad")
    Bpk_    = Bpk+dB 
    return Spk*matrix([[sin(Bpk_.rad)*cos(Apk.rad)],
                       [sin(Bpk_.rad)*sin(Apk.rad)],
                       [cos(Bpk_.rad)]])
"""
# LAxyz Control
LAcur_P            = matrix([[angle("40 27 21.6348" ,"dms")],
                             [angle("34 52 32.1786" ,"dms")],
                             [587.218]])

LAxyz_P            = matrix([[3987699.9846],
                             [2778937.0697],
                             [4117697.4796]])
    
Apk, Bpk, Spk      = angle("74 32 56.4819", "dms"), angle("88 53 25.6000", "dms"), 13794.682

LAxyz_K            = LAxyz(Apk, Bpk, Spk, LAcur_P[0,0], k = .13, ell = 'Hayford')
"""
def LAxyz2ATxyz(LAcur_P, LAxyz_P, LAxyz_K):
    """
    | LAT_ |        ----------+
    | LON_ |                  |
    |  H_  |LA_P              |
                              |   
                              |      
                              |      
                              |      
                              |      
                              |                      |  X   |
    |  X_  |        ----------+---+-------------->   |  Y   |
    |  Y_  |                      |                  |  Z   |AT_K(Obs. Point)
    |  Z_  |LA_P                  |                  
                                  |   
                                  |  
                                  |                 
                                |  X_  |      
                                |  Y_  |
                                |  Z_  |LA_K
    """              
    # given====================================================================
    # LAcur_P   ---->   Cur.Coor. (latLA, lonLA)[2x1] of
    #            ||     Station Point(P) on Local Ast. Coor. Sys.
    #            ||
    #            |----lat_a ----> Astronomic latitude on Local Ast.
    #            |                Coor. Sys.
    #            |
    #            -----lon_a ----> Astronomic longtitude on Local Ast.
    #                             Coor. Sys.
    # LAxyz_P   ---->   Car.Coor. (xLA, yLA, zLA)[3x1] of
    #                   Station Point(P) on Local Ast. Coor. Sys.
    # LAxyz_K   ---->   Car.Coor. (xLA, yLA, zLA)[3x1] of
    #                   Obs. Point(K) on Local Ast. Coor. Sys.
    #
    # asked====================================================================
    # ATxyz     ---->   Car.Coor. (xLA, yLA, zLA)[3x1] of
    #                   Obs. Point(K) on Avg. Terr. Coor. Sys.
    #
    return LAxyz_P+R3(angle(pi, "rad")-LAcur_P[1,0])*R2(angle(pi/2, "rad")-LAcur_P[0,0])*P2*LAxyz_K
"""
#LAxyz2ATxyz control

LAcur_P          = matrix([[angle("40 27 21.6348" ,"dms")],
                           [angle("34 52 32.1786" ,"dms")],
                           [587.218]])

LAxyz_P          = matrix([[3987699.9846],
                           [2778937.0697],
                           [4117697.4796]])
    
LAxyz_K          = matrix([[ 3674.41435092],
                           [13293.66792893],
                           [  265.18607453]])
ATxyz_K = LAxyz2ATxyz(LAcur_P, LAxyz_P,  LAxyz_K)
"""
def ATxyz2LAxyz(ATcur_P, ATxyz_P, ATxyz_K):
    """
    |  X   |        ----------+
    |  Y   |                  |
    |  Z   |AT_P(Sta. P)      |
                              |   
                              |      
                              |      
                              |      
                              |      
                              |                      |  X_  |
    |  X   |        ----------+---+-------------->   |  Y_  |
    |  Y   |                      |                  |  Z_  |LA_K(Obs. P)
    |  Z   |AT_K(Obs. P)          |                  
                                  |   
                                  |  
                                  |                 
                                | LAT  |      
                                | LON  |
                                |  H   |AT_P
   """
    # given====================================================================
    # ATcur_P   ---->   Cur.Coor. (lat_a, lon_a, H)[3x1] of
    #           |||     Station Point(P) on Avg. Terr. Coor. Sys.
    #           |||
    #           ||----lat_a ----> Astronomic latitude on Avg. Terr.
    #           ||                Coor. Sys.
    #           ||
    #           |-----lon_a ----> Astronomic longtitude on Avg. Terr.
    #           |                 Coor. Sys.
    #           |
    #           ----- H     ----> Orthometric height on Avg. Terr.
    #                             Coor. Sys.
    # ATxyz_P   ---->   Car.Coor. (xLA, yLA, zLA)[3x1] of
    #                   Station Point(P) on Avg. Terr. Coor. Sys.
    # ATxyz_K   ---->   Car.Coor. (xLA, yLA, zLA)[3x1] of
    #                   Obs. Point(K) on Avg. Terr. Coor. Sys.
    #
    # asked====================================================================
    # LAxyz_K   ---->   Car.Coor. (xLA, yLA, zLA)[3x1] of
    #                   Obs. Point(K) on Local Ast. Coor. Sys.
    #
    return (R3(angle(pi, "rad")-ATcur_P[1,0])*R2(angle(pi/2, "rad")-ATcur_P[0,0])*P2).transpose()*(ATxyz_K-ATxyz_P)
"""
#ATxyz2LAxyz control

ATcur_P     = matrix([[angle("40 27 21.6348","dms")],
                      [angle("34 52 32.1786","dms")],
                      [587.218]])
    
ATxyz_P     = matrix([[3987699.9846],
                      [2778937.0697],
                      [4117697.4796]])
    
ATxyz_K     = matrix([[3978308.2677],
                      [2788595.2368],
                      [4120665.4274]])
    
LAxyz_K     = ATxyz2LAxyz(ATcur_P, ATxyz_P, ATxyz_K)   
"""
def ASBI(LAxyz_K):
    """    
    |  X_  |                                 Astr. Azimuth , Zenith   , Slope
    |  Y_  |        -----------------------> from P to K     distance   distance
    |  Z_  |LA_K (Observed from              (Apk)           (Bpk)      (Spk)
                  station point P)
   """
    # given====================================================================
    # ATxyz_K   ---->   Car.Coor. (xLA, yLA, zLA)[3x1] of
    #                   Obs. Point(K) calculated by using 
    #                   measurement values from station 
    #                   point P on Local Ast. Coor. Sys.
    #
    # asked====================================================================
    # Apk       ---->   Astronomic azimuth
    # Spk       ---->   Distance between sta. and obs. point
    # Bpk       ---->   Zenith distance
    #
    Apk     = atan(LAxyz_K[1,0], LAxyz_K[0,0])
    Spk     = (LAxyz_K[0,0]**2+LAxyz_K[1,0]**2+LAxyz_K[2,0]**2) ** .5
    Bpk     = acos(LAxyz_K[2,0]/Spk)
    return angle(Apk, "rad"), Spk, angle(Bpk, "rad")
"""
#ASBI control

ATcur_P     = matrix([[angle("40 27 21.6348","dms")],
                      [angle("34 52 32.1786","dms")],
                      [587.218]])
    
ATxyz_P     = matrix([[3987699.9846],
                      [2778937.0697],
                      [4117697.4796]])
    
ATxyz_K     = matrix([[3978308.2677],
                      [2788595.2368],
                      [4120665.4274]])
    
LAxyz_K     = ATxyz2LAxyz(ATcur_P, ATxyz_P, ATxyz_K)
Apk, Spk, Bpk = ASBI(LAxyz_K)  
"""
def ASBII(ATcur_P, ATxyz_P, ATxyz_K):
    """    
    |  X   |                                 Astr. Azimuth , Zenith   , Slope
    |  Y   |        ---------------+-------> from P to K     distance   distance
    |  Z   |AT_P                   |         (Apk)           (Bpk)      (Spk)
                                   |
                                   |
                                   |
                                   |                | LAT  |   
                                   +--------------- | LON  |
                                   |                |  H   |AT_P
                                   |
    |  X   |                       |
    |  Y   |        ---------------+
    |  Z   |AT_K                            
   """
    # given====================================================================
    # ATcur_P   ---->   Cur.Coor. (LAT, LON, H)[3x1] of
    #           |||     Station Point(P) on Avg. Terr. Coor. Sys.
    #           |||
    #           ||----LAT ----> Astronomic latitude on Avg. Terr.
    #           ||              Coor. Sys.
    #           ||
    #           |-----LON ----> Astronomic longtitude on Avg. Terr.
    #           |               Coor. Sys.
    #           |
    #           ----- H   ----> Orthometric height on Avg. Terr.
    #                           Coor. Sys.
    # ATxyz_P   ---->   Car.Coor. (X, Y, Z)[3x1] of
    #                   Station Point(P) on Avg. Terr. Coor. Sys.
    # ATxyz_K   ---->   Car.Coor. (X, Y, Z)[3x1] of
    #                   Obs. Point(K) on Avg. Terr. Coor. Sys.
    #
    # asked====================================================================
    # Apk       ---->   Astronomic azimuth
    # Spk       ---->   Distance between sta. and obs. point
    # Bpk       ---->   Zenith distance
    #
    d       = ATxyz_K-ATxyz_P
    Apk     = atan(-d[0,0]*sin(ATcur_P[1,0].rad)+d[1,0]*cos(ATcur_P[1,0].rad), -d[0,0]*sin(ATcur_P[0,0].rad)*cos(ATcur_P[1,0].rad)-d[1,0]*sin(ATcur_P[0,0].rad)*sin(ATcur_P[1,0].rad)+d[2,0]*cos(ATcur_P[0,0].rad))
    Spk     = (d[0,0]**2+d[1,0]**2+d[2,0]**2)**.5
    Bpk     = acos((d[0,0]*cos(ATcur_P[0,0].rad)*cos(ATcur_P[1,0].rad)+d[1,0]*cos(ATcur_P[0,0].rad)*sin(ATcur_P[1,0].rad)+d[2,0]*sin(ATcur_P[0,0].rad))/Spk)
    return angle(Apk, "rad"), Spk, angle(Bpk, "rad")
"""
#ASBII control

ATcur_P     = matrix([[angle("40 27 21.6348","dms")],
                      [angle("34 52 32.1786","dms")],
                      [587.218]])
    
ATxyz_P     = matrix([[3987699.9846],
                      [2778937.0697],
                      [4117697.4796]])
    
ATxyz_K     = matrix([[3978308.2677],
                      [2788595.2368],
                      [4120665.4274]])

Apk, Spk, Bpk = ASBII(ATcur_P, ATxyz_P, ATxyz_K)
"""
def LGxyz(apk, Bpk, Spk, e, n):
    """
    Ellipsoidal     ----------+     Deflection of
    Azimuth                   |     vertices    
    (apk),                    |     (e, n)
                              |     |
                              |     |
                              |     |               |  x  |
    Zenith          ----------+-----+----------->   |  y  |
    Angle                     |                     |  z  |LG_K(Obs. Point)
    w/o corr.                 |
    (Bpk),                    |
                              |  
    Slope           ----------+              
    Distance                     
    (Spk)                        
                                 
   +-------------------+         
  *| values are        |         
   | measured from     |         
   | sta. point(P) on  |         
   | local geodetic    |     
  *| coor. sys.        |*
   +-------------------+
   """
    # given====================================================================
    # Spk       ---->   Slope distance between point P and K
    # apk       ---->   Astronomic azimuth angle from point P to K
    # Bpk       ---->   Zenith angle (w/o any corrections)
    # e, n      ---->   Deflection of vertices
    #
    # asked====================================================================
    # LGxyz     ---->   Car. Coor. (xG, yG, zG)[3x1] of 
    #                   Obs. Point(K) on Local Geo. Coor. Sys.
    #
    Bpk_    = Bpk+e*cos(apk.rad)+n*sin(apk.rad)
    
    return Spk*matrix([[sin(Bpk_.rad)*cos(apk.rad)],
                       [sin(Bpk_.rad)*sin(apk.rad)],
                       [cos(Bpk_.rad)]])
"""
#ASBII control
apk         = angle("74 32 50.9614", "dms")
Bpk         = angle("88 53 54.56307", "dms")
Spk         = 13794.682
e, n        = angle("0 0 -14.2393", "dms"), angle("0 0 6.4728", "dms")

LGxyz_K     = LGxyz(apk, Bpk, Spk, e, n)
"""
def LGxyz2Gxyz(Gcur_P, Gxyz_P, LGxyz_K):
    """    
    |  x   |                              |  x   |
    |  y   |        ----------+---------> |  y   |
    |  z   |G_P               |           |  z   |G_P
                              |
                              |
                              |
                              |           | lat  |   
                              +---------- | lon  |
                              |           |  h   |G_P
                              |
    |  x_  |                  |
    |  y_  |        ----------+
    |  z_  |LG_K                              
                                
                                
   """
    # given====================================================================
    # Gcur_P    ---->   Cur.Coor. (lat, lon, h)[3x1] of
    #           |||     Station Point(P) on Geo. Coor. Sys.
    #           |||
    #           ||----lat ----> Astronomic latitude on Geo.
    #           ||              Coor. Sys.
    #           ||
    #           |-----lon ----> Astronomic longtitude on Geo.
    #           |               Coor. Sys.
    #           |
    #           ------h   ----> Orthometric height on Geo.
    #                           Coor. Sys.
    # Gxyz_P    ---->   Car.Coor. (x, y, z)[3x1] of
    #                   Station Point(P) on Geo. Coor. Sys.
    # LGxyz_K   ---->   Car.Coor. (x_, y_, z_)[3x1] of
    #                   Obs. Point(K) on Local Geo. Coor. Sys.
    #
    # asked====================================================================
    # Gxyz_K    ---->   Car.Coor. (x, y, z)[3x1] of
    #                   Obs. Point(K) on  Geo. Coor. Sys.
    #
    return Gxyz_P+R3(angle(pi, "rad")-Gcur_P[1,0])*R2(angle(pi/2, "rad")-Gcur_P[0,0])*P2*LGxyz_K
"""
#LGxyz2Gxyz control
apk         = angle("74 32 50.9614", "dms")
Bpk         = angle("88 53 54.56307", "dms")
Spk         = 13794.682
e, n        = angle("0 0 -14.2393", "dms"), angle("0 0 6.4728", "dms")

Gcur_P      = matrix([[angle("40 27 35.8741", "dms")],
                      [angle("34 52 23.6714", "dms")],
                      [623.845]])
    
Gxyz_P      = matrix([[3987577.612],
                      [2779005.035],
                      [4117452.563]])
    
LGxyz_K     = LGxyz(apk, Bpk, Spk, e, n)

Gxyz_K      = LGxyz2Gxyz(Gcur_P, Gxyz_P, LGxyz_K)
"""
def Gxyz2LGxyz(Gcur_P, Gxyz_P, Gxyz_K):
    """
    |  x   |        ----------+
    |  y   |                  |
    |  z   |G_P(Sta. P)       |
                              |   
                              |      
                              |      
                              |      
                              |      
                              |                      |  x_  |
    |  x   |        ----------+---+-------------->   |  y_  |
    |  y   |                      |                  |  z_  |LA_K
    |  z   |G_K(Obs. P)           |                  
                                  |   
                                  |  
                                  |                 
                                | lat  |      
                                | lon  |
                                |  h   |G_P
    """
    # given====================================================================
    # Gcur_P    ---->   Cur.Coor. (lat, lon, h)[3x1] of
    #           |||     Station Point(P) on Geo. Coor. Sys.
    #           |||
    #           ||----lat ----> Astronomic latitude on Geo.
    #           ||              Coor. Sys.
    #           ||
    #           |-----lon ----> Astronomic longtitude on Geo.
    #           |               Coor. Sys.
    #           |
    #           ------h   ----> Orthometric height on Geo.
    #                           Coor. Sys.
    # Gxyz_P    ---->   Car.Coor. (x, y, z)[3x1] of
    #                   Station Point(P) on Geo. Coor. Sys.
    # Gxyz_K    ---->   Car.Coor. (x, y, z)[3x1] of
    #                   Obs. Point(K) on Geo. Coor. Sys.
    #
    # asked====================================================================
    # LGxyz_K    ---->   Car.Coor. (x_, y_, z_)[3x1] of
    #                   Obs. Point(K) on  Local Geo. Coor. Sys.
    #
    return (R3(angle(pi, "rad")-Gcur_P[1,0])*R2(angle(pi/2, "rad")-Gcur_P[0,0])*P2).transpose()*(Gxyz_K-Gxyz_P)
"""
#Gxyz2LGxyz  control
Gxyz_P      = matrix([[3987577.6120],
                      [2779005.0350],
                      [4117452.5630]])
    
Gxyz_K      = matrix([[3978185.8880],
                      [2788663.1907],
                      [4120420.5253]])
    
Gcur_PI     = Gxyz2GcurI (Gxyz_P, ell = 'Hayford')
Gcur_PII    = Gxyz2GcurII(Gxyz_P, ell = 'Hayford')

Gcur_KI     = Gxyz2GcurI (Gxyz_K, ell = 'Hayford')
Gcur_KII    = Gxyz2GcurII(Gxyz_K, ell = 'Hayford')

LGxyz_KI     = Gxyz2LGxyz(Gcur_PI, Gxyz_P, Gxyz_K)
LGxyz_KII     = Gxyz2LGxyz(Gcur_PII, Gxyz_P, Gxyz_K)
"""
def LAxyz2LGxyz(LAxyz_K, e, n, dA):
    """
    |  X_ |         ----------+     Deflection of
    |  Y_ |                   |     vertices    
    |  Z_ |LA_K               |     (e, n)
                              |     |
                              |     |
                              |     |               |  x_ |
                              +-----+----------->   |  y_ |
                              |                     |  z_ |LG_K(Obs. Point)
                              |
                              |
                            Diff. between
                            Astr. and Geo.(Ell.)
                            azimuths
                            (Apk - apk)[dA]
    """
    # given====================================================================
    # LAxyz_K   ---->   Car.Coor. (X_, Y_, Z_)[3x1] of
    #                   Station Point(P) on Local Astr. Coor. Sys.
    # e, n      ---->   Deflection of vertices
    # dA        ---->   Diff. between Astr. and Geo.(Ell.) azimuths
    #
    # asked====================================================================
    # LGxyz_K   ---->   Car.Coor. (x_, y_, z_)[3x1] of
    #                   Obs. Point(K) on a Local Geo. 
    #                   Coor. Sys. established at point P
    #
    return R(1, n, -e, dA)*LAxyz_K
"""
#LAxyz2LGxyz  control
LAxyz_K     = matrix([[ 3674.4144],
                      [13293.6679],
                      [  265.1862]]) 
dA          = angle("0 0 5.5205", "dms")
e, n        = angle("0 0 -14.2393", "dms"), angle("0 0 6.4728", "dms")
"""
def LGxyz2LAxyz(LAxyz_K, e, n, dA):
    """
    |  x_ |         ----------+     Deflection of
    |  y_ |                   |     vertices    
    |  z_ |LG_K               |     (e, n)
                              |     |
                              |     |
                              |     |               |  X_ |
                              +-----+----------->   |  Y_ |
                              |                     |  Z_ |LA_K(Obs. Point)
                              |
                              |
                            Diff. between
                            Astr. and Geo.(Ell.)
                            azimuths
                            (Apk - apk)[dA]
    """
    # given====================================================================
    # LGxyz_K   ---->   Car.Coor. (x_, y_, z_)[3x1] of
    #                   Station Point(P) on Local Geo. Coor. Sys.
    # e, n      ---->   Deflection of vertices
    # dA        ---->   Diff. between Astr. and Geo.(Ell.) azimuths
    #
    # asked====================================================================
    # LAxyz_K   ---->   Car.Coor. (X_, Y_, Z_)[3x1] of
    #                   Obs. Point(K) on a Local Astr. 
    #                   Coor. Sys. established at point P
    #
    return R(1, n, -e, dA).I*LAxyz_K
# -----------------------------------------------------------------------------
def findFunc(params):
    """                              *
    *         *                   *       *
  * All those  *             *    A full list  *
    parameters                 *  of funcs  *
 *  you don't know -----+  +----> you might  *
    where to use*       |  |      be looking  
          *             |  |      for  *
      *                 |  |         *  *
                      * v  ^ *
                  *    MAGIC   *
                 *   DONE HERE   *
                    *  *    *
    """
    funcs = {
                "ITxyz2ATxyz"   : ITxyz2ATxyz.__code__.co_varnames  ,
                "ATxyz2ITxyz"   : ATxyz2ITxyz.__code__.co_varnames  ,
                "ITcur2ATcur"   : ITcur2ATcur.__code__.co_varnames  ,
                "Gcur2Gxyz"     : Gcur2Gxyz.__code__.co_varnames    ,
                "Gxyz2GcurI"    : Gxyz2GcurI.__code__.co_varnames   ,
                "Gxyz2GcurII"   : Gxyz2GcurII.__code__.co_varnames  ,
                "Gxyz2ATxyz"    : Gxyz2ATxyz.__code__.co_varnames   ,
                "ATxyz2Gxyz"    : ATxyz2Gxyz.__code__.co_varnames   ,
                "DOV"           : DOV.__code__.co_varnames          ,
                "LAxyz"         : LAxyz.__code__.co_varnames        ,
                "LAxyz2ATxyz"   : LAxyz2ATxyz.__code__.co_varnames  ,
                "ATxyz2LAxyz"   : ATxyz2LAxyz.__code__.co_varnames  ,
                "ASBI"          : ASBI.__code__.co_varnames         ,
                "ASBII"         : ASBII.__code__.co_varnames        ,
                "LGxyz"         : LGxyz.__code__.co_varnames        ,
                "LGxyz2Gxyz"    : LGxyz2Gxyz.__code__.co_varnames   ,
                "Gxyz2LGxyz"    : Gxyz2LGxyz.__code__.co_varnames   ,
                "LAxyz2LGxyz"   : LAxyz2LGxyz.__code__.co_varnames  ,
                "LGxyz2LAxyz"   : LGxyz2LAxyz.__code__.co_varnames  ,
                "Just Wandering": ()
            }
    coverage = []
    for func in funcs:
        diff = seq(None,params,funcs[func])
        coverage.append([func, diff.ratio()]) if diff.ratio() != 0 else 1
        coverage.sort(key = lambda l: l[1], reverse = True)
    print("You might be looking for...")
    for l in coverage:
        print("\t[{:3.2f}%] {:12s} ".format(l[1]*100, l[0]))

