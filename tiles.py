import math
import numpy as np

# map ant settings begin
UpperLeftX = -76000
UpperLeftY = 7818000
LowerRightX= 739460
LowerRightY=6474000
SizeX=556687
SizeY=917504
# map and settings end
# ETRS-TM35FIN begin
a = 6378137.0
b = 6356752.314245
fm = 298.257223563
f = 1/ fm
n = (a-b)/(a+b)
A1= (a/(1+n))*(1 + n**2 / 4 + n**4 / 64)
e2 = 2*f - f**2
e=math.sqrt(e2)
edot2 = e2/(1-e2)
mylambda0=27.0

h1=n/2 -2/3*n**2 + 37/96*n**3 - 1/360*n**4
h2=     1/48*n**2 + 1/15*n**3 - 437/1440*n**4
h3=                17/480*n**3 - 37/840*n**4
h4=                            4397/161280*n**4

hdot1 = 1/2*n -2/3*n**2 + 5/16*n**3 + 41/180*n**4
hdot2 =       13/48*n**2 - 3/5*n**3 + 557/1440*n**4
hdot3 =                  61/240*n**3 - 103/140*n**4
hdot4 =                             49561/161280*n**4

# http://
# docs.jhs-suositukset.fi/jhs-suositukset/JHS197_liite1/JHS197_liite1.html
k0 = 0.9996 # and
E0 = 500000
#http://docs.jhs-suositukset.fi/jhs-suositukset/JHS197_liite3/JHS197_liite3.html
# ETRS-TM35FIN end


class Tiles():
    """ 
    Stuff related to mapant tiles of Finland 
    x runs from west to east, westmost longitude is
    y runs from NORTH to south, northmost latitude is 

    https://en.wikipedia.org/wiki/Tiled_web_map
    http://wiki.openstreetmap.org/wiki/Slippy_map_tilenames#Python

    https://en.wikipedia.org/wiki/Transverse_Mercator_projection
    (CHAPTER Spherical normal Mercator revisited)

    mapant tiles http://www.mapant.se/fi/tiles/Z/X/Y.png
    mapant starts from west 

    x runs from 0 to 2^(Z+2)
    y runs from 0 to 2^(Z+2)+some extra numbers (it is not square)

    X in meters runs from -76000 to 
        <UpperLeftX>-76000</UpperLeftX>
        <UpperLeftY>7818000</UpperLeftY>
        <LowerRightX>739460</LowerRightX>
        <LowerRightY>6474000</LowerRightY>
                <SizeX>556687</SizeX>
                <SizeY>917504</SizeY>

    suomen itäisin ETRS-TM35FIN  P 6983865  I: 733028
                                 lat 62.54.6110 31.34.8166
                                     62.908514  31.587112

    """
    


    def deg2NE(self, north_deg, east_deg, ):
        """
        Finnish ETRS-TM35FIN coodinate system
http://docs.jhs-suositukset.fi/jhs-suositukset/JHS197_liite2/JHS197_liite2.html
        """
        
        phi = math.radians(north_deg)
        mylambda = math.radians(east_deg)

        Q_dot = math.asinh(math.tan(phi))
        Q_star = math.atanh(e*math.sin(phi))
        Q = Q_dot - e * Q_star
        l = mylambda - math.radians(mylambda0)
        beta = math.atan(math.sinh(Q))
        eta_dot = math.atanh(math.cos(beta) * math.sin(l))
        ksi_dot = math.asin(math.sin(beta) * math.cosh(eta_dot))
        ksi1 = hdot1*math.sin(2*ksi_dot)*math.cosh(2*eta_dot)
        ksi2 = hdot2*math.sin(4*ksi_dot)*math.cosh(4*eta_dot)
        ksi3 = hdot3*math.sin(6*ksi_dot)*math.cosh(6*eta_dot)
        ksi4 = hdot4*math.sin(8*ksi_dot)*math.cosh(8*eta_dot)

        eta1 = hdot1*math.cos(2*ksi_dot)*math.sinh(2*eta_dot)
        eta2 = hdot2*math.cos(4*ksi_dot)*math.sinh(4*eta_dot)
        eta3 = hdot3*math.cos(6*ksi_dot)*math.sinh(6*eta_dot)
        eta4 = hdot4*math.cos(8*ksi_dot)*math.sinh(8*eta_dot)

        ksi = ksi_dot + ksi1 + ksi2 + ksi3 + ksi4
        eta = eta_dot + eta1 + eta2 + eta3 + eta4

        N = A1 * ksi * k0
        E = A1 * eta * k0 + E0

        #print("phi, mylambda: ", phi, mylambda)
        #print("Q_dot, Q_star, Q, l", Q_dot, Q_star, Q, l)
        #print("beta, eta_dot, ksi_dot", beta, eta_dot, ksi_dot)
        #print("ksi1,ksi2,ksi3,ksi4", ksi1,ksi2,ksi3,ksi4)
        
        return (N, E)

    def xy2tile(self, y, x, zoom):
        #0..1
        #xscaled = (x-0.0)/(LowerRightX-UpperLeftX)
        x = x-UpperLeftX # in meters
        resolution = 2**(zoom+2)/3000
        xp = x * resolution
        xtile = int(xp/256) # tiles are 256X256 pixels
        
        
        #xscaled = (x-UpperLeftX) / SizeX
        #xtile = int(xscaled * 2.0 ** (zoom+2))
        print("xtile: ", xtile, 2.0 ** (zoom+2))
    
# http://
# docs.jhs-suositukset.fi/jhs-suositukset/JHS197_liite3/JHS197_liite3.html
# python3 tiles.py 60.3851068722 19.84813676944 0

# Koitere
# python3 tiles.py 62.9460858 30.7148771 7
# map ant gives xtile=509
if __name__ == '__main__':
    import sys
    test = Tiles()
    lon_deg=float(sys.argv[1])
    lat_deg=float(sys.argv[2])
    zoom   =float(sys.argv[3])
    N, E = test.deg2NE(lon_deg, lat_deg)
    print("NE", N, E)
    test.xy2tile(N,E,zoom)
    #test.get_west()
