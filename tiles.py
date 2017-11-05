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
edot2 = e**2/(1-e**2)
mylambda0=27.0


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

    python3 getTile.py 31.587112 62.8152 9
    """
    


    def deg2xy(self, north_deg, east_deg, ):
        """
        Finnish ETRS-TM35FIN coodinate system
http://docs.jhs-suositukset.fi/jhs-suositukset/JHS197_liite2/JHS197_liite2.html
        """
        
        phi = math.radians(north_deg)
        mylambda = math.radians(east_deg)

        Q_dot = math.asinh(math.tan(lon_rad))
        Q_star = math.atanh(e*math.sin(lon_rad))
        Q = Q_dot - e * Q_star
        l = mylambda - math.radians(mylambda0)
        beta = math.atan(math.sinh(Q))
        eta_dot = math.atanh(math.cos(beta) * sin(l))
        ksi_dot = math.asin(math.sin(beta) * math.cosh(eta_dot))

        secant=1/math.cos(lat_rad)
        earth_circumference=2*math.pi*R/secant
        x = (lon_rad)*R/secant
        r=earth_circumference/(2*math.pi)
        flatx=2*r*math.sin(lon_rad/2)
        print("x, flatx", x, flatx)
        x = x + west_shift
        #n = 2.0 ** (zoom+2)
        #x = lon_rad 
        #xtile = int((lon_deg + 180.0) / 360.0 * n)
        #ytile = int((1.0 - math.log(math.tan(lat_rad) + (1 / math.cos(lat_rad))) / math.pi) / 2.0 * n)
        y=0
        return (x, y)

    def xy2tile(self, x, y, zoom):
        #0..1
        #xscaled = (x-0.0)/(LowerRightX-UpperLeftX)
        xscaled = (x-UpperLeftX)/(LowerRightX-UpperLeftX)
        #xscaled = (x-UpperLeftX) / SizeX
        xtile = int(xscaled * 2.0 ** (zoom+2))
        print("xtile: ", xtile, 2.0 ** (zoom+2))
    
    def get_west(self):

        lat_deg=60
        lat_rad = math.radians(lat_deg)
        secant=1/math.cos(lat_rad)
        for lon_deg in np.arange(10.34,10.38,0.0001):
            lon_rad = math.radians(lon_deg)
            print(500000 - lon_rad*R/secant, lon_deg)

        #center is 27deg = 500000
        #500000-lon_rad*R/secant=-76000
        # yields 10.36deg
        # so their east long is 16.64 e at arigin

# python3 getTile.py 25 60 0
if __name__ == '__main__':
    import sys
    test = Tiles()
    lon_deg=float(sys.argv[1])
    lat_deg=float(sys.argv[2])
    zoom   =float(sys.argv[3])
    x, y = test.deg2xy(lon_deg, lat_deg)
    print("XY", x, y)
    test.xy2tile(x,y,zoom)
    #test.get_west()
