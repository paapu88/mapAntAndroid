import java.nio.channels.Pipe.SourceChannel;


public class Tiles {
	/*
	 * 
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
	 */
	
    private double[] NE = {0.0d, 0.0d}; //result: NE in ETRS-TM35FIN
    private int[] tiles = {0, 0}; // result: tiles (X, Y)

	
	private final int UpperLeftX = -76000;
	private final int UpperLeftY = 7818000;
	private final int LowerRightX= 739460;
	private final int LowerRightY=6474000;
	private final int SizeX=556687;
	private final int SizeY=917504;
	// map and settings end
	// ETRS-TM35FIN begin
	private final double a = 6378137.0d;
	private final double b = 6356752.314245d;
	private final double fm = 298.257223563d;
	// http://
	//docs.jhs-suositukset.fi/jhs-suositukset/JHS197_liite1/JHS197_liite1.html
	private final double k0 = 0.9996d; // and
	private final double E0 = 500000d; // and
	//http://docs.jhs-suositukset.fi/jhs-suositukset/JHS197_liite3/JHS197_liite3.html
	// ETRS-TM35FIN end
	private double f = 1.0d/ fm;
	private double n = (a-b)/(a+b);
	private double A1= (a/(1.0d+n))*(1.0d + Math.pow(n,2) / 4.0d + Math.pow(n,4) / 64.0d);
	private double e2 = 2.0d*f - f*f;
	private double e=Math.sqrt(e2);
	private double edot2 = e2/(1.0d-e2);
	private final double mylambda0=27.0d;

	private double h1=n/2 -2/3*Math.pow(n, 2) + 37/96*Math.pow(n, 3) - 1/360*Math.pow(n, 4);
	private double h2=     1/48*Math.pow(n, 2) + 1/15*Math.pow(n, 3) - 437/1440*Math.pow(n, 4);
	private double h3=                17/480*Math.pow(n, 3) - 37/840*Math.pow(n, 4);
	private double h4=                            4397/161280*Math.pow(n, 4);

	
	private double hdot1 = 1.0d/2.0d*n -2.0d/3.0d*Math.pow(n, 2) + 5.0d/16.0d*Math.pow(n, 3) + 41.0d/180.0d*Math.pow(n, 4);
	private double hdot2 =       13.0d/48.0d*Math.pow(n, 2)- 3.0d/5.0d*Math.pow(n, 3) + 557.0d/1440.0d*Math.pow(n, 4);
	private double hdot3 =                  61.0d/240.0d*Math.pow(n, 3) - 103.0d/140.0d*Math.pow(n, 4);
	private double hdot4 =                             49561.0d/161280.0d*Math.pow(n, 4);

	
	static double asinh(double x)
	{
	return Math.log(x + Math.sqrt(x*x + 1.0d));
	}

	static double acosh(double x)
	{
	return Math.log(x + Math.sqrt(x*x - 1.0d));
	}

	/*
	static double atanh(double x)
	{
	return 0.5d*Math.log( (x + 1.0d) / (x - 1.0d) );
	} 
	*/
	
	
	static double atanh(double x){
		return (x + Math.pow(x, 3)/3.0d +
				Math.pow(x, 5)/5.0d +
				Math.pow(x, 7)/7.0d +
				Math.pow(x, 9)/9.0d +
				Math.pow(x, 11)/11.0d); 
	}
	
	public double[] deg2NE(double north_deg, double east_deg){
        double phi = Math.toRadians(north_deg);
        double mylambda = Math.toRadians(east_deg);

    	System.out.println("N, HDOT1 "+ n + " " + hdot1);

        
        System.out.println("PHI, LAMBDA:"+phi + " " + mylambda);
        double Q_dot = asinh(Math.tan(phi));
        double Q_star = atanh(e*Math.sin(phi));
        double Q = Q_dot - e * Q_star;
        System.out.println("QDOT, QSTAR"+Q_dot+" "+Q_star);

        double l = mylambda - Math.toRadians(mylambda0);
        double beta = Math.atan(Math.sinh(Q));
        double eta_dot = atanh(Math.cos(beta) * Math.sin(l));
        double ksi_dot = Math.asin(Math.sin(beta) * Math.cosh(eta_dot));
        System.out.println(eta_dot + " " + ksi_dot);

        
        double ksi1 = hdot1*Math.sin(2.0d*ksi_dot)*Math.cosh(2.0d*eta_dot);
        double ksi2 = hdot2*Math.sin(4.0d*ksi_dot)*Math.cosh(4.0d*eta_dot);
        double ksi3 = hdot3*Math.sin(6.0d*ksi_dot)*Math.cosh(6.0d*eta_dot);
        double ksi4 = hdot4*Math.sin(8.0d*ksi_dot)*Math.cosh(8.0d*eta_dot);
        System.out.println("KSI1, HDOT1 "+ksi1+" "+hdot1);
        
        double eta1 = hdot1*Math.cos(2*ksi_dot)*Math.sinh(2*eta_dot);
        double eta2 = hdot2*Math.cos(4*ksi_dot)*Math.sinh(4*eta_dot);
        double eta3 = hdot3*Math.cos(6*ksi_dot)*Math.sinh(6*eta_dot);
        double eta4 = hdot4*Math.cos(8*ksi_dot)*Math.sinh(8*eta_dot);
        
        double ksi = ksi_dot + ksi1 + ksi2 + ksi3 + ksi4;
        double eta = eta_dot + eta1 + eta2 + eta3 + eta4;
        
        System.out.println("KSI, ETA "+ ksi + " " + eta);

        double N = A1 * ksi * k0;
        double E = A1 * eta * k0 + E0;
        
        NE[0] = N;
        NE[1] = E;
        return NE;
	}
	
	
	public int[] xy2tile(double y, double x, int zoom){
        //0..1
        x = x - UpperLeftX; // in meters, running east
        y = UpperLeftY - y; // in meters, running south
        double resolution = Math.pow(2.0d,(zoom+2))/3000.0d;
        double xp = x * resolution;
        double yp = y * resolution;
        int xtile = (int) Math.floor(xp/256.0d); //tiles are 256X256 pixels
        int ytile = (int) Math.floor(yp/256.0d); //tiles are 256X256 pixels
        tiles[0] = xtile;
        tiles[1] = ytile;
        return tiles;
	}   
	
	
	public static void main(String[] args) {
		/*
		 # http://
		 # docs.jhs-suositukset.fi/jhs-suositukset/JHS197_liite3/JHS197_liite3.html
	     # python3 tiles.py 60.3851068722 19.84813676944 0
		 */
		   int zoom = 1;
		   double north_deg = 60.3851068722d;
		   double east_deg = 19.84813676944d;
		   
		   Tiles ti = new Tiles();
		   double[] outNE = {0.0, 0.0};
		   outNE = ti.deg2NE(north_deg, east_deg);
		   System.out.println(outNE[0] + " " + outNE[1]);
		   int[] outTILES = {0, 0};
		   outTILES = ti.xy2tile(outNE[0], outNE[1], 7);
		   System.out.println("TILES:" + outTILES[0] + " " + outTILES[1]);
	
	}
}

