package piedpipers.group4;

import java.util.*;
import java.awt.Rectangle;
import java.awt.geom.Line2D;

import piedpipers.sim.Point;

public class Player extends piedpipers.sim.Player { 

	static int npipers;
	static Point gate;
	static double pspeed = 0.49;
	static double mpspeed = 0.09;

	private Rat[] rats;
	private int[] thetas;

	private Point myLocation;
	private double mySpeed;
	private Rectangle myPartition;

	//which strategy?
	private boolean s1 = true;

	private boolean s2 = false;
	private Rat targetRat = null;
	private Point targetPoint = null;

	private boolean s3 = false;

	//STRATEGY 4
	private boolean s4 = false;
	private boolean s4init = true;

	private boolean player = false;
	private int ratToFix = -1;
	private Point ratToFixPoint = null;
	//private int blink = 0;

	private int nplayers = 1;
	private int ncollectors = npipers - 1;
	private int reducedWidth = 50;
	private int reducedHeight = 100;
	private Rectangle reducedField;
	private Rectangle reducedFieldWithBuffer;
	//

	private boolean initialize = true;
	private boolean sprintToGate = true;

	private class Rat {
		Point currentLocation = new Point();			
		Point lastLocation = new Point();		
		Point velocity = new Point();				
		boolean captured = false;			
		boolean returned = false; //if returned = true, rat is on the left side of the fence 			
	}

	public Player() {
		super();
	}

	public void init() {
	}

	public void initialize(Point[] piperPositions, Point[] ratPositions, boolean[] music, int[] thetas) {
		initRats(ratPositions);
		gate = new Point(dimension/2, dimension/2);
		npipers = piperPositions.length;

		Rectangle[] partitions = createPartitions(dimension/2, dimension, npipers);
		myPartition = partitions[id];

		initialize = false;
	}

	public Rectangle[] createPartitions(int fieldWidth, int fieldHeight, int npipers) {
		Rectangle[] partitions = new Rectangle[npipers];
		int corner_x = dimension/2;
		int corner_y;
		for (int i = 0; i < npipers; i++) {
			corner_y = (dimension/2 - fieldHeight/2) + fieldHeight / npipers * i;
			partitions[i] = new Rectangle(corner_x, corner_y, fieldWidth, fieldHeight / npipers);
			//System.out.println(corner_x + ", " + corner_y + " - " + fieldHeight + " x " + fieldHeight/npipers);
		}
		return partitions;
	}

	public void refresh(Point[] piperPositions, Point[] ratPositions, boolean[] music, int[] thetas) {
		refreshRats(piperPositions, ratPositions, music);
		this.thetas = thetas;

		myLocation = piperPositions[id];
		if (this.music) {
			mySpeed = mpspeed;
		}
		else {
			mySpeed = pspeed;
		}

		updateStrategy();
	}

	public void updateStrategy() {
		if (myLocation.x < gate.x && !noRatsRemain()) {
			s1 = true;
			s2 = s3 = s4 = false;
		}

		else if (noRatsRemain()) {
			s3 = true;
			s1 = s2 = s4 = false;
		}

		else {
			s2 = true;
			s1 = s3 = false;
			//s4 = true;
			//s1 = s2 = s3 = false;
		}

	}


	public Point move(Point[] piperPositions, Point[] ratPositions, boolean[] music, int[] thetas) {
		if (initialize) {
			initialize(piperPositions, ratPositions, music, thetas);
		}

		refresh(piperPositions, ratPositions, music, thetas);

		if (this.s1) {
			return performS1();
		}
		else if (this.s2) {
			return performS2();
		}
		else if (this.s3) {
			return performS3();
		}
		else if (this.s4) {
			return performS4(ratPositions);
		}
		else {
			return performS1();
		}
		
	}

	public Point performS1() {
		this.music = false;
		double dist = distance(myLocation, gate);
		double ox = (gate.x - myLocation.x) / dist * mySpeed;
		double oy = (gate.y - myLocation.y) / dist * mySpeed;

		return new Point(myLocation.x + ox, myLocation.y + oy);
	}

	public Point performS2() {
		this.music = true;
		if (targetRat == null) {
			targetRat = findFastestRatToCatchInPartition();
			if (targetRat == null) {
				targetRat = findFastestRatToCatch();
			}
			targetPoint = findBestApproach(targetRat);
		}

		else if (Math.abs(myLocation.x - targetPoint.x) < 1 && Math.abs(myLocation.y - targetPoint.y) < 1) {
			targetRat = findFastestRatToCatchInPartition();
			if (targetRat == null) {
				targetRat = findFastestRatToCatch();
			}
			targetPoint = findBestApproach(targetRat);
		}


		double dist = distance(myLocation, targetPoint);
		double ox = (targetPoint.x - myLocation.x) / dist * mySpeed;
		double oy = (targetPoint.y - myLocation.y) / dist * mySpeed;
		return new Point(myLocation.x + ox, myLocation.y + oy);
	}

	public Point performS3() {
		Point dropOff = new Point(gate.x - 10, gate.y);
		this.music = true;
		double dist;
		double ox;
		double oy;
		if (myLocation.x >= gate.x) {
			dist = distance(myLocation, gate);
			ox = (gate.x - myLocation.x) / dist * mySpeed;
			oy = (gate.y - myLocation.y) / dist * mySpeed;
		}
		else {
			dist = distance(myLocation, dropOff);
			ox = (dropOff.x - myLocation.x) / dist * mySpeed;
			oy = (dropOff.y - myLocation.y) / dist * mySpeed;
		}
		return new Point(myLocation.x + ox, myLocation.y + oy);
	}

	public Point performS4(Point[] ratPositions) {
		//this assumes more than 1 piper

		if (s4init) {
			nplayers = 1;
			ncollectors = npipers - 1;
			reducedWidth = 50;
			reducedHeight = 100;

			reducedField = new Rectangle(dimension/2, dimension/2 - reducedHeight/2, reducedWidth, reducedHeight);
			reducedFieldWithBuffer = new Rectangle(dimension/2, dimension/2 - reducedHeight/2 - 11, reducedWidth + 11, reducedHeight + 22);

			if (id == npipers - 1) {
				this.player = true;
				this.music = false;
			}

			Rectangle[] partitions = createPartitions(reducedWidth, reducedHeight, ncollectors);
			if (!this.player) {
				myPartition = partitions[id];
			}

			s4init = false;
		}
	
		if (this.player) {
			return performS4P(ratPositions); //player
		}

		else {
			return performS4C(); //collector
		}

	}

	public Point performS4P(Point[] ratPositions) {
		if (this.music && reducedField.contains(myLocation.x, myLocation.y)) {
			System.out.println("WHY AM I HERE?!?!?!");
			System.out.println(myLocation.x + ", " + myLocation.y + " - seeking " + ratToFix + " at " + ratToFixPoint.x + ", " + ratToFixPoint.y);

		}

		if (ratToFix < 0 || isAcceptableAngle(ratToFix)) {
			this.music = false;

			ArrayList<Integer> ratsToFix = new ArrayList<Integer>();
			Rat rat;
			//for (int i = 0; i < rats.length; i++) {
			for (int i = 0; i < ratPositions.length; i++){
				//rat = rats[i];			
				//dimension/2, dimension/2 - reducedHeight/2, reducedWidth, reducedHeight
				//dimension/2, dimension/2 - reducedHeight/2 - 11, reducedWidth + 11, reducedHeight + 22
				//if (!rat.returned && !isAcceptableAngle(i) && !(rat.currentLocation.x >= dimension/2 && rat.currentLocation.x <=dimension/2 +reducedWidth) && !(rat.currentLocation.y >= dimension/2 - reducedHeight/2 && rat.currentLocation.y <= dimension/2 - reducedHeight/2 + reducedHeight)) {
				if (!(ratPositions[i].x < dimension/2) && !isAcceptableAngle(i) && !reducedField.contains(ratPositions[i].x, ratPositions[i].y)){
					Point intercept = findBestApproach(rats[i]);
					System.out.println("Intercept Point..!!:--" + intercept.x + intercept.y);
					if (!reducedField.contains(intercept.x, intercept.y)) {
					//if (!((intercept.x >= dimension/2 && intercept.x <=dimension/2 +reducedWidth) && (intercept.y >= dimension/2 - reducedHeight/2 && intercept.y <= dimension/2 - reducedHeight/2 + reducedHeight))) {
						ratsToFix.add(i);
					}
				} 
			}

			int minTimeRat = 0;
			double currentMinTime = Double.MAX_VALUE;
			int ratId;
			for (int i = 0; i < ratsToFix.size(); i++) {
				ratId = ratsToFix.get(i);

				double timeToCatch = timeToCatch(rats[ratId]);
				if (timeToCatch < currentMinTime) {
					currentMinTime = timeToCatch;
					minTimeRat = i;
				}
			}

			ratToFix = minTimeRat;
			ratToFixPoint = findBestApproach(rats[ratToFix]);
		}

		if (Math.abs(myLocation.x - ratToFixPoint.x) > 1 || Math.abs(myLocation.y - ratToFixPoint.y) > 1) {
			double dist = distance(myLocation, ratToFixPoint);
			double ox = (ratToFixPoint.x - myLocation.x) / dist * mySpeed;
			double oy = (ratToFixPoint.y - myLocation.y) / dist * mySpeed;
			return new Point(myLocation.x + ox, myLocation.y + oy);
		}

		//else reached rat!
		/*
		if (Math.abs(rats[ratToFix].currentLocation.x - myLocation.x) > 2 || Math.abs(rats[ratToFix].currentLocation.y - myLocation.y) > 2 || 
			reducedFieldWithBuffer.contains(myLocation.x, myLocation.y)) {
			System.out.println("Nothing here");
			ratToFix = -1;
			return myLocation;
		}
		*/


		blink();
		return myLocation;
		
		
		
		//return myLocation;
	}

	public void blink() {
		if (this.music) {
			this.music = false;
		}
		else {
			this.music = true;
		}
	}

	public boolean isAcceptableAngle(int ratId) {
		Rat rat = rats[ratId];
		int theta = thetas[ratId];
		Line2D trajectory = new Line2D.Double(rat.currentLocation.x, 
												rat.currentLocation.y, 
												(dimension*Math.sin(theta * Math.PI / 180) + rat.currentLocation.x), 
												(dimension*Math.cos(theta * Math.PI / 180) + rat.currentLocation.y));
		return trajectory.intersects(dimension/2, dimension/2 - reducedHeight/2, reducedWidth, reducedHeight);
	}

	public Point performS4C() {
		this.music = true;
		if (targetRat == null) {
			targetRat = findFastestRatToCatchInPartition();
			if (targetRat == null) {
				targetRat = findFastestRatToCatch();
				if (targetRat == null) {
					System.out.println("Can't find shit");
					return myLocation;
				}
			}
			targetPoint = findBestApproach(targetRat);
		}

		else if (Math.abs(myLocation.x - targetPoint.x) < 1 && Math.abs(myLocation.y - targetPoint.y) < 1) {
			targetRat = findFastestRatToCatchInPartition();
			if (targetRat == null) {
				targetRat = findFastestRatToCatchInReducedField();
				if (targetRat == null) {
					System.out.println("Can't find shit");
					return myLocation;
				}
			}
			targetPoint = findBestApproach(targetRat);
		}


		double dist = distance(myLocation, targetPoint);
		double ox = (targetPoint.x - myLocation.x) / dist * mySpeed;
		double oy = (targetPoint.y - myLocation.y) / dist * mySpeed;
		return new Point(myLocation.x + ox, myLocation.y + oy);

	}

	public Rat findFastestRatToCatch() {
		Rat minTimeRat = null;
		double currentMinTime = Double.MAX_VALUE;

		Rat rat;
		for (int i = 0; i < rats.length; i++) {
			rat = rats[i];
			if (!rat.captured && !rat.returned) {
				double timeToCatch = timeToCatch(rat);
				if (timeToCatch < currentMinTime) {

					currentMinTime = timeToCatch;
					minTimeRat = rat;

				}
			}
		}
		return minTimeRat;
	}

	public Rat findFastestRatToCatchInPartition() {
		Rat minTimeRat = null;
		double currentMinTime = Double.MAX_VALUE;

		Rat rat;
		for (int i = 0; i < rats.length; i++) {
			rat = rats[i];
			if (!rat.captured && !rat.returned) {
				double timeToCatch = timeToCatch(rat);
				if (timeToCatch < currentMinTime) {

					if (canCatchInMyPartition(rat)) {
						currentMinTime = timeToCatch;
						minTimeRat = rat;
					}

				}
			}
		}
		return minTimeRat;
	}

	public boolean canCatchInMyPartition(Rat rat) {
		Point intercept = findBestApproach(rat);
		if (myPartition.contains(intercept.x, intercept.y)) {
			return true;
		}
		//System.out.println(id + " can't catch " + intercept.x + " , " + intercept.y);
		return false;
	}

	public Rat findFastestRatToCatchInReducedField() {
		Rat minTimeRat = null;
		double currentMinTime = Double.MAX_VALUE;

		Rat rat;
		for (int i = 0; i < rats.length; i++) {
			rat = rats[i];
			if (!rat.captured && !rat.returned) {
				double timeToCatch = timeToCatch(rat);
				if (timeToCatch < currentMinTime) {

					if (canCatchInReducedField(rat)) {
						currentMinTime = timeToCatch;
						minTimeRat = rat;
					}

				}
			}
		}
		return minTimeRat;
	}

	public boolean canCatchInReducedField(Rat rat) {
		Point intercept = findBestApproach(rat);
		if (reducedField.contains(intercept.x, intercept.y)) {
			return true;
		}
		//System.out.println(id + " can't catch " + intercept.x + " , " + intercept.y);
		return false;
	}


	//RAT UPDATE STUFF
	//INITIALIZE AND REFRESH RAT LIST, CHECK RAT LIST TO SEE IF ALL CAPTURED
	public void initRats(Point[] ratPositions) {
		rats = new Rat[ratPositions.length];

		for (int i = 0; i < rats.length; i++) {
			rats[i] = new Rat();
		}
		
		initialize = false;
	}
	
	public void refreshRats(Point[] piperPositions, Point[] ratPositions, boolean[] music) {
		for (int i = 0; i < rats.length; i++) {
			Rat rat = rats[i];
			rat.lastLocation = rat.currentLocation;
			rat.currentLocation = ratPositions[i];
			double vx = rat.currentLocation.x - rat.lastLocation.x;
			double vy = rat.currentLocation.y - rat.lastLocation.y;
			rat.velocity = new Point(vx, vy);
			
			if (rat.currentLocation.x < gate.x) { 
				rat.returned = true;
			}
			else {
				rat.returned = false;
			}

			rat.captured = false;
			for (int j = 0; j < piperPositions.length; j++) {
				if (music[j]) {
					if (distance(piperPositions[j], rat.currentLocation) <= 10) {
						rat.captured = true;
						
					}
				}
			}
		}
	}

	public boolean noRatsRemain()
	{
		if(rats.length > 0)
		{
			for(int i = 0; i < rats.length; i++)
			{
				if (!rats[i].captured && !rats[i].returned)
					return false;
			}
			return true;			
		}
		else
			return false;
	}

	//MATH
	
	static double distance(Point a, Point b) {
		return Math.sqrt((a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y));
	}

	private double[] solveQuad(double a, double b, double c) {
		double[] sol = {0,0};
		if (a == 0) {
			if (b == 0) {
				if (c != 0)
					sol = null;
			} else {
				sol[0] = -c/b;
				sol[1] = -c/b;
			}
		} else {
			double disc = b*b - 4*a*c;
			if (disc >= 0) {
				disc = Math.sqrt(disc);
				a = 2*a;
				sol[0] = (-b-disc)/a;
				sol[1] = (-b+disc)/a;
			}
		}
		return sol;
	}

	//RAT INTERCEPTION STUFF

	/*
	private boolean isValidIntercept(Point intercept) {
		Rectangle field = new Rectangle(dimension/2, 0, dimension/2, dimension);
		if (field.contains(intercept.x, intercept.y)) {
			return true;
		}
		return false;
	}
	*/
	

	//RAT INTERCEPTION STUFF

	private boolean isValidIntercept(Point intercept) {
		if (intercept.x <= dimension && 
			intercept.x >= gate.x &&
			intercept.y <= dimension && 
			intercept.y >= 0) {
			return true;
		}
		return false;
	}

	
	// Returns the point the current piper should aim at if attempting to intercept an object at point p moving at velocity v, assuming no wall reflections
	private Point findIntercept(Point p, Point v, double tPast) {
		// Adapted from http://stackoverflow.com/questions/2248876/2d-game-fire-at-a-moving-target-by-predicting-intersection-of-projectile-and-u
		double s = mySpeed;
		
		double tx, ty, tvx, tvy;
		tx =  (p.x - v.x*tPast) - myLocation.x;
		ty =  (p.y - v.y*tPast) - myLocation.y;
		tvx = v.x;
		tvy = v.y;
		
		double a = tvx*tvx + tvy*tvy - s*s;
		double b = 2 * (tvx * tx + tvy * ty);
		double c = tx*tx + ty*ty;
		
		double[] ts = solveQuad(a, b, c);
		
		Point sol = null;
		if (ts != null) {
			double t0 = ts[0];
			double t1 = ts[1];
			double t = Math.min(t0, t1);
			
			if (t < tPast) t = Math.max(t0, t1);
			if (t > tPast) {
				sol = new Point((p.x - v.x*tPast) + tvx*t, (p.y - v.y*tPast) + tvy*t);
			}
		}
		
		return sol;
	}
	
	// Returns the point that the current piper should aim at to most quickly approach the given rat, including the effect of bouncing off walls
	private Point findBestApproach(Rat rat) {
		Point rat_loc = new Point(rat.currentLocation.x, rat.currentLocation.y);
		Point rat_v = new Point(rat.velocity.x, rat.velocity.y);
		double [] next;
		Point intercept = null;
		Point sol = null;
		double t = 0;
		double tIntercept;
		double tBest = 99999999;
		
		int numReflections = 0;
		
		// Keep reflecting off walls until we find an approach vector that works
		// Maximum number of tries
		do {
			intercept = findIntercept(rat_loc, rat_v, t);		
			if (intercept != null) {
				tIntercept = distance(myLocation, intercept)/mySpeed;
				if (tIntercept < tBest && isValidIntercept(intercept)) {
					tBest = tIntercept;
					sol = intercept;
				}
			}
			
			next = intersectWall(rat_loc, rat_v);
			rat_loc = new Point(next[0], next[1]);
			rat_v = new Point(next[3], next[4]);
			t += next[2];
			numReflections++;
		} while (numReflections < 5);

		
		if (sol == null)
			sol = new Point(rat.currentLocation.x, rat.currentLocation.y);
		
		return sol;
	}

	private double timeToCatch(Rat rat) {
		Point rat_loc = new Point(rat.currentLocation.x, rat.currentLocation.y);
		Point rat_v = new Point(rat.velocity.x, rat.velocity.y);
		double [] next;
		Point intercept = null;
		Point sol = null;
		double t = 0;
		double tIntercept;
		double tBest = 99999999;
		
		int numReflections = 0;
		
		// Keep reflecting off walls until we find an approach vector that works
		// Maximum number of tries
		do {
			intercept = findIntercept(rat_loc, rat_v, t);		
			if (intercept != null) {
				tIntercept = distance(myLocation, intercept)/mySpeed;
				if (tIntercept < tBest && isValidIntercept(intercept)) {
					tBest = tIntercept;
					sol = intercept;
				}
			}
			
			next = intersectWall(rat_loc, rat_v);
			rat_loc = new Point(next[0], next[1]);
			rat_v = new Point(next[3], next[4]);
			t += next[2];
			numReflections++;
		} while (numReflections < 5);
		
		if (sol == null)
			sol = new Point(rat.currentLocation.x, rat.currentLocation.y);
		
		return tBest;
	}	
	
	// Given a point and velocity, returns the point, time, and new heading at which it will hit a wall
	private double[] intersectWall(Point p, Point v) {
		// Try it for each wall and return the one with the smallest positive time
		double t  = 999999;
		double tc;
		double sx = 1;
		double sy = 1;
		
		// Right wall: x = dimension
		tc = (dimension - p.x)/v.x;
		if (tc > 0 && tc < t) {
			t = tc;
			sx = -1;
			sy = 1;
		}
		
		// Top wall: y = 0
		tc = (0 - p.y)/v.y;
		if (tc > 0 && tc < t) {
			t = tc;
			sx = 1;
			sy = -1;
		}
		
		// Bottom wall: y = dimension
		tc = (dimension - p.y)/v.y;
		if (tc > 0 && tc < t) {
			t = tc;
			sx = 1;
			sy = -1;
		}
		
		// Left wall: x = 0
		tc = (0 - p.x)/v.x;
		if (tc > 0 && tc < t) {
			t = tc;
			sx = -1;
			sy = 1;
		}
		
		// Center wall: x = dimension/2 with a 2m hole
		tc = (dimension/2 - p.x)/v.x;
		if (tc > 0) {
			Point pc = new Point(p.x + v.x*tc, p.y +v.y*tc);
			if ((pc.y < dimension/2 - 1 || pc.y > dimension/2 + 1) && tc < t) {
				t = Math.min(tc, t);
				sx = -1;
				sy = 1;				
			}
		}
		
		double[] sol = {p.x + v.x*t, p.y + v.y*t, t, v.x*sx, v.y*sy};
		
		return sol;
	}

}