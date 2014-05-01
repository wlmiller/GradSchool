/* aspherical_coexistence.c
 *
 * Performs a Monte Carlo of a box containing one-half fcc crystalline and 
 * one-half liquid-phase hard aspherical particles. Particles are constructed
 * as described in [1], by placing some number of balls in an overlapping 
 * configuration.  The number of balls and degree of overlap are adjustable,
 * allowing for control of the degree of asphericity.
 *
 * The number of crystalline particles is tracked as a function of Monte Carlo
 * step.  Coexstence is occurs at the pressure at which this quantity stabilizes,
 * neither increasing to encompass the whole system nor decreasing to near zero.
 * See [2] for more details.
 *
 * Usage:
 *      ./aspherical_coexistence nBall sigma pressure
 *          * nBall = the number of balls of which the particles are constructed
 *          * sigma = how far from the particle center balls are placed in particle construction
 *          * pressure = the system pressure
 *
 * [1] W. L. Miller, B. Bozorgui, and A. Cacciuto. Crystallization of hard aspherical particles.
 *      J. Chem. Phys. 132, 134901 (2010).
 * [2] W. L. Miller and A. Cacciuto. On the phase behavior of hard aspherical particles.
 *      J. Chem. Phys. 133, 234903 (2010).*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define NLIQ 1008
#define EPSILON 1e-12

typedef struct {
	double x, y, z;
} vec;                      // a three-dimensional vector

typedef struct {
	double w, x, y, z;
} quat;                     // a quaternion

typedef struct {
	vec pos;				// particle position
	vec vpos;				// particle "Verlet position"
	int nVerl;				// number of particles in Verlet list 
	int vList[100];			// Verlet list
	int nBond;				// number of bonds to other particles  
	int bList[100];			// list of bonded particles
	int nClus;				// number of particles in cluster
	int cList[100];			// list of particles in cluster
	vec ballpos[12];		// "balls" making up a single particle
	double balldiam;		// diameter of the "balls"
	quat orient;			// particle orienttion
} PARTICLE;

typedef struct {
	double vdiam, vgap;		// Verlet diameter and gap	
	double bdiam, cdiam;	// maximum diameter of a bond and cluster
} SYSTEM;

typedef struct{
	int     ncont;             // number of neighboring particles
	double  q6q6treshold;      // minimum cutoff for "crystallinity"
	int     ncluster,nsolid;   // number of clusters and crystals in it
	double  Y6coef[14];        // coeff spherical harmonics
	int     p_c[2*NLIQ+1];
} ORDER_PARAMETER;

typedef struct{
	double q6Re[13];			// real part of the order parameter
	double q6Im[13];			// Imaginary part of the order parameter
	int nconn;					// number of connecting neighbors
	int state;					// liquid/crystal (0/1)
} BONDS;

SYSTEM sys;
PARTICLE* part;
ORDER_PARAMETER xtal;
BONDS *bond;
FILE *outfile, *outfile2;
int solid_list[2*NLIQ+1];

void buildPart(int);
void checkVerlet();
void setVerlet();
void setVerletClust();
void melt();
void melt2();
void MCMove(int);
void MCRot(int);
void MCVol(int);
void write(int);
double energyTotal();
double energy(int);
void set_Order_Parameter();
void q6();
int q6dot();
vec setRandVec();
void Rotate(int, double, vec);

double sx, sy, sz;
double sx2, sy2, sz2;
double barrier = 0;
int hard = 0;
int NPART = NLIQ;
double pressure, sigma;
int NBALL;

main(int argc, char* argv[]) {
	int i, j, k, l = 0;
	double n;
	double x, y, z;
	double dense, hsdens = 0.545;
	double scale;
	double dx, dy, dz;

	srand48(111111);

	NBALL = atoi(argv[1]);          // the number of balls making up 1 particle
	sigma = atof(argv[2]);          // how far from the center of the particle the balls are placed
	pressure = atof(argv[3]);       // the system pressure

	printf("%d %lf %lf\n",NBALL,sigma,pressure);	

	outfile = fopen("config.xyz","w");
	outfile2 = fopen("config.out","w");

	part = (PARTICLE *) malloc((2016)*sizeof(PARTICLE));
	bond = (BONDS *) malloc((12*2016)*sizeof(PARTICLE));

	sys.vdiam = 1.8;
	sys.bdiam = 1.5;
	sys.cdiam = 1.5;
	sys.vgap = 0.5*(sys.vdiam - 1.0);

	set_Order_Parameter();

	dense = 4./3.*M_PI*(0.5*0.5*0.5)*1008/(7.5*(7*sqrt(3.)/2. + 1)*(17.*sqrt(3.)/2. + 1));
	scale = pow(dense/hsdens,1./3.);
	sx = scale*7.;

	sy = scale*(7.*sqrt(3.)/2.+1);
	sz = 4./3.*M_PI*(0.5*0.5*0.5)*1008/(0.494*sx*sy) - 1.5;

	for (i = 0; i < NLIQ; i++) {
		// build a particle
		buildPart(i);

		retry:
			part[i].pos.x = sx*(0.5-drand48());
			part[i].pos.y = sy*(0.5-drand48());
			part[i].pos.z = (sz)*(0.5-drand48());

		// randomly place the particle, ensuring it's not too close to any existing ones
		for (j = l; j < i; j++) {
				dx = part[i].pos.x - part[j].pos.x;
				dy = part[i].pos.y - part[j].pos.y;
				dz = part[i].pos.z - part[j].pos.z;
			
				if (fabs(dx) > sx/2)
					dx -= sx*fabs(dx)/dx;
				if (fabs(dy) > sy/2)
					dy -= sy*fabs(dy)/dy;
				if (fabs(dz) > sz/2)
					dz -= sz*fabs(dz)/dz;

				if (sqrt(dx*dx + dy*dy + dz*dz) < 0.85) goto retry;
		}
	}

	// scale the particle positions by the system size
	for (i = 0; i < NLIQ; i++) {
		part[i].pos.x /= sx;
		part[i].pos.y /= sy;
		part[i].pos.z /= sz;
	}
	
	printf("%lf\n",1008*4./3.*M_PI*(0.5*0.5*0.5)/(sx*sy*(sz+1)));

	// allow the liquid particles to relax so there are no overlaps
	melt();

	// scale the positions back out for addition of the solid
	for (i = 0; i < NLIQ; i++) {
		part[i].pos.x *= sx;
		part[i].pos.y *= sy;
		part[i].pos.z *= sz;
	}

	NPART = 2*NLIQ;
	l = NLIQ;

	// set the size of the "crystalline" size of the box
	sx2 = sx;
	sy2 = sy;
	sz2 = scale*(17.*sqrt(3.)/2. + 1);

	// place crystalline particles
	for (i = 0; i < 14; i++) {
		for (j = 0; j < 4; j++) {
			for (k = 0; k < 18; k++) {
				buildPart(l);
				part[l].pos.x = scale*(-3.5 + 0.5*(double)(i) + 0.5*(k%2));
				part[l].pos.y = scale*(-4.*sqrt(3.)/2. + 0.25 + sqrt(3.)/2.*(i%2)+sqrt(3.)*(double)(j)+sqrt(3.)/4.*(k%2));
				part[l].pos.z = scale*(-9.*sqrt(3.)/2. + 0.5 + k*sqrt(3.)/2.);
				l++;
			}
		}
	}

	for (i = NLIQ; i < NPART; i++) {
		if (fabs(part[i].pos.y) > sy/2.) printf("%d y\n",i);
		if (fabs(part[i].pos.z) > sz/2.) printf("%d z\n",i);
	}	

	// scale the z positions of the liquid particles to fit in the new box correctly
	for (i = 0; i < NLIQ; i++) {
		if (part[i].pos.z < 0) part[i].pos.z = sz2/2. + (sz+1.5)/2. + part[i].pos.z;
		else part[i].pos.z = -sz2/2. - (sz+1.5)/2. + part[i].pos.z;
	}

	sz = sz + sz2+1.5;

	// scale everything back down
	for (i = 0; i < NPART; i++) {
		part[i].pos.x /= sx;
		part[i].pos.y /= sy;
		part[i].pos.z /= sz;
	}

	// allow the crytalline particles to relax so there aren't overlaps 
	melt2();

	hard = 1;
	
	write(0);

	for (i = 0; i < NPART; i++) {
		if (energy(i) > 0) printf("%d\n",i);
	}

	// freeze the crystal, relax the liquid
	for (i = 1; i <= 250000; i++) {
		for (j = 0; j < 2*NLIQ; j++) {
			k = (int)(2*NLIQ*drand48());
			if (k < NLIQ) MCMove(k);
			else MCRot(k-NLIQ);
		}
		if (i % 25000 == 0) write(i);
	}

	// freeze the liquid, relax the crystal
	for (i = 250001; i <= 500000; i++) {
		for (j = 0; j < 2*(NPART-NLIQ); j++) {
			k = NLIQ + (int)((NPART-NLIQ)*2*drand48());
			if (k < NPART) MCMove(k);
			else MCRot(k-(NPART-NLIQ));
		}
		if (i % 25000 == 0) write(i);
	}

	// run Monte Carlo on the entire system
	for (i = 500001; i <= 10000000; i++) {
		for (j = 0; j < 2*NPART+1; j++) {
			k = (int)((2*NPART+1)*drand48());
			if (k < NPART) MCMove(k);
			else if (k < 2*NPART) MCRot(k-NPART);
			else {
				k = (int)(3*drand48());
				MCVol(k);
			}
		}
		if (i % 25000 == 0) write(i);
	}
}

// Build a particle out of a set of smaller "balls." 
void buildPart(int k) {
	int i, j;
	double x, y, z, dx, dy, dz;
	int inside = 0, outside = 0;
	int flag = 1;
	double scale;

	for (i = 0; i < NBALL; i++) {
		part[k].ballpos[i] = setRandVec();
		part[k].ballpos[i].x *= sigma;
		part[k].ballpos[i].y *= sigma;
		part[k].ballpos[i].z *= sigma;
	}
	
	for (i = 0; i < 1000000; i++) {	
		flag = 1;

		x = (sigma+0.5)*2.*(drand48()-0.5);
		y = (sigma+0.5)*2.*(drand48()-0.5);
		z = (sigma+0.5)*2.*(drand48()-0.5);
	
		for (j = 0; j < NBALL; j++) {
			dx = x - part[k].ballpos[j].x;
			dy = y - part[k].ballpos[j].y;
			dz = z - part[k].ballpos[j].z;

			if (sqrt(dx*dx+dy*dy+dz*dz) < 0.5) {
				inside++;
				flag = 0;
				break;
			}
		}
		outside += flag;
	}

	// scale the balls so that the total volume of the particle is the same as the volume of a hard sphere
    // with radius 0.5.
    scale = pow(1./6.*M_PI*(double)(inside+outside)/((1+2*sigma)*(1+2*sigma)*(1+2*sigma)*(double)(inside)),1./3.);
	part[k].balldiam = scale;
	
	for (i = 0; i < NBALL; i++) {
		part[k].ballpos[i].x *= scale;
		part[k].ballpos[i].y *= scale;
		part[k].ballpos[i].z *= scale;
	}

	part[k].orient.w = 1.;
	part[k].orient.x = 0.;
	part[k].orient.y = 0.;
	part[k].orient.z = 0.;
}

// Check if any particles have moved too far from the position at which
// the Verlet list was constructed.  If so, rebuild it.
void checkVerlet(int k) {
	double dx, dy, dz, dsq;

	dx = part[k].pos.x - part[k].vpos.x;
	dy = part[k].pos.y - part[k].vpos.y;
	dz = part[k].pos.z - part[k].vpos.z;

	if (fabs(dx) > 0.5)
		dx -= fabs(dx)/dx;
	if (fabs(dy) > 0.5)
		dy -= fabs(dy)/dy;
	if (fabs(dz) > 0.5)
		dz -= fabs(dz)/dz;

	dsq = (dx*dx*sx*sx + dy*dy*sy*sy + dz*dz*sz*sz);
	if (dsq > sys.vgap*sys.vgap) setVerlet();
}

// Build Verlet lists
void setVerlet() {
	int i, j;
	double dx, dy, dz, dsq;

	for (i = 0; i < NPART; i++) {
		part[i].nVerl = 0;
		part[i].vpos = part[i].pos;
		part[i].nBond = 0;
	}

	for (i = 0; i < NPART-1; i++) {
		for (j = i+1; j < NPART; j++) {
			dx = part[i].pos.x-part[j].pos.x;
			dy = part[i].pos.y-part[j].pos.y;
			dz = part[i].pos.z-part[j].pos.z;

			if (fabs(dx) > 0.5)
				dx -= fabs(dx)/dx;
			if (fabs(dy) > 0.5)
				dy -= fabs(dy)/dy;
			if (fabs(dz) > 0.5)
				dz -= fabs(dz)/dz;

			dsq = (dx*dx*sx*sx + dy*dy*sy*sy + dz*dz*sz*sz);

			if (dsq < sys.vdiam*sys.vdiam) {
				part[i].vList[part[i].nVerl] = j;
				part[i].nVerl++;
				part[j].vList[part[j].nVerl] = i;
				part[j].nVerl++;
			}

			if (dsq < sys.bdiam*sys.bdiam) {
				part[i].bList[part[i].nBond] = j;
				part[i].nBond++;
				part[j].bList[part[j].nBond] = i;
				part[j].nBond++;
			}
		}
	}
}

// Same as above, but for clusters
void setVerletClust(int nsol){  
	int p,t,i,k;
	double r2;
	double dx,dy,dz;

	for (i=0;i<=NPART;i++){
		part[i].nClus=0;
	}

	for (p=1;p<=nsol-1;p++){
		for (t=p+1;t<=nsol;t++){
			i=solid_list[p];
			k=solid_list[t];

			dx=part[i].pos.x-part[k].pos.x;
			dy=part[i].pos.y-part[k].pos.y;
			dz=part[i].pos.z-part[k].pos.z;

			if (dx>.5) dx-=1.;
			if (dy>.5) dy-=1.;
			if (dz>.5) dz-=1.;

			if (dx<-.5) dx+=1.;
			if (dy<-.5) dy+=1.;
			if (dz<-.5) dz+=1.;

			r2=(dx*dx*sx*sx + dy*dy*sy*sy + dz*dz*sz*sz);
			if (r2<sys.cdiam*sys.cdiam){
				part[k].cList[ part[k].nClus ]=i;
				part[i].cList[ part[i].nClus ]=k;
				part[i].nClus++; part[k].nClus++;
			}
		}
	}   // end while loop
}

// Allows the liquid particles to relax
void melt() {
	int i, j, k, l;
	vec oldpos;
	double olden, newen;
	quat oldorient;
	vec oldball[NBALL];

	setVerlet();

	barrier++;

	outfile = fopen("scramble.xyz","w");
	fclose(outfile);

	barrier++;

	i = 0;
	while (!(energyTotal() == 0.000000000)) {	
		i++;
		for (j = 0; j < 2.*NLIQ; j++) {
			k = (int)(2*NLIQ*drand48());

			if (k < NLIQ) {
				checkVerlet(k);

				oldpos = part[k].pos;

				olden = energy(k);

				part[k].pos.x += 0.01*(drand48()-0.5)/sx;
				part[k].pos.y += 0.01*(drand48()-0.5)/sy;
				part[k].pos.z += 0.01*(drand48()-0.5)/sz;

				if (fabs(part[k].pos.x) > 0.5)
					part[k].pos.x -= fabs(part[k].pos.x)/part[k].pos.x;
				if (fabs(part[k].pos.y) > 0.5)
					part[k].pos.y -= fabs(part[k].pos.y)/part[k].pos.y;
				if (fabs(part[k].pos.z) > 0.5)
					part[k].pos.z -= fabs(part[k].pos.z)/part[k].pos.z;

				checkVerlet(k);

				newen = energy(k);

				if (drand48() > exp(-(newen - olden))) {
					part[k].pos = oldpos;
				}
				
			} else {
				k -= NLIQ;

				oldorient = part[k].orient;
				olden = energy(k);

				for (l = 0; l < NBALL; l++) {
					oldball[l] = part[k].ballpos[l];
				}

				Rotate(k, 2.*2.*(drand48()-0.5), setRandVec());

				newen = energy(k);

				if (drand48() > exp(-(newen - olden))) {
					part[k].orient = oldorient;
					for (l = 0; l < NBALL; l++) {
						part[k].ballpos[l] = oldball[l];
					}
				}
			}
		}

		if (i % 10 == 0) {
			printf("%d %lf\n",i,energyTotal());
			barrier++;
		}
	}
}

// Allow the crystalline particles to relax
void melt2() {
	int i, j, k, l;
	vec oldpos;
	double olden, newen;
	quat oldorient;
	vec oldball[NBALL];

	setVerlet();

	barrier++;

	outfile = fopen("scramble.xyz","w");
	fclose(outfile);

	barrier++;

	i = 0;
	while (!(energyTotal() == 0.000000000)) {	
		i++;
		for (j = 0; j < 2.*(NPART-NLIQ); j++) {
			k = NLIQ + (int)(2*(NPART-NLIQ)*drand48());

			if (k < NPART) {
				checkVerlet(k);

				oldpos = part[k].pos;

				olden = energy(k);

				part[k].pos.x += 0.01*(drand48()-0.5)/sx;
				part[k].pos.y += 0.01*(drand48()-0.5)/sy;
				part[k].pos.z += 0.01*(drand48()-0.5)/sz;

				if (fabs(part[k].pos.x) > 0.5)
					part[k].pos.x -= fabs(part[k].pos.x)/part[k].pos.x;
				if (fabs(part[k].pos.y) > 0.5)
					part[k].pos.y -= fabs(part[k].pos.y)/part[k].pos.y;
				if (fabs(part[k].pos.z) > 0.5)
					part[k].pos.z -= fabs(part[k].pos.z)/part[k].pos.z;

				checkVerlet(k);

				newen = energy(k);

				if (drand48() > exp(-(newen - olden))) {
					//	if (newen == 999999) {
					part[k].pos = oldpos;
				}
				
			} else {
				k -= (NPART-NLIQ);

				oldorient = part[k].orient;
				olden = energy(k);

				for (l = 0; l < NBALL; l++) {
					oldball[l] = part[k].ballpos[l];
				}

				Rotate(k, 2.*2.*(drand48()-0.5), setRandVec());

				newen = energy(k);

				if (drand48() > exp(-(newen - olden))) {
					part[k].orient = oldorient;
					for (l = 0; l < NBALL; l++) {
						part[k].ballpos[l] = oldball[l];
					}
				}
			}
		}

		if (i % 10 == 0) {
			printf("%d %lf\n",i,energyTotal());
			barrier++;
		}
	}
}

// Attempt a Monte Carlo move of a particle
void MCMove(int k) {
	int i, j;
	vec oldpos;
	double olden, newen;

	oldpos = part[k].pos;

	part[k].pos.x += 0.01*(drand48()-0.5)/sx;
	part[k].pos.y += 0.01*(drand48()-0.5)/sy;
	part[k].pos.z += 0.01*(drand48()-0.5)/sz;

	if (fabs(part[k].pos.x) > 0.5)
		part[k].pos.x -= fabs(part[k].pos.x)/part[k].pos.x;
	if (fabs(part[k].pos.y) > 0.5)
		part[k].pos.y -= fabs(part[k].pos.y)/part[k].pos.y;
	if (fabs(part[k].pos.z) > 0.5)
		part[k].pos.z -= fabs(part[k].pos.z)/part[k].pos.z;

	checkVerlet(k);

	newen = energy(k);

	if (newen == 999999) {
		part[k].pos = oldpos;
	}
}

// Attempt a Monte Carlo particle rotation
void MCRot(int k) {
	int i;
	vec oldball[NBALL];
	quat oldorient;
	double newen;

	for (i = 0; i < NBALL; i++) {
		oldball[i] = part[k].ballpos[i];
	}

	oldorient = part[k].orient;

	Rotate(k, 2.*2.*(0.5-drand48()), setRandVec());

	newen = energy(k);

	if (newen == 999999) {
		for (i = 0; i < NBALL; i++) {
			part[k].ballpos[i] = oldball[i];
		}

		part[k].orient = oldorient;
	}
}

// Attempt a Monte Carlo volume change
void MCVol(int k) {
	int i, j;
	double oldsx, oldsy, oldsz;
	double olden, newen;
	double oldvol, newvol;

	oldsx = sx; 
	oldsy = sy; 
	oldsz = sz; 

	oldvol = oldsx*oldsy*oldsz;
	newvol = exp(log(oldvol) + 2*(drand48()-0.5)*0.01);

	switch (k) {
		case 0:
			sx = newvol/(sy*sz);
			break;
		case 1:
			sy = newvol/(sx*sz);
			break;
		case 2:
			sz = newvol/(sx*sy);
			break;
	}

	newen = energyTotal();

	if (newen > 0) {
		sx = oldsx;
		sy = oldsy;
		sz = oldsz;
		return;
	}

	if (drand48() > exp(-(pressure*(newvol-oldvol)-(NPART+1)*log(newvol/oldvol)))) {
		sx = oldsx;
		sy = oldsy;
		sz = oldsz;
	}
}

// Write the system configuration.
write(int k) {
	int i;
	int q6;

	q6 = q6dot();
	outfile = fopen("config.xyz","a");
	fprintf(outfile,"%d\n\n",NPART);
	
	// write particle positions
    for (i = 0; i < NPART; i++) {
		if (bond[i].state) fprintf(outfile,"He ");
		else fprintf(outfile, "H ");
		fprintf(outfile,"%lf %lf %lf\n",part[i].pos.x*sx,part[i].pos.y*sy,part[i].pos.z*sz);
	}
	fclose(outfile);

	outfile2 = fopen("config.out","a");
	// write order parameters - crystallinity and density
    fprintf(outfile2,"%d %d %lf %lf\n",k,q6,4./3.*M_PI*(1./8.)*NPART/(sx*sy*sz),energyTotal());
	printf("%d %d %lf %lf\n",k,q6,4./3.*M_PI*(1./8.)*NPART/(sx*sy*sz),energyTotal());
	fclose(outfile2);
}

// Calculate the potential energy of a particle.
// For hard balls, this is just checking for overlaps.
double energy (int k) {
	int i, j, l, m;
	double dx, dy, dz, d;
	double en = 0;
	double diam;

	for (i = 0; i < part[k].nVerl; i++) {
		l = part[k].vList[i];
		diam = (0.5+sigma)*(part[k].balldiam + part[l].balldiam);
		
		dx = part[k].pos.x - part[l].pos.x;
		dy = part[k].pos.y - part[l].pos.y;
		dz = part[k].pos.z - part[l].pos.z;

		if (fabs(dx) > 0.5)
			dx -= fabs(dx)/dx;
		if (fabs(dy) > 0.5)
			dy -= fabs(dy)/dy;
		if (fabs(dz) > 0.5)
			dz -= fabs(dz)/dz;

		d = sqrt(dx*dx*sx*sx + dy*dy*sy*sy + dz*dz*sz*sz);

		if (d < diam) { 
			diam = 0.5*(part[k].balldiam+part[l].balldiam);
			
			for (j = 0; j < NBALL; j++) {
				for (m = 0; m < NBALL; m++) {
					dx = (part[k].pos.x - part[l].pos.x) + (part[k].ballpos[j].x - part[l].ballpos[m].x)/sx;
					dy = (part[k].pos.y - part[l].pos.y) + (part[k].ballpos[j].y - part[l].ballpos[m].y)/sy;
					dz = (part[k].pos.z - part[l].pos.z) + (part[k].ballpos[j].z - part[l].ballpos[m].z)/sz;

					if (fabs(dx) > 0.5)
						dx -= fabs(dx)/dx;
					if (fabs(dy) > 0.5)
						dy -= fabs(dy)/dy;
					if (fabs(dz) > 0.5)
						dz -= fabs(dz)/dz;

					d = sqrt(dx*dx*sx*sx + dy*dy*sy*sy + dz*dz*sz*sz);

					if (d < diam) {
						if (hard) {
							return 999999;
						} else {
							en += diam*barrier/d;
						}
					}
				}
			}
		}
	}
	return en;
}

// Calculate the potential energy of the whole system.
double energyTotal() {
	int i;
	double entot=0;

	if (!hard) {
		for (i = 0; i < NPART; i++) {
			entot += energy(i);
		}

		return entot/2.;
	} else {
		for (i = 0; i < NPART; i++) {
			if (energy(i) == 999999) return 999999;
		}

		return 0;
	}
}

// Set up the order paramter coefficients
void set_Order_Parameter(void){
	int i;
	for (i=1;i<=NPART;i++) xtal.p_c[i]=0;
	xtal.q6q6treshold=.5;
	xtal.ncont=6;
	xtal.ncluster=0; 
	xtal.nsolid=0;
	xtal.Y6coef[0]=1./32.*sqrt(13./M_PI);
	xtal.Y6coef[1]=1./16.*sqrt(273./(2*M_PI));
	xtal.Y6coef[2]=1./64.*sqrt(1365./M_PI);
	xtal.Y6coef[3]=1./32.*sqrt(1365./M_PI);
	xtal.Y6coef[4]=3./32.*sqrt(91./(2*M_PI)); 
	xtal.Y6coef[5]=3./32.*sqrt(1001./M_PI);
	xtal.Y6coef[6]=1./64.*sqrt(3003./M_PI);
}

// Calculate the q6 order parameter over the system.
void q6(void){
	int i,j,k,o;
	double Q6[2][13],B,q6n,Q6n[13];
	double dr,dr2,dx,dy,dz;
	double y6[2][13];
	double two_cosphi,cosphi[7],sinphi[7];
	double costh,costh2,costh4,costh6;
	double sinth,sinth2,sinth3,sinth4,sinth5,sinth6;

	for (i=0;i<NPART;i++){

		for (k=0;k<=6;k++){
			Q6[0][k]=0.;  //  Real  part
			Q6[1][k]=0.;  //  Imag. part    
		}

		for (j=0;j<part[i].nBond;j++){

			o=part[i].bList[j];
			dx=part[i].pos.x-part[o].pos.x;
			dy=part[i].pos.y-part[o].pos.y;
			dz=part[i].pos.z-part[o].pos.z;

			if (dx>.5) dx-=1.;
			if (dy>.5) dy-=1.;
			if (dz>.5) dz-=1.;

			if (dx<-.5) dx+=1.;
			if (dy<-.5) dy+=1.;
			if (dz<-.5) dz+=1.;

			dx*=sx;dy*=sy;dz*=sz;

			dr2=(dx*dx+dy*dy+dz*dz);
			dr=sqrt(dr2);

			costh=dz/dr;
			costh2=costh*costh;
			costh4=costh2*costh2;
			costh6=costh4*costh2;

			sinth2=1.-costh2;
			sinth=sqrt(sinth2);
			sinth3=sinth2*sinth;
			sinth4=sinth2*sinth2;
			sinth5=sinth3*sinth2;
			sinth6=sinth3*sinth3;

			if (fabs(sinth)<EPSILON){
				cosphi[1]=0.;
				sinphi[1]=0.;	
			}else{
				cosphi[1]=dx/(dr*sinth);
				sinphi[1]=dy/(dr*sinth);	
			}

			two_cosphi=2.0*cosphi[1];
			cosphi[2]=two_cosphi*cosphi[1]-1.0;
			sinphi[2]=two_cosphi*sinphi[1];
			cosphi[3]=two_cosphi*cosphi[2]-cosphi[1];
			sinphi[3]=two_cosphi*sinphi[2]-sinphi[1];
			cosphi[4]=two_cosphi*cosphi[3]-cosphi[2];
			sinphi[4]=two_cosphi*sinphi[3]-sinphi[2];
			cosphi[5]=two_cosphi*cosphi[4]-cosphi[3];
			sinphi[5]=two_cosphi*sinphi[4]-sinphi[3];
			cosphi[6]=two_cosphi*cosphi[5]-cosphi[4];
			sinphi[6]=two_cosphi*sinphi[5]-sinphi[4];

			y6[0][0]=xtal.Y6coef[0]*(-5.+105.*costh2-315.*costh4+231.*costh6);
			y6[1][0]=0.;

			B=xtal.Y6coef[1]*costh*(5.-30.*costh2+33.*costh4)*sinth;
			y6[0][1]=-B*cosphi[1];
			y6[1][1]=-B*sinphi[1];

			B=xtal.Y6coef[2]*(1.-18.*costh2+33.*costh4)*sinth2;
			y6[0][2]=B*cosphi[2];
			y6[1][2]=B*sinphi[2];

			B=xtal.Y6coef[3]*costh*(-3.+11.*costh2)*sinth3;
			y6[0][3]=-B*cosphi[3];
			y6[1][3]=-B*sinphi[3];

			B=xtal.Y6coef[4]*(-1.+11.*costh2)*sinth4;
			y6[0][4]=B*cosphi[4];
			y6[1][4]=B*sinphi[4];

			B=xtal.Y6coef[5]*costh*sinth5;
			y6[0][5]=-B*cosphi[5];
			y6[1][5]=-B*sinphi[5];

			B=xtal.Y6coef[6]*sinth6;
			y6[0][6]=B*cosphi[6];
			y6[1][6]=B*sinphi[6];

			for (k=0;k<=6;k++){
				Q6[0][k]=Q6[0][k]+y6[0][k];
				Q6[1][k]=Q6[1][k]+y6[1][k];
			}	   


		}

		// ADDED NORMALIZATION 

		for (k=0;k<=6;k++){
			Q6n[k]=Q6[0][k]*Q6[0][k]+Q6[1][k]*Q6[1][k];
		}

		q6n=0.;
		for (k=1;k<=6;k++){
			q6n=q6n+2*Q6n[k];
		}  
		q6n = q6n + Q6n[0];
		q6n=sqrt(q6n);
		// END  ADDED NORMALIZATION (set q6n=1 if not needed)

		for (k=0;k<=6;k++){
			bond[i].q6Re[k]=Q6[0][k]/q6n;
			bond[i].q6Im[k]=Q6[1][k]/q6n;
		}    

	} 
}

// Calculate the dot product q6*q6 for pairs of particles in clusters
// throughout the system.
int q6dot(void){
	int i,j,k,o,tempor,nc,l,p;
	int uu,s,flag_no,number;
	int flag[NPART+1],which_c=0;
	int keep[NPART+1],keep2[NPART+1];
	int cluster[NPART+1],nClus;
	double q6q6;

	setVerlet();  
	q6();
	tempor=0;  
	nClus=0;

	for (i=0;i<NPART;i++){
		bond[i].nconn=0;
		bond[i].state=0;    
		solid_list[i]=0;
		cluster[i]=0;
		flag[i]=0;
	}  

	xtal.nsolid=0;

	for (i=0;i<NPART;i++){
		for (j=0;j<part[i].nBond;j++){
			q6q6=0.; 
			o=part[i].bList[j];
			for (k=1;k<=4;k++){
				q6q6+=bond[i].q6Re[k]*bond[o].q6Re[k]+
					bond[i].q6Im[k]*bond[o].q6Im[k];
			}
			q6q6=2.*q6q6+bond[i].q6Re[0]*bond[o].q6Re[0];
			if (q6q6>xtal.q6q6treshold) bond[i].nconn++;
		}

		// if the particle has xtal.ncont neighbors with a high order paramter...
		if (bond[i].nconn>=xtal.ncont){
			xtal.nsolid++;
			bond[i].state=1; // solid-like particle
			solid_list[xtal.nsolid]=i;
		}else bond[i].state=0;
	} 

	/* NOW MAKE CLUSTERS */

	nc=0;
	if (xtal.nsolid==0) {return 0;}
	setVerletClust(xtal.nsolid);   
	nc=1; uu=0;

	for (k=1;k<=xtal.nsolid;k++){
		flag_no=0;
		l=0;   
		s=solid_list[k];
		if (flag[s]==1){ flag_no=1; goto stopcluster;}
		uu=1; flag[s]=1;    
		keep[uu]=s;

		for (p=0;p<part[s].nClus;p++){
			o=part[s].cList[p];
			if (flag[o]==0){	
				uu++;l++;flag[o]=1;keep[l]=o;
			}
		}
		if (l==0) goto stopcluster;
		number=l;

		while (l!=0){
			l=0;
			for (j=1;j<=number;j++){
				for (p=1;p<=part[ keep[j] ].nClus;p++){
					o=part[keep[j]].cList[p];
					if (flag[o]==0){ // inside neigh cells
						uu++;l++;flag[o]=1;keep2[l]=o;
					}
				}
			}
			number=l;for (p=1;p<=number;p++) keep[p]=keep2[p];
		} 


	stopcluster: 
		if (uu>0 && flag_no!=1) {
			nClus=uu;
			if (nClus>=tempor) {
				tempor=nClus;
				which_c=nc;nc++;
			} 
		}
	}

	skip:  
		return tempor;  
}

// Choose a random vector.
vec setRandVec() {
	double x, y, z, s = 2.;
	vec e;

	while (s > 1.) {
		x = 2.*drand48() - 1.;
		y = 2.*drand48() - 1.;
		s = x*x + y*y;
	}

	e.z = 1 - 2.*s;
	s = 2.*sqrt(1.-s);
	e.x = s*x;
	e.y = s*y;

	return e;
}

// Rotate particle i by angle theta around vec vector,
// using quaternion multiplication.
void Rotate(int i, double theta, vec vector) {
	int j;
	double cn, sn;
	double w, x, y, z;
	double xx, yy, zz, xy, xz, yz, wx, wy, wz;
	vec u;
	quat neworient;
	double diff;
	double rotmat[3+1][3+1];

	cn = cos(0.5*theta*M_PI/180.);
	sn = sin(0.5*theta*M_PI/180.);

	w = cn;
	x = vector.x*sn;
	y = vector.y*sn;
	z = vector.z*sn;

	neworient.w = w*part[i].orient.w - x*part[i].orient.x - y*part[i].orient.y - z*part[i].orient.z;
	neworient.x = w*part[i].orient.x + x*part[i].orient.w + y*part[i].orient.z - z*part[i].orient.y;
	neworient.y = w*part[i].orient.y - x*part[i].orient.z + y*part[i].orient.w + z*part[i].orient.x;
	neworient.z = w*part[i].orient.z + x*part[i].orient.y - y*part[i].orient.x + z*part[i].orient.w;

	part[i].orient = neworient;

	xx = x*x;
	yy = y*y;
	zz = z*z;
	xy = x*y;
	xz = x*z;
	yz = y*z;
	wx = w*x;
	wy = w*y;
	wz = w*z;

	rotmat[1][1] = 1. - 2.*(yy + zz);
	rotmat[2][1] = 2.*(xy - wz);
	rotmat[3][1] = 2.*(xz + wy);

	rotmat[1][2] = 2.*(xy + wz);
	rotmat[2][2] = 1. - 2.*(xx + zz);
	rotmat[3][2] = 2.*(yz - wx);

	rotmat[1][3] = 2.*(xz - wy);
	rotmat[2][3] = 2.*(yz + wx);
	rotmat[3][3] = 1. - 2.*(xx + yy);

	for (j = 0; j < NBALL; j++) {
		u.x = rotmat[1][1]*part[i].ballpos[j].x	+ rotmat[2][1]*part[i].ballpos[j].y	+ rotmat[3][1]*part[i].ballpos[j].z;
		u.y = rotmat[1][2]*part[i].ballpos[j].x	+ rotmat[2][2]*part[i].ballpos[j].y	+ rotmat[3][2]*part[i].ballpos[j].z;
		u.z = rotmat[1][3]*part[i].ballpos[j].x	+ rotmat[2][3]*part[i].ballpos[j].y	+ rotmat[3][3]*part[i].ballpos[j].z;

		part[i].ballpos[j].x = u.x;
		part[i].ballpos[j].y = u.y;
		part[i].ballpos[j].z = u.z;
	}
}
