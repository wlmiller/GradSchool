#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/*	MD simulation using a Langevin thermostat and "leap frog" integration
	Rotations are performed by manipulating the rotation matrix R		
    Particles are generalizations of "Janus" particles, which attract
    on one hemisphere and repel on the other. Particles in this simulation
    are permitted to have some fraction different than 50% attractive. 
    
    For more information, see:
        W. L. Miller and A. Cacciuto. Hierarchical self-assembly of 
            asymmetric amphiphatic spherical colloidal particles. 
            Phys. Rev. E 80, 021404 (2009).
*/

#define NMOL 1000			    //Particle number
#define NGROUP 1			    //Number of patches per molecule
#define MDSTEPS 10000000		//Number of MD steps
#define WRITE 100000			//Writing interval
#define EPSILON 1e-12			//Small number

typedef struct {
	double x, y, z;
} vec3D;				//3D vector

typedef struct {
	double w, i, j, k;
} quat;					//quaternion data structure	

typedef struct {
	double id[9+1];
} matrix;				//Rotational matrix	

typedef struct {
	double etaC;
	double p1diam, p1diamsq;	//Particle diameter
	double patchdiam, patchdiamsq;	//Patch particle diameter
	double side, iside, halfside;	//System side length	
	double vol, ivol;		//System volume 
	double vdiam, vdiamsq, vgap;	//Verlet parameters
	double rcut, rcutsq;		//Potential cutoff radius for patches
	double rtcut, rtcutsq;		//Interaction range of full molecules 
	double rpcut, rpcutsq;		//Interaction range of molecules via patch
	double taillength;	
	double temp, itemp;		//Temperature
	double dt, dtsq, idt, halfdt;	//Time step
	double kv, kr;			//Velocity and radial kinetic energies
	double nbar;			//Range of the angular potential
	double thetamax, thetamin;	//Angular potential cutoff
	double cosmax, cosmin;		//cos(theta*)
	double e0;			//Minimum attraction energy
	double ecut;			//Energy shift due to potential cutoff
	double range, rangesq;
} SYSTEM;

typedef struct {
	vec3D pos;			//Position
	vec3D vel;			//Velocity
	vec3D force;			//Force
	vec3D avel;			//Angular velocity
	vec3D torq;			//Torque
	vec3D angle;			//Theta, phi, psy orientation
	vec3D verl;			//Verlet coordinates
	int vlist[300];			//Verlet list
	int clist[10];			//Cluster list
	int nVerl;			//Length of Verlet list
	int nClust;
	vec3D patch[NGROUP+1];		//Positons of patches with respect to the molecule;
	vec3D bond[NGROUP+1];		//Coordinates of particle bonds (patches)
} MOLECULE;

matrix *RT;				//Transposed R matrix for each molecule
SYSTEM sys;
MOLECULE *part;

void setAll();			//Set system parameters
void setParts();		//Set initial particle coordinates and velocites
void setPatchesandAngles(int);	//Set patch positions and angular orientations
void setVerlet();		//Set Verlet lists of all particles
void checkVerlet();		//Check if the Verlet lists need to be rebuilt
void force();			//Calculate the force on all particles
void bondForces(int, int);	//Calculate forces between bonds
void integrate();		//Make a MD step of time dt
void reorientBonds();		//Apply rotation matrix to find bond coordinates
matrix updateRotMat(vec3D, int);	//Update the rotation matrix
matrix getRotMat(quat, int);	//Gnerate rotation matrix initially
vec3D rotate(matrix, vec3D);	//Rotate a vector
vec3D rotateT(matrix, vec3D);	//Rotate using the transpose of R
matrix matMult(matrix, matrix);	//Multiply two matrices together
vec3D gasdev();			//Generate a random number from a normal distribution
vec3D setRandVec();		//Set a random vector
void tempcheck();		//Verify that the velocity distribution matches the temperature
void makeMovie();

double en = 0., ecut;
double Gamma, rGamma;
double gfric, rgfric;
double noise, rnoise;
double VMAG;

FILE *moviefile, *outfile;
char filename[100], outfilename[100];
double inputvar1, inputvar2, inputvar3, inputvar4;

/****************************************/

main(int argc, char **argv) {
	int i;
	inputvar2 = 0;              // The minimum angle of interaction
	inputvar1 = atof(argv[1]);  // The maximum angle of interaction
	inputvar3 = atof(argv[2]);
	inputvar4 = atof(argv[3]);
	
	setAll();
	sys.etaC = 4./3.*M_PI*(double)NMOL/8./sys.vol;

	printf("	%lf          %lf\n",inputvar2, inputvar1);	
	printf("	%lf          %lf\n",sys.cosmin, sys.cosmax);	
	
	for (i = 1; i <= MDSTEPS; i++) {
		integrate();
		force();

		if (i % WRITE == 0) {
			tempcheck();
			printf("%d %1.2lf Kv = %lf Kr = %lf tot = %lf\n",i,180./(double)(M_PI)*(sys.thetamin+sys.thetamax)/2.,sys.kv/(double)(2*NMOL),
				sys.kr/(double)(2*NMOL),(sys.kv + sys.kr)/(double)(2*NMOL));
			makeMovie();
		}
	}
}

/****************************************/

void setAll() {
	srand48(111111);
	
	sys.taillength = 10.*(double)(M_PI)/180.;

	sys.dt = 0.0001;
	sys.dtsq = sys.dt * sys.dt;
	sys.idt = 1./sys.dt;
	sys.halfdt = 0.5*sys.dt;

	Gamma = 5.0;		//Friction coefficient
	rGamma = 2.0;		//Radial friction coefficient

	sys.temp = 1.0;
	sys.itemp = 1./sys.temp;

	VMAG = sqrt(3.*(1. - 1./(double)(NMOL))*sys.temp);

	noise = sqrt(6. * Gamma * sys.temp/sys.dt);
	rnoise = sqrt(6. * rGamma * sys.temp/sys.dt);
	gfric = 1. - Gamma*sys.dt;
	rgfric = 1. - rGamma*sys.dt;

	sys.thetamax = inputvar1*(double)(M_PI)/180.;
	if (sys.thetamax + sys.taillength <= M_PI) {
		sys.cosmax = cos(sys.thetamax + sys.taillength);
	} else {
		sys.cosmax = -1.0;
	}

	sys.thetamin = inputvar2*(double)(M_PI)/180.;
	if (sys.thetamin - sys.taillength >= 0) {
		sys.cosmin = cos(sys.thetamin - sys.taillength);
	} else {
		sys.cosmin = 1.0;
	}

	sys.etaC = 0.01;
	sys.vol = 4./3. * M_PI * (double)(NMOL)/8./sys.etaC;
	sys.ivol = 1./sys.vol;
		
	sys.side = pow(sys.vol,1./3.);
	sys.halfside = sys.side/2.;
	sys.iside = 1./sys.side;

	sys.p1diam = 1.0;
	sys.p1diamsq = sys.p1diam * sys.p1diam;

	sys.range = inputvar4*sys.p1diam;					//range of attraction
	sys.rangesq = sys.range*sys.range;

	sys.rcut = pow(2,1./6.)*sys.range;					//bottom of attractive well
	sys.rcutsq = sys.rcut*sys.rcut;

	sys.rtcut = pow(2,1./6.) * sys.p1diam;				//bottom of repulsive well (attraction will be shifted so these overlap)
	sys.rtcutsq = sys.rtcut * sys.rtcut;

	sys.rpcut = sys.p1diam + 2.5*sys.range;				//attractive potential cutoff
	sys.rpcutsq = sys.rpcut * sys.rpcut;

	sys.vdiam = 2.8 * (sys.rpcut);
	sys.vdiamsq = sys.vdiam * sys.vdiam;
	
	sys.e0 = inputvar3;
	sys.ecut = -4. * sys.e0 * (pow(sys.range/(sys.rpcut - sys.rtcut + sys.rcut),12) - pow(sys.range/(sys.rpcut - sys.rtcut + sys.rcut),6));

	sys.vgap = 0.5 * (sys.vdiam - sys.rpcut);

	part = (MOLECULE *) malloc((NMOL+1)*sizeof(MOLECULE));
	RT = (matrix *) malloc((NMOL+1)*sizeof(matrix));

	sprintf(filename,"theta_min_%1.0lf_max_%1.0lf_e0_%1.2lf.xyz", inputvar2, inputvar1, sys.e0);
	
	moviefile = fopen(filename,"w");
	fclose(moviefile);

	setParts();

	setVerlet();
}

/****************************************/

void setParts() {
	int i, j, k = 0;
	double dx, dy, dz, diffsq;	//Particle distance
	double Vsq;			//Velocity squared
	double Vx, Vy, Vz;		//Velocity
	double scale;			//Scaling factor
	double tt;
	vec3D e, tmp;

	double ll = sqrt(NMOL);

	for (i = 1; i <= NMOL; i++) {
		part[i].force.x = 0.;
		part[i].force.y = 0.;
		part[i].force.z = 0.;

		part[i].torq.x = 0.;
		part[i].torq.y = 0.;
		part[i].torq.z = 0.;

		again:
			part[i].pos.x = (0.5 - drand48()) * sys.side;
			part[i].pos.y = (0.5 - drand48()) * sys.side;
			part[i].pos.z = (0.5 - drand48()) * sys.side;

			dx = part[i].pos.x;
			dy = part[i].pos.y;
			dz = part[i].pos.z;

			diffsq = dx*dx + dy*dy + dz*dz;

			for (j = 1; j < i; j++) {
				dx = part[i].pos.x - part[j].pos.x;
				dy = part[i].pos.y - part[j].pos.y;
				dz = part[i].pos.z - part[j].pos.z;
				
				if (dx > sys.halfside) dx -= sys.side;
				if (dy > sys.halfside) dy -= sys.side;
				if (dz > sys.halfside) dz -= sys.side;

				if (dx < -sys.halfside) dx += sys.side;
				if (dy < -sys.halfside) dy += sys.side;
				if (dz < -sys.halfside) dz += sys.side;

				diffsq = dx*dx + dy*dy + dz*dz;

				if (diffsq < sys.p1diamsq) goto again;
			}
		setPatchesandAngles(i);
		reorientBonds(i);
	}

	/*Rescale to be consistent with the temperature*/

	Vx = 0.;
	Vy = 0.;
	Vz = 0.;

	for (i = 1; i <= NMOL; i++) {
		part[i].vel = gasdev();
		
		part[i].vel.x *= VMAG;
		part[i].vel.y *= VMAG;
		part[i].vel.z *= VMAG;
		
		Vx += part[i].vel.x;
		Vy += part[i].vel.y;
		Vz += part[i].vel.z;
	}

	Vsq = 0.;

	for (i = 1; i <= NMOL; i++) {
		part[i].vel.x = part[i].vel.x - Vx/(double)(NMOL);
		part[i].vel.y = part[i].vel.y - Vy/(double)(NMOL);
		part[i].vel.z = part[i].vel.z - Vz/(double)(NMOL);

		Vsq += (part[i].vel.x*part[i].vel.x + 
			part[i].vel.y*part[i].vel.y + 
			part[i].vel.z*part[i].vel.z);
	}

	for (i = 1; i <= NMOL; i++) {
		e = setRandVec();

		tt = e.x*e.x + e.y*e.y + e.z*e.z;

		tt = VMAG/sqrt(tt);

		tmp.x = tt*e.x;
		tmp.y = tt*e.y;
		tmp.z = tt*e.z;
		
		part[i].avel = rotate(RT[i],tmp);
	}
	
	force();

	for (i = 1; i <= NMOL; i++) {	//For the Leap Frog algorithm, velocity should be dt/2 behind position
		part[i].vel.x = part[i].vel.x - 0.5*sys.dt*part[i].force.x;
		part[i].vel.y = part[i].vel.y - 0.5*sys.dt*part[i].force.y;
		part[i].vel.z = part[i].vel.z - 0.5*sys.dt*part[i].force.z;

		part[i].avel.x = part[i].avel.x - 0.5*sys.dt*part[i].torq.x;
		part[i].avel.y = part[i].avel.y - 0.5*sys.dt*part[i].torq.y;
		part[i].avel.z = part[i].avel.z - 0.5*sys.dt*part[i].torq.z;
	}

	makeMovie();
	tempcheck();
	
	printf("	Kv = %lf Kr = %lf\n",sys.kv/(double)(2*NMOL),sys.kr/(double)(2*NMOL));
}

/****************************************/

void setPatchesandAngles(int k) {	//Define patches in particle reference frame
	double s, x, y;
	double a1, a2, a3;
	vec3D e;
	quat q;

	part[k].patch[1].x = 0.;
	part[k].patch[1].y = 0.;
	part[k].patch[1].z = 0.5*sys.p1diam;

	e = setRandVec();

	part[k].angle.x = atan2(e.x,e.y); 	//phi
	part[k].angle.y = acos(e.z);		//theta
	part[k].angle.z = 2 * M_PI * drand48();	//psi

	a1 = 0.5 * part[k].angle.y;
	a2 = 0.5 * (part[k].angle.x - part[k].angle.z);
	a3 = 0.5 * (part[k].angle.x + part[k].angle.z);

	/*Set quaternion holding angular coordinates - for fast rotation only*/

	q.w = cos(a1) * cos(a3);
	q.i = sin(a1) * cos(a2);
	q.j = sin(a1) * sin(a2);
	q.k = cos(a1) * sin(a3);

	/*Reorient patches accordingly*/
	RT[k] = getRotMat(q,1);
}

/****************************************/

void setVerlet() {

	int i, j;
	double dx, dy, dz, diffsq;

	for (i = 1; i <= NMOL; i++) {
		part[i].nVerl = 0;
		part[i].verl.x = part[i].pos.x;
		part[i].verl.y = part[i].pos.y;
		part[i].verl.z = part[i].pos.z;
	}
	
	for (i = 1; i < NMOL; i++) {
		for (j = i + 1; j <= NMOL; j++) {
			dx = part[i].pos.x - part[j].pos.x;
			dy = part[i].pos.y - part[j].pos.y;
			dz = part[i].pos.z - part[j].pos.z;

			if (dx > sys.halfside) dx -= sys.side;
			if (dy > sys.halfside) dy -= sys.side;
			if (dz > sys.halfside) dz -= sys.side;

			if (dx < -sys.halfside) dx += sys.side;
			if (dy < -sys.halfside) dy += sys.side;
			if (dz < -sys.halfside) dz += sys.side;
			
			diffsq = (dx*dx + dy*dy + dz*dz);

			if (diffsq < sys.vdiamsq) {
				part[i].nVerl++;
				part[i].vlist[part[i].nVerl] = j;

				part[j].nVerl++;
				part[j].vlist[part[j].nVerl] = i;
			}
		}
	}
}

/****************************************/

void checkVerlet(int k) {
	double dx, dy, dz, diffsq, diff;

	dx = part[k].pos.x - part[k].verl.x;
	dy = part[k].pos.y - part[k].verl.y;
	dz = part[k].pos.z - part[k].verl.z;

	if (dx > sys.halfside) dx -= sys.side;
	if (dy > sys.halfside) dy -= sys.side;
	if (dz > sys.halfside) dz -= sys.side;

	if (dx < -sys.halfside) dx += sys.side;
	if (dy < -sys.halfside) dy += sys.side;
	if (dz < -sys.halfside) dz += sys.side;

	diffsq = dx*dx + dy*dy + dz*dz;
	diff = sqrt(diffsq);

	if (diff > sys.vgap) setVerlet();
}	
	 
/****************************************/

void force() {
	int i, j;
	double dx, dy, dz, diffsq;
	double force;
	double idiffsq, idiff6;
	vec3D n1, n2, nn;

	for (i = 1; i <= NMOL; i++) {

        /*Thermostat - adds random fluctuations to force*/
		nn = gasdev();
		part[i].force.x = noise*nn.x;
		part[i].force.y = noise*nn.y;
		part[i].force.z = noise*nn.z;

		n1 = gasdev();
		part[i].torq.x = n1.x*rnoise;
		part[i].torq.y = n1.y*rnoise;
		part[i].torq.z = n1.z*rnoise;

		for (j = 1; j <= part[i].nVerl; j++) {
			dx = part[i].pos.x - part[part[i].vlist[j]].pos.x;
			dy = part[i].pos.y - part[part[i].vlist[j]].pos.y;
			dz = part[i].pos.z - part[part[i].vlist[j]].pos.z;

			if (dx > sys.halfside) dx -= sys.side;
			if (dy > sys.halfside) dy -= sys.side;
			if (dz > sys.halfside) dz -= sys.side;

			if (dx < -sys.halfside) dx += sys.side;
			if (dy < -sys.halfside) dy += sys.side;
			if (dz < -sys.halfside) dz += sys.side;

			diffsq = dx*dx + dy*dy + dz*dz;

			/*If the particles are close enough for their patches to interact, calculate the force*/	
            if (diffsq <= sys.rtcutsq) {
				idiffsq = sys.p1diamsq/diffsq;
				idiff6 = idiffsq * idiffsq * idiffsq;

				force = sys.e0*48./diffsq * idiff6 * (idiff6 - 0.5);

				part[i].force.x += force * dx;
				part[i].force.y += force * dy;
				part[i].force.z += force * dz;
			} 
			if (diffsq <= sys.rpcutsq) {
				bondForces(i,part[i].vlist[j]);
			}
		}
	}	
}

/****************************************/
/*bondForces calculates attractions and
 * torques due to patches.*/
void bondForces (int k, int j) {
	int s,t;
	double dx, dy, dz, diffsq, diff, tmp, diffsq2, diff2;
	double force, force1;
	double idiffsq, idiff6, idiffsqb, idiff6b;
	vec3D n1, n2, nn, b1, b2, cross, fA, FT;
	double b1b2, vtheta1, vtheta2, theta1, theta2, ntheta1, ntheta2, b1v, b2v, n2theta1, n2theta2;
	double ntheta;
	int outside = 1;

	for (s = 1; s <= NGROUP; s++) {
		FT.x = 0.;
		FT.y = 0.;
		FT.z = 0.;

		b1.x = part[k].bond[s].x - part[k].pos.x;
		b1.y = part[k].bond[s].y - part[k].pos.y;
		b1.z = part[k].bond[s].z - part[k].pos.z;

		dx = part[k].pos.x - part[j].pos.x;
		dy = part[k].pos.y - part[j].pos.y;
		dz = part[k].pos.z - part[j].pos.z;

		if (dx > sys.halfside) dx -= sys.side;
		if (dy > sys.halfside) dy -= sys.side;
		if (dz > sys.halfside) dz -= sys.side;

		if (dx < -sys.halfside) dx += sys.side;
		if (dy < -sys.halfside) dy += sys.side;
		if (dz < -sys.halfside) dz += sys.side;
		
		diffsq = dx*dx + dy*dy + dz*dz;
		diff = sqrt(diffsq);

		for (t = 1; t <= NGROUP; t++) {
			b2.x = (part[j].bond[t].x - part[j].pos.x);
			b2.y = (part[j].bond[t].y - part[j].pos.y);
			b2.z = (part[j].bond[t].z - part[j].pos.z);
			
			b1v = 2./diff * -(b1.x*dx + b1.y*dy + b1.z*dz);
			b2v = 2./diff * (b2.x*dx + b2.y*dy + b2.z*dz);
			b1b2 = 4. * (b1.x*b2.x + b1.y*b2.y + b1.z*b2.z);	//Times 4 since |b| = 0.5sigma

			if ((b1v >= sys.cosmax) && (b2v >= sys.cosmax) && (b1v <= sys.cosmin) && (b2v <= sys.cosmin)) {
				diff2 = diff - sys.rtcut + sys.rcut;
				diffsq2 = diff2*diff2;
					
				if (1) {
					idiffsq = sys.rangesq/diffsq2;
					idiff6 = idiffsq * idiffsq * idiffsq;
					if (diffsq > sys.rtcutsq) {
						force = sys.e0 * 48./diffsq2 * idiff6 * (idiff6 - 0.5);
					} else force = 0;

					if (b1v > (1 - EPSILON)) b1v = (1-EPSILON);
					if (b1v < (-1 + EPSILON)) b1v = (-1+EPSILON);
					if (b2v > (1 - EPSILON)) b2v = (1-EPSILON);
					if (b2v < (-1 + EPSILON)) b2v = (-1+EPSILON);

					theta1 = acos(b1v);
					theta2 = acos(b2v);

					ntheta1 = (theta1 - sys.thetamax)*((double)(M_PI)/(2*sys.taillength));
					ntheta2 = (theta2 - sys.thetamax)*((double)(M_PI)/(2*sys.taillength));

					n2theta1 = (sys.thetamin - theta1)*((double)(M_PI)/(2*sys.taillength));
					n2theta2 = (sys.thetamin - theta2)*((double)(M_PI)/(2*sys.taillength));

					if (theta1 > sys.thetamax) {
						vtheta1 = cos(ntheta1)*cos(ntheta1);
					} else if (theta1 < sys.thetamin) {
						vtheta1 = cos(n2theta1)*cos(n2theta1);
					} else {
						vtheta1 = 1;
					}	

					if (theta2 > sys.thetamax) {
						vtheta2 = cos(ntheta2)*cos(ntheta2);
					} else if (theta2 < sys.thetamin) {
						vtheta2 = cos(n2theta2)*cos(n2theta2);
					} else {
						vtheta2 = 1;
					}

					fA.x = force*dx*vtheta1*vtheta2;
					fA.y = force*dy*vtheta1*vtheta2;
					fA.z = force*dz*vtheta1*vtheta2;
					
					if (diffsq > sys.rtcutsq) {
						force1 = (4. * sys.e0 * idiff6 * (idiff6 - 1.)) + sys.ecut;
					} else {
						idiffsqb = sys.p1diam/diffsq;
						idiff6b = idiffsqb*idiffsqb*idiffsqb;
						force1 = (4. * sys.e0 * idiff6b * (idiff6b - 1.)) + sys.ecut;
					}

					if (theta1 > sys.thetamax) {
						fA.x += force1*b1.x*cos(ntheta1)*sin(ntheta1)/sin(theta1)*(double)(M_PI)/sys.taillength*2./diff*vtheta2;
						fA.y += force1*b1.y*cos(ntheta1)*sin(ntheta1)/sin(theta1)*(double)(M_PI)/sys.taillength*2./diff*vtheta2;
						fA.z += force1*b1.z*cos(ntheta1)*sin(ntheta1)/sin(theta1)*(double)(M_PI)/sys.taillength*2./diff*vtheta2;
					} else if (theta1 < sys.thetamin) {
						fA.x += -force1*b1.x*cos(n2theta1)*sin(n2theta1)/sin(theta1)*(double)(M_PI)/sys.taillength*2./diff*vtheta2;
						fA.y += -force1*b1.y*cos(n2theta1)*sin(n2theta1)/sin(theta1)*(double)(M_PI)/sys.taillength*2./diff*vtheta2;
						fA.z += -force1*b1.z*cos(n2theta1)*sin(n2theta1)/sin(theta1)*(double)(M_PI)/sys.taillength*2./diff*vtheta2;
					}

					if (theta2 > sys.thetamax) {
						fA.x += -force1*b2.x*cos(ntheta2)*sin(ntheta2)/sin(theta2)*(double)(M_PI)/sys.taillength*2./diff*vtheta1;
						fA.y += -force1*b2.y*cos(ntheta2)*sin(ntheta2)/sin(theta2)*(double)(M_PI)/sys.taillength*2./diff*vtheta1;
						fA.z += -force1*b2.z*cos(ntheta2)*sin(ntheta2)/sin(theta2)*(double)(M_PI)/sys.taillength*2./diff*vtheta1;
					} else if (theta2 < sys.thetamin) {
						fA.x += force1*b2.x*cos(n2theta2)*sin(n2theta2)/sin(theta2)*(double)(M_PI)/sys.taillength*2./diff*vtheta1;
						fA.y += force1*b2.y*cos(n2theta2)*sin(n2theta2)/sin(theta2)*(double)(M_PI)/sys.taillength*2./diff*vtheta1;
						fA.z += force1*b2.z*cos(n2theta2)*sin(n2theta2)/sin(theta2)*(double)(M_PI)/sys.taillength*2./diff*vtheta1;
					}


					part[k].force.x += fA.x;
					part[k].force.y += fA.y;
					part[k].force.z += fA.z;

					FT.x += fA.x;
					FT.y += fA.y;
					FT.z += fA.z;

					if (theta1 <= sys.thetamax && theta1 >= sys.thetamin) {
						tmp = 0.;
					} else if (theta1 > sys.thetamax) {
						tmp = -(double)(M_PI)/(sys.taillength)*cos(ntheta1)*sin(ntheta1)/sin(theta1);
					} else if (theta1 < sys.thetamin) {
						tmp = (double)(M_PI)/(sys.taillength)*cos(n2theta1)*sin(n2theta1)/sin(theta1);
					}

					if (theta2 > sys.thetamax) {
						tmp = tmp*cos(ntheta2)*cos(ntheta2);
					} else if (theta2 < sys.thetamin) {
						tmp = tmp*cos(n2theta2)*cos(n2theta2);
					}

					cross.x = 2./diff*(b1.y*dz - b1.z*dy);
					cross.y = 2./diff*(b1.z*dx - b1.x*dz);
					cross.z = 2./diff*(b1.x*dy - b1.y*dx);
					
					part[k].torq.x += -force1 * tmp * cross.x;
					part[k].torq.y += -force1 * tmp * cross.y;
					part[k].torq.z += -force1 * tmp * cross.z;
				}
			}	
		}

	}

}

/****************************************/
/* Perform a molecular dynamics step    */
void integrate() {
	int i;
	
	for (i = 1; i <= NMOL; i++) {
		part[i].avel.x = part[i].avel.x*rgfric + sys.dt*part[i].torq.x;
		part[i].avel.y = part[i].avel.y*rgfric + sys.dt*part[i].torq.y;
		part[i].avel.z = part[i].avel.z*rgfric + sys.dt*part[i].torq.z;

		RT[i] = updateRotMat(part[i].avel,i);

		part[i].vel.x = part[i].vel.x*gfric + sys.dt*part[i].force.x;
		part[i].vel.y = part[i].vel.y*gfric + sys.dt*part[i].force.y;
		part[i].vel.z = part[i].vel.z*gfric + sys.dt*part[i].force.z;

		part[i].pos.x += part[i].vel.x*sys.dt;
		part[i].pos.y += part[i].vel.y*sys.dt;
		part[i].pos.z += part[i].vel.z*sys.dt;

		if (part[i].pos.x > sys.halfside) part[i].pos.x -= sys.side;
		if (part[i].pos.y > sys.halfside) part[i].pos.y -= sys.side;
		if (part[i].pos.z > sys.halfside) part[i].pos.z -= sys.side;

		if (part[i].pos.x < -sys.halfside) part[i].pos.x += sys.side;
		if (part[i].pos.y < -sys.halfside) part[i].pos.y += sys.side;
		if (part[i].pos.z < -sys.halfside) part[i].pos.z += sys.side;

		reorientBonds(i);
		checkVerlet(i);
	}	
}

/****************************************/

void reorientBonds(int k) {
	int j;
	vec3D shift;

	for (j = 1; j <= NGROUP; j++) {
		shift = rotate(RT[k],part[k].patch[j]);

		part[k].bond[j].x = part[k].pos.x + shift.x;
		part[k].bond[j].y = part[k].pos.y + shift.y;
		part[k].bond[j].z = part[k].pos.z + shift.z;
	}	
}

/****************************************/

matrix updateRotMat(vec3D angle, int k) {		//Optimized
	double c[3], s[3];
	double t;
	double c0s2, c0c2, s0c2, s0s2;
	matrix M1;					//UxUyUz
	matrix M2;					//UzUyUx
	matrix M1M2;					//UxUyUzUzUyUx
	matrix Rnew;					//Updated R matrix

	angle.x = angle.x * sys.halfdt;
	angle.y = angle.y * sys.halfdt;
	angle.z = angle.z * sys.halfdt;

	t = .25 * angle.x * angle.x;
	c[0] = (1. - t)/(1. + t);
	s[0] = angle.x/(1. + t);

	t = .25 * angle.y * angle.y;
	c[1] = (1. - t)/(1. + t);
	s[1] = angle.y/(1. + t);

	t = .25 * angle.z * angle.z;
	c[2] = (1. - t)/(1. + t);
	s[2] = angle.z/(1. + t);

	c0c2 = c[0]*c[2];
	c0s2 = c[0]*s[2];
	s0c2 = s[0]*c[2];
	s0s2 = s[0]*s[2];

	/*Matrix is numbered down columns and then across!*/

	M1.id[1] = c[1]*c[2];
	M1.id[2] = s0c2*s[1] + c0s2;
	M1.id[3] = -c0c2*s[1] + s0s2;

	M1.id[4] = -c[1]*s[2];
	M1.id[5] = -s0s2*s[1] + c0c2;
	M1.id[6] = c0s2*s[1] + s0c2;

	M1.id[7] = s[1];
	M1.id[8] = -s[0]*c[1];
	M1.id[9] = c[0]*c[1];

	M2.id[1] = M1.id[1];
	M2.id[2] = -M1.id[4];
	M2.id[3] = -M1.id[7];

	M2.id[4] = s0c2*s[1] - c0s2;
	M2.id[5] = s0s2*s[1] + c0c2;
	M2.id[6] = -M1.id[8];

	M2.id[7] = c0c2*s[1] + s0s2;
	M2.id[8] = c0s2*s[1] - s0c2;
	M2.id[9] = M1.id[9];

	M1M2 = matMult(M1,M2);
	Rnew = matMult(M1M2,RT[k]);

	return Rnew;
}

/****************************************/

matrix getRotMat(quat Q, int transp) {				//Optimized
	int s;
	double qiqi, qiqj, qiqk, qiqw;
	double qjqj, qjqk, qjqw;
	double qkqk, qkqw;
	double qwqw;
	matrix RQ;

	if (transp == 1) s = 1;					//transpose the matrix
	else s = -1;						//don't transpose it

	qiqi = Q.i*Q.i;
	qiqj = Q.i*Q.j;
	qiqk = Q.i*Q.k;
	qiqw = Q.i*Q.w;

	qjqj = Q.j*Q.j;
	qjqk = Q.j*Q.k;
	qjqw = Q.j*Q.w;

	qkqk = Q.k*Q.k;
	qkqw = Q.k*Q.w;

	qwqw = Q.w*Q.w;

	RQ.id[1] = 2.*(qiqi + qwqw) - 1.;
	RQ.id[2] = 2.*(qiqj + s*qkqw);
	RQ.id[3] = 2.*(qiqk - s*qjqw);

	RQ.id[4] = 2.*(qiqj - s*qkqw);
	RQ.id[5] = 2.*(qjqj + qwqw) - 1.;
	RQ.id[6] = 2.*(qjqk + s*qiqw);

	RQ.id[7] = 2.*(qiqk + s*qjqw);
	RQ.id[8] = 2.*(qjqk - s*qiqw);
	RQ.id[9] = 2.*(qkqk + qwqw) - 1.;

	return RQ;
}	

/****************************************/

vec3D rotate(matrix R, vec3D s) {
	vec3D r;

	r.x = R.id[1]*s.x + R.id[4]*s.y + R.id[7]*s.z;
	r.y = R.id[2]*s.x + R.id[5]*s.y + R.id[8]*s.z;
	r.z = R.id[3]*s.x + R.id[6]*s.y + R.id[9]*s.z;

	return r;
}

/****************************************/

vec3D rotateT(matrix R, vec3D s) {
	vec3D r;

	r.x = R.id[1]*s.x + R.id[2]*s.y + R.id[3]*s.z;
	r.y = R.id[4]*s.x + R.id[5]*s.y + R.id[6]*s.z;
	r.z = R.id[7]*s.x + R.id[8]*s.y + R.id[9]*s.z;

	return r;
}

/****************************************/

matrix matMult(matrix A, matrix B) {
	matrix U;

	U.id[1] = A.id[1]*B.id[1] + A.id[4]*B.id[2] + A.id[7]*B.id[3];
	U.id[2] = A.id[2]*B.id[1] + A.id[5]*B.id[2] + A.id[8]*B.id[3];
	U.id[3] = A.id[3]*B.id[1] + A.id[6]*B.id[2] + A.id[9]*B.id[3];

	U.id[4] = A.id[1]*B.id[4] + A.id[4]*B.id[5] + A.id[7]*B.id[6];
	U.id[5] = A.id[2]*B.id[4] + A.id[5]*B.id[5] + A.id[8]*B.id[6];
	U.id[6] = A.id[3]*B.id[4] + A.id[6]*B.id[5] + A.id[9]*B.id[6];
	
	U.id[7] = A.id[1]*B.id[7] + A.id[4]*B.id[8] + A.id[7]*B.id[9];
	U.id[8] = A.id[2]*B.id[7] + A.id[5]*B.id[8] + A.id[8]*B.id[9];
	U.id[9] = A.id[3]*B.id[7] + A.id[6]*B.id[8] + A.id[9]*B.id[9];

	return U;
}

/****************************************/

vec3D gasdev() {		//"Box-Muller" Gaussian random number generator
	double s, x, y;
	vec3D p;

	s = 2;
	while (s > 1.) {
		x = 2.*drand48() - 1;
		y = 2.*drand48() - 1;
		s = x*x + y*y;
	}

	p.z = 1. - 2.*s;
	s = 2.*sqrt(1.-s);
	p.x = s*x;
	p.y = s*y;

	return p;
}

/****************************************/

vec3D setRandVec() {
	double x, y, s = 2.;
	vec3D e;

	while (s > 1.) {
		x = 2. * drand48() - 1.;
		y = 2. * drand48() - 1.;
		s = x*x + y*y;
	}

	e.z = 1. - 2.*s;
	s = 2.*sqrt(1.-s);
	e.x = s * x;
	e.y = s * y;

	return e;
}

/****************************************/

void tempcheck() {
	int i;

	sys.kv = 0.;
	sys.kr = 0.;

	for (i = 1; i <= NMOL; i++) {
		sys.kv += (part[i].vel.x*part[i].vel.x + 
			part[i].vel.y*part[i].vel.y + 
			part[i].vel.z*part[i].vel.z);

		sys.kr += (part[i].avel.x*part[i].avel.x + 
			part[i].avel.y*part[i].avel.y + 
			part[i].avel.z*part[i].avel.z);
	}
}	

/****************************************/

void makeMovie(void){
	int i, j;
	vec3D angle, RR;
 
	moviefile = fopen(filename,"a");
 
	fprintf(moviefile,"%d\n\n",(NGROUP * NMOL) + NMOL);  
	for (i = 1; i <= NMOL; i++) {
		fprintf(moviefile,"H %lf %lf %lf\n", part[i].pos.x, part[i].pos.y, part[i].pos.z);

		for (j = 1; j <= NGROUP; j++) {
			fprintf(moviefile, "HZ %lf %lf %lf\n",part[i].bond[j].x,
					part[i].bond[j].y, part[i].bond[j].z);
		}			
	}
  
	fclose(moviefile);  
}
