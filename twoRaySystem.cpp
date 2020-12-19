/* fastOPDTilt.cpp, R Miyakawa, 11/2013
 *
 * OPD = fastOPD(Xcoords, Ycoords, z1_um, z2_um, lambda_um, T_um, orders, alpha (tiltX), beta (tiltY));
 * 
 * Computes the OPD from the self-interference of a wave incident on a 
 * defocused grating in an LSI experiment.  This function can be used to study
 * the systematic aberrations associated with using a grating at high numerical
 * apertures.  While this result is not perfectly analytical, it can be 
 * solved numerically much more quickly than running a simulation.
 *
 * Here Xcoords and Ycoords are arrays containing the detector coordinates.
 *
 * Orders is a 1x2 array specifying which 2 grating orders are interfering
 * e.g., [1, -1] specifies that its the +/- 1 orders that are interfering.
 *
 */


#include <stdio.h>
#include <mex.h>
#include <math.h>


// Recursively computes the binomial coefficient nCk:
double nChooseK(double n, double k) {
    if (k == 0 && n >= 0)
        return 1;
    else if (n == 0 && k > 0)
        return 0;
    else
        return nChooseK(n - 1, k - 1) + nChooseK(n - 1, k);
}

// Retrieves the value the Zernike Polynomial order J at the coordinate [r, th]:
double zgenpt(int j, double r, double th) {
    
    // Get dual index [n,m] from j:
    int smct = 2;
    int ct = 0;
    int numInSmct = 3;
    while (j > 0) {
        smct += 2;
        j -= numInSmct;
        numInSmct = smct + 1;
        ct++;
    }
    
    int n = 2 * ct;
    int m = 0;
    
    for (int q = 1; q <= abs(j); q++) {
        if (q % 2 == 1) {
            n--;
            m++;
        }
    }
    if ((abs(j)) % 2 == 1) {
        m *= -1;
    }
    
    // Compute value of radial function:
    int p = (n - abs(m)) / 2;
    double kParity = 1;
    double Rv = 0;
    
    for (int k = 0; k <= p; k += 1) {
        if (k % 2 == 1)
            kParity = -1;
        else
            kParity = 1;
        
        Rv += kParity*nChooseK((double)(n - k), (double)k)
            *nChooseK((double)(n - 2 * k), (double)(p - k))
            *pow(r, (double)(n - 2 * k));
    }

    // Compute value of azimuthal function:
    double Av = 0;
    if (m > 0)
        Av = cos(m*th);
    else if (m < 0)
        Av = -sin(m*th);
    else
        Av = 1;
    
    return Rv*Av;
}

// Objective function to solve for theta as a function of detector coordinate
double getThetaObjective(double th, double X, double Z, double alpha, double z1, 
                         double z2, double lambda, double T, int order){
    
    double s1       = z1 * sin(th)/cos(th - alpha);
    double sinphi   = sin(th-alpha) + order*lambda/T;
    double tanphi   = tan(asin(sinphi));
            
    //printf("tanphi: %0.5f\n", tanphi);
    return (z1 - s1 * sin(alpha)) * tan(th) - (z2 - s1 * sin(alpha)+Z)*tanphi - X;
    
}



double secantSolveTheta(double xm2, double xm1, double X,double Z, double alpha, double z1,
                double z2, double lambda, double T, int order,
                    int currentIter, int maxIter, double tolX){
    
    
    // Compute the value of f(xm1) and f(xm2)
    double fxm1 = getThetaObjective(xm1, X, Z, alpha, z1, z2, lambda, T, order);
    double fxm2 = getThetaObjective(xm2, X, Z, alpha, z1, z2, lambda, T, order);
    
    // Compute new guess
    double x0 = xm1 - fxm1*(xm1 - xm2)/(fxm1 - fxm2);
                    
    // Exit because iterations have exceeded limit
    if (currentIter >= maxIter){
        printf("*******MAXIMUM ITERATIONS REACHED");
        return x0;
    }
     
    // exit because function value is within tolerance
    if (sqrt((x0 - xm1)*(x0 - xm1)) < tolX ){
        //printf("Total iterations: %d\n", currentIter);
        //printf("x0 = %e, xm1 = %e, diff = %e, tolX = %e\n", x0, xm1,sqrt((x0 - xm1)*(x0 - xm1)), tolX);
        return x0;
    }
    
    // Call recursively with new guess
    return secantSolveTheta(xm1, x0, X, Z, alpha, z1, z2, lambda, T, order, 
                        currentIter + 1, maxIter, tolX);
   
}




// Main function: generates an array of coordinates of trapezoids representing
// the desired zoneplate. A pointer to the 2D array is returned
double * computeOPD(double * X, double * Y, double * Z, double alpha, double beta, double z1,
                    double z2, double lambda,
                    double T, double NA, double * orders, long numPts, double * zernikeCouples, int nz, bool isX){
    
    /* Split computation into three parts
    
     * 1: Based on detector coords, compute thetas of rays coming to virtual 
     *    focus.  Use these to compute pathlength from virutal focus to 
     *    grating.
     *
     * 2: Compute phase from grating itself.  Make sure to apply opposite
     *    signs to phase from + and - orders.
     *
     * 3: Compute pathlength from grating to detector
     * 
    */
    
    double plp1, plp2, plp3, plp4, plm1, plm2, plm3, plm4, thp, thm, thpy, thmy, sXp, sXm, sYp, sYm, thisOPD;
    double grx_p, grx_m;
    double s1p, s1m, s1py, s1my;
    double thisX, thisY ,thisZ;
    double thgxp, thgxm;
    double ppx, pmx, ppy, pmy, ppr, ppth, pmr, pmth;
    
    double * OPD;
    
    OPD = new double[numPts];
    
    for (long k = 0; k < numPts; k++){
        
        thisX   = X[k];
        thisY   = Y[k];
        thisZ   = Z[k];

        thgxp = 0;
        thgxm = 0;
        
        // Part 1:
        if (isX){ // X
            thp     = secantSolveTheta(thgxp, thgxp + .01, thisX, thisZ, alpha, z1, z2, lambda, T, orders[0], 0, 100, 1e-6);
            thm     = secantSolveTheta(thgxm, thgxm + .01, thisX, thisZ, alpha, z1, z2, lambda, T, orders[1], 0, 100, 1e-6);
            thpy     = secantSolveTheta(thgxp, thgxp + .01, thisY, thisZ, beta, z1, z2, lambda, T, orders[0], 0, 100, 1e-6);
            thmy     = secantSolveTheta(thgxm, thgxm + .01, thisY, thisZ, beta, z1, z2, lambda, T, orders[0], 0, 100, 1e-6);
        } 
        else { // Y
            thp     = secantSolveTheta(thgxp, thgxp + .01, thisX, thisZ, alpha, z1, z2, lambda, T, orders[0], 0, 100, 1e-6);
            thm     = secantSolveTheta(thgxm, thgxm + .01, thisX, thisZ, alpha, z1, z2, lambda, T, orders[0], 0, 100, 1e-6);
            thpy     = secantSolveTheta(thgxp, thgxp + .01, thisY, thisZ, beta, z1, z2, lambda, T, orders[0], 0, 100, 1e-6);
            thmy     = secantSolveTheta(thgxm, thgxm + .01, thisY, thisZ, beta, z1, z2, lambda, T, orders[1], 0, 100, 1e-6);
        }
        

        
        // Generate pupil coordinates projected onto plane
        
        ppx = tan(thp)/tan(asin(NA));
        pmx = tan(thm)/tan(asin(NA));
        
        ppy = tan(thpy)/tan(asin(NA));
        pmy = tan(thmy)/tan(asin(NA));
        
        // Transform to polar coordinates
        ppr = sqrt(ppx*ppx + ppy*ppy);
        ppth = atan2(ppy, ppx);
        pmr = sqrt(pmx*pmx + pmy*pmy);
        pmth = atan2(pmy, pmx);
        
        
        s1p     = z1 * sin(thp)/cos(thp - alpha);
        s1m     = z1 * sin(thm)/cos(thm - alpha);
        s1py    = z1 * sin(thpy)/cos(thpy - beta);
        s1my     = z1 * sin(thmy)/cos(thmy - beta);
        //printf("th1: %0.5f, th2: %0.5f\n", thp, thm);

        // Coordinates of ray at grating plane
        sXp     = (z1 - s1p*sin(alpha))*tan(thp);//???
        sXm     = (z1 - s1m*sin(alpha))*tan(thm);
        sYp     = (z1 - s1py*sin(beta))*tan(thpy);//???
        sYm     = (z1 - s1my*sin(beta))*tan(thmy);
        //sYp     = thisY * z1/(z2 + z1);
        //sYm     = thisY * z1/(z2 + z1);

        plp1    = -sqrt(sXp*sXp + sYp*sYp + (z1 - s1p*sin(alpha) - s1py*sin(beta))*(z1 - s1p*sin(alpha) - s1py*sin(beta)));//+(z1 - s1py*sin(beta))*(z1 - s1py*sin(beta)));
        plm1    = -sqrt(sXm*sXm + sYm*sYm + (z1 - s1m*sin(alpha) - s1my*sin(beta))*(z1 - s1m*sin(alpha) - s1my*sin(beta)));//+(z1 - s1my*sin(beta))*(z1 - s1my*sin(beta)));
        

        // Part 2:
        if (isX){
            plp2 = -orders[0]*lambda * s1p / T;
            plm2 = -orders[1]*lambda * s1m / T;
        } else {
            plp2 = -orders[0]*lambda * s1py / T;
            plm2 = -orders[1]*lambda * s1my / T;
        }
        
        // Part 3:
        plp3 = sqrt( (sXp - thisX)*(sXp - thisX) + (sYp - thisY)*(sYp - thisY)
                            + (z2 - s1p * sin(alpha) - s1py*sin(beta)- thisZ)*(z2 - s1p * sin(alpha) - s1py*sin(beta)- thisZ)  );
        plm3 = sqrt( (sXm - thisX)*(sXm - thisX) + (sYm - thisY)*(sYm - thisY)
                            + (z2 - s1m * sin(alpha) - s1my*sin(beta)- thisZ)*(z2 - s1m * sin(alpha) - s1my*sin(beta)- thisZ)  );

        // Part 4:
        plp4 = 0;
        plm4 = 0;
        for (int n = 0; n < nz; n++){
            plp4 += zernikeCouples[2*n + 1] * zgenpt((int)zernikeCouples[2*n], ppr, ppth) * lambda; 
            plm4 += zernikeCouples[2*n + 1] * zgenpt((int)zernikeCouples[2*n], pmr, pmth) * lambda; 
        }
        
        thisOPD = (plp1 + plp2 + plp3 + plp4) - (plm1 + plm2 + plm3 + plm4);
        OPD[k] = thisOPD;
    }
    
    return OPD;
    
    
}

// Gateway function for MATLAB interface
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) 
{
  
    mwSize rowLen, colLen , rowAber;
    double *detectorCoordsX; 
    double *detectorCoordsY;
    double *detectorCoordsZ;
    double *zernikeCouples;
    double *OPDx, *OPDy;
    double *orders;
    double alpha, beta;
    
    
    detectorCoordsX = mxGetPr(prhs[0]);
    detectorCoordsY = mxGetPr(prhs[1]);
    detectorCoordsZ = mxGetPr(prhs[2]);
    rowLen          = mxGetM(prhs[0]);
    colLen          = mxGetN(prhs[0]);
    

    
    double z1       = (double)(mxGetScalar(prhs[3]));
    double z2       = (double)(mxGetScalar(prhs[4]));
    double lambda   = (double)(mxGetScalar(prhs[5]));
    double T        = (double)(mxGetScalar(prhs[6]));
    
    orders          = mxGetPr(prhs[7]);
    
    alpha           = (double)(mxGetScalar(prhs[8]));
    beta            = (double)(mxGetScalar(prhs[9]));
    zernikeCouples  =  mxGetPr(prhs[10]);
    rowAber          = mxGetN(prhs[10]);
    
    double NA        = (double)(mxGetScalar(prhs[11]));
    
    OPDx = computeOPD(detectorCoordsX, detectorCoordsY, detectorCoordsZ, alpha, beta, z1, z2, lambda, T, NA,
                        orders, ((long)rowLen) * ((long) colLen), zernikeCouples, rowAber, true);
    OPDy = computeOPD(detectorCoordsX, detectorCoordsY, detectorCoordsZ, alpha, beta, z1, z2, lambda, T, NA,
                        orders, ((long)rowLen) * ((long) colLen), zernikeCouples, rowAber, false);
    
    plhs[0]         = mxCreateDoubleMatrix(rowLen, colLen, mxREAL);
    plhs[1]         = mxCreateDoubleMatrix(rowLen, colLen, mxREAL);
    
    double * OPDPtrx = mxGetPr(plhs[0]);
    double * OPDPtry = mxGetPr(plhs[1]);

    // Reassign computed array to matlab-allocated array
    for (long k = 0 ; k < rowLen * colLen; k++)
       OPDPtrx[k] = OPDx[k];
    
    // Reassign computed array to matlab-allocated array
    for (long k = 0 ; k < rowLen * colLen; k++)
       OPDPtry[k] = OPDy[k];
    
}
