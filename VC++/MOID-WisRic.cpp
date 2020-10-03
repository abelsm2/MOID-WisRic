// MOID-WisRic.cpp
//
// Original Research Paper: http://moid.cbk.waw.pl/orbity/static/MOID.pdf
//
// Implementation of the algorithm published by Wisniowski and Rickman for calculating MOID.  This is
// just a conversion from their Fortran code to visual c++ and is represnted by the function moid_wisric().
// main() just includes some code to run the test cases published by the author to just the correctness of
// the conversion as well as to judge the speed of the algorithm.


#define _USE_MATH_DEFINES

#include <iostream>
#include <math.h>
#include <iomanip>
#include <cassert>
#include <chrono>  // for high_resolution_clock
#include <vector>
#include <algorithm>
#include <random>

// declarations
std::vector<long double> moid_wisric(const std::vector<long double>&, const std::vector<long double>&, const std::vector<long double>&, const std::vector<long double>&, const std::vector<long double>&);

// random number generator for methods like shuffle
static std::mt19937_64 urng(std::chrono::high_resolution_clock::now().time_since_epoch().count());

int main() {

    // This is the number of cases that will be run
    const int numCases(1000);  

    // uniform random number distribution for perturbing variables
    std::uniform_real_distribution<double> rnd(0.8, 1.2);

    std::vector<std::vector<long double>> moid(numCases,std::vector<long double> (3,0));
    std::vector<std::vector<long double>> semiMajorAxis(numCases, std::vector<long double>(2,0));
    std::vector<std::vector<long double>> eccentricity(numCases, std::vector<long double>(2, 0));
    std::vector<std::vector<long double>> inclination(numCases, std::vector<long double>(2, 0));
    std::vector<std::vector<long double>> RAAN(numCases, std::vector<long double>(2, 0));
    std::vector<std::vector<long double>> argPeriapsis(numCases, std::vector<long double>(2, 0));

    //# Results from Research Paper - Angles were kept in degrees here for ease of entry
    //# http://moid.cbk.waw.pl/orbity/static/MOID.pdf
    //                                                 radius-peri      e        i(deg)   Omega(deg)   w(deg)      Actual MOID     D - min
    std::vector<std::vector<long double>> Asteroids { {2.55343183, 0.0777898,  10.58785,  80.35052,  72.14554, 0.13455874348909, 0.1346719 },
                                                      {2.12995319, 0.2313469,  34.84268, 173.12520, 310.03850, 0.00289925623680, 0.0151316 },
                                                      {1.98948966, 0.2552218,  12.97943, 169.90317, 248.22602, 0.07817951779390, 0.0784614 },
                                                      {2.15354370, 0.0882196,   7.13426, 103.89537, 150.08873, 0.08735595371552, 0.0878071 },
                                                      {2.08388391, 0.1905003,   5.36719, 141.60955, 358.80654, 0.14532630925408, 0.1455354 },
                                                      {2.48391159, 0.9543470, 119.29902,  39.00301, 357.90012, 0.26938418933051, 0.2709789 },
                                                      {2.36382356, 0.9006860, 160.41316, 297.34820, 102.45000, 0.54491059333263, 0.5479893 },
                                                      {0.13964163, 0.8901393,  22.23224, 265.28749, 322.11933, 0.70855959609279, 0.7152397 },
                                                      {0.35420623, 0.8363753,  11.68912,  28.13011, 208.66724, 0.03943927946198, 0.1120682 },
                                                      {0.52469070, 0.7715449,  12.56792,   7.25167, 122.30952, 0.18225709092897, 0.2183000 },
                                                      {2.74144856, 0.1153501,   0.00431, 272.90217, 251.43828, 0.14766834758223, 0.1483890 },
                                                      {2.50571901, 0.1924270,   0.01522,  94.14405, 304.71343, 0.00010493251317, 0.0102467 },
                                                      {2.11312640, 0.1215091,   0.02244, 321.26045, 109.96758, 0.00030783183432, 0.0141584 },
                                                      {2.09876663, 0.1543590,   0.02731,  88.64817,  67.91991, 0.00098583168214, 0.0084652 },
                                                      {2.67112178, 0.1328536,   0.02809,  41.39822, 274.65080, 0.20707625146740, 0.2087313 },
                                                      {1.99601821, 0.1875129,   1.26622, 238.06043,  31.32645, 0.00000003815330, 0.0069443 },
                                                      {2.03086844, 0.1653922,   0.66023, 339.21518,  89.47548, 0.00000419348257, 0.0011561 },
                                                      {1.77550824, 0.1928808,   3.43901, 140.55651, 216.20834, 0.00000627704688, 0.0402762 },
                                                      {1.96745453, 0.1837814,   3.69269,  98.95749, 227.52626, 0.00000785853673, 0.0043333 },
                                                      {2.15731280, 0.1007470,   2.91058, 138.77805, 231.93187, 0.00001189165231, 0.0100350 } };
                    

    // target defined in research paper linked above
    std::vector<long double> Target { 2.036, .164, 0, 0, (250.227 * M_PI / 180.0), (2.036 / (1 - .164)) };

    // convert all the angles to radians and calculate the semi-major axis
    for (int i = 0; i < 20; i++) {
        Asteroids[i][2] = Asteroids[i][2] * M_PI / 180.0;
        Asteroids[i][3] = Asteroids[i][3] * M_PI / 180.0;
        Asteroids[i][4] = Asteroids[i][4] * M_PI / 180.0;

        // calculate semi-major axis from the data provided
        Asteroids[i].push_back(Asteroids[i][0] / (1 - Asteroids[i][1]));
    }

    // fill in the 1st 20 test cases with the data above and the rest with slightly random data
    for (int i = 0; i < numCases; i++) {

        // load up our test data with known results
        if (i < 20) {
            semiMajorAxis[i] = { Target[5], Asteroids[i][7] };
            eccentricity[i] = { Target[1], Asteroids[i][1] };
            inclination[i] = { Target[2],  Asteroids[i][2] };
            RAAN[i] = { Target[3], Asteroids[i][3] };
            argPeriapsis[i] = { Target[4], Asteroids[i][4] };
        }

        // this is randomized data, if numTest is higher than 20
        else {

            semiMajorAxis[i] = { Target[5], Asteroids[i % 20][7] * rnd(urng) };
            eccentricity[i] = { Target[1], Asteroids[i % 20][1] * rnd(urng) };
            inclination[i] = { Target[2],  Asteroids[i % 20][2] * rnd(urng) };
            RAAN[i] = { Target[3], Asteroids[i % 20][3] * rnd(urng) };
            argPeriapsis[i] = { Target[4], Asteroids[i % 20][4] * rnd(urng) };
        }
    }

    // Record start time
    auto start = std::chrono::high_resolution_clock::now();

    for (int i = 0; i < numCases; i++) {

        moid[i] = moid_wisric(semiMajorAxis[i], eccentricity[i], argPeriapsis[i], RAAN[i], inclination[i]);
    }

    // Record end time
    auto finish = std::chrono::high_resolution_clock::now();
    
    for (int i = 0; i < 20; i++) {

        // for a debug build, verify that we are getting the correct results to within a small error
        assert((moid[i][0] - Asteroids[i][5]) < 1e-8);

        // for all builds, print out the test cases that match the data above
        std::cout << std::fixed << std::setprecision(14) << "MOID: " << moid[i][0] << "\t" << Asteroids[i][5] 
                                << "\tv(A): " << std::setprecision(5) << std::setw(8) << moid[i][1] 
                                << "\tv(B): " << std::setw(8) << moid[i][2] << std::endl;
    }

    std::chrono::duration<double, std::micro> averageElapsed = (finish - start) / numCases;
    std::cout << std::setprecision(3) << "Average Execution Time: " << averageElapsed.count() << " microseconds\n";
}

// From the initial authors website: http://moid.cbk.waw.pl/orbity/static/MOID.F
//
// current version: 4.0
// version 1.0(July 31th, 2013) 
// version 2.0(August 23th,2013)   works faster, some unnecessary caclulations removed 
// version 3.0(September 5th,2014) handles extremely eccentric orbits, corrections 
//                                 done thanks to comments submitted by Zong-Fu,Sie NCU Taiwan
// version 4.0(July 30th,2017)     removes important bug injected in version 3.0,  
//                                 corrections done thanks to comments submitted by
//                                 Robert Jedicke,University of Hawaii
//--------------------------------
// This program calculates the MOID between two asteroids - Lutetia and Magdalena.
// It uses all ideas and solutions described with details in the paper by 
// T.Wisniowski and H.Rickman "A Fast, Geometric Method for Calculating Accurate
// Minimum Orbit Intersection Distances (MOIDs)" published in 2013 in Acta Astronomica.
// The program is free and may be used without limits as the core of any other program.
// The authors will appreciate for mentioning them, when the program will appear to be useful.

std::vector<long double> moid_wisric(const std::vector<long double>& semiMajorAxis,
                                     const std::vector<long double>& eccentricity, 
                                     const std::vector<long double>& argPeriapsis, 
                                     const std::vector<long double>& RAAN, 
                                     const std::vector<long double>& inclination) {

    // variables used in loops but that need scope outside of those loops
    long double rA, rA2, rB, sintmp, costmp, Bz_sq ;
    long double tmp1, tmp2, step, threshold;
    long double dist, dist_min, dist_oo, trueB_o, trueB_m;
    long double longit, longit_m, longit_o;
    long double tmptrueB[11]{ 0 }, tmplongit[11]{ 0 }, tmpmoid[11]{ 0 };
    long double rAt[4]{ 0 }, rBt[4]{ 0 }, Axt[4]{ 0 }, Ayt[4]{ 0 }, Bxt[4]{ 0 }, Byt[4]{ 0 }, Bzt[4]{ 0 };
    long double moidStart = 1e50;
    int nmax;

    //....orbital parameters of body A
    long double saxisA = semiMajorAxis[0];
    long double eccenA = eccentricity[0];
    long double argpeA = argPeriapsis[0];  // radians
    long double omegaA = RAAN[0];          // radians
    long double incliA = inclination[0];   // radians
    
    //....orbital parameters of body B
    long double saxisB = semiMajorAxis[1];
    long double eccenB = eccentricity[1];
    long double argpeB = argPeriapsis[1];  // radians
    long double omegaB = RAAN[1];          // radians
    long double incliB = inclination[1];   // radians
    
    //  Parameters of the program   
    //  These steps are optimized with respect to speed and reliability
    long double cstep = 0.12;     // scanning step of true anomaly / longitude in radians
    long double stepini = 0.07;   // initial step of first tuning in radians
    long double steptresh = 1e-5; // final step of first tunig(to choose the MOID) in radians
    
    //  This step depends on expected final accuracy
    long double stepmin = 1e-14;  // threshold step of second tuning in radians

    /////////////////////////////// START OF PREPARING THE ORBITS /////////////////////////////

    //....computing parameters of transition matrix c11...c33
    long double c11 = cos(omegaA) * cos(argpeA) - sin(omegaA) * cos(incliA) * sin(argpeA);
    long double c12 = sin(omegaA) * cos(argpeA) + cos(omegaA) * cos(incliA) * sin(argpeA);
    long double c13 = sin(incliA) * sin(argpeA);
    long double c21 = -cos(omegaA) * sin(argpeA) - sin(omegaA) * cos(incliA) * cos(argpeA);
    long double c22 = -sin(omegaA) * sin(argpeA) + cos(omegaA) * cos(incliA) * cos(argpeA);
    long double c23 = sin(incliA) * cos(argpeA);
    long double c31 = sin(incliA) * sin(omegaA);
    long double c32 = -sin(incliA) * cos(omegaA);
    long double c33 = cos(incliA);

    //....calculating new values of Euler angles using transition matrix 
    long double sintmpi = sin(incliB);
    long double costmpi = cos(incliB);
    long double costmpo = cos(omegaB);
    long double sintmpo = sin(omegaB);
    long double costmpa = cos(argpeB);
    long double sintmpa = sin(argpeB);
    long double x1 = costmpo * costmpa - sintmpo * costmpi * sintmpa;
    long double x2 = sintmpo * costmpa + costmpo * costmpi * sintmpa;
    long double x3 = sintmpi * sintmpa;
    long double y1 = -costmpo * sintmpa - sintmpo * costmpi * costmpa;
    long double y2 = -sintmpo * sintmpa + costmpo * costmpi * costmpa;
    long double y3 = sintmpi * costmpa;
    long double z1 = sintmpi * sintmpo;
    long double z2 = -sintmpi * costmpo;
    long double z3 = costmpi;
    long double z1n = c11 * z1 + c12 * z2 + c13 * z3;
    long double z2n = c21 * z1 + c22 * z2 + c23 * z3;
    long double z3n = c31 * z1 + c32 * z2 + c33 * z3;
    long double y3n = c31 * y1 + c32 * y2 + c33 * y3;
    long double x3n = c31 * x1 + c32 * x2 + c33 * x3;
    incliB = atan2(sqrt(z1n * z1n + z2n * z2n), z3n);
    omegaB = -atan2(z1n, -z2n);
    argpeB = -atan2(x3n, y3n);

    //....helpful precalculated values
    costmpo = cos(omegaB);
    sintmpo = sin(omegaB);
    sintmpi = sin(incliB);  // = z1n / sintmpo
    costmpi = z3n;          // = cos(incliB)
    long double sint = sintmpo * costmpi;
    long double cost = costmpo * costmpi;
    long double radA = saxisA * (1 - eccenA * eccenA);
    long double radB = saxisB * (1 - eccenB * eccenB);
    /////////////////////////////// END OF PREPARING THE ORBITS /////////////////////////////

    ///////////////////////////////      START OF SCANNING      /////////////////////////////
    // This tool yields a preliminary approach to the minima of the distance function.
    // By scanning one full revolution of meridional plane we look for the local minima. 
    
    //......initial parameters                          
    long double trueB = -2 * cstep;
    long double moid = moidStart;
    long double dist_o = moidStart;  // something big
    tmpmoid[1] = moidStart;
    tmpmoid[2] = moidStart;
    tmpmoid[3] = moidStart;
    tmpmoid[4] = moidStart;
        
    //.....Looking for the minima with rotating meridional plane
    
    //.......a)at first we calculate the coordinates of two additional positions of the plane to create first triplet
    for (int iii = 1; iii <= 2; iii++) {

        rB = radB / (1 + eccenB * cos(trueB)); // compute the radius for B
        sintmp = sin(trueB + argpeB);
        costmp = cos(trueB + argpeB);
        Bz_sq = sintmpi * sintmp;
        Bz_sq = Bz_sq * Bz_sq; // square of Z - coordinate for B
        longit = atan2(sintmpo * costmp + sintmp * cost, costmpo * costmp - sintmp * sint); // compute the longitude for A
        tmp2 = eccenA * cos(longit); // temporary value
        rA = radA / (1 + tmp2); // compute the radius for A(two possibilities)
        rA2 = radA / (1 - tmp2);
        tmp1 = rB * sqrt(1 - Bz_sq); // temporary value
        
        if (abs(tmp1 - rA) > abs(tmp1 + rA2)) {
            rA = rA2;
            longit = longit - M_PI; // the second possibility gives smaller distance
            tmp1 = tmp1 + rA2;
        }

        else {
            tmp1 = tmp1 - rA;
        }
         
        dist = rB * rB * Bz_sq + tmp1 * tmp1; // square of the distance A - B

        if (iii == 1) {
            dist_oo = dist;
        }

        else {
            dist_o = dist;
            trueB_o = trueB;
            longit_o = longit;
        }
        
        trueB = trueB + cstep;
    }
    
    //.......b)now we scan doing one full revolution of the meridional plane  
    nmax = 0; // counts the minima
    dist_min = dist;
    
    // loop for true anomaly of B
    while (trueB < (2 * M_PI + cstep)) {                    

        rB = radB / (1 + eccenB * cos(trueB)); // compute the radius for B
        sintmp = sin(trueB + argpeB);
        costmp = cos(trueB + argpeB);
        Bz_sq = sintmpi * sintmp;
        Bz_sq = Bz_sq * Bz_sq; // square of Z - coordinate for B
        longit = atan2(sintmpo * costmp + sintmp * cost, costmpo * costmp - sintmp * sint); // compute the longitude for A
        tmp2 = eccenA * cos(longit); // temporary value
        rA = radA / (1 + tmp2); // compute the radius for A(two possibilities)
        rA2 = radA / (1 - tmp2);
        tmp1 = rB * sqrt(1 - Bz_sq); // temporary value

        if (abs(tmp1 - rA) > abs(tmp1 + rA2)) {
            rA = rA2;
            longit = longit - M_PI; // the second possibility gives smaller distance
            tmp1 = tmp1 + rA2;
        }

        else {
            tmp1 = tmp1 - rA;
        }
            
        dist = rB * rB * Bz_sq + tmp1 * tmp1; // square of the distance A - B

        // the minimum was found
        if ((dist_o <= dist) && (dist_o <= dist_oo)) {    
            nmax += 1;
            tmptrueB[nmax] = trueB_o;
            tmplongit[nmax] = longit_o;
            tmpmoid[nmax] = dist_o;
        }
        
        if (dist_min > dist) {
            dist_min = dist;
        }

        dist_oo = dist_o;
        trueB_o = trueB;
        longit_o = longit;
        dist_o = dist;
        trueB = trueB + cstep;
    }
    
    ///////////////////////////////      END OF SCANNING      /////////////////////////////

    // Enclosed the water procedure and tuning procedure in a while loop to maintain
    // the flow control of the original Fortran code but without the GOTO statements
    // so that we can run the water procedure in the case that two minima are very 
    // close to each other during the tuning phase.

    bool done = false; // flow control to handle the GOTO statement in the original Fortran code

    while (!done) {

        done = true; // the only time this will become false is when checking if two minimum are very close together

        //////////////////////////// "WATER" PROCEDURE" //////////////////////////////////
        // In case only one minimum is detected we take a special care
        // to avoid the risk of missing the minima. 
        // Instead of starting the tuning with one detected minimum,
        // we start it with four positions of the meridional plane, evenly
        // distributed along the inclined orbit. The later tuning procedure
        // moves the points along the orbits similarly as water droplets
        // are leading by the gravity to the points of minimum height, so
        // we called this phase "water procedure". It slightly slows down
        // the calculations, but minimizes the risk of missing the MOID.
        // With "water procedure" the speed is 9-12 seconds per 100,000 MOIDs, risk of missing <1E-6
        // Without "water procedure" the speed is <9 seconds per 100,000 MOIDs, risk of missing about 3E-5
        // (speed measured on fast single CPU core)

        // only one minimum was detected
        if (nmax < 2) {
            nmax = 4;

            for (int iii = 1; iii <= 4; iii++) {
                tmptrueB[iii] = (.25 + .5 * iii) * M_PI;  // evenly distributed points
                sintmp = sin(tmptrueB[iii] + argpeB);
                costmp = cos(tmptrueB[iii] + argpeB);
                tmplongit[iii] = atan2(sintmpo * costmp + sintmp * cost, costmpo * costmp - sintmp * sint); // compute the longitude for A
                tmpmoid[iii] = moidStart; // something big
            }
        }
        ////////////////////////// END OF "WATER" PROCEDURE" ///////////////////////////////

        ///////////////////////////////     START OF PARALLEL TUNING      /////////////////////////////
        // After the scanning phase, we typically have a few minima on a meridional plane. 
        // The goal of the tuning procedure is to move objects separately along their orbits
        // in order to find the smallest possible distance between them, which is no longer a meridional distance.      
          
        for (int jjj = 1; jjj <= nmax + 1; jjj++) {

            if (jjj <= nmax) {
                moid = tmpmoid[jjj];
                trueB_m = tmptrueB[jjj];
                longit_m = tmplongit[jjj];
                step = stepini;
                threshold = steptresh;
            }

            else {

                if (nmax == 2) {

                    // in case of two minima are very close to each other(<1E-4) run the "water procedure"
                    if (abs(tmpmoid[1] - tmpmoid[2]) < 1e-4) {
                        nmax = 1;// make sure we can get into the water procedure
                        done = false;// set false so we start at the top of the loop again
                        break; // get out of the for - loop and start over
                    }

                    else if (tmpmoid[1] < moid) {
                        moid = tmpmoid[1];
                        trueB_m = tmptrueB[1];
                        longit_m = tmplongit[1];
                    }
                }

                else {

                    for (int iii = 1; iii <= nmax - 1; iii++) { // the choice of moids for final tuning

                        if (tmpmoid[iii] < moid) {
                            moid = tmpmoid[iii];
                            trueB_m = tmptrueB[iii];
                            longit_m = tmplongit[iii];
                        }
                    }
                }

                step = 2 * stepini;  // initial step for final tuning
                threshold = stepmin;  // terminal step for final tuning
            }
                
            rBt[2] = radB / (1 + eccenB * cos(trueB_m));
            sintmp = sin(trueB_m + argpeB);
            costmp = cos(trueB_m + argpeB);
            Bxt[2] = costmpo * costmp - sintmp * sint;
            Byt[2] = sintmpo * costmp + sintmp * cost;
            Bzt[2] = sintmpi * sintmp;
            rAt[2] = radA / (1 + eccenA * cos(longit_m));
            Axt[2] = cos(longit_m);
            Ayt[2] = sin(longit_m);
            bool aleft{ true };
            bool aright{ true };
            bool bleft{ true };
            bool bright{ true };

            while (step >= threshold) {
                int lpoints{ 0 };
                int j1min{ 1 };
                int j1max{ 3 };
                int i1min{ 1 };
                int i1max{ 3 };
                bool calc1{ false };
                bool calc2{ false };
                bool calc3{ false };
                bool calc4{ false };

                if (bleft) {
                    rBt[1] = radB / (1 + eccenB * cos(trueB_m - step));
                    sintmp = sin(trueB_m - step + argpeB);
                    costmp = cos(trueB_m - step + argpeB);
                    Bxt[1] = costmpo * costmp - sintmp * sint;
                    Byt[1] = sintmpo * costmp + sintmp * cost;
                    Bzt[1] = sintmpi * sintmp;
                    lpoints = lpoints + 1;
                }

                if (bright) {
                    rBt[3] = radB / (1 + eccenB * cos(trueB_m + step));
                    sintmp = sin(trueB_m + step + argpeB);
                    costmp = cos(trueB_m + step + argpeB);
                    Bxt[3] = costmpo * costmp - sintmp * sint;
                    Byt[3] = sintmpo * costmp + sintmp * cost;
                    Bzt[3] = sintmpi * sintmp;
                    lpoints = lpoints + 1;
                }

                if (aleft) {
                    rAt[1] = radA / (1 + eccenA * cos(longit_m - step));
                    Axt[1] = cos(longit_m - step);
                    Ayt[1] = sin(longit_m - step);
                    lpoints = lpoints + 1;
                }

                if (aright) {
                    rAt[3] = radA / (1 + eccenA * cos(longit_m + step));
                    Axt[3] = cos(longit_m + step);
                    Ayt[3] = sin(longit_m + step);
                    lpoints = lpoints + 1;
                }

                int j1_t = 2;
                int i1_t = 2;

                if (lpoints == 1) {

                    if (aleft) i1max = 1;
                    if (aright) i1min = 3;
                    if (bleft) j1max = 1;
                    if (bright) j1min = 3;
                }

                if (lpoints == 2) {

                    if (aleft && bright) calc1 = true;
                    if (aleft && bleft) calc2 = true;
                    if (aright && bright) calc3 = true;
                    if (aright && bleft) calc4 = true;
                }

                for (int j1 = j1min; j1 <= j1max; j1++) {

                    for (int i1 = i1min; i1 <= i1max; i1++) {

                        if (lpoints == 2) {

                            if (i1 != 1) {
                                if (((j1 != 3) && calc1) || ((j1 != 1) && calc2)) continue;
                            }

                            if (i1 != 3) {
                                if (((j1 != 3) && calc3) || ((j1 != 1) && calc4)) continue;
                            }
                        }

                        if ((i1 == 2) && (j1 == 2)) continue;

                        long double Dx = rBt[j1] * Bxt[j1] - rAt[i1] * Axt[i1];
                        long double Dy = rBt[j1] * Byt[j1] - rAt[i1] * Ayt[i1];
                        long double Dz = rBt[j1] * Bzt[j1];
                        dist = (Dx * Dx + Dy * Dy + Dz * Dz);

                        if (dist < moid) {
                            moid = dist;
                            j1_t = j1;
                            i1_t = i1;
                        }
                    }
                }


                if ((j1_t != 2) || (i1_t != 2)) {

                    aleft = false;
                    aright = false;
                    bleft = false;
                    bright = false;

                    if (i1_t != 2) {

                        if (i1_t == 1) {
                            aleft = true;
                            longit_m = longit_m - step;
                            rAt[3] = rAt[2];
                            Axt[3] = Axt[2];
                            Ayt[3] = Ayt[2];
                            rAt[2] = rAt[1];
                            Axt[2] = Axt[1];
                            Ayt[2] = Ayt[1];
                        }

                        else {
                            aright = true;
                            longit_m = longit_m + step;
                            rAt[1] = rAt[2];
                            Axt[1] = Axt[2];
                            Ayt[1] = Ayt[2];
                            rAt[2] = rAt[3];
                            Axt[2] = Axt[3];
                            Ayt[2] = Ayt[3];
                        }
                    }

                    if (j1_t != 2) {

                        if (j1_t == 1) {
                            bleft = true;
                            trueB_m = trueB_m - step;
                            rBt[3] = rBt[2];
                            Bxt[3] = Bxt[2];
                            Byt[3] = Byt[2];
                            Bzt[3] = Bzt[2];
                            rBt[2] = rBt[1];
                            Bxt[2] = Bxt[1];
                            Byt[2] = Byt[1];
                            Bzt[2] = Bzt[1];
                        }

                        else {
                            bright = true;
                            trueB_m = trueB_m + step;
                            rBt[1] = rBt[2];
                            Bxt[1] = Bxt[2];
                            Byt[1] = Byt[2];
                            Bzt[1] = Bzt[2];
                            rBt[2] = rBt[3];
                            Bxt[2] = Bxt[3];
                            Byt[2] = Byt[3];
                            Bzt[2] = Bzt[3];
                        }
                    }
                }

                else {
                    aleft = true;
                    aright = true;
                    bleft = true;
                    bright = true;
                    step = step * .15; // 0.15 is optimal value
                }
            }

            if (jjj <= nmax) {
                tmpmoid[jjj] = moid;
                tmptrueB[jjj] = trueB_m;
                tmplongit[jjj] = longit_m;
            }
        }
    //###### END OF PARALLEL TUNING ##########################################
    }

    // added ability to return the true anomaly where the distance is minimum in the range -pi to +pi
    long double trueAnomalyA, trueAnomalyB;
    
    if (longit_m > M_PI) {
        trueAnomalyA = 2 * M_PI - longit_m;
    }
     
    else if (longit_m < -M_PI) {
        trueAnomalyA = 2 * M_PI + longit_m;
    }
    
    else {
        trueAnomalyA = -longit_m;
    }
            
    if (trueB_m > M_PI) {
        trueAnomalyB = 2 * M_PI - trueB_m;
    }
     
    else if (trueB_m < -M_PI) {
        trueAnomalyB = -2 * M_PI + trueB_m;
    }

    else {
        trueAnomalyB = -trueB_m;
    }
        
    return std::vector<long double> {sqrt(moid), trueAnomalyA, trueAnomalyB}; // we dealt with squares
}



