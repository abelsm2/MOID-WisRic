## function [moid, dist_min_sq] = MOID_WisRic (semiMajorAxis, eccentricity, argPeriapsis, RAAN, inclination)
##
##    Note on Inputs:  All inputs are a column vector that represents the orbital parameters for two orbits.  The main
##                     body (Target) orbit is in the first row and the object of interest is in the second row.
##
##    Input
##        semiMajorAxis   Column vector containing semi-major axes for both orbits using any units desired.
##        eccentricity    Column vector containing the eccentricity for both orbits.
##        argPeriapsis    Column vector containing the Argument of Periapsis for both orbits, in degrees.
##        RAAN            Column vector containing the Right Ascension of the Ascending Node for both orbits, in degrees.
##        inclination     Column vector containing the inclination of both orbits, in degrees.
##
##    Returns
##        moid            The minimum orbit intersection distance.  Units are the same as those used for SemiMajorAxis
##        dist_min_sq     This is the SQARE of the minimum distance found during the scanning phase of the algorithm.
##                        Units are the same as those used for SemiMajor Axis.  This is largely only used for initial
##                        error checking to compare results with those of the original authors.
##
## Created for Octave: 2020-09-25
##
## From the initial authors website: http://moid.cbk.waw.pl/orbity/static/MOID.F
##
## current version: 4.0
## version 1.0(July 31th, 2013) 
## version 2.0(August 23th,2013)   works faster, some unnecessary caclulations removed 
## version 3.0(September 5th,2014) handles extremely eccentric orbits, corrections 
##                                 done thanks to comments submitted by Zong-Fu,Sie NCU Taiwan
## version 4.0(July 30th,2017)     removes important bug injected in version 3.0,  
##                                 corrections done thanks to comments submitted by
##                                 Robert Jedicke,University of Hawaii
##--------------------------------
## This program calculates the MOID between two asteroids - Lutetia and Magdalena.
## It uses all ideas and solutions described with details in the paper by 
## T.Wisniowski and H.Rickman "A Fast, Geometric Method for Calculating Accurate
## Minimum Orbit Intersection Distances (MOIDs)" published in 2013 in Acta Astronomica.
## The program is free and may be used without limits as the core of any other program.
## The authors will appreciate for mentioning them, when the program will appear to be useful.
##

function [moid, dist_min_sq] = MOID_WisRic (semiMajorAxis, eccentricity, argPeriapsis, RAAN, inclination)

##....orbital parameters of body A - asteroid Magdalena, angles in rads
    saxisA = semiMajorAxis(1); % semiaxis   
    eccenA = eccentricity(1);  % eccenticity
    argpeA = argPeriapsis(1);  % argument of perihelion
    omegaA = RAAN(1);          % longitude of ascending node
    incliA = inclination(1);   % inclination
      
##....orbital parameters of body B - asteroid Lutetia, angles in rads
    saxisB = semiMajorAxis(2); % semiaxis
    eccenB = eccentricity(2);  % eccenticity
    argpeB = argPeriapsis(2);  % argument of perihelion
    omegaB = RAAN(2);          % longitude of ascending node
    incliB = inclination(2);   % inclination

##  Parameters of the program   
##  These steps are optimized with respect to speed and reliability
    cstep     = 0.12;  % scanning step of true anomaly/longitude in radians
    stepini   = 0.07;  % initial step of first tuning in radians
    steptresh = 1e-5;  % final step of first tunig (to choose the MOID) in radians

##  This step depends on expected final accuracy
    stepmin   = 1e-14; % threshold step of second tuning in radians
    
################### START OF PREPARING THE ORBITS ##########################

##....computing parameters of transition matrix c11...c33
    c11 = cos(omegaA) * cos(argpeA) - sin(omegaA) * cos(incliA) * sin(argpeA);
    c12 = sin(omegaA) * cos(argpeA) + cos(omegaA) * cos(incliA) * sin(argpeA);
    c13 = sin(incliA) * sin(argpeA);
    c21 = -cos(omegaA) * sin(argpeA) - sin(omegaA) * cos(incliA) * cos(argpeA);
    c22 = -sin(omegaA) * sin(argpeA) + cos(omegaA) * cos(incliA) * cos(argpeA);
    c23 = sin(incliA) * cos(argpeA);
    c31 = sin(incliA) * sin(omegaA);
    c32 = -sin(incliA) * cos(omegaA);
    c33 = cos(incliA);
    
##....calculating new values of Euler angles using transition matrix 
    sintmpi = sin(incliB);
    costmpi = cos(incliB);
    costmpo = cos(omegaB);
    sintmpo = sin(omegaB);
    costmpa = cos(argpeB);
    sintmpa = sin(argpeB);
    x1 = costmpo * costmpa - sintmpo * costmpi * sintmpa;
    x2 = sintmpo * costmpa + costmpo * costmpi * sintmpa;
    x3 = sintmpi * sintmpa;
    y1 = -costmpo * sintmpa - sintmpo * costmpi * costmpa;
    y2 = -sintmpo * sintmpa + costmpo * costmpi * costmpa;
    y3 = sintmpi * costmpa;
    z1 = sintmpi * sintmpo;
    z2 = -sintmpi * costmpo;
    z3 = costmpi;
    z1n = c11 * z1 + c12 * z2 + c13 * z3;
    z2n = c21 * z1 + c22 * z2 + c23 * z3;
    z3n = c31 * z1 + c32 * z2 + c33 * z3;
    y3n = c31 * y1 + c32 * y2 + c33 * y3;
    x3n = c31 * x1 + c32 * x2 + c33 * x3;
    incliB = atan2(sqrt(z1n * z1n + z2n * z2n), z3n);
    omegaB = -atan2(z1n, -z2n);
    argpeB = -atan2(x3n, y3n);
    
##....helpful precalculated values
    costmpo = cos(omegaB);
    sintmpo = sin(omegaB);
    sintmpi = sin(incliB); % = z1n / sintmpo
    costmpi = z3n;         % = cos(incliB)
    sint = sintmpo * costmpi;
    cost = costmpo * costmpi;
    radA = saxisA * (1 - eccenA * eccenA);
    radB = saxisB * (1 - eccenB * eccenB);
###################### END OF PREPARING THE ORBITS #########################


########################### START OF SCANNING ##############################
## This tool yields a preliminary approach to the minima of the distance function.
## By scanning one full revolution of meridional plane we look for the local minima. 

##......initial parameters                          
    trueB = -2 * cstep;
    moid = 1e6;
    dist_o = 1E6;  # something big
    tmpmoid(1) = 1e6;
    tmpmoid(2) = 1e6;
    tmpmoid(3) = 1e6;
    tmpmoid(4) = 1e6;
    
##.....Looking for the minima with rotating meridional plane

##.......a)at first we calculate the coordinates of two additional positions of the plane to create first triplet
    for iii = 1 : 2 
        rB = radB / (1 + eccenB * cos(trueB));      #compute the radius for B
        sintmp = sin(trueB + argpeB);
        costmp = cos(trueB + argpeB);
        Bz_sq = sintmpi * sintmp;
        Bz_sq = Bz_sq * Bz_sq;                      # square of Z-coordinate for B
        longit = atan2(sintmpo * costmp + sintmp * cost, costmpo * costmp - sintmp * sint); # compute the longitude for A        
        tmp2 = eccenA * cos(longit);                # temporary value         
        rA = radA / (1 + tmp2);                     # compute the radius for A (two possibilities)
        rA2 = radA / (1 - tmp2); 
        tmp1 = rB * sqrt(1 - Bz_sq);                # temporary value
        
        if (abs(tmp1 - rA) > abs(tmp1 + rA2))
            rA = rA2;
            longit = longit - pi;                   # the second possibility gives smaller distance
            tmp1 = tmp1 + rA2;
            
        else  
            tmp1 = tmp1 - rA;
            
        endif
        
        dist = rB * rB * Bz_sq + tmp1 * tmp1;       # square of the distance A-B
        
        if (iii == 1)
          dist_oo = dist;
          
        else
          dist_o = dist;
          trueB_o = trueB;
          longit_o = longit;
          
        endif
        
        trueB = trueB + cstep;
    end
  
##.......b)now we scan doing one full revolution of the meridional plane  
    nmax = 0; # counts the minima
    dist_min_sq = dist;
    
    while (trueB < (2 * pi + cstep))                    # loop for true anomaly of B           
        rB = radB / (1 + eccenB * cos(trueB));          # compute the radius for B
        sintmp = sin(trueB + argpeB);
        costmp = cos(trueB + argpeB);
        Bz_sq = sintmpi * sintmp;
        Bz_sq = Bz_sq * Bz_sq;                          # square of Z-coordinate for B
        longit = atan2(sintmpo * costmp + sintmp * cost, costmpo * costmp - sintmp * sint);     # compute the longitude for A        
        tmp2 = eccenA * cos(longit);                    # temporary value         
        rA = radA / (1 + tmp2);                         # compute the radius for A (two possibilities)
        rA2 = radA / (1 - tmp2); 
        tmp1 = rB * sqrt(1 - Bz_sq);                    # temporary value
        
        if (abs(tmp1 - rA) > abs(tmp1 + rA2))
            rA = rA2;
            longit = longit - pi;                       # the second possibility gives smaller distance
            tmp1 = tmp1 + rA2;
            
        else
            tmp1 = tmp1 - rA;
            
        endif
        
        dist = rB * rB * Bz_sq + tmp1 * tmp1;           # square of the distance A-B            
        
        if ((dist_o <= dist) && (dist_o <= dist_oo))    # the minimum was found        
            nmax = nmax + 1;
            tmptrueB(nmax) = trueB_o;
            tmplongit(nmax) = longit_o;
            tmpmoid(nmax) = dist_o;
            
        endif
        
        if (dist_min_sq > dist) 
            dist_min_sq = dist;
            
        endif
            
        dist_oo = dist_o;
        trueB_o = trueB;
        longit_o = longit;
        dist_o = dist;
        trueB = trueB + cstep;
    endwhile
########################### END OF SCANNING ##############################


    ## Enclosed the water procedure and tuning procedure in a while loop to maintain
    ## the flow control of the original Fortran code but without the GOTO statements
    ## so that we can run the water procedure in the case that two minima are very 
    ## close to each other during the tuning phase.
    
    done = false; % flow control
    
    while (!done)
                
        done = true; % the only time this will become false if when checking if two minimum are very close together

        ### "WATER" PROCEDURE ####################################################
        ### In case only one minimum is detected we take a special care
        ### to avoid the risk of missing the minima. 
        ### Instead of starting the tuning with one detected minimum,
        ### we start it with four positions of the meridional plane, evenly
        ### distributed along the inclined orbit. The later tuning procedure
        ### moves the points along the orbits similarly as water droplets
        ### are leading by the gravity to the points of minimum height, so
        ### we called this phase "water procedure". It slightly slows down
        ### the calculations, but minimizes the risk of missing the MOID.
        ### With "water procedure" the speed is 9-12 seconds per 100,000 MOIDs, risk of missing <1E-6
        ### Without "water procedure" the speed is <9 seconds per 100,000 MOIDs, risk of missing about 3E-5
        ### (speed measured on fast single CPU core)

        if (nmax < 2) # only one minimum was detected
            nmax = 4;
            
            for iii = 1:4
                tmptrueB(iii) = (.25 + .5 * iii) * pi;  # evenly distributed points      
                sintmp = sin(tmptrueB(iii) + argpeB);
                costmp = cos(tmptrueB(iii) + argpeB);
                tmplongit(iii) = atan2(sintmpo * costmp + sintmp * cost, costmpo * costmp - sintmp * sint); # compute the longitude for A
                tmpmoid(iii) = 1e6; # something big
            end
        endif
        ###### END OF "WATER" PROCEDURE ##########################################

        ###### START OF PARALLEL TUNING ##########################################  
        ### After the scanning phase, we typically have a few minima on a meridional plane. 
        ### The goal of the tuning procedure is to move objects separately along their orbits
        ### in order to find the smallest possible distance between them, which is no longer a meridional distance.      
        for jjj = 1 : nmax + 1    
          
            if (jjj <= nmax)
                moid = tmpmoid(jjj);
                trueB_m = tmptrueB(jjj);
                longit_m = tmplongit(jjj);
                step = stepini;
                threshold = steptresh;
                
            else
                
                if (nmax == 2)
                    
                    ### in case of two minima are very close to each other(<1E-4) run the "water procedure"
                    if (abs(tmpmoid(1) - tmpmoid(2)) < 1e-4) 
                        nmax = 1; % make sure we can get into the water procedure
                        done = false; % set false so we start at the top of the loop again
                        break; % get out of the for-loop and start over
            
                    else
                   
                        if (tmpmoid(1) < moid)
                            moid = tmpmoid(1);
                            trueB_m = tmptrueB(1);
                            longit_m = tmplongit(1);
                        endif
                    endif
                 
                else  
                    for iii = 1 : nmax - 1;  # the choice of moids for final tuning
                        if (tmpmoid(iii) < moid)
                            moid = tmpmoid(iii);
                            trueB_m = tmptrueB(iii);
                            longit_m = tmplongit(iii);
                        endif
                    end
                endif
                
                step = 2 * stepini;  # initial step for final tuning
                threshold = stepmin;  # terminal step for final tuning
            endif
            
            rBt(2) = radB / (1 + eccenB * cos(trueB_m));
            sintmp = sin(trueB_m + argpeB);
            costmp = cos(trueB_m + argpeB);
            Bxt(2) = costmpo * costmp - sintmp * sint;
            Byt(2) = sintmpo * costmp + sintmp * cost;
            Bzt(2) = sintmpi * sintmp;
            rAt(2) = radA / (1 + eccenA * cos(longit_m));
            Axt(2) = cos(longit_m);
            Ayt(2) = sin(longit_m);
            aleft = true;
            aright = true;
            bleft = true;
            bright = true;
            
            while (step >= threshold)
                lpoints = 0;
                j1min = 1;
                j1max = 3;
                i1min = 1;
                i1max = 3;
                calc1 = false;
                calc2 = false; 
                calc3 = false;
                calc4 = false;
                
                if (bleft)
                    rBt(1) = radB / (1 + eccenB * cos(trueB_m - step));
                    sintmp = sin(trueB_m - step + argpeB);
                    costmp = cos(trueB_m - step + argpeB);
                    Bxt(1) = costmpo * costmp - sintmp * sint;
                    Byt(1) = sintmpo * costmp + sintmp * cost;
                    Bzt(1) = sintmpi * sintmp;
                    lpoints = lpoints + 1;
                endif  
                
                if (bright)
                    rBt(3) = radB / (1 + eccenB * cos(trueB_m + step));
                    sintmp = sin(trueB_m + step + argpeB);
                    costmp = cos(trueB_m + step + argpeB);
                    Bxt(3) = costmpo * costmp - sintmp * sint;
                    Byt(3) = sintmpo * costmp + sintmp * cost;
                    Bzt(3) = sintmpi * sintmp;
                    lpoints = lpoints + 1;
                endif  
                
                if (aleft)
                    rAt(1) = radA / (1 + eccenA * cos(longit_m - step));
                    Axt(1) = cos(longit_m - step);
                    Ayt(1) = sin(longit_m - step);
                    lpoints = lpoints + 1;
                endif  
                
                if (aright)
                    rAt(3) = radA / (1 + eccenA * cos(longit_m + step));
                    Axt(3) = cos(longit_m + step);
                    Ayt(3) = sin(longit_m + step);
                    lpoints = lpoints + 1;
                endif
                
                j1_t = 2;
                i1_t = 2;
                
                if (lpoints == 1)
                    
                    if (aleft) 
                        i1max = 1;
                    endif         
                    
                    if (aright) 
                        i1min = 3;
                    endif
                      
                    if (bleft) 
                        j1max = 1;
                    endif
                    
                    if (bright) 
                        j1min = 3;
                    endif
                endif  
                
                if (lpoints == 2)
                    
                    if (aleft && bright)
                        calc1 = true;
                    endif

                    if (aleft && bleft)
                        calc2 = true;
                    endif

                    if (aright && bright)
                        calc3 = true;
                    endif

                    if (aright && bleft) 
                        calc4 = true;
                    endif
                endif  
                
                for j1 = j1min : j1max  
                    for i1 = i1min : i1max
                        
                        if (lpoints == 2)
                            if (i1 != 1)
                                if (((j1 != 3) && calc1) || ((j1 != 1) && calc2)) 
                                    continue;
                                endif
                            endif
                            
                            if (i1 != 3)
                                if (((j1 != 3) && calc3) || ((j1 != 1) && calc4)) 
                                    continue;
                                endif
                            endif    
                        endif  
                        
                        if ((i1 == 2) && (j1 == 2)) 
                            continue;
                        endif
                        
                        Dx = rBt(j1) * Bxt(j1) - rAt(i1) * Axt(i1);
                        Dy = rBt(j1) * Byt(j1) - rAt(i1) * Ayt(i1);
                        Dz = rBt(j1) * Bzt(j1);
                        dist = (Dx * Dx + Dy * Dy + Dz * Dz);

                        if (dist < moid)
                            moid = dist;
                            j1_t = j1;
                            i1_t = i1;
                        endif
                    end
                end
                
                
                if ((j1_t != 2) || (i1_t != 2))
                    
                    aleft = false;
                    aright = false;
                    bleft = false;
                    bright = false;
                    
                    if (i1_t != 2)
                        if (i1_t == 1)
                            aleft = true;
                            longit_m = longit_m - step;
                            rAt(3) = rAt(2);
                            Axt(3) = Axt(2);
                            Ayt(3) = Ayt(2);
                            rAt(2) = rAt(1);
                            Axt(2) = Axt(1);
                            Ayt(2) = Ayt(1);
                        else
                            aright = true;
                            longit_m = longit_m + step;
                            rAt(1) = rAt(2);
                            Axt(1) = Axt(2);
                            Ayt(1) = Ayt(2);
                            rAt(2) = rAt(3);
                            Axt(2) = Axt(3);
                            Ayt(2) = Ayt(3);
                        endif
                    endif  
                    
                  
                    if (j1_t != 2)
                        if (j1_t == 1)
                            bleft = true;
                            trueB_m = trueB_m - step;
                            rBt(3) = rBt(2);
                            Bxt(3) = Bxt(2);
                            Byt(3) = Byt(2);
                            Bzt(3) = Bzt(2);
                            rBt(2) = rBt(1);
                            Bxt(2) = Bxt(1);
                            Byt(2) = Byt(1);
                            Bzt(2) = Bzt(1);
                        else                                       
                            bright = true;
                            trueB_m = trueB_m + step;
                            rBt(1) = rBt(2);
                            Bxt(1) = Bxt(2);
                            Byt(1) = Byt(2);
                            Bzt(1) = Bzt(2);
                            rBt(2) = rBt(3);
                            Bxt(2) = Bxt(3);
                            Byt(2) = Byt(3);
                            Bzt(2) = Bzt(3);
                        endif 
                    endif
                  
                else
                    aleft = true;
                    aright = true;
                    bleft = true;
                    bright = true;
                    step = step * .15; # 0.15 is optimal value
                endif  
            endwhile
            
            if (jjj <= nmax)
                tmpmoid(jjj) = moid;
                tmptrueB(jjj) = trueB_m;
                tmplongit(jjj) = longit_m;
            endif
        end
        ###### END OF PARALLEL TUNING ##########################################
    endwhile
    
    moid = sqrt(moid); # we dealt with squares
    
endfunction
