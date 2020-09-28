# MOID
Implementation of Minimum Orbit Intersection Distance (MOID) Algorithm created by T. Wisniowski and H. Rickman.

## BACKGROUND
I found some Fortran code created by various researchers for the calculation of MOID.  This repo is simply a conversion of those codes into other languages that I found useful along with tests to show that the code is working properly.  As you can see below, I was able to get pretty decent matches to the values published by the original authors with the differences possibly due to a more precise value of Pi being used within Octave and C++.

I did not make any real attempt to optimize the code or use more modern programming style aside from what it took to convert the fortran over to working code so it is implemented as a single function simialr to the original.

* Original research paper: http://moid.cbk.waw.pl/orbity/static/MOID.pdf
* Original Fotran Code: http://moid.cbk.waw.pl/orbity/static/MOID.F

## GNU OCTAVE RESULTS
Running the MOID_Test script on Octave will execute 20 tests which should result in the table below.  This has not been tested in Matlab but should probably work there as well.  On a Ryzen R7 3800XT the time Elapsed for 20 test cases: 0.789262 seconds (Average of ~40 ms per iteration).

Test  | MOID Calc | MOID Expected | Dmin Calc | Dmin Expected 
 ------|----------------|-------------------------|------------------|------------------
1 |      0.13455874619444 |       0.13455874348909  |      0.13467195  |    0.13467190
2 |      0.00289925626282 |       0.00289925623680  |      0.01513159   |   0.01513160
3 |      0.07817951806849 |       0.07817951779390  |      0.07846140   |   0.07846140
4 |      0.08735595327857 |       0.08735595371552  |      0.08780712   |   0.08780710
5 |      0.14532630845989 |       0.14532630925408  |      0.14553541  |    0.14553540
6 |      0.26938418767873 |       0.26938418933051  |      0.27097889  |    0.27097890
7 |      0.54491059218717 |       0.54491059333263  |      0.54798929  |    0.54798930
8 |      0.70855958463834 |       0.70855959609279  |      0.71523969  |    0.71523970
9 |      0.03943927452247 |       0.03943927946198  |      0.11206812  |    0.11206820
10|      0.18225709316049 |       0.18225709092897  |      0.21829996  |    0.21830000
11|      0.14766834353602 |       0.14766834758223  |      0.14838896  |    0.14838900
12|      0.00010493251424 |       0.00010493251317  |      0.01024673  |    0.01024670
13|      0.00030783183885 |       0.00030783183432  |      0.01415839  |    0.01415840
14|      0.00098583168085 |       0.00098583168214  |      0.00846524  |    0.00846520
15|      0.20707624718093 |       0.20707625146740  |      0.20873132  |    0.20873130
16 |     0.00000003860552 |       0.00000003815330  |      0.00694430  |    0.00694430
17 |     0.00000419364072 |       0.00000419348257  |      0.00115612  |    0.00115610
18 |     0.00000627750835 |       0.00000627704688  |      0.04027618  |    0.04027620
19 |     0.00000785937722 |       0.00000785853673  |      0.00433327  |    0.00433330
20 |     0.00001189234779 |       0.00001189165231  |      0.01003500  |    0.01003500

## VC++ RESULTS

On a Ryzen R7 3800XT the average Execution Time is ~14 microseconds.  This was tested by running 10e6 cases and perturbing the original data slightly to generate 10e6 unique orbits similar to the ones provided in the research paper.  This also had the effect of testing hyperbolic orbits when the eccentricity went over 1.0.  While this was done to test the speed of the algorithm, it is unknown if the actual calculated valued represent the correct value.

Using the same 20 test cases as those contained within the research paper, the results are presented below:

MOID Calc | MOID Expected |
 ------|---------------
0.13455874619444      |  0.13455874348909
0.00289925626282      |  0.00289925623680
0.07817951806849      |  0.07817951779390
0.08735595327857      |  0.08735595371552
0.14532630845989      |  0.14532630925408
0.26938418767873      |  0.26938418933051
0.54491059218717      |  0.54491059333263
0.70855958463834      |  0.70855959609279
0.03943927452247      |  0.03943927946198
0.18225709316049      |  0.18225709092897
0.14766834353602      |  0.14766834758223
0.00010493251424      |  0.00010493251317
0.00030783183885      |  0.00030783183432
0.00098583168085      |  0.00098583168214
0.20707624718093      |  0.20707625146740
0.00000003860552      |  0.00000003815330
0.00000419364072      |  0.00000419348257
0.00000627750835      |  0.00000627704688
0.00000785937722      |  0.00000785853673
0.00001189234779      |  0.00001189165231


