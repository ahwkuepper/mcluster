## McLuster - a tool to make a star cluster

 * Developer: Andreas H.W. Kuepper at [AIfA Bonn, Germany](https://astro.uni-bonn.de/en)
 * Collaborators: [Pavel Kroupa](https://astro.uni-bonn.de/~pavel/) &
   [Holger Baumgardt](https://people.smp.uq.edu.au/HolgerBaumgardt/)
 * Date: 25 March 2013

Detailed information on the code and the parameters can 
be found in the following paper:

> Kuepper A.H.W., Maschberger T., Kroupa P., Baumgardt H., 2011, MNRAS, 417, 2300
> [Mass segregation and fractal substructure in young massive clusters:
> (I) the McLuster code and method calibration](https://arxiv.org/abs/1107.2395)

This reference should be used for any publications that
utilize McLuster.

For comments, questions and bug reports contact:
ahwkuepper@gmail.com

### Compiling the code

For compiling the basic McLuster version use the MAKEFILE by typing
`make mcluster` or simply type:
```bash
cc -o mcluster main.c -O2 -fopenmp -lm
```

This should work on most computers on which a recent
C-compiler is installed.

For the SSE-version of McLuster you have to use the MAKEFILE.
Try compiling by typing `make mcluster_sse`.

You may have to change the C or FORTRAN compiler in the
MAKEFILE by editing the variable CC or FC, respectively.
Furthermore, you may have to give the correct directory of
your gfortran installation to the C compiler. This can be
done by changing the path in the CFLAGS variable `-L/usr/lib/`
to the correct path of your gfortran installation.

Note that some warnings like the following may appear:
ld warning: alignment lost in merging tentative definition
which can be ignored.

### Display available commands

Get a quick help by typing
```bash
mcluster -h
```
for the basic version, or
```bash
mcluster_sse -h
```
respectively for the SSE version.

### How to use it

McLuster can be run from the command line by passing arguments
to the code which specify the desired cluster, e.g.,
```bash
mcluster -N1000 -R2.0 -P0
```
would give a cluster with a Plummer profile and a half-mass
radius of 2 pc, consisting of 1000 stars. For parameters
which are not specified default values are taken, which can be
changed within the code. Go to the top of the 'main' routine
for a complete list of available parameters with descriptions.

By typing
```bash
mcluster -h
```
you get the allowed arguments which can be passed to McLuster
which will be described in detail in the following:

* *-N <number> (number of stars)*
   
   The number of stars may vary from 3 to ~10^6. Remember though
   that some procedures within the code require ~N^2 computational
   steps which may take incredibly long. Furthermore, processes
   like mass segregation and fractalization may temporarily
   need a lot of memory/or may temporarily generate more stars than
   the required number. The parameter NMAX within the main
   routine may have to be set to a higher value then.

* *-M <value> (mass of cluster; specify either N or M)*

   If you specify the total mass of your cluster then McLuster
   will generate stars until the total mass reaches M.

* *-P <0|1|2|3|-1> (density profile; 0= Plummer, 1= King (1966), 2= Subr et al. (2007) mass-segregated, 3= 2-dimensional EFF/Nuker template,-1= no density gradient)*

   The density profile is used to generate the positions &
   velocities of the stars. The final cluster is then scaled to
   exactly match the half-mass radius which you specified. In case
   you chose the EFF/Nuker profile McLuster, will integrate the density
   distribution and determine the expected half-mass radius of
   that specific distribution. The final cluster is then scaled
   to this expectation value.

   For the generation of the Subr et al. profile the routine
   PLUMIX is used, whose most recent version was incorporated
   in McLuster.

   Note that some choices for the profile necessitate further
   parameters to be specified.

* *-W <1-12> (W0 parameter for King model)*

   In case you choose the King profile you"ll have to specify its
   concentration which is done by fixing its W0 value.

* *-R <value> (half-mass radius [pc], ignored for P = 3)*

   The half-mass radius of the cluster has to be specified for all
   density profiles except for the EFF/Nuker profile. For the latter it
   gets determined automatically according to your choice of the
   scale radius and the cut-off radius. All radii are set in pc!
   By setting the half-mass radius to -1, McLuster will use the
   cluster mass-half-mass radius relation from Marks & Kroupa (2012)
   to estimate an initial half-mass radius according to your choice
   of mass.

* *-r <value> (scale radius of EFF/Nuker template [pc])*

   This specifies the scale radius, see Elson, Fall & Freeman (1987)
   and or Lauer et al. (1995) for Nuker template.

* *-c <value> (cut-off radius of EFF/Nuker template [pc])*

   Since the EFF and Nuker profiles extend infinitely they have to be
   cut off at some radius.

* *-g <value> (power-law slope(s) of EFF/Nuker template; use once for EFF template; use three times for Nuker template (outer slope, inner slope, transition)*
   
   With this option you can set the power-law slopes of the EFF and Nuker
   profiles. Use once for EFF to set the outer slope. Use three times for
   outer slope, inner slope and the transition parameter. Note that the EFF
   profile is a special case of the Nuker profile with the inner slope set
   to 0.0 and the transition parameter set to 2.0.

* *-S <0.0-1.0> (degree of mass segregation; 0.0= no segregation)*

   Mass segregation can be added to any choice of profile. In case of
   the Subr et al. profile the S value should lie between 0.0 and 0.5
   in order to produce reasonable clusters (see Subr, Kroupa & Baumgardt
   2007 for further details). All other density profiles allow values
   from 0.0 to 1.0, where 1.0 corresponds to total segregation, i.e.
   the most massive star sits most deeply within the cluster potential,
   followed by the second most massive star, etc. The segregation of
   stars is done following Baumgardt, de Marchi & Kroupa (2008).

* *-D <1.6-3.0> (fractal dimension; 3.0= no fractalization)*

   Any density profile can be fractalized. The fractal dimension
   specifies the degree of fractalization. Values of 2.6 are
   reasonable, and values lower than 1.6 should be omitted. The
   fractalization is realized similar as described in Whitworth &
   Goodwin (2004). If you choose P = -1 you can get a fractal
   distribution of stars without explicit density gradient.

* *-T <value> (tcrit in N-body units)*

   This specifies the time after which the N-body calculations
   should be stopped. For details see Sverre Aarseth's NbodyX
   handbook.

* *-Q <value> (virial ratio)*

   The virial ratio can be any positive value. 0.5 means virial
   equlibrium, values >0.5 mean initial expansion, and values
   <0.5 result in collapsing clusters.

* *-C <0|1|3|4> (code; 0= Nbody6, 1= Nbody4, 3= table of stars, 4= Nbody7)*

   If you want to use the output of McLuster as input for your
   N-body computations you should specify the corresponding code.
   You can also generate an ascii table of stars instead.

* *-A <value> (dtadj in N-body units)*

   Sets the adjustment time-step for Sverre's NbodyX codes.

* *-O <value> (deltat in N-body units)*

   Sets the output interval of NbodyX.

* *-G <0|1> (GPU usage; 0= no GPU, 1= use GPU)*

   Specify whether you want to use Nbody6 with a GPU.

* *-o <name> (output name of cluster model)*

   Change the name of the output file(s).

* *-f <0|1|2|3|4> (IMF; 0= no IMF, 1= Kroupa (2001), 2= user defined, 3= equal to "=1" but with optimal sampling, 4= L3 IMF from Maschberger 2012)*

   Set f to 0 for a cluster consisting of single-mass stars. For f
   set to 1 the canonical two-part power-law Kroupa IMF is used from
   0.08 Msun to 100 Msun. Those default values can be changed within
   the main routine of the code or with the -m option. The first -m
   will specify the lower mass limit, the second one will set the
   upper mass limit.
   
   If f is set to 2 you can specify the limiting masses and slopes
   between those limits for up to 10 power-law slopes and 11 mass limits.
   If you need more, change the parameters MAX_AN and MAX_MN within
   main.h accordingly.
   
   If f is set to 4 you can also use the -m option to set upper and
   lower mass limits, and you can also use the -a option to specify
   the power law slopes and the mu parameter of the L3 IMF. First -a
   will set alpha, the second will set beta and the third is for mu.

* *-a <value> (IMF slope for user defined IMF, may be used  multiple times, from low mass to high mass; or for parameters alpha, beta and mu of L3 IMF)*
* *-m <value> (IMF mass limits for the IMF; for user defined IMF it may be used multiple times; always set limits from low mass to high mass [Msun])*

   Specify the mass limits and the power-law slopes of the IMF as follows

   ```bash
   mcluster -f2 -m 0.08 -a -1.3 -m 0.5 -a -2.3 -m 150.0
   ```
   This would yield a Kroupa IMF ranging from 0.08 Msun to 150 Msun.
   ```bash
   mcluster -f4 -m 0.01 -m 100.0 -a 2.3 -a 1.4 -a 0.2
   ```
   This would give a L3 IMF from 0.01 to 100 Msun with alpha = 2.3,
   beta = 1.4 and mu = 0.2.

* *-B <number> (number of binary systems)*
* *-b <value> (binary fraction, specify either B or b)*

   For the binary content you can specify the fraction of stars
   in binaries, *0.0 <= b <= 1.0*, or the number of binary systems,
   *0 <= B <= N/2*. Note that the total mass of the cluster may exceed
   the mass you specified initially when using high binary fractions
   and eigenevolution. This is due to the form of the eigenevolution
   algorithm. Also, do not use eigenevolution when using BSE!

* *-p <0|1|2> (binary pairing, 0= random, 1= ordered for M>5Msun, 2= random but separate pairing for M>5Msun)*

   If p is 0 then the components of each binary are randomly drawn
   from the stars of the cluster. If p is 1 then first all stars above
   5 Msun are put in binaries, where the most massive star is put
   together with the second most massive star, the third with the
   forth, and so on until the desired binary fraction/number of
   binaries is reached, or until no more stars above 5 Msun are
   available. In this case the following stars are again drawn
   randomly. The limit between random and ordered sampling can be
   set within the main routine by changin the variable 'msort'. In
   addition, binaries with primary masses above msort can be chosen
   to follow a period distribution as given in Sana & Evans (2011) by
   setting the parameter 'OBperiods' in the main routine to 1 (recommended).
   If p is set to 2 then stars above 5 Msun are paired randomly among each
   other and the rest is also paired randomly.

* *-s <number> (seed for randomization; 0= randomize by timer)*

   If s is set to 0 then each call of McLuster will generate an
   independent realization of the specified cluster parameters (random
   number generator is seeded with the local time). Any other integer
   number will result in exactly the same realization (except for rounding
   differences).

* *-t <1|2|3> (tidal field; 1= near-field, 2= point-mass, 3= Milky-Way potential)*

   The tidal field can be specified according to the NbodyX manual.
   By setting the cluster's galactocentric radius and its orbital
   velocity, McLuster calculates the theoretical tidal radius and
   cuts off the Plummer profile at this radius (if P is set to 0).

* *-e <value> (epoch for stellar evolution [Myr])*

   This parameter can be used to set the cluster age but only with
   mcluster_sse. McLuster will use the SSE routines by Hurley et al.
   to evolve the cluster stars from their initial masses to the
   masses at the specified age. It also calculates the according
   luminosites, stellar radii, etc. If you are interested in a
   cluster of age zero but also need those stellar parameters, then
   use mcluster_sse and set e to a small value such as 0.1. Within
   the main.c file you can also activate BSE for binary evolution of
  primordial binaries. Therefore, set BSE = 1.

* *-Z <value> (metallicity [0.0001-0.03, 0.02 = solar])*

   Use Z to set the cluster metallicity. Alternatively, you can set
   the parameter FeH within the main routine to your desired [Fe/H]
   value such that McLuster will calculate the according Z value.

* *-X <value> (galactocentric radius vector, use 3x, [pc])*

   Use X three times to specify the position of the cluster within
   the galaxy, e.g.,
   ```bash
   mcluster -X 8500 -X 0 -X 0
   ```
   to set the cluster to a position of x = 8.5 kpc, y = z = 0 kpc.

* *-V <value> (cluster velocity vector, use 3x, [km/s])*

   With V you can specify the cluster's orbital velocity, e.g.,
   ```bash
   mcluster -V 0 -V 220 -V 0
   ```
   for a cluster orbital velocity of 220 km/s within the x-y-
   plane.

* *-u <0|1> (output units; 0= Nbody, 1= astrophysical)*

   When you generate a table of stars you may be interested in
   astrophysical units rather than Nbody units.

* *-h (display this help)*
* *-? (display this help)*

   Get help with those options.

Further parameters can be specified within the code. Note
that you have to compile the code again after any change.
Therefore, first type
```bash
make clean
```
ignore any error message, and compile the code again.
