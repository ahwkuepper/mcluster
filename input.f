      subroutine input
*
*
*       read input parameters
*       ----------------
*
      include 'common.h'
*
      integer i,ixx1,ixx2,ixx3,verboseTmp, tmpInt
*
      real*8 xxx4, tmp, bodyMin, bodyMax
      real*4 xxx32
      CHARACTER(len=32) :: arg
      CHARACTER(len=100) :: key
      CHARACTER(len=1000) :: conf
      character*1 st1
      character*2 st2

      ixx1 = 0
      ixx2 = 0
      ixx3 = 0
      xxx4 = 0.0d0

      conf = "mcluster.ini"


      ! check input parameters
      DO i = 1, iargc()
        CALL getarg(i, arg)
          print*,"INFO: verbose mode enabled"
        if (arg.eq."--help".OR.arg.eq."-h") then
      print*,"MOCCA = MOnte Carlo Cluster simulAtor"
      print*,""
      print*,"  --help    | -h   prints this help"
      print*,"  --conf    | -c   specify configuration file with the "
      print*,"                   initial parameters. If file is not "
      print*,"                   specified then file mocca.ini is read"
      print*,"                   by default. If mocca.ini file is not "
      print*,"                   present then all default values are "
      print*,"                   used"
          stop
        else if (arg.eq."--conf".OR.arg.eq."-c") then
          CALL getarg(i + 1, conf)
          print*,"INFO: reading configuration from ",conf
        endif
!        WRITE (*,*) arg
      END DO

      ! check if there are unrecognized parameters in the mocca.ini file
      call config_validate(conf)
     
      key = "Mcluster:n"
      call config_getstr(npopchar, "100000", key, conf);
      call char_to_arrayint(npopchar, npop);

      numpop = 0
      do i = 1, 10
        numpop = numpop + 1
        ixx2 = ixx2 + npop(i)
        if (npop(i).EQ.0) then
          numpop = numpop - 1
          exit
        endif
      end do

      print*,"Total number of objects ",ixx2

      key = "Mcluster:fracb"
      call config_getstr(fracbchar, "0.2", key, conf);
      call char_to_arraydouble(fracbchar, fracb);
      do i = 1, numpop
      if (fracb(i).lt.0.0d0.OR.fracb(i).gt.1.0d0) then
      print*,"ERROR Mcluster:fracb ",i,"-th is wrong ",
     &        "file. It should be in the range [0,1]"
        stop
      endif
      end do

      key = "Mcluster:initialModel"
      call config_getstr(initmodelchar, "1", key, conf);
      call char_to_arrayint(initmodelchar, initmodel);
      do i = 1, numpop
      if (initmodel(i).lt.0.OR.initmodel(i).gt.4) then
      print*,"ERROR Mcluster:initialModel ",i,"-th",
     &     " is wrong. It should be 0, 1, 2, 3 or 4"
         stop
      endif 
      end do

      key = "Mcluster:w0"
      call config_getstr(w0char, "1", key, conf);
      call char_to_arraydouble(w0char, w0);
      do i = 1, numpop
       if (w0(i).le.1.0d0.OR.w0(i).gt.12.0d0) then
      print*,"ERROR Mcluster:w0 ", i, "th",
     &   " is wrong. It should be between 1.0 and 12.0"
        stop
      endif
      end do
      
      key = "Mcluster:S"
      call config_getstr(Segchar, "0.0", key, conf);
      call char_to_arraydouble(Segchar, Seg);
      do i = 1, numpop
       if (Seg(i).lt.0.0d0.OR.Seg(i).gt.1.0d0) then
      print*,"ERROR Mcluster:S ",i,"th",
     &   " is wrong. It should be between 0 and 1"
         stop
      endif
      end do
       
      key = "Mcluster:fractal"
      call config_getstr(fractalchar, "3.0", key, conf);
      call char_to_arraydouble(fractalchar, fractal);
      do i = 1, numpop
       if (fractal(i).lt.1.6d0.OR.fractal(i).gt.3.0d0) then
      print*,"ERROR Mcluster:fractal ",i,"th",
     &   " is wrong. It should be between 1.6 and 3.0"
         stop
      endif
      end do

      key = "Mcluster:qvir"
      call config_getdouble(qvir, 0.5d0, key, conf);

      key = "Mcluster:mfunc"
      call config_getstr(imfgchar, "1.0", key, conf);
      call char_to_arrayint(imfgchar, imfg);

      key = "Mcluster:single_mass"
      call config_getstr(equalmasschar, "1.0", key, conf);
      call char_to_arraydouble(equalmasschar, equalmass);

      key = "Mcluster:mlow"
      call config_getstr(mlowchar, "0.08", key, conf);
      call char_to_arraydouble(mlowchar, mlow);

      key = "Mcluster:mup"
      call config_getstr(mupchar, "100.0", key, conf);
      call char_to_arraydouble(mupchar, mup);

      key = "Mcluster:alpha_imf"
      call config_getstr(alphaimfchar,
     &      "-1.35, -2.35, -2.7, 0.0, 0.0", key, conf);

      key = "Mcluster:mlim_imf"
      call config_getstr(mlimimfchar,
     &    "0.08, 0.5, 4.0, 100.0, 0.0, 0.0", key, conf);

      key = "Mcluster:alpha_L3"
      call config_getstr(alpha_L3char, "2.3", key, conf);
      call char_to_arraydouble(alpha_L3char, alpha_L3);

      key = "Mcluster:beta_L3"
      call config_getstr(beta_L3char, "1.4", key, conf);
      call char_to_arraydouble(beta_L3char, beta_L3);

      key = "Mcluster:mu_L3"
      call config_getstr(mu_L3char, "0.2", key, conf);
      call char_to_arraydouble(mu_L3char, mu_L3);

      key = "Mcluster:pairing"
      call config_getstr(pairingchar, "3", key, conf);
      call char_to_arrayint(pairingchar, pairing);
      
      key = "Mcluster:adis"
      call config_getstr(adischar, "3", key, conf);
      call char_to_arrayint(adischar, adis);
      do i = 1, numpop
       if (adis(i).eq.6.AND.fracb(i).gt.0.5d0) then
      print*,"ERROR Mcluster:fracb ",i,"th",
     &   " is wrong. It should not be greater",
     &   " then 0.5 if  adis = 6"
         stop
      endif
      end do

      key = "Mcluster:eigen"
      call config_getstr(eigenchar, "0", key, conf);
      call char_to_arrayint(eigenchar, eigen);

      key = "Mcluster:amin"
      call config_getstr(aminchar, "-1.0", key, conf);
      call char_to_arraydouble(aminchar, amin);

      key = "Mcluster:amax"
      call config_getstr(amaxchar, "10747.0", key, conf);
      call char_to_arraydouble(amaxchar, amax);

      key = "Mcluster:tf"
      call config_getint(tf, 0, key, conf);

      key = "Mcluster:rbar"
      call config_getdouble(rbar, 35.8d0, key, conf);

      key = "Mcluster:rh_mcl"
      call config_getdouble(rh_mcl, 1.0d0, key, conf);

      key = "Mcluster:conc_pop"
      call config_getstr(conc_popchar, "0.5", key, conf);
      call char_to_arraydouble(conc_popchar, conc_pop);

      key = "Mcluster:potential_energy"
      call config_getint(potential_energy, 1, key, conf);

      key = "Mcluster:epoch"
      call config_getstr(epochchar, "0.0", key, conf);
      call char_to_arraydouble(epochchar, epoch_pop);

      key = "Mcluster:zini"
      call config_getstr(zinichar, "0.001", key, conf);
      call char_to_arraydouble(zinichar, zini_pop);
      do i = 1, numpop
          if (zini_pop(i).le.0.0d0) then
            print*,"ERROR StarCluster:zini ",i,"-th",
     &        " is wrong.  It should be > 0.0"
              stop
       endif
      end do
      zini = zini_pop(1);

      key = "Mcluster:seedmc"
      call config_getint(seedmc, 0, key, conf);

      key = "Mcluster:outputf"
      call config_getint(outputf, 0, key, conf);

      key = "Mcluster:check_en"
      call config_getint(check_en, 1, key, conf);

*
      return
*
      end
*
