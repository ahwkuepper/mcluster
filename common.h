*
      implicit none

*
* mcluster parameters
*
      integer, dimension(1:10) :: npop, initmodel, imfg, pairing,
     &         adis, eigen
      integer potential_energy, tf, mclusteron, seedmc,numpop,
     &         outputf

      REAL(KIND=8) qvir, rbar, zini,rh_mcl
      REAL(KIND=8), dimension(1:10) :: fracb, w0, conc_pop, Seg,
     &       equalmass, mlow, mup, alpha_L3,
     &       beta_L3, mu_L3, amin, amax, epoch_pop,
     &       zini_pop, fractal

      character(len=1000) :: npopchar,fracbchar, initmodelchar, w0char,
     &       conc_popchar, Segchar, fractalchar,
     &       imfgchar, equalmasschar, mlowchar, mupchar,
     &       alpha_L3char, beta_L3char, mu_L3char, pairingchar,
     &       adischar, eigenchar, aminchar, amaxchar,
     &       zinichar, epochchar
      character(len=100000) :: mlimimfchar, alphaimfchar

      common /mclusterarri/  npop, initmodel, imfg, pairing, adis,
     &        eigen

      common /mclusterarrd/ fracb, w0, conc_pop, Seg, fractal,
     &       equalmass, mlow, mup, alpha_L3,
     &       beta_L3, mu_L3, amin, amax, epoch_pop,
     &       zini_pop

      common /mclusteri/ potential_energy, tf, mclusteron, seedmc,
     &       outputf

      common /mclusterd/ qvir, rbar, rh_mcl
      common /mclusterchar/ alphaimfchar, mlimimfchar
