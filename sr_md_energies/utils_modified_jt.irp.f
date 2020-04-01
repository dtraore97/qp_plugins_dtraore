 ! routine that helps in building the x/c potentials on the AO basis for a GGA functional with a short-range interaction
 ! From Emmanuel's plugins: dft_utils_one_e/utils.irp.f
 !
 !-----------------------------------------------------------------------------------------------------------------------------------
 ! Parameter             : Definition                                                                   ; Source
 ! ----------------------------------------------------------------------------------------------------------------------------------
 ! a, b and c            :                                                                              ; from (2), eq. ??
 ! B1, C1, D1, E1, F1    :                                                                              ; from (2)
 ! beta_rho              :                                                                              ; from (2), eq. ??
 ! delta_rho             :                                                                              ; from (2), eq. ??
 ! dbeta_drho            : d(beta)/dn                                                                   ;
 ! ddeltamu_drho         : d(delta*mu)/dn                                                               ;
 ! dgamma_drho           : dgamma/dn                                                                    ;
 ! dgammamu2_drho        : d(gamma*mu**2)/dn                                                            ;
 ! dg0_drho              : dg0/dn                                                                       ; ask Emmanuel
 ! dg0_drho              : dg0/dn                                                                       ; from (5), eq. 12
 ! dg0_drs               : dg0/drs                                                                      ; from (5), eq. 14
 ! decPBE_drho           : decPBE/dn                                                                    ;
 ! dexPBE_drho           : dexPBE/dn                                                                    ;
 ! ddexPBE_ddrho         : d2exPBE/dn2                                                                  ;
 ! ec                    : Usual density of correlation energy                                          ; from (1), already done in QP
 ! ex                    : Usual density of exchange energy                                             ; from (1), already done in QP
 ! gamma_rho             :                                                                              ; from (2), eq. ??
 ! g0                    : On-top pair-distribution function of the spin-unpolarized UEG                ;
 ! g0_UEG_mu_inf         :                                                                              ; rsdft_ecmd/ueg_on_top.irp.f
 ! grad_rho              : gradient of density                                                          ; 
 ! n2_UEG                : On-top pair density of the uniform electron gas                              ; from (2), eq. 51
 ! n2xc_UEG              : On-top exchange-correlation pair density of the uniform electron gas         ; from (2), below eq. 55
 ! rho                   : rho_a + rho_b (both densities of spins alpha and beta)                       ;
 ! rs                    : Seitz radius                                                                 ; rsdft_ecmd/ueg_on_top.irp.f
 ! vc_rho_sr_pbe_md      : derivative of ec_md^sr^mu^PBE with respect to rho                            ; from (3) and (4)
 ! vc_grad_rho_sr_pbe_md : derivative of ec_md^sr^mu^PBE with respect to grad_rho                       ; from (3) and (4)
 ! vx_rho_sr_pbe_md      : derivative of ex_md^sr^mu^PBE with respect to rho                            ; from (3) and (4)
 ! vx_grad_rho_sr_pbe_md : derivative of ex_md^sr^mu^PBE with respect to grad_rho                       ; from (3) and (4)
 !-----------------------------------------------------------------------------------------------------------------------------------
 ! SOURCES
 !-----------------------------------------------------------------------------------------------------------------------------------
 ! (1) : Generalized Gradient Approximation Made Simple - J. P. Perdew, K. Burke, M. Ernzerhof - PRL(77),18 (28/10/96)
 ! (2) : Short-range exchange and correlation density functionnals - J. Toulouse (Odense-Paris collaboration)
 ! (3) : A. Ferté's Thesis
 ! (4) : Developpement of eq. to be done
 ! (5) : Supplementary Materials for 'A density-based basis-set incomleteness correction for GW Methods' 
 !       - P.-F. Loos, B. Pradines, A. Scemama, E. Giner, J. Toulouse - ???????????
 !----------------------------------------------------------------------------------------------------------------------------------- 

subroutine excmdsrPBE(mu,rho_a,rho_b,grad_rho_a_2,grad_rho_b_2,grad_rho_a_b, &
       ex_srmuPBE,dexdrho_a,dexdrho_b,dexdgrad_rho_a_2,dexdgrad_rho_b_2,dexdgrad_rho_a_b, &
       ec_srmuPBE,decdrho_a,decdrho_b,decdgrad_rho_a_2,decdgrad_rho_b_2,decdgrad_rho_a_b)
 
 implicit none
 BEGIN_DOC
 ! Give exchange and correlation energies and chemical potentials
 ! Use qp_plugins_dtraore/sr_md_energies/utils_modified_jt.irp.f's plugins : exmdsrPBE and ecmdsrPBE
 END_DOC
 double precision, intent(in)  :: mu
 double precision, intent(in)  :: rho_a,rho_b,grad_rho_a_2,grad_rho_b_2,grad_rho_a_b
 double precision, intent(out) :: ex_srmuPBE,dexdrho_a,dexdrho_b,dexdgrad_rho_a_2,dexdgrad_rho_b_2,dexdgrad_rho_a_b
 double precision, intent(out) :: ec_srmuPBE,decdrho_a,decdrho_b,decdgrad_rho_a_2,decdgrad_rho_b_2,decdgrad_rho_a_b

 call exmdsrPBE(mu,rho_a,rho_b,grad_rho_a_2,grad_rho_b_2,grad_rho_a_b,ex_srmuPBE,dexdrho_a,dexdrho_b,dexdgrad_rho_a_2,dexdgrad_rho_b_2,dexdgrad_rho_a_b)
 call ecmdsrPBE(mu,rho_a,rho_b,grad_rho_a_2,grad_rho_b_2,grad_rho_a_b,ec_srmuPBE,decdrho_a,decdrho_b,decdgrad_rho_a_2,decdgrad_rho_b_2,decdgrad_rho_a_b)

 end subroutine excmdsrPBE

!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine g0_dg0(rho, rho_a, rho_b, g0, dg0drho)
  
  implicit none
  BEGIN_DOC
  ! Give the on-top pair distribution function g0 and its derivative according to rho dg0drho
  END_DOC

  double precision, intent (in) :: rho, rho_a, rho_b
  double precision, intent (out) :: g0, dg0drho
  double precision :: pi
  double precision :: g0_UEG_mu_inf, dg0drs
  double precision :: C1, F1, D1, E1, B1, rs

  pi = dacos(-1.d0)
  C1 = 0.0819306d0
  F1 = 0.752411d0
  D1 = -0.0127713d0
  E1 = 0.00185898d0
  B1 = 0.7317d0 - F1
  rs = (3.d0 / (4.d0*pi*rho))**(1.d0/3.d0) 
   
  g0 = g0_UEG_mu_inf(rho_a, rho_b)
  dg0drs = 0.5d0*((-B1 + 2.d0*C1*rs + 3.d0*D1*rs**2 + 4.d0*E1*rs**3)-F1*(1.d0 - B1*rs + C1*rs**2 + D1*rs**3 + E1*rs**4))*exp(-F1*rs)
  dg0drho = -((6.d0*dsqrt(pi)*rho**2)**(-2.d0/3.d0))*dg0drs
 
  end subroutine g0_dg0 

!-------------------------------------------------------------------------------------------------------------------------------------------
  subroutine exmdsrPBE(mu,rho_a,rho_b,grad_rho_a_2,grad_rho_b_2,grad_rho_a_b,ex_srmuPBE,dexdrho_a,dexdrho_b,dexdgrad_rho_a_2,dexdgrad_rho_b_2,dexdgrad_rho_a_b)
  
  implicit none
  BEGIN_DOC
  ! Calculation of exchange energy and chemical potential in PBE approximation using multideterminantal wave function (short-range part)
  END_DOC
  double precision, intent(in)  :: mu
  double precision, intent(in)  :: rho_a,rho_b,grad_rho_a_2,grad_rho_b_2,grad_rho_a_b
  double precision, intent(out) :: ex_srmuPBE,dexdrho_a,dexdrho_b,dexdgrad_rho_a_2,dexdgrad_rho_b_2,dexdgrad_rho_a_b
  double precision              :: exPBE,dexPBEdrho_a,dexPBEdrho_b, dexPBEdrho, dexPBEdgrad_rho_a_2,dexPBEdgrad_rho_b_2,dexPBEdgrad_rho_a_b
  double precision              :: dexdrho, gamma, dgammadrho, delta, ddeltadrho, denom, ddenomdrho
  double precision              :: pi, a, b
  double precision              :: rho, m  
  double precision              :: n2_UEG, dn2_UEGdrho, n2xc_UEG, dn2xc_UEGdrho, g0, dg0drho
! 
! n = na + nb  | m = na - nb
! na = (n+m)/2 | nb = (n-m)/2
! gn2 = ga2 + gb2 + 2 gab | ga2 = (gn2 + gm2 + 2gnm)/4
! gm2 = ga2 + gb2 - 2 gab | gb2 = (gn2 + gm2 - 2gnm)/4
! gnm = ga2 - gb2         | gab = (gn2 - gm2)/4

! dedn = dedna * dnadn + dednb * dnbdn = (1/2) * (dedna + dednb)
! dedgn2 = dedga2 * dga2/dgn2 + dedgb2 * dgb2dgn2 + dedgab * dgabdgn2 = 

  pi = dacos(-1.d0)
  rho = rho_a + rho_b
  m = rho_a - rho_b

! exchange PBE standard
  call ex_pbe_sr(1.d-12,rho_a,rho_b,grad_rho_a_2,grad_rho_b_2,grad_rho_a_b,exPBE,dexPBEdrho_a,dexPBEdrho_b,dexPBEdgrad_rho_a_2,dexPBEdgrad_rho_b_2,dexPBEdgrad_rho_a_b)
  call g0_dg0(rho, rho_a, rho_b, g0, dg0drho)

! calculation of energy
  a = pi / 2.d0
  b = 2*dsqrt(pi)*(2*dsqrt(2.d0) - 1.d0)/3.d0   
  n2_UEG = (rho**2)*g0
  n2xc_UEG = n2_UEG - rho**2
  gamma = exPBE / (a*n2xc_UEG)
  delta = -b*n2_UEG*gamma**2 / exPBE
  denom = 1.d0 + delta*mu + gamma*(mu**2)

  ex_srmuPBE=exPBE/denom

! calculation of derivatives
  dexPBEdrho = 0.5d0 *(dexPBEdrho_a + dexPBEdrho_b)
  dn2_UEGdrho = 2.d0*rho*g0 + (rho**2)*dg0drho
  dn2xc_UEGdrho = dn2_UEGdrho - 2.d0*rho

  dgammadrho = (1.d0/(a*n2xc_UEG))*dexPBEdrho  + (exPBE/(a*n2xc_UEG**2))*dn2xc_UEGdrho
  ddeltadrho = -((b*gamma**2)/exPBE)*dn2_UEGdrho -(b*n2_UEG/exPBE)*2.d0*gamma*dgammadrho + (b*n2_UEG*gamma**2)/(exPBE**2)*dexPBEdrho

  ddenomdrho = ddeltadrho * mu + dgammadrho * (mu**2)
  dexdrho = dexPBEdrho/denom - exPBE*ddenomdrho/denom**2
 
  end subroutine exmdsrPBE
!---------------------------------------------------------------------------------------------------------------------------------------------
  subroutine ecmdsrPBE(mu,rho_a,rho_b,grad_rho_a_2,grad_rho_b_2,grad_rho_a_b,ec_srmuPBE,decdrho_a,decdrho_b,decdgrad_rho_a_2,decdgrad_rho_b_2,decdgrad_rho_a_b)

  implicit none
  BEGIN_DOC
  ! Calculation of correlation energy and chemical potential in PBE approximation using multideterminantal wave function (short-range part)
  END_DOC
 
  double precision, intent(in)  :: mu
  double precision, intent(in)  :: rho_a,rho_b,grad_rho_a_2,grad_rho_b_2,grad_rho_a_b
  double precision, intent(out) :: ec_srmuPBE,decdrho_a,decdrho_b,decdgrad_rho_a_2,decdgrad_rho_b_2,decdgrad_rho_a_b
  double precision              :: ecPBE,decPBEdrho_a,decPBEdrho_b, decPBEdrho, decPBEdgrad_rho_a_2,decPBEdgrad_rho_b_2,decPBEdgrad_rho_a_b 
  double precision              :: rho_c, rho_o,grad_rho_c_2,grad_rho_o_2,grad_rho_c_o,decPBEdrho_c,decPBEdrho_o,decPBEdgrad_rho_c_2,decPBEdgrad_rho_o_2, decPBEdgrad_rho_c_o
  double precision              :: decdrho, beta, dbetadrho, denom, ddenomdrho
  double precision              :: pi, c
  double precision              :: rho, m  
  double precision              :: n2_UEG, dn2_UEGdrho, n2xc_UEG, dn2xc_UEGdrho, g0, dg0drho
 
  pi = dacos(-1.d0)
  rho = rho_a + rho_b
  m = rho_a - rho_b

! correlation PBE standard
  call ex_pbe_sr(1.d-12,rho_c,rho_o,grad_rho_c_2,grad_rho_o_2,grad_rho_c_o,ecPBE,decPBEdrho_c,decPBEdrho_o,decPBEdgrad_rho_c_2,decPBEdgrad_rho_o_2, decPBEdgrad_rho_c_o)
  call v_rho_oc_to_v_rho_ab(decPBEdrho_o, decPBEdrho_c, decPBEdrho_a, decPBEdrho_b)
  call v_grad_rho_oc_to_v_grad_rho_ab(decPBEdgrad_rho_o_2, decPBEdgrad_rho_c_2, decPBEdgrad_rho_c_o, decPBEdgrad_rho_a_2, decPBEdgrad_rho_b_2, decPBEdgrad_rho_a_b)
  call g0_dg0(rho, rho_a, rho_b, g0, dg0drho)

! calculation of energy
  c = 2*dsqrt(pi)*(1.d0 - dsqrt(2.d0))/3.d0
  n2_UEG = rho**2*g0
  beta = ecPBE/(c*n2_UEG)
  denom = 1.d0 + beta*mu**3
  ec_srmuPBE=ecPBE/denom

! calculation of derivatives 
  decPBEdrho = 0.5d0 *(decPBEdrho_a + decPBEdrho_b)
  dn2_UEGdrho = 2.d0*rho*g0 + (rho**2)*dg0drho
  dbetadrho = decPBEdrho/(c*n2_UEG) - (ecPBE/(c*n2_UEG**2))*dn2_UEGdrho
  ddenomdrho = dbetadrho*mu**3
  decdrho = decPBEdrho/denom - ecPBE*ddenomdrho/(denom**2)
 
  end subroutine ecmdsrPBE


