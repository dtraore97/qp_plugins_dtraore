
subroutine GGA_sr_type_functionals_mu(mu_usual,mu,r,rho_a,rho_b,grad_rho_a_2,grad_rho_b_2,grad_rho_a_b, &
       ex,vx_rho_a,vx_rho_b,vx_grad_rho_a_2,vx_grad_rho_b_2,vx_grad_rho_a_b, vx_rho_sr_pbe_md1, vx_rho_sr_pbe_md2, vx_rho_sr_pbe_md3, &
       ec,vc_rho_a,vc_rho_b,vc_grad_rho_a_2,vc_grad_rho_b_2,vc_grad_rho_a_b, vc_rho_sr_pbe_md, vc_grad_rho_sr_pbe_md)
 
 implicit none
 BEGIN_DOC
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
 END_DOC
 double precision, intent(in)  :: r(3),rho_a(N_states),rho_b(N_states),grad_rho_a_2(N_states),grad_rho_b_2(N_states),grad_rho_a_b(N_states),mu, mu_usual
 double precision, intent(out) :: ex(N_states),vx_rho_a(N_states),vx_rho_b(N_states),vx_grad_rho_a_2(N_states),vx_grad_rho_b_2(N_states),vx_grad_rho_a_b(N_states), vx_rho_sr_pbe_md1(N_states), vx_rho_sr_pbe_md2(N_states), vx_rho_sr_pbe_md3(N_states)
 double precision, intent(out) :: ec(N_states),vc_rho_a(N_states),vc_rho_b(N_states),vc_grad_rho_a_2(N_states),vc_grad_rho_b_2(N_states),vc_grad_rho_a_b(N_states), vc_rho_sr_pbe_md(N_states), vc_grad_rho_sr_pbe_md(N_states) 
 integer          :: istate
 double precision :: r2(3),dr2(3), local_potential,r12,dx2, pi, a, b, c, rho, grad_rho, grad_rho_a, grad_rho_b, ex_PBE, ec_PBE, dexPBE_drho, ddexPBE_ddrho, decPBE_drho, n2_UEG, n2xc_UEG, g0, dg0_drho, g0_UEG_mu_inf, beta_rho, dbeta_drho, dbeta_decPBE, gamma_rho, delta_rho, dgammamu2_drho, dgamma_drho, ddeltamu_drho, C1, F1, D1, E1, B1, rs, dg0_drs 
 
 do istate = 1, N_states
 !------------------Usual PBE--------------------
  call ex_pbe_sr(mu_usual,rho_a(istate),rho_b(istate),grad_rho_a_2(istate),grad_rho_b_2(istate),grad_rho_a_b(istate),ex(istate),vx_rho_a(istate),vx_rho_b(istate),vx_grad_rho_a_2(istate),vx_grad_rho_b_2(istate),vx_grad_rho_a_b(istate))

  double precision :: rhoc,rhoo,sigmacc,sigmaco,sigmaoo,vrhoc,vrhoo,vsigmacc,vsigmaco,vsigmaoo
  ! convertion from (alpha,beta) formalism to (closed, open) formalism
   call rho_ab_to_rho_oc(rho_a(istate),rho_b(istate),rhoo,rhoc)
   call grad_rho_ab_to_grad_rho_oc(grad_rho_a_2(istate),grad_rho_b_2(istate),grad_rho_a_b(istate),sigmaoo,sigmacc,sigmaco)

   call ec_pbe_sr(mu_usual,rhoc,rhoo,sigmacc,sigmaco,sigmaoo,ec(istate),vrhoc,vrhoo,vsigmacc,vsigmaco,vsigmaoo)

 !--------------------------Convert--------------------------
   call v_rho_oc_to_v_rho_ab(vrhoo,vrhoc,vc_rho_a(istate),vc_rho_b(istate)) ! conversion closed/open -> alpha/beta
   call v_grad_rho_oc_to_v_grad_rho_ab(vsigmaoo,vsigmacc,vsigmaco,vc_grad_rho_a_2(istate),vc_grad_rho_b_2(istate),vc_grad_rho_a_b(istate))

 !-----------------------------------------------------------
 ! CONSTANTS
 pi = dacos(-1.d0)                                                                          ! OK
 a  = pi/2.d0                                                                               ! OK
 b  = dsqrt(pi) + 2.d0*dsqrt(pi)*(2.d0*dsqrt(2.d0)-1.d0)/3.d0                               ! OK
 c  = 2.d0*dsqrt(pi)*(1.d0 - dsqrt(2.d0))/3.d0                                              ! OK
 C1 = 0.0819306d0                                                                           ! OK
 F1 = 0.752411d0                                                                            ! OK
 D1 = -0.0127713d0                                                                          ! OK
 E1 = 0.00185898d0                                                                          ! OK
 B1 = 0.7317d0 - F1                                                                         ! OK

 ! RHO
 rho        = rho_a(istate) + rho_b(istate)                             ! OK
 grad_rho_a = dsqrt(grad_rho_a_2(istate))                               ! Pas sûre
 grad_rho_b = dsqrt(grad_rho_b_2(istate))                               ! Pas sûre
 grad_rho   = grad_rho_a + grad_rho_b                                   ! Pas sûre

 ! ENERGIES DENSITIES / DERIVATIVES
 ex_PBE        = ex(istate)                                                                                ! OK
 ec_PBE        = ec(istate)                                                                                ! OK
 dexPBE_drho   = vx_rho_a(istate) + vx_rho_b(istate)                                                       ! OK
 ddexPBE_ddrho = 2.d0*vx_grad_rho_a_2(istate)*grad_rho_a + 2.d0*vx_grad_rho_b_2(istate)*grad_rho_b         ! Ferte's notes
 decPBE_drho   = vc_rho_a(istate) + vc_rho_b(istate)
 
 ! ON-TOP PAIR DENSITIES 
 g0       = g0_UEG_mu_inf(rho_a(istate), rho_b(istate))                                                    ! OK
 rs       = (3.d0 / (4.d0*pi*rho))**(1.d0/3.d0)                                                            ! OK
 dg0_drs  = 0.5d0*((-B1 + 2.d0*C1*rs + 3.d0*D1*rs**2 + 4.d0*E1*rs**3)-F1*(1.d0 - B1*rs + C1*rs**2 + D1*rs**3 + E1*rs**4))*exp(-F1*rs)    ! OK
 dg0_drho = -((6.d0*dsqrt(pi)*rho**2)**(-2.d0/3.d0))*dg0_drs                                               ! OK 
 n2_UEG   = (rho**2)*g0                                             ! OK (Attention cas spin = 0) sinon n^2(1-zeta^2)*g0
 n2xc_UEG = n2_UEG - rho**2                                         ! OK         

 ! ALPHA/BETA/GAMMA AND DERIVATIVES
 beta_rho       = ec_PBE / (c*n2_UEG)                                                                            ! OK
 dbeta_drho     = (1.d0/c)*( (decPBE_drho/n2_UEG) - (ec_PBE/(n2_UEG**2))*(2*rho*g0 + dg0_drho*rho**2) )           ! OK
 dbeta_decPBE   = (c/(n2_UEG))                                                                                   ! OK
 gamma_rho      = ex_PBE/(a*n2xc_UEG)                                                                            ! OK
 delta_rho      = -(b*n2_UEG*gamma_rho**2)/ex_PBE                                                                ! OK
 dgammamu2_drho = ((mu**2)/a)*(dexPBE_drho - (1.d0/(n2xc_UEG**2))*(2.d0*g0*rho + (rho**2)*dg0_drho - 2.d0*rho )) ! OK
 dgamma_drho    = dgammamu2_drho/(mu**2)                                                                         ! OK
 ddeltamu_drho  = -mu*b*( (2*g0*rho*gamma_rho**2)/(ex_PBE) +  ((gamma_rho**2)*(rho**2)/(ex_PBE))*dg0_drho + ((2.d0*gamma_rho*n2_UEG)/ex_PBE)*dgamma_drho - (n2_UEG*gamma_rho**2/(ex_PBE**2))*dexPBE_drho )                                                     ! Copy OK 

 ! EXCHANGE POTENTIAL
 vx_rho_sr_pbe_md1(istate)      = dexPBE_drho*(1.d0/(1.d0 + delta_rho*mu + gamma_rho*mu**2)) - ex_PBE*((ddeltamu_drho + dgammamu2_drho)/((1.d0 + delta_rho*mu + gamma_rho*mu**2)**2))                                ! NON
 vx_rho_sr_pbe_md2(istate)      = (mu*ex_PBE/((1.d0 + gamma_rho*mu**2 + delta_rho*mu)**2)) * ( -mu*(dexPBE_drho*(1.d0/(a*n2xc_UEG)) - (ex_PBE/a)*((2.d0*rho*(g0 - 1.d0))/(n2xc_UEG**2))) + b*((2.d0*rho*g0*gamma_rho**2)/(ex_PBE) + ((2.d0*gamma_rho*n2_UEG)/ex_PBE)*dgamma_drho) - ((n2_UEG*gamma_rho**2)/(ex_PBE**2))*dexPBE_drho ) + dexPBE_drho*(1.d0/(1.d0 + gamma_rho*mu**2 + delta_rho*mu)) ! EN COURS DE TEST
 
 double precision :: Part1, Part2, Part21, Part22
  Part1                         = -(mu/((1.d0 + gamma_rho*mu**2 + delta_rho*mu)**2))*( (dexPBE_drho/(a*n2xc_UEG)) + (ex_PBE/(a*n2xc_UEG**2))*(2.d0*rho*(g0-1.d0) + (rho**2)*dg0_drho) )
  Part21                        = -(b*n2_UEG/ex_PBE)*2.d0*gamma_rho*(dexPBE_drho/(a*n2xc_UEG) - (ex_PBE/(a*n2xc_UEG**2))*(2.d0*rho*(g0-1.d0)+(rho**2)*dg0_drho))
  Part22                        = ((b*n2_UEG*gamma_rho**2)/(ex_PBE**2))*dexPBE_drho - ((b*gamma_rho**2)/ex_PBE)*(2.d0*rho*g0 + (rho**2)*dg0_drho)
  Part2                         = (mu/((1.d0 + gamma_rho*mu**2 + delta_rho*mu)**2))*( Part21 + Part22  )
 vx_rho_sr_pbe_md3(istate)       = dexPBE_drho*(1.d0/(1.d0 + gamma_rho*mu**2 + delta_rho*mu)) + ex_PBE*(Part1 + Part2 )
! vx_grad_rho_sr_pbe_md(istate) = 2.d0*grad_rho*ddexPBE_ddrho*(1.d0/(1.d0 + delta_rho*mu + gamma_rho**2))*(1.d0 -(mu/(a*n2xc_UEG))  - (b*n2_UEG*gamma_rho**2)/((ex_PBE**2)*(1.d0 + delta_rho*mu + gamma_rho*mu**2)))  ! NON

 ! CORRELATION POTENTIAL
 vc_rho_sr_pbe_md(istate)      = (1.d0 /(1.d0+beta_rho*mu**3))*(decPBE_drho) - ((ec_PBE*mu**3)/((1.d0 +beta_rho*mu**3)**2))*(dbeta_drho)
 vc_grad_rho_sr_pbe_md(istate) = ((1.d0)/(1.d0 + beta_rho*mu**3)) - ((ec_PBE*mu**3)/((1.d0 + beta_rho*mu**3)**2))*dbeta_decPBE

 enddo
end

subroutine GGA_sr_type_functionals_mu_ok(mu,r,rho_a,rho_b,grad_rho_a_2,grad_rho_b_2,grad_rho_a_b, &
                                ex,vx_rho_a,vx_rho_b,vx_grad_rho_a_2,vx_grad_rho_b_2,vx_grad_rho_a_b, &
                                ec,vc_rho_a,vc_rho_b,vc_grad_rho_a_2,vc_grad_rho_b_2,vc_grad_rho_a_b )
 implicit none
 BEGIN_DOC
 ! routine that helps in building the x/c potentials on the AO basis for a GGA functional with a short-range interaction
 END_DOC
 double precision, intent(in)  :: r(3),rho_a(N_states),rho_b(N_states),grad_rho_a_2(N_states),grad_rho_b_2(N_states),grad_rho_a_b(N_states)
 double precision, intent(out) :: ex(N_states),vx_rho_a(N_states),vx_rho_b(N_states),vx_grad_rho_a_2(N_states),vx_grad_rho_b_2(N_states),vx_grad_rho_a_b(N_states)
 double precision, intent(out) :: ec(N_states),vc_rho_a(N_states),vc_rho_b(N_states),vc_grad_rho_a_2(N_states),vc_grad_rho_b_2(N_states),vc_grad_rho_a_b(N_states)
 integer          :: istate
 double precision :: r2(3),dr2(3), local_potential,r12,dx2,mu
 do istate = 1, N_states
  call ex_pbe_sr(mu,rho_a(istate),rho_b(istate),grad_rho_a_2(istate),grad_rho_b_2(istate),grad_rho_a_b(istate),ex(istate),vx_rho_a(istate),vx_rho_b(istate),vx_grad_rho_a_2(istate),vx_grad_rho_b_2(istate),vx_grad_rho_a_b(istate))

  double precision :: rhoc,rhoo,sigmacc,sigmaco,sigmaoo,vrhoc,vrhoo,vsigmacc,vsigmaco,vsigmaoo
  ! convertion from (alpha,beta) formalism to (closed, open) formalism
   call rho_ab_to_rho_oc(rho_a(istate),rho_b(istate),rhoo,rhoc)
   call grad_rho_ab_to_grad_rho_oc(grad_rho_a_2(istate),grad_rho_b_2(istate),grad_rho_a_b(istate),sigmaoo,sigmacc,sigmaco)

   call ec_pbe_sr(mu,rhoc,rhoo,sigmacc,sigmaco,sigmaoo,ec(istate),vrhoc,vrhoo,vsigmacc,vsigmaco,vsigmaoo)

   call v_rho_oc_to_v_rho_ab(vrhoo,vrhoc,vc_rho_a(istate),vc_rho_b(istate))
   call v_grad_rho_oc_to_v_grad_rho_ab(vsigmaoo,vsigmacc,vsigmaco,vc_grad_rho_a_2(istate),vc_grad_rho_b_2(istate),vc_grad_rho_a_b(istate))
 enddo
end


