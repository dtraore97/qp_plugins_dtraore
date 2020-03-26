
subroutine GGA_sr_type_functionals_mu(mu,r,rho_a,rho_b,grad_rho_a_2,grad_rho_b_2,grad_rho_a_b, &
                                ex,vx_rho_a,vx_rho_b,vx_grad_rho_a_2,vx_grad_rho_b_2,vx_grad_rho_a_b, &
                                ec,vc_rho_a,vc_rho_b,vc_grad_rho_a_2,vc_grad_rho_b_2,vc_grad_rho_a_b ) 
 
 implicit none
 BEGIN_DOC
 ! routine that helps in building the x/c potentials on the AO basis for a GGA functional with a short-range interaction
 ! From Emmanuel's plugins: dft_utils_one_e/utils.irp.f
 !
 !-----------------------------------------------------------------------------------------------------------------------------------
 ! Parameter             : Definition                                                                   ; Source
 ! ----------------------------------------------------------------------------------------------------------------------------------
 ! a                     :                                                                              ; from (2), eq. ??
 ! b                     :                                                                              ; from (2), eq. ??
 ! beta_rho              :                                                                              ; from (2), eq. ??
 ! delta_rho             :                                                                              ; from (2), eq. ??
 ! dbeta_drho            : d(beta)/dn                                                                   ;
 ! ddeltamu_drho         : d(delta*mu)/dn                                                               ;
 ! dgamma_drho           : dgamma/dn                                                                    ;
 ! dgammamu2_drho        : d(gamma*mu**2)/dn                                                            ;
 ! dg0_drho              : dg0/dn                                                                       ; ask Emmanuel
 ! decPBE_drho           : decPBE/dn                                                                    ;
 ! dexPBE_drho           : dexPBE/dn                                                                    ;
 ! ddexPBE_ddrho         : d2exPBE/dn2                                                                  ;
 ! ec                    : Usual density of correlation energy                                          ; from (1)
 ! ex                    : Usual density of exchange energy                                             ; from (1)
 ! gamma_rho             :                                                                              ; from (2), eq. ??
 ! g0                    : On-top pair-distribution function of the spin-unpolarized UEG                ;
 ! g0_UEG_mu_inf         :                                                                              ; rsdft_ecmd/ueg_on_top.irp.f
 ! grad_rho              : gradient of density                                                          ; 
 ! n2_UEG                : On-top pair density of the uniform electron gas                              ;
 ! n2xc_UEG              : On-top exchange-correlation pair density of the uniform electron gas         ;
 ! rho                   : rho_a + rho_b (both densities of spins alpha and beta)                       ;
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
 !----------------------------------------------------------------------------------------------------------------------------------- 
 END_DOC
 double precision, intent(in)  :: r(3),rho_a(N_states),rho_b(N_states),grad_rho_a_2(N_states),grad_rho_b_2(N_states),grad_rho_a_b(N_states),mu
 double precision, intent(out) :: ex(N_states),vx_rho_a(N_states),vx_rho_b(N_states),vx_rho_pbe_a(N_states),vx_rho_pbe_b(N_states),vx_grad_rho_pbe_a_2(N_states),vx_grad_rho_pbe_b_2(N_states),vx_grad_rho_pbe_a_b(N_states),vx_rho_sr_pbe_md_a(N_states),vx_rho_sr_pbe_md_b(N_states),vx_grad_rho_sr_pbe_md_a_2(N_states),vx_grad_rho_sr_pbe_md_b_2(N_states),vx_grad_rho_sr_pbe_md_a_b(N_states)
 double precision, intent(out) :: ec(N_states),vc_rho_a(N_states),vc_rho_b(N_states),vc_rho_pbe_a(N_states),vc_rho_pbe_b(N_states),vc_grad_rho_pbe_a_2(N_states),vc_grad_rho_pbe_b_2(N_states),vc_grad_rho_pbe_a_b(N_states), vc_rho_sr_pbe_md_a(N_states),vc_rho_sr_pbe_md_b(N_states),vc_grad_rho_sr_pbe_md_a_2(N_states),vc_grad_rho_sr_pbe_md_b_2(N_states),vc_grad_rho_sr_pbe_md_a_b(N_states)
 integer          :: istate
 double precision :: r2(3),dr2(3), local_potential,r12,dx2, a, b, n2_UEG, n2xc_UEG, g0, g0_UEG_mu_inf
 
 do istate = 1, N_states
 !------------------Usual PBE--------------------
  call ex_pbe_sr(mu,rho_a(istate),rho_b(istate),grad_rho_a_2(istate),grad_rho_b_2(istate),grad_rho_a_b(istate),ex(istate),vx_rho_a(istate),vx_rho_b(istate),vx_grad_rho_a_2(istate),vx_grad_rho_b_2(istate),vx_grad_rho_a_b(istate))

  double precision :: rhoc,rhoo,sigmacc,sigmaco,sigmaoo,vrhoc,vrhoo,vsigmacc,vsigmaco,vsigmaoo
  ! convertion from (alpha,beta) formalism to (closed, open) formalism
   call rho_ab_to_rho_oc(rho_a(istate),rho_b(istate),rhoo,rhoc)
   call grad_rho_ab_to_grad_rho_oc(grad_rho_a_2(istate),grad_rho_b_2(istate),grad_rho_a_b(istate),sigmaoo,sigmacc,sigmaco)

   call ec_pbe_sr(mu,rhoc,rhoo,sigmacc,sigmaco,sigmaoo,ec(istate),vrhoc,vrhoo,vsigmacc,vsigmaco,vsigmaoo)

 !-------------------Multideterminant sr-PBE-----------------
 ! Constants
 ! 

 ! EVERYTHING ABOUT RHO HERE
 rho = rho_a + rho_b                                                                       ! OK
 grad_rho = dsqrt(grad_rho_a_2) + dsqrt(grad_rho_b_2)                                      ! ????
 
 ex_PBE = ex(istate)                                                                       ! OK
 ec_PBE = ec(istate)                                                                       ! OK
 dexPBE_drho= vx_rho_a(istate) + vx_rho_b(istate)                                          ! Not sure
 
 pi = dacos(-1.d0)                                                                         ! OK
 a = pi/2.d0                                                                               ! OK
 b = dsqrt(pi) + 2.d0*dsqrt(pi)*(2.d0*dsqrt(2.d0)-1.d0)/3.d0                               ! OK   
 c = 2.d0*dsqrt(pi)*(1.d0 - dsqrt(2.d0))/3.d0                                              ! OK
 n2_UEG =
 n2xc_UEG =
 g0 = g0_UEG_mu_inf(rho_a(istate), rho_b(istate))
 dg0_drho = ! Demander à Emmanuel où se trouve la dérivée de g0 dans QP
 beta_rho  = ec_PBE / (c*n2_UEG)
 dbeta_drho = (1.d0/c)*( (decPBE_drho/n2_UEG) + ec_PBE(2.d0*rho*g0 + (n**2)*dg0_drho)   )
 dbeta_decPBE = (c/(n2_UEG))
 delta_rho =
 gamma_rho =
 
 dexPBE_drho =
 ddexPBE_ddrho =
 
 decPBE_drho =

 d_gammamu = ((mu**2)/a)(d_exPBE - (1.d0/(n2xc_UEG**2))*(2.d0*g0*rho + (rho**2)*d_g0 - 2.d0*rho )
 d_gamma = d_gammamu/(mu**2)
 d_deltamu = - mu*b((2*g0*gamma_rho**2)/(exPBE) + (((gamma_rho**2)*(rho**2))/exPBE)*d_g0 + ((2*gamma_rho*n2_UEG)/exPBE)*d_gamma - ((n2_UEG*gamma_rho**2)/(exPBE**2))*d_exPBE)
 decsrmdpbedgradrho = !Eq 15   !!! Find an other label

 ! Exchange potential
 vx_rho_sr_pbe_md(istate) = d_exPBE*(1.d0/(1.d0 + delta_rho*mu + gamma_rho*mu**2)) - exPBE*((d_deltamu + d_gammamu)/((1.d0 + delta_rho*mu + gamma_rho*mu**2)**2))
 vx_grad_rho_sr_pbe_md(istate) = 2.d0*grad_rho*dd_exPBE*(1.d0/(1.d0 + delta_rho*mu + gamma_mu**2))*(1.d0 -(mu/(a*n2xc_UEG))  - (b*n2_UEG*gamma_rho**2)/((ex_PBE**2)*(1.d0 + delta_rho*mu + gamma_rho*mu**2)))

 ! Correlation potential
 vc_rho_sr_pbe_md(istate) = (1.d0 /(1.d0+beta_rho*mu**3))*(d_ecPBE) - ((ec_PBE*mu**3)/((1.d0 +beta_rho*mu**3)**2))*(d_beta)
 vc_grad_rho_sr_pbe_md(istate) = ((1.d0)/(1.d0 + beta_rho*mu**3)) - ((ecPBE*mu**3)/((1.d0 + beta_rho*mu**3)**2))*d_beta_decPBE
 !--------------------------Convert--------------------------
   call v_rho_oc_to_v_rho_ab(vrhoo,vrhoc,vc_rho_a(istate),vc_rho_b(istate)) ! conversion closed/open -> alpha/beta
   call v_grad_rho_oc_to_v_grad_rho_ab(vsigmaoo,vsigmacc,vsigmaco,vc_grad_rho_a_2(istate),vc_grad_rho_b_2(istate),vc_grad_rho_a_b(istate))
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


