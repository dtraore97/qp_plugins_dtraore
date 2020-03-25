
!subroutine GGA_sr_type_functionals_mu(mu,r,rho_a,rho_b,grad_rho_a_2,grad_rho_b_2,grad_rho_a_b, &
!                                ex,vx_rho_a,vx_rho_b,vx_grad_rho_a_2,vx_grad_rho_b_2,vx_grad_rho_a_b, &
!                                ec,vc_rho_a,vc_rho_b,vc_grad_rho_a_2,vc_grad_rho_b_2,vc_grad_rho_a_b ) 
 
! implicit none
! BEGIN_DOC
! ! routine that helps in building the x/c potentials on the AO basis for a GGA functional with a short-range interaction
 END_DOC
! double precision, intent(in)  :: r(3),rho_a(N_states),rho_b(N_states),grad_rho_a_2(N_states),grad_rho_b_2(N_states),grad_rho_a_b(N_states),mu
! double precision, intent(out) :: ex(N_states),vx_rho_a(N_states),vx_rho_b(N_states),vx_rho_pbe_a(N_states),vx_rho_pbe_b(N_states),vx_grad_rho_pbe_a_2(N_states),vx_grad_rho_pbe_b_2(N_states),vx_grad_rho_pbe_a_b(N_states),vx_rho_sr_pbe_md_a(N_states),vx_rho_sr_pbe_md_b(N_states),vx_grad_rho_sr_pbe_md_a_2(N_states),vx_grad_rho_sr_pbe_md_b_2(N_states),vx_grad_rho_sr_pbe_md_a_b(N_states)
! double precision, intent(out) :: ec(N_states),vc_rho_a(N_states),vc_rho_b(N_states),vc_rho_pbe_a(N_states),vc_rho_pbe_b(N_states),vc_grad_rho_pbe_a_2(N_states),vc_grad_rho_pbe_b_2(N_states),vc_grad_rho_pbe_a_b(N_states), vc_rho_sr_pbe_md_a(N_states),vc_rho_sr_pbe_md_b(N_states),vc_grad_rho_sr_pbe_md_a_2(N_states),vc_grad_rho_sr_pbe_md_b_2(N_states),vc_grad_rho_sr_pbe_md_a_b(N_states)
! integer          :: istate
! double precision :: r2(3),dr2(3), local_potential,r12,dx2
 
! do istate = 1, N_states
! !------------------Usual PBE--------------------
!  call ex_pbe_sr(mu,rho_a(istate),rho_b(istate),grad_rho_a_2(istate),grad_rho_b_2(istate),grad_rho_a_b(istate),ex(istate),vx_rho_a(istate),vx_rho_b(istate),vx_grad_rho_a_2(istate),vx_grad_rho_b_2(istate),vx_grad_rho_a_b(istate))

!  double precision :: rhoc,rhoo,sigmacc,sigmaco,sigmaoo,vrhoc,vrhoo,vsigmacc,vsigmaco,vsigmaoo
!  ! convertion from (alpha,beta) formalism to (closed, open) formalism
!   call rho_ab_to_rho_oc(rho_a(istate),rho_b(istate),rhoo,rhoc)
!   call grad_rho_ab_to_grad_rho_oc(grad_rho_a_2(istate),grad_rho_b_2(istate),grad_rho_a_b(istate),sigmaoo,sigmacc,sigmaco)

!   call ec_pbe_sr(mu,rhoc,rhoo,sigmacc,sigmaco,sigmaoo,ec(istate),vrhoc,vrhoo,vsigmacc,vsigmaco,vsigmaoo)

 !-------------------Multideterminant sr-PBE-----------------
! ! Constants
! pi = dacos(-1.d0)
! beta_rho_delta_rho_m = 3.d0 /(2.d0*dsqrt(pi)*(1.d0 - dsqrt(2.d0)))
! rho = rho_a(istate) + rho_b(istate)
! dbetadrho =  !Eq. 11          !!! Find an other label
! decsrmdpbedgradrho = !Eq 15   !!! Find an other label
!! g0 = g0_UEG_inf(rho_a(istate), rho_b(istate))
!! n2_UEG = (rho**2)*g0

! ! Exchange potential
! vx_rho_sr_pbe_md_a(istate) = 
! vx_rho_sr_pbe_md_b(istate)
! vx_grad_rho_sr_pbe_md_a_2(istate)
! vx_grad_rho_sr_pbe_md_b_2(istate)
! vx_grad_rho_sr_pbe_md_a_b(istate)

 ! Correlation potential
! vc_rho_sr_pbe_md_a(istate) = (1.d0 /(1.d0+beta_rho_delta_rho_m*mu**3))*('decpbe/dn') - ('ecpbe'*mu**3)/((1.d0 +beta_rho_delta_rho*mu**3)**2)*('dBeta/dn') !Derivative of ex^(sr,PBE)_(md) with respect to the density eq. 10 : Supplementary Materials for 'A Density-Based Basis-set Incompleteness Correction for GW Methods'
! vc_rho_sr_pbe_md_b(istate)
! vc_grad_rho_sr_pbe_md_a_2(istate)
! vc_grad_rho_sr_pbe_md_b_2(istate)
! vc_grad_rho_sr_pbe_md_a_b(istate)
 !--------------------------Convert--------------------------
!   call v_rho_oc_to_v_rho_ab(vrhoo,vrhoc,vc_rho_a(istate),vc_rho_b(istate)) ! conversion closed/open -> alpha/beta
!   call v_grad_rho_oc_to_v_grad_rho_ab(vsigmaoo,vsigmacc,vsigmaco,vc_grad_rho_a_2(istate),vc_grad_rho_b_2(istate),vc_grad_rho_a_b(istate))
! enddo
!end

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


