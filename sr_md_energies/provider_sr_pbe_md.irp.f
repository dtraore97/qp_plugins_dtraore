BEGIN_PROVIDER[double precision, energy_c_sr_pbe_md_copy2, (N_states) ]
&BEGIN_PROVIDER[double precision, energy_x_sr_pbe_md_copy2, (N_states) ]
 implicit none
 BEGIN_DOC
! Last modification : Mon., March 23th
! exchange/correlation energy with the short range pbe functional, multideterminantal form (VERSION 2)
! 11/03/20 Modified version of qp_plugins_eginer/stable/rsdft_cipsi/functionals/sr_pbe.irp
! ------------------------------------------PARAMETERS-------------------------------------------------
! a, b and c                   : Constants                                               eq. 45 and 47 
! beta_n_m_grad_n              :                                                         eq. 50
! contrib_grad                 : 
! delta_n_m_grad_n             :                                                         eq. 57
! energy_c_sr_pbe_md_copy2     : Sum(ec_prime*weight)                                    eq. 48
! energy_x_sr_pbe_md_copy2     : Sum(ex_prime*weight)                                    eq. 53
! ex & ec                      : ec_PBE & ex_PBE                                         eq. 49 and 54
! ex_prime & ec_prime          : exchange and correlation energies perelementaty volume  eq. 54 and 49
! gamma_n_m_grad_n             :                                                         eq. 55
! GGA_sr_type_functionals      : Emmanuel's                                              dft_utils_one_e/utils.irp.f
! GGA_sr_type_functionals_mu   : Modified version of GGA_sr_type_functionals where mu is intent(in) so we can have Ecsr value for a given mu
! g0_UEG_mu_inf  (g0)              : from rsdft_ecmd/ueg_on_top.irp.f
! grad_rho_a & grad_rho_b      : gradient of the densities of spin alpha and beta
! grad_rho_a_2 & grad_rho_b_2  : square root of grad_rho_a & grad_rho_b
! grad_rho_a_b                 : grad_rho_a*grad_rho_b
! m_spin                       : rho_a - rho_b
! mu                           : interaction parameter (constant for the moment)
! n2_UEG                       : On top pair density of the uniform electron gas         eq. 51
! n2_xc_UEG                    : on top exchange/correlation pair density of the UEG     eq. 55
! rho                          : rho_a + rho_b
! rho_a & rho_b                : densities of spin alpha and beta
! vc_rho_a & vc_rho_b          :
! vx_rho_a & vx_rho_b          :
! weight                       : dr
! zeta_m_n                     : m/n
! ----------------------------------------------------------------------------------------------------
 END_DOC
 integer :: istate,i,j,m, i_unit_output1, i_unit_output2, i_unit_output3, i_unit_output4, getUnitAndOpen
 double precision :: r(3), pi, a, b, c, g0_UEG_mu_inf, thr, r_norm, grad_rho_2, rho, m_spin, zeta_m_n, g0, n2_UEG, n2_xc_UEG, ec_prime, ex_prime, beta_n_m_delta_n, delta_n_m_delta_n, gamma_n_m_delta_n
 double precision :: mu,weight
 double precision, allocatable :: ex(:), ec(:)
 double precision, allocatable :: rho_a(:),rho_b(:),grad_rho_a(:,:),grad_rho_b(:,:),grad_rho_a_2(:),grad_rho_b_2(:),grad_rho_a_b(:)
 double precision, allocatable :: contrib_grad_xa(:,:),contrib_grad_xb(:,:),contrib_grad_ca(:,:),contrib_grad_cb(:,:)
 double precision, allocatable :: vc_rho_a(:), vc_rho_b(:), vx_rho_a(:), vx_rho_b(:)
 double precision, allocatable :: vx_grad_rho_a_2(:), vx_grad_rho_b_2(:), vx_grad_rho_a_b(:), vc_grad_rho_a_2(:), vc_grad_rho_b_2(:), vc_grad_rho_a_b(:)

 character*(128) :: output1, output2, output3, output4, type1 

 allocate(vc_rho_a(N_states), vc_rho_b(N_states), vx_rho_a(N_states), vx_rho_b(N_states))
 allocate(vx_grad_rho_a_2(N_states), vx_grad_rho_b_2(N_states), vx_grad_rho_a_b(N_states), vc_grad_rho_a_2(N_states), vc_grad_rho_b_2(N_states), vc_grad_rho_a_b(N_states))
 allocate(rho_a(N_states), rho_b(N_states),grad_rho_a(3,N_states),grad_rho_b(3,N_states))
 allocate(grad_rho_a_2(N_states),grad_rho_b_2(N_states),grad_rho_a_b(N_states), ex(N_states), ec(N_states))
 

!----------------Constantes------------------ 
 pi = dacos(-1.d0)
 a = pi/2.d0
 b = 2*dsqrt(pi)*(2*dsqrt(2.d0) - 1.d0)/3.d0  
 c = 2*dsqrt(pi)*(1.d0 - dsqrt(2.d0))/3.d0
 mu = mu_erf_dft
 thr = 1.d-12

 energy_x_sr_pbe_md_copy2 = 0.d0
 energy_c_sr_pbe_md_copy2 = 0.d0
!--------------------------------------------

do i = 1, n_points_final_grid
 r(1) = final_grid_points(1,i)
 r(2) = final_grid_points(2,i)
 r(3) = final_grid_points(3,i)
 weight = final_weight_at_r_vector(i)

 do istate = 1, N_states
  rho_a(istate) =  one_e_dm_and_grad_alpha_in_r(4,i,istate)
  rho_b(istate) =  one_e_dm_and_grad_beta_in_r(4,i,istate)
  grad_rho_a(1:3,istate) =  one_e_dm_and_grad_alpha_in_r(1:3,i,istate)
  grad_rho_b(1:3,istate) =  one_e_dm_and_grad_beta_in_r(1:3,i,istate)
  grad_rho_a_2 = 0.d0
  grad_rho_b_2 = 0.d0
  grad_rho_a_b = 0.d0
   
  do m = 1, 3
    grad_rho_a_2(istate) += grad_rho_a(m,istate) * grad_rho_a(m,istate)
    grad_rho_b_2(istate) += grad_rho_b(m,istate) * grad_rho_b(m,istate)
    grad_rho_a_b(istate) += grad_rho_a(m,istate) * grad_rho_b(m,istate)
  enddo
 enddo   
                           ! inputs
 call GGA_sr_type_functionals_mu(1.d-12,r,rho_a,rho_b,grad_rho_a_2,grad_rho_b_2,grad_rho_a_b,         &  ! outputs exchange      
                             ex,vx_rho_a,vx_rho_b,vx_grad_rho_a_2,vx_grad_rho_b_2,vx_grad_rho_a_b,   &  ! outputs correlation
                             ec,vc_rho_a,vc_rho_b,vc_grad_rho_a_2,vc_grad_rho_b_2,vc_grad_rho_a_b  ) 
 
 do istate = 1, N_states 
  if(dabs(ex(istate)).lt.thr)then
   ex(istate) = 1.d-12
  endif

  rho = rho_a(istate) + rho_b(istate)
  if(dabs(rho).lt.thr)then
   rho = 1.d-12
  endif 
  
  m_spin = rho_a(istate) - rho_b(istate)
  if(dabs(m_spin).lt.thr)then
   m_spin = 1.d-12
  endif
      
  zeta_m_n = m_spin/rho
  if(dabs(zeta_m_n).lt.thr)then
   zeta_m_n = 1.d-12
  endif

  g0 = g0_UEG_mu_inf(rho_a(istate), rho_b(istate))
  if(dabs(g0).lt.thr)then
   g0 = 1.d-12
  endif

  n2_UEG = (rho**2)*(1.0d0 - zeta_m_n**2)*g0
  if(dabs(n2_UEG).lt.thr)then
   n2_UEG = 1.d-12
  endif

  n2_xc_UEG = n2_UEG - rho**2
  if(dabs(n2_xc_UEG).lt.thr)then
   n2_xc_UEG = 1.d-12
  endif

  ! a, b and c defined before the loop
  
  beta_n_m_delta_n = ec(istate) / (c*n2_UEG)
  if(dabs(beta_n_m_delta_n).lt.thr)then
   beta_n_m_delta_n = 1.d-12
  endif

  gamma_n_m_delta_n = ex(istate) / (a*n2_xc_UEG)
  if(dabs(gamma_n_m_delta_n).lt.thr)then
   gamma_n_m_delta_n = 1.d-12
  endif

  delta_n_m_delta_n = - b*n2_UEG*gamma_n_m_delta_n**2 / ex(istate)
  if(dabs(delta_n_m_delta_n).lt.thr)then
   delta_n_m_delta_n = 1.d-12
  endif

  ec_prime = ec(istate) / (1.0d0 + beta_n_m_delta_n*mu**3)
  ex_prime = ex(istate) / (1.0d0 + delta_n_m_delta_n*mu + gamma_n_m_delta_n*mu**2)
 
  energy_c_sr_pbe_md_copy2(istate) += ec_prime * weight
  energy_x_sr_pbe_md_copy2(istate) += ex_prime * weight
  enddo
 enddo
END_PROVIDER
