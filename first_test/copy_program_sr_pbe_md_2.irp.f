program energy_x_c_md_test_2

 implicit none
 BEGIN_DOC
! exchange/correlation energy with the short range pbe functional, multideterminantal form (VERSION 2)
! 11/03/20 Modified version of qp_plugins_eginer/stable/rsdft_cipsi/functionals/sr_pbe.irp.f
! This program aim to test the providers from sr_pbe_version_2.irp.f
! ------------------------------------------PARAMETERS-------------------------------------------------
! a, b and c                   : Constants                                               eq. 45 and 47 
! beta_n_m_grad_n              :                                                         eq. 50
! contrib_grad                 : 
! delta_n_m_grad_n             :                                                         eq. 57
! energy_c_sr_pbe_md_copy      : Sum(ec_prime*weight)                                    eq. 48
! energy_x_sr_pbe_md_copy      : Sum(ex_prime*weight)                                    eq. 53
! ex & ec                      : ec_PBE & ex_PBE                                         eq. 49 and 54
! ex_prime & ec_prime          : exchange and correlation energies perelementaty volume  eq. 54 and 49
! gamma_n_m_grad_n             :                                                         eq. 55
! g0_UEG_mu  (g0)              : from rsdft_ecmd/ueg_on_top.irp.f
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
 integer :: istate,i,j,m
 double precision :: r(3)
 double precision :: mu,weight
 double precision, allocatable :: ex(:), ec(:)
 double precision, allocatable :: rho_a(:),rho_b(:),grad_rho_a(:,:),grad_rho_b(:,:),grad_rho_a_2(:),grad_rho_b_2(:),grad_rho_a_b(:)
 double precision, allocatable :: contrib_grad_xa(:,:),contrib_grad_xb(:,:),contrib_grad_ca(:,:),contrib_grad_cb(:,:)
 double precision, allocatable :: vc_rho_a(:), vc_rho_b(:), vx_rho_a(:), vx_rho_b(:)
 double precision, allocatable :: vx_grad_rho_a_2(:), vx_grad_rho_b_2(:), vx_grad_rho_a_b(:), vc_grad_rho_a_2(:), vc_grad_rho_b_2(:), vc_grad_rho_a_b(:)
!-----------Added-----------
 double precision :: pi, a, b, c, g0_UEG_mu
 double precision, allocatable :: rho(:), m_spin(:), zeta_m_n(:), g0(:), n2_UEG(:), n2_xc_UEG(:), ec_prime(:), ex_prime(:), beta_n_m_delta_n(:), delta_n_m_delta_n(:), gamma_n_m_delta_n(:), energy_x_sr_pbe_md_copy(:), energy_c_sr_pbe_md_copy(:)
!---------------------------

 allocate(vc_rho_a(N_states), vc_rho_b(N_states), vx_rho_a(N_states), vx_rho_b(N_states))
 allocate(vx_grad_rho_a_2(N_states), vx_grad_rho_b_2(N_states), vx_grad_rho_a_b(N_states), vc_grad_rho_a_2(N_states), vc_grad_rho_b_2(N_states), vc_grad_rho_a_b(N_states))
!-----------Added------------
 allocate(rho(N_states), m_spin(N_states), zeta_m_n(N_states), g0(N_states), n2_UEG(N_states), n2_xc_UEG(N_states), ec_prime(N_states), ex_prime(N_states), beta_n_m_delta_n(N_states), delta_n_m_delta_n(N_states), gamma_n_m_delta_n(N_states), energy_x_sr_pbe_md_copy(N_states), energy_c_sr_pbe_md_copy(N_states))
!----------------------------

 allocate(rho_a(N_states), rho_b(N_states),grad_rho_a(3,N_states),grad_rho_b(3,N_states))
 allocate(grad_rho_a_2(N_states),grad_rho_b_2(N_states),grad_rho_a_b(N_states), ex(N_states), ec(N_states))
 
 energy_x_sr_pbe_md_copy = 0.d0
 energy_c_sr_pbe_md_copy = 0.d0

!----------------Constantes------------------ 
 pi = dacos(-1.d0)
 a = pi/2.d0
 b = 2*dsqrt(pi)*(2*dsqrt(2.d0) - 1.d0)/3.d0  
 c = 2*dsqrt(pi)*(1.d0 - dsqrt(2.d0))/3.d0
 mu = 0.5d0
!--------------------------------------------

 do istate = 1, N_states
  do i = 1, n_points_final_grid
   r(1) = final_grid_points(1,i)
   r(2) = final_grid_points(2,i)
   r(3) = final_grid_points(3,i)
   weight = final_weight_at_r_vector(i)
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

                             ! inputs
   call GGA_sr_type_functionals(r,rho_a,rho_b,grad_rho_a_2,grad_rho_b_2,grad_rho_a_b,                 &  ! outputs exchange
                             ex,vx_rho_a,vx_rho_b,vx_grad_rho_a_2,vx_grad_rho_b_2,vx_grad_rho_a_b, &  ! outputs correlation
                             ec,vc_rho_a,vc_rho_b,vc_grad_rho_a_2,vc_grad_rho_b_2,vc_grad_rho_a_b  )
   
   ! Define the variables below as arrows of size N_states
   rho = rho_a + rho_b
   m_spin = rho_a - rho_b
   zeta_m_n = m_spin/rho
   g0 = g0_UEG_mu(mu, rho_a, rho_b)
   n2_UEG = (rho**2)*(1.0d0 - zeta_m_n**2)*g0
   n2_xc_UEG = n2_UEG - rho**2
   
   ! a, b, c defined before the loop
   beta_n_m_delta_n = ec / (c*n2_UEG)
   delta_n_m_delta_n = - b*n2_UEG**2 / ex
   gamma_n_m_delta_n = ex / (a*n2_xc_UEG)
 
   ec_prime = ec / (1.0d0 + beta_n_m_delta_n*mu**3)
   ex_prime = ex / (1.0d0 + delta_n_m_delta_n*mu + gamma_n_m_delta_n*mu**2)
   
   energy_c_sr_pbe_md_copy += ec_prime * weight
   energy_x_sr_pbe_md_copy += ex_prime * weight
  
   !--------------- Print : Tested on H2O cc-pvdz ----------------

   !rho               : rho_copy_program_sr_pbe_md_2.dat
   !m_spin            : m_spin_copy_program_sr_pbe_md_2.dat
   !zeta_m_n          : zeta_copy_program_sr_pbe_md_2.dat
   !g0                : g0_copy_program_sr_pbe_md_2.dat
   !n2_UEG            : n2_UEG_copy_program_sr_pbe_md_2.dat
   !n2_xc_UEG         : n2xcUEG_copy_program_sr_pbe_md_2.dat
   !beta_n_m_delta_n  : beta_copy_program_sr_pbe_md_2.dat   !!!!
   !delta_n_m_delta_n : delta_copy_program_sr_pbe_md_2.dat  !!!!
   !gamma_n_m_delta_n : gamma_copy_program_sr_pbe_md_2.dat  
   !ec_prime          : ecprime_copy_program_sr_pbe_md_2.dat
   !ex_prime          : ex_prime_copy_program_sr_pbe_md_2.dat
  enddo
 enddo

end program
!END_PROVIDER
