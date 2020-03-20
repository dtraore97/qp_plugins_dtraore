program energy_x_c_md_test_2
 implicit none
 BEGIN_DOC
! Last modification : Thurs., March 19th, 6.25pm
! WARNING: This version contains some mistakes that will be checked tomorrow
! 
! Wedn., March 18th : arrays for r and rho : check with Julien's Mathematica's program
! exchange/correlation energy with the short range pbe functional, multideterminantal form (VERSION 2)
! 11/03/20 Modified version of qp_plugins_eginer/stable/rsdft_cipsi/functionals/sr_pbe.irp.f
! This program aim to test the providers from sr_pbe_version_2.irp.f
! ------------------------------------------PARAMETERS-------------------------------------------------
! a, b and c                   : Constants                                               eq. 45 and 47 
! beta_n_m_grad_n              :                                                         eq. 50
! contrib_grad                 : 
! delta_n_m_grad_n             :                                                         eq. 57
! energy_c_sr_pbe_copy         : Non md case
! energy_x_sr_pbe_copy         : Non md case 
! energy_c_sr_pbe_md_copy      : Sum(ec_prime*weight)                                    eq. 48
! energy_x_sr_pbe_md_copy      : Sum(ex_prime*weight)                                    eq. 53
! ex & ec                      : ec_PBE & ex_PBE                                         eq. 49 and 54
! ex_prime & ec_prime          : exchange and correlation energies perelementaty volume  eq. 54 and 49
! gamma_n_m_grad_n             :                                                         eq. 55
! GGA_sr_type_functionals      : Emmanuel's                                              dft_utils_one_e/utils.irp.f
! GGA_sr_type_functionals_mu   : Modified version of GGA_sr_type_functionals where mu is intent(in) so we can have Ecsr value for a given mu
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
 integer :: test_r, k, p
 double precision :: pi, a, b, c, g0_UEG_mu, thr, r_norm
! double precision, allocatable :: rho(:), m_spin(:), zeta_m_n(:), g0(:), n2_UEG(:), n2_xc_UEG(:), ec_prime(:), ex_prime(:), beta_n_m_delta_n(:), delta_n_m_delta_n(:), gamma_n_m_delta_n(:), energy_x_sr_pbe_md_copy(:), energy_c_sr_pbe_md_copy(:)
 double precision :: rho, m_spin, zeta_m_n, g0, n2_UEG, n2_xc_UEG, ec_prime, ex_prime, beta_n_m_delta_n, delta_n_m_delta_n, gamma_n_m_delta_n
 double precision, allocatable :: energy_x_sr_pbe_md_copy(:), energy_c_sr_pbe_md_copy(:), energy_x_sr_pbe_copy(:), energy_c_sr_pbe_copy(:), r_norm_prec(:), mu_array(:), ex_pbe(:), ec_pbe(:)
!---------------------------

 allocate(vc_rho_a(N_states), vc_rho_b(N_states), vx_rho_a(N_states), vx_rho_b(N_states))
 allocate(vx_grad_rho_a_2(N_states), vx_grad_rho_b_2(N_states), vx_grad_rho_a_b(N_states), vc_grad_rho_a_2(N_states), vc_grad_rho_b_2(N_states), vc_grad_rho_a_b(N_states))
!-----------Added------------
 allocate(energy_x_sr_pbe_md_copy(N_states), energy_c_sr_pbe_md_copy(N_states), energy_x_sr_pbe_copy(N_states), energy_c_sr_pbe_copy(N_states),r_norm_prec(n_points_final_grid), mu_array(20), ex_pbe(N_states), ec_pbe(N_states))
!----------------------------

 allocate(rho_a(N_states), rho_b(N_states),grad_rho_a(3,N_states),grad_rho_b(3,N_states))
 allocate(grad_rho_a_2(N_states),grad_rho_b_2(N_states),grad_rho_a_b(N_states), ex(N_states), ec(N_states))
 
 energy_x_sr_pbe_md_copy = 0.d0
 energy_c_sr_pbe_md_copy = 0.d0
 energy_x_sr_pbe_copy = 0.d0
 energy_c_sr_pbe_copy = 0.d0


 r_norm_prec = 0.d0
!----------------Constantes------------------ 
 pi = dacos(-1.d0)
 a = pi/2.d0
 b = 2*dsqrt(pi)*(2*dsqrt(2.d0) - 1.d0)/3.d0  
 c = 2*dsqrt(pi)*(1.d0 - dsqrt(2.d0))/3.d0
 !mu = 1.d12
 mu_array = (/ 0.d0, 0.125d0, 0.25d0, 0.375d0, 0.5d0, 0.625d0, 0.75d0, 0.875d0, 1.d0, 1.5d0, 2.d0, 2.5d0, 3.d0, 4.d0, 5.d0, 6.d0, 7.d0, 8.d0, 9.d0, 10.d0 /)
 thr = 1.d-12
 r_norm_prec = 1.d-12
!--------------------------------------------
!-------------------FILES--------------------
! rho(r) file for a given mu
OPEN(UNIT=14, FILE='rho_r_he_ccpvdz.dat')
! With multideterminantal extension
OPEN(UNIT=10, FILE='mu_ecsr_md_he_ccpvdz20.dat')
OPEN(UNIT=11, FILE='mu_exsr_md_he_ccpvdz20.dat')
! Without md extension
OPEN(UNIT=12, FILE='mu_ecsr_he_ccpvdz20.dat')
OPEN(UNIT=13, FILE='mu_exsr_he_ccpvdz20.dat')
!--------------------------------------------


do p = 1, 20  ! loop over mu_array
 mu = mu_array(p)
 do i = 1, n_points_final_grid ! do1
  do istate = 1, N_states ! do2
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
   do m = 1, 3 ! do3
    grad_rho_a_2(istate) += grad_rho_a(m,istate) * grad_rho_a(m,istate)
    grad_rho_b_2(istate) += grad_rho_b(m,istate) * grad_rho_b(m,istate)
    grad_rho_a_b(istate) += grad_rho_a(m,istate) * grad_rho_b(m,istate)
   enddo ! do3
  enddo !enddo2
                             ! inputs
   call GGA_sr_type_functionals_mu(mu,r,rho_a,rho_b,grad_rho_a_2,grad_rho_b_2,grad_rho_a_b,                 &  ! outputs exchange
                             ex_pbe,vx_rho_a,vx_rho_b,vx_grad_rho_a_2,vx_grad_rho_b_2,vx_grad_rho_a_b, &  ! outputs correlation
                             ec_pbe,vc_rho_a,vc_rho_b,vc_grad_rho_a_2,vc_grad_rho_b_2,vc_grad_rho_a_b  )
   
   call GGA_sr_type_functionals(r,rho_a,rho_b,grad_rho_a_2,grad_rho_b_2,grad_rho_a_b,                 &  ! outputs exchange
                             ex,vx_rho_a,vx_rho_b,vx_grad_rho_a_2,vx_grad_rho_b_2,vx_grad_rho_a_b, &  ! outputs correlation
                             ec,vc_rho_a,vc_rho_b,vc_grad_rho_a_2,vc_grad_rho_b_2,vc_grad_rho_a_b  )

   do istate = 1, N_states ! do4
     
    if(dabs(ex(istate)).lt.thr)then
       ex(istate) = 1.d-12
     endif

     rho = rho_a(istate) + rho_b(istate)
     if(dabs(rho).lt.thr)then
       rho = 1.d-12
     endif
     
     m_spin = rho_a(istate) - rho_b(istate) !0.0000
     if(dabs(m_spin).lt.thr)then
       m_spin = 1.d-12
     endif
      
     zeta_m_n = m_spin/rho !0.00000
     if(dabs(zeta_m_n).lt.thr)then
       zeta_m_n = 1.d-12
     endif

     g0 = g0_UEG_mu(mu, rho_a(istate), rho_b(istate))
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
     
     energy_c_sr_pbe_copy(istate) += ec_pbe(istate) * weight
     energy_x_sr_pbe_copy(istate) += ex_pbe(istate) * weight 
     energy_c_sr_pbe_md_copy(istate) += ec_prime * weight
     energy_x_sr_pbe_md_copy(istate) += ex_prime * weight
   enddo !do4
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

!-------------------- rho(r) bloc ------------------------
   r_norm = dsqrt(r(1)**2 + r(2)**2 + r(3)**2)
   
   test_r=0

   do k=1, i
     if(r_norm == r_norm_prec(k))then
       test_r = 1
       exit
     endif

     if(dabs(r_norm - r_norm_prec(i-1)) < 1.d-10)then
       test_r = 1
     endif
   enddo
   
   r_norm_prec(i) = r_norm
   
   if(test_r==0)then
     write(14,*) mu,' ', r_norm, ' ',rho  
     print*, r_norm, rho
   endif
!--------------------------------------------------------

 enddo
 
 write(10, *) mu, ' ', energy_c_sr_pbe_md_copy(1)
 write(11, *) mu, ' ', energy_x_sr_pbe_md_copy(1)
 write(12, *) mu, ' ', energy_c_sr_pbe_copy(1)
 write(13, *) mu, ' ', energy_x_sr_pbe_copy(1)

enddo
end program
!END_PROVIDER
