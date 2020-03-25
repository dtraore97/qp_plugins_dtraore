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
! g0_UEG_mu_inf  (g0)          : from rsdft_ecmd/ueg_on_top.irp.f
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
 call GGA_sr_type_functionals_mu_ok(1.d-12,r,rho_a,rho_b,grad_rho_a_2,grad_rho_b_2,grad_rho_a_b,         &  ! outputs exchange      
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

 BEGIN_PROVIDER[double precision, aos_sr_vc_alpha_sr_pbe_md_w  , (ao_num,n_points_final_grid,N_states)]
&BEGIN_PROVIDER[double precision, aos_sr_vc_beta_sr_pbe_md_w   , (ao_num,n_points_final_grid,N_states)]
&BEGIN_PROVIDER[double precision, aos_sr_vx_alpha_sr_pbe_md_w  , (ao_num,n_points_final_grid,N_states)]
&BEGIN_PROVIDER[double precision, aos_sr_vx_beta_sr_pbe_md_w   , (ao_num,n_points_final_grid,N_states)]
&BEGIN_PROVIDER[double precision, aos_dsr_vc_alpha_sr_pbe_md_w  , (ao_num,n_points_final_grid,N_states)]
&BEGIN_PROVIDER[double precision, aos_dsr_vc_beta_sr_pbe_md_w   ,  (ao_num,n_points_final_grid,N_states)]
&BEGIN_PROVIDER[double precision, aos_dsr_vx_alpha_sr_pbe_md_w  , (ao_num,n_points_final_grid,N_states)]
&BEGIN_PROVIDER[double precision, aos_dsr_vx_beta_sr_pbe_md_w   ,  (ao_num,n_points_final_grid,N_states)]
 implicit none
 BEGIN_DOC
! aos_sr_vxc_alpha_sr_pbe_md_w(j,i) = ao_i(r_j) * (v^x_alpha(r_j) + v^c_alpha(r_j)) * W(r_j)
 END_DOC
 integer :: istate,i,j,m
 double precision :: r(3)
 double precision :: mu,weight
 double precision, allocatable :: ex(:), ec(:)
 double precision, allocatable :: rho_a(:),rho_b(:),grad_rho_a(:,:),grad_rho_b(:,:),grad_rho_a_2(:),grad_rho_b_2(:),grad_rho_a_b(:)
 double precision, allocatable :: contrib_grad_xa(:,:),contrib_grad_xb(:,:),contrib_grad_ca(:,:),contrib_grad_cb(:,:)
 double precision, allocatable :: vc_rho_a(:), vc_rho_b(:), vx_rho_a(:), vx_rho_b(:)
 double precision, allocatable :: vx_grad_rho_a_2(:), vx_grad_rho_b_2(:), vx_grad_rho_a_b(:), vc_grad_rho_a_2(:), vc_grad_rho_b_2(:), vc_grad_rho_a_b(:)
 allocate(vc_rho_a(N_states), vc_rho_b(N_states), vx_rho_a(N_states), vx_rho_b(N_states))
 allocate(vx_grad_rho_a_2(N_states), vx_grad_rho_b_2(N_states), vx_grad_rho_a_b(N_states), vc_grad_rho_a_2(N_states), vc_grad_rho_b_2(N_states), vc_grad_rho_a_b(N_states))


 allocate(rho_a(N_states), rho_b(N_states),grad_rho_a(3,N_states),grad_rho_b(3,N_states))
 allocate(grad_rho_a_2(N_states),grad_rho_b_2(N_states),grad_rho_a_b(N_states), ex(N_states), ec(N_states))
 allocate(contrib_grad_xa(3,N_states),contrib_grad_xb(3,N_states),contrib_grad_ca(3,N_states),contrib_grad_cb(3,N_states))
 aos_dsr_vc_alpha_sr_pbe_md_w= 0.d0
 aos_dsr_vc_beta_sr_pbe_md_w = 0.d0
 aos_dsr_vx_alpha_sr_pbe_md_w= 0.d0
 aos_dsr_vx_beta_sr_pbe_md_w = 0.d0
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
   call GGA_sr_type_functionals_mu_ok(1.d-12,r,rho_a,rho_b,grad_rho_a_2,grad_rho_b_2,grad_rho_a_b,                 &  ! outputs exchange
                             ex,vx_rho_a,vx_rho_b,vx_grad_rho_a_2,vx_grad_rho_b_2,vx_grad_rho_a_b, &  ! outputs correlation
                             ec,vc_rho_a,vc_rho_b,vc_grad_rho_a_2,vc_grad_rho_b_2,vc_grad_rho_a_b  )
   vx_rho_a(istate) *= weight
   vc_rho_a(istate) *= weight
   vx_rho_b(istate) *= weight
   vc_rho_b(istate) *= weight
   do m= 1,3
    contrib_grad_ca(m,istate) = weight * (2.d0 * vc_grad_rho_a_2(istate) *  grad_rho_a(m,istate) + vc_grad_rho_a_b(istate)  * grad_rho_b(m,istate))
    contrib_grad_xa(m,istate) = weight * (2.d0 * vx_grad_rho_a_2(istate) *  grad_rho_a(m,istate) + vx_grad_rho_a_b(istate)  * grad_rho_b(m,istate))
    contrib_grad_cb(m,istate) = weight * (2.d0 * vc_grad_rho_b_2(istate) *  grad_rho_b(m,istate) + vc_grad_rho_a_b(istate)  * grad_rho_a(m,istate))
    contrib_grad_xb(m,istate) = weight * (2.d0 * vx_grad_rho_b_2(istate) *  grad_rho_b(m,istate) + vx_grad_rho_a_b(istate)  * grad_rho_a(m,istate))
   enddo
   do j = 1, ao_num
    aos_sr_vc_alpha_sr_pbe_md_w(j,i,istate) = vc_rho_a(istate) * aos_in_r_array(j,i)
    aos_sr_vc_beta_sr_pbe_md_w (j,i,istate) = vc_rho_b(istate) * aos_in_r_array(j,i)
    aos_sr_vx_alpha_sr_pbe_md_w(j,i,istate) = vx_rho_a(istate) * aos_in_r_array(j,i)
    aos_sr_vx_beta_sr_pbe_md_w (j,i,istate) = vx_rho_b(istate) * aos_in_r_array(j,i)
   enddo
   do j = 1, ao_num
    do m = 1,3
     aos_dsr_vc_alpha_sr_pbe_md_w(j,i,istate) += contrib_grad_ca(m,istate) * aos_grad_in_r_array_transp_xyz(m,j,i)
     aos_dsr_vc_beta_sr_pbe_md_w (j,i,istate) += contrib_grad_cb(m,istate) * aos_grad_in_r_array_transp_xyz(m,j,i)
     aos_dsr_vx_alpha_sr_pbe_md_w(j,i,istate) += contrib_grad_xa(m,istate) * aos_grad_in_r_array_transp_xyz(m,j,i)
     aos_dsr_vx_beta_sr_pbe_md_w (j,i,istate) += contrib_grad_xb(m,istate) * aos_grad_in_r_array_transp_xyz(m,j,i)
    enddo
   enddo
  enddo
 enddo

 END_PROVIDER


 BEGIN_PROVIDER [double precision, pot_sr_scal_x_alpha_ao_sr_pbe_md, (ao_num,ao_num,N_states)]
&BEGIN_PROVIDER [double precision, pot_sr_scal_c_alpha_ao_sr_pbe_md, (ao_num,ao_num,N_states)]
&BEGIN_PROVIDER [double precision, pot_sr_scal_x_beta_ao_sr_pbe_md, (ao_num,ao_num,N_states)]
&BEGIN_PROVIDER [double precision, pot_sr_scal_c_beta_ao_sr_pbe_md, (ao_num,ao_num,N_states)]
 implicit none
 integer                        :: istate
   BEGIN_DOC
   ! intermediate quantity for the calculation of the vxc potentials for the GGA functionals  related to the scalar part of the potential 
   END_DOC
   pot_sr_scal_c_alpha_ao_sr_pbe_md = 0.d0
   pot_sr_scal_x_alpha_ao_sr_pbe_md = 0.d0
   pot_sr_scal_c_beta_ao_sr_pbe_md = 0.d0
   pot_sr_scal_x_beta_ao_sr_pbe_md = 0.d0
   double precision               :: wall_1,wall_2
   call wall_time(wall_1)
   do istate = 1, N_states
     ! correlation alpha
     call dgemm('N','T',ao_num,ao_num,n_points_final_grid,1.d0,                                                                                       &
                 aos_sr_vc_alpha_sr_pbe_md_w(1,1,istate),size(aos_sr_vc_alpha_sr_pbe_md_w,1),                                                                   &
                 aos_in_r_array,size(aos_in_r_array,1),1.d0,                                                                                          &
                 pot_sr_scal_c_alpha_ao_sr_pbe_md(1,1,istate),size(pot_sr_scal_c_alpha_ao_sr_pbe_md,1))
     ! correlation beta
     call dgemm('N','T',ao_num,ao_num,n_points_final_grid,1.d0,                                                                                         &
                 aos_sr_vc_beta_sr_pbe_md_w(1,1,istate),size(aos_sr_vc_beta_sr_pbe_md_w,1),                                                                       &
                 aos_in_r_array,size(aos_in_r_array,1),1.d0,                                                                                            &
                 pot_sr_scal_c_beta_ao_sr_pbe_md(1,1,istate),size(pot_sr_scal_c_beta_ao_sr_pbe_md,1))
     ! exchange alpha
     call dgemm('N','T',ao_num,ao_num,n_points_final_grid,1.d0,                                                                                         &
                 aos_sr_vx_alpha_sr_pbe_md_w(1,1,istate),size(aos_sr_vx_alpha_sr_pbe_md_w,1),                                                                     &
                 aos_in_r_array,size(aos_in_r_array,1),1.d0,                                                                                            &
                 pot_sr_scal_x_alpha_ao_sr_pbe_md(1,1,istate),size(pot_sr_scal_x_alpha_ao_sr_pbe_md,1))
     ! exchange beta
     call dgemm('N','T',ao_num,ao_num,n_points_final_grid,1.d0,                                                                                            &
                 aos_sr_vx_beta_sr_pbe_md_w(1,1,istate),size(aos_sr_vx_beta_sr_pbe_md_w,1),                                                                          &
                 aos_in_r_array,size(aos_in_r_array,1),1.d0,                                                                                               &
                 pot_sr_scal_x_beta_ao_sr_pbe_md(1,1,istate), size(pot_sr_scal_x_beta_ao_sr_pbe_md,1))
 
   enddo
 call wall_time(wall_2)

END_PROVIDER 


 BEGIN_PROVIDER [double precision, pot_sr_grad_x_alpha_ao_sr_pbe_md,(ao_num,ao_num,N_states)]
&BEGIN_PROVIDER [double precision, pot_sr_grad_x_beta_ao_sr_pbe_md,(ao_num,ao_num,N_states)]
&BEGIN_PROVIDER [double precision, pot_sr_grad_c_alpha_ao_sr_pbe_md,(ao_num,ao_num,N_states)]
&BEGIN_PROVIDER [double precision, pot_sr_grad_c_beta_ao_sr_pbe_md,(ao_num,ao_num,N_states)]
   implicit none
   BEGIN_DOC
   ! intermediate quantity for the calculation of the vxc potentials for the GGA functionals  related to the gradienst of the density and orbitals 
   END_DOC
   integer                        :: istate
   double precision               :: wall_1,wall_2
   call wall_time(wall_1)
   pot_sr_grad_c_alpha_ao_sr_pbe_md = 0.d0
   pot_sr_grad_x_alpha_ao_sr_pbe_md = 0.d0
   pot_sr_grad_c_beta_ao_sr_pbe_md = 0.d0
   pot_sr_grad_x_beta_ao_sr_pbe_md = 0.d0
   do istate = 1, N_states
       ! correlation alpha
       call dgemm('N','N',ao_num,ao_num,n_points_final_grid,1.d0,                                                                                             &
                  aos_dsr_vc_alpha_sr_pbe_md_w(1,1,istate),size(aos_dsr_vc_alpha_sr_pbe_md_w,1),                                                                      &
                  aos_in_r_array_transp,size(aos_in_r_array_transp,1),1.d0,                                                                                &
                  pot_sr_grad_c_alpha_ao_sr_pbe_md(1,1,istate),size(pot_sr_grad_c_alpha_ao_sr_pbe_md,1))
       ! correlation beta
       call dgemm('N','N',ao_num,ao_num,n_points_final_grid,1.d0,                                                                                             &
                  aos_dsr_vc_beta_sr_pbe_md_w(1,1,istate),size(aos_dsr_vc_beta_sr_pbe_md_w,1),                                                                      &
                  aos_in_r_array_transp,size(aos_in_r_array_transp,1),1.d0,                                                                                &
                  pot_sr_grad_c_beta_ao_sr_pbe_md(1,1,istate),size(pot_sr_grad_c_beta_ao_sr_pbe_md,1))
       ! exchange alpha
       call dgemm('N','N',ao_num,ao_num,n_points_final_grid,1.d0,                                                                                             &
                  aos_dsr_vx_alpha_sr_pbe_md_w(1,1,istate),size(aos_dsr_vx_alpha_sr_pbe_md_w,1),                                                                      &
                  aos_in_r_array_transp,size(aos_in_r_array_transp,1),1.d0,                                                                                &
                  pot_sr_grad_x_alpha_ao_sr_pbe_md(1,1,istate),size(pot_sr_grad_x_alpha_ao_sr_pbe_md,1))
       ! exchange beta
       call dgemm('N','N',ao_num,ao_num,n_points_final_grid,1.d0,                                                                                             &
                  aos_dsr_vx_beta_sr_pbe_md_w(1,1,istate),size(aos_dsr_vx_beta_sr_pbe_md_w,1),                                                                      &
                  aos_in_r_array_transp,size(aos_in_r_array_transp,1),1.d0,                                                                                &
                  pot_sr_grad_x_beta_ao_sr_pbe_md(1,1,istate),size(pot_sr_grad_x_beta_ao_sr_pbe_md,1))
   enddo
   
 call wall_time(wall_2)

END_PROVIDER

 BEGIN_PROVIDER [double precision, potential_x_alpha_ao_sr_pbe_md,(ao_num,ao_num,N_states)]
&BEGIN_PROVIDER [double precision, potential_x_beta_ao_sr_pbe_md,(ao_num,ao_num,N_states)]
&BEGIN_PROVIDER [double precision, potential_c_alpha_ao_sr_pbe_md,(ao_num,ao_num,N_states)]
&BEGIN_PROVIDER [double precision, potential_c_beta_ao_sr_pbe_md,(ao_num,ao_num,N_states)]
   implicit none
 BEGIN_DOC
 ! exchange / correlation potential for alpha / beta electrons  with the Perdew-Burke-Ernzerhof GGA functional 
 END_DOC 
   integer :: i,j,istate
   do istate = 1, n_states 
    do i = 1, ao_num
     do j = 1, ao_num
      potential_x_alpha_ao_sr_pbe_md(j,i,istate) = pot_sr_scal_x_alpha_ao_sr_pbe_md(j,i,istate) + pot_sr_grad_x_alpha_ao_sr_pbe_md(j,i,istate) + pot_sr_grad_x_alpha_ao_sr_pbe_md(i,j,istate)
      potential_x_beta_ao_sr_pbe_md(j,i,istate) = pot_sr_scal_x_beta_ao_sr_pbe_md(j,i,istate) + pot_sr_grad_x_beta_ao_sr_pbe_md(j,i,istate) + pot_sr_grad_x_beta_ao_sr_pbe_md(i,j,istate)

      potential_c_alpha_ao_sr_pbe_md(j,i,istate) = pot_sr_scal_c_alpha_ao_sr_pbe_md(j,i,istate) + pot_sr_grad_c_alpha_ao_sr_pbe_md(j,i,istate) + pot_sr_grad_c_alpha_ao_sr_pbe_md(i,j,istate)
      potential_c_beta_ao_sr_pbe_md(j,i,istate) = pot_sr_scal_c_beta_ao_sr_pbe_md(j,i,istate) + pot_sr_grad_c_beta_ao_sr_pbe_md(j,i,istate) + pot_sr_grad_c_beta_ao_sr_pbe_md(i,j,istate)
     enddo
    enddo
   enddo

END_PROVIDER 


 BEGIN_PROVIDER[double precision, aos_sr_vxc_alpha_sr_pbe_md_w  , (ao_num,n_points_final_grid,N_states)]
&BEGIN_PROVIDER[double precision, aos_sr_vxc_beta_sr_pbe_md_w   , (ao_num,n_points_final_grid,N_states)]
&BEGIN_PROVIDER[double precision, aos_dsr_vxc_alpha_sr_pbe_md_w  , (ao_num,n_points_final_grid,N_states)]
&BEGIN_PROVIDER[double precision, aos_dsr_vxc_beta_sr_pbe_md_w   ,  (ao_num,n_points_final_grid,N_states)]
 implicit none
 BEGIN_DOC
! aos_sr_vxc_alpha_pbe_w(j,i) = ao_i(r_j) * (v^x_alpha(r_j) + v^c_alpha(r_j)) * W(r_j)
 END_DOC
 integer :: istate,i,j,m
 double precision :: r(3)
 double precision :: mu,weight
 double precision, allocatable :: ex(:), ec(:)
 double precision, allocatable :: rho_a(:),rho_b(:),grad_rho_a(:,:),grad_rho_b(:,:),grad_rho_a_2(:),grad_rho_b_2(:),grad_rho_a_b(:)
 double precision, allocatable :: contrib_grad_xa(:,:),contrib_grad_xb(:,:),contrib_grad_ca(:,:),contrib_grad_cb(:,:)
 double precision, allocatable :: vc_rho_a(:), vc_rho_b(:), vx_rho_a(:), vx_rho_b(:)
 double precision, allocatable :: vx_grad_rho_a_2(:), vx_grad_rho_b_2(:), vx_grad_rho_a_b(:), vc_grad_rho_a_2(:), vc_grad_rho_b_2(:), vc_grad_rho_a_b(:)
 allocate(vc_rho_a(N_states), vc_rho_b(N_states), vx_rho_a(N_states), vx_rho_b(N_states))
 allocate(vx_grad_rho_a_2(N_states), vx_grad_rho_b_2(N_states), vx_grad_rho_a_b(N_states), vc_grad_rho_a_2(N_states), vc_grad_rho_b_2(N_states), vc_grad_rho_a_b(N_states))


 allocate(rho_a(N_states), rho_b(N_states),grad_rho_a(3,N_states),grad_rho_b(3,N_states))
 allocate(grad_rho_a_2(N_states),grad_rho_b_2(N_states),grad_rho_a_b(N_states), ex(N_states), ec(N_states))
 allocate(contrib_grad_xa(3,N_states),contrib_grad_xb(3,N_states),contrib_grad_ca(3,N_states),contrib_grad_cb(3,N_states))

 aos_dsr_vxc_alpha_sr_pbe_md_w = 0.d0
 aos_dsr_vxc_beta_sr_pbe_md_w = 0.d0

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
   vx_rho_a(istate) *= weight
   vc_rho_a(istate) *= weight
   vx_rho_b(istate) *= weight
   vc_rho_b(istate) *= weight
   do m= 1,3
    contrib_grad_ca(m,istate) = weight * (2.d0 * vc_grad_rho_a_2(istate) *  grad_rho_a(m,istate) + vc_grad_rho_a_b(istate)  * grad_rho_b(m,istate))
    contrib_grad_xa(m,istate) = weight * (2.d0 * vx_grad_rho_a_2(istate) *  grad_rho_a(m,istate) + vx_grad_rho_a_b(istate)  * grad_rho_b(m,istate))
    contrib_grad_cb(m,istate) = weight * (2.d0 * vc_grad_rho_b_2(istate) *  grad_rho_b(m,istate) + vc_grad_rho_a_b(istate)  * grad_rho_a(m,istate))
    contrib_grad_xb(m,istate) = weight * (2.d0 * vx_grad_rho_b_2(istate) *  grad_rho_b(m,istate) + vx_grad_rho_a_b(istate)  * grad_rho_a(m,istate))
   enddo
   do j = 1, ao_num
    aos_sr_vxc_alpha_sr_pbe_md_w(j,i,istate) = ( vc_rho_a(istate) + vx_rho_a(istate) ) * aos_in_r_array(j,i)
    aos_sr_vxc_beta_sr_pbe_md_w (j,i,istate) = ( vc_rho_b(istate) + vx_rho_b(istate) ) * aos_in_r_array(j,i)
   enddo
   do j = 1, ao_num
    do m = 1,3
     aos_dsr_vxc_alpha_sr_pbe_md_w(j,i,istate) += ( contrib_grad_ca(m,istate) + contrib_grad_xa(m,istate) ) * aos_grad_in_r_array_transp_xyz(m,j,i)
     aos_dsr_vxc_beta_sr_pbe_md_w (j,i,istate) += ( contrib_grad_cb(m,istate) + contrib_grad_xb(m,istate) ) * aos_grad_in_r_array_transp_xyz(m,j,i)
    enddo
   enddo
  enddo
 enddo

 END_PROVIDER


 BEGIN_PROVIDER [double precision, pot_sr_scal_xc_alpha_ao_sr_pbe_md, (ao_num,ao_num,N_states)]
&BEGIN_PROVIDER [double precision, pot_sr_scal_xc_beta_ao_sr_pbe_md, (ao_num,ao_num,N_states)]
 implicit none
 integer                        :: istate
   BEGIN_DOC
   ! intermediate quantity for the calculation of the vxc potentials for the GGA functionals  related to the scalar part of the potential 
   END_DOC
   pot_sr_scal_xc_alpha_ao_sr_pbe_md = 0.d0
   pot_sr_scal_xc_beta_ao_sr_pbe_md = 0.d0
   double precision               :: wall_1,wall_2
   call wall_time(wall_1)
   do istate = 1, N_states
     ! exchange - correlation alpha
     call dgemm('N','T',ao_num,ao_num,n_points_final_grid,1.d0,                                                                                       &
                 aos_sr_vxc_alpha_sr_pbe_md_w(1,1,istate),size(aos_sr_vxc_alpha_sr_pbe_md_w,1),                                                                   &
                 aos_in_r_array,size(aos_in_r_array,1),1.d0,                                                                                          &
                 pot_sr_scal_xc_alpha_ao_sr_pbe_md(1,1,istate),size(pot_sr_scal_xc_alpha_ao_sr_pbe_md,1))
     ! exchange - correlation beta
     call dgemm('N','T',ao_num,ao_num,n_points_final_grid,1.d0,                                                                                         &
                 aos_sr_vxc_beta_sr_pbe_md_w(1,1,istate),size(aos_sr_vxc_beta_sr_pbe_md_w,1),                                                                       &
                 aos_in_r_array,size(aos_in_r_array,1),1.d0,                                                                                            &
                 pot_sr_scal_xc_beta_ao_sr_pbe_md(1,1,istate),size(pot_sr_scal_xc_beta_ao_sr_pbe_md,1))
   enddo
 call wall_time(wall_2)

END_PROVIDER 


 BEGIN_PROVIDER [double precision, pot_sr_grad_xc_alpha_ao_sr_pbe_md,(ao_num,ao_num,N_states)]
&BEGIN_PROVIDER [double precision, pot_sr_grad_xc_beta_ao_sr_pbe_md,(ao_num,ao_num,N_states)]
   implicit none
   BEGIN_DOC
   ! intermediate quantity for the calculation of the vxc potentials for the GGA functionals  related to the gradienst of the density and orbitals 
   END_DOC
   integer                        :: istate
   double precision               :: wall_1,wall_2
   call wall_time(wall_1)
   pot_sr_grad_xc_alpha_ao_sr_pbe_md = 0.d0
   pot_sr_grad_xc_beta_ao_sr_pbe_md = 0.d0
   do istate = 1, N_states
       ! exchange - correlation alpha
       call dgemm('N','N',ao_num,ao_num,n_points_final_grid,1.d0,                                                                                             &
                  aos_dsr_vxc_alpha_sr_pbe_md_w(1,1,istate),size(aos_dsr_vxc_alpha_sr_pbe_md_w,1),                                                                      &
                  aos_in_r_array_transp,size(aos_in_r_array_transp,1),1.d0,                                                                                &
                  pot_sr_grad_xc_alpha_ao_sr_pbe_md(1,1,istate),size(pot_sr_grad_xc_alpha_ao_sr_pbe_md,1))
       ! exchange - correlation beta
       call dgemm('N','N',ao_num,ao_num,n_points_final_grid,1.d0,                                                                                             &
                  aos_dsr_vxc_beta_sr_pbe_md_w(1,1,istate),size(aos_dsr_vxc_beta_sr_pbe_md_w,1),                                                                      &
                  aos_in_r_array_transp,size(aos_in_r_array_transp,1),1.d0,                                                                                &
                  pot_sr_grad_xc_beta_ao_sr_pbe_md(1,1,istate),size(pot_sr_grad_xc_beta_ao_sr_pbe_md,1))
   enddo
   
 call wall_time(wall_2)

END_PROVIDER

 BEGIN_PROVIDER [double precision, potential_xc_alpha_ao_sr_pbe_md,(ao_num,ao_num,N_states)]
&BEGIN_PROVIDER [double precision, potential_xc_beta_ao_sr_pbe_md,(ao_num,ao_num,N_states)]
   implicit none
 BEGIN_DOC
 ! exchange / correlation potential for alpha / beta electrons  with the Perdew-Burke-Ernzerhof GGA functional 
 END_DOC 
   integer :: i,j,istate
   do istate = 1, n_states 
    do i = 1, ao_num
     do j = 1, ao_num
      potential_xc_alpha_ao_sr_pbe_md(j,i,istate) = pot_sr_scal_xc_alpha_ao_sr_pbe_md(j,i,istate) + pot_sr_grad_xc_alpha_ao_sr_pbe_md(j,i,istate) + pot_sr_grad_xc_alpha_ao_sr_pbe_md(i,j,istate)
      potential_xc_beta_ao_sr_pbe_md(j,i,istate)  = pot_sr_scal_xc_beta_ao_sr_pbe_md(j,i,istate)  + pot_sr_grad_xc_beta_ao_sr_pbe_md(j,i,istate)  + pot_sr_grad_xc_beta_ao_sr_pbe_md(i,j,istate)
     enddo
    enddo
   enddo

END_PROVIDER 
