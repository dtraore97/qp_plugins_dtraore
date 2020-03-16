program energy_x_c_md_test
! BEGIN_PROVIDER[double precision, energy_x_sr_pbe_md, (N_states) ]
!&BEGIN_PROVIDER[double precision, energy_c_sr_pbe_md, (N_states) ]
 implicit none
 BEGIN_DOC
! exchange/correlation energy with the short range pbe functional, multideterminantal form

! 11/03/20 Modified version of qp_plugins_eginer/stable/rsdft_cipsi/functionals/sr_pbe.irp.f
!-------- Parameters (alphabetical order)--------
! contrib_grad_xa(xb, ca,cb)
! energy_x/c_sr_pbe_md_copy : "Short-range exchange correlation density functionals" - Odense-Paris collaboration (11/03/20 version)
! g0                   : "Short-range exchange correlation density functionals" - Odense-Paris collaboration (11/03/20 version) -> from rsdft_ecmd/ueg_on_top.irp.f
! grad_rho_a(b)        : gradient of alpha(beta) spin density per state
! grad_rho_a(b)_2      : square of the gradient of alpha(beta) spin density per state
! grad_rho_a_b         : grad_rho_a*grad_rho_b
! rho_a(b)             : density of alpha(beta) spin per state
!------------------------
! Result for H2O : NaN -> parameters to check
! in sr_pbe.irp.f (Emmanuel's version) ex is used as 'ex' but it doesn't work in my code: I have to use 'ex(istate)' -> i should define delta, gamma and beta as arrows...
 END_DOC
 integer :: istate,i,j,m
 double precision :: r(3)
 !double precision :: mu,weight, delta_n_m_grad_n, gamma_n_m_grad_n, beta_n_m_grad_n, n2_xc_ueg, n2_ueg, a, b, c, g0, zeta_var, pi, g0_UEG_mu, g0_UEG_mu_inf
 double precision :: mu,weight, n2_xc_ueg, n2_ueg, a, b, c, g0, zeta_var, pi, g0_UEG_mu, g0_UEG_mu_inf 
 double precision, allocatable :: ex(:), ec(:), delta_n_m_grad_n(:), gamma_n_m_grad_n(:), beta_n_m_grad_n(:), energy_x_sr_pbe_md_copy(:), energy_c_sr_pbe_md_copy(:)
 double precision, allocatable :: rho_a(:),rho_b(:),grad_rho_a(:,:),grad_rho_b(:,:),grad_rho_a_2(:),grad_rho_b_2(:),grad_rho_a_b(:)
 double precision, allocatable :: contrib_grad_xa(:,:),contrib_grad_xb(:,:),contrib_grad_ca(:,:),contrib_grad_cb(:,:)
 double precision, allocatable :: vc_rho_a(:), vc_rho_b(:), vx_rho_a(:), vx_rho_b(:)
 double precision, allocatable :: vx_grad_rho_a_2(:), vx_grad_rho_b_2(:), vx_grad_rho_a_b(:), vc_grad_rho_a_2(:), vc_grad_rho_b_2(:), vc_grad_rho_a_b(:)
 allocate(vc_rho_a(N_states), vc_rho_b(N_states), vx_rho_a(N_states), vx_rho_b(N_states))
 allocate(vx_grad_rho_a_2(N_states), vx_grad_rho_b_2(N_states), vx_grad_rho_a_b(N_states), vc_grad_rho_a_2(N_states), vc_grad_rho_b_2(N_states), vc_grad_rho_a_b(N_states), energy_x_sr_pbe_md_copy(N_states), energy_c_sr_pbe_md_copy(N_states))


 allocate(rho_a(N_states), rho_b(N_states),grad_rho_a(3,N_states),grad_rho_b(3,N_states))
 allocate(grad_rho_a_2(N_states),grad_rho_b_2(N_states),grad_rho_a_b(N_states), ex(N_states), ec(N_states), delta_n_m_grad_n(N_states), gamma_n_m_grad_n(N_states), beta_n_m_grad_n(N_states))
 
 energy_x_sr_pbe_md_copy = 0.d0
 energy_c_sr_pbe_md_copy = 0.d0
 pi = dacos(-1.d0)
 a = pi/2.d0
 b = 2.d0*dsqrt(pi)*(2.d0*dsqrt(2.d0) - 1.d0)/3.d0
 c = 2.d0*dsqrt(pi)*(1.d0 - dsqrt(2.d0))/3.d0
 mu = 0.5d0

 print*, 'n-points_final_grid= ', n_points_final_grid
 print*, 'N_states= ', N_states
 
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
!GGA_sr_type_functionals from dft_utils_one_e/utils.irp.f
! -> ex_pbe_sr from dft_utils_one_e/exc_sr_pbe.irp.f
                             ! inputs
   call GGA_sr_type_functionals(r,rho_a,rho_b,grad_rho_a_2,grad_rho_b_2,grad_rho_a_b,                 &  ! outputs exchange
                             ex,vx_rho_a,vx_rho_b,vx_grad_rho_a_2,vx_grad_rho_b_2,vx_grad_rho_a_b, &  ! outputs correlation
                             ec,vc_rho_a,vc_rho_b,vc_grad_rho_a_2,vc_grad_rho_b_2,vc_grad_rho_a_b  )

   zeta_var = (rho_a(istate) - rho_b(istate))/(rho_a(istate) + rho_b(istate)) ! = 0 for spin multiplicity = 1
   ! zeta_var = mu / ( ((4d0/(9d0*pi))**(1d0/3d0)) * ((3d0 / (4d0*pi*(rho_a(istate)+rho_b(istate))))**(1d0/3d0)) )**(1d0)
  !  print*, 'zeta_var=', zeta_var
   g0 = g0_UEG_mu(mu, rho_a(istate),rho_b(istate)) !ordre de grandeur ok
   ! print*, 'g0=',g0
   !g0 = g0_UEG_mu_inf(rho_a,rho_b)
   n2_ueg = ((rho_a(istate) + rho_b(istate))**2)*(1d0 - zeta_var**2)*g0
   !print*, 'n2_ueg=', n2_ueg
   n2_xc_ueg = n2_ueg - (rho_a(istate) + rho_b(istate))**2
   !print*, 'n2_xc_ueg=', n2_xc_ueg
   gamma_n_m_grad_n(istate) = ex(istate)/(a*n2_xc_ueg)
   !print*, 'gamma_n_m_grad=', gamma_n_m_grad_n(istate)
   delta_n_m_grad_n(istate) = -(b*n2_ueg*gamma_n_m_grad_n(istate)**2)/ex(istate)
   beta_n_m_grad_n(istate)  = ec(istate)/(c*n2_ueg)
   
   energy_x_sr_pbe_md_copy(istate) += (ex(istate)/(1.d0 + delta_n_m_grad_n(istate)*mu + gamma_n_m_grad_n(istate)*mu**2)) * weight 
   energy_c_sr_pbe_md_copy(istate) += (ec(istate)/(1.d0 + beta_n_m_grad_n(istate)*mu**3)) * weight
   
  ! print*, 'energy_x_sr_pbe_md_copy(istate=',istate,'i=',i,')=',(ex(istate)/(1.d0 + delta_n_m_grad_n(istate)*mu + gamma_n_m_grad_n(istate)*mu**2)) * weight,' Sum is now: ', energy_x_sr_pbe_md_copy(istate)
  enddo
 enddo


end program