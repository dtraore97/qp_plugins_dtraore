
 BEGIN_PROVIDER[double precision, energy_x_sr_pbe_2, (N_states) ]
&BEGIN_PROVIDER[double precision, energy_c_sr_pbe_2, (N_states) ]
 implicit none
 BEGIN_DOC
! exchange/correlation energy with the short range pbe functional
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
 double precision :: pi, a, b, c
 double precision, allocatable :: n(:), m_spin(:), zeta_m_n(:), g0(:), n2_UEG(:), n2_xc_UEG(:), ec_prime(:), ex_prime(:), beta_n_m_delta_n(:), delta_n_m_delta_n(:), gamma_n_m_delta_n(:) 
!---------------------------

 allocate(vc_rho_a(N_states), vc_rho_b(N_states), vx_rho_a(N_states), vx_rho_b(N_states))
 allocate(vx_grad_rho_a_2(N_states), vx_grad_rho_b_2(N_states), vx_grad_rho_a_b(N_states), vc_grad_rho_a_2(N_states), vc_grad_rho_b_2(N_states), vc_grad_rho_a_b(N_states))
!-----------Added------------
 !allocate()
!----------------------------

 allocate(rho_a(N_states), rho_b(N_states),grad_rho_a(3,N_states),grad_rho_b(3,N_states))
 allocate(grad_rho_a_2(N_states),grad_rho_b_2(N_states),grad_rho_a_b(N_states), ex(N_states), ec(N_states))
 
 energy_x_sr_pbe = 0.d0
 energy_c_sr_pbe = 0.d0

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
   
   ! Define the variables below as arrows of size Nstates
   n = rho_a + rho_b
   m_spin = rho_a - rho_b
   zeta_m_n = m_spin/n
   g0 = g0_mu_UEG(mu, rho_a, rho_b) ! g0 is arrow of size n (???)
   n2_UEG = (n**2)*(1.0d0 - zeta_m_n**2)*g0
   n2_xc_UEG = n2_UEG - n**2
   
   ! a, b, c defined before the loop
   beta_n_m_delta_n = ec / (c*n2_UEG)
   delta_n_m_delta_n = - b*n2_UEG**2 / ex
   gamma_n_m_delta_n = ex / (a*n2_xc_UEG)
 
   ec_prime = ec / (1.0d0 + beta_n_m_delta_n*mu**3)
   ex_prime = ex / (1.0d0 + delta_n_m_delta_n*mu + gamma_n_m_delta_n*mu**2)
   
   energy_c_sr_pbe_2 += ec_prime * weight
   energy_x_sr_pbe_2 += ex_prime * weight
  enddo
 enddo


END_PROVIDER


