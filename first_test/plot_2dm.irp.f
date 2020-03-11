program plot_2_dm
 implicit none
 read_wf = .True.
 touch read_wf
 double precision :: r(3), r1(3), r2(3)
 double precision, allocatable :: mos_array(:),dm_a(:),dm_b(:),e_pbe(:)
 allocate(mos_array(mo_num),dm_a(N_states),dm_b(N_states),e_pbe(N_states))
 integer :: i,nx, ntheta,m 
 double precision :: two_dm_in_r_diata , two_dm, r_initial

! double precision :: dx , xmas
! xmax = 5.d0
! nx = 1000
! dx = xmax/dble(nx)
! print*,'xmax,nx,dx',xmax,nx,dx
! r(1) = 0.d0
! r(2) = 0.d0
! r(3) = 0.d0
 
! Coordonnées cylindriques:  
 double precision :: dtheta, dthetamax, theta, alpha, pi
 pi = dacos(-1.d0)
 dthetamax = 2.d0*pi
 ntheta = 500
 dtheta = dthetamax / dble(ntheta)
 r_initial = 0.50
 theta = 0.d0
 ! Travelling electron r1
 r1(1) = r_initial*dsin(theta)
 r1(2) = r_initial*dcos(theta)
 r1(3) = 0.d0
 ! Fixed electron r2
 r2(1) = 0.d0
 r2(2) = r_initial
 r2(3) = 0.d0
 ! alpha = pi/2.d0 - acos((r_initial/2.d0)*((cos(theta)*cos(theta)) + (1.d0 + sin(theta)*(1.d0 + sin(theta))))) 
 integer :: istate
 istate = 1
! do i = 1, nx
!  call give_all_mos_at_r(r,mos_array)
!  call dm_dft_alpha_beta_at_r(r,dm_a,dm_b)
!  two_dm = two_dm_in_r_diata(r,r,istate)
!  write(33,*)r(1),dm_a(istate)*dm_b(istate),two_dm
!  r(1) += dx
! enddo


 double precision :: mu
 mu = 0.d0
 write(34, *) 'theta   ', 'r1(1)   ', 'r1(2)   ', 'dsqrt(r1(1)**2 + r1(2)**2)    ', 'dma*dmb    ', 'n2     ', 'mos(1)**4    ', 'e_pbe'
 do i = 1, ntheta
  call give_all_mos_at_r(r1,mos_array)
  call dm_dft_alpha_beta_at_r(r1,dm_a,dm_b)
  two_dm = two_dm_in_r_diata(r1,r2,istate)
  call ec_pbe_compact(r1,mu,e_pbe)
  ! write(34,*)r1(1), r1(2), r2(1), r2(2), theta , dm_a(istate)*dm_b(istate) , two_dm
  write(34,'(100(F16.10,X))')theta, r1(1), r1(2), dsqrt(r1(1)**2 + r1(2)**2), dm_a(istate)*dm_b(istate) , two_dm,mos_array(1)**4.d0,e_pbe(istate)
 ! write(34,*)i, r(1), r(2), theta
  theta += dtheta
  ! alpha = pi/2.d0 - acos((r_initial/2.d0)*((cos(theta)*cos(theta)) + (1.d0 + sin(theta)*(1.d0 + sin(theta))))) 
  r1(1) = r_initial * dsin(theta)
  r1(2) = r_initial * dcos(theta)
 enddo

! Number of electrons calculation
 double precision :: accu, weight
 !integer :: n_points_final_grid
 !n_points_final_grid = 500
 accu = 0.d0
 do i = 1, n_points_final_grid
  do m = 1, 3
   r1(m) = final_grid_points(m,i) ! ------ ????? ------
  enddo 
  weight = final_weight_at_r_vector(i)

  call dm_dft_alpha_beta_at_r(r1,dm_a,dm_b)
  accu += weight * (dm_a(istate)+dm_b(istate))
 enddo
 print*,'accu = ',accu

! ------- Correction calculation ----------
! double precision, allocatable :: w_B(:), beta(:), mo_r1(:), mo_r2(:) ! w_B: distribution - Signification à revoir (eq. A8 - article 149 2018), beta: eq. 14b - article J. Phys. Chem. 2019, mo_ri: orbitales en ri // eq. A6 opérateur creations et annihilations?? 
! double precision :: c, Sum_w ! Sum_: somme 
! integer :: mo_num_i, mo_num_j, mo_num_k, mo_num_l ! Numéro des orbitales
! allocate(w_B(n_points_final_grid), beta(n_points_final_grid), mo_r1(n_points_final_grid), mo_r2(n_points_final_grid))
 
! c = 3.d0 / (2.d0 * dsqrt(pi)*(1.d0 - dsqrt(2.d0)))
! Sum_w = 0
! mo_num_i = 1
! mo_num_j = 1 
! mo_num_k = 1 
! mo_num_l = 1 

! do i = 1, n_points_final_grid
!  call ec_pbe_compact(r1,mu,e_pbe)
!  call beta_from_on_top_grad_on_top_e_pbe(r1, mu, Beta)
!  call give_all_mos_at_r(r1,mo_num_i)
!  call give_all_mos_at_r(r1,mo_num_j)  
!  call give_all_mos_at_r(r1,mo_num_k)  
!  call give_all_mos_at_r(r1,mo_num_l)   
  
!  w_B(i) = e_pbe(i)**(2) * mu**(3) / (1.d0 + Beta(i)*mu**(3))**(2) 
!  Sum_w += final_weight_at_r_vector(i) * mo_r1(mo_num_i)*mo_r1(mo_num_j)*mo_r1(mo_num_k)*mo_r1(mo_num_l)*w_B(i) 
  
! enddo
! print *,'sum = ', Sum_w
! print*,'youpi'
end
