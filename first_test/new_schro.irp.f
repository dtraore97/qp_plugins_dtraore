program New_Schro
 implicit none
 read_wf = .True.
 touch read_wf
 
 BEGIN_DOC
 ! All integrals for sr_dft E0 calculations without lda_UEG approximation
 ! ---------- EHxc ---------
 ! 1) WB_sr:
 !    ok, tested
 !
 ! 2) VB_sr:
 !    ok, to test
 !
 ! ---------- Elr ----------- (looking for PsiMu)
 ! 3) T_mu:
 !    See : mo_one_e_ints/mo_one_e_ints.irp.f ---->  mo_one_e_integrals
 !
 ! 4) Wee_lr:
 !    I'm looking for this integral in the plugins
 !    Need to find erf.
 !
 ! 6) ------- To add ! Vne -----------
 !    Check if One_e_integrals contains Vne AND T  
 ! 
 ! 5) File with integrals
 END_DOC

 double precision, allocatable :: w_B(:), beta(:), mo_r1(:), mo_r2(:), r1(:), e_pbe(:) ! w_B: distribution - Signification à revoir (eq. A8 - article 149 2018), beta: eq. 14b - article J. Phys. Chem. 2019, mo_ri: orbitales en ri // eq. A6 opérateur creations et annihilations?? 
 double precision :: c_1, Sum_w, Sum_w_lim, pi, mu, theta, dtheta, thetamax, r_initial, mu_max, two_dm, two_dm_in_r_diata, get_two_e_integral, bielec_integral ! Sum_w: somme 
 integer :: mo_num_i, mo_num_j, mo_num_k, mo_num_l, i, m, istate, j ! Numéro des orbitales
 
 allocate(w_B(n_points_final_grid), v_B(n_points_final_grid), beta(n_points_final_grid), mo_r1(mo_num), mo_r2(mo_num), r1(3), e_pbe(N_states))

 
 pi = dacos(-1.d0)
 ! mu = 1.d0
 c_1 = 3.d0 / (2.d0 * dsqrt(pi)*(1.d0 - dsqrt(2.d0)))
 mu = 0.d0
 mo_num_i = 1
 mo_num_j = 1
 mo_num_k = 1 
 mo_num_l = 1 
 istate = 1
 j = 37
 print*, 'N_states= ', N_states, 'n_points_final_grid= ', n_points_final_grid

! -------------- WB = Sum(space) Psi_r(i)Psi_r(j) wB Psi_r(k)Psi_r(l)----------------
!do while (mu < 10.d0)
 Sum_w = 0.d0
 Sum_w_lim = 0.d0
 do i = 1, n_points_final_grid
  do m = 1, 3
   r1(m) = final_grid_points(m,i) ! ------ ????? ------
  enddo

  call ec_pbe_compact(r1,mu,e_pbe)
  call beta_from_on_top_grad_on_top_e_pbe(r1, mu, beta)
  call give_all_mos_at_r(r1,mo_r1)
  two_dm = two_dm_in_r_diata(r1,r1,istate)
  write(37,*)mu, i, Beta(istate)
 
  !w_B(i) = e_pbe(istate)**(2) * mu**(3) * c_1 / (two_dm**(2) *(1.d0 + beta(istate)*mu**(3))**(2)) ! ATTENTION e_pbe : tableau sur le nombre d'états
  mu = mu_of_r_vector(i)
  w_B(i) = e_pbe(istate)**(2) * mu**(3) * c_1 / (two_dm**(2) *(1.d0 + beta(istate)*mu**(3))**(2))
  Sum_w += final_weight_at_r_vector(i)*mo_r1(mo_num_i)*mo_r1(mo_num_j)*mo_r2(mo_num_k)*mo_r2(mo_num_l)*w_B(i) ! Final_weight_at_r_vector
  !print*,'final_weight_at_r_vector(',i,')= ', final_weight_at_r_vector(i)
  Sum_w_lim += final_weight_at_r_vector(i)*mo_r1(mo_num_i)*mo_r1(mo_num_j)*mo_r1(mo_num_k)*mo_r1(mo_num_l)/(c_1*mu**3)
   if (mu < 1.5d0 .AND. mu > 0.7d0) then
  !  write(j,*)mu, r1(1), r1(2), r1(3), Sum_w, Sum_w_lim
   end if
 enddo
 j+=1
 !print *,'sum = ', Sum_w
  bielec_integral = get_two_e_integral(1,1,1,1,mo_integrals_map)
! write(36,*)mu, Beta(1), Sum_w, Sum_w_lim, bielec_integral
 write(36,*) Sum_w, Sum_w_lim, bielec_integral
!  mu += 0.1d0

 ! hmono = 2.d0 * mo_one_e_integrals(2,2)
  !print*,'e = ', hmono + bielec_integral
!enddo

! ---------------- VB _ short range ----------------
double precision, allocatable :: VB_n(:),dm_a(:),dm_b(:)
double precision :: Sum_v, cx
Sum_v = 0.d0
cx = 
allocate(VB_n(n_points_final_grid), dm_a(N_states), dm_b(N_states))
do i = 1, n_points_final_grid
 do m = 1, 3
  r1(m) = final_grid_points(m,i)
 enddo
  call dm_dft_alpha_beta_at_r(r1,dm_a,dm_b)
  call give_all_mos_at_r(r1,mo_r2)

VB_n(i) = 4.d0/3.d0 * cx * (dm_a(mo_num_j) + dm_a(mo_num_j))**(1.d0/3.d0)
Sum_v += finale_weight_at_r_vector(i)*VB_n(i)*mo_r1(mo_num_j)*mo_r1(mo_num_j)
enddo

!--------------- 1 electron sum (kinetic energy) ---------
double precision :: Sum_Ek
double precision, allocatable :: Ek_i(:)
allocate(Ek_i(n_points_final_grid))

Sum_Ek = 0.d0

do i = 1, n_points_final_grid
 do m = 1, 3
  r1(m) = final_grid_points(m,i)
 enddo

 Ek_i(i) = mo_one_e_integrals(i,j) 
 Sum_Ek += Ek_i(i)*mo_r1(mo_num_i)*mo_r1(mo_num_j)

enddo

!------------  Long-range Coulomb Energy --------------
double precision, allocatable :: Wee_lr(i)
double precision :: Sum_We_lr
allocate(Wee_lr(n_points_final_grid))

do i = 1, n_points_final_grid
 do m = 1, 3
  r1(m) = final_grid_points(m,i)
 enddo

Wee_lr = 
Sum_We_lr += final_weight_at_r_vector(i)*Wee_lr(i)*mo_r1(mo_num_i)*mo_r1(mo_num_j)*mo_r2(mo_num_k)*mo_r2(mo_num_l)

enddo

!----------- File with integrals ----------------
 write(36,*) Sum_w, Sum_Ek, Sum_Weelr, Sum_v !Sum_w = Wee_sr

end
