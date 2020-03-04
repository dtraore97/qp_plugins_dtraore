program effectiv_potential
 implicit none
 read_wf = .True.
 touch read_wf

 BEGIN_DOC
 ! Here we compute the 'complementary part' of the effectiv interaction !
 ! ---- Variables ----
 ! -- Dimension > 1 --
 ! - w_B(n_points_final_grid) = Delta(EbarB)/Delta(on top pair density) ! - Beta = function of e_pbe and on top pair density
 ! - mo_r1(mo_num) and mo_r2(mo_num) = molecular orbitals (or just orbitals as overlaps can occur) at r1 and r2 ; initialize with 'give_all_mos_at_r(r,mo_r)
 ! - r1(3) : coordinates 
 ! - e_pbe : 'e^pbe_c', from 'ec_pbe_compact(r, mu, e_pbe)
 ! 
 ! -- Constant --
 ! - Sum_w = sum_over_the_grid_k' ( weight*phi_i(rk')*phi_j(rk')*phi_k(rk')*phi_l(rk')*w_B(rk') ) 
 ! - Sum_w_lim = sum_over_the_grid_k' ( weight*phi_i(rk')*phi_j(rk')*phi_k(rk')*phi_l(rk')*(cst/mu**3) ) : to check the high mu behaviour (lim mu tend to infinity)
 ! - c = [2*sqrt(pi)*(1-sqrt(2))]/3
 ! ------------------- 
 END_DOC

 double precision, allocatable :: w_B(:), beta(:), mo_r1(:), mo_r2(:), r1(:), e_pbe(:), mu_table(:)
 double precision :: c, Sum_w, Sum_w_lim, pi, mu, two_dm, two_dm_in_r_diata, pas, j_dp
 integer :: i, j, istate, m, mo_num_i, mo_num_j, mo_num_k, mo_num_l
 allocate(w_B(n_points_final_grid), beta(n_points_final_grid), mo_r1(mo_num), mo_r2(mo_num), r1(3), e_pbe(N_states), mu_table(10))

 
 pi = dacos(-1.d0)
 c = (2.d0 * dsqrt(pi)*(1.d0 - dsqrt(2.d0))) / 3.d0
 Sum_w = 0.d0
 Sum_w_lim = 0.d0
 istate = 1
 mo_num_i = 1
 mo_num_j = 1 
 mo_num_k = 1
 mo_num_l = 1
 r1(1) = 0.d0
 r1(2) = 0.d0
 r1(3) = 0.d0
 pas = 0.5
 !do i = 1, n_points_final_grid

 j_dp = 0.d0
 do i = 1, 10
   mu_table(i) = j_dp
   j_dp += 1.d0
 enddo

 do j = 40, 50
 mu = mu_table(j)
 r1(1) = 0.d0
  do i = 1, 500 
   !do m = 1, 3
   !  r1(m) = final_grid_points(m,i)
   ! enddo 
   ! ---- mu = constante then -------
   ! mu = 0.d0
   ! ---- mu = funtion of r then ----
   ! mu = mu_of_r_vector(i)
 
   two_dm = two_dm_in_r_diata(r1,r1, istate) ! on_top pair density
   call beta_from_on_top_grad_on_top_e_pbe(r1, mu, beta)
   call ec_pbe_compact(r1,mu, e_pbe)
   call give_all_mos_at_r(r1,mo_r1)

 
   ! ---- w_B calculation -----------
   w_B(i) = e_pbe(istate)**(2) * mu**(3) / ( two_dm**(2) * (1.d0 + beta(istate)*mu**(3))**(2) * c)
   ! ---- Sums ----------------------
   Sum_w += final_weight_at_r_vector(i)*mo_r1(mo_num_i)*mo_r1(mo_num_j)*mo_r1(mo_num_k)*mo_r1(mo_num_l)*w_B(i)
   Sum_w_lim += final_weight_at_r_vector(i)*mo_r1(mo_num_i)*mo_r1(mo_num_j)*mo_r1(mo_num_k)*mo_r1(mo_num_l)*c/(mu**3) 

   r1(1) += pas
   write(j,*) r1(1), w_B(i), Sum_w, Sum_w_lim
  enddo
 enddo
 ! write(0403,*) Sum_w, Sum_w_lim

end
