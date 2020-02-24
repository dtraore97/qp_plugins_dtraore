double precision function two_dm_in_r_diata(r1,r2,istate)
 implicit none
 BEGIN_DOC
 ! two body density evaluated ar two points in real space 
 END_DOC
 integer, intent(in) :: istate
 double precision, intent(in) :: r1(3),r2(3)
 integer :: i,j,k,l
 double precision, allocatable :: mos_array_r1(:), mos_array_r2(:)
 allocate(mos_array_r2(mo_num), mos_array_r1(mo_num))
 call give_all_mos_at_r(r1,mos_array_r1)
 call give_all_mos_at_r(r2,mos_array_r2)
 two_dm_in_r_diata = 0.d0
 do l = 1, mo_num
  do k = 1, mo_num
    do j = 1, mo_num
     do i = 1, mo_num
     !                                                   1 2 1 2 
     two_dm_in_r_diata += all_states_act_two_rdm_alpha_beta_mo(i,j,k,l,istate) * mos_array_r1(i) * mos_array_r1(k) * mos_array_r2(j) * mos_array_r2(l)
    enddo
   enddo
  enddo
 enddo
 two_dm_in_r_diata = max(two_dm_in_r_diata,1.d-15)
end



subroutine ec_pbe_compact(r,mu,e_pbe)
 implicit none
 double precision, intent(in) :: mu
 double precision, intent(in) :: r(3) 
 double precision, intent(out) :: e_pbe(N_states)
 BEGIN_DOC
 ! routine that computes the pbe correlation energy for a given mu at a given point in r
 !
 ! it uses the density of the current wave function
 END_DOC
 double precision :: aos_array(ao_num), grad_aos_array(3,ao_num)
 double precision :: rho_a(N_states),rho_b(N_states)
 double precision :: grad_rho_a(3,N_states),grad_rho_b(3,N_states)
 double precision :: grad_rho_a_2(N_states),grad_rho_b_2(N_states),grad_rho_a_b(N_states)
 double precision :: rhoc,rhoo,sigmacc,sigmaco,sigmaoo,vrhoc,vrhoo,vsigmacc,vsigmaco,vsigmaoo
 integer :: m, istate
 call density_and_grad_alpha_beta_and_all_aos_and_grad_aos_at_r(r,rho_a,rho_b, grad_rho_a, grad_rho_b, aos_array, grad_aos_array)
 grad_rho_a_2 = 0.d0
 grad_rho_b_2 = 0.d0
 grad_rho_a_b = 0.d0
 do istate = 1, N_states
  do m = 1, 3
   grad_rho_a_2(istate) += grad_rho_a(m,istate)*grad_rho_a(m,istate)
   grad_rho_b_2(istate) += grad_rho_b(m,istate)*grad_rho_b(m,istate)
   grad_rho_a_b(istate) += grad_rho_a(m,istate)*grad_rho_b(m,istate)
  enddo
  ! convertion from (alpha,beta) formalism to (closed, open) formalism
  call rho_ab_to_rho_oc(rho_a(istate),rho_b(istate),rhoo,rhoc)
  sigmaco = 0.d0
  sigmaoo = 0.d0
  call grad_rho_ab_to_grad_rho_oc(grad_rho_a_2(istate),grad_rho_b_2(istate),grad_rho_a_b(istate),sigmaoo,sigmacc,sigmaco)
  double precision :: delta,two_dm_corr,rhoo_2
  call ec_pbe_only(mu,rhoc,rhoo,sigmacc,sigmaco,sigmaoo,e_PBE(istate))
 enddo

end

subroutine beta_from_on_top_grad_on_top_e_pbe(r, mu, Beta)  !!  A COMPLETER
 implicit none
 double precision, intent(in) :: r(3), mu
 double precision, intent (out) :: Beta(N_states)
 BEGIN_DOC
 ! routine that computes the beta parameter in WB
 END_DOC

 double precision :: c_1, pi, two_dm, two_dm_in_r_diata, e_pbe(N_states)
 integer :: i
 pi = dacos(-1.d0)
 c_1  = 3.d0 / (2.d0 * dsqrt(pi)*(1.d0 - dsqrt(2.d0)))
 call ec_pbe_compact(r,mu,e_pbe) 
 ! Somme sur les états
 do i = 1, N_states
  ! call ec_pbe_compact(r,mu,e_pbe) 
  two_dm = two_dm_in_r_diata(r,r,i)  
  Beta(i) = c_1*e_pbe(i)/two_dm
 enddo
 !print*,'Beta(1)= ', Beta(1), 'Beta(2)= ', Beta(2), 'Beta(3)= ', Beta(3)
end 
