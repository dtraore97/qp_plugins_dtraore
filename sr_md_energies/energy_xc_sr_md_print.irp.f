program energy_md
implicit none
 BEGIN_DOC
 ! 
 END_DOC
 print*,'Energy_x_sr_pbe_md= ', energy_x_sr_pbe_md_copy2,' , Energy_x_sr_pbe= ', energy_x_sr_pbe
 print*,'Energy_c_sr_pbe_md= ', energy_c_sr_pbe_md_copy2,' , Energy_c_sr_pbe= ', energy_c_sr_pbe

 integer           :: istate, i
 double precision  :: mu
 double precision  :: r(3), r_norm, weight
 double precision  :: rho_a,rho_b,grad_rho_a_2,grad_rho_b_2,grad_rho_a_b
 double precision  :: grad_rho_ax, grad_rho_ay, grad_rho_az
 double precision  :: grad_rho_bx, grad_rho_by, grad_rho_bz
 double precision  :: ex_srmuPBE,dexdrho_a,dexdrho_b,dexdgrad_rho_a_2,dexdgrad_rho_b_2,dexdgrad_rho_a_b
 double precision  :: ec_srmuPBE,decdrho_a,decdrho_b,decdgrad_rho_a_2,decdgrad_rho_b_2,decdgrad_rho_a_b
 double precision  :: energy_c_sr_pbe_md(N_states), energy_x_sr_pbe_md(N_states)
 mu = 0.5d0

do i=1, n_points_final_grid
 r(1) = final_grid_points(1,i)
 r(2) = final_grid_points(2,i)
 r(3) = final_grid_points(3,i)
 r_norm = dsqrt(r(1)**2 + r(2)**2 + r(3)**2)
 weight = final_weight_at_r_vector(i)

 do istate=1, N_states
   rho_a =  one_e_dm_and_grad_alpha_in_r(4,i,istate)
   rho_b =  one_e_dm_and_grad_beta_in_r(4,i,istate)
   grad_rho_ax = one_e_dm_and_grad_alpha_in_r(1,i,istate)
   grad_rho_ay = one_e_dm_and_grad_alpha_in_r(2,i,istate)
   grad_rho_az = one_e_dm_and_grad_alpha_in_r(3,i,istate)
   grad_rho_bx = one_e_dm_and_grad_alpha_in_r(1,i,istate)
   grad_rho_by = one_e_dm_and_grad_alpha_in_r(2,i,istate)
   grad_rho_bz = one_e_dm_and_grad_alpha_in_r(3,i,istate)  
   grad_rho_a_2 = 0.d0
   grad_rho_b_2 = 0.d0
   grad_rho_a_b = 0.d0

   grad_rho_a_2 = grad_rho_ax * grad_rho_ax + grad_rho_ay * grad_rho_ay + grad_rho_az * grad_rho_az
   grad_rho_b_2 = grad_rho_bx * grad_rho_bx + grad_rho_by * grad_rho_by + grad_rho_bz * grad_rho_bz
   grad_rho_a_b = grad_rho_ax * grad_rho_bx + grad_rho_ay * grad_rho_by + grad_rho_az * grad_rho_bz
  
 call exmdsrPBE(mu,rho_a,rho_b,grad_rho_a_2,grad_rho_b_2,grad_rho_a_b,ex_srmuPBE,dexdrho_a,dexdrho_b,dexdgrad_rho_a_2,dexdgrad_rho_b_2,dexdgrad_rho_a_b)
 call ecmdsrPBE(mu,rho_a,rho_b,grad_rho_a_2,grad_rho_b_2,grad_rho_a_b,ec_srmuPBE,decdrho_a,decdrho_b,decdgrad_rho_a_2,decdgrad_rho_b_2,decdgrad_rho_a_b)
 
 energy_c_sr_pbe_md(istate) += ec_srmuPBE*weight
 energy_x_sr_pbe_md(istate) += ex_srmuPBE*weight
 
 enddo
 print*,'Ec_sr_pbe_md=', energy_c_sr_pbe_md(1)
 print*,'Ex_sr_pbe_md=', energy_x_sr_pbe_md(1) 
enddo 
end program

