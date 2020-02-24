program first_test
  implicit none
  use bitmasks ! you need to include the bitmasks_module.f90 features
  BEGIN_DOC
! We want to read the wave function and look at it
  END_DOC
  read_wf = .True. 
  touch read_wf 
  print*,'N_det = ',N_det
  integer :: i,j
  integer(bit_kind), allocatable :: det_i(:,:),det_j(:,:)
  double precision               :: hij
  double precision, allocatable :: hmat(:,:)
  double precision, allocatable :: hmat_diag(:,:)
  allocate(hmat(N_det,N_det))
  allocate(hmat_diag(N_det,N_det))


  allocate(det_i(N_int,2),det_j(N_int,2))
  do i = 1, N_det 
   print*,''
   print*,''
   print*,'i',i,psi_coef(i,1)
   !
   det_i(1,1) = psi_det(1,1,i) ! alpha component of ith determinant 
   det_i(1,2) = psi_det(1,2,i) ! beta  component of ith determinant 
   call print_det(det_i,N_int)
   do j = 1, N_det
    det_j(1,1) = psi_det(1,1,j) ! alpha component of ith determinant 
    det_j(1,2) = psi_det(1,2,j) ! beta  component of ith determinant 
    call i_H_j(det_i,det_j,N_int,hij)
    hmat(i,j) = hij
   enddo
   print*,''
  enddo

  do i = 1, N_det
   hmat_diag(i,i) = hmat(i,i)
   ! hmat_diag(i,i) = 0. 
    do j = 1, N_det
     if (j /= i) then
     hmat_diag(i,i) = hmat_diag(i,i) + hmat(i,j)*hmat(j,i)/(hmat(i,i) - hmat(j,j)) 
     end if
    enddo 
  enddo

  do i = 1, N_det
   do j = 1, N_det
    if (i /= j) then
     hmat_diag(i,j) = 0.
    end if 
   enddo
  enddo

  print*,''
  print*,'Hamiltonian matrix'
  do i = 1, N_det
   write(*,'(100(F10.5,X))')hmat(i,:)
  enddo
  print*,''
  print*,''

  print*,''
  print*,'Diagonalized Hamiltonian'
  do i = 1, N_det
   write(*,'(100(F10.5,X))')hmat_diag(i,:)
  enddo
  print*,''
  print*,''
  
  double precision :: get_two_e_integral, bielec_integral, hmono
  double precision, allocatable :: hmat_calc(:,:)
  allocate(hmat_calc(N_det,N_det))

  print*,'Hamiltonian from integral calculations'
  ! Hamiltonian from integrals calculations
  ! a b c d
  !   e f g
  !     h i
  !       j

  !------- Diagonal terms --------
  ! a
  !                                    <ij|kl>  
  bielec_integral = get_two_e_integral(1,1,1,1,mo_integrals_map)
  hmono = 2.d0 * mo_one_e_integrals(1,1)
  print*,'a = ', hmono + bielec_integral
  hmat_calc(1,1) = hmono + bielec_integral

  ! e
  bielec_integral = get_two_e_integral(2,2,2,2,mo_integrals_map)
  hmono = 2.d0 * mo_one_e_integrals(2,2)
  print*,'e = ', hmono + bielec_integral
  hmat_calc(2,2) = hmono + bielec_integral
  ! h and j 
  bielec_integral = get_two_e_integral(1,2,1,2,mo_integrals_map)
  hmono = 1.d0 * mo_one_e_integrals(2,2) + 1.d0 * mo_one_e_integrals(1,1)
  print*,'h = j =', hmono + bielec_integral
  hmat_calc(3,3) = hmono + bielec_integral
  hmat_calc(4,4) = hmono + bielec_integral

  !------ Anti-diagonal terms ------
  ! b and i
  bielec_integral = get_two_e_integral(2,2,1,1,mo_integrals_map) 
  print*,'b = i =', bielec_integral  
  hmat_calc(1,2) = bielec_integral
  hmat_calc(2,1) = bielec_integral
  hmat_calc(3,4) = bielec_integral
  hmat_calc(4,3) = bielec_integral

  ! c and d
  bielec_integral = get_two_e_integral(1,1,1,2,mo_integrals_map) 
  hmono = 1.d0 * mo_one_e_integrals(2,1)
  print*,'c = d =', hmono + bielec_integral
  hmat_calc(3,1) = hmono + bielec_integral
  hmat_calc(1,3) = hmono + bielec_integral
  hmat_calc(4,1) = hmono + bielec_integral
  hmat_calc(1,4) = hmono + bielec_integral

  ! f and g
  bielec_integral = get_two_e_integral(2,2,1,2,mo_integrals_map) 
  hmono = 1.d0 * mo_one_e_integrals(1,2)
  print*,'f = g =', hmono + bielec_integral 
  hmat_calc(3,2) = hmono + bielec_integral
  hmat_calc(2,3) = hmono + bielec_integral
  hmat_calc(4,2) = hmono + bielec_integral
  hmat_calc(2,4) = hmono + bielec_integral

  print*,''
  do i = 1, N_det
   write(*,'(100(F10.5,X))')hmat_calc(i,:)
  enddo
  print*,''
  print*,'' 

end
