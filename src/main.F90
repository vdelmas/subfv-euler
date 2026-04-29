program main
  use mpi
  use subfv_mpi_module
  use subfv_precision_module
  use global_data_module
  use subfv_mesh_module
  use subfv_mesh_geometry_module
  use subfv_mesh_reading_module
  use subfv_mesh_connectivity_module
  use io_module
  use euler_module
  use subfv_sparse_csr_linear_module
  use subfv_sparse_bcsr_linear_module

  implicit none

  integer(kind=ENTIER) :: me, num_procs, mpi_ierr

  integer(kind=ENTIER) :: iaff = 1, iter = 0, id_mesh
  real(kind=DOUBLE) :: t = 0., dt = 1., dt2
  real(kind=DOUBLE), dimension(:), allocatable :: t_sol_aff

  real(kind=DOUBLE) :: tstart, tend, tstartio, tendio, ttotio
  real(kind=DOUBLE) :: t1, t2

  type :: simulation_type
    type(mesh_type) :: mesh

    type(mpi_send_recv_type) :: mpi_send_recv

    real(kind=DOUBLE), dimension(:, :), allocatable :: sol, new_sol, sol1
    real(kind=DOUBLE), dimension(:), allocatable :: residu
    real(kind=DOUBLE), dimension(:), allocatable :: sum_lambda

    real(kind=DOUBLE), dimension(:, :), allocatable :: euler_flux_sum

    !Second order variables
    real(kind=DOUBLE), dimension(:, :, :), allocatable :: grad
    real(kind=DOUBLE), dimension(:, :), allocatable :: divF, u_nodal
    real(kind=DOUBLE), dimension(:), allocatable :: alpha_mcfl
  end type simulation_type

  type(simulation_type), dimension(:), allocatable :: simu

  integer :: fn_residual_error
  real(kind=DOUBLE), dimension(10) :: dummy_read
  logical :: res_exists

  integer(kind=ENTIER) :: ngp

  call MPI_INIT(mpi_ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, num_procs, mpi_ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, me, mpi_ierr)

  call read_input_parameters("input_data.f")


  allocate (simu(n_meshes))

  do id_mesh = 1, n_meshes
  call read_mesh_msh(simu(id_mesh)%mesh, meshfile_path, meshfile(id_mesh), &
    n_bc, bc_name, me, num_procs, simu(id_mesh)%mpi_send_recv)

  if (use_cylinder_map) call map_cylinder(simu(id_mesh)%mesh, dim_cylinder_map)
  !if (odd_even_perturb) call odd_even(simu(id_mesh)%mesh, odd_even_perturb_l)
  if (odd_even_perturb) call odd_even(simu(id_mesh)%mesh)

  call build_mesh(simu(id_mesh)%mesh, num_procs, simu(id_mesh)%mpi_send_recv, .true., boundary_2d)

  call compute_geometry_mesh(simu(id_mesh)%mesh, .true., boundary_2d)

  allocate (simu(id_mesh)%sol(5, simu(id_mesh)%mesh%n_elems))
  simu(id_mesh)%sol(:, :) = 0.0_DOUBLE
  allocate (simu(id_mesh)%sum_lambda(simu(id_mesh)%mesh%n_elems))
  simu(id_mesh)%sum_lambda(:) = 0.0_DOUBLE
  allocate (simu(id_mesh)%new_sol(5, simu(id_mesh)%mesh%n_elems))
  simu(id_mesh)%new_sol(:, :) = 0.0_DOUBLE
  allocate (simu(id_mesh)%residu(simu(id_mesh)%mesh%n_elems))
  simu(id_mesh)%residu(:) = 0.0_DOUBLE

  allocate (simu(id_mesh)%euler_flux_sum(5, simu(id_mesh)%mesh%n_elems))
  simu(id_mesh)%euler_flux_sum(:, :) = 0.0_DOUBLE

  allocate (simu(id_mesh)%u_nodal(3, simu(id_mesh)%mesh%n_vert))
  simu(id_mesh)%u_nodal = 0.0_DOUBLE

  allocate (simu(id_mesh)%alpha_mcfl(simu(id_mesh)%mesh%n_vert))
  simu(id_mesh)%alpha_mcfl = 0.0_DOUBLE

  if (timedisc == "rk2") then
    allocate (simu(id_mesh)%sol1(5, simu(id_mesh)%mesh%n_elems))
    simu(id_mesh)%sol1(:, :) = 0.0_DOUBLE
  else if (timedisc /= "euler") then
    print *, "Bad timedisc !"
    error stop
  end if

  if (second_order .or. write_cell_size ) then
    allocate (simu(id_mesh)%grad(3, 5, simu(id_mesh)%mesh%n_elems))
    simu(id_mesh)%grad(:, :, :) = 0.0_DOUBLE
    allocate (simu(id_mesh)%divF(5, simu(id_mesh)%mesh%n_elems))
    simu(id_mesh)%divF = 0.0_DOUBLE
  end if
  end do

  t = 0.0_DOUBLE
  allocate (t_sol_aff(n_sol_aff))
  call init_t_sol_aff(t, t_max(n_meshes), n_sol_aff, t_sol_aff)

  if (write_residual) then
    inquire(file="residual.dat", exist=res_exists)
    open(newunit=fn_residual, file="residual.dat")
    if (init_restart .and. res_exists) then
      do
      read(fn_residual, *, iostat = fn_residual_error) restart_iter, restart_time, &
        dummy_read(1:3), restart_cpu_time
      if ( fn_residual_error == iostat_end ) exit
      end do
      close(fn_residual)
      open(newunit=fn_residual, file="residual.dat", position="append")
    end if
  end if

  ttotio = 0.0_DOUBLE
  iter = 0
  call cpu_time(tstart)
  call MPI_ALLREDUCE(MPI_IN_PLACE, tstart, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD, mpi_ierr)

  do id_mesh = 1, n_meshes

  if (id_mesh == 1) call init_sol(simu(1)%mesh, simu(1)%sol, me, num_procs, simu(1)%mpi_send_recv)

  if (second_order) then
    if(method == 0) then
      simu(id_mesh)%grad(:, :, :) = 0.0_DOUBLE
    else if(method == 1) then
      call compute_grad_least_squares(simu(id_mesh)%mesh, simu(id_mesh)%sol, &
        simu(id_mesh)%grad)
    else if(method == 2 .or. method == 3 .or. method == 4 &
      .or. method == 5 .or. method == 6 .or. method == 7) then
      call compute_cell_grad_from_nodal_grad(simu(id_mesh)%mesh, simu(id_mesh)%sol, &
        simu(id_mesh)%grad, method)
    else
      print*,"Unknown method for second order reconstruction!"
      error stop
    end if
    if (num_procs > 1) call mpi_memory_exchange(simu(id_mesh)%mpi_send_recv, &
      simu(id_mesh)%mesh%n_elems, 3*5, simu(id_mesh)%grad)
  endif
  call check_write_sol(simu(id_mesh)%mesh, simu(id_mesh)%sol, &
    simu(id_mesh)%grad, simu(id_mesh)%residu, &
    simu(id_mesh)%alpha_mcfl, &
    t, dt, t_sol_aff, iaff, me, num_procs, simu(id_mesh)%mpi_send_recv)

  do while (t < t_max(id_mesh) .and. iter < n_max_iter)

  if (timedisc == "euler") then

    call compute_dt_a_priori(simu(id_mesh)%mesh, t, dt, t_max(id_mesh), cfl, simu(id_mesh)%sol, num_procs)

    if( second_order ) then
      if(method == 0) then
        simu(id_mesh)%grad(:, :, :) = 0.0_DOUBLE
      else if(method == 1) then
        call compute_grad_least_squares(simu(id_mesh)%mesh, simu(id_mesh)%sol, &
          simu(id_mesh)%grad)
      else if( method == 2 .or. method == 3 .or. method == 4 &
        .or. method == 5 .or. method == 6 .or. method == 7) then
        call compute_cell_grad_from_nodal_grad(simu(id_mesh)%mesh, simu(id_mesh)%sol, &
          simu(id_mesh)%grad, method)
      end if
      if (num_procs > 1) call mpi_memory_exchange(simu(id_mesh)%mpi_send_recv, &
        simu(id_mesh)%mesh%n_elems, 3*5, simu(id_mesh)%grad)
      !Compute divF
      call compute_divF(simu(id_mesh)%mesh, simu(id_mesh)%sol, &
        simu(id_mesh)%grad, simu(id_mesh)%divF)
    end if

    call compute_euler_flux_mcfl(simu(id_mesh)%mesh, simu(id_mesh)%sol, &
      simu(id_mesh)%euler_flux_sum, simu(id_mesh)%grad, t, dt, simu(id_mesh)%divF,&
      simu(id_mesh)%u_nodal, simu(id_mesh)%sum_lambda, simu(id_mesh)%alpha_mcfl)
    call update_new_sol_with_flux(simu(id_mesh)%mesh, simu(id_mesh)%new_sol, simu(id_mesh)%sol, &
      simu(id_mesh)%euler_flux_sum, &
      simu(id_mesh)%sum_lambda, simu(id_mesh)%residu, dt, num_procs)

    t = t + dt
  else if (timedisc == "rk2") then
    if( second_order ) then
      if(method == 0) then
        simu(id_mesh)%grad(:, :, :) = 0.0_DOUBLE
      else if(method == 1) then
        call compute_grad_least_squares(simu(id_mesh)%mesh, simu(id_mesh)%sol, &
          simu(id_mesh)%grad)
      else if( method == 2 .or. method == 3 .or. method == 4 &
        .or. method == 5 .or. method == 6 .or. method == 7) then
        call compute_cell_grad_from_nodal_grad(simu(id_mesh)%mesh, simu(id_mesh)%sol, &
          simu(id_mesh)%grad, method)
      end if
      if (num_procs > 1) call mpi_memory_exchange(simu(id_mesh)%mpi_send_recv, &
        simu(id_mesh)%mesh%n_elems, 3*5, simu(id_mesh)%grad)
    end if
    call compute_dt_a_priori(simu(id_mesh)%mesh, t, dt, t_max(id_mesh), cfl, simu(id_mesh)%sol, num_procs)
    call compute_euler_flux_mcfl(simu(id_mesh)%mesh, simu(id_mesh)%sol, &
      simu(id_mesh)%euler_flux_sum, simu(id_mesh)%grad, t, dt, simu(id_mesh)%divF, &
      simu(id_mesh)%u_nodal, simu(id_mesh)%sum_lambda, simu(id_mesh)%alpha_mcfl)
    dt2 = 0.5_DOUBLE*dt
    call update_new_sol_with_flux(simu(id_mesh)%mesh, simu(id_mesh)%sol1, simu(id_mesh)%sol, &
      simu(id_mesh)%euler_flux_sum, &
      simu(id_mesh)%sum_lambda, simu(id_mesh)%residu, dt2, num_procs)
    t = t + dt2
    if (num_procs > 1) call mpi_memory_exchange(simu(id_mesh)%mpi_send_recv, &
      simu(id_mesh)%mesh%n_elems, 5, simu(id_mesh)%sol1)

    if( second_order ) then
      if(method == 0) then
        simu(id_mesh)%grad(:, :, :) = 0.0_DOUBLE
      else if(method == 1) then
        call compute_grad_least_squares(simu(id_mesh)%mesh, simu(id_mesh)%sol1, &
          simu(id_mesh)%grad)
      else if( method == 2 .or. method == 3 .or. method == 4 &
        .or. method == 5 .or. method == 6 .or. method == 7) then
        call compute_cell_grad_from_nodal_grad(simu(id_mesh)%mesh, simu(id_mesh)%sol1, &
          simu(id_mesh)%grad, method)
      end if
      if (num_procs > 1) call mpi_memory_exchange(simu(id_mesh)%mpi_send_recv, &
        simu(id_mesh)%mesh%n_elems, 3*5, simu(id_mesh)%grad)
    end if

    call compute_euler_flux_mcfl(simu(id_mesh)%mesh, simu(id_mesh)%sol1, &
      simu(id_mesh)%euler_flux_sum, simu(id_mesh)%grad, t, dt, simu(id_mesh)%divF, &
      simu(id_mesh)%u_nodal, simu(id_mesh)%sum_lambda, simu(id_mesh)%alpha_mcfl)
    call update_new_sol_with_flux(simu(id_mesh)%mesh, simu(id_mesh)%new_sol, simu(id_mesh)%sol, &
      simu(id_mesh)%euler_flux_sum, &
      simu(id_mesh)%sum_lambda, simu(id_mesh)%residu, dt, num_procs)
    t = t + 0.5_DOUBLE*dt
  else
    print *, "Bad timedisc !"
    error stop
  end if

  if (write_residual .and. mod(iter, n_iter_residual) == 0) call write_dat_residual(simu(id_mesh)%mesh, &
    simu(id_mesh)%sol, simu(id_mesh)%new_sol, &
    simu(id_mesh)%euler_flux_sum, &
    iter, t, dt, tstart, fn_residual, num_procs, simu(id_mesh)%mpi_send_recv)

  simu(id_mesh)%sol = simu(id_mesh)%new_sol
  if (num_procs > 1) call mpi_memory_exchange(simu(id_mesh)%mpi_send_recv, &
    simu(id_mesh)%mesh%n_elems, 5, simu(id_mesh)%sol)

  iter = iter + 1
  call cpu_time(tstartio)
  call check_write_sol(simu(id_mesh)%mesh, simu(id_mesh)%sol, &
    simu(id_mesh)%grad, simu(id_mesh)%residu, &
    simu(id_mesh)%alpha_mcfl, &
    t, dt, t_sol_aff, iaff, me, num_procs, simu(id_mesh)%mpi_send_recv)
  call check_print_status(iter, n_iter_aff, t, dt, me)
  if (write_dt) call write_dat_dt(iter, dt)
  call cpu_time(tendio)
  ttotio = ttotio + tendio - tstartio
  end do

  !Project solution onto next mesh
  if (id_mesh < n_meshes) then
    call cpu_time(t1)
    if( num_procs > 1 ) then
      call mpi_project_sol_box(5, simu(id_mesh)%mesh, simu(id_mesh)%sol, &
        simu(id_mesh + 1)%mesh, simu(id_mesh + 1)%sol)
      if (num_procs > 1) call mpi_memory_exchange(simu(id_mesh+1)%mpi_send_recv, &
        simu(id_mesh+1)%mesh%n_elems, 5, simu(id_mesh+1)%sol)
    else
      call project_sol_box(5, simu(id_mesh)%mesh, simu(id_mesh)%sol, &
        simu(id_mesh + 1)%mesh, simu(id_mesh + 1)%sol)
    end if
    call cpu_time(t2)
    print *, ""//achar(27)//"[33m[*] Time for projection onto next mesh :"//achar(27)//"[0m", t2 - t1
  end if
  end do

  call cpu_time(tend)

  id_mesh = n_meshes

  if (write_radial_sol) call write_dat_radial_sol(simu(id_mesh)%mesh, simu(id_mesh)%sol)
  if (write_radial_cp) call write_dat_radial_cp(simu(id_mesh)%mesh, simu(id_mesh)%sol)

  if (init_isentropic_vortex) call compute_error_isentropic(simu(id_mesh)%mesh, simu(id_mesh)%sol, t)
  if (error_shear) call compute_error_shear(simu(id_mesh)%mesh, simu(id_mesh)%sol)

  print *, ""//achar(27)//"[33m[*] Time with io     "//achar(27)//"[0m", tend - tstart
  print *, ""//achar(27)//"[33m[*] Time without io  "//achar(27)//"[0m", tend - tstart - ttotio

  if (write_residual) then
    close(fn_residual)
  end if

  call MPI_FINALIZE(mpi_ierr)
end program main
