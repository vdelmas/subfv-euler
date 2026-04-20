module global_data_module
  use precision_module
  implicit none

  integer, parameter :: mnmesh = 10 !Maximum number of sucessive meshes
  integer, parameter :: mnbc = 10 !Maximum number of boundaries
  real(kind=DOUBLE), parameter :: PI = 4.0_DOUBLE*datan(1.0_DOUBLE)
  real(kind=DOUBLE) :: gamma = 1.4_DOUBLE
  real(kind=DOUBLE) :: Cv_p = 720.19471_DOUBLE

  logical :: mesh_adapt = .false.

  !Mesh
  integer(kind=ENTIER) :: n_meshes = 1
  character(len=255) :: meshfile_path, meshfile(mnmesh)
  logical :: mesh_wasilij = .FALSE.
  logical :: use_cylinder_map = .FALSE.
  logical :: odd_even_perturb = .FALSE.
  real(kind=DOUBLE) :: odd_even_perturb_l=1e-3
  integer(kind=ENTIER) :: dim_cylinder_map = 2

  !Boundary conditions
  integer(kind=ENTIER) :: n_bc = 0
  character(len=255), dimension(mnbc) :: bc_name, bc_type
  real(kind=DOUBLE), dimension(5, mnbc) :: bc_val

  !Time
  character(len=255) :: timedisc = "euler"
  real(kind=DOUBLE) :: t_max(mnmesh), cfl
  logical :: dt_is_fixed = .false.
  real(kind=DOUBLE) :: fixed_dt=1e-2
  integer(kind=ENTIER) :: n_sol_aff = 0, n_iter_aff = 0
  logical :: local_time_step=.false.
  real(kind=DOUBLE) :: t_max_simu = 86400.
  integer(kind=ENTIER) :: n_max_iter = huge(1_ENTIER)

  !Space
  logical :: second_order = .FALSE.
  integer(kind=ENTIER) :: method = 2

  !Scheme
  character(len=255) :: scheme = ""
  logical :: lambda_l_eq_lambda_r = .FALSE.
  logical :: hybrid_mp2p = .FALSE.
  logical :: low_dissip_switch = .FALSE.
  logical :: exclude_bound_vert = .FALSE.
  real(kind=DOUBLE) :: mach_threshold = 0.3_DOUBLE
  integer(kind=ENTIER) :: bc_style = 1
  logical :: boundary_2d = .FALSE.

  !Init
  logical :: init_uniform = .FALSE.
  real(kind=DOUBLE), dimension(5) :: sol_uniform
  logical :: init_sod = .FALSE.
  real(kind=DOUBLE) :: x_sod = 0.
  real(kind=DOUBLE) :: sod_angle = 0.
  real(kind=DOUBLE), dimension(5) :: sol_left, sol_right
  logical :: init_hexa = .FALSE.
  logical :: init_ellig = .FALSE.
  logical :: init_sedov = .FALSE.
  real(kind=DOUBLE) :: r_sedov = 0.
  logical :: init_noh = .FALSE.
  logical :: init_isentropic_vortex = .FALSE.
  logical :: error_2d = .FALSE.
  real(kind=DOUBLE) :: error_2d_h = 0.025_DOUBLE
  logical :: init_gresho = .FALSE.
  logical :: init_kelvin = .FALSE.
  logical :: init_double_mach = .FALSE.
  logical :: init_grad_p = .FALSE.
  logical :: init_spherical_sod = .FALSE.
  real(kind=DOUBLE) :: r_spherical_sod = 0.
  real(kind=DOUBLE), dimension(5) :: spherical_sod_sol_left
  real(kind=DOUBLE), dimension(5) :: spherical_sod_sol_right
  logical :: init_restart = .FALSE.
  character(len=255) :: restart_file
  integer(kind=ENTIER) :: id_vtk_restart = 0
  integer(kind=ENTIER) :: restart_iter=0
  real(kind=DOUBLE) :: restart_time=0.0, restart_cpu_time=0.0
  logical :: init_double_vortex_rings = .FALSE.
  logical :: init_vortex_ring_generator = .FALSE.
  logical :: cyl_bound_wasilij = .FALSE.
  real(kind=DOUBLE) :: cyl_bound_wasilij_r1 = 0.5
  real(kind=DOUBLE) :: cyl_bound_wasilij_r2 = 1.25
  logical :: error_shear

  !Output
  logical :: write_cell_size = .FALSE.
  logical :: write_dt = .FALSE.
  logical :: write_mass = .FALSE.
  logical :: write_radial_sol = .FALSE.
  logical :: write_radial_cp = .FALSE.
  integer(kind=ENTIER) :: cp_tag = 0
  real(kind=DOUBLE), dimension(3) :: radial_center = 0.
  logical :: compute_cp = .FALSE.
  logical :: compute_coeffs = .FALSE.
  integer(kind=ENTIER) :: coeffs_surf=1
  real(kind=DOUBLE) :: pinf = 0., rhoinf = 0., vinf = 0.
  logical :: plot_solution_dat = .FALSE.
  real(kind=DOUBLE) :: xmin_dat, xmax_dat
  real(kind=DOUBLE) :: ymin_dat, ymax_dat
  real(kind=DOUBLE) :: zmin_dat, zmax_dat
  logical :: print_lambdas = .FALSE.
  integer(kind=ENTIER) :: fn_lambdas
  logical :: write_residual = .FALSE.
  integer(kind=ENTIER) :: fn_residual, n_iter_residual=100

contains
  subroutine read_input_parameters(filename)
    implicit none

    character(len=*), intent(in) :: filename

    integer(kind=ENTIER) :: funit

    namelist /INPUT_PARAM/ &
      n_meshes, meshfile_path, meshfile, mesh_wasilij, cyl_bound_wasilij, &
      cyl_bound_wasilij_r1, cyl_bound_wasilij_r2, &
      use_cylinder_map, dim_cylinder_map, &
      odd_even_perturb, odd_even_perturb_l, t_max_simu, &
      n_bc, bc_name, bc_type, bc_val, &
      timedisc, t_max, cfl, n_sol_aff, n_iter_aff, &
      second_order, method, &
      scheme, lambda_l_eq_lambda_r, exclude_bound_vert, hybrid_mp2p, &
      local_time_step, &
      dt_is_fixed, fixed_dt, &
      mach_threshold, bc_style, boundary_2d, low_dissip_switch, &
      init_uniform, sol_uniform, &
      init_sod, x_sod, sol_left, sol_right, sod_angle, &
      init_hexa, &
      init_ellig, &
      init_sedov, r_sedov, &
      init_noh, &
      init_isentropic_vortex, error_2d, error_2d_h, &
      init_gresho, &
      init_kelvin, &
      init_double_mach, &
      init_grad_p, &
      init_spherical_sod, r_spherical_sod, &
      spherical_sod_sol_left, spherical_sod_sol_right, &
      init_restart, restart_file, id_vtk_restart, &
      init_double_vortex_rings, &
      init_vortex_ring_generator, &
      write_dt, &
      write_mass, &
      write_radial_sol, write_radial_cp, cp_tag, radial_center, &
      compute_cp, pinf, rhoinf, vinf, &
      compute_coeffs, coeffs_surf, &
      write_cell_size, &
      plot_solution_dat, &
      xmin_dat, xmax_dat, &
      ymin_dat, ymax_dat, &
      zmin_dat, zmax_dat, &
      mesh_adapt, &
      error_shear, &
      n_max_iter, &
      write_residual, n_iter_residual


    open (newunit=funit, file=trim(adjustl(filename)))
    read (nml=INPUT_PARAM, unit=funit)
    close (unit=funit)
  end subroutine read_input_parameters
end module global_data_module
