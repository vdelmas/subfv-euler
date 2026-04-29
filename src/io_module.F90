module io_module
  use subfv_precision_module
  use subfv_mesh_module
  use subfv_mpi_module
  implicit none

contains

  subroutine check_print_status(iter, n_iter_aff, t, dt, me)
    implicit none

    integer(kind=ENTIER), intent(in) :: iter, n_iter_aff, me
    real(kind=DOUBLE), intent(in) :: t, dt

    if (mod(iter, n_iter_aff) == 0) then
      if( me==0 ) write (*, '(a, i10, a, g10.3, a, g10.3, a)') &
        achar(27)//"[33m [*] Iteration : ", iter, &
        " Time : ", t, " dt : ", dt, &
        ""//achar(27)//"[0m"
    end if
  end subroutine check_print_status

  subroutine check_write_sol(mesh, sol, grad, residu, &
      alpha_mcfl, t, dt, t_sol_aff, iaff, me, num_procs, mpi_send_recv)
    use global_data_module
    implicit none

    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(in) :: sol
    real(kind=DOUBLE), intent(in) :: t, dt
    integer(kind=ENTIER), intent(in) :: me, num_procs
    integer(kind=ENTIER), intent(inout) :: iaff
    real(kind=DOUBLE), dimension(n_sol_aff), intent(in) :: t_sol_aff
    real(kind=DOUBLE), dimension(mesh%n_vert), intent(in) :: alpha_mcfl
    type(mpi_send_recv_type), intent(inout) :: mpi_send_recv

    real(kind=DOUBLE), dimension(3, 5, mesh%n_elems), intent(inout) :: grad
    real(kind=DOUBLE), dimension(mesh%n_elems), intent(in) :: residu

    character(len=255) :: iaff_char, me_str

    if (iaff <= n_sol_aff) then
      if (t >= t_sol_aff(iaff) - 0.6_DOUBLE*dt) then
        if (init_restart) then
          write (iaff_char, *) iaff - 1 + id_vtk_restart
        else
          write (iaff_char, *) iaff - 1
        end if

        write (me_str, *) me
        call write_sol_vtu(mesh, trim(adjustl(me_str))//"_output_"//trim(adjustl(iaff_char))//".vtu", &
          sol, grad, residu, alpha_mcfl)
        call write_sol_meta_pvtu("output_"//trim(adjustl(iaff_char))//".pvtu", iaff_char)
        !call write_vert_coord(mesh, "vert_coord_"//trim(adjustl(iaff_char))//".dat")
        !call write_sol_simple(mesh, sol, "sol_"//trim(adjustl(iaff_char))//".dat")
        if (plot_solution_dat) call write_sol_dat(mesh, &
          trim(adjustl(me_str))//"_output_"//trim(adjustl(iaff_char))//".dat", sol)
        if (compute_coeffs) then
          call write_sol_meta_pvtu_coeffs("coeffs_"//trim(adjustl(iaff_char))//".pvtu", iaff_char)
          call write_coeffs_vtu(mesh, &
            trim(adjustl(me_str))//"_coeffs_"//trim(adjustl(iaff_char))//".vtu", sol)
        end if

        if( write_cell_size ) call write_cell_size_dat("cell_sizes_"//trim(adjustl(iaff_char))//".dat", &
          mesh, sol, grad, mpi_send_recv)
        if( me == 0 ) print *, ""//achar(27)//"[32m[+] Sol written "//achar(27)//"[0m", iaff - 1, t
        iaff = iaff + 1
      end if
    end if
  end subroutine check_write_sol

  subroutine write_vert_coord(mesh, filename)
    implicit none

    type(mesh_type), intent(in) :: mesh
    character(len=*), intent(in) :: filename


    integer(kind=ENTIER) :: fn, i

    open(newunit=fn, file=filename)
    do i=1, mesh%n_vert
      write(fn, *) mesh%vert(i)%coord
    end do
    close(fn)

  end subroutine write_vert_coord

  subroutine write_sol_simple(mesh, sol, filename)
    use euler_module, only: conserv_to_primit
    implicit none

    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(in) :: sol
    character(len=*), intent(in) :: filename


    integer(kind=ENTIER) :: fn, i

    open(newunit=fn, file=filename)
    do i=1, mesh%n_elems
      write(fn, *) sol(:, i), conserv_to_primit(sol(:, i))
    end do
    close(fn)

  end subroutine write_sol_simple

  subroutine write_sol_vtu(mesh, filename, &
      sol, grad, residu, alpha_mcfl)
    use mpi
    use global_data_module
    use euler_module, only: sol_isentropic_vortex, conserv_to_primit, compute_nodal_grad
    use subfv_linear_solver_module
    implicit none

    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(in) :: sol
    real(kind=DOUBLE), dimension(:, :, :), intent(in) :: grad
    real(kind=DOUBLE), dimension(mesh%n_elems), intent(in) :: residu
    real(kind=DOUBLE), dimension(mesh%n_vert), intent(in) :: alpha_mcfl

    character(len=*), intent(in) :: filename


    integer(kind=ENTIER) :: me, num_procs, mpi_ierr
    integer(kind=ENTIER) :: n_interior_elems, n_interior_vert
    integer(kind=ENTIER) :: fn, i, j, s, id_face, size_tot, k
    real(kind=DOUBLE) :: dmin, dmax

    real(kind=DOUBLE), dimension(5) :: w
    real(kind=DOUBLE), dimension(3, 5) :: nodal_grad
    integer(kind=ENTIER), dimension(:), allocatable :: id_vert_no_ghost
    integer(kind=ENTIER), dimension(:), allocatable :: local_id_vert_no_ghost

    call MPI_COMM_SIZE(MPI_COMM_WORLD, num_procs, mpi_ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, me, mpi_ierr)

    allocate(id_vert_no_ghost(mesh%n_vert))
    id_vert_no_ghost = 0
    n_interior_vert = 0
    do i=1, mesh%n_vert
      if ( .not. mesh%vert(i)%is_ghost ) then
        n_interior_vert = n_interior_vert + 1
        id_vert_no_ghost(i) = n_interior_vert
      end if
    end do

    n_interior_elems = 0
    do i=1, mesh%n_elems
      if ( .not. mesh%elem(i)%is_ghost ) n_interior_elems = n_interior_elems + 1
    end do

    open(newunit=fn, file=trim(adjustl(filename)))

    write(fn, *) "<VTKFile type='UnstructuredGrid' version='1.0' &
      &byte_order='LittleEndian' header_type='UInt64'>"
    write(fn, *) "<UnstructuredGrid>\n<Piece NumberOfPoints='", n_interior_vert, &
      &"' NumberOfCells='", n_interior_elems, "'>"

    dmin = 1e100_DOUBLE
    dmax = -1e100_DOUBLE
    do i=1, mesh%n_vert
      if( .not. mesh%vert(i)%is_ghost ) then
        do j=1, 3
          if( mesh%vert(i)%coord(j) > dmax ) dmax = mesh%vert(i)%coord(j)
          if( mesh%vert(i)%coord(j) < dmin ) dmin = mesh%vert(i)%coord(j)
        end do
      end if
    end do

    write(fn, *) "<Points>"
    write(fn, *) " <DataArray type='Float64' Name='Points' NumberOfComponents='3' &
      & format='ascii' RangeMin='", dmin, "' RangeMax='", dmax, "'>"

    do i=1, mesh%n_vert
      if( .not. mesh%vert(i)%is_ghost ) then
        write(fn, *) mesh%vert(i)%coord(:)
      end if
    end do

    write(fn, *) "</DataArray>"
    write(fn, *) "</Points>"
    write(fn, *) "<Cells>"
    write(fn, *) "<DataArray type='Int64' Name='connectivity' format='ascii' &
      &RangeMin='0' RangeMax='", mesh%n_vert-1, "'>"
    do i=1, mesh%n_elems
      if( .not. mesh%elem(i)%is_ghost ) then
        write(fn, *) id_vert_no_ghost(mesh%elem(i)%vert) - 1
      end if
    end do
    write(fn, *) "</DataArray>"

    s = 0
    do i=1, mesh%n_elems
      if( .not. mesh%elem(i)%is_ghost ) then
        s = s + mesh%elem(i)%n_vert
      end if
    end do

    write(fn, *) "<DataArray type='Int64' Name='offsets' format='ascii' &
      &RangeMin='", 0, "' RangeMax='", s, "'>"

    s = 0
    do i=1, mesh%n_elems
      if( .not. mesh%elem(i)%is_ghost ) then
        s = s + mesh%elem(i)%n_vert
        write(fn, *) s
      end if
    end do

    write(fn, *) "</DataArray>"

    write(fn, *) "<DataArray type='UInt8' Name='types' format='ascii' &
      &RangeMin='42' RangeMax='42'>"

    do i=1, mesh%n_elems
      if( .not. mesh%elem(i)%is_ghost ) then
        write(fn, *) 42 !Type for polyhedra
      end if
    end do

    write(fn, *) "</DataArray>"
    write(fn, *) "<DataArray type='Int64' Name='faces' format='ascii' &
      &RangeMin='0' RangeMax='", mesh%n_vert-1, "'>"

    do i=1, mesh%n_elems
      if( .not. mesh%elem(i)%is_ghost ) then
        write(fn, *) mesh%elem(i)%n_faces
        do j=1, mesh%elem(i)%n_faces
          id_face = mesh%elem(i)%face(j)
          allocate(local_id_vert_no_ghost(mesh%face(id_face)%n_vert))
          do k=1, mesh%face(id_face)%n_vert
            local_id_vert_no_ghost(k) = id_vert_no_ghost(mesh%face(id_face)%vert(k))
          end do
          !write(fn, *) mesh%face(id_face)%n_vert, mesh%face(id_face)%vert(:) - 1
          write(fn, *) mesh%face(id_face)%n_vert, local_id_vert_no_ghost - 1
          deallocate(local_id_vert_no_ghost)
        end do
      end if
    end do

    write(fn, *) "</DataArray>"

    s = 0
    do i=1, mesh%n_elems
      if( .not. mesh%elem(i)%is_ghost ) then
        do j=1, mesh%elem(i)%n_faces
          id_face = mesh%elem(i)%face(j)
          s = s + mesh%face(id_face)%n_vert + 1
        end do
        s = s + 1
      end if
    end do

    write(fn, *) "<DataArray type='Int64' Name='faceoffsets' format='ascii' &
      &RangeMin='", 0, "' RangeMax='", s, "'>"

    s = 0
    do i=1, mesh%n_elems
      if( .not. mesh%elem(i)%is_ghost ) then
        do j=1, mesh%elem(i)%n_faces
          id_face = mesh%elem(i)%face(j)
          s = s + mesh%face(id_face)%n_vert + 1
        end do
        s = s + 1
        write(fn, *) s
      end if
    end do

    write(fn, *) "</DataArray>"
    write(fn, *) "</Cells>"

    !Point data
    write(fn, *) "<PointData>"

    write(fn, *) "<DataArray type='Float64' Name='Nodal_Grad_Density' format='ascii' NumberOfComponents='3'>"
    do i=1, mesh%n_vert
      if( .not. mesh%vert(i)%is_ghost ) then
        call compute_nodal_grad(mesh, i, sol, nodal_grad)
        write (fn, *) nodal_grad(:, 1)
      end if
    end do
    write(fn, *) "</DataArray>"

    write(fn, *) "<DataArray type='Float64' Name='Nodal_Grad_Velocity_X' format='ascii' NumberOfComponents='3'>"
    do i=1, mesh%n_vert
      if( .not. mesh%vert(i)%is_ghost ) then
        call compute_nodal_grad(mesh, i, sol, nodal_grad)
        write (fn, *) nodal_grad(:, 2)
      end if
    end do
    write(fn, *) "</DataArray>"

    write(fn, *) "<DataArray type='Float64' Name='Nodal_Grad_Velocity_Y' format='ascii' NumberOfComponents='3'>"
    do i=1, mesh%n_vert
      if( .not. mesh%vert(i)%is_ghost ) then
        call compute_nodal_grad(mesh, i, sol, nodal_grad)
        write (fn, *) nodal_grad(:, 3)
      end if
    end do
    write(fn, *) "</DataArray>"

    write(fn, *) "<DataArray type='Float64' Name='Nodal_Grad_Velocity_Z' format='ascii' NumberOfComponents='3'>"
    do i=1, mesh%n_vert
      if( .not. mesh%vert(i)%is_ghost ) then
        call compute_nodal_grad(mesh, i, sol, nodal_grad)
        write (fn, *) nodal_grad(:, 4)
      end if
    end do
    write(fn, *) "</DataArray>"

    write(fn, *) "<DataArray type='Float64' Name='Nodal_Grad_Pressure' format='ascii' NumberOfComponents='3'>"
    do i=1, mesh%n_vert
      if( .not. mesh%vert(i)%is_ghost ) then
        call compute_nodal_grad(mesh, i, sol, nodal_grad)
        write (fn, *) nodal_grad(:, 5)
      end if
    end do
    write(fn, *) "</DataArray>"

    write(fn, *) "<DataArray type='Float64' Name='Alpha_mcfl' format='ascii' NumberOfComponents='1'>"
    do i=1, mesh%n_vert
      if( .not. mesh%vert(i)%is_ghost ) then
        write (fn, *) alpha_mcfl(i)
      end if
    end do
    write(fn, *) "</DataArray>"

    write(fn, *) "</PointData>"

    !Cell data
    write(fn, *) "<CellData>"

    write(fn, *) "<DataArray type='Float64' Name='Centroid' format='ascii' NumberOfComponents='3'>"
    do i=1, mesh%n_elems
      if( .not. mesh%elem(i)%is_ghost ) then
        write(fn, *) mesh%elem(i)%coord
      end if
    end do
    write(fn, *) "</DataArray>"

    write(fn, *) "<DataArray type='Int32' Name='Tag' format='ascii' NumberOfComponents='1'>"
    do i=1, mesh%n_elems
      if( .not. mesh%elem(i)%is_ghost ) then
        write (fn, *) mesh%elem(i)%tag
      end if
    end do
    write(fn, *) "</DataArray>"

    write(fn, *) "<DataArray type='Int32' Name='MPI_Color' format='ascii' NumberOfComponents='1'>"
    do i=1, mesh%n_elems
      if( .not. mesh%elem(i)%is_ghost ) then
        write(fn, *) me
      end if
    end do
    write(fn, *) "</DataArray>"

    write(fn, *) "<DataArray type='Float64' Name='Density' format='ascii' NumberOfComponents='1'>"
    do i=1, mesh%n_elems
      if( .not. mesh%elem(i)%is_ghost ) then
        write(fn, *) sol(1, i)
      end if
    end do
    write(fn, *) "</DataArray>"

    write(fn, *) "<DataArray type='Float64' Name='Momentum' format='ascii' NumberOfComponents='3'>"
    do i=1, mesh%n_elems
      if( .not. mesh%elem(i)%is_ghost ) then
        write(fn, *) sol(2:4, i)
      end if
    end do
    write(fn, *) "</DataArray>"

    write(fn, *) "<DataArray type='Float64' Name='Energy' format='ascii' NumberOfComponents='1'>"
    do i=1, mesh%n_elems
      if( .not. mesh%elem(i)%is_ghost ) then
        write(fn, *) sol(5, i)/sol(1, i)
      end if
    end do
    write(fn, *) "</DataArray>"

    write(fn, *) "<DataArray type='Float64' Name='Velocity' format='ascii' NumberOfComponents='3'>"
    do i=1, mesh%n_elems
      if( .not. mesh%elem(i)%is_ghost ) then
        write(fn, *) sol(2:4, i)/sol(1, i)
      end if
    end do
    write(fn, *) "</DataArray>"

    write(fn, *) "<DataArray type='Float64' Name='Pressure' format='ascii' NumberOfComponents='1'>"
    do i=1, mesh%n_elems
      if( .not. mesh%elem(i)%is_ghost ) then
        w = conserv_to_primit(sol(:, i))
        write(fn, *) w(5)
      end if
    end do
    write(fn, *) "</DataArray>"

    write(fn, *) "<DataArray type='Float64' Name='Internal_energy' format='ascii' NumberOfComponents='1'>"
    do i=1, mesh%n_elems
      if( .not. mesh%elem(i)%is_ghost ) then
        write (fn, *) sol(5, i)/sol(1, i) - 0.5_DOUBLE*norm2(sol(2:4, i)/sol(1, i))**2
      end if
    end do
    write(fn, *) "</DataArray>"

    write(fn, *) "<DataArray type='Float64' Name='H' format='ascii' NumberOfComponents='1'>"
    do i=1, mesh%n_elems
      if( .not. mesh%elem(i)%is_ghost ) then
        w = conserv_to_primit(sol(:, i))
        write (fn, *) sol(5, i)/sol(1, i) + w(5)/sol(1, i)
      end if
    end do
    write(fn, *) "</DataArray>"

    write(fn, *) "<DataArray type='Float64' Name='Mach' format='ascii' NumberOfComponents='1'>"
    do i=1, mesh%n_elems
      if( .not. mesh%elem(i)%is_ghost ) then
        w = conserv_to_primit(sol(:, i))
        write (fn, *) norm2(sol(2:4, i)/sol(1, i))/(sqrt(gamma*w(5)/w(1)))
      end if
    end do
    write(fn, *) "</DataArray>"

    write(fn, *) "<DataArray type='Float64' Name='Residu' format='ascii' NumberOfComponents='1'>"
    do i=1, mesh%n_elems
      if( .not. mesh%elem(i)%is_ghost ) then
        write (fn, *) residu(i)
      end if
    end do
    write(fn, *) "</DataArray>"

    write(fn, *) "<DataArray type='Float64' Name='Temperature' format='ascii' NumberOfComponents='1'>"
    do i=1, mesh%n_elems
      if( .not. mesh%elem(i)%is_ghost ) then
        write (fn, *) (sol(5, i)/sol(1, i) - 0.5_DOUBLE*norm2(sol(2:4, i)/sol(1, i))**2)/Cv_p
      end if
    end do
    write(fn, *) "</DataArray>"

    write(fn, *) "<DataArray type='Float64' Name='r' format='ascii' NumberOfComponents='1'>"
    do i=1, mesh%n_elems
      if( .not. mesh%elem(i)%is_ghost ) then
        write (fn, *) norm2(mesh%elem(i)%coord)
      end if
    end do
    write(fn, *) "</DataArray>"

    write(fn, *) "<DataArray type='Float64' Name='rxy' format='ascii' NumberOfComponents='1'>"
    do i=1, mesh%n_elems
      if( .not. mesh%elem(i)%is_ghost ) then
        write (fn, *) norm2(mesh%elem(i)%coord(1:2))
      end if
    end do
    write(fn, *) "</DataArray>"

    write(fn, *) "<DataArray type='Float64' Name='Theta3' format='ascii' NumberOfComponents='1'>"
    do i=1, mesh%n_elems
      if( .not. mesh%elem(i)%is_ghost ) then
        write (fn, *) atan2(-mesh%elem(i)%coord(2), -mesh%elem(i)%coord(1))
      end if
    end do
    write(fn, *) "</DataArray>"

    if( second_order ) then
      write(fn, *) "<DataArray type='Float64' Name='Density_Grad' format='ascii' NumberOfComponents='3'>"
      do i=1, mesh%n_elems
        if( .not. mesh%elem(i)%is_ghost ) then
          write (fn, *) grad(:, 1, i)
        end if
      end do
      write(fn, *) "</DataArray>"

      write(fn, *) "<DataArray type='Float64' Name='Velocity_X_Grad' format='ascii' NumberOfComponents='3'>"
      do i=1, mesh%n_elems
        if( .not. mesh%elem(i)%is_ghost ) then
          write (fn, *) grad(:, 2, i)
        end if
      end do
      write(fn, *) "</DataArray>"

      write(fn, *) "<DataArray type='Float64' Name='Velocity_Y_Grad' format='ascii' NumberOfComponents='3'>"
      do i=1, mesh%n_elems
        if( .not. mesh%elem(i)%is_ghost ) then
          write (fn, *) grad(:, 3, i)
        end if
      end do
      write(fn, *) "</DataArray>"

      write(fn, *) "<DataArray type='Float64' Name='Velocity_Z_Grad' format='ascii' NumberOfComponents='3'>"
      do i=1, mesh%n_elems
        if( .not. mesh%elem(i)%is_ghost ) then
          write (fn, *) grad(:, 4, i)
        end if
      end do
      write(fn, *) "</DataArray>"

      write(fn, *) "<DataArray type='Float64' Name='Pressure_Grad' format='ascii' NumberOfComponents='3'>"
      do i=1, mesh%n_elems
        if( .not. mesh%elem(i)%is_ghost ) then
          write (fn, *) grad(:, 5, i)
        end if
      end do
      write(fn, *) "</DataArray>"
    end if

    write(fn, *) "</CellData>"

    write(fn, *) "</Piece>"
    write(fn, *) "</UnstructuredGrid>"
    write(fn, *) "</VTKFile>"

    close(fn)
  end subroutine write_sol_vtu

  subroutine write_sol_dat(mesh, filename, sol)
    use global_data_module
    use euler_module, only: sol_isentropic_vortex, conserv_to_primit
    implicit none

    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(in) :: sol

    character(len=*), intent(in) :: filename

    integer(kind=ENTIER) :: i, vtkout
    real(kind=DOUBLE), dimension(5) :: w

    open (newunit=vtkout, file=trim(filename), status="unknown")
    do i = 1, mesh%n_elems
      if (.not. mesh%elem(i)%is_ghost) then
        if ((mesh%elem(i)%coord(1) < xmax_dat .and. &
          mesh%elem(i)%coord(1) > xmin_dat) .and. &
          (mesh%elem(i)%coord(2) < ymax_dat .and. &
          mesh%elem(i)%coord(2) > ymin_dat) .and. &
          (mesh%elem(i)%coord(3) < zmax_dat .and. &
          mesh%elem(i)%coord(3) > zmin_dat)) then

          w = conserv_to_primit(sol(:, i))
          write (vtkout, *) mesh%elem(i)%coord, sol(:, i), w(:)
        end if
      end if
    end do
    close (vtkout)
  end subroutine write_sol_dat

  subroutine write_sol_vtk(mesh, filename, t, &
      sol, grad, residu, u_nodal)
    use global_data_module
    use euler_module, only: sol_isentropic_vortex, conserv_to_primit, compute_nodal_grad
    use subfv_linear_solver_module
    implicit none

    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), intent(in) :: t
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(in) :: sol
    real(kind=DOUBLE), dimension(:, :, :), intent(in) :: grad
    real(kind=DOUBLE), dimension(mesh%n_elems), intent(in) :: residu
    real(kind=DOUBLE), dimension(3, mesh%n_vert), intent(in) :: u_nodal

    character(len=*), intent(in) :: filename

    integer(kind=ENTIER) :: i, j, k, size_tot, vtkout
    integer(kind=ENTIER) :: id_vert, id_face, id_sub_face, id_sub_elem, id_elem
    real(kind=DOUBLE) :: cpmax, r
    real(kind=DOUBLE), dimension(5) :: w
    real(kind=DOUBLE), dimension(3) :: vort, norm
    real(kind=DOUBLE), dimension(3, 5) :: nodal_grad

    open (newunit=vtkout, file=trim(filename), status="unknown")

    write (vtkout, '(A)') '# vtk DataFile Version 3.0'
    write (vtkout, '(A)') filename
    write (vtkout, '(A)') 'ASCII'
    write (vtkout, '(A)') 'DATASET UNSTRUCTURED_GRID'

    write (vtkout, *) 'POINTS ', mesh%n_vert, ' double'
    do i = 1, mesh%n_vert
      write (vtkout, *) mesh%vert(i)%coord(1), mesh%vert(i)%coord(2), &
        mesh%vert(i)%coord(3)
    end do

    size_tot = 0
    do i = 1, mesh%n_elems
      if (.not. mesh%elem(i)%is_ghost) then
        size_tot = size_tot + mesh%elem(i)%n_vert + 1
      end if
    end do

    write (vtkout, *) 'CELLS ', mesh%n_interior_elems, size_tot
    do i = 1, mesh%n_elems
      if (.not. mesh%elem(i)%is_ghost) then
        write (vtkout, '(i12)', advance='no') mesh%elem(i)%n_vert
        do j = 1, mesh%elem(i)%n_vert - 1
          write (vtkout, '(i12)', advance='no') mesh%elem(i)%vert(j) - 1
        end do
        write (vtkout, '(i12)') mesh%elem(i)%vert(mesh%elem(i)%n_vert) - 1
      end if
    end do

    write (vtkout, '(a,i11)') 'CELL_TYPES ', mesh%n_interior_elems
    do i = 1, mesh%n_elems
      if (.not. mesh%elem(i)%is_ghost) then
        if (mesh%elem(i)%elem_kind == 4) then
          write (vtkout, '(i12)') 10
        else if (mesh%elem(i)%elem_kind == 5) then
          write (vtkout, '(i12)') 12
        else if (mesh%elem(i)%elem_kind == 6) then
          write (vtkout, '(i12)') 13
        else if (mesh%elem(i)%elem_kind == 7) then
          write (vtkout, '(i12)') 14
        end if
      end if
    end do

    write (vtkout, '(a,i11)') 'CELL_DATA ', mesh%n_interior_elems
    write (vtkout, *) 'SCALARS Centroid double', 3
    write (vtkout, *) 'LOOKUP_TABLE default'
    do i = 1, mesh%n_elems
      if (.not. mesh%elem(i)%is_ghost) then
        write (vtkout, *) mesh%elem(i)%coord
      end if
    end do

    write (vtkout, *) 'SCALARS r double', 1
    write (vtkout, *) 'LOOKUP_TABLE default'
    do i = 1, mesh%n_elems
      if (.not. mesh%elem(i)%is_ghost) then
        write (vtkout, *) sqrt(mesh%elem(i)%coord(1)**2 + mesh%elem(i)%coord(2)**2)
      end if
    end do

    write (vtkout, *) 'SCALARS Theta1 double', 1
    write (vtkout, *) 'LOOKUP_TABLE default'
    do i = 1, mesh%n_elems
      if (.not. mesh%elem(i)%is_ghost) then
        write (vtkout, *) atan2(mesh%elem(i)%coord(2), mesh%elem(i)%coord(1))
      end if
    end do

    write (vtkout, *) 'SCALARS Theta2 double', 1
    write (vtkout, *) 'LOOKUP_TABLE default'
    do i = 1, mesh%n_elems
      if (.not. mesh%elem(i)%is_ghost) then
        write (vtkout, *) atan2(-mesh%elem(i)%coord(2), mesh%elem(i)%coord(1))
      end if
    end do

    write (vtkout, *) 'SCALARS Theta3 double', 1
    write (vtkout, *) 'LOOKUP_TABLE default'
    do i = 1, mesh%n_elems
      if (.not. mesh%elem(i)%is_ghost) then
        write (vtkout, *) atan2(-mesh%elem(i)%coord(2), -mesh%elem(i)%coord(1))
      end if
    end do

    write (vtkout, *) 'SCALARS Theta4 double', 1
    write (vtkout, *) 'LOOKUP_TABLE default'
    do i = 1, mesh%n_elems
      if (.not. mesh%elem(i)%is_ghost) then
        write (vtkout, *) atan2(mesh%elem(i)%coord(2), -mesh%elem(i)%coord(1))
      end if
    end do

    write (vtkout, *) 'SCALARS Density double', 1
    write (vtkout, *) 'LOOKUP_TABLE default'
    do i = 1, mesh%n_elems
      if (.not. mesh%elem(i)%is_ghost) then
        write (vtkout, *) sol(1, i)
      end if
    end do

    write (vtkout, *) 'VECTORS Momentum double'
    do i = 1, mesh%n_elems
      if (.not. mesh%elem(i)%is_ghost) then
        write (vtkout, *) sol(2, i), sol(3, i), sol(4, i)
      end if
    end do

    write (vtkout, *) 'SCALARS Energy double', 1
    write (vtkout, *) 'LOOKUP_TABLE default'
    do i = 1, mesh%n_elems
      if (.not. mesh%elem(i)%is_ghost) then
        write (vtkout, *) sol(5, i)
      end if
    end do

    write (vtkout, *) 'VECTORS Velocity double'
    do i = 1, mesh%n_elems
      if (.not. mesh%elem(i)%is_ghost) then
        write (vtkout, *) sol(2, i)/sol(1, i), &
          sol(3, i)/sol(1, i), sol(4, i)/sol(1, i)
      end if
    end do

    write (vtkout, *) 'SCALARS Internal_energy double', 1
    write (vtkout, *) 'LOOKUP_TABLE default'
    do i = 1, mesh%n_elems
      if (.not. mesh%elem(i)%is_ghost) then
        write (vtkout, *) sol(5, i)/sol(1, i) - 0.5_DOUBLE*norm2(sol(2:4, i)/sol(1, i))**2
      end if
    end do

    write (vtkout, *) 'SCALARS Pressure double', 1
    write (vtkout, *) 'LOOKUP_TABLE default'
    do i = 1, mesh%n_elems
      if (.not. mesh%elem(i)%is_ghost) then
        w = conserv_to_primit(sol(:, i))
        write (vtkout, *) w(5)
      end if
    end do

    write (vtkout, *) 'SCALARS Mach double', 1
    write (vtkout, *) 'LOOKUP_TABLE default'
    do i = 1, mesh%n_elems
      if (.not. mesh%elem(i)%is_ghost) then
        w = conserv_to_primit(sol(:, i))
        write (vtkout, *) norm2(w(2:4))/sqrt(gamma*w(5)/w(1))
      end if
    end do

    write (vtkout, *) 'SCALARS H double', 1
    write (vtkout, *) 'LOOKUP_TABLE default'
    do i = 1, mesh%n_elems
      if (.not. mesh%elem(i)%is_ghost) then
        w = conserv_to_primit(sol(:, i))
        write (vtkout, *) (gamma*w(5)/w(1))/(gamma - 1) + 0.5_DOUBLE*norm2(w(2:4))**2
      end if
    end do

    write (vtkout, *) 'SCALARS Residu double', 1
    write (vtkout, *) 'LOOKUP_TABLE default'
    do i = 1, mesh%n_elems
      if (.not. mesh%elem(i)%is_ghost) then
        write (vtkout, *) residu(i)
      end if
    end do

    if (init_noh) then
      write (vtkout, *) 'SCALARS Density_exact_noh double', 1
      write (vtkout, *) 'LOOKUP_TABLE default'
      do i = 1, mesh%n_elems
        if (.not. mesh%elem(i)%is_ghost) then
          r = sqrt(mesh%elem(i)%coord(1)**2 + mesh%elem(i)%coord(2)**2)
          if (r < t/3.0_DOUBLE) then
            write (vtkout, *) ((gamma + 1)/(gamma - 1))**2
          else
            write (vtkout, *) (1.0_DOUBLE + t/r)
          end if
        end if
      end do

      write (vtkout, *) 'SCALARS Pressure_exact_noh double', 1
      write (vtkout, *) 'LOOKUP_TABLE default'
      do i = 1, mesh%n_elems
        if (.not. mesh%elem(i)%is_ghost) then
          r = sqrt(mesh%elem(i)%coord(1)**2 + mesh%elem(i)%coord(2)**2)
          if (r < t/3.0_DOUBLE) then
            write (vtkout, *) 16.0_DOUBLE/3.0_DOUBLE
          else
            write (vtkout, *) 1e-6
          end if
        end if
      end do

      write (vtkout, *) 'SCALARS Internal_energy_exact_noh double', 1
      write (vtkout, *) 'LOOKUP_TABLE default'
      do i = 1, mesh%n_elems
        if (.not. mesh%elem(i)%is_ghost) then
          r = sqrt(mesh%elem(i)%coord(1)**2 + mesh%elem(i)%coord(2)**2)
          if (r < t/3.0_DOUBLE) then
            write (vtkout, *) 1.0_DOUBLE/2.0_DOUBLE
          else
            write (vtkout, *) 1e-6
          end if
        end if
      end do

      write (vtkout, *) 'VECTORS Velocity_Mexact_noh double'
      do i = 1, mesh%n_elems
        if (.not. mesh%elem(i)%is_ghost) then
          r = sqrt(mesh%elem(i)%coord(1)**2 + mesh%elem(i)%coord(2)**2)
          if (r < t/3.0_DOUBLE) then
            write (vtkout, *) 0.0_DOUBLE, 0.0_DOUBLE, 0.0_DOUBLE
          else
            write (vtkout, *) - mesh%elem(i)%coord(1)/r, -mesh%elem(i)%coord(2)/r, 0.0_DOUBLE
          end if
        end if
      end do
    end if

    if (init_isentropic_vortex) then
      call compute_error_isentropic(mesh, sol, t)
      write (vtkout, *) 'SCALARS Density_exact_isentropic double', 1
      write (vtkout, *) 'LOOKUP_TABLE default'
      do i = 1, mesh%n_elems
        if (.not. mesh%elem(i)%is_ghost) then
          call sol_isentropic_vortex(mesh%elem(i)%coord, w, t)
          write (vtkout, *) w(1)
        end if
      end do
    end if

    if (compute_cp) then
      write (vtkout, *) 'SCALARS Cp double', 1
      write (vtkout, *) 'LOOKUP_TABLE default'
      do i = 1, mesh%n_elems
        if (.not. mesh%elem(i)%is_ghost) then
          w = conserv_to_primit(sol(:, i))
          write (vtkout, *) (w(5) - pinf)/(0.5_DOUBLE*rhoinf*vinf**2)
        end if
      end do

      cpmax = 2.0_DOUBLE/(gamma*vinf**2)*( &
        ((gamma + 1)**2*vinf**2/(4*gamma*vinf**2 - 2*(gamma - 1)))**(gamma/(gamma - 1)) &
        *(1 - gamma + 2*gamma*vinf**2)/(gamma + 1) - 1)

      write (vtkout, *) 'SCALARS Modified_Newtonian_Law double', 4
      write (vtkout, *) 'LOOKUP_TABLE default'
      do i = 1, mesh%n_elems
        if (.not. mesh%elem(i)%is_ghost) then
          write (vtkout, *) cpmax*sin(pi/2.0_DOUBLE - atan2(mesh%elem(i)%coord(2), mesh%elem(i)%coord(1)))**2, &
            cpmax*sin(pi/2.0_DOUBLE - atan2(-mesh%elem(i)%coord(2), mesh%elem(i)%coord(1)))**2, &
            cpmax*sin(pi/2.0_DOUBLE - atan2(-mesh%elem(i)%coord(2), -mesh%elem(i)%coord(1)))**2, &
            cpmax*sin(pi/2.0_DOUBLE - atan2(mesh%elem(i)%coord(2), -mesh%elem(i)%coord(1)))**2
        end if
      end do

    end if

    write (vtkout, *) 'SCALARS Tag int', 1
    write (vtkout, *) 'LOOKUP_TABLE default'
    do i = 1, mesh%n_elems
      if (.not. mesh%elem(i)%is_ghost) then
        write (vtkout, *) mesh%elem(i)%tag
      end if
    end do

    if (second_order) then
      write (vtkout, *) 'VECTORS Density_grad double'
      do i = 1, mesh%n_elems
        if (.not. mesh%elem(i)%is_ghost) then
          write (vtkout, *) grad(1, 1, i), grad(2, 1, i), grad(3, 1, i)
        end if
      end do

      write (vtkout, *) 'VECTORS Pressure_grad double'
      do i = 1, mesh%n_elems
        if (.not. mesh%elem(i)%is_ghost) then
          write (vtkout, *) grad(1, 5, i), grad(2, 5, i), grad(3, 5, i)
        end if
      end do
    end if

    write (vtkout, '(a,i11)') 'POINT_DATA ', mesh%n_vert

    write (vtkout, *) 'SCALARS is_bound double', 1
    write (vtkout, *) 'LOOKUP_TABLE default'
    do i = 1, mesh%n_vert
      if( mesh%vert(i)%is_bound ) then
        write (vtkout, *) 1.0_DOUBLE
      else
        write (vtkout, *) 0.0_DOUBLE
      end if
    end do

    write (vtkout, *) 'VECTORS u_nodal double'
    do i = 1, mesh%n_vert
      write (vtkout, *) u_nodal(1, i), u_nodal(2, i), u_nodal(3, i)
    end do

    write (vtkout, *) 'VECTORS Nodal_Grad_Density double'
    do i = 1, mesh%n_vert
      call compute_nodal_grad(mesh, i, sol, nodal_grad)
      write (vtkout, *) nodal_grad(:, 1)
    end do

    write (vtkout, *) 'VECTORS Nodal_Grad_Pressure double'
    do i = 1, mesh%n_vert
      call compute_nodal_grad(mesh, i, sol, nodal_grad)
      write (vtkout, *) nodal_grad(:, 5)
    end do

    write (vtkout, *) 'VECTORS Nodal_Grad_Velocity_X double'
    do i = 1, mesh%n_vert
      call compute_nodal_grad(mesh, i, sol, nodal_grad)
      write (vtkout, *) nodal_grad(:, 2)
    end do

    write (vtkout, *) 'VECTORS Nodal_Grad_Velocity_Y double'
    do i = 1, mesh%n_vert
      call compute_nodal_grad(mesh, i, sol, nodal_grad)
      write (vtkout, *) nodal_grad(:, 3)
    end do

    write (vtkout, *) 'VECTORS Nodal_Grad_Velocity_Z double'
    do i = 1, mesh%n_vert
      call compute_nodal_grad(mesh, i, sol, nodal_grad)
      write (vtkout, *) nodal_grad(:, 4)
    end do

    write (vtkout, *) 'VECTORS Vorticity_operator double'
    do id_vert = 1, mesh%n_vert
      vort = 0.0_DOUBLE
      do j = 1, mesh%vert(id_vert)%n_sub_elems_neigh
        id_sub_elem = mesh%vert(id_vert)%sub_elem_neigh(j)
        id_elem = mesh%sub_elem(id_sub_elem)%mesh_elem
        w = conserv_to_primit(sol(:, id_elem))
        do k = 1, mesh%sub_elem(id_sub_elem)%n_sub_faces
          id_sub_face = mesh%sub_elem(id_sub_elem)%sub_face(k)
          id_face = mesh%sub_face(id_sub_face)%mesh_face
          if(mesh%sub_face(id_sub_face)%left_elem_neigh == id_elem) then
            norm = mesh%face(id_face)%norm
          else
            norm = -mesh%face(id_face)%norm
          end if
          if(boundary_2d) norm(3) = 0.0_DOUBLE
          vort = vort - mesh%sub_face(id_sub_face)%area/mesh%sub_elem(id_sub_elem)%volume&
            *cross_product(norm, w(2:4))
        end do
      end do
      write (vtkout, *) vort(1), vort(2), vort(3)
    end do

    close (vtkout)
  end subroutine write_sol_vtk

  subroutine print_mat(mat)
    implicit none

    integer(kind=ENTIER) :: i, j
    real(kind=DOUBLE), intent(in) :: mat(:, :)

    print *, ""
    do i = 1, size(mat(:, 1))
      do j = 1, size(mat(i, :))
        if (mat(i, j) > 1e-12_DOUBLE) then
          write (*, '(a,f8.4,a)', advance="no") ""//achar(27)//"[36m", mat(i, j), ""//achar(27)//"[0m"
        else if (mat(i, j) < -1e-12_DOUBLE) then
          write (*, '(a,f8.4,a)', advance="no") ""//achar(27)//"[34m", mat(i, j), ""//achar(27)//"[0m"
        else
          write (*, '(f8.4)', advance="no") mat(i, j)
        end if
      end do
      write (*, *) " "
    end do
    print *, ""
  end subroutine print_mat

  pure subroutine compute_error_sod(mesh, sol, t, x_sod, sol_left, sol_right, gamma)
    implicit none

    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), dimension(5), intent(in) :: sol_left, sol_right
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(inout) :: sol
    real(kind=DOUBLE), intent(in) :: t, gamma, x_sod

    integer(kind=ENTIER) :: n_iter, n_iter_max, i
    real(kind=DOUBLE) :: x1, x2, x3, x4
    real(kind=DOUBLE) :: s1, s2, s3, s4
    real(kind=DOUBLE) :: p1, p2, p3, p4, p5
    real(kind=DOUBLE) :: u1, u2, u3, u4, u5
    real(kind=DOUBLE) :: rho1, rho2, rho3, rho4, rho5
    real(kind=DOUBLE) :: al, ar
    real(kind=DOUBLE) :: eps, res, gamma_2
    real(kind=DOUBLE) :: sqrt1, sqrt2, beta
    real(kind=DOUBLE) :: error, tot_volume

    rho1 = sol_left(1)
    u1 = sol_left(2)
    p1 = sol_left(5)

    rho5 = sol_right(1)
    u5 = sol_right(2)
    p5 = sol_right(5)

    al = sqrt(gamma*p1/rho1)
    ar = sqrt(gamma*p5/rho5)

    s1 = u1 - al
    x1 = x_sod + (u1 - al)*t

    gamma_2 = (gamma - 1.0_DOUBLE)/(gamma + 1.0_DOUBLE)
    beta = (gamma - 1.0_DOUBLE)/(2.0_DOUBLE*gamma)

    eps = 1e-8_DOUBLE
    p3 = 0.5_DOUBLE*(p1 + p5)
    res = 42.0_DOUBLE
    n_iter_max = 10000
    n_iter = 0
    do while (abs(res) > eps .and. n_iter < n_iter_max)
      sqrt1 = sqrt((1.0_DOUBLE - gamma_2)/(sol_right(1)*(p3 + gamma_2*p5)))
      sqrt2 = sqrt((1.0_DOUBLE - gamma_2**2)*p1**(1.0_DOUBLE/gamma)/(gamma_2**2*sol_left(1)))
      res = (4.0_DOUBLE*p3 + p5 + (p1**beta - p3**beta)*sqrt2/sqrt1)/5.0_DOUBLE - p3
      p3 = p3 + res
      u3 = (p3 - p5)*sqrt1
      u4 = (p1**beta - p3**beta)*sqrt2
      n_iter = n_iter + 1
    end do

    p4 = p3

    rho3 = rho1*(p3/p1)**(1.0_DOUBLE/gamma)
    rho4 = rho5*(p4 + gamma_2*p5)/(p5 + gamma_2*p4)

    s4 = (rho5*u5 - rho4*u4)/(rho5 - rho4)
    x4 = x_sod + s4*t

    s2 = u3 - sqrt(gamma*p3/(rho3))
    x2 = x_sod + s2*t

    s3 = u3
    x3 = x_sod + s3*t

    error = 0.0_DOUBLE
    tot_volume = 0.0_DOUBLE
    do i = 1, mesh%n_elems

      if (mesh%elem(i)%coord(1) < x1) then
        sol(1, i) = rho1
        sol(2, i) = u1
        sol(3, i) = 0.0_DOUBLE
        sol(4, i) = 0.0_DOUBLE
        sol(5, i) = p1
      else if (mesh%elem(i)%coord(1) > x1 .and. mesh%elem(i)%coord(1) < x2) then
        u2 = 2.0_DOUBLE/(gamma + 1.0_DOUBLE)*(al + (mesh%elem(i)%coord(1) - x_sod)/t)
        rho2 = rho1*(1 - (gamma - 1)/2.0_DOUBLE*u2/al)**(2.0_DOUBLE/(gamma - 1))
        p2 = p1*(1 - (gamma - 1)/2.0_DOUBLE*u2/al)**(2.0_DOUBLE*gamma/(gamma - 1))
        sol(1, i) = rho2
        sol(2, i) = u2
        sol(3, i) = 0.0_DOUBLE
        sol(4, i) = 0.0_DOUBLE
        sol(5, i) = p2
      else if (mesh%elem(i)%coord(1) > x2 .and. mesh%elem(i)%coord(1) < x3) then
        sol(1, i) = rho3
        sol(2, i) = u3
        sol(3, i) = 0.0_DOUBLE
        sol(4, i) = 0.0_DOUBLE
        sol(5, i) = p3
      else if (mesh%elem(i)%coord(1) > x3 .and. mesh%elem(i)%coord(1) < x4) then
        sol(1, i) = rho4
        sol(2, i) = u4
        sol(3, i) = 0.0_DOUBLE
        sol(4, i) = 0.0_DOUBLE
        sol(5, i) = p4
      else if (mesh%elem(i)%coord(1) > x4) then
        sol(1, i) = rho5
        sol(2, i) = u5
        sol(3, i) = 0.0_DOUBLE
        sol(4, i) = 0.0_DOUBLE
        sol(5, i) = p5
      end if

      error = error + mesh%elem(i)%volume*(sol(1, i) - sol(1, i))**2
      tot_volume = tot_volume + mesh%elem(i)%volume
    end do
  end subroutine compute_error_sod

  subroutine compute_error_isentropic(mesh, sol, t)
    use global_data_module, only: error_2d, error_2d_h
    use euler_module, only: sol_isentropic_vortex, conserv_to_primit
    implicit none

    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(in) :: sol
    real(kind=DOUBLE), intent(in) :: t

    integer(kind=ENTIER) :: i
    real(kind=DOUBLE) :: error, volume
    real(kind=DOUBLE), dimension(5) :: wexact, wsol

    error = 0.0_DOUBLE
    volume = 0.0_DOUBLE
    do i = 1, mesh%n_elems
      !if( norm2(mesh%elem(i)%coord) < 3.0_DOUBLE ) then
        call sol_isentropic_vortex(mesh%elem(i)%coord, wexact, t)
        wsol = conserv_to_primit(sol(:, i))
        error = error + mesh%elem(i)%volume*(wsol(1) - wexact(1))**2
        volume = volume + mesh%elem(i)%volume
      !end if
    end do

    if (error_2d) then
      print *, "Error Vortex: ", sqrt((volume/error_2d_h)/mesh%n_elems), sqrt(error)
    else
      print *, "Error Vortex: ", (volume/mesh%n_elems)**(1.0_DOUBLE/3.0_DOUBLE), sqrt(error)
    end if
  end subroutine compute_error_isentropic

  subroutine print_send_recv_info(mpi_send_recv)
    use subfv_mpi_module
    implicit none

    type(mpi_send_recv_type), intent(in) :: mpi_send_recv

    integer(kind=ENTIER) :: i

    print *, "Send neigh", mpi_send_recv%n_mpi_send_neigh
    do i = 1, mpi_send_recv%n_mpi_send_neigh
      print *, "Send ", mpi_send_recv%mpi_send_neigh(i)%n_elems, " to ", &
        mpi_send_recv%mpi_send_neigh(i)%partition_id
    end do

    print *, "Recv neigh", mpi_send_recv%n_mpi_recv_neigh
    do i = 1, mpi_send_recv%n_mpi_recv_neigh
      print *, "Receive ", mpi_send_recv%mpi_recv_neigh(i)%n_elems, " from ", &
        mpi_send_recv%mpi_recv_neigh(i)%partition_id
    end do
  end subroutine print_send_recv_info

  subroutine write_basic_sol_vtk(mesh, filename, sol)
    implicit none

    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(in) :: sol
    character(len=*), intent(in) :: filename

    integer(kind=ENTIER) :: i, j, size_tot, vtkout

    open (newunit=vtkout, file=trim(filename), status="unknown")

    write (vtkout, '(A)') '# vtk DataFile Version 3.0'
    write (vtkout, '(A)') filename
    write (vtkout, '(A)') 'ASCII'
    write (vtkout, '(A)') 'DATASET UNSTRUCTURED_GRID'

    write (vtkout, *) 'POINTS ', mesh%n_vert, ' double'
    do i = 1, mesh%n_vert
      write (vtkout, *) mesh%vert(i)%coord(1), mesh%vert(i)%coord(2), &
        mesh%vert(i)%coord(3)
    end do

    size_tot = 0
    do i = 1, mesh%n_elems
      if (.not. mesh%elem(i)%is_ghost) then
        size_tot = size_tot + mesh%elem(i)%n_vert + 1
      end if
    end do

    write (vtkout, *) 'CELLS ', mesh%n_interior_elems, size_tot
    do i = 1, mesh%n_elems
      if (.not. mesh%elem(i)%is_ghost) then
        write (vtkout, '(i12)', advance='no') mesh%elem(i)%n_vert
        do j = 1, mesh%elem(i)%n_vert - 1
          write (vtkout, '(i12)', advance='no') mesh%elem(i)%vert(j) - 1
        end do
        write (vtkout, '(i12)') mesh%elem(i)%vert(mesh%elem(i)%n_vert) - 1
      end if
    end do

    write (vtkout, '(a,i11)') 'CELL_TYPES ', mesh%n_interior_elems
    do i = 1, mesh%n_elems
      if (.not. mesh%elem(i)%is_ghost) then
        if (mesh%elem(i)%elem_kind == 4) then
          write (vtkout, '(i12)') 10
        else if (mesh%elem(i)%elem_kind == 5) then
          write (vtkout, '(i12)') 12
        else if (mesh%elem(i)%elem_kind == 6) then
          write (vtkout, '(i12)') 13
        else if (mesh%elem(i)%elem_kind == 7) then
          write (vtkout, '(i12)') 14
        end if
      end if
    end do

    write (vtkout, '(a,i11)') 'CELL_DATA ', mesh%n_interior_elems
    write (vtkout, *) 'SCALARS Centroid double', 3
    write (vtkout, *) 'LOOKUP_TABLE default'
    do i = 1, mesh%n_elems
      if (.not. mesh%elem(i)%is_ghost) then
        write (vtkout, *) mesh%elem(i)%coord
      end if
    end do

    write (vtkout, *) 'SCALARS Sol double', 5
    write (vtkout, *) 'LOOKUP_TABLE default'
    do i = 1, mesh%n_elems
      if (.not. mesh%elem(i)%is_ghost) then
        write (vtkout, *) sol(:, i)
      end if
    end do

    close (vtkout)
  end subroutine write_basic_sol_vtk

  subroutine write_sol_vtk_wasilij(mesh, filename, t, &
      sol, grad, residu)
    use global_data_module
    use euler_module, only: sol_isentropic_vortex, conserv_to_primit
    use subfv_linear_solver_module
    implicit none

    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), intent(in) :: t
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(in) :: sol
    real(kind=DOUBLE), dimension(:, :, :), intent(in) :: grad
    real(kind=DOUBLE), dimension(mesh%n_elems), intent(in) :: residu

    character(len=*), intent(in) :: filename

    integer(kind=ENTIER) :: i, j, size_tot, vtkout
    real(kind=DOUBLE) :: r
    real(kind=DOUBLE), dimension(5) :: w

    open (newunit=vtkout, file=trim(filename), status="unknown")

    write (vtkout, '(A)') '# vtk DataFile Version 3.0'
    write (vtkout, '(A)') filename
    write (vtkout, '(A)') 'ASCII'
    write (vtkout, '(A)') 'DATASET UNSTRUCTURED_GRID'

    write (vtkout, *) 'POINTS ', mesh%n_vert/2, ' double'
    do i = 1, mesh%n_vert/2
      write (vtkout, *) mesh%vert(i)%coord(1), mesh%vert(i)%coord(2), &
        mesh%vert(i)%coord(3)
    end do

    size_tot = 0
    do i = 1, mesh%n_elems
      if (.not. mesh%elem(i)%is_ghost) then
        size_tot = size_tot + mesh%elem(i)%n_vert/2 + 1
      end if
    end do

    write (vtkout, *) 'CELLS ', mesh%n_interior_elems, size_tot
    do i = 1, mesh%n_elems
      if (.not. mesh%elem(i)%is_ghost) then
        write (vtkout, '(i12)', advance='no') mesh%elem(i)%n_vert/2
        do j = 1, mesh%elem(i)%n_vert/2 - 1
          write (vtkout, '(i12)', advance='no') mesh%elem(i)%vert(j) - 1
        end do
        write (vtkout, '(i12)') mesh%elem(i)%vert(mesh%elem(i)%n_vert/2) - 1
      end if
    end do

    write (vtkout, '(a,i11)') 'CELL_TYPES ', mesh%n_interior_elems
    do i = 1, mesh%n_elems
      if (.not. mesh%elem(i)%is_ghost) then
        write (vtkout, '(i12)') 7
      end if
    end do

    write (vtkout, '(a,i11)') 'CELL_DATA ', mesh%n_interior_elems
    write (vtkout, *) 'SCALARS Centroid double', 3
    write (vtkout, *) 'LOOKUP_TABLE default'
    do i = 1, mesh%n_elems
      if (.not. mesh%elem(i)%is_ghost) then
        write (vtkout, *) mesh%elem(i)%coord
      end if
    end do

    write (vtkout, *) 'SCALARS r double', 1
    write (vtkout, *) 'LOOKUP_TABLE default'
    do i = 1, mesh%n_elems
      if (.not. mesh%elem(i)%is_ghost) then
        write (vtkout, *) sqrt(mesh%elem(i)%coord(1)**2 + mesh%elem(i)%coord(2)**2)
      end if
    end do

    write (vtkout, *) 'SCALARS Theta1 double', 1
    write (vtkout, *) 'LOOKUP_TABLE default'
    do i = 1, mesh%n_elems
      if (.not. mesh%elem(i)%is_ghost) then
        write (vtkout, *) atan2(mesh%elem(i)%coord(2), mesh%elem(i)%coord(1))
      end if
    end do

    write (vtkout, *) 'SCALARS Theta2 double', 1
    write (vtkout, *) 'LOOKUP_TABLE default'
    do i = 1, mesh%n_elems
      if (.not. mesh%elem(i)%is_ghost) then
        write (vtkout, *) atan2(-mesh%elem(i)%coord(2), mesh%elem(i)%coord(1))
      end if
    end do

    write (vtkout, *) 'SCALARS Theta3 double', 1
    write (vtkout, *) 'LOOKUP_TABLE default'
    do i = 1, mesh%n_elems
      if (.not. mesh%elem(i)%is_ghost) then
        write (vtkout, *) atan2(-mesh%elem(i)%coord(2), -mesh%elem(i)%coord(1))
      end if
    end do

    write (vtkout, *) 'SCALARS Theta4 double', 1
    write (vtkout, *) 'LOOKUP_TABLE default'
    do i = 1, mesh%n_elems
      if (.not. mesh%elem(i)%is_ghost) then
        write (vtkout, *) atan2(mesh%elem(i)%coord(2), -mesh%elem(i)%coord(1))
      end if
    end do

    write (vtkout, *) 'SCALARS Density double', 1
    write (vtkout, *) 'LOOKUP_TABLE default'
    do i = 1, mesh%n_elems
      if (.not. mesh%elem(i)%is_ghost) then
        write (vtkout, *) sol(1, i)
      end if
    end do

    write (vtkout, *) 'VECTORS Momentum double'
    do i = 1, mesh%n_elems
      if (.not. mesh%elem(i)%is_ghost) then
        write (vtkout, *) sol(2, i), sol(3, i), sol(4, i)
      end if
    end do

    write (vtkout, *) 'SCALARS Energy double', 1
    write (vtkout, *) 'LOOKUP_TABLE default'
    do i = 1, mesh%n_elems
      if (.not. mesh%elem(i)%is_ghost) then
        write (vtkout, *) sol(5, i)
      end if
    end do

    write (vtkout, *) 'VECTORS Velocity double'
    do i = 1, mesh%n_elems
      if (.not. mesh%elem(i)%is_ghost) then
        write (vtkout, *) sol(2, i)/sol(1, i), &
          sol(3, i)/sol(1, i), sol(4, i)/sol(1, i)
      end if
    end do

    write (vtkout, *) 'SCALARS Internal_energy double', 1
    write (vtkout, *) 'LOOKUP_TABLE default'
    do i = 1, mesh%n_elems
      if (.not. mesh%elem(i)%is_ghost) then
        write (vtkout, *) sol(5, i)/sol(1, i) - 0.5_DOUBLE*norm2(sol(2:4, i)/sol(1, i))**2
      end if
    end do

    write (vtkout, *) 'SCALARS Pressure double', 1
    write (vtkout, *) 'LOOKUP_TABLE default'
    do i = 1, mesh%n_elems
      if (.not. mesh%elem(i)%is_ghost) then
        w = conserv_to_primit(sol(:, i))
        write (vtkout, *) w(5)
      end if
    end do

    write (vtkout, *) 'SCALARS Mach double', 1
    write (vtkout, *) 'LOOKUP_TABLE default'
    do i = 1, mesh%n_elems
      if (.not. mesh%elem(i)%is_ghost) then
        w = conserv_to_primit(sol(:, i))
        write (vtkout, *) norm2(w(2:4))/sqrt(gamma*w(5)/w(1))
      end if
    end do

    write (vtkout, *) 'SCALARS H double', 1
    write (vtkout, *) 'LOOKUP_TABLE default'
    do i = 1, mesh%n_elems
      if (.not. mesh%elem(i)%is_ghost) then
        w = conserv_to_primit(sol(:, i))
        write (vtkout, *) (gamma*w(5)/w(1))/(gamma - 1) + 0.5_DOUBLE*norm2(w(2:4))**2
      end if
    end do

    write (vtkout, *) 'SCALARS Residu double', 1
    write (vtkout, *) 'LOOKUP_TABLE default'
    do i = 1, mesh%n_elems
      if (.not. mesh%elem(i)%is_ghost) then
        write (vtkout, *) residu(i)
      end if
    end do

    if (init_noh) then
      write (vtkout, *) 'SCALARS Density_exact_noh double', 1
      write (vtkout, *) 'LOOKUP_TABLE default'
      do i = 1, mesh%n_elems
        if (.not. mesh%elem(i)%is_ghost) then
          r = sqrt(mesh%elem(i)%coord(1)**2 + mesh%elem(i)%coord(2)**2)
          if (r < t/3.0_DOUBLE) then
            write (vtkout, *) ((gamma + 1)/(gamma - 1))**2
          else
            write (vtkout, *) (1.0_DOUBLE + t/r)
          end if
        end if
      end do

      write (vtkout, *) 'SCALARS Pressure_exact_noh double', 1
      write (vtkout, *) 'LOOKUP_TABLE default'
      do i = 1, mesh%n_elems
        if (.not. mesh%elem(i)%is_ghost) then
          r = sqrt(mesh%elem(i)%coord(1)**2 + mesh%elem(i)%coord(2)**2)
          if (r < t/3.0_DOUBLE) then
            write (vtkout, *) 16.0_DOUBLE/3.0_DOUBLE
          else
            write (vtkout, *) 1e-6
          end if
        end if
      end do

      write (vtkout, *) 'SCALARS Internal_energy_exact_noh double', 1
      write (vtkout, *) 'LOOKUP_TABLE default'
      do i = 1, mesh%n_elems
        if (.not. mesh%elem(i)%is_ghost) then
          r = sqrt(mesh%elem(i)%coord(1)**2 + mesh%elem(i)%coord(2)**2)
          if (r < t/3.0_DOUBLE) then
            write (vtkout, *) 1.0_DOUBLE/2.0_DOUBLE
          else
            write (vtkout, *) 1e-6
          end if
        end if
      end do

      write (vtkout, *) 'VECTORS Velocity_Mexact_noh double'
      do i = 1, mesh%n_elems
        if (.not. mesh%elem(i)%is_ghost) then
          r = sqrt(mesh%elem(i)%coord(1)**2 + mesh%elem(i)%coord(2)**2)
          if (r < t/3.0_DOUBLE) then
            write (vtkout, *) 0.0_DOUBLE, 0.0_DOUBLE, 0.0_DOUBLE
          else
            write (vtkout, *) - mesh%elem(i)%coord(1)/r, -mesh%elem(i)%coord(2)/r, 0.0_DOUBLE
          end if
        end if
      end do
    end if

    if (init_isentropic_vortex) then
      call compute_error_isentropic(mesh, sol, t)
      write (vtkout, *) 'SCALARS Density_exact_isentropic double', 1
      write (vtkout, *) 'LOOKUP_TABLE default'
      do i = 1, mesh%n_elems
        if (.not. mesh%elem(i)%is_ghost) then
          call sol_isentropic_vortex(mesh%elem(i)%coord, w, t)
          write (vtkout, *) w(1)
        end if
      end do
    end if

    if (compute_cp) then
      write (vtkout, *) 'SCALARS Cp double', 1
      write (vtkout, *) 'LOOKUP_TABLE default'
      do i = 1, mesh%n_elems
        if (.not. mesh%elem(i)%is_ghost) then
          w = conserv_to_primit(sol(:, i))
          write (vtkout, *) (w(5) - pinf)/(0.5_DOUBLE*rhoinf*vinf**2)
        end if
      end do
    end if

    write (vtkout, *) 'SCALARS Tag int', 1
    write (vtkout, *) 'LOOKUP_TABLE default'
    do i = 1, mesh%n_elems
      if (.not. mesh%elem(i)%is_ghost) then
        write (vtkout, *) mesh%elem(i)%tag
      end if
    end do

    if (second_order) then
      write (vtkout, *) 'VECTORS Density_grad double'
      do i = 1, mesh%n_elems
        if (.not. mesh%elem(i)%is_ghost) then
          write (vtkout, *) grad(1, 1, i), grad(2, 1, i), grad(3, 1, i)
        end if
      end do

      write (vtkout, *) 'VECTORS Pressure_grad double'
      do i = 1, mesh%n_elems
        if (.not. mesh%elem(i)%is_ghost) then
          write (vtkout, *) grad(1, 5, i), grad(2, 5, i), grad(3, 5, i)
        end if
      end do
    end if

    close (vtkout)
  end subroutine write_sol_vtk_wasilij

  subroutine write_sol_meta_pvtu(filename, iaff_char)
    use mpi
    use subfv_mpi_module
    use global_data_module, only: second_order
    implicit none

    character(len=*), intent(in) :: filename, iaff_char

    integer(kind=ENTIER) :: fn, i, num_procs, me, mpi_ierr
    character(len=255) :: i_char

    call MPI_COMM_SIZE(MPI_COMM_WORLD, num_procs, mpi_ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, me, mpi_ierr)

    if( me == 0 ) then
      open(newunit=fn, file=trim(adjustl(filename)))

      write(fn, *) "<VTKFile type='PUnstructuredGrid' version='0.1' byte_order='LittleEndian' header_type='UInt64'>"
      write(fn, *) "<PUnstructuredGrid GhostLevel='0'>"
      write(fn, *) "<PPoints>"
      write(fn, *) "<PDataArray type='Float64' Name='Points' NumberOfComponents='3'/>"
      write(fn, *) "</PPoints>"
      write(fn, *) "<PCells>"
      write(fn, *) "<PDataArray type='Int64' Name='connectivity' NumberOfComponents='1'/>"
      write(fn, *) "<PDataArray type='Int64' Name='offsets'      NumberOfComponents='1'/>"
      write(fn, *) "<PDataArray type='UInt8' Name='types'        NumberOfComponents='1'/>"
      write(fn, *) "</PCells>"
      write(fn, *) "<PPointData>"
      write(fn, *) "<PDataArray type='Float64' Name='Alpha_mcfl' NumberOfComponents='1'/>"
      write(fn, *) "<PDataArray type='Float64' Name='Nodal_Grad_Density' NumberOfComponents='3'/>"
      write(fn, *) "<PDataArray type='Float64' Name='Nodal_Grad_Velocity_X' NumberOfComponents='3'/>"
      write(fn, *) "<PDataArray type='Float64' Name='Nodal_Grad_Velocity_Y' NumberOfComponents='3'/>"
      write(fn, *) "<PDataArray type='Float64' Name='Nodal_Grad_Velocity_Z' NumberOfComponents='3'/>"
      write(fn, *) "<PDataArray type='Float64' Name='Nodal_Grad_Pressure' NumberOfComponents='3'/>"
      write(fn, *) "</PPointData>"
      write(fn, *) "<PCellData>"
      write(fn, *) "<PDataArray type='Float64' Name='Centroid' NumberOfComponents='3'/>"
      write(fn, *) "<PDataArray type='Int32' Name='Tag' NumberOfComponents='1'/>"
      write(fn, *) "<PDataArray type='Int32' Name='MPI_Color' NumberOfComponents='1'/>"
      write(fn, *) "<PDataArray type='Float64' Name='Density' NumberOfComponents='1'/>"
      write(fn, *) "<PDataArray type='Float64' Name='Momentum' NumberOfComponents='3'/>"
      write(fn, *) "<PDataArray type='Float64' Name='Energy' NumberOfComponents='1'/>"
      write(fn, *) "<PDataArray type='Float64' Name='Velocity' NumberOfComponents='3'/>"
      write(fn, *) "<PDataArray type='Float64' Name='Pressure' NumberOfComponents='1'/>"
      write(fn, *) "<PDataArray type='Float64' Name='Internal_energy' NumberOfComponents='1'/>"
      write(fn, *) "<PDataArray type='Float64' Name='H' NumberOfComponents='1'/>"
      write(fn, *) "<PDataArray type='Float64' Name='Mach' NumberOfComponents='1'/>"
      write(fn, *) "<PDataArray type='Float64' Name='Residu' NumberOfComponents='1'/>"
      write(fn, *) "<PDataArray type='Float64' Name='Temperature' NumberOfComponents='1'/>"
      write(fn, *) "<PDataArray type='Float64' Name='r' NumberOfComponents='1'/>"
      write(fn, *) "<PDataArray type='Float64' Name='rxy' NumberOfComponents='1'/>"
      write(fn, *) "<PDataArray type='Float64' Name='Theta3' NumberOfComponents='1'/>"
      if( second_order ) then
        write(fn, *) "<PDataArray type='Float64' Name='Density_Grad' NumberOfComponents='3'/>"
        write(fn, *) "<PDataArray type='Float64' Name='Velocity_X_Grad' NumberOfComponents='3'/>"
        write(fn, *) "<PDataArray type='Float64' Name='Velocity_Y_Grad' NumberOfComponents='3'/>"
        write(fn, *) "<PDataArray type='Float64' Name='Velocity_Z_Grad' NumberOfComponents='3'/>"
        write(fn, *) "<PDataArray type='Float64' Name='Pressure_Grad' NumberOfComponents='3'/>"
      end if
      write(fn, *) "</PCellData>"

      do i=0, num_procs-1
        write(i_char, *) i
        write(fn, *) "<Piece Source='"//trim(adjustl(i_char))//"_output_"//&
          trim(adjustl(iaff_char))//".vtu'/>"
      end do

      write(fn, *) "</PUnstructuredGrid>"
      write(fn, *) "</VTKFile>"

      close(fn)
    end if

    call mpi_barrier(mpi_comm_world, mpi_ierr)
  end subroutine write_sol_meta_pvtu

  subroutine write_coeffs_vtu(mesh, filename, sol)
    use global_data_module, only: compute_coeffs, coeffs_surf, &
      pinf, rhoinf, vinf, Cv_p, gamma
    use mpi
    use subfv_mpi_module
    use euler_module, only: conserv_to_primit
    use subfv_linear_solver_module, only: lu_inverse, tensor_product
    implicit none

    type(mesh_type), intent(in) :: mesh
    character(len=*), intent(in) :: filename
    real(kind=DOUBLE), dimension(5, mesh%n_elems) :: sol

    integer(kind=ENTIER) :: size_tot, fn
    integer(kind=ENTIER) :: i, j, k
    integer(kind=ENTIER) :: id_elem, id_face, id_sub_face, id_vert
    integer(kind=ENTIER) :: id_sub_elem, id_sub_elem_loc, id_sub_face_k_loc
    integer(kind=ENTIER) :: id_elem_k, id_face_k, id_sub_face_k, id_sub_elem_k
    integer(kind=ENTIER) :: n_bdy_faces
    real(kind=DOUBLE), dimension(5) :: w
    integer(kind=ENTIER), dimension(:), allocatable :: face_ids
    real(kind=DOUBLE), dimension(:), allocatable :: cp

    integer(kind=ENTIER), dimension(:), allocatable :: id_vert_no_ghost
    integer(kind=ENTIER), dimension(:), allocatable :: local_id_vert_no_ghost
    real(kind=DOUBLE) :: dmin, dmax
    integer(kind=ENTIER) :: me, num_procs, mpi_ierr, sm
    integer(kind=ENTIER) :: n_interior_elems, n_interior_vert

    !Get the face ids from the boundary coeffs_surf
    n_bdy_faces = 0
    do i = 1, mesh%n_faces
      if (-mesh%face(i)%right_neigh == coeffs_surf) n_bdy_faces = n_bdy_faces + 1
    end do
    allocate (face_ids(n_bdy_faces))
    n_bdy_faces = 0
    do i = 1, mesh%n_faces
      if (-mesh%face(i)%right_neigh == coeffs_surf) then
        n_bdy_faces = n_bdy_faces + 1
        face_ids(n_bdy_faces) = i
      end if
    end do

    allocate (cp(n_bdy_faces))
    !Compute Cp for the face
    do i = 1, n_bdy_faces
      id_face = face_ids(i)
      id_elem = mesh%face(id_face)%left_neigh
      w = conserv_to_primit(sol(:, id_elem))
      cp(i) = (w(5) - pinf)/(0.5_DOUBLE*rhoinf*vinf**2)
    end do

    !Write VTU Data
    call MPI_COMM_SIZE(MPI_COMM_WORLD, num_procs, mpi_ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, me, mpi_ierr)

    allocate(id_vert_no_ghost(mesh%n_vert))
    id_vert_no_ghost = 0
    n_interior_vert = 0
    do i=1, mesh%n_vert
      if ( .not. mesh%vert(i)%is_ghost ) then
        n_interior_vert = n_interior_vert + 1
        id_vert_no_ghost(i) = n_interior_vert
      end if
    end do

    n_interior_elems = 0
    do i = 1, n_bdy_faces
      id_face = face_ids(i)
      if( .not. mesh%elem(mesh%face(id_face)%left_neigh)%is_ghost ) n_interior_elems = n_interior_elems + 1
    end do

    open(newunit=fn, file=trim(adjustl(filename)))

    write(fn, *) "<VTKFile type='UnstructuredGrid' version='1.0' &
      &byte_order='LittleEndian' header_type='UInt64'>"
    write(fn, *) "<UnstructuredGrid>\n<Piece NumberOfPoints='", n_interior_vert, &
      &"' NumberOfCells='", n_interior_elems, "'>"

    dmin = 1e100_DOUBLE
    dmax = -1e100_DOUBLE
    do i=1, mesh%n_vert
      if( .not. mesh%vert(i)%is_ghost ) then
        do j=1, 3
          if( mesh%vert(i)%coord(j) > dmax ) dmax = mesh%vert(i)%coord(j)
          if( mesh%vert(i)%coord(j) < dmin ) dmin = mesh%vert(i)%coord(j)
        end do
      end if
    end do

    write(fn, *) "<Points>"
    write(fn, *) " <DataArray type='Float64' Name='Points' NumberOfComponents='3' &
      & format='ascii' RangeMin='", dmin, "' RangeMax='", dmax, "'>"

    do i=1, mesh%n_vert
      if( .not. mesh%vert(i)%is_ghost ) then
        write(fn, *) mesh%vert(i)%coord(:)
      end if
    end do

    write(fn, *) "</DataArray>"
    write(fn, *) "</Points>"

    write(fn, *) "<Cells>"
    write(fn, *) "<DataArray type='Int64' Name='connectivity' format='ascii' &
      &RangeMin='0' RangeMax='", mesh%n_vert-1, "'>"
    do i = 1, n_bdy_faces
      id_face = face_ids(i)
      if( .not. mesh%elem(mesh%face(id_face)%left_neigh)%is_ghost ) then
        write(fn, *) id_vert_no_ghost(mesh%face(id_face)%vert) - 1
      end if
    end do
    write(fn, *) "</DataArray>"

    sm = 0
    do i = 1, n_bdy_faces
      id_face = face_ids(i)
      if( .not. mesh%elem(mesh%face(id_face)%left_neigh)%is_ghost ) then
        sm = sm + mesh%face(id_face)%n_vert
      end if
    end do

    write(fn, *) "<DataArray type='Int64' Name='offsets' format='ascii' &
      &RangeMin='", 0, "' RangeMax='", sm, "'>"

    sm = 0
    do i = 1, n_bdy_faces
      id_face = face_ids(i)
      if( .not. mesh%elem(mesh%face(id_face)%left_neigh)%is_ghost ) then
        sm = sm + mesh%face(id_face)%n_vert
        write(fn, *) sm
      end if
    end do

    write(fn, *) "</DataArray>"

    write(fn, *) "<DataArray type='UInt8' Name='types' format='ascii' &
      &RangeMin='7' RangeMax='7'>"

    do i = 1, n_bdy_faces
      id_face = face_ids(i)
      if( .not. mesh%elem(mesh%face(id_face)%left_neigh)%is_ghost ) then
        write(fn, *) 7 !Type for polygon
      end if
    end do

    write(fn, *) "</DataArray>"
    write(fn, *) "</Cells>"

    write(fn, *) "<CellData>"

    write(fn, *) "<DataArray type='Float64' Name='Centroid' format='ascii' NumberOfComponents='3'>"
    do i = 1, n_bdy_faces
      id_face = face_ids(i)
      if( .not. mesh%elem(mesh%face(id_face)%left_neigh)%is_ghost ) then
        write (fn, *) mesh%face(id_face)%coord
      end if
    end do
    write(fn, *) "</DataArray>"

    write(fn, *) "<DataArray type='Float64' Name='Cp' format='ascii' NumberOfComponents='1'>"
    do i = 1, n_bdy_faces
      id_face = face_ids(i)
      if( .not. mesh%elem(mesh%face(id_face)%left_neigh)%is_ghost ) then
        write (fn, *) cp(i)
      end if
    end do
    write(fn, *) "</DataArray>"

    write(fn, *) "<DataArray type='Float64' Name='Density' format='ascii' NumberOfComponents='1'>"
    do i = 1, n_bdy_faces
      id_face = face_ids(i)
      if( .not. mesh%elem(mesh%face(id_face)%left_neigh)%is_ghost ) then
        id_elem = mesh%face(id_face)%left_neigh
        write (fn, *) sol(1, id_elem)
      end if
    end do
    write(fn, *) "</DataArray>"

    write(fn, *) "<DataArray type='Float64' Name='Temperature' format='ascii' NumberOfComponents='1'>"
    do i = 1, n_bdy_faces
      id_face = face_ids(i)
      if( .not. mesh%elem(mesh%face(id_face)%left_neigh)%is_ghost ) then
        id_elem = mesh%face(id_face)%left_neigh
        write (fn, *) (sol(5, id_elem)/sol(1, id_elem) &
          - 0.5_DOUBLE*norm2(sol(2:4, id_elem)/sol(1, id_elem))**2)/Cv_p
      end if
    end do
    write(fn, *) "</DataArray>"

    write(fn, *) "<DataArray type='Float64' Name='Pressure' format='ascii' NumberOfComponents='1'>"
    do i = 1, n_bdy_faces
      id_face = face_ids(i)
      if( .not. mesh%elem(mesh%face(id_face)%left_neigh)%is_ghost ) then
        id_elem = mesh%face(id_face)%left_neigh
        w = conserv_to_primit(sol(:, id_elem))
        write (fn, *) w(5)
      end if
    end do
    write(fn, *) "</DataArray>"

    write(fn, *) "<DataArray type='Float64' Name='Theta3' format='ascii' NumberOfComponents='1'>"
    do i = 1, n_bdy_faces
      id_face = face_ids(i)
      if( .not. mesh%elem(mesh%face(id_face)%left_neigh)%is_ghost ) then
        id_elem = mesh%face(id_face)%left_neigh
        write (fn, *) atan2(-mesh%elem(id_elem)%coord(2), -mesh%elem(id_elem)%coord(1))
      end if
    end do
    write(fn, *) "</DataArray>"

    write(fn, *) "</CellData>"

    write(fn, *) "</Piece>"
    write(fn, *) "</UnstructuredGrid>"
    write(fn, *) "</VTKFile>"
    close(fn)
  end subroutine write_coeffs_vtu

  subroutine write_sol_meta_pvtu_coeffs(filename, iaff_char)
    use mpi
    use subfv_mpi_module
    use global_data_module, only: second_order
    implicit none

    character(len=*), intent(in) :: filename, iaff_char

    integer(kind=ENTIER) :: fn, i, num_procs, me, mpi_ierr
    character(len=255) :: i_char

    call MPI_COMM_SIZE(MPI_COMM_WORLD, num_procs, mpi_ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, me, mpi_ierr)

    if( me == 0 ) then
      open(newunit=fn, file=trim(adjustl(filename)))

      write(fn, *) "<VTKFile type='PUnstructuredGrid' version='0.1' byte_order='LittleEndian' header_type='UInt64'>"
      write(fn, *) "<PUnstructuredGrid GhostLevel='0'>"
      write(fn, *) "<PPoints>"
      write(fn, *) "<PDataArray type='Float64' Name='Points' NumberOfComponents='3'/>"
      write(fn, *) "</PPoints>"
      write(fn, *) "<PCells>"
      write(fn, *) "<PDataArray type='Int64' Name='connectivity' NumberOfComponents='1'/>"
      write(fn, *) "<PDataArray type='Int64' Name='offsets'      NumberOfComponents='1'/>"
      write(fn, *) "<PDataArray type='UInt8' Name='types'        NumberOfComponents='1'/>"
      write(fn, *) "</PCells>"
      write(fn, *) "<PCellData>"
      write(fn, *) "<PDataArray type='Float64' Name='Centroid' NumberOfComponents='3'/>"
      write(fn, *) "<PDataArray type='Float64' Name='Density' NumberOfComponents='1'/>"
      write(fn, *) "<PDataArray type='Float64' Name='Pressure' NumberOfComponents='1'/>"
      write(fn, *) "<PDataArray type='Float64' Name='Temperature' NumberOfComponents='1'/>"
      write(fn, *) "<PDataArray type='Float64' Name='Theta3' NumberOfComponents='1'/>"
      write(fn, *) "<PDataArray type='Float64' Name='Cp' NumberOfComponents='1'/>"
      write(fn, *) "</PCellData>"

      do i=0, num_procs-1
        write(i_char, *) i
        write(fn, *) "<Piece Source='"//trim(adjustl(i_char))//"_coeffs_"//&
          trim(adjustl(iaff_char))//".vtu'/>"
      end do

      write(fn, *) "</PUnstructuredGrid>"
      write(fn, *) "</VTKFile>"

      close(fn)
    end if

    call mpi_barrier(mpi_comm_world, mpi_ierr)
  end subroutine write_sol_meta_pvtu_coeffs

  subroutine write_cell_size_dat(filename, mesh, sol, grad, mpi_send_recv)
    use mpi
    use subfv_mpi_module
    use global_data_module, only: second_order
    use euler_module, only: compute_cell_grad_from_nodal_grad
    implicit none

    character(len=*), intent(in) :: filename

    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(in) :: sol 
    real(kind=DOUBLE), dimension(3, 5, mesh%n_elems), intent(inout) :: grad
    type(mpi_send_recv_type), intent(inout) :: mpi_send_recv
    integer(kind=ENTIER) :: fn, i, num_procs, me, mpi_ierr, p, n_cells

    real(kind=DOUBLE), dimension(mesh%n_vert) :: cell_size


    call MPI_COMM_SIZE(MPI_COMM_WORLD, num_procs, mpi_ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, me, mpi_ierr)

    call compute_cell_grad_from_nodal_grad(mesh, sol, grad, 2)
    if (num_procs > 1) call mpi_memory_exchange(mpi_send_recv, &
      mesh%n_elems, 3*5, grad)
    call compute_cell_size(mesh, sol, grad, cell_size)

    do p=0, num_procs-1

      n_cells = 0
      do i=1, mesh%n_vert
        if( cell_size(i) > 0.0_DOUBLE) then
          n_cells = n_cells + 1
        end if
      end do
      call MPI_ALLREDUCE(MPI_IN_PLACE, n_cells, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD, mpi_ierr)

      if( me == p ) then
        if( p==0 ) then
          open(newunit=fn, file=trim(adjustl(filename)))
          write(fn, *) n_cells
        else
          open(newunit=fn, file=trim(adjustl(filename)), position='append', status='old')
        end if

        n_cells = 0
        do i=1, mesh%n_vert
          if( cell_size(i) > 0.0_DOUBLE) then
            n_cells = n_cells + 1
          end if
        end do

        do i=1, mesh%n_vert
          if( cell_size(i) > 0.0_DOUBLE) then
            write(fn, *) mesh%vert(i)%coord, cell_size(i)
          end if
        end do

        if( p==num_procs-1 ) write(fn, *) 0
        close(fn)
      end if
      call mpi_barrier(mpi_comm_world, mpi_ierr)
    end do
  end subroutine write_cell_size_dat

  subroutine compute_cell_size(mesh, sol, grad, cell_size)
    use global_data_module, only: gamma, Cv_p
    use euler_module, only: conserv_to_primit, compute_nodal_grad
    implicit none

    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), dimension(5, mesh%n_elems) :: sol
    real(kind=DOUBLE), dimension(3, 5, mesh%n_elems) :: grad
    real(kind=DOUBLE), dimension(mesh%n_vert) :: cell_size

    real(kind=DOUBLE), parameter :: reduc_ratio_1 = 1.0_DOUBLE
    real(kind=DOUBLE), parameter :: reduc_ratio_2 = 20.0_DOUBLE
    real(kind=DOUBLE), parameter :: pressure_grad_thresh = 0.1_DOUBLE
    real(kind=DOUBLE), parameter :: min_cell_size = 0.001_DOUBLE
    real(kind=DOUBLE), parameter :: max_cell_size = 10._DOUBLE
    !real(kind=DOUBLE), parameter :: delta_target_mach = 0.2_DOUBLE
    real(kind=DOUBLE), parameter :: delta_target_mach = 1.0_DOUBLE
    real(kind=DOUBLE), parameter :: delta_target_temperature = 0.02_DOUBLE

    integer(kind=ENTIER) :: i, j, id_elem
    real(kind=DOUBLE) :: ratio_temperature, ratio_mach, wp, a
    real(kind=DOUBLE), dimension(5) :: sol_p, sol_p_w
    real(kind=DOUBLE), dimension(3, 5) :: nodal_grad
    real(kind=DOUBLE) :: grad_mach, grad_temperature

    cell_size = 0.0_DOUBLE
    do i=1, mesh%n_vert

      wp = 0.0_DOUBLE
      sol_p = 0.0_DOUBLE
      do j = 1, mesh%vert(i)%n_elems_neigh
        id_elem = mesh%vert(i)%elem_neigh(j)
        wp = wp + mesh%elem(id_elem)%volume
        sol_p = sol_p + mesh%elem(id_elem)%volume*sol(:, id_elem)
      end do
      sol_p = sol_p/wp
      sol_p_w = conserv_to_primit(sol_p)
      wp = wp/mesh%vert(i)%n_elems_neigh


      call compute_nodal_grad(mesh, i, sol, nodal_grad)

      grad_temperature = abs(norm2(nodal_grad(:, 5))/sol_p_w(1) &
        - sol_p_w(5)/sol_p_w(1)**2 * norm2(nodal_grad(:, 1)))&
        /((gamma-1.0_DOUBLE)*Cv_p)

      a = sqrt(gamma*sol_p_w(5)/sol_p_w(1))
      grad_mach = abs(norm2(nodal_grad(:, 2:4))/a &
        - 0.5_DOUBLE*norm2(sol_p_w(2:4))/a**3*gamma*&
        abs(norm2(nodal_grad(:, 5))/sol_p_w(1) - sol_p_w(5)*norm2(nodal_grad(:, 1))/sol_p_w(1)**2))
      !grad_mach = norm2(nodal_grad(:, 5))

      ratio_temperature = grad_temperature*wp**(1.0_DOUBLE/3.0_DOUBLE)/delta_target_temperature + 1e-8_DOUBLE
      ratio_mach = grad_mach*wp**(1.0_DOUBLE/3.0_DOUBLE)/delta_target_mach + 1e-8_DOUBLE

      cell_size(i) = min(&
        max(wp**(1.0_DOUBLE/3.0_DOUBLE)/max(ratio_mach, ratio_temperature), min_cell_size), &
        max_cell_size)

    end do

    !Smooth
    !do i=1, mesh%n_vert
    !end do

    !Clip
    do i=1, mesh%n_vert
      if( cell_size(i) >= (1.0_DOUBLE-1e-8)*max_cell_size) cell_size(i) = -1.0_DOUBLE
    end do
  end subroutine compute_cell_size
end module io_module
