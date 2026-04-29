module euler_module
  use subfv_precision_module
  use subfv_mesh_module
  implicit none
contains
  subroutine compute_divF(mesh, sol, grad, divF)
    implicit none

    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(in) :: sol
    real(kind=DOUBLE), dimension(3, 5, mesh%n_elems), intent(in) :: grad
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(inout) :: divF

    integer(kind=ENTIER) :: i, k, p, idse, idv, idsf, idf
    real(kind=DOUBLE), dimension(3) :: norm
    real(kind=DOUBLE), dimension(5) :: sol_w_idv, sol_idv

    !$OMP PARALLEL DO PRIVATE(i, k, p, idse, idv, idsf, idf, &
    !$OMP norm, sol_w_idv, sol_idv)
    do i=1, mesh%n_elems
      divF(:, i) = 0.0_DOUBLE
      do k=1, mesh%elem(i)%n_sub_elems
        idse = mesh%elem(i)%sub_elem(k)
        idv = mesh%sub_elem(idse)%mesh_vert
        sol_w_idv = conserv_to_primit(sol(:, i))
        call reconstruct_with_grad(sol_w_idv, grad(:, :, i), &
          mesh%vert(idv)%coord - mesh%elem(i)%coord)
        sol_idv = primit_to_conserv(sol_w_idv)
        do p=1,mesh%sub_elem(idse)%n_sub_faces
          idsf = mesh%sub_elem(idse)%sub_face(p)
          idf = mesh%sub_face(idsf)%mesh_face
          if( mesh%sub_face(idsf)%left_elem_neigh == i ) then
            norm = mesh%face(idf)%norm
          else
            norm = -mesh%face(idf)%norm
          end if
          divF(:, i) = divF(:, i) + mesh%sub_face(idsf)%area*matmul(euler_flux(sol_idv), norm)
        end do
      end do
      divF(:, i) = divF(:, i)/mesh%elem(i)%volume
    end do
  end subroutine compute_divF

  subroutine compute_cell_grad_from_nodal_grad(mesh, sol, grad, method)
    implicit none

    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(in) :: sol
    real(kind=DOUBLE), dimension(3, 5, mesh%n_elems), intent(inout) :: grad
    !real(kind=DOUBLE), dimension(3, 5, mesh%n_vert), intent(inout) :: all_nodal_grad
    integer(kind=ENTIER), intent(in) :: method

    integer(kind=ENTIER) :: i, var
    real(kind=DOUBLE), dimension(5, mesh%n_elems) :: weight_sum

    grad = 0.0_DOUBLE
    weight_sum = 0.0_DOUBLE
    !$OMP PARALLEL DO
    do i=1, mesh%n_vert
      call compute_nodal_grad_cell_grad(mesh, i, sol, grad, weight_sum)
    end do

    !$OMP PARALLEL DO PRIVATE(var)
    do i=1, mesh%n_elems
      do var=1, 5
        grad(:, var, i) = grad(:, var, i) / weight_sum(var, i)
      end do
      call physical_slope_limiter(mesh, sol, i, grad)
      call relaxed_maximum_principle(mesh, sol, i, grad)
    end do
  end subroutine compute_cell_grad_from_nodal_grad

  subroutine relaxed_maximum_principle(mesh, sol, i, grad)
    use mpi
    implicit none

    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(in) :: sol
    integer(kind=ENTIER), intent(in) :: i
    real(kind=DOUBLE), dimension(3, 5, mesh%n_elems), intent(inout) :: grad

    integer(kind=ENTIER) :: j, id_vert, k, id_elem
    real(kind=DOUBLE) :: omega
    real(kind=DOUBLE), dimension(5) :: sol_w, delta_sol_w

    real(kind=DOUBLE) :: rho_min, rho_max, p_min, p_max

    real(kind=DOUBLE), parameter :: tol = 1e-8_DOUBLE

    integer(kind=ENTIER) :: me, mpi_ierr

    call MPI_COMM_RANK(MPI_COMM_WORLD, me, mpi_ierr)

    !Compute min max
    rho_min = 1e100_DOUBLE
    rho_max = 0.0_DOUBLE
    p_min = 1e100_DOUBLE
    p_max = 0.0_DOUBLE
    do j=1, mesh%elem(i)%n_vert
      id_vert = mesh%elem(i)%vert(j)
      do k=1, mesh%vert(id_vert)%n_elems_neigh
        id_elem = mesh%vert(id_vert)%elem_neigh(k)
        sol_w = conserv_to_primit(sol(:, id_elem))
        p_min = min(p_min, sol_w(5))
        p_max = max(p_max, sol_w(5))
        rho_min = min(rho_min, sol_w(1))
        rho_max = max(rho_max, sol_w(1))
      end do
    end do

    !Relax
    p_min = (1.0_DOUBLE-1e-2_DOUBLE)*p_min
    rho_min = (1.0_DOUBLE-1e-2_DOUBLE)*rho_min
    p_max = (1.0_DOUBLE+1e-2_DOUBLE)*p_max
    rho_max = (1.0_DOUBLE+1e-2_DOUBLE)*rho_max

    sol_w = conserv_to_primit(sol(:, i))
    omega = 1.0_DOUBLE
    do j=1, mesh%elem(i)%n_vert
      id_vert = mesh%elem(i)%vert(j)
      delta_sol_w = matmul(transpose(grad(:, :, i)), &
        mesh%vert(id_vert)%coord - mesh%elem(i)%coord)
      if( delta_sol_w(1) < 0.0_DOUBLE ) then
        omega = min(omega, &
          max(0.0_DOUBLE, (sol_w(1) - rho_min)/(tol - delta_sol_w(1))))
      else 
        omega = min(omega, &
          max(0.0_DOUBLE, -(sol_w(1) - rho_max)/(tol + delta_sol_w(1))))
      end if
      if( delta_sol_w(5) < 0.0_DOUBLE ) then
        omega = min(omega, &
          max(0.0_DOUBLE, (sol_w(5) - p_min)/(tol - delta_sol_w(5))))
      else
        omega = min(omega, &
          max(0.0_DOUBLE, -(sol_w(5) - p_max)/(tol + delta_sol_w(5))))
      end if
    end do
    grad(:, :, i) = omega * grad(:, :, i)
  end subroutine relaxed_maximum_principle

  subroutine physical_slope_limiter(mesh, sol, i, grad)
    use mpi
    implicit none

    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(in) :: sol
    integer(kind=ENTIER), intent(in) :: i
    real(kind=DOUBLE), dimension(3, 5, mesh%n_elems), intent(inout) :: grad

    integer(kind=ENTIER) :: j, id_vert, k, id_elem
    real(kind=DOUBLE) :: omega
    real(kind=DOUBLE), dimension(5) :: sol_w, delta_sol_w

    real(kind=DOUBLE), parameter :: tol = 1e-8_DOUBLE

    integer(kind=ENTIER) :: me, mpi_ierr

    call MPI_COMM_RANK(MPI_COMM_WORLD, me, mpi_ierr)

    sol_w = conserv_to_primit(sol(:, i))
    omega = 1.0_DOUBLE
    do j=1, mesh%elem(i)%n_vert
      id_vert = mesh%elem(i)%vert(j)
      delta_sol_w = matmul(transpose(grad(:, :, i)), &
        mesh%vert(id_vert)%coord - mesh%elem(i)%coord)
      if( delta_sol_w(1) < 0.0_DOUBLE ) then
        omega = min(omega, &
          max(0.0_DOUBLE, sol_w(1)/(tol - delta_sol_w(1))))
      end if
      if( delta_sol_w(5) < 0.0_DOUBLE ) then
        omega = min(omega, &
          max(0.0_DOUBLE, sol_w(5)/(tol - delta_sol_w(5))))
      end if
    end do
    grad(:, :, i) = omega * grad(:, :, i)
  end subroutine physical_slope_limiter

  subroutine compute_cell_grad_method_2(mesh, all_nodal_grad, id_elem, grad)
    implicit none

    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), dimension(3, 5, mesh%n_vert), intent(in) :: all_nodal_grad
    integer(kind=ENTIER), intent(in) :: id_elem
    real(kind=DOUBLE), dimension(3, 5), intent(inout) :: grad

    integer(kind=ENTIER) :: j, id_vert
    real(kind=DOUBLE) :: dual_vol_sum

    grad(:, :) = 0.0_DOUBLE
    dual_vol_sum = 0.0_DOUBLE
    do j=1, mesh%elem(id_elem)%n_vert
      id_vert = mesh%elem(id_elem)%vert(j)
      grad(:, :) = grad(:, :) &
        + mesh%vert(id_vert)%volume*all_nodal_grad(:, :, id_vert)
      dual_vol_sum = dual_vol_sum &
        + mesh%vert(id_vert)%volume
    end do
    grad = grad/dual_vol_sum
  end subroutine compute_cell_grad_method_2

  subroutine compute_cell_grad_method_3(mesh, all_nodal_grad, id_elem, grad)
    implicit none

    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), dimension(3, 5, mesh%n_vert), intent(in) :: all_nodal_grad
    integer(kind=ENTIER), intent(in) :: id_elem
    real(kind=DOUBLE), dimension(3, 5), intent(inout) :: grad

    integer(kind=ENTIER) :: j, id_vert, var

    id_vert = mesh%elem(id_elem)%vert(1)
    grad(:, :) = all_nodal_grad(:, :, id_vert)
    do j=2, mesh%elem(id_elem)%n_vert
      id_vert = mesh%elem(id_elem)%vert(j)
      do var=1, 5
        if(OI(all_nodal_grad(:, var, id_vert)) < OI(grad(:, var))) then
          grad(:, var) = all_nodal_grad(:, var, id_vert)
        end if
      end do
    end do
  end subroutine compute_cell_grad_method_3

  subroutine compute_cell_grad_method_4(mesh, all_nodal_grad, id_elem, grad)
    implicit none

    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), dimension(3, 5, mesh%n_vert), intent(in) :: all_nodal_grad
    integer(kind=ENTIER), intent(in) :: id_elem
    real(kind=DOUBLE), dimension(3, 5), intent(inout) :: grad

    real(kind=DOUBLE), parameter :: eps=1e-3_DOUBLE

    integer(kind=ENTIER) :: j, id_vert, var
    real(kind=DOUBLE), dimension(5) :: weight_sum

    !Non-linear combination with alpha=1
    grad(:, :) = 0.0_DOUBLE
    weight_sum(:) = 0.0_DOUBLE
    do j=1, mesh%elem(id_elem)%n_vert
      id_vert = mesh%elem(id_elem)%vert(j)
      do var=1, 5
        grad(:, var) = grad(:, var) &
          + mesh%vert(id_vert)%volume/(eps+OI(all_nodal_grad(:, var, id_vert))) &
          * all_nodal_grad(:, var, id_vert)
        weight_sum(var) = weight_sum(var) &
          + mesh%vert(id_vert)%volume/(eps+OI(all_nodal_grad(:, var, id_vert)))
      end do
    end do

    do var=1, 5
      grad(:, var) = grad(:, var)/weight_sum(var)
    end do
  end subroutine compute_cell_grad_method_4

  subroutine compute_cell_grad_method_5(mesh, all_nodal_grad, id_elem, grad)
    implicit none

    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), dimension(3, 5, mesh%n_vert), intent(in) :: all_nodal_grad
    integer(kind=ENTIER), intent(in) :: id_elem
    real(kind=DOUBLE), dimension(3, 5), intent(inout) :: grad

    real(kind=DOUBLE), parameter :: eps=1e-3_DOUBLE

    integer(kind=ENTIER) :: j, id_vert, var
    real(kind=DOUBLE), dimension(5) :: weight_sum

    !Non-linear combination with alpha=1
    grad(:, :) = 0.0_DOUBLE
    weight_sum(:) = 0.0_DOUBLE
    do j=1, mesh%elem(id_elem)%n_vert
      id_vert = mesh%elem(id_elem)%vert(j)
      do var=1, 5
        grad(:, var) = grad(:, var) &
          + mesh%vert(id_vert)%volume/(eps+OI(all_nodal_grad(:, var, id_vert)))**2 &
          * all_nodal_grad(:, var, id_vert)
        weight_sum(var) = weight_sum(var) &
          + mesh%vert(id_vert)%volume/(eps+OI(all_nodal_grad(:, var, id_vert)))**2
      end do
    end do

    do var=1, 5
      grad(:, var) = grad(:, var)/weight_sum(var)
    end do
  end subroutine compute_cell_grad_method_5

  subroutine compute_cell_grad_method_6(mesh, all_nodal_grad, id_elem, grad)
    implicit none

    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), dimension(3, 5, mesh%n_vert), intent(in) :: all_nodal_grad
    integer(kind=ENTIER), intent(in) :: id_elem
    real(kind=DOUBLE), dimension(3, 5), intent(inout) :: grad

    real(kind=DOUBLE), parameter :: eps=1e-3_DOUBLE

    integer(kind=ENTIER) :: j, id_vert, var
    real(kind=DOUBLE), dimension(5) :: weight_sum
    real(kind=DOUBLE) :: central_weight_sum
    real(kind=DOUBLE), dimension(3, 5) :: central_grad

    !Non-linear combination with alpha=1
    grad(:, :) = 0.0_DOUBLE
    weight_sum(:) = 0.0_DOUBLE
    central_grad(:, :) = 0.0_DOUBLE
    central_weight_sum = 0.0_DOUBLE
    do j=1, mesh%elem(id_elem)%n_vert
      id_vert = mesh%elem(id_elem)%vert(j)
      do var=1, 5
        grad(:, var) = grad(:, var) &
          + mesh%vert(id_vert)%volume/(eps+OI(all_nodal_grad(:, var, id_vert)))**2 &
          * all_nodal_grad(:, var, id_vert)
        weight_sum(var) = weight_sum(var) &
          + mesh%vert(id_vert)%volume/(eps+OI(all_nodal_grad(:, var, id_vert)))**2
      end do
      central_grad(:, :) = central_grad(:, :) &
        + mesh%vert(id_vert)%volume*all_nodal_grad(:, :, id_vert)
      central_weight_sum = central_weight_sum &
        + mesh%vert(id_vert)%volume
    end do

    central_grad = central_grad/central_weight_sum
    grad(:, :) = grad(:, :) &
      + central_weight_sum/(eps+OI(central_grad(:, :)))**2 &
      * central_grad(:, :)
    weight_sum(:) = weight_sum(:) &
      + central_weight_sum/(eps+OI(central_grad(:, :)))**2

    do var=1, 5
      grad(:, var) = grad(:, var)/weight_sum(var)
    end do
  end subroutine compute_cell_grad_method_6

  pure function OI(grad)
    implicit none

    real(kind=DOUBLE), dimension(3), intent(in) :: grad
    real(kind=DOUBLE) :: OI

    OI = norm2(grad(:))
  end function OI

  subroutine three_wave(sol_l, sol_r, sol_w_l, sol_w_r, lr_flux, sl, sr)
    use global_data_module, only: gamma
    implicit none

    real(kind=DOUBLE), dimension(5), intent(in) :: sol_l, sol_r
    real(kind=DOUBLE), dimension(5), intent(in) :: sol_w_l, sol_w_r
    real(kind=DOUBLE), dimension(5, 2), intent(inout) :: lr_flux
    real(kind=DOUBLE), intent(inout) :: sl, sr

    real(kind=DOUBLE) :: rhol, rhor, ul, vl, wl, ur, vr, wr
    real(kind=DOUBLE) :: pl, pr, rhol_et, rhor_et
    real(kind=DOUBLE) :: ul_et, vl_et, wl_et, ur_et, vr_et, wr_et, u_et
    real(kind=DOUBLE) :: el, er, al, ar
    real(kind=DOUBLE) :: u_bar, pl_bar, pr_bar
    real(kind=DOUBLE) :: lambda_l, lambda_r
    real(kind=DOUBLE), dimension(5) :: fl, fr
    real(kind=DOUBLE), dimension(5) :: sol_l_et, sol_r_et

    rhol = sol_w_l(1)
    ul = sol_w_l(2)
    vl = sol_w_l(3)
    wl = sol_w_l(4)
    pl = sol_w_l(5)
    el = sol_l(5)/rhol
    al = sqrt(gamma*pl/rhol)

    rhor = sol_w_r(1)
    ur = sol_w_r(2)
    vr = sol_w_r(3)
    wr = sol_w_r(4)
    pr = sol_w_r(5)
    er = sol_r(5)/rhor
    ar = sqrt(gamma*pr/rhor)

    fl(:) = ul*sol_l(:) + (/0.0_DOUBLE, pl, 0.0_DOUBLE, 0.0_DOUBLE, pl*ul/)
    fr(:) = ur*sol_r(:) + (/0.0_DOUBLE, pr, 0.0_DOUBLE, 0.0_DOUBLE, pr*ur/)

    !Compute lambdas from two-point positivity conditions
    lambda_l = max(al*rhol, sqrt(rhol*max(0.0_DOUBLE, pr - pl)), -rhol*(ur - ul))
    lambda_r = max(ar*rhor, sqrt(rhor*max(0.0_DOUBLE, pl - pr)), -rhor*(ur - ul))
    u_bar = (lambda_l*ul + lambda_r*ur - (pr - pl))/(lambda_r + lambda_l)

    u_et = u_bar

    rhol_et = 1.0_DOUBLE/(1.0_DOUBLE/rhol + (u_et - ul)/lambda_l)
    ul_et = u_et
    vl_et = vl
    wl_et = wl
    pl_bar = pl - lambda_l*(u_et - ul)

    sol_l_et(1) = rhol_et
    sol_l_et(2) = rhol_et*ul_et
    sol_l_et(3) = rhol_et*vl_et
    sol_l_et(4) = rhol_et*wl_et
    sol_l_et(5) = rhol_et*(el + (pl*ul - pl_bar*ul_et)/lambda_l)

    rhor_et = 1.0_DOUBLE/(1.0_DOUBLE/rhor + (ur - u_et)/lambda_r)
    ur_et = u_et
    vr_et = vr
    wr_et = wr
    pr_bar = pr + lambda_r*(u_et - ur)

    sol_r_et(1) = rhor_et
    sol_r_et(2) = rhor_et*ur_et
    sol_r_et(3) = rhor_et*vr_et
    sol_r_et(4) = rhor_et*wr_et
    sol_r_et(5) = rhor_et*(er + (pr_bar*ur_et - pr*ur)/lambda_r)

    sl = ul - lambda_l/rhol
    sr = ur + lambda_r/rhor

    lr_flux(:, 1) = 0.5_DOUBLE*(fl + fr) - 0.5_DOUBLE* &
      (abs(sl)*(sol_l_et - sol_l) + &
      abs(u_et)*(sol_r_et - sol_l_et) + &
      abs(sr)*(sol_r - sol_r_et))

    lr_flux(:, 2) = -lr_flux(:, 1)
  end subroutine three_wave

  subroutine modified_three_wave(sol_l, sol_r, sol_w_l, sol_w_r, lr_flux, sl, sr)
    use global_data_module, only: gamma
    implicit none

    real(kind=DOUBLE), dimension(5), intent(in) :: sol_l, sol_r
    real(kind=DOUBLE), dimension(5), intent(in) :: sol_w_l, sol_w_r
    real(kind=DOUBLE), dimension(5, 2), intent(inout) :: lr_flux
    real(kind=DOUBLE), intent(inout) :: sl, sr

    real(kind=DOUBLE) :: rhol, rhor, ul, vl, wl, ur, vr, wr
    real(kind=DOUBLE) :: pl, pr, rhol_et, rhor_et
    real(kind=DOUBLE) :: ul_et, vl_et, wl_et, ur_et, vr_et, wr_et, u_et
    real(kind=DOUBLE) :: el, er, al, ar
    real(kind=DOUBLE) :: piv_et, v_et, piw_et, w_et
    real(kind=DOUBLE) :: u_bar, pl_bar, pr_bar
    real(kind=DOUBLE) :: lambda_l, lambda_r
    real(kind=DOUBLE), dimension(5) :: fl, fr
    real(kind=DOUBLE), dimension(5) :: sol_l_et, sol_r_et

    rhol = sol_w_l(1)
    ul = sol_w_l(2)
    vl = sol_w_l(3)
    wl = sol_w_l(4)
    pl = sol_w_l(5)
    el = sol_l(5)/rhol
    al = sqrt(gamma*pl/rhol)

    rhor = sol_w_r(1)
    ur = sol_w_r(2)
    vr = sol_w_r(3)
    wr = sol_w_r(4)
    pr = sol_w_r(5)
    er = sol_r(5)/rhor
    ar = sqrt(gamma*pr/rhor)

    fl(:) = ul*sol_l(:) + (/0.0_DOUBLE, pl, 0.0_DOUBLE, 0.0_DOUBLE, pl*ul/)
    fr(:) = ur*sol_r(:) + (/0.0_DOUBLE, pr, 0.0_DOUBLE, 0.0_DOUBLE, pr*ur/)

    !Compute lambdas from two-point positivity conditions
    lambda_l = max(al*rhol, sqrt(rhol*max(0.0_DOUBLE, pr - pl)), -rhol*(ur - ul))
    lambda_r = max(ar*rhor, sqrt(rhor*max(0.0_DOUBLE, pl - pr)), -rhor*(ur - ul))
    u_bar = (lambda_l*ul + lambda_r*ur - (pr - pl))/(lambda_r + lambda_l)

    u_et = u_bar

    piv_et = - lambda_l*lambda_r/(lambda_l+lambda_r) * (vr - vl)
    piw_et = - lambda_l*lambda_r/(lambda_l+lambda_r) * (wr - wl)

    v_et = (lambda_l * vl + lambda_r * vr)/(lambda_l+lambda_r)
    w_et = (lambda_l * wl + lambda_r * wr)/(lambda_l+lambda_r)

    rhol_et = 1.0_DOUBLE/(1.0_DOUBLE/rhol + (u_et - ul)/lambda_l)
    ul_et = u_et
    vl_et = v_et
    wl_et = w_et
    pl_bar = pl - lambda_l*(u_et - ul)

    sol_l_et(1) = rhol_et
    sol_l_et(2) = rhol_et*ul_et
    sol_l_et(3) = rhol_et*vl_et
    sol_l_et(4) = rhol_et*wl_et
    sol_l_et(5) = rhol_et*(el + (pl*ul - pl_bar*ul_et &
      - piv_et*v_et - piw_et*w_et)/lambda_l)

    rhor_et = 1.0_DOUBLE/(1.0_DOUBLE/rhor + (ur - u_et)/lambda_r)
    ur_et = u_et
    vr_et = v_et
    wr_et = w_et
    pr_bar = pr + lambda_r*(u_et - ur)

    sol_r_et(1) = rhor_et
    sol_r_et(2) = rhor_et*ur_et
    sol_r_et(3) = rhor_et*vr_et
    sol_r_et(4) = rhor_et*wr_et
    sol_r_et(5) = rhor_et*(er + (pr_bar*ur_et - pr*ur&
      + piv_et*v_et + piw_et*w_et)/lambda_r)

    sl = ul - lambda_l/rhol
    sr = ur + lambda_r/rhor

    lr_flux(:, 1) = 0.5_DOUBLE*(fl + fr) - 0.5_DOUBLE* &
      (abs(sl)*(sol_l_et - sol_l) + &
      abs(u_et)*(sol_r_et - sol_l_et) + &
      abs(sr)*(sol_r - sol_r_et))

    lr_flux(:, 2) = -lr_flux(:, 1)
  end subroutine modified_three_wave

  subroutine two_wave(sol_l, sol_r, sol_w_l, sol_w_r, lr_flux, sl, sr)
    use global_data_module, only: gamma
    implicit none

    real(kind=DOUBLE), dimension(5), intent(in) :: sol_l, sol_r
    real(kind=DOUBLE), dimension(5), intent(in) :: sol_w_l, sol_w_r
    real(kind=DOUBLE), dimension(5, 2), intent(inout) :: lr_flux
    real(kind=DOUBLE), intent(inout) :: sl, sr

    real(kind=DOUBLE) :: rhol, rhor, ul, vl, wl, ur, vr, wr
    real(kind=DOUBLE) :: pl, pr, lambda_l, lambda_r
    real(kind=DOUBLE) :: al, ar
    real(kind=DOUBLE), dimension(5) :: fl, fr

    rhol = sol_w_l(1)
    ul = sol_w_l(2)
    vl = sol_w_l(3)
    wl = sol_w_l(4)
    pl = sol_w_l(5)
    al = sqrt(gamma*pl/rhol)

    rhor = sol_w_r(1)
    ur = sol_w_r(2)
    vr = sol_w_r(3)
    wr = sol_w_r(4)
    pr = sol_w_r(5)
    ar = sqrt(gamma*pr/rhor)

    fl(:) = ul*sol_l(:) + (/0.0_DOUBLE, pl, 0.0_DOUBLE, 0.0_DOUBLE, pl*ul/)
    fr(:) = ur*sol_r(:) + (/0.0_DOUBLE, pr, 0.0_DOUBLE, 0.0_DOUBLE, pr*ur/)

    lambda_l = max(al*rhol, sqrt(rhol*max(0.0_DOUBLE, pr - pl)), -rhol*(ur - ul))
    lambda_r = max(ar*rhor, sqrt(rhor*max(0.0_DOUBLE, pl - pr)), -rhor*(ur - ul))

    sl = ul - lambda_l/rhol
    sr = ur + lambda_r/rhor

    lr_flux(:, 1) = (sr*fl - sl*fr + sr*sl*(sol_r - sol_l))/(sr - sl)
    lr_flux(:, 2) = -lr_flux(:, 1)
  end subroutine two_wave

  subroutine modified_two_wave(sol_l, sol_r, sol_w_l, sol_w_r, lr_flux, sl, sr)
    use global_data_module, only: gamma
    implicit none

    real(kind=DOUBLE), dimension(5), intent(in) :: sol_l, sol_r
    real(kind=DOUBLE), dimension(5), intent(in) :: sol_w_l, sol_w_r
    real(kind=DOUBLE), dimension(5, 2), intent(inout) :: lr_flux
    real(kind=DOUBLE), intent(inout) :: sl, sr

    real(kind=DOUBLE) :: rhol, rhor, ul, vl, wl, ur, vr, wr, pl, pr
    real(kind=DOUBLE) :: el, er, al, ar, taul, taur
    real(kind=DOUBLE) :: lambda_l, lambda_r
    real(kind=DOUBLE), dimension(5) :: fl, fr
    real(kind=DOUBLE) :: u_et, u_bar, rho_et, e_et, p_bar, pv_bar, tau_et
    real(kind=DOUBLE), dimension(5) :: sol_l_et, sol_r_et

    rhol = sol_w_l(1)
    taul = 1.0_DOUBLE/rhol
    ul = sol_w_l(2)
    vl = sol_w_l(3)
    wl = sol_w_l(4)
    pl = sol_w_l(5)
    el = sol_l(5)/rhol
    al = sqrt(gamma*pl/rhol)

    rhor = sol_w_r(1)
    taur = 1.0_DOUBLE/rhor
    ur = sol_w_r(2)
    vr = sol_w_r(3)
    wr = sol_w_r(4)
    pr = sol_w_r(5)
    er = sol_r(5)/rhor
    ar = sqrt(gamma*pr/rhor)

    fl(:) = ul*sol_l(:) + (/0.0_DOUBLE, pl, 0.0_DOUBLE, 0.0_DOUBLE, pl*ul/)
    fr(:) = ur*sol_r(:) + (/0.0_DOUBLE, pr, 0.0_DOUBLE, 0.0_DOUBLE, pr*ur/)

    !Compute lambdas from two-point positivity conditions
    lambda_l = max(al*rhol, sqrt(rhol*max(0.0_DOUBLE, pr - pl)), -rhol*(ur - ul))
    lambda_r = max(ar*rhor, sqrt(rhor*max(0.0_DOUBLE, pl - pr)), -rhor*(ur - ul))

    u_bar = (lambda_l*ur + lambda_r*ul + lambda_l*lambda_r*(taur - taul))/(lambda_r+lambda_l)
    u_et = (lambda_l*ul + lambda_r*ur - (pr - pl))/(lambda_r+lambda_l)
    tau_et = (lambda_l*taul+lambda_r*taur + ur-ul)/(lambda_l+lambda_r)
    rho_et = 1.0_DOUBLE/tau_et
    e_et = (lambda_l*el + lambda_r*er - (pr*ur - pl*ul))/(lambda_l+lambda_r)

    sol_l_et(1) = rho_et
    sol_l_et(2) = rho_et*u_et
    sol_l_et(3) = rho_et*vl
    sol_l_et(4) = rho_et*wl
    sol_l_et(5) = rho_et*e_et

    sol_r_et(1) = rho_et
    sol_r_et(2) = rho_et*u_et
    sol_r_et(3) = rho_et*vr
    sol_r_et(4) = rho_et*wr
    sol_r_et(5) = rho_et*e_et

    sl = ul - lambda_l/rhol
    sr = ur + lambda_r/rhor

    lr_flux(:, 1) = 0.5_DOUBLE*(fl + fr) - 0.5_DOUBLE* &
      (abs(sl)*(sol_l_et - sol_l) + &
      abs(u_bar)*(sol_r_et - sol_l_et) + &
      abs(sr)*(sol_r - sol_r_et))

    lr_flux(:, 2) = -lr_flux(:, 1)
  end subroutine modified_two_wave

  subroutine multi_point(sol_l, sol_r, sol_w_l, sol_w_r, lr_flux, u_nodal, &
      lambda_l, lambda_r, sl, sr)
    implicit none

    real(kind=DOUBLE), intent(in) :: u_nodal
    real(kind=DOUBLE), dimension(5), intent(in) :: sol_l, sol_r
    real(kind=DOUBLE), dimension(5), intent(in) :: sol_w_l, sol_w_r
    real(kind=DOUBLE), intent(inout) :: lambda_l, lambda_r
    real(kind=DOUBLE), dimension(5, 2), intent(inout) :: lr_flux
    real(kind=DOUBLE), intent(inout) :: sl, sr

    real(kind=DOUBLE) :: rhol, rhor, ul, vl, wl, ur, vr, wr
    real(kind=DOUBLE) :: pl, pr, rhol_et, rhor_et
    real(kind=DOUBLE) :: ul_et, vl_et, wl_et, ur_et, vr_et, wr_et, u_et
    real(kind=DOUBLE) :: el, er
    real(kind=DOUBLE) :: u_bar, pl_bar, pr_bar
    real(kind=DOUBLE), dimension(5) :: fl, fr
    real(kind=DOUBLE), dimension(5) :: sol_l_et, sol_r_et

    rhol = sol_w_l(1)
    ul = sol_w_l(2)
    vl = sol_w_l(3)
    wl = sol_w_l(4)
    pl = sol_w_l(5)
    el = sol_l(5)/rhol

    rhor = sol_w_r(1)
    ur = sol_w_r(2)
    vr = sol_w_r(3)
    wr = sol_w_r(4)
    pr = sol_w_r(5)
    er = sol_r(5)/rhor

    fl(:) = ul*sol_l(:) + (/0.0_DOUBLE, pl, 0.0_DOUBLE, 0.0_DOUBLE, pl*ul/)
    fr(:) = ur*sol_r(:) + (/0.0_DOUBLE, pr, 0.0_DOUBLE, 0.0_DOUBLE, pr*ur/)

    u_bar = (lambda_l*ul + lambda_r*ur)/(lambda_l + lambda_r) &
      - (pr - pl)/(lambda_r + lambda_l)

    u_et = u_nodal

    rhol_et = 1.0_DOUBLE/(1.0_DOUBLE/rhol + (u_et - ul)/lambda_l)
    ul_et = u_et
    vl_et = vl
    wl_et = wl
    pl_bar = pl - lambda_l*(u_et - ul)

    sol_l_et(1) = rhol_et
    sol_l_et(2) = rhol_et*ul_et
    sol_l_et(3) = rhol_et*vl_et
    sol_l_et(4) = rhol_et*wl_et
    sol_l_et(5) = rhol_et*(el + (pl*ul - pl_bar*u_et)/lambda_l)

    rhor_et = 1.0_DOUBLE/(1.0_DOUBLE/rhor + (ur - u_et)/lambda_r)
    ur_et = u_et
    vr_et = vr
    wr_et = wr
    pr_bar = pr + lambda_r*(u_et - ur)

    sol_r_et(1) = rhor_et
    sol_r_et(2) = rhor_et*ur_et
    sol_r_et(3) = rhor_et*vr_et
    sol_r_et(4) = rhor_et*wr_et
    sol_r_et(5) = rhor_et*(er + (pr_bar*u_et - pr*ur)/lambda_r)

    if (rhol_et < 0.0_DOUBLE .or. rhor_et < 0.0_DOUBLE) then
      print *, "Negative specific volume MPCC !", rhol_et, rhor_et
      print*, rhol, rhor
      print*, lambda_l, lambda_r
      print*, ul_et, vl_et, wl_et
      print*, ur_et, vr_et, wr_et
      print*,"R"
      print*, 1.0_DOUBLE/rhor, (ur - u_et)/lambda_r
      error stop
    end if

    sl = ul - lambda_l/rhol
    sr = ur + lambda_r/rhor

    lr_flux(:, 1) = 0.5_DOUBLE*(fl + fr) - 0.5_DOUBLE* &
      (abs(sl)*(sol_l_et - sol_l) + &
      abs(u_et)*(sol_r_et - sol_l_et) + &
      abs(sr)*(sol_r - sol_r_et)) &
      - 0.5_DOUBLE*(pr_bar - pl_bar)* &
      (/0.0_DOUBLE, 1.0_DOUBLE, 0.0_DOUBLE, 0.0_DOUBLE, u_et/)

    lr_flux(:, 2) = 0.5_DOUBLE*(fl + fr) - 0.5_DOUBLE* &
      (abs(sl)*(sol_l_et - sol_l) + &
      abs(u_et)*(sol_r_et - sol_l_et) + &
      abs(sr)*(sol_r - sol_r_et)) &
      + 0.5_DOUBLE*(pr_bar - pl_bar)* &
      (/0.0_DOUBLE, 1.0_DOUBLE, 0.0_DOUBLE, 0.0_DOUBLE, u_et/)

    lr_flux(:, 2) = -lr_flux(:, 2)
  end subroutine multi_point

  subroutine compute_dt_a_priori(mesh, t, dt, t_max, cfl, sol, num_procs)
    use global_data_module, only: gamma, dt_is_fixed, fixed_dt
    use mpi
    implicit none

    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), intent(in) :: t, t_max
    real(kind=DOUBLE), intent(inout) :: dt
    real(kind=DOUBLE), intent(in) :: cfl
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(in) :: sol
    integer(kind=ENTIER), intent(in) :: num_procs

    integer(kind=ENTIER) :: i, j, id_face, mpi_ierr
    real(kind=DOUBLE) :: eigen_value, sum_eig
    real(kind=DOUBLE), dimension(5) :: sol_w

    dt = 1e16
    !$OMP PARALLEL DO REDUCTION(min:dt) &
    !$OMP& PRIVATE(j, id_face, eigen_value, sol_w, sum_eig)
    do i = 1, mesh%n_elems
      if (.not. mesh%elem(i)%is_ghost) then
        sum_eig = 0.0_DOUBLE
        sol_w = conserv_to_primit(sol(:, i))
        do j=1, mesh%elem(i)%n_faces
          id_face = mesh%elem(i)%face(j)
          eigen_value = abs(dot_product(sol_w(2:4), mesh%face(id_face)%norm)) &
            + sqrt(gamma*sol_w(5)/sol_w(1))
          !sum_eig = sum_eig + eigen_value*mesh%face(id_face)%area
          dt = min(dt, cfl*mesh%elem(i)%volume/(eigen_value*mesh%face(id_face)%area))
        end do
        !dt = min(dt, cfl*mesh%elem(i)%volume/sum_eig)
      end if
    end do
    dt = cfl * dt

    if( dt_is_fixed ) dt = fixed_dt

    if( t + dt > t_max ) dt = t_max - t
    if (num_procs > 1) then
      call MPI_ALLREDUCE(MPI_IN_PLACE, dt, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD, mpi_ierr)
    end if
  end subroutine compute_dt_a_priori

  subroutine compute_dt(mesh, dt, cfl, sum_lambda, num_procs)
    use mpi
    implicit none

    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), dimension(mesh%n_elems), intent(in) :: sum_lambda
    real(kind=DOUBLE), intent(inout) :: dt
    real(kind=DOUBLE), intent(in) :: cfl
    integer(kind=ENTIER), intent(in) :: num_procs

    integer(kind=ENTIER) :: i, mpi_ierr

    dt = 1e10
    !$OMP PARALLEL DO REDUCTION(min:dt)
    do i = 1, mesh%n_elems
      if (.not. mesh%elem(i)%is_ghost) then
        dt = min(dt, cfl*mesh%elem(i)%volume/sum_lambda(i))
      end if
    end do

    if (num_procs > 1) then
      call MPI_ALLREDUCE(MPI_IN_PLACE, dt, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD, mpi_ierr)
    end if
  end subroutine compute_dt

  subroutine init_sol(mesh, sol, me, num_procs, mpi_send_recv)
    use global_data_module
    use mpi
    use subfv_mpi_module
    implicit none

    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(inout) :: sol
    integer(kind=ENTIER), intent(in) :: me, num_procs
    type(mpi_send_recv_type) :: mpi_send_recv

    integer(kind=ENTIER) :: i
    real(kind=DOUBLE) :: tot_vol, r
    real(kind=DOUBLE) :: H_helm, M_helm, r_helm, delta_helm
    real(kind=DOUBLE), dimension(5) :: w
    real(kind=DOUBLE) :: w_gresho, m_gresho, p_gresho
    real(kind=DOUBLE), dimension(3) :: coord2

    if (init_uniform) then
      do i = 1, mesh%n_elems
        sol(:, i) = primit_to_conserv(sol_uniform)
      end do
    else if (init_grad_p) then
      do i = 1, mesh%n_elems
        sol_uniform = &
          (/1.0_DOUBLE, 0.0_DOUBLE, 0.0_DOUBLE, 0.0_DOUBLE, 100.0_DOUBLE + mesh%elem(i)%coord(1)/)
        sol(:, i) = primit_to_conserv(sol_uniform)
      end do
    else if (init_sod) then
      do i = 1, mesh%n_elems
        w = sol_primit_sod_angle(mesh%elem(i)%coord, &
          x_sod, sod_angle, sol_left, sol_right)
        sol(:, i) = primit_to_conserv(w)
      end do
    else if (init_double_mach) then
      sol_left = (/1.4_DOUBLE, 0.0_DOUBLE, 0.0_DOUBLE, 0.0_DOUBLE, 1.0_DOUBLE/)
      sol_right = (/8.0_DOUBLE, 7.145_DOUBLE, -4.125_DOUBLE, 0.0_DOUBLE, 116.5_DOUBLE/)
      do i = 1, mesh%n_elems
        if (mesh%elem(i)%coord(2) < 1.732_DOUBLE*(mesh%elem(i)%coord(1) - 0.1667_DOUBLE)) then
          sol(:, i) = primit_to_conserv(sol_left)
        else
          sol(:, i) = primit_to_conserv(sol_right)
        end if
      end do
    else if (init_spherical_sod) then
      do i = 1, mesh%n_elems
        if (norm2(mesh%elem(i)%coord) < r_spherical_sod) then
          sol(:, i) = primit_to_conserv(spherical_sod_sol_left)
        else
          sol(:, i) = primit_to_conserv(spherical_sod_sol_right)
        end if
      end do
    else if (init_ellig) then
      do i = 1, mesh%n_elems
        if (mesh%elem(i)%coord(1) < 50.0_DOUBLE) then
          if (mesh%elem(i)%coord(2) > 19.75_DOUBLE .and. &
            mesh%elem(i)%coord(2) < 20.25_DOUBLE) then
            sol_left = (/1.0_double, 0.0_double, &
              0.0_double, 0.0_double, 5.0_double/63.0_double/)
            sol(:, i) = primit_to_conserv(sol_left)
          else
            sol_left = (/1.0_DOUBLE, 1.0_DOUBLE, &
              0.0_DOUBLE, 0.0_DOUBLE, 5.0_DOUBLE/63.0_DOUBLE/)
            sol(:, i) = primit_to_conserv(sol_left)
          end if
        else
          sol_right = (/27.0_DOUBLE/7.0_DOUBLE, 7.0_DOUBLE/27.0_DOUBLE, &
            0.0_DOUBLE, 0.0_DOUBLE, 155.0_DOUBLE/189.0_DOUBLE/)
          sol(:, i) = primit_to_conserv(sol_right)
        end if
      end do

    else if (init_sedov) then
      tot_vol = 0.0_DOUBLE
      do i = 1, mesh%n_elems
        ! if (sqrt(mesh%elem(i)%coord(1)**2 + mesh%elem(i)%coord(2)**2) < r_sedov) then
        if (norm2(mesh%elem(i)%coord) < r_sedov) then
          tot_vol = tot_vol + mesh%elem(i)%volume
        end if
      end do

      if (tot_vol < 1e-10_DOUBLE) then
        print *, "[-] r_sedov too small, tot_vol_init:", tot_vol
        error stop
      end if

      do i = 1, mesh%n_elems
        if (norm2(mesh%elem(i)%coord) < r_sedov) then
          sol(1, i) = 1.0_DOUBLE
          sol(2, i) = 0.0_DOUBLE
          sol(3, i) = 0.0_DOUBLE
          sol(4, i) = 0.0_DOUBLE
          ! sol(5,i) = 0.311357_DOUBLE / tot_vol
          ! sol(5, i) = 0.244816_DOUBLE/tot_vol
          sol(5, i) = 0.851072_DOUBLE/tot_vol
        else
          sol(1, i) = 1.0_DOUBLE
          sol(2, i) = 0.0_DOUBLE
          sol(3, i) = 0.0_DOUBLE
          sol(4, i) = 0.0_DOUBLE
          sol(5, i) = 1e-12_DOUBLE/(gamma - 1.0_DOUBLE)
        end if
      end do
    else if (init_noh) then
      gamma = 5.0_DOUBLE/3.0_DOUBLE
      do i = 1, mesh%n_elems
        w(1) = 1.0_DOUBLE
        w(2) = -mesh%elem(i)%coord(1)/sqrt(mesh%elem(i)%coord(1)**2 + mesh%elem(i)%coord(2)**2)
        w(3) = -mesh%elem(i)%coord(2)/sqrt(mesh%elem(i)%coord(1)**2 + mesh%elem(i)%coord(2)**2)
        w(4) = 0.0_DOUBLE
        w(5) = 1e-6_DOUBLE

        sol(:, i) = primit_to_conserv(w)
      end do
    else if (init_isentropic_vortex) then
      do i = 1, mesh%n_elems
        call sol_isentropic_vortex(mesh%elem(i)%coord, w, 0.0_DOUBLE)
        sol(:, i) = primit_to_conserv(w)
      end do
    else if (init_gresho) then
      w_gresho = 0.2_DOUBLE
      m_gresho = 1e-3
      p_gresho = 1.0_DOUBLE/(gamma*m_gresho**2)
      do i = 1, mesh%n_elems
        coord2(:) = mesh%elem(i)%coord(:) - (/0.5_DOUBLE, 0.5_DOUBLE, 0.5_DOUBLE/)
        r = norm2(coord2(:2))

        if (r < w_gresho) then
          w(1) = 1.0_DOUBLE
          w(2) = -5.0_DOUBLE*coord2(2)
          w(3) = 5.0_DOUBLE*coord2(1)
          w(4) = 0.0_DOUBLE
          w(5) = p_gresho + 12.5*r**2
        else if (r > 1.0_DOUBLE*w_gresho .and. r < 2.0_DOUBLE*w_gresho) then
          w(1) = 1.0_DOUBLE
          w(2) = (5.0_DOUBLE - 2.0_DOUBLE/r)*coord2(2)
          w(3) = (-5.0_DOUBLE + 2.0_DOUBLE/r)*coord2(1)
          w(4) = 0.0_DOUBLE
          w(5) = p_gresho + 12.5_DOUBLE*r**2 + &
            4 - 20.0_DOUBLE*r + 4.0_DOUBLE*log(5.0_DOUBLE*r)
        else
          w(1) = 1.0_DOUBLE
          w(2) = 0.0_DOUBLE
          w(3) = 0.0_DOUBLE
          w(4) = 0.0_DOUBLE
          w(5) = p_gresho - 2.0_DOUBLE + 4.0_DOUBLE*log(2.0_DOUBLE)
        end if

        sol(:, i) = primit_to_conserv(w)
      end do
    else if (init_kelvin) then
      M_Helm = 1e-2_DOUBLE
      r_helm = 1e-3_DOUBLE
      delta_helm = 0.1_DOUBLE
      t_max = 0.8_DOUBLE/M_helm
      print *, "T_max modified to", t_max
      do i = 1, mesh%n_elems
        H_helm = kelvin_helmholtz(mesh%elem(i)%coord(2))
        w(1) = gamma + r_helm*H_helm
        w(2) = M_helm*H_helm
        w(3) = delta_helm*M_helm*sin(2.0_DOUBLE*pi*mesh%elem(i)%coord(1))
        w(4) = 0.0_DOUBLE
        w(5) = 1.0_DOUBLE
        sol(:, i) = primit_to_conserv(w)
      end do
    else if (init_restart) then
      call restart_from_vtu_file(mesh, sol, restart_file, me, num_procs)
    else if (init_double_vortex_rings) then
      sol(1, :) = 1.0_DOUBLE
      sol(2:4, :) = 0.0_DOUBLE
      sol(5, :) = 1.0_DOUBLE
      call put_vortex_ring(mesh, sol, (/0.0_DOUBLE, 0.0_DOUBLE, 0.5_DOUBLE/), &
        (/0.0_DOUBLE, 0.0_DOUBLE, 0.0_DOUBLE/))
      do i = 1, mesh%n_elems
        sol(:, i) = primit_to_conserv(sol(:,i))
      end do
    else if (init_vortex_ring_generator) then
      do i = 1, mesh%n_elems
        if( mesh%elem(i)%coord(3) < 0.25 .or. mesh%elem(i)%coord(3) > 3.75) then
          sol(:, i) = primit_to_conserv((/1.0_double, 0.0_double, &
            0.0_double, 0.0_double, 1.1_double/))
        else
          sol(:, i) = primit_to_conserv((/1.0_double, 0.0_double, &
            0.0_double, 0.0_double, 1.0_double/))
        end if
      end do
    end if

    if (num_procs > 1) call mpi_memory_exchange(mpi_send_recv, &
      mesh%n_elems, 5, sol)
  end subroutine init_sol

  subroutine compute_euler_flux_mcfl(mesh, sol, euler_flux_sum, grad, t, dt, &
      divF, u_nodal, sum_lambda, alpha_mcfl)
    use global_data_module, &
      only: bc_style, mach_threshold, scheme, second_order, &
      exclude_bound_vert, low_dissip_switch, timedisc, local_time_step
    use subfv_linear_solver_module
    use mpi
    implicit none

    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(in) :: sol
    real(kind=DOUBLE), dimension(:, :, :), intent(in) :: grad
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(inout) :: euler_flux_sum
    real(kind=DOUBLE), dimension(:, :), intent(in) :: divF
    real(kind=DOUBLE), dimension(3, mesh%n_vert), intent(inout) :: u_nodal
    real(kind=DOUBLE), dimension(mesh%n_vert), intent(inout) :: alpha_mcfl
    real(kind=DOUBLE), dimension(mesh%n_elems), intent(inout) :: sum_lambda
    real(kind=DOUBLE), intent(in) :: t, dt

    integer(kind=ENTIER) :: j, k, idse, ide, id_vert, nsfn, nsen
    real(kind=DOUBLE) :: alpha
    real(kind=DOUBLE), dimension(3) :: u_vert
    real(kind=DOUBLE), dimension(:), allocatable :: sum_lambda_vert_o1
    real(kind=DOUBLE), dimension(:), allocatable :: sum_lambda_vert_o2
    real(kind=DOUBLE), dimension(:, :), allocatable :: flux_sum_vert_o1
    real(kind=DOUBLE), dimension(:, :), allocatable :: flux_sum_vert_o2

    euler_flux_sum = 0.0_DOUBLE
    u_nodal = 0.0_DOUBLE
    sum_lambda = 0.0_DOUBLE

    !$OMP PARALLEL DO PRIVATE(j, k, idse, ide, &
    !$OMP id_vert, nsfn, nsen, alpha, u_vert, sum_lambda_vert_o1, &
    !$OMP sum_lambda_vert_o2, flux_sum_vert_o1, flux_sum_vert_o2)
    do id_vert = 1, mesh%n_vert
      if (.not. mesh%vert(id_vert)%is_ghost) then
        nsfn = mesh%vert(id_vert)%n_sub_faces_neigh
        nsen = mesh%vert(id_vert)%n_sub_elems_neigh
        if( second_order ) then
          allocate(sum_lambda_vert_o1(nsen))
          sum_lambda_vert_o1 = 0.0_DOUBLE
          allocate(flux_sum_vert_o1(5, nsen))
          flux_sum_vert_o1 = 0.0_DOUBLE

          allocate(sum_lambda_vert_o2(nsen))
          sum_lambda_vert_o2 = 0.0_DOUBLE
          allocate(flux_sum_vert_o2(5, nsen))
          flux_sum_vert_o2 = 0.0_DOUBLE

          call compute_rhs_around_vert(mesh, sol, grad, divF, &
            t, dt, nsfn, nsen, u_vert, sum_lambda_vert_o1, &
            flux_sum_vert_o1, .false., id_vert)

          call compute_rhs_around_vert(mesh, sol, grad, divF, &
            t, dt, nsfn, nsen, u_vert, sum_lambda_vert_o2, &
            flux_sum_vert_o2, .true., id_vert)

          !Compute limiting coeff alpha
          alpha = 1.0_DOUBLE
          do j=1, mesh%vert(id_vert)%n_sub_elems_neigh
            idse = mesh%vert(id_vert)%sub_elem_neigh(j)
            ide = mesh%sub_elem(idse)%mesh_elem
            alpha = min(alpha, max_omega_for_positivity(sol(:, ide), &
              -dt/mesh%sub_elem(idse)%volume&
              *(flux_sum_vert_o2(:, j) - flux_sum_vert_o1(:, j))))
          end do

          alpha_mcfl(id_vert) = alpha
          u_nodal(:, id_vert) = u_vert
          do j=1, mesh%vert(id_vert)%n_sub_elems_neigh
            idse = mesh%vert(id_vert)%sub_elem_neigh(j)
            ide = mesh%sub_elem(idse)%mesh_elem
            do k=1, 5
              !$OMP ATOMIC UPDATE
              euler_flux_sum(k, ide) = euler_flux_sum(k, ide) &
                + flux_sum_vert_o1(k, j) &
                + alpha * (flux_sum_vert_o2(k, j) - flux_sum_vert_o1(k, j))
            end do
            sum_lambda(ide) = sum_lambda(ide) + sum_lambda_vert_o1(j)
          end do
          deallocate(sum_lambda_vert_o1, flux_sum_vert_o1)
          deallocate(sum_lambda_vert_o2, flux_sum_vert_o2)
        else !First order
          allocate(sum_lambda_vert_o1(nsen))
          sum_lambda_vert_o1 = 0.0_DOUBLE
          allocate(flux_sum_vert_o1(5, nsen))
          flux_sum_vert_o1 = 0.0_DOUBLE

          call compute_rhs_around_vert(mesh, sol, grad, divF, &
            t, dt, nsfn, nsen, u_vert, sum_lambda_vert_o1, &
            flux_sum_vert_o1, .false., id_vert)

          u_nodal(:, id_vert) = u_vert
          do j=1, mesh%vert(id_vert)%n_sub_elems_neigh
            idse = mesh%vert(id_vert)%sub_elem_neigh(j)
            ide = mesh%sub_elem(idse)%mesh_elem
            do k=1, 5
              !$OMP ATOMIC UPDATE
              euler_flux_sum(k, ide) = euler_flux_sum(k, ide) &
                + flux_sum_vert_o1(k, j)
            end do
            sum_lambda(ide) = sum_lambda(ide) + sum_lambda_vert_o1(j)
          end do
          deallocate(sum_lambda_vert_o1, flux_sum_vert_o1)
        end if
      end if
    end do
  end subroutine compute_euler_flux_mcfl

  subroutine compute_rhs_around_vert(mesh, sol, grad, divF, &
      t, dt, nsfn, nsen, u_vert, sum_lambda_vert, flux_sum_vert, second_order, id_vert)
    use global_data_module, &
      only: bc_style, scheme, &
      exclude_bound_vert, timedisc, local_time_step
    use subfv_linear_solver_module
    use mpi
    implicit none

    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(in) :: sol
    real(kind=DOUBLE), dimension(:, :, :), intent(in) :: grad
    real(kind=DOUBLE), dimension(:, :), intent(in) :: divF
    real(kind=DOUBLE), intent(in) :: t, dt
    integer(kind=ENTIER), intent(in) :: nsfn, nsen
    real(kind=DOUBLE), dimension(3), intent(inout) :: u_vert
    real(kind=DOUBLE), dimension(nsen), intent(inout) :: sum_lambda_vert
    real(kind=DOUBLE), dimension(5, nsen), intent(inout) :: flux_sum_vert
    logical, intent(in) :: second_order
    integer(kind=ENTIER), intent(in) :: id_vert

    integer(kind=ENTIER) :: k, j, id_sub_face, id_face
    integer(kind=ENTIER) :: le, re, lse, rse, lse_loc, rse_loc
    real(kind=DOUBLE) :: p_bound, omega, p_nodal, sl, sr
    real(kind=DOUBLE), dimension(5, 2, nsfn) :: sol_lr, sol_w_lr
    real(kind=DOUBLE), dimension(5, 2, nsfn) :: lr_flux, lr_flux_tw
    real(kind=DOUBLE), dimension(2, nsfn) :: lambda
    real(kind=DOUBLE), dimension(nsfn) :: warea
    real(kind=DOUBLE), dimension(nsfn) :: u_bars, p_bars

    flux_sum_vert = 0.0_DOUBLE
    u_vert = 0.0_DOUBLE
    sum_lambda_vert = 0.0_DOUBLE
    warea = 0.0_DOUBLE

    do j = 1, mesh%vert(id_vert)%n_sub_faces_neigh
      id_sub_face = mesh%vert(id_vert)%sub_face_neigh(j)
      id_face = mesh%sub_face(id_sub_face)%mesh_face

      le = mesh%sub_face(id_sub_face)%left_elem_neigh
      re = mesh%sub_face(id_sub_face)%right_elem_neigh

      warea(j) = mesh%sub_face(id_sub_face)%area

      sol_lr(:, 1, j) = sol(:, le)
      sol_w_lr(:, 1, j) = conserv_to_primit(sol_lr(:, 1, j))

      if (second_order) then
        call reconstruct_with_grad(sol_w_lr(:, 1, j), grad(:, :, le), &
          mesh%vert(id_vert)%coord - mesh%elem(le)%coord)
        sol_lr(:, 1, j) = primit_to_conserv(sol_w_lr(:, 1, j))
        !Hancock : progress hafl a time step
        if(timedisc == "euler" .and. .not. local_time_step) then
          omega = min(1.0_DOUBLE, 0.9_DOUBLE*max_omega_for_positivity(sol_lr(:, 1, j), - 0.5_DOUBLE*dt*divF(:, le)))
          sol_lr(:, 1, j) = sol_lr(:, 1, j) - 0.5_DOUBLE*omega*dt*divF(:, le)
          sol_w_lr(:, 1, j) = conserv_to_primit(sol_lr(:, 1, j))
        end if
      end if

      if( re > 0 ) then
        sol_lr(:, 2, j) = sol(:, re)
        sol_w_lr(:, 2, j) = conserv_to_primit(sol_lr(:, 2, j))

        if (second_order) then
          call reconstruct_with_grad(sol_w_lr(:, 2, j), grad(:, :, re), &
            mesh%vert(id_vert)%coord - mesh%elem(re)%coord)
          sol_lr(:, 2, j) = primit_to_conserv(sol_w_lr(:, 2, j))
          !Hancock : progress hafl a time step
          if(timedisc == "euler" .and. .not. local_time_step) then
            omega = min(1.0_DOUBLE, &
              max_omega_for_positivity(sol_lr(:, 2, j), - 0.5_DOUBLE*dt*divF(:, re)))
            sol_lr(:, 2, j) = sol_lr(:, 2, j) - 0.5_DOUBLE*omega*dt*divF(:, re)
            sol_w_lr(:, 2, j) = conserv_to_primit(sol_lr(:, 2, j))
          end if
        end if
      else
        call compute_right_state(mesh, id_face, le, re, sol_lr(:, 1, j), sol_lr(:, 2, j), t)
        sol_w_lr(:, 2, j) = conserv_to_primit(sol_lr(:, 2, j))
      end if

      call base_change(sol_lr(2:4, 1, j), sol_lr(2:4, 2, j), mesh%face(id_face)%norm, .TRUE.)
      call base_change(sol_w_lr(2:4, 1, j), sol_w_lr(2:4, 2, j), mesh%face(id_face)%norm, .TRUE.)
    end do

    warea = minval(warea)/warea

    !Solve nodal system if multi_point scheme
    if (scheme == "multi_point" ) then
      lambda(:, :) = 0.0_DOUBLE
      call compute_lambdas_and_solve_nodal_velocity(mesh, id_vert, &
        sol_w_lr, lambda, u_bars, u_vert, p_bound)
    else if (scheme == "multi_point_iso" ) then
      lambda(:, :) = 0.0_DOUBLE
      call compute_lambdas_and_solve_nodal_velocity_iso(mesh, id_vert, &
        sol_w_lr, lambda, u_bars, u_vert, p_bound)
    else if(scheme =="multi_point_pressure") then
      lambda(:, :) = 0.0_DOUBLE
      call compute_lambdas_and_solve_nodal_pressure(mesh, id_vert, &
        sol_w_lr, lambda, p_bars, p_nodal)
    end if

    !!Compute flux across each sub_face
    do j = 1, mesh%vert(id_vert)%n_sub_faces_neigh
      id_sub_face = mesh%vert(id_vert)%sub_face_neigh(j)
      id_face = mesh%sub_face(id_sub_face)%mesh_face
      le = mesh%sub_face(id_sub_face)%left_elem_neigh
      re = mesh%sub_face(id_sub_face)%right_elem_neigh
      lse = mesh%sub_face(id_sub_face)%left_sub_elem_neigh
      lse_loc = mesh%sub_elem(lse)%id_loc_around_node
      rse = mesh%sub_face(id_sub_face)%right_sub_elem_neigh
      if( rse > 0 ) then
        rse_loc = mesh%sub_elem(rse)%id_loc_around_node
      end if

      if (scheme == "multi_point") then
        if( mesh%vert(id_vert)%is_bound .and. exclude_bound_vert) then
          call two_wave(sol_lr(:, 1, j), sol_lr(:, 2, j), &
            sol_w_lr(:, 1, j), sol_w_lr(:, 2, j), lr_flux(:, :, j), sl, sr)
        else
          if (is_wall(re) .and. bc_style == 1) then
            call multi_point(sol_lr(:, 1, j), sol_lr(:, 2, j), &
              sol_w_lr(:, 1, j), sol_w_lr(:, 2, j), &
              lr_flux(:, :, j), 0.0_DOUBLE, &
              lambda(1, j), lambda(2, j), sl, sr)
          else
            call multi_point(sol_lr(:, 1, j), sol_lr(:, 2, j), &
              sol_w_lr(:, 1, j), sol_w_lr(:, 2, j), &
              lr_flux(:, :, j), &
              dot_product(u_vert, mesh%face(id_face)%norm), &
              lambda(1, j), lambda(2, j), sl, sr)
          end if
        end if
      else if (scheme == "multi_point_iso") then
        call three_wave(sol_lr(:, 1, j), sol_lr(:, 2, j), &
          sol_w_lr(:, 1, j), sol_w_lr(:, 2, j), lr_flux_tw(:, :, j), sl, sr)
        if( mesh%vert(id_vert)%is_bound .and. exclude_bound_vert) then
          call two_wave(sol_lr(:, 1, j), sol_lr(:, 2, j), &
            sol_w_lr(:, 1, j), sol_w_lr(:, 2, j), lr_flux(:, :, j), sl, sr)
        else
          if (is_wall(re) .and. bc_style == 1) then
            call multi_point(sol_lr(:, 1, j), sol_lr(:, 2, j), &
              sol_w_lr(:, 1, j), sol_w_lr(:, 2, j), &
              lr_flux(:, :, j), 0.0_DOUBLE, &
              lambda(1, j), lambda(2, j), sl, sr)
          else
            call multi_point(sol_lr(:, 1, j), sol_lr(:, 2, j), &
              sol_w_lr(:, 1, j), sol_w_lr(:, 2, j), &
              lr_flux(:, :, j), &
              dot_product(u_vert, mesh%face(id_face)%norm), &
              lambda(1, j), lambda(2, j), sl, sr)
          end if
        end if
        lr_flux(:, :, j) = warea(j) * lr_flux(:, :, j) &
          + (1.0_DOUBLE - warea(j)) * lr_flux_tw(:, :, j)
      else if(scheme=="three_wave") then
        call three_wave(sol_lr(:, 1, j), sol_lr(:, 2, j), &
          sol_w_lr(:, 1, j), sol_w_lr(:, 2, j), lr_flux(:, :, j), sl, sr)
      else if(scheme == "two_wave") then
        call two_wave(sol_lr(:, 1, j), sol_lr(:, 2, j), &
          sol_w_lr(:, 1, j), sol_w_lr(:, 2, j), lr_flux(:, :, j), sl, sr)
      else if(scheme == "modified_three_wave") then
        call modified_three_wave(sol_lr(:, 1, j), sol_lr(:, 2, j), &
          sol_w_lr(:, 1, j), sol_w_lr(:, 2, j), lr_flux(:, :, j), sl, sr)
      else
        print*,"Unknown scheme"
        error stop
      end if

      !Rotate flux back
      call base_change(lr_flux(2:4, 1, j), lr_flux(2:4, 2, j), mesh%face(id_face)%norm, .FALSE.)

      !$OMP ATOMIC UPDATE
      sum_lambda_vert(lse_loc) = sum_lambda_vert(lse_loc) + mesh%sub_face(id_sub_face)%area*max(0.0_DOUBLE, -sl)
      do k=1, 5
        !$OMP ATOMIC UPDATE
        flux_sum_vert(k, lse_loc) = flux_sum_vert(k, lse_loc) + mesh%sub_face(id_sub_face)%area*lr_flux(k, 1, j)
      end do

      if (re > 0) then
        !$OMP ATOMIC UPDATE
        sum_lambda_vert(rse_loc) = sum_lambda_vert(rse_loc) + mesh%sub_face(id_sub_face)%area*max(0.0_DOUBLE, sr)
        do k=1, 5
          !$OMP ATOMIC UPDATE
          flux_sum_vert(k, rse_loc) = flux_sum_vert(k, rse_loc) + mesh%sub_face(id_sub_face)%area*lr_flux(k, 2, j)
        end do
      end if
    end do
  end subroutine compute_rhs_around_vert

  subroutine compute_right_state(mesh, id_face, le, re, sol_l, sol_r, t)
    use global_data_module, only: gamma, bc_type, bc_val, sol_left, sol_right, x_sod, pi, &
      sod_angle
    implicit none

    type(mesh_type), intent(in) :: mesh
    integer(kind=ENTIER), intent(in) :: id_face, le, re
    real(kind=DOUBLE), dimension(5), intent(in) :: sol_l
    real(kind=DOUBLE), dimension(5), intent(inout) :: sol_r
    real(kind=DOUBLE), intent(in) :: t

    real(kind=DOUBLE) :: r
    real(kind=DOUBLE), dimension(3) :: coord_re
    real(kind=DOUBLE), dimension(5) :: w

    if (re > 0) then
      print*, "re > 0, bad use of compute right sate !"
      error stop
    else if (re == 0) then !Default option wall
      sol_r(:) = sol_l
      sol_r(2:4) = sol_r(2:4) &
        - 2.0_DOUBLE*dot_product(sol_r(2:4), mesh%face(id_face)%norm)*mesh%face(id_face)%norm
    else !!Boundary
      if (trim(adjustl(bc_type(-re))) == &
        "outflowsupersonic") then
        sol_r = sol_l
      else if (trim(adjustl(bc_type(-re))) == &
        "freestream") then
        sol_r = primit_to_conserv(bc_val(:, -re))
      else if (trim(adjustl(bc_type(-re))) == &
        "sod") then
        sol_r = primit_to_conserv(sol_primit_sod_angle(&
          mesh%elem(le)%coord + 10*(mesh%face(id_face)%coord - mesh%elem(le)%coord), &
          x_sod, sod_angle, sol_left, sol_right))
      else if (trim(adjustl(bc_type(-re))) == &
        "inflow_pond") then
        sol_r = primit_to_conserv(bc_val(:, -re))
        if (dot_product(sol_r(2:4), mesh%face(id_face)%norm) > 0.0_DOUBLE) then
          sol_r = sol_l
        end if
      else if (trim(adjustl(bc_type(-re))) == &
        "inflow_pond_sod") then
        coord_re = mesh%face(id_face)%coord + (mesh%face(id_face)%coord - mesh%elem(le)%coord)
        sol_r = primit_to_conserv(bc_val(:, -re))
        if (dot_product(sol_r(2:4), mesh%face(id_face)%norm) > 0.0_DOUBLE) then
          sol_r = sol_l
        else
          if (coord_re(1) < x_sod - sin(sod_angle/360._DOUBLE*2._DOUBLE*pi)*coord_re(2)) then
            sol_r = primit_to_conserv(sol_left)
          else
            sol_r = primit_to_conserv(sol_right)
          end if
        end if
      else if (trim(adjustl(bc_type(-re))) == &
        "inout_double_mach") then

        coord_re = mesh%face(id_face)%coord + (mesh%face(id_face)%coord - mesh%elem(le)%coord)
        if (coord_re(2) > 1.732_DOUBLE*(coord_re(1) - 0.1667_DOUBLE - 10.0_DOUBLE*t)) then
          sol_r = primit_to_conserv((/8.0_DOUBLE, 7.145_DOUBLE, -4.125_DOUBLE, 0.0_DOUBLE, 116.5_DOUBLE/))
          if (sol_r(2) > 0.0_DOUBLE) then
            sol_r = sol_l
          end if
        else
          sol_r = primit_to_conserv((/1.4_DOUBLE, 0.0_DOUBLE, 0.0_DOUBLE, 0.0_DOUBLE, 1.0_DOUBLE/))
        end if

      else if (trim(adjustl(bc_type(-re))) == &
        "noh") then
        coord_re = mesh%face(id_face)%coord + (mesh%face(id_face)%coord - mesh%elem(le)%coord)
        r = sqrt(coord_re(1)**2 + coord_re(2)**2)
        w(1) = (1.0_DOUBLE + t/r)
        w(2) = -coord_re(1)/r
        w(3) = -coord_re(2)/r
        w(4) = 0.0_DOUBLE
        w(5) = w(1)*(gamma - 1)*1e-6
        sol_r = primit_to_conserv(w)
      else if ((trim(adjustl(bc_type(-re))) == "wall") .or. &
        (trim(adjustl(bc_type(-re))) == "")) then
        sol_r(:) = sol_l
        sol_r(2:4) = sol_r(2:4) &
          - 2.0_DOUBLE*dot_product(sol_r(2:4), mesh%face(id_face)%norm)*mesh%face(id_face)%norm
      else
        print *, "BC TYPE NOT RECOGNIZED !"
        error stop
      end if
    end if
  end subroutine compute_right_state

  subroutine write_dat_radial_cp(mesh, sol)
    use global_data_module
    implicit none

    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(in) :: sol

    integer(kind=ENTIER) :: i, funit
    real(kind=DOUBLE) :: r, cp, cpmax, r2d
    real(kind=DOUBLE) :: t1, t2, t3, t4, mnl1, mnl2, mnl3, mnl4
    real(kind=DOUBLE), dimension(5) :: w

    open (newunit=funit, file="radial_cp.dat", status="unknown")

    do i = 1, mesh%n_elems
      if (mesh%elem(i)%tag == cp_tag) then
        r = sqrt((mesh%elem(i)%coord(1) - radial_center(1))**2 &
          + (mesh%elem(i)%coord(2) - radial_center(2))**2 &
          + (mesh%elem(i)%coord(3) - radial_center(3))**2)

        w = conserv_to_primit(sol(:, i))
        cp = (w(5) - pinf)/(0.5_DOUBLE*rhoinf*vinf**2)
        cpmax = 2.0_DOUBLE/(gamma*vinf**2)*( &
          ((gamma + 1)**2*vinf**2/(4*gamma*vinf**2 - 2*(gamma - 1)))**(gamma/(gamma - 1)) &
          *(1 - gamma + 2*gamma*vinf**2)/(gamma + 1) - 1)

        r2d = sqrt(mesh%elem(i)%coord(2)**2 + mesh%elem(i)%coord(3)**2)
        t1 = atan2(r2d, mesh%elem(i)%coord(1))
        mnl1 = cpmax*sin(pi/2.0_DOUBLE - t1)**2
        t2 = atan2(r2d, -mesh%elem(i)%coord(1))
        mnl2 = cpmax*sin(pi/2.0_DOUBLE - t2)**2
        t3 = atan2(-r2d, -mesh%elem(i)%coord(1))
        mnl3 = cpmax*sin(pi/2.0_DOUBLE - t3)**2
        t4 = atan2(-r2d, mesh%elem(i)%coord(1))
        mnl4 = cpmax*sin(pi/2.0_DOUBLE - t4)**2

        write (funit, *) mesh%elem(i)%coord(1), r2d, &
          t1, mnl1, t2, mnl2, t3, mnl3, t4, mnl4, cp, cpmax
      end if
    end do

    close (funit)
  end subroutine write_dat_radial_cp

  subroutine write_dat_radial_sol(mesh, sol)
    use global_data_module
    implicit none

    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(in) :: sol

    integer(kind=ENTIER) :: i, funit
    real(kind=DOUBLE) :: r

    open (newunit=funit, file="radial_sol.dat", status="unknown")

    do i = 1, mesh%n_elems
      r = sqrt((mesh%elem(i)%coord(1) - radial_center(1))**2 &
        + (mesh%elem(i)%coord(2) - radial_center(2))**2 &
        + (mesh%elem(i)%coord(3) - radial_center(3))**2)
      write (funit, *) r, sol(1, i), sol(2, i), sol(3, i), sol(4, i)
    end do

    close (funit)
  end subroutine write_dat_radial_sol

  subroutine write_dat_dt(iter, dt)
    implicit none

    integer(kind=ENTIER), intent(in) :: iter
    real(kind=DOUBLE), intent(in) :: dt

    integer(kind=ENTIER) :: funit

    if (iter == 1) then
      open (newunit=funit, file="dt_iter.dat", status="unknown")
      write (funit, *) iter, dt
      close (funit)
    else
      open (newunit=funit, file="dt_iter.dat", status="old", position='append')
      write (funit, *) iter, dt
      close (funit)
    end if
  end subroutine write_dat_dt

  pure subroutine init_t_sol_aff(t, t_max, n_sol_aff, t_sol_aff)
    implicit none

    integer(kind=ENTIER), intent(in) :: n_sol_aff
    real(kind=DOUBLE), intent(in) :: t, t_max
    real(kind=DOUBLE), dimension(n_sol_aff), intent(inout) :: t_sol_aff

    integer(kind=ENTIER) :: i

    do i = 1, n_sol_aff
      t_sol_aff(i) = t + (i - 1)*(t_max - t)/(n_sol_aff - 1)
    end do
  end subroutine init_t_sol_aff

  pure function is_wall(re)
    use global_data_module, only: bc_type
    implicit none

    integer(kind=ENTIER), intent(in) :: re
    logical :: is_wall

    if (re > 0) then
      is_wall = .FALSE.
    else if (re == 0) then
      is_wall = .TRUE.
    else if (re < 0) then
      if ((trim(adjustl(bc_type(-re))) == "wall") &
        .or. (trim(adjustl(bc_type(-re))) == "")) then
        is_wall = .TRUE.
      else
        is_wall = .FALSE.
      end if
    end if
  end function is_wall

  subroutine compute_lambdas_and_solve_nodal_velocity(mesh, id_vert, sol_w_lr, lambda, u_bars, u_node, p_bound)
    use global_data_module, only: gamma, lambda_l_eq_lambda_r, bc_style, boundary_2d
    use subfv_linear_solver_module
    implicit none

    type(mesh_type), intent(in) :: mesh
    integer(kind=ENTIER), intent(in) :: id_vert
    real(kind=DOUBLE), dimension(5, 2, mesh%vert(id_vert)%n_sub_faces_neigh), &
      intent(in) :: sol_w_lr
    real(kind=DOUBLE), dimension(2, mesh%vert(id_vert)%n_sub_faces_neigh), &
      intent(inout) :: lambda
    real(kind=DOUBLE), dimension(mesh%vert(id_vert)%n_sub_faces_neigh), &
      intent(inout) :: u_bars
    real(kind=DOUBLE), dimension(3), intent(inout) :: u_node
    real(kind=DOUBLE), intent(inout) :: p_bound

    integer(kind=ENTIER) :: iter, n_iter
    integer(kind=ENTIER) :: j, id_sub_face, id_face, le, re
    real(kind=DOUBLE) :: rhol, ul, pl, al, lambda_l
    real(kind=DOUBLE) :: rhor, ur, pr, ar, lambda_r
    real(kind=DOUBLE) :: u_et
    real(kind=DOUBLE), dimension(3) :: Rp, Bp, u_node_old
    real(kind=DOUBLE), dimension(3, 3) :: mat, mat_inv

    n_iter = 20
    u_node = 0.0_DOUBLE
    u_node_old = 1.0_DOUBLE
    iter = 0
    do while ( maxval(abs(u_node_old - u_node)) > 1e-12_DOUBLE &
        .and. iter < n_iter )

      iter = iter + 1
      u_node_old = u_node

      mat(:, :) = 0.0_DOUBLE
      mat_inv(:, :) = 0.0_DOUBLE
      Rp(:) = 0.0_DOUBLE

      !Build Bp
      Bp(:) = 0.0_DOUBLE
      do j = 1, mesh%vert(id_vert)%n_sub_faces_neigh
        id_sub_face = mesh%vert(id_vert)%sub_face_neigh(j)
        id_face = mesh%sub_face(id_sub_face)%mesh_face

        le = mesh%sub_face(id_sub_face)%left_elem_neigh
        re = mesh%sub_face(id_sub_face)%right_elem_neigh

        if(is_wall(re)) then
          if(boundary_2d) then
            if( abs(dot_product(mesh%face(id_face)%norm, (/0.0_DOUBLE, 0.0_DOUBLE, 1.0_DOUBLE/))) &
              < 1e-8_DOUBLE) then
              Bp = Bp + mesh%sub_face(id_sub_face)%area*mesh%face(id_face)%norm
            end if
          else
            Bp = Bp + mesh%sub_face(id_sub_face)%area*mesh%face(id_face)%norm
          end if
        end if
      end do

      !Compute initial lambdas
      do j = 1, mesh%vert(id_vert)%n_sub_faces_neigh
        id_sub_face = mesh%vert(id_vert)%sub_face_neigh(j)
        id_face = mesh%sub_face(id_sub_face)%mesh_face

        le = mesh%sub_face(id_sub_face)%left_elem_neigh
        re = mesh%sub_face(id_sub_face)%right_elem_neigh

        rhol = sol_w_lr(1, 1, j)
        ul = sol_w_lr(2, 1, j)
        pl = sol_w_lr(5, 1, j)
        al = sqrt(gamma*pl/rhol)

        rhor = sol_w_lr(1, 2, j)
        ur = sol_w_lr(2, 2, j)
        pr = sol_w_lr(5, 2, j)
        ar = sqrt(gamma*pr/rhor)

        lambda_l = lambda(1, j)
        lambda_r = lambda(2, j)

        if (iter == 1) then
          lambda_l = max(1e-6_DOUBLE, lambda_l, al*rhol, 1.01_DOUBLE*sqrt(rhol*max(0.0_DOUBLE, pr - pl)), &
            -1.01_DOUBLE*rhol*(ur - ul))
          lambda_r = max(1e-6_DOUBLE, lambda_r, ar*rhor, 1.01_DOUBLE*sqrt(rhor*max(0.0_DOUBLE, pl - pr)), &
            -1.01_DOUBLE*rhor*(ur - ul))

          u_bars(j) = (lambda_l*ul + lambda_r*ur - (pr - pl))/(lambda_r + lambda_l)

          lambda_l = max(lambda_l, al*rhol*(1.0_DOUBLE + 1.5_DOUBLE*max(0.0_DOUBLE, -(u_bars(j) - ul)/al)))
          lambda_r = max(lambda_r, ar*rhor*(1.0_DOUBLE + 1.5_DOUBLE*max(0.0_DOUBLE, (u_bars(j) - ur)/ar)))
        else
          if (is_wall(re) .and. bc_style == 1) then
            lambda_l = max(lambda_l, al*rhol*(1.0_DOUBLE + 1.5_DOUBLE*max(0.0_DOUBLE, -(u_bars(j) - ul)/al)))
            lambda_r = max(lambda_r, ar*rhor*(1.0_DOUBLE + 1.5_DOUBLE*max(0.0_DOUBLE, (u_bars(j) - ur)/ar)))
          else
            u_et = dot_product(u_node, mesh%face(id_face)%norm)
            lambda_l = max(lambda_l, al*rhol*(1.0_DOUBLE + 1.5_DOUBLE*max(0.0_DOUBLE, -(u_et - ul)/al)))
            lambda_r = max(lambda_r, ar*rhor*(1.0_DOUBLE + 1.5_DOUBLE*max(0.0_DOUBLE, (u_et - ur)/ar)))
          end if
        end if

        if (lambda_l_eq_lambda_r) then
          lambda_l = max(lambda_l, lambda_r)
          lambda_r = lambda_l
        end if

        u_bars(j) = (lambda_l*ul + lambda_r*ur - (pr - pl))/(lambda_r + lambda_l)

        lambda(1, j) = max(lambda(1, j), lambda_l)
        lambda(2, j) = max(lambda(2, j), lambda_r)

        lambda_l = lambda(1, j)
        lambda_r = lambda(2, j)

        if (is_wall(re)) then
          if (bc_style == 0 .or. bc_style == 2) then
            mat = mat + (lambda_r + lambda_l)*mesh%sub_face(id_sub_face)%area &
              *tensor_product(mesh%face(id_face)%norm,mesh%face(id_face)%norm)

            if(bc_style == 0) then
              Rp = Rp + mesh%sub_face(id_sub_face)%area* &
                (lambda_r + lambda_l)*u_bars(j)*mesh%face(id_face)%norm
            else if(bc_style == 2) then
              if(boundary_2d) then
                if( abs(dot_product(mesh%face(id_face)%norm, (/0.0_DOUBLE, 0.0_DOUBLE, 1.0_DOUBLE/))) &
                  < 1e-8_DOUBLE) then
                  Rp = Rp + mesh%sub_face(id_sub_face)%area*pl*mesh%face(id_face)%norm
                end if
              else
                Rp = Rp + mesh%sub_face(id_sub_face)%area*pl*mesh%face(id_face)%norm
              end if
            end if
          end if
        else
          mat = mat + (lambda_r + lambda_l)*mesh%sub_face(id_sub_face)%area &
            *tensor_product(mesh%face(id_face)%norm,mesh%face(id_face)%norm)
          Rp = Rp + mesh%sub_face(id_sub_face)%area* &
            (lambda_r + lambda_l)*u_bars(j)*mesh%face(id_face)%norm
        end if
      end do

      if(bc_style == 2)then
        if(maxval(abs(Bp)) > 1e-8_DOUBLE) then
          call inverse_3_by_3(mat, mat_inv)
          p_bound = dot_product(matmul(mat_inv, Rp),Bp)/dot_product(matmul(mat_inv, Bp),Bp)
          u_node = matmul(mat_inv, Rp - p_bound*Bp)

          if( abs(dot_product(u_node, Bp)) > 1e-12_DOUBLE) then
            print*, "Boundary velocity not orthogonal to average boundary normal !!!"
            print*, Bp
            print*, u_node
            error stop
          end if
        else
          p_bound = -1.0_DOUBLE
          call qr_solve(3, 3, mat, u_node, Rp)
        end if
      else
        call qr_solve(3, 3, mat, u_node, Rp)
      end if
    end do

    if( iter >= n_iter ) then
      print*, "Multi-point nodal solver iter maxiter reached", iter, n_iter
      error stop
    end if
  end subroutine compute_lambdas_and_solve_nodal_velocity

  subroutine compute_lambdas_and_solve_nodal_velocity_iso(mesh, id_vert, sol_w_lr, lambda, u_bars, u_node, p_bound)
    use global_data_module, only: gamma, lambda_l_eq_lambda_r, bc_style, boundary_2d
    use subfv_linear_solver_module
    implicit none

    type(mesh_type), intent(in) :: mesh
    integer(kind=ENTIER), intent(in) :: id_vert
    real(kind=DOUBLE), dimension(5, 2, mesh%vert(id_vert)%n_sub_faces_neigh), &
      intent(in) :: sol_w_lr
    real(kind=DOUBLE), dimension(2, mesh%vert(id_vert)%n_sub_faces_neigh), &
      intent(inout) :: lambda
    real(kind=DOUBLE), dimension(mesh%vert(id_vert)%n_sub_faces_neigh), &
      intent(inout) :: u_bars
    real(kind=DOUBLE), dimension(3), intent(inout) :: u_node
    real(kind=DOUBLE), intent(inout) :: p_bound

    integer(kind=ENTIER) :: iter, n_iter
    integer(kind=ENTIER) :: j, id_sub_face, id_face, le, re
    real(kind=DOUBLE) :: rhol, ul, pl, al, lambda_l
    real(kind=DOUBLE) :: rhor, ur, pr, ar, lambda_r
    real(kind=DOUBLE) :: u_et
    real(kind=DOUBLE), dimension(3) :: Rp, Bp, u_node_old
    real(kind=DOUBLE), dimension(3, 3) :: mat, mat_inv

    n_iter = 20
    u_node = 0.0_DOUBLE
    u_node_old = 1.0_DOUBLE
    iter = 0
    do while ( maxval(abs(u_node_old - u_node)) > 1e-12_DOUBLE &
        .and. iter < n_iter )

      iter = iter + 1
      u_node_old = u_node

      mat(:, :) = 0.0_DOUBLE
      mat_inv(:, :) = 0.0_DOUBLE
      Rp(:) = 0.0_DOUBLE

      !Build Bp
      Bp(:) = 0.0_DOUBLE
      do j = 1, mesh%vert(id_vert)%n_sub_faces_neigh
        id_sub_face = mesh%vert(id_vert)%sub_face_neigh(j)
        id_face = mesh%sub_face(id_sub_face)%mesh_face

        le = mesh%sub_face(id_sub_face)%left_elem_neigh
        re = mesh%sub_face(id_sub_face)%right_elem_neigh

        if(is_wall(re)) then
          if(boundary_2d) then
            if( abs(dot_product(mesh%face(id_face)%norm, (/0.0_DOUBLE, 0.0_DOUBLE, 1.0_DOUBLE/))) &
              < 1e-8_DOUBLE) then
              Bp = Bp + mesh%sub_face(id_sub_face)%area*mesh%face(id_face)%norm
            end if
          else
            Bp = Bp + mesh%sub_face(id_sub_face)%area*mesh%face(id_face)%norm
          end if
        end if
      end do

      !Compute initial lambdas
      do j = 1, mesh%vert(id_vert)%n_sub_faces_neigh
        id_sub_face = mesh%vert(id_vert)%sub_face_neigh(j)
        id_face = mesh%sub_face(id_sub_face)%mesh_face

        le = mesh%sub_face(id_sub_face)%left_elem_neigh
        re = mesh%sub_face(id_sub_face)%right_elem_neigh

        rhol = sol_w_lr(1, 1, j)
        ul = sol_w_lr(2, 1, j)
        pl = sol_w_lr(5, 1, j)
        al = sqrt(gamma*pl/rhol)

        rhor = sol_w_lr(1, 2, j)
        ur = sol_w_lr(2, 2, j)
        pr = sol_w_lr(5, 2, j)
        ar = sqrt(gamma*pr/rhor)

        lambda_l = lambda(1, j)
        lambda_r = lambda(2, j)

        if (iter == 1) then
          lambda_l = max(1e-6_DOUBLE, lambda_l, al*rhol, 1.01_DOUBLE*sqrt(rhol*max(0.0_DOUBLE, pr - pl)), &
            -1.01_DOUBLE*rhol*(ur - ul))
          lambda_r = max(1e-6_DOUBLE, lambda_r, ar*rhor, 1.01_DOUBLE*sqrt(rhor*max(0.0_DOUBLE, pl - pr)), &
            -1.01_DOUBLE*rhor*(ur - ul))

          u_bars(j) = (lambda_l*ul + lambda_r*ur - (pr - pl))/(lambda_r + lambda_l)

          lambda_l = max(lambda_l, al*rhol*(1.0_DOUBLE + 1.5_DOUBLE*max(0.0_DOUBLE, -(u_bars(j) - ul)/al)))
          lambda_r = max(lambda_r, ar*rhor*(1.0_DOUBLE + 1.5_DOUBLE*max(0.0_DOUBLE, (u_bars(j) - ur)/ar)))
        else
          if (is_wall(re) .and. bc_style == 1) then
            lambda_l = max(lambda_l, al*rhol*(1.0_DOUBLE + 1.5_DOUBLE*max(0.0_DOUBLE, -(u_bars(j) - ul)/al)))
            lambda_r = max(lambda_r, ar*rhor*(1.0_DOUBLE + 1.5_DOUBLE*max(0.0_DOUBLE, (u_bars(j) - ur)/ar)))
          else
            u_et = dot_product(u_node, mesh%face(id_face)%norm)
            lambda_l = max(lambda_l, al*rhol*(1.0_DOUBLE + 1.5_DOUBLE*max(0.0_DOUBLE, -(u_et - ul)/al)))
            lambda_r = max(lambda_r, ar*rhor*(1.0_DOUBLE + 1.5_DOUBLE*max(0.0_DOUBLE, (u_et - ur)/ar)))
          end if
        end if

        if (lambda_l_eq_lambda_r) then
          lambda_l = max(lambda_l, lambda_r)
          lambda_r = lambda_l
        end if

        u_bars(j) = (lambda_l*ul + lambda_r*ur - (pr - pl))/(lambda_r + lambda_l)

        lambda(1, j) = max(lambda(1, j), lambda_l)
        lambda(2, j) = max(lambda(2, j), lambda_r)

        lambda_l = lambda(1, j)
        lambda_r = lambda(2, j)

        if (is_wall(re)) then
          if (bc_style == 0 .or. bc_style == 2) then
            mat = mat + (lambda_r + lambda_l)&
              *tensor_product(mesh%face(id_face)%norm,mesh%face(id_face)%norm)

            if(bc_style == 0) then
              Rp = Rp &
                + (lambda_r + lambda_l)*u_bars(j)*mesh%face(id_face)%norm
            else if(bc_style == 2) then
              if(boundary_2d) then
                if( abs(mesh%face(id_face)%norm(3)) < 1e-8_DOUBLE) then
                  Rp = Rp + pl*mesh%face(id_face)%norm
                end if
              else
                Rp = Rp + pl*mesh%face(id_face)%norm
              end if
            end if
          end if
        else
          mat = mat &
            + (lambda_r + lambda_l)*tensor_product(mesh%face(id_face)%norm,mesh%face(id_face)%norm)
          Rp = Rp + (lambda_r + lambda_l)*u_bars(j)*mesh%face(id_face)%norm
        end if
      end do

      if(bc_style == 2)then
        if(maxval(abs(Bp)) > 1e-8_DOUBLE) then
          call inverse_3_by_3(mat, mat_inv)
          p_bound = dot_product(matmul(mat_inv, Rp),Bp)/dot_product(matmul(mat_inv, Bp),Bp)
          u_node = matmul(mat_inv, Rp - p_bound*Bp)

          if( abs(dot_product(u_node, Bp)) > 1e-12_DOUBLE) then
            print*, "Boundary velocity not orthogonal to average boundary normal !!!"
            print*, Bp
            print*, u_node
            error stop
          end if
        else
          p_bound = -1.0_DOUBLE
          call qr_solve(3, 3, mat, u_node, Rp)
        end if
      else
        call qr_solve(3, 3, mat, u_node, Rp)
      end if
    end do

    if( iter >= n_iter ) then
      print*, "Multi-point nodal solver iter maxiter reached", iter, n_iter
      error stop
    end if
  end subroutine compute_lambdas_and_solve_nodal_velocity_iso

  pure function primit_to_conserv(w) result(u)
    use global_data_module, only: gamma
    implicit none

    real(kind=DOUBLE), dimension(5), intent(in) :: w
    real(kind=DOUBLE), dimension(5) :: u

    u(1) = w(1)
    u(2) = w(2)*w(1)
    u(3) = w(3)*w(1)
    u(4) = w(4)*w(1)
    u(5) = w(5)/(gamma - 1) + 0.5_DOUBLE*w(1)*(w(2)**2 + w(3)**2 + w(4)**2)
  end function primit_to_conserv

  function conserv_to_primit(u) result(w)
    use global_data_module, only: gamma
    implicit none

    real(kind=DOUBLE), dimension(5), intent(in) :: u
    real(kind=DOUBLE), dimension(5) :: w

    w(1) = u(1)
    w(2) = u(2)/u(1)
    w(3) = u(3)/u(1)
    w(4) = u(4)/u(1)
    w(5) = (gamma - 1)*(u(5) - 0.5_DOUBLE*w(1)*(w(2)**2 + w(3)**2 + w(4)**2))

    if (u(1) < 0.0_DOUBLE .or. u(5) < 0.0_DOUBLE) then
      print *, "error negatif conserv primit input"
      print *, u
      error stop
    else if (w(1) < 0.0_DOUBLE .or. w(5) < 0.0_DOUBLE) then
      print *, "error in conserv to primit output"
      print *, u
      print *, w
      print *, (u(5) - 0.5_DOUBLE*w(1)*(w(2)**2 + w(3)**2 + w(4)**2))
      print *, (u(5) - 0.5_DOUBLE*w(1)*(w(2)**2 + w(3)**2 + w(4)**2))
      error stop
    end if
  end function conserv_to_primit

  pure subroutine sol_isentropic_vortex(coord, w, t)
    use global_data_module
    implicit none

    real(kind=DOUBLE), intent(in) :: t
    real(kind=DOUBLE), dimension(3), intent(in) :: coord
    real(kind=DOUBLE), dimension(5), intent(inout) :: w

    real(kind=DOUBLE), dimension(3) :: center_coord, vel
    real(kind=DOUBLE) :: beta, r

    beta = 5.0_DOUBLE
    vel(:) = (/0.0_DOUBLE, 0.0_DOUBLE, 0.0_DOUBLE/)
    center_coord(:) = (/0.0_DOUBLE, 0.0_DOUBLE, 0.0_DOUBLE/) + vel*t
    r = (coord(1) - center_coord(1))**2 + (coord(2) - center_coord(2))**2
    w(1) = 1.0_DOUBLE*(1.0_DOUBLE - ((gamma - 1)*beta**2)/(8.0_DOUBLE*gamma*pi**2)* &
      exp(1.0_DOUBLE - r))**(1.0_DOUBLE/(gamma - 1.0_DOUBLE))
    w(2) = vel(1) &
      - (coord(2) - center_coord(2))*beta/(2.0_DOUBLE*pi)*exp(0.5_DOUBLE*(1.0_DOUBLE - r))
    w(3) = vel(2) &
      + (coord(1) - center_coord(1))*beta/(2.0_DOUBLE*pi)*exp(0.5_DOUBLE*(1.0_DOUBLE - r))
    w(4) = 0.0_DOUBLE
    w(5) = w(1)**gamma
  end subroutine sol_isentropic_vortex

  subroutine reconstruct_with_grad(sol, grad, v)
    implicit none

    real(kind=DOUBLE), dimension(5), intent(inout) :: sol
    real(kind=DOUBLE), dimension(3), intent(in) :: v
    real(kind=DOUBLE), dimension(3, 5), intent(in) :: grad

    integer(kind=ENTIER) :: var

    do var = 1, 5
      sol(var) = sol(var) + dot_product(grad(:, var), v)
    end do

    if (sol(1) < 0.0_DOUBLE .or. sol(5) < 0.0_DOUBLE) then
      print *, "Bad reconstruction !"
      print *, "recons sol", sol
      print *, "grad 1", grad(:3, 1)
      print *, "grad 5", grad(:3, 5)
      error stop
    end if
  end subroutine reconstruct_with_grad

  function kelvin_helmholtz(y) result(v)
    use global_data_module, only: pi
    implicit none

    real(kind=DOUBLE), intent(in) :: y
    real(kind=DOUBLE) :: v

    real(kind=DOUBLE) :: w = 1.0_DOUBLE/16.0_DOUBLE

    if (y > -1.0_DOUBLE/4.0_DOUBLE - w/2.0_DOUBLE .and. &
      y < -1.0_DOUBLE/4.0_DOUBLE + w/2.0_DOUBLE) then
      v = -sin(pi*(y + 1.0_DOUBLE/4.0_DOUBLE)/w)
    else if (y > -1.0_DOUBLE/4.0_DOUBLE + w/2.0_DOUBLE .and. &
      y < 1.0_DOUBLE/4.0_DOUBLE - w/2.0_DOUBLE) then
      v = -1.0_DOUBLE
    else if (y > 1.0_DOUBLE/4.0_DOUBLE - w/2.0_DOUBLE .and. &
      y < 1.0_DOUBLE/4.0_DOUBLE + w/2.0_DOUBLE) then
      v = sin(pi*(y - 1.0_DOUBLE/4.0_DOUBLE)/w)
    else
      v = 1.0_DOUBLE
    end if
  end function kelvin_helmholtz

  pure subroutine base_change(v1, v2, n1, forward)
    implicit none

    real(kind=DOUBLE), dimension(3), intent(in) :: n1
    real(kind=DOUBLE), dimension(3), intent(inout) :: v1
    real(kind=DOUBLE), dimension(3), intent(inout), optional :: v2
    logical, intent(in) :: forward

    real(kind=DOUBLE) :: n2n, n3n
    real(kind=DOUBLE), dimension(3) :: n2, n3
    real(kind=DOUBLE), dimension(3) :: v1tmp, v2tmp

    n2 = (/n1(3), 0.0_DOUBLE, -n1(1)/)
    n2n = norm2(n2)
    if (n2n > 1e-8) then
      n2 = n2/n2n
    else
      n2 = (/0.0_DOUBLE, -n1(3), n1(2)/)
      n2n = norm2(n2)
      n2 = n2/n2n
    end if

    n3(1) = n1(2)*n2(3) - n1(3)*n2(2)
    n3(2) = -(n1(1)*n2(3) - n1(3)*n2(1))
    n3(3) = n1(1)*n2(2) - n1(2)*n2(1)

    n3n = norm2(n3)
    n3 = n3/n3n

    if (forward) then
      v1tmp(1) = v1(1)*n1(1) + v1(2)*n1(2) + v1(3)*n1(3)
      v1tmp(2) = v1(1)*n2(1) + v1(2)*n2(2) + v1(3)*n2(3)
      v1tmp(3) = v1(1)*n3(1) + v1(2)*n3(2) + v1(3)*n3(3)

      if (present(v2)) then
        v2tmp(1) = v2(1)*n1(1) + v2(2)*n1(2) + v2(3)*n1(3)
        v2tmp(2) = v2(1)*n2(1) + v2(2)*n2(2) + v2(3)*n2(3)
        v2tmp(3) = v2(1)*n3(1) + v2(2)*n3(2) + v2(3)*n3(3)
      end if
    else
      v1tmp(1) = v1(1)*n1(1) + v1(2)*n2(1) + v1(3)*n3(1)
      v1tmp(2) = v1(1)*n1(2) + v1(2)*n2(2) + v1(3)*n3(2)
      v1tmp(3) = v1(1)*n1(3) + v1(2)*n2(3) + v1(3)*n3(3)

      if (present(v2)) then
        v2tmp(1) = v2(1)*n1(1) + v2(2)*n2(1) + v2(3)*n3(1)
        v2tmp(2) = v2(1)*n1(2) + v2(2)*n2(2) + v2(3)*n3(2)
        v2tmp(3) = v2(1)*n1(3) + v2(2)*n2(3) + v2(3)*n3(3)
      end if
    end if

    v1 = v1tmp
    if (present(v2)) v2 = v2tmp
  end subroutine base_change

  subroutine compute_nodal_grad(mesh, id_vert, sol, nodal_grad)
    use subfv_linear_solver_module, only: tensor_product, &
      tensor_product_3,&
      pseudo_inverse_inplace_lapack, print_mat
    implicit none

    type(mesh_type), intent(in) :: mesh
    integer(kind=ENTIER), intent(in) :: id_vert
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(in) :: sol
    real(kind=DOUBLE), dimension(3, 5), intent(inout) :: nodal_grad

    integer(kind=ENTIER) :: i, id_sub_face, id_face, le, re
    real(kind=DOUBLE), dimension(5) :: sol_w_l, sol_w_r
    real(kind=DOUBLE), dimension(3, 3) :: mat

    nodal_grad = 0.0_DOUBLE
    mat = 0.0_DOUBLE
    do i=1, mesh%vert(id_vert)%n_sub_faces_neigh
      id_sub_face = mesh%vert(id_vert)%sub_face_neigh(i)
      id_face = mesh%sub_face(id_sub_face)%mesh_face
      le = mesh%face(id_face)%left_neigh
      re = mesh%face(id_face)%right_neigh
      if( re > 0 ) then
        sol_w_l = conserv_to_primit(sol(:, le))
        sol_w_r = conserv_to_primit(sol(:, re))
        nodal_grad = nodal_grad &
          + mesh%sub_face(id_sub_face)%area&
          * tensor_product(mesh%face(id_face)%norm, &
          sol_w_r - sol_w_l)
        mat = mat &
          + mesh%sub_face(id_sub_face)%area&
          * tensor_product_3(mesh%face(id_face)%norm, &
          mesh%elem(re)%coord - mesh%elem(le)%coord)
      end if
    end do
    call pseudo_inverse_inplace_lapack(3, mat)
    nodal_grad = matmul(mat, nodal_grad)
    if( mesh%vert(id_vert)%is_bound ) nodal_grad = 0.0_DOUBLE
  end subroutine compute_nodal_grad

  subroutine compute_nodal_grad_cell_grad(mesh, id_vert, sol, grad, weight_sum)
    use subfv_linear_solver_module, only: tensor_product, &
      tensor_product_3,&
      pseudo_inverse_inplace_lapack, print_mat
    use global_data_module, only: method
    implicit none

    type(mesh_type), intent(in) :: mesh
    integer(kind=ENTIER), intent(in) :: id_vert
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(in) :: sol
    real(kind=DOUBLE), dimension(3, 5, mesh%n_elems), intent(inout) :: grad
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(inout) :: weight_sum

    integer(kind=ENTIER) :: i, id_sub_face, id_face, le, re
    integer(kind=ENTIER) :: id_elem, id_sub_elem, var, k
    real(kind=DOUBLE) :: weight
    real(kind=DOUBLE), dimension(5) :: sol_w_l, sol_w_r
    real(kind=DOUBLE), dimension(3, 3) :: mat
    real(kind=DOUBLE), dimension(3, 5) :: nodal_grad

    real(kind=DOUBLE), parameter :: eps=1e-3_DOUBLE

    !Compute nodal grad
    nodal_grad = 0.0_DOUBLE
    mat = 0.0_DOUBLE
    do i=1, mesh%vert(id_vert)%n_sub_faces_neigh
      id_sub_face = mesh%vert(id_vert)%sub_face_neigh(i)
      id_face = mesh%sub_face(id_sub_face)%mesh_face
      le = mesh%face(id_face)%left_neigh
      re = mesh%face(id_face)%right_neigh
      if( re > 0 ) then
        sol_w_l = conserv_to_primit(sol(:, le))
        sol_w_r = conserv_to_primit(sol(:, re))
        nodal_grad = nodal_grad &
          + mesh%sub_face(id_sub_face)%area&
          * tensor_product(mesh%face(id_face)%norm, &
          sol_w_r - sol_w_l)
        mat = mat &
          + mesh%sub_face(id_sub_face)%area&
          * tensor_product_3(mesh%face(id_face)%norm, &
          mesh%elem(re)%coord - mesh%elem(le)%coord)
      end if
    end do
    call pseudo_inverse_inplace_lapack(3, mat)
    nodal_grad = matmul(mat, nodal_grad)
    if( mesh%vert(id_vert)%is_bound ) nodal_grad = 0.0_DOUBLE

    !Accumulate 
    do i=1, mesh%vert(id_vert)%n_sub_elems_neigh
      id_sub_elem = mesh%vert(id_vert)%sub_elem_neigh(i)
      id_elem = mesh%sub_elem(id_sub_elem)%mesh_elem

      if( method == 2 ) then
        weight = mesh%sub_elem(id_sub_elem)%volume
        do var=1, 5
          do k=1, 3
            !$OMP ATOMIC UPDATE
            grad(k, var, id_elem) = grad(k, var, id_elem) &
              + weight * nodal_grad(k, var)
          end do
          !$OMP ATOMIC UPDATE
          weight_sum(var, id_elem) = weight_sum(var, id_elem) &
            + weight
        end do

      else if( method == 3 ) then
        !$OMP CRITICAL
        if( weight_sum(1, id_elem) < 1e-2_DOUBLE ) then
          do var = 1, 5
            do k=1, 3
              grad(k, var, id_elem) = nodal_grad(k, var)
            end do
            weight_sum(var, id_elem) = 1.0_DOUBLE
          end do
        else
          do var = 1, 5
            if( OI(nodal_grad(:, var)) < OI(grad(:, var, id_elem)) ) then
              do k=1, 3
                grad(k, var, id_elem) = nodal_grad(k, var)
              end do
              weight_sum(var, id_elem) = 1.0_DOUBLE
            end if
          end do
        end if
        !$OMP END CRITICAL
      else if( method == 4 ) then
        do var=1, 5
          weight = mesh%sub_elem(id_sub_elem)%volume/(eps + OI(nodal_grad(:, var)))
          do k=1, 3
            !$OMP ATOMIC UPDATE
            grad(k, var, id_elem) = grad(k, var, id_elem) &
              + weight * nodal_grad(k, var)
          end do
          !$OMP ATOMIC UPDATE
          weight_sum(var, id_elem) = weight_sum(var, id_elem) &
            + weight
        end do

      else if( method == 5 ) then
        do var=1, 5
          weight = mesh%sub_elem(id_sub_elem)%volume/(eps + OI(nodal_grad(:, var)))**2
          do k=1, 3
            !$OMP ATOMIC UPDATE
            grad(k, var, id_elem) = grad(k, var, id_elem) &
              + weight * nodal_grad(k, var)
          end do
          !$OMP ATOMIC UPDATE
          weight_sum(var, id_elem) = weight_sum(var, id_elem) &
            + weight
        end do
      else
        print*, "Error unknown second order method nodal_grad_cell_grad"
        error stop
      end if
    end do
  end subroutine compute_nodal_grad_cell_grad

  function max_omega_for_positivity(U, deltaU) result(omega)
    implicit none

    real(kind=DOUBLE), dimension(5), intent(in) :: U, deltaU
    real(kind=DOUBLE) :: omega

    real(kind=DOUBLE), parameter :: tol = 1e-16_DOUBLE
    real(kind=DOUBLE), parameter :: safe = 0.9_DOUBLE
    real(kind=DOUBLE) :: a, b, c, discrim, x1, x2

    omega = 1.0_DOUBLE
    if( maxval(abs(deltaU)) < tol ) return
    if( U(1)+omega*deltaU(1) < 0.0_DOUBLE ) then
      omega = min(1.0_DOUBLE, -safe*U(1)/deltaU(1))
    end if

    a = deltaU(1)*deltaU(5) - 0.5_DOUBLE*dot_product(deltaU(2:4), deltaU(2:4))
    b = U(5)*deltaU(1) + U(1)*deltaU(5) - dot_product(U(2:4), deltaU(2:4))
    c = U(1)*U(5) - 0.5_DOUBLE*dot_product(U(2:4), U(2:4))

    if( a*omega**2+b*omega+c > 0.0_DOUBLE ) return
    if( c < 0.0_DOUBLE ) then
      print*," ERROR INPUT OMEGA POS "
      error stop
    end if
    if( abs(a) < tol .and. abs(b) > tol ) then
      omega = min(omega, -safe*c/b)
      return
    end if

    discrim = b**2 - 4.0_DOUBLE*a*c
    if( discrim < 0.0_DOUBLE ) then
      return
    else
      x1 = (-b+sqrt(discrim))/(2.0_DOUBLE*a)
      if(x1 < 0.0_DOUBLE) x1 = 10.0_DOUBLE
      x2 = (-b-sqrt(discrim))/(2.0_DOUBLE*a)
      if(x2 < 0.0_DOUBLE) x2 = 10.0_DOUBLE
      omega = min(omega, safe*x1, safe*x2)
    end if
  end function max_omega_for_positivity

  subroutine compute_grad_least_squares(mesh, sol, grad)
    use global_data_module, only: scheme
    use subfv_linear_solver_module
    implicit none

    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(in) :: sol
    real(kind=DOUBLE), dimension(:, :, :), intent(inout) :: grad

    integer(kind=ENTIER) :: i

    !$OMP PARALLEL DO 
    do i = 1, mesh%n_elems
      call compute_cell_grad_least_squares(mesh, sol, i, grad)
      call physical_slope_limiter(mesh, sol, i, grad)
      call relaxed_maximum_principle(mesh, sol, i, grad)
    end do
  end subroutine compute_grad_least_squares

  subroutine compute_cell_grad_least_squares(mesh, sol, i, grad)
    use subfv_linear_solver_module, only: tensor_product, &
      pseudo_inverse_inplace_lapack
    implicit none

    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(in) :: sol
    integer(kind=ENTIER), intent(in) :: i
    real(kind=DOUBLE), dimension(:, :, :), intent(inout) :: grad

    integer(kind=ENTIER) :: j, k, var, id_vert
    real(kind=DOUBLE), dimension(3, 3) :: mat
    real(kind=DOUBLE), dimension(3, 5) :: rhs
    real(kind=DOUBLE), dimension(5) :: sol_wi, sol_wj
    real(kind=DOUBLE), dimension(3) :: v

    rhs = 0.0_DOUBLE
    mat = 0.0_DOUBLE
    sol_wi = conserv_to_primit(sol(:, i))
    do k=1, mesh%elem(i)%n_neigh_by_vert
      j = mesh%elem(i)%neigh_by_vert(k)
      sol_wj = conserv_to_primit(sol(:, j))
      v = mesh%elem(j)%coord - mesh%elem(i)%coord
      do var=1, 5
        rhs(:, var) = rhs(:, var) + (sol_wj(var) - sol_wi(var))*v
      end do
      mat = mat + tensor_product(v, v)
    end do
    call pseudo_inverse_inplace_lapack(3, mat)
    do var=1, 5
      grad(:, var, i) = matmul(mat, rhs(:, var))
    end do

    do k=1, mesh%elem(i)%n_vert
      id_vert = mesh%elem(i)%vert(k)
      if( mesh%vert(id_vert)%is_bound ) grad(:, :, i) = 0.0_DOUBLE
    end do
  end subroutine compute_cell_grad_least_squares

  function euler_flux(U)
    use subfv_linear_solver_module, only: tensor_product
    implicit none

    real(kind=DOUBLE), dimension(5), intent(in) :: U
    real(kind=DOUBLE), dimension(5,3) :: euler_flux

    integer(kind=ENTIER) :: i
    real(kind=DOUBLE), dimension(5) :: W
    real(kind=DOUBLE), dimension(3,3) :: identity

    identity = 0.0_DOUBLE
    do i=1, 3
      identity(i, i) = 1.0_DOUBLE
    end do

    W = conserv_to_primit(U)

    euler_flux(1, :) = U(2:4)
    euler_flux(2:4, :) = tensor_product(U(2:4), W(2:4)) + W(5)*identity
    euler_flux(5, :) = (U(5) + W(5))*W(2:4)
  end function euler_flux

  subroutine print_mat(mat)
    implicit none

    integer(kind=ENTIER) :: i, j
    real(kind=DOUBLE), intent(in) :: mat(:, :)

    print *, ""
    do i = 1, size(mat, dim=1)
      do j = 1, size(mat, dim=2)
        if (mat(i, j) > 1e-12_DOUBLE) then
          write (*, '(a,f16.8,a)', advance="no") ""//achar(27)//"[36m", mat(i, j), ""//achar(27)//"[0m"
        else if (mat(i, j) < -1e-12_DOUBLE) then
          write (*, '(a,f16.8,a)', advance="no") ""//achar(27)//"[34m", mat(i, j), ""//achar(27)//"[0m"
        else
          write (*, '(f16.8)', advance="no") mat(i, j)
        end if
      end do
      write (*, *) " "
    end do
    print *, ""
  end subroutine print_mat

  subroutine write_dat_residual(mesh, sol, new_sol, flux_sum, iter, t, dt, tstart, funit, &
      num_procs, mpi_send_recv)
    use global_data_module, only: restart_iter, restart_time, restart_cpu_time
    use mpi
    use subfv_mpi_module
    use subfv_sparse_csr_linear_module, only: mpi_norm2
    implicit none

    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(in) :: sol, new_sol, flux_sum
    integer(kind=ENTIER), intent(in) :: iter, funit, num_procs
    real(kind=DOUBLE), intent(in) :: t, tstart, dt
    type(mpi_send_recv_type) :: mpi_send_recv

    integer(kind=ENTIER) :: i, j, mpi_ierr
    real(kind=DOUBLE) :: r1, r2, r3, tcur, r4

    if( num_procs > 1 ) then 
      r1 = mpi_norm2(mesh%n_elems, 5, new_sol - sol, mpi_send_recv)
      r4 = mpi_norm2(mesh%n_elems, 5, flux_sum, mpi_send_recv)
    else
      r1 = norm2(new_sol - sol)
      r4 = norm2(flux_sum)
    end if

    r2 = 0.0_DOUBLE
    r3 = 0.0_DOUBLE
    do i=1, mesh%n_elems
      if( .not. mesh%elem(i)%is_ghost ) then
        r2 = max(r2, abs((new_sol(1, i) - sol(1, i))/sol(1, i)))
        do j=1, 5
          r3 = max(r3, abs((new_sol(j, i) - sol(j, i)/sol(j, i))))
        end do
      end if
    end do

    if( num_procs > 1) then
      call MPI_ALLREDUCE(MPI_IN_PLACE, r2, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD, mpi_ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE, r3, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD, mpi_ierr)
    end if

    call cpu_time(tcur)
    write(funit, *) iter+restart_iter, t+restart_time, r1, r2, r3, &
      restart_cpu_time+tcur-tstart, r4, dt
  end subroutine write_dat_residual

  subroutine update_new_sol_with_flux(mesh, new_sol, sol, rhs, sum_lambda, residu, dt, &
      num_procs)
    use mpi
    use global_data_module, only: gamma, local_time_step, second_order
    implicit none

    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(inout) :: new_sol
    real(kind=DOUBLE), dimension(mesh%n_elems), intent(inout) :: residu
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(in) :: sol
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(inout) :: rhs
    real(kind=DOUBLE), dimension(mesh%n_elems), intent(inout) :: sum_lambda
    real(kind=DOUBLE), intent(inout) :: dt
    integer(kind=ENTIER), intent(in) :: num_procs

    integer(kind=ENTIER) :: i

    !$OMP PARALLEL DO
    do i = 1, mesh%n_elems
      new_sol(:, i) = sol(:, i) - rhs(:, i)*dt/mesh%elem(i)%volume
    end do
  end subroutine update_new_sol_with_flux

  subroutine multi_point_pressure(sol_l, sol_r, sol_w_l, sol_w_r, lr_flux, p_nodal, &
      lambda_l, lambda_r, sl, sr)
    implicit none

    real(kind=DOUBLE), intent(in) :: p_nodal
    real(kind=DOUBLE), dimension(5), intent(in) :: sol_l, sol_r
    real(kind=DOUBLE), dimension(5), intent(in) :: sol_w_l, sol_w_r
    real(kind=DOUBLE), intent(inout) :: lambda_l, lambda_r
    real(kind=DOUBLE), dimension(5, 2), intent(inout) :: lr_flux
    real(kind=DOUBLE), intent(inout) :: sl, sr

    real(kind=DOUBLE) :: rhol, rhor, ul, vl, wl, ur, vr, wr
    real(kind=DOUBLE) :: pl, pr, rhol_et, rhor_et
    real(kind=DOUBLE) :: ul_et, vl_et, wl_et, ur_et, vr_et, wr_et
    real(kind=DOUBLE) :: el, er
    real(kind=DOUBLE) :: u_bar, p_bar, pl_bar, pr_bar
    real(kind=DOUBLE), dimension(5) :: fl, fr
    real(kind=DOUBLE), dimension(5) :: sol_l_et, sol_r_et
    real(kind=DOUBLE) :: eig_val

    rhol = sol_w_l(1)
    ul = sol_w_l(2)
    vl = sol_w_l(3)
    wl = sol_w_l(4)
    pl = sol_w_l(5)
    el = sol_l(5)/rhol

    rhor = sol_w_r(1)
    ur = sol_w_r(2)
    vr = sol_w_r(3)
    wr = sol_w_r(4)
    pr = sol_w_r(5)
    er = sol_r(5)/rhor

    fl(:) = ul*sol_l(:) + (/0.0_DOUBLE, pl, 0.0_DOUBLE, 0.0_DOUBLE, pl*ul/)
    fr(:) = ur*sol_r(:) + (/0.0_DOUBLE, pr, 0.0_DOUBLE, 0.0_DOUBLE, pr*ur/)

    p_bar = (lambda_r*pl + lambda_l*pr - lambda_r*lambda_l*(ur - ul))/(lambda_r + lambda_l)
    u_bar = (lambda_l*ul + lambda_r*ur)/(lambda_l + lambda_r) &
      - (pr - pl)/(lambda_r + lambda_l)

    ul_et = ul - (p_nodal - pl)/lambda_l
    rhol_et = 1.0_DOUBLE/(1.0_DOUBLE/rhol + (u_bar - ul)/lambda_l)
    vl_et = vl
    wl_et = wl
    pl_bar = p_nodal

    sol_l_et(1) = rhol_et
    sol_l_et(2) = rhol_et*ul_et
    sol_l_et(3) = rhol_et*vl_et
    sol_l_et(4) = rhol_et*wl_et
    sol_l_et(5) = rhol_et*(el + (pl*ul - pl_bar*u_bar)/lambda_l)

    ur_et = ur + (p_nodal - pr)/lambda_r
    rhor_et = 1.0_DOUBLE/(1.0_DOUBLE/rhor + (ur - u_bar)/lambda_r)
    vr_et = vr
    wr_et = wr
    pr_bar = p_nodal

    sol_r_et(1) = rhor_et
    sol_r_et(2) = rhor_et*ur_et
    sol_r_et(3) = rhor_et*vr_et
    sol_r_et(4) = rhor_et*wr_et
    sol_r_et(5) = rhor_et*(er + (pr_bar*u_bar - pr*ur)/lambda_r)

    if (rhol_et < 0.0_DOUBLE .or. rhor_et < 0.0_DOUBLE) then
      print *, "Negative specific volume !", rhol_et, rhor_et
      error stop
    end if

    sl = ul - lambda_l/rhol
    sr = ur + lambda_r/rhor

    eig_val = max(eig_val, abs(sl), abs(sr))

    lr_flux(:, 1) = 0.5_DOUBLE*(fl + fr) - 0.5_DOUBLE* &
      (abs(sl)*(sol_l_et - sol_l) + &
      abs(u_bar)*(sol_r_et - sol_l_et) + &
      abs(sr)*(sol_r - sol_r_et))

    lr_flux(:, 2) = 0.5_DOUBLE*(fl + fr) - 0.5_DOUBLE* &
      (abs(sl)*(sol_l_et - sol_l) + &
      abs(u_bar)*(sol_r_et - sol_l_et) + &
      abs(sr)*(sol_r - sol_r_et))

    lr_flux(:, 2) = -lr_flux(:, 2)
  end subroutine multi_point_pressure

  subroutine compute_lambdas_and_solve_nodal_pressure(mesh, id_vert, sol_w_lr, lambda, p_bars, p_node)
    use global_data_module, only: gamma, lambda_l_eq_lambda_r, bc_style
    use subfv_linear_solver_module
    implicit none

    type(mesh_type), intent(in) :: mesh
    integer(kind=ENTIER), intent(in) :: id_vert
    real(kind=DOUBLE), dimension(5, 2, mesh%vert(id_vert)%n_sub_faces_neigh), &
      intent(in) :: sol_w_lr
    real(kind=DOUBLE), dimension(2, mesh%vert(id_vert)%n_sub_faces_neigh), &
      intent(inout) :: lambda
    real(kind=DOUBLE), dimension(mesh%vert(id_vert)%n_sub_faces_neigh), &
      intent(inout) :: p_bars
    real(kind=DOUBLE), intent(inout) :: p_node

    integer(kind=ENTIER) :: j, id_sub_face, id_face, le, re
    integer(kind=ENTIER) :: iter, n_iter
    real(kind=DOUBLE) :: rhol, ul, pl, al, lambda_l
    real(kind=DOUBLE) :: rhor, ur, pr, ar, lambda_r
    real(kind=DOUBLE) :: s1, s2

    n_iter = 10
    p_node = 0.0_DOUBLE
    do iter = 1, n_iter
      s1 = 0.0_DOUBLE
      s2 = 0.0_DOUBLE

      !Compute initial lambdas
      do j = 1, mesh%vert(id_vert)%n_sub_faces_neigh
        id_sub_face = mesh%vert(id_vert)%sub_face_neigh(j)
        id_face = mesh%sub_face(id_sub_face)%mesh_face

        le = mesh%sub_face(id_sub_face)%left_elem_neigh
        re = mesh%sub_face(id_sub_face)%right_elem_neigh

        rhol = sol_w_lr(1, 1, j)
        ul = sol_w_lr(2, 1, j)
        pl = sol_w_lr(5, 1, j)
        al = sqrt(gamma*pl/rhol)

        rhor = sol_w_lr(1, 2, j)
        ur = sol_w_lr(2, 2, j)
        pr = sol_w_lr(5, 2, j)
        ar = sqrt(gamma*pr/rhor)

        lambda_l = lambda(1, j)
        lambda_r = lambda(2, j)

        if (iter == 1) then
          lambda_l = max(lambda_l, al*rhol, sqrt(rhol*max(0.0_DOUBLE, pr - pl)), -rhol*(ur - ul))
          lambda_r = max(lambda_r, ar*rhor, sqrt(rhor*max(0.0_DOUBLE, pl - pr)), -rhor*(ur - ul))

          p_bars(j) = (lambda_r*pl + lambda_l*pr - lambda_r*lambda_l*(ur - ul))/(lambda_r + lambda_l)
          lambda_l = max(lambda_l, (al*rhol + 1.5_DOUBLE*sqrt(max(0.0_DOUBLE, (p_bars(j) - pl)*rhol))))
          lambda_r = max(lambda_r, (ar*rhor + 1.5_DOUBLE*sqrt(max(0.0_DOUBLE, (p_bars(j) - pr)*rhor))))

          ! u_bar = (lambda_l*ul + lambda_r*ur - (pr - pl))/(lambda_r + lambda_l)
          ! lambda_l = max(lambda_l, al*rhol*(1.0_DOUBLE + 1.5_DOUBLE*max(0.0_DOUBLE, -(u_bar - ul)/al)))
          ! lambda_r = max(lambda_r, ar*rhor*(1.0_DOUBLE + 1.5_DOUBLE*max(0.0_DOUBLE, (u_bar - ur)/ar)))
        else
          if (bc_style == 0) then
            lambda_l = max(lambda_l, al*rhol + 1.5_DOUBLE*sqrt(max(0.0_DOUBLE, (p_node - pl)*rhol)))
            lambda_r = max(lambda_r, ar*rhor + 1.5_DOUBLE*sqrt(max(0.0_DOUBLE, (p_node - pr)*rhor)))
          else if (bc_style == 1) then
            if (is_wall(re)) then
              lambda_l = max(lambda_l, (al*rhol + 1.5_DOUBLE*sqrt(max(0.0_DOUBLE, (p_bars(j) - pl)*rhol))))
              lambda_r = max(lambda_r, (ar*rhor + 1.5_DOUBLE*sqrt(max(0.0_DOUBLE, (p_bars(j) - pr)*rhor))))
            else
              lambda_l = max(lambda_l, al*rhol + 1.5_DOUBLE*sqrt(max(0.0_DOUBLE, (p_node - pl)*rhol)))
              lambda_r = max(lambda_r, ar*rhor + 1.5_DOUBLE*sqrt(max(0.0_DOUBLE, (p_node - pr)*rhor)))
            end if
          else
            print *, "bad bc style"
            error stop
          end if
        end if

        if (lambda_l_eq_lambda_r) then
          lambda_l = max(lambda_l, lambda_r)
          lambda_r = lambda_l
        end if

        p_bars(j) = (lambda_r*pl + lambda_l*pr - lambda_r*lambda_l*(ur - ul))/(lambda_r + lambda_l)

        if (bc_style == 0) then
          s1 = s1 + (1.0_DOUBLE/lambda_r + 1.0_DOUBLE/lambda_l)*mesh%sub_face(id_sub_face)%area*p_bars(j)
          s2 = s2 + (1.0_DOUBLE/lambda_r + 1.0_DOUBLE/lambda_l)*mesh%sub_face(id_sub_face)%area
        else if (bc_style == 1) then
          if (is_wall(re)) then
            ! s1 = s1 + (1.0_DOUBLE/lambda_r + 1.0_DOUBLE/lambda_l)*mesh%sub_face(id_sub_face)%area*p_bars(j)
            s1 = s1 + (1.0_DOUBLE/lambda_r + 1.0_DOUBLE/lambda_l)*mesh%sub_face(id_sub_face)%area &
              *(lambda_r*pl + lambda_r*lambda_l*ul)/(lambda_r + lambda_l)
            s2 = s2 + (1.0_DOUBLE/lambda_l)*mesh%sub_face(id_sub_face)%area
          else
            s1 = s1 + (1.0_DOUBLE/lambda_r + 1.0_DOUBLE/lambda_l)*mesh%sub_face(id_sub_face)%area*p_bars(j)
            s2 = s2 + (1.0_DOUBLE/lambda_r + 1.0_DOUBLE/lambda_l)*mesh%sub_face(id_sub_face)%area
          end if
        else
          print *, "Bad bc_style !"
          error stop
        end if

        lambda(1, j) = lambda_l
        lambda(2, j) = lambda_r
      end do

      p_node = s1/s2
    end do
  end subroutine compute_lambdas_and_solve_nodal_pressure

  subroutine put_vortex_ring(mesh, sol, c, v)
    use global_data_module, only : pi
    implicit none

    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(inout) :: sol
    real(kind=DOUBLE), dimension(3) :: c, v

    integer(kind=ENTIER) :: i
    real(kind=DOUBLE) :: w, r, theta, velocity, phi, distXY, largeR

    w = 0.1_DOUBLE
    largeR = 0.25_DOUBLE
    do i=1, mesh%n_elems
      distXY = sqrt((mesh%elem(i)%coord(1)-c(1))**2 + &
        (mesh%elem(i)%coord(2)-c(2))**2) + 1e-14
      r = sqrt((distXY-largeR)**2 + (mesh%elem(i)%coord(3)-c(3))**2)
      phi = atan2(mesh%elem(i)%coord(2)-c(2), mesh%elem(i)%coord(1)-c(1)+1e-14)
      theta = asin((mesh%elem(i)%coord(3)-c(3))/(r+1e-14_DOUBLE))

      if(distXY < largeR) theta = pi - theta

      if(r < w) then
        velocity = 10.0_DOUBLE*r
      else 
        velocity = max(0.0_DOUBLE, 2.0_DOUBLE - 10.0_DOUBLE*r)
      end if

      if(distXY > 1e-12) then
        sol(2, i) = sol(2, i) - sin(theta)*cos(phi)*velocity/distXY
        sol(3, i) = sol(3, i) - sin(theta)*sin(phi)*velocity/distXY
        sol(4, i) = sol(4, i) + cos(theta)*velocity/distXY
      end if
    end do
  end subroutine put_vortex_ring

  subroutine restart_from_vtu_file(mesh, sol, restart_file, me, num_procs)
    implicit none

    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), dimension(5, mesh%n_elems), intent(inout) :: sol
    integer(kind=ENTIER), intent(in) :: me, num_procs
    character(len=*), intent(in) :: restart_file

    integer(kind=ENTIER) :: i, vtkin
    character(len=255) :: text, me_str, mpi_restart_file

    if (num_procs > 1) then
      write (me_str, *) me
      mpi_restart_file = trim(adjustl(me_str))//"_"//trim(adjustl(restart_file))
      open (newunit=vtkin, file=trim(mpi_restart_file), status="old")
    else
      open (newunit=vtkin, file=trim(restart_file), status="old")
    end if

    read (vtkin, '(a)') text
    do while (" <DataArray type='Float64' Name='Density' format='ascii' NumberOfComponents='1'>" /= text(:80))
      read (vtkin, '(a)') text
    end do
    do i = 1, mesh%n_elems
      if (.not. mesh%elem(i)%is_ghost) then
        read (vtkin, *) sol(1, i)
      end if
    end do

    do while (" <DataArray type='Float64' Name='Momentum' format='ascii' NumberOfComponents='3'>" /= text(:81))
      read (vtkin, '(a)') text
    end do
    do i = 1, mesh%n_elems
      if (.not. mesh%elem(i)%is_ghost) then
        read (vtkin, *) sol(2, i), sol(3, i), sol(4, i)
      end if
    end do

    do while (" <DataArray type='Float64' Name='Energy' format='ascii' NumberOfComponents='1'>" /= text(:79))
      read (vtkin, '(a)') text
    end do
    do i = 1, mesh%n_elems
      if (.not. mesh%elem(i)%is_ghost) then
        read (vtkin, *) sol(5, i)
        sol(5, i) = sol(5, i) * sol(1, i)
      end if
    end do

    close (vtkin)
  end subroutine restart_from_vtu_file

  function sol_primit_sod_angle(coord, x_sod, sod_angle, sol_left, sol_right) result(w)
    use global_data_module, only: pi
    implicit none

    real(kind=DOUBLE), dimension(3), intent(in) :: coord
    real(kind=DOUBLE), intent(in) :: sod_angle, x_sod
    real(kind=DOUBLE), dimension(5), intent(in) :: sol_left, sol_right
    real(kind=DOUBLE), dimension(5) :: w

    real(kind=DOUBLE) :: sod_angle_rad
    real(kind=DOUBLE), dimension(3,3) :: mat_rot
    real(kind=DOUBLE), dimension(3) :: coord_rot

    sod_angle_rad = sod_angle/360._DOUBLE*2._DOUBLE*pi
    mat_rot = 0.0_DOUBLE
    mat_rot(1, 1) = cos(sod_angle_rad)
    mat_rot(2, 1) = -sin(sod_angle_rad)
    mat_rot(1, 2) = sin(sod_angle_rad)
    mat_rot(2, 2) = cos(sod_angle_rad)
    mat_rot(3, 3) = 1.0_DOUBLE
    coord_rot = matmul(mat_rot, coord)

    if (coord_rot(1) < x_sod) then
      w = sol_left
    else
      w = sol_right
    end if
    w(2:4) = matmul(transpose(mat_rot), w(2:4))
  end function sol_primit_sod_angle

  subroutine compute_error_shear(mesh, sol)
    use global_data_module, only: x_sod, sod_angle, &
      sol_left, sol_right, error_2d, error_2d_h
    implicit none

    type(mesh_type), intent(in) :: mesh
    real(kind=DOUBLE), dimension(5, mesh%n_elems) :: sol

    integer(kind=ENTIER) :: i
    real(kind=DOUBLE) :: error, volume
    real(kind=DOUBLE), dimension(5) :: wexact, wsol

    error = 0.0_DOUBLE
    volume = 0.0_DOUBLE
    do i = 1, mesh%n_elems
      if ( abs(mesh%elem(i)%coord(1)) < 0.8 .and. &
        abs(mesh%elem(i)%coord(2)) < 0.8 ) then
        wexact = sol_primit_sod_angle(mesh%elem(i)%coord, x_sod, sod_angle, &
          sol_left, sol_right)
        wsol = conserv_to_primit(sol(:, i))
        error = error &
          + mesh%elem(i)%volume*maxval(abs(wsol(2:3) - wexact(2:3)))
        volume = volume + mesh%elem(i)%volume
      end if
    end do

    if( error_2d ) then
      print *, "Error Shear: ", sqrt((volume/error_2d_h)/mesh%n_elems), &
        ",", error
    else
      print*, "Error 2D not set"
      error stop
    end if
  end subroutine compute_error_shear
end module euler_module
