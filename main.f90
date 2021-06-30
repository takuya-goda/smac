! 11.SMAC法
! 角柱流れ

program fortran

    implicit none

! 変数
    integer, parameter :: nxc = 200
    integer, parameter :: nyc = 100
    integer, parameter :: nxd = nxc + 1
    integer, parameter :: nyd = nyc + 1
    real(8), parameter :: lx = 0.2d0
    real(8), parameter :: ly = 0.1d0
    real(8), parameter :: dx = 0.001d0
    real(8), parameter :: dy = 0.001d0
    real(8), parameter :: uin = 0.01d0
    real(8), parameter :: lt = 5d0

    real(8), parameter :: Re = 1000d0 !=uin*lx/visc
    real(8), parameter :: Co = 0.1d0 !クーラン数

    real(8) :: dt
    real(8) :: dens, visc, arufa, err_total
    real(8) :: start_time, end_time
    integer :: i, j, time, nt, outputstep
    integer :: ic, jc
    real(8) :: x, y
    real(8), dimension(0:nxd + 1, 0:nyc + 1) :: u
    real(8), dimension(0:nxc + 1, 0:nyd + 1) :: v
    real(8), dimension(0:nxd + 1, 0:nyc + 1) :: u_aux
    real(8), dimension(0:nxc + 1, 0:nyd + 1) :: v_aux
    real(8), dimension(0:nxc + 1, 0:nyc + 1) :: p
    real(8), dimension(0:nxc + 1, 0:nyc + 1) :: phi
    real(8), dimension(0:nxc + 1, 0:nyc + 1) :: dp
    real(8), dimension(nxc, nyc) :: theta

    integer, parameter :: prizm_left = 90
    integer, parameter :: prizm_right = 110
    integer, parameter :: prizm_bottom = 40
    integer, parameter :: prizm_top = 60

! 初期設定
    dt = Co*dx/uin
    nt = nint(lt/dt)

    dens = 1.0e+3
! visc = 1d0/1e+6
    visc = uin*lx/Re !0.01*0.2/1000
    err_total = 1d0/1e+12
    arufa = 1.99d0

    u(:, :) = uin
    v(:, :) = 0d0
    u_aux(:, :) = 0d0
    v_aux(:, :) = 0d0
    p(:, :) = 0d0
    phi(:, :) = 0d0
    theta(:, :) = 0d0

    nt = 75000
    outputstep = 100

    ! u(0,:) = uin*2d0 - u(1,:)
    ! u(:, 0) = uin*2d0 - u(:, 1)
    ! u(:, nyc + 1) = uin*2d0 - u(:, nyc)
    u(0, :) = uin
    u(:, 0) = u(:, 1)
    u(:, nyc + 1) = u(:, nyc)

    v(0, :) = v(1, :)
    v(:, 0) = 0d0
    v(:, nyd) = 0d0

    call setBoundary(u, v, uin, prizm_left, prizm_right, prizm_top, prizm_bottom)

! メインルーチン
    call cpu_time(start_time)
    do time = 1, nt
        ! do time = 1, 30
        write (*, *) time
        call computeAuxiallyVelocity(u_aux, v_aux, u, v, p, dx, dy, dt, visc, uin, &
                                     prizm_left, prizm_right, prizm_top, prizm_bottom)
        call computeDivergenceAuxiallyVelocity(theta, u_aux, v_aux, dx, dy, dt, &
                                               prizm_left, prizm_right, prizm_top, prizm_bottom)
        call computePressurePoisson(p, phi, dp, u, v, dx, dy, theta, err_total, arufa, dens, visc, dt, &
                                    prizm_left, prizm_right, prizm_top, prizm_bottom)
        call computeVelocity(u, v, u_aux, v_aux, dt, dx, dy, dens, phi, uin, &
                             prizm_left, prizm_right, prizm_top, prizm_bottom)
        if (mod(time, outputstep) == 0) then
            call output(u, v, p, phi, time, theta)
            ! elseif (time < 110 .and. mod(time, 10) == 0) then
            !     call output(u, v, p, phi, time, theta)
        elseif (time > 4900 .and. time < 5000 .and. mod(time, 10) == 0) then
            call output(u, v, p, phi, time, theta)
        end if
    end do
    call cpu_time(end_time)
    write (*, *) "elapsed time = ", end_time - start_time

    open (unit=30, file="elapsed.txt")
    write (30, *) end_time - start_time
    close (30)

contains

    subroutine output(u, v, p, phi, time, theta)
        implicit none
        real(8), intent(in) :: u(0:, 0:)
        real(8), intent(in) :: v(0:, 0:)
        real(8), intent(in) :: p(0:, 0:)
        real(8), intent(in) :: phi(0:, 0:)
        integer, intent(in) :: time
        real(8), intent(in) :: theta(:, :)

        integer :: i, j
        real(8) :: uout, vout
        character*128 filename1, filename2, filename3

        write (filename1, '("pressure/pres",i5.5,".txt")') time
        write (filename2, '("phi/phi",i5.5,".txt")') time
        write (filename3, '("theta/theta",i5.5,".txt")') time
        open (unit=10, file=filename1)
        open (unit=20, file=filename2)
        open (unit=30, file=filename3)
        do j = 1, nyc
        do i = 1, nxc
            uout = (u(i, j) + u(i + 1, j))/2d0
            vout = (v(i, j) + v(i, j + 1))/2d0
            if (prizm_left <= i .and. i <= prizm_right .and. prizm_bottom <= j .and. j <= prizm_top) then
                uout = 0d0
                vout = 0d0
            elseif (i == prizm_right .and. prizm_bottom <= j .and. j <= prizm_top) then
                uout = 0d0
            elseif (prizm_left <= i .and. i <= prizm_right .and. j == prizm_top) then
                vout = 0d0
            end if
            write (10, '(2i8,3f23.14)') i, j, p(i, j), uout, vout
            write (20, '(2i8,3f23.14)') i, j, phi(i, j), uout, vout
            write (30, '(2i8,3f23.14)') i, j, theta(i, j), uout, vout
        end do
        end do
        close (10)
        close (20)
        close (30)

    end subroutine

    subroutine setBoundary(u, v, uin, prizm_left, prizm_right, prizm_top, prizm_bottom)
        implicit none
        real(8), intent(inout) :: u(0:, 0:)
        real(8), intent(inout) :: v(0:, 0:)
        real(8), intent(in) :: uin
        integer, intent(in) :: prizm_left, prizm_right, prizm_top, prizm_bottom

        ! u(0, :) = uin*2d0 - u(1, :)
        ! u(:, 0) = uin*2d0 - u(:, 1)
        ! u(:, nyc + 1) = uin*2d0 - u(:, nyc)
        !&<
        u(0, :) = uin
        u(:, 0) = u(:, 1)
        u(:, nyc + 1) = u(:, nyc)

        ! v(0, :) = 0d0
        v(0, :) = v(1, :)
        v(:, 1) = 0d0
        v(:, nyd) = 0d0

        ! 流出境界条件
        u(nxd + 1, :) = u(nxd, :)
        v(nxc + 1, :) = v(nxc, :)

        ! 角柱表面
        u(prizm_left + 1:prizm_right, prizm_bottom) = -u(prizm_left + 1:prizm_right, prizm_bottom - 1)
        u(prizm_left + 1:prizm_right, prizm_top   ) = -u(prizm_left + 1:prizm_right, prizm_top    + 1)
        u(prizm_left + 1            , prizm_bottom:prizm_top) = -u(prizm_left  - 1, prizm_bottom:prizm_top)
        u(prizm_right               , prizm_bottom:prizm_top) = -u(prizm_right + 2, prizm_bottom:prizm_top)

        v(prizm_left            , prizm_bottom + 1:prizm_top) = -v(prizm_left  - 1, prizm_bottom + 1:prizm_top)
        v(prizm_right           , prizm_bottom + 1:prizm_top) = -v(prizm_right + 1, prizm_bottom + 1:prizm_top)
        v(prizm_left:prizm_right, prizm_bottom + 1) = -v(prizm_left:prizm_right, prizm_bottom - 1)
        v(prizm_left:prizm_right, prizm_top       ) = -v(prizm_left:prizm_right, prizm_top    + 2)

        ! 角柱と流体の境界
        u(prizm_left     , prizm_bottom:prizm_top) = 0d0
        u(prizm_right + 1, prizm_bottom:prizm_top) = 0d0
        v(prizm_left:prizm_right, prizm_bottom ) = 0d0
        v(prizm_left:prizm_right, prizm_top + 1) = 0d0

        ! 角柱内部
        u(prizm_left + 2:prizm_right - 1, prizm_bottom + 1:prizm_top - 1) = 0d0
        v(prizm_left + 1:prizm_right - 1, prizm_bottom + 2:prizm_top - 1) = 0d0
        !&>

    end subroutine

    subroutine computeAuxiallyVelocity(u_aux, v_aux, u, v, p, dx, dy, dt, visc, uin, &
                                       prizm_left, prizm_right, prizm_top, prizm_bottom)
        implicit none
        real(8), intent(inout) :: u_aux(0:, 0:)
        real(8), intent(inout) :: v_aux(0:, 0:)
        real(8), intent(in) :: u(0:, 0:)
        real(8), intent(in) :: v(0:, 0:)
        real(8), intent(in) :: p(0:, 0:)
        real(8), intent(in) :: dx, dy, dt, visc, uin
        integer, intent(in) :: prizm_left, prizm_right, prizm_top, prizm_bottom

        !&<
        do jc = 1, nyc
        do i = 1, nxd
            j = jc
            ic = i
            if (prizm_left <= i .and. i <= prizm_right .and. prizm_bottom <= jc .and. jc <= prizm_top) cycle
            if (i == prizm_right + 1 .and. prizm_bottom <= jc .and. jc <= prizm_top) cycle !角柱の右辺上
            u_aux(i, jc) = u(i, jc) - dt*((u(i  , jc ) + u(i-1 , jc ))/2d0*(u(i  , jc  ) - u(i-1, jc  ))/dx &
                                        + (u(i+1, jc ) + u(i   , jc ))/2d0*(u(i+1, jc  ) - u(i  , jc  ))/dx &
                                        + (v(ic , j+1) + v(ic-1, j+1))/2d0*(u(i  , jc+1) - u(i  , jc  ))/dy &
                                        + (v(ic , j  ) + v(ic-1, j  ))/2d0*(u(i  , jc  ) - u(i  , jc-1))/dy)/2d0 &
                         + dt*visc*((u(i+1, jc  ) - 2d0*u(i, jc) + u(i-1, jc  ))/dx**2 &
                                  + (u(i  , jc+1) - 2d0*u(i, jc) + u(i  , jc-1))/dy**2) &
                         - dt*(p(ic, jc) - p(ic - 1, jc))/dens/dx
        end do
        end do

        do j = 1, nyd
        do ic = 1, nxc
            jc = j
            i = ic
            if (prizm_left <= ic .and. ic <= prizm_right .and. prizm_bottom <= j .and. j <= prizm_top) cycle
            if (prizm_left <= ic .and. ic <= prizm_right .and. j == prizm_top + 1) cycle !角柱の上辺上
            v_aux(ic, j) = v(ic, j) - dt*((u(i  , jc ) + u(i  , jc-1))/2d0*(v(ic  , j  ) - v(ic-1, j  ))/dx &
                                        + (u(i+1, jc ) + u(i+1, jc-1))/2d0*(v(ic+1, j  ) - v(ic  , j  ))/dx &
                                        + (v(ic , j+1) + v(ic , j   ))/2d0*(v(ic  , j+1) - v(ic  , j  ))/dy &
                                        + (v(ic , j  ) + v(ic , j-1 ))/2d0*(v(ic  , j  ) - v(ic  , j-1))/dy)/2d0 &
                         + dt*visc*((v(ic+1, j  ) - 2d0*v(ic, j) + v(ic-1, j  ))/dx**2 &
                                  + (v(ic  , j+1) - 2d0*v(ic, j) + v(ic  , j-1))/dy**2) &
                         - dt*(p(ic, jc) - p(ic, jc - 1))/dens/dy
        end do
        end do
        !&>

        call setBoundary(u_aux, v_aux, uin, prizm_left, prizm_right, prizm_top, prizm_bottom)

    end subroutine

    subroutine computeDivergenceAuxiallyVelocity(theta, u_aux, v_aux, dx, dy, dt, prizm_left, prizm_right, prizm_top, prizm_bottom)
        implicit none
        real(8), intent(inout) :: theta(:, :)
        real(8), intent(in) :: u_aux(0:, 0:)
        real(8), intent(in) :: v_aux(0:, 0:)

        real(8), intent(in) :: dt, dx, dy
        integer, intent(in) :: prizm_left, prizm_right, prizm_top, prizm_bottom

        do jc = 1, nyc
        do ic = 1, nxc
            j = jc
            i = ic
            ! if (prizm_left <= ic .and. ic <= prizm_right .and. prizm_bottom <= jc .and. jc <= prizm_top) cycle ! 修正ポイント3
            theta(ic, jc) = (u_aux(i + 1, jc) - u_aux(i, jc))/dx + (v_aux(ic, j + 1) - v_aux(ic, j))/dy
        end do
        end do

    end subroutine

    subroutine computePressurePoisson(p, phi, dp, u, v, dx, dy, theta, err_total, arufa, dens, visc, dt, &
                                      prizm_left, prizm_right, prizm_top, prizm_bottom)
        implicit none
        real(8), intent(inout) :: p(0:, 0:)
        real(8), intent(inout) :: phi(0:, 0:)
        real(8), intent(inout) :: dp(0:, 0:)
        real(8), intent(in) :: u(0:, 0:)
        real(8), intent(in) :: v(0:, 0:)
        real(8), intent(in) :: theta(:, :)

        real(8), intent(in) :: dx, dy, err_total, arufa, dens, visc, dt
        integer, intent(in) :: prizm_left, prizm_right, prizm_top, prizm_bottom
        real(8) :: err, err_phi, err_dp, err_p, test_sum
        integer count

        dp(:, :) = 0d0
        ! phi(0:, 0:) = 0d0
        phi(:, :) = 0d0 !修正ポイント1
        err = 1d0
        count = 0
        test_sum = 0d0

        !&<
        do while (err > err_total)
            ! count = count + 1
            err_p = 0d0
            err_dp = 0d0
            err_phi = 0d0
            do jc = 1, nyc
            do ic = 1, nxc
                j = jc
                i = ic
                if (prizm_left <= ic .and. ic <= prizm_right .and. prizm_bottom <= jc .and. jc <= prizm_top) cycle
                dp(ic, jc) = ((dy**2d0*(phi(ic + 1, jc) + phi(ic - 1, jc)) + &
                               dx**2d0*(phi(ic, jc + 1) + phi(ic, jc - 1)) - &
                               dx**2d0*dy**2d0*theta(ic, jc)*dens/dt)/(2d0*(dx**2d0 + dy**2d0))) - phi(ic, jc)
                phi(ic, jc) = phi(ic, jc) + arufa*dp(ic, jc)
                err_dp = err_dp + abs(dp(ic, jc))
                err_phi = err_phi + abs(phi(ic, jc))
            end do
            end do
            ! 境界条件　角柱
            phi(prizm_left , prizm_bottom:prizm_top) = phi(prizm_left  - 1, prizm_bottom:prizm_top)
            phi(prizm_right, prizm_bottom:prizm_top) = phi(prizm_right + 1, prizm_bottom:prizm_top)
            phi(prizm_left:prizm_right, prizm_bottom) = phi(prizm_left:prizm_right, prizm_bottom - 1)
            phi(prizm_left:prizm_right, prizm_top   ) = phi(prizm_left:prizm_right, prizm_top    + 1)

            ! 境界条件　流入流出　圧力勾配０
            ! phi(0    , :    ) = phi(1  , :  )
            ! phi(nxc+1, :    ) = phi(nxc, :  )
            ! phi(:    , 0    ) = phi(:  , 1  )
            ! phi(:    , nyc+1) = phi(:  , nyc)

            if (err_phi < 1.0d-20) then
                err_phi = 1d0
            end if
            err = err_dp/err_phi

            ! print *,err,err_dp,err_phi

        end do
        !&>

        ! p(0:, 0:) = p(0:, 0:) + phi(0:, 0:)
        p(:, :) = p(:, :) + phi(:, :) !修正ポイント2
        ! print *, sum(p), sum(phi), test_sum, err_phi
        ! print *, theta(89, 60), theta(89, 50)

    end subroutine

    subroutine computeVelocity(u, v, u_aux, v_aux, dt, dx, dy, dens, phi, uin, prizm_left, prizm_right, prizm_top, prizm_bottom)
        implicit none
        real(8), intent(inout) :: u(0:, 0:)
        real(8), intent(inout) :: v(0:, 0:)
        real(8), intent(inout) :: phi(0:, 0:)
        real(8), intent(in) :: u_aux(0:, 0:)
        real(8), intent(in) :: v_aux(0:, 0:)

        integer, intent(in) :: prizm_left, prizm_right, prizm_top, prizm_bottom
        real(8), intent(in) :: dt, dx, dy, dens, uin

        do jc = 1, nyc
        do i = 1, nxd
            j = jc
            ic = i
            if (prizm_left <= i .and. i <= prizm_right .and. prizm_bottom <= jc .and. jc <= prizm_top) cycle
            if (prizm_bottom <= jc .and. jc <= prizm_top .and. i == prizm_right + 1) cycle
            u(i, jc) = u_aux(i, jc) - dt*(phi(ic, jc) - phi(ic - 1, jc))/dx/dens
        end do
        end do

        do j = 1, nyd
        do ic = 1, nxc
            jc = j
            i = ic
            if (prizm_left <= ic .and. ic <= prizm_right .and. prizm_bottom <= j .and. j <= prizm_top) cycle
            if (prizm_left <= ic .and. ic <= prizm_right .and. j == prizm_top + 1) cycle
            v(ic, j) = v_aux(ic, j) - dt*(phi(ic, jc) - phi(ic, jc - 1))/dy/dens
        end do
        end do

        call setBoundary(u, v, uin, prizm_left, prizm_right, prizm_top, prizm_bottom)

    end subroutine

end program
