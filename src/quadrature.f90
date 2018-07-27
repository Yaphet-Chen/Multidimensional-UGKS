module quadrature
contains
    subroutine gaussian_maxwell(npoints, abscissas, weight) ! For Maxwell polynomial exp(-x^2)*x^p (p=0 only, half-space).
        !---------------------------------------------------------------------------------------------------------------------------
        ! References:
        ! 1. A Gaussian quadrature procedure for use in the solution of the Boltzmann equation and related problems
        ! 2. Spectral Methods in Chemistry and Physics: Applications to Kinetic Theory and Quantum Mechanics
        ! 3. Half-Range Generalized Hermite Polynomials and the Related Gaussian Quadratures
        !---------------------------------------------------------------------------------------------------------------------------
        use iso_fortran_env, only : REAL64
        implicit none
        integer, intent(in) :: npoints ! Number of points for half-space.
        real(kind=REAL64), allocatable, dimension(:), intent(out) :: abscissas
        real(kind=REAL64), allocatable, dimension(:), intent(out) :: weight
        real(kind=REAL64), allocatable, dimension(:,:) :: eigenvector
        real(kind=REAL64), allocatable, dimension(:) :: workspace
        real(kind=REAL64), allocatable, dimension(:) :: alpha
        real(kind=REAL64), allocatable, dimension(:) :: beta
        real(kind=REAL64), allocatable, dimension(:) :: gn
        real(kind=REAL64), parameter :: PI = 4.0*atan(1.0_REAL64)
        real(kind=REAL64) :: residual
        real(kind=REAL64) :: gn_old
        real(kind=REAL64) :: y
        integer :: iter
        integer :: err
        integer :: p
        integer :: i

        allocate(abscissas(npoints))
        allocate(weight(npoints))
        allocate(alpha(npoints))
        allocate(beta(npoints))
        allocate(gn(max(51,npoints+1)))
        allocate(eigenvector(npoints, npoints))
        allocate(workspace(max(1,2*npoints-2)))

        !---------------------------------------------------------------------------------------------------------------------------
        ! Initialize alpha and beta in the tridiagonal matrix, see Ref.1.
        !---------------------------------------------------------------------------------------------------------------------------
        p = 0
        alpha(1) = 1.0_REAL64/sqrt(PI)
        alpha(2) = 2.0_REAL64/sqrt(PI)/(PI-2.0)
        beta(1) = 0.0
        beta(2) = 0.0+(p+1)/2.0_REAL64-alpha(1)**2-beta(1)

        !---------------------------------------------------------------------------------------------------------------------------
        ! Compute alpha and beta in the tridiagonal matrix, see Ref.3. For other methods, see Ref.2.
        !---------------------------------------------------------------------------------------------------------------------------
        ! Starting values.
        gn(1) = beta(1)-(0.0+p/2.0_REAL64)/6.0
        gn(2) = beta(2)-(1.0+p/2.0_REAL64)/6.0

        ! Initial estimation of gn.
        do i = 3, min(npoints, 9)
            y = 2*(i-2)+p
            gn(i) = (y+1)/3.0_REAL64-gn(i-1)-&
                    ((y/6.0_REAL64-gn(i-1))**2-p**2/16.0_REAL64)**2/((y-1)/3.0_REAL64-gn(i-1)-gn(i-2))/(y/12.0_REAL64+gn(i-1))**2
        end do
        do i = min(npoints, 9)+1, size(gn)
            gn(i) = compute_gn(i-1)
        end do

        ! Fix point iteration of gn (accuracy improvement). Newton's method has faster convergence speed.
        iter = 0
        do while(.true.)
            residual = 0.0
            do i = 3, size(gn)-1
                gn_old = gn(i)
                gn(i) = iterate_gn(i)
                residual = max(residual, abs(gn_old-gn(i)))
            end do
            iter = iter+1
            if (residual <= 1.0E-16 .or. iter >= 50) then
                exit
            end if
        end do

        ! Compute alpha and sqrt(beta).
        do i = 1, npoints
            alpha(i) = sqrt((2*(i-1)+p+1)/3.0_REAL64-gn(i+1)-gn(i))
            beta(i) = sqrt(gn(i+1)+(i+p/2.0_REAL64)/6.0) ! Start from beta(2) since beta(1) is not used in the matrix.
        end do

        !---------------------------------------------------------------------------------------------------------------------------
        ! Compute abscissas and weights, see Ref.1.
        !---------------------------------------------------------------------------------------------------------------------------
        call dstev('V', npoints, alpha, beta, eigenvector, npoints, workspace, err) ! Eigenvalue and eigenvector by lapack.

        abscissas = alpha
        weight = gamma((p+1)/2.0_REAL64)/2.0_REAL64*eigenvector(1,:)**2
    contains
        elemental function compute_gn(n) result(g)
            integer, intent(in) :: n
            real(kind=REAL64) :: g
            real(kind=REAL64) :: c0, c1, c2, c3, y

            c0 = 1.0_REAL64/36.0-p**2/8.0_REAL64
            c1 = 23.0_REAL64/432.0-11.0_REAL64/48.0*p**2+3.0_REAL64/32.0*p**4
            c2 = 1189.0_REAL64/2592.0-409.0_REAL64/192.0*p**2+75.0_REAL64/64.0*p**4+9.0_REAL64/64.0*p**6
            c3 = 196057.0_REAL64/20736.0&
                -153559.0_REAL64/3456.0*p**2+7111.0_REAL64/256.0*p**4+639.0_REAL64/128.0*p**6+135.0_REAL64/512.0*p**8
            y = 2*n+p
            g = c0/y+c1/y**3+c2/y**5+c3/y**7
        end function compute_gn

        function iterate_gn(n) result(g)
            integer, intent(in) :: n
            real(kind=REAL64) :: g
            real(kind=REAL64) :: y, a, b, c, d, e

            y = 2*(n-1)+p
            a = gn(n+1)-(y+1)/3.0_REAL64
            b = gn(n-1)-(y-1)/3.0_REAL64
            c = y/12.0_REAL64
            d = y/6.0_REAL64
            e = p**2/16.0_REAL64
            g = gn(n)**3*(-4*d-2*c-(a+b))+gn(i)**2*(6*d**2-2*e**2-c**2-(a+b)*2*c-a*b)-2*(e*d)**2+e**4+d**4-a*b*c**2
            g = g/((a+b)*c**2+2*a*b*c+4*d**3-4*e**2*d)
        end function iterate_gn
    end subroutine gaussian_maxwell
end module quadrature

program main
    use iso_fortran_env, only : REAL64
    use quadrature
    implicit none
    real(kind=REAL64), allocatable, dimension(:) :: abscissas
    real(kind=REAL64), allocatable, dimension(:) :: weight
    character(len=:), allocatable :: command
    integer :: length
    integer :: npoint
    integer :: i

    call get_command_argument(1, length=length)
    allocate(character(len=length) :: command)
    call get_command_argument(1, command)
    read(command, *) npoint

    call gaussian_maxwell(npoint, abscissas, weight)
    do i = 1, size(abscissas)
        write(*,*) abscissas(i), weight(i)
    end do
end program main

! vim: set ft=fortran ff=unix tw=132:
