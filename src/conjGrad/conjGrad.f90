! Small example detailing a more reallistic application of openACC
! to a conjugate gradient solver with a "dense" matrix. Compares 2
! different implementations of the CG algorithm against the OpenMP
! version. Missing a version using shared memory.
! Suggestions:
!     1. Test V0 and V1 without ACC (remove -gpu:* and -acc compilation options from Makefile)
!     2. Remove ACC from sections of the code to see the effects.
!     3. On V1, remove the atomic update and see the effects.
!     4. Test different vector/matrix sizes.
!     5. Find out why the first run of V0 is slower than the rest.
!     6. On V0, try substituting the ACC kernels by the DIY version.
!     7. Check the nsys report after running "make run".

module solver
    implicit none
        PUBLIC :: conjGrad_OMP, conjGrad_V0, conjGrad_V1
        PRIVATE
        integer(4), parameter :: maxIter = 100
        real(4)   , parameter :: tol = 1.0e-8
        integer(4)            :: iter
        integer(8)            :: i, j
        real(4)               :: alpha, beta, rho, err, aux
        contains
            subroutine conjGrad_V0(A,b,x,n)
                implicit none
                    integer(8), intent(in)  :: n
                    real(4)   , intent(in)  :: A(n,n), b(n)
                    real(4)   , intent(out) :: x(n)
                    real(4)                 :: r(n), p(n), Ap(n)
                    ! Initialize auxliary variables
                    alpha = 0.0
                    beta  = 0.0
                    rho   = 0.0
                    err   = 0.0
                    aux   = 0.0
                    ! Initialize auxilliary arrays
                    !$acc parallel loop
                    do i = 1, n
                        x(i) = 0.0
                        r(i) = 0.0
                        p(i) = 0.0
                        Ap(i) = 0.0
                    end do
                    !$acc end parallel loop
                    ! Compute initial residual
                    !$acc kernels
                    do i = 1, n
                        aux = 0.0
                        do j = 1, n
                            aux = aux + A(i,j)*x(j)
                        end do
                        r(i) = b(i) - aux
                    end do
                    !$acc end kernels
                    print*, "V0, r0*r0 := ", dot_product(r,r)
                    ! Start conjugate gradient iterations
                    outer:do iter = 1,maxIter
                        !$acc parallel loop
                        do i = 1, n
                            p(i) = r(i) + beta*p(i)
                        end do
                        !$acc end parallel loop
                        rho = 0.0
                        !$acc parallel loop reduction(+:rho)
                        do i = 1, n
                            rho = rho + r(i)*r(i)
                        end do
                        !$acc end parallel loop
                        !$acc kernels
                        do i = 1, n
                            aux = 0.0
                            do j = 1,n
                                aux = aux + A(i,j)*p(j)
                            end do
                            Ap(i) = aux
                        end do
                        !$acc end kernels
                        aux = 0.0
                        !$acc parallel loop reduction(+:aux)
                        do i = 1, n
                            aux = aux + p(i)*Ap(i)
                        end do
                        !$acc end parallel loop
                        alpha = rho/aux
                        !$acc parallel loop
                        do i = 1, n
                            x(i) = x(i) + alpha*p(i)
                            r(i) = r(i) - alpha*Ap(i)
                        end do
                        !$acc end parallel loop
                        err = 0.0
                        !$acc parallel loop reduction(+:err)
                        do i = 1, n
                            err = err + r(i)*r(i)
                        end do
                        !$acc end parallel loop
                        err = sqrt(err)
                        if (err < tol) then
                            write(*,*) "Conjugate gradient converged in ", iter, " iterations"
                            write(*,*) "Residual norm: ", err
                            exit outer
                        end if
                        write(*,*) "Iteration: ", iter, " Residual norm: ", err
                        beta = (err**2)/rho
                    end do outer
            end subroutine conjGrad_V0

            subroutine conjGrad_V1(A,b,x,n)
                implicit none
                    integer(8), intent(in)  :: n
                    real(4)   , intent(in)  :: A(n,n), b(n)
                    real(4)   , intent(out) :: x(n)
                    real(4)                 :: r(n), p(n), Ap(n)
                    ! Initialize auxliary variables
                    alpha = 0.0
                    beta  = 0.0
                    rho   = 0.0
                    err   = 0.0
                    aux   = 0.0
                    ! Initialize auxilliary arrays
                    !$acc parallel loop
                    do i = 1, n
                        x(i) = 0.0
                        r(i) = 0.0
                        p(i) = 0.0
                        Ap(i) = 0.0
                    end do
                    !$acc end parallel loop
                    ! Compute initial residual
                    !$acc parallel loop
                    do j = 1, n
                        !$acc loop vector
                        do i = 1, n
                            !$acc atomic update
                            r(i) = r(i) + (A(i,j)*x(j))
                            !$acc end atomic
                        end do
                    end do
                    !$acc end parallel loop
                    !$acc parallel loop
                    do i = 1, n
                        r(i) = b(i) - r(i)
                    end do
                    !$acc end parallel loop
                    print*, "V1, r0*r0 := ", dot_product(r,r)
                    ! Start conjugate gradient iterations
                    outer:do iter = 1, maxIter
                        !$acc parallel loop
                        do i = 1, n
                            p(i) = r(i) + beta*p(i)
                        end do
                        !$acc end parallel loop
                        rho = 0.0
                        !$acc parallel loop reduction(+:rho)
                        do i = 1, n
                            rho = rho + r(i)*r(i)
                        end do
                        !$acc end parallel loop
                        !$acc parallel loop
                        do i = 1, n
                            Ap(i) = 0.0
                        end do
                        !$acc end parallel loop
                        !$acc parallel loop
                        do j = 1, n
                            !$acc loop vector
                            do i = 1, n
                                !$acc atomic update
                                Ap(i) = Ap(i) + A(i,j)*p(j)
                                !$acc end atomic
                            end do
                        end do
                        !$acc end parallel loop
                        aux = 0.0
                        !$acc parallel loop reduction(+:aux)
                        do i = 1, n
                            aux = aux + p(i)*Ap(i)
                        end do
                        !$acc end parallel loop
                        alpha = rho/aux
                        !$acc parallel loop
                        do i = 1, n
                            x(i) = x(i) + alpha*p(i)
                            r(i) = r(i) - alpha*Ap(i)
                        end do
                        !$acc end parallel loop
                        err = 0.0
                        !$acc parallel loop reduction(+:err)
                        do i = 1, n
                            err = err + r(i)*r(i)
                        end do
                        !$acc end parallel loop
                        err = sqrt(err)
                        if (err < tol) then
                            write(*,*) "Conjugate gradient converged in ", iter, " iterations"
                            write(*,*) "Residual norm: ", err
                            exit outer
                        end if
                        write(*,*) "Iteration: ", iter, " Residual norm: ", err
                        beta = (err**2)/rho
                    end do outer
            end subroutine conjGrad_V1

            subroutine conjGrad_OMP(A,b,x,n)
                implicit none
                    integer(8), intent(in)  :: n
                    real(4)   , intent(in)  :: A(n,n), b(n)
                    real(4)   , intent(out) :: x(n)
                    real(4)                 :: r(n), p(n), Ap(n)
                    ! Initialize auxliary variables
                    alpha = 0.0
                    beta  = 0.0
                    rho   = 0.0
                    err   = 0.0
                    aux   = 0.0
                    ! Initialize auxilliary arrays
                    !$omp parallel do
                    do i = 1, n
                        x(i) = 0.0
                        r(i) = 0.0
                        p(i) = 0.0
                        Ap(i) = 0.0
                    end do
                    !$omp end parallel do
                    ! Compute initial residual
                    !$omp parallel do private(aux)
                    do i = 1, n
                        aux = 0.0
                        !$omp simd reduction(+:aux)
                        do j = 1, n
                            aux = aux + A(i,j)*x(j)
                        end do
                        r(i) = b(i) - aux
                    end do
                    !$omp end parallel do
                    ! Start conjugate gradient iterations
                    outer:do iter = 1,maxIter
                        !$omp parallel do
                        do i = 1, n
                            p(i) = r(i) + beta*p(i)
                        end do
                        !$omp end parallel do
                        rho = 0.0
                        !$omp parallel do reduction(+:rho)
                        do i = 1, n
                            rho = rho + r(i)*r(i)
                        end do
                        !$omp end parallel do
                        !$omp parallel do private(aux)
                        do i = 1, n
                            aux = 0.0
                            !$omp simd reduction(+:aux)
                            do j = 1,n
                                aux = aux + A(i,j)*p(j)
                            end do
                            Ap(i) = aux
                        end do
                        !$omp end parallel do
                        aux = 0.0
                        !$omp parallel do reduction(+:aux)
                        do i = 1, n
                            aux = aux + p(i)*Ap(i)
                        end do
                        !$omp end parallel do
                        alpha = rho/aux
                        !$omp parallel do
                        do i = 1, n
                            x(i) = x(i) + alpha*p(i)
                            r(i) = r(i) - alpha*Ap(i)
                        end do
                        !$omp end parallel do
                        err = 0.0
                        !$omp parallel do reduction(+:err)
                        do i = 1, n
                            err = err + r(i)*r(i)
                        end do
                        !$omp end parallel do
                        err = sqrt(err)
                        if (err < tol) then
                            write(*,*) "Conjugate gradient converged in ", iter, " iterations"
                            write(*,*) "Residual norm: ", err
                            exit outer
                        end if
                        write(*,*) "Iteration: ", iter, " Residual norm: ", err
                        beta = (err**2)/rho
                    end do outer
            end subroutine conjGrad_OMP

            subroutine conjGrad_OMPGPU(A,b,x,n)
                implicit none
                    integer(8), intent(in)  :: n
                    real(4)   , intent(in)  :: A(n,n), b(n)
                    real(4)   , intent(out) :: x(n)
                    real(4)                 :: r(n), p(n), Ap(n)
                    ! Initialize auxliary variables
                    alpha = 0.0
                    beta  = 0.0
                    rho   = 0.0
                    err   = 0.0
                    aux   = 0.0
                    ! Initialize auxilliary arrays
                    !$omp target teams distribute parallel do
                    do i = 1, n
                        x(i) = 0.0
                        r(i) = 0.0
                        p(i) = 0.0
                        Ap(i) = 0.0
                    end do
                    !$omp end target teams distribute parallel do
                    ! Compute initial residual
                    !$omp parallel do private(aux)
                    do i = 1, n
                        aux = 0.0
                        !$omp simd reduction(+:aux)
                        do j = 1, n
                            aux = aux + A(i,j)*x(j)
                        end do
                        r(i) = b(i) - aux
                    end do
                    !$omp end parallel do
                    ! Start conjugate gradient iterations
                    outer:do iter = 1,maxIter
                        !$omp parallel do
                        do i = 1, n
                            p(i) = r(i) + beta*p(i)
                        end do
                        !$omp end parallel do
                        rho = 0.0
                        !$omp parallel do reduction(+:rho)
                        do i = 1, n
                            rho = rho + r(i)*r(i)
                        end do
                        !$omp end parallel do
                        !$omp parallel do private(aux)
                        do i = 1, n
                            aux = 0.0
                            !$omp simd reduction(+:aux)
                            do j = 1,n
                                aux = aux + A(i,j)*p(j)
                            end do
                            Ap(i) = aux
                        end do
                        !$omp end parallel do
                        aux = 0.0
                        !$omp parallel do reduction(+:aux)
                        do i = 1, n
                            aux = aux + p(i)*Ap(i)
                        end do
                        !$omp end parallel do
                        alpha = rho/aux
                        !$omp parallel do
                        do i = 1, n
                            x(i) = x(i) + alpha*p(i)
                            r(i) = r(i) - alpha*Ap(i)
                        end do
                        !$omp end parallel do
                        err = 0.0
                        !$omp parallel do reduction(+:err)
                        do i = 1, n
                            err = err + r(i)*r(i)
                        end do
                        !$omp end parallel do
                        err = sqrt(err)
                        if (err < tol) then
                            write(*,*) "Conjugate gradient converged in ", iter, " iterations"
                            write(*,*) "Residual norm: ", err
                            exit outer
                        end if
                        write(*,*) "Iteration: ", iter, " Residual norm: ", err
                        beta = (err**2)/rho
                    end do outer
            end subroutine conjGrad_OMPGPU
end module solver

program test
#ifdef USE_OMP
    use omp_lib
#endif
#ifdef USE_ACC
    use openacc
#endif
    use solver
        implicit none
            integer(8), parameter   :: n =5000 ! Matrix and array size (change at will)
            integer(8), parameter   :: nTimes = 10
            real(4)   , parameter   :: pi = 3.141592653589793238462643383279502884197169399375105820974944592307816406286
            real(4)   , allocatable :: A(:,:), b(:), x(:), xOMP(:), xV0(:), xV1(:)
            integer(8)              :: i, j, iTimes
            real(4)                 :: t0, t1, avgTime
#ifdef USE_OMP
            print*, "Number of threads: ", omp_get_max_threads()
#endif
#ifdef USE_ACC
            print*, "Number of devices: ", acc_get_num_devices(acc_device_nvidia)
#endif
            print*, "Memory required: ", (4.0*n*n + 4.0*2*n + 4.0*3*n)/1024.0/1024.0/1024.0, " GB"
            ! Allocate arrays
            allocate(A(n,n), b(n), x(n))
            allocate(xOMP(n), xV0(n), xV1(n))
            ! Generate the matrix A
            A(:,:) = 0.0
            A(1,n) = ((1.0/real(n,4))/6.0)*1.0
            A(n,1) = ((1.0/real(n,4))/6.0)*1.0
            !$omp parallel do private(i)
            do i = 1, n
                A(i,i)   = ((1.0/real(n,4))/6.0)*4.0
            end do
            !$omp end parallel do
            !$omp parallel do private(i)
            do i = 1, n-1
                A(i,i+1) = ((1.0/n)/6.0)*1.0
                A(i+1,i) = ((1.0/n)/6.0)*1.0
            end do
            !$omp end parallel do
            ! Generate array b
            do i = 1, n
                b(i) = 20.0*sin(2.0*pi*i/n)
            end do
            ! Solve the system nTimes
            open(1,file="cg_omp.txt",status="replace")
            do iTimes = 1, nTimes
                ! Solve the system
                call CPU_TIME(t0)
                call conjGrad_OMP(A,b,x,n)
                call CPU_TIME(t1)
                write(*,'(a,f8.4,a)') "OMP Time: ", t1-t0, "s"
                write(1,'(i0,f8.4)') iTimes, t1-t0
            end do
            close(1)
            xOMP(:) = x(:)
            ! Solve the system nTimes
            open(1,file="cg_V0.txt",status="replace")
            do iTimes = 1, nTimes
                ! Solve the system
                call CPU_TIME(t0)
                call conjGrad_V0(A,b,x,n)
                call CPU_TIME(t1)
                write(*,'(a,f8.4,a)') "V0 Time: ", t1-t0, "s"
                write(1,'(i0,f8.4)') iTimes, t1-t0
            end do
            close(1)
            xV0(:) = x(:)
            ! Solve the system nTimes
            open(1,file="cg_V1.txt",status="replace")
            do iTimes = 1, nTimes
                ! Solve the system
                call CPU_TIME(t0)
                call conjGrad_V1(A,b,x,n)
                call CPU_TIME(t1)
                write(*,'(a,f8.4,a)') "V1 Time: ", t1-t0, "s"
                write(1,'(i0,f8.4)') iTimes, t1-t0
            end do
            close(1)
            xV1(:) = x(:)

            ! Write the results
            open(1,file="cg_results.txt",status="replace")
            do i = 1, n
                write(1,*) 2*pi*i/n, xOMP(i), xV0(i), xV1(i)
            end do
            close(1)
end program test