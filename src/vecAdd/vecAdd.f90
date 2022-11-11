module vecFunctions
    implicit none
        real :: t0,t1
        contains
            ! Vector addition routine using OpenACC
            subroutine vecAdd_ACC(n,a,b,c)
                implicit none
                    integer(8), intent(in) :: n
                    real(4), intent(in)    :: a(n), b(n)
                    real(4), intent(out)   :: c(n)
                    integer(8)             :: i
                    call CPU_TIME(t0)
                    ! Add appropriate acc directive
                    do i = 1, n
                        c(i) = a(i) + b(i)
                    end do
                    ! End acc region
                    call CPU_TIME(t1)
            end subroutine vecAdd_ACC
            ! Vector addition routine using OpenMP
            subroutine vecAdd_OMP(n,a,b,c)
                implicit none
                    integer(8), intent(in) :: n
                    real(4), intent(in)    :: a(n), b(n)
                    real(4), intent(out)   :: c(n)
                    integer(8)             :: i
                    call CPU_TIME(t0)
                    !$omp parallel do
                    do i = 1, n
                        c(i) = a(i) + b(i)
                    end do
                    !$omp end parallel do
                    call CPU_TIME(t1)
            end subroutine vecAdd_OMP
end module vecFunctions

! Driver program to test the routines
program test
    use omp_lib
    use openacc
    use vecFunctions
    implicit none
        integer(8), parameter   :: nIter = 10
        integer(8), parameter   :: n = 4*(512**3)   ! Array size
        integer(8)              :: i
        real(4)                 :: resOMP, resACC   ! Results
        real(4),    allocatable :: a(:), b(:), c(:) ! Arrays
        print*, "Nummber of running threads := ",omp_get_max_threads()
        print*, "Number of available devices := ",acc_get_num_devices(acc_device_nvidia)
        print*, "Required memory := ",(n*4*3)/1024.0/1024.0/1024.0," GB"
        allocate(a(n), b(n), c(n))
        !$omp parallel do
        do i = 1, n
            a(i) = 1.0
            b(i) = 2.0
            c(i) = 0.0
        end do
        !$omp end parallel do
        ! Call and time OMP routine
        do i = 1, nIter
            call vecAdd_OMP(n, a, b, c)
            print*, "OMP time = ", t1-t0, "s"
        end do
        resOMP = c(1)
        ! Call and time ACC routine
        a(:) = 1.0
        b(:) = 2.0
        c(:) = 0.0
        do i = 1, nIter
            call vecAdd_ACC(n, a, b, c)
            print*, "ACC time = ", t1-t0, "s"
        end do
        resACC = c(1)
        ! Compare results
        if (resOMP == resACC) then
            print*, "Results match"
        else
            print*, "Results do not match"
            stop 1
        end if
end program test