! Small test program to compare OpenMP and OpenACC, and introduce the use of
! OpenACC in Fortran. The program performs a vector addition of two arrays,
! and compares the results, running the kernels several times to get an
! average execution time.
!
! Suggestions:
!     1. Compare usage of different ACC directives (magic, DIY)
!     2. Compare execution times for different array sizes
!     3. Compare execution times for different number of threads
!     4. Comment lines 88-92 and compare execution times
!     5. Compile and run without "managed" option (see Makefile)
!     6. Make kernel operations more complex (e.g. c(i) = cos(a(i))**2 + sin(b(i))**2)


module vecFunctions
    implicit none
        real :: t0, t1, tAvg_omp, tAvg_acc
        contains
            ! Vector addition routine using OpenACC
            subroutine vecAdd_ACC(n,a,b,c)
                implicit none
                    integer(8), intent(in) :: n
                    real(4), intent(in)    :: a(n), b(n)
                    real(4), intent(out)   :: c(n)
                    integer(8)             :: i
                    call CPU_TIME(t0)
                    ! Add appropriate acc directive (magic or DIY)
                    do i = 1, n
                        c(i) = a(i) + b(i)
                    end do
                    ! End acc region (!$acc end [same as above])
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

    ! Modules and libraries to use
    use omp_lib      ! Library of OpenMP instructions (optional)
    use openacc      ! Library of OpenACC instructions (optional)
    use vecFunctions ! The module (above) containing the functions

    ! Variable declaration
    implicit none
    integer(8), parameter   :: nIter = 10       ! Number of kernel repetitions (Modify at will)
    integer(8), parameter   :: n = 4*(512**3)   ! Array size (Modify at will)
    integer(8)              :: i                ! Array index
    real(4)                 :: resOMP, resACC   ! Results
    real(4),    allocatable :: a(:), b(:), c(:) ! Arrays (must be allocated)

    ! Program starts here
    print*, "Nummber of running threads := ",omp_get_max_threads()                   ! How many omp threads are active
    print*, "Number of available devices := ",acc_get_num_devices(acc_device_nvidia) ! change to acc_device_radeon for AMD
    print*, "Required memory := ",(n*4*3)/1024.0/1024.0/1024.0," GB"                 ! (3 arrrrays of float type)

    ! Generate arrays
    allocate(a(n), b(n), c(n)) ! Allocate (reserve) memory for arrays (a,b,c) with size n
    !$omp parallel do
    do i = 1, n
        a(i) = 1.0
        b(i) = 2.0
        c(i) = 0.0
    end do
    !$omp end parallel do

    ! Call and time OMP routine
    tAvg_omp = 0.0
    do i = 1, nIter
        call vecAdd_OMP(n, a, b, c)
        print*, "OMP time = ", t1-t0, "s"
        tAvg_omp = tAvg_omp + (t1-t0)
    end do
    resOMP = c(1)
    tAvg_omp = tAvg_omp/real(nIter,4)

    ! Call and time ACC routine
    !$acc kernels
    a(:) = 1.0
    b(:) = 2.0
    c(:) = 0.0
    !$acc end kernels
    tAvg_acc = 0.0
    do i = 1, nIter
        call vecAdd_ACC(n, a, b, c)
        print*, "ACC time = ", t1-t0, "s"
        tAvg_acc = tAvg_acc + (t1-t0)
    end do
    resACC = c(1)
    tAvg_acc = tAvg_acc/real(nIter,4)

    ! Compare results
    if (resOMP == resACC) then
        print*, "Results match"
        print*, "Average OMP time = ", tAvg_omp, "s"
        print*, "Average ACC time = ", tAvg_acc, "s"
    else
        print*, "Results do not match"
        stop 1
    end if

end program test