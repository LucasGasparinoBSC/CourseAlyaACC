module mesh
	implicit none
		public :: generateMesh
		contains
			pure subroutine generateMesh(nElem,pOrder,nNodes,nPoints,listConnec,xyz,hElem)
				implicit none
					integer(4), intent(in)  :: pOrder, nNodes
					integer(8), intent(in)  :: nElem, nPoints
					integer(8), intent(out) :: listConnec(nElem,nNodes)
					integer(4)              :: iNode
					integer(8)              :: iElem, iPoint
					real(4)   , intent(out) :: xyz(nPoints)
					real(4)   , intent(out) :: hElem
					real(4)   , parameter   :: xMin = 0.0, xMax = 1.0
					real(4)                 :: aux, x0
					! Initialize the output arrays to 0
					listConnec(:,:) = 0
					xyz(:)          = 0.0
					! Compute eelement size hElem
					hElem = (xMax - xMin) / real(nElem,8)
					! Generate linear mesh data: connectivity and coordinates
					! Connectivity
					do iElem = 1, nElem
						do iNode = 1, 2
							listConnec(iElem,iNode) = iElem + iNode - 1
						end do
					end do
					! Coordinates
					do iElem = 1, nElem
						do iNode = 1, 2
							iPoint = listConnec(iElem,iNode)
							xyz(iPoint) = xMin + (iPoint-1)*((xMax-xMin)/nElem)
						end do
					end do
					! If order > 1, generate extra points and connectivity
					if (pOrder > 1) then
						! Connectivity
						do iElem = 1, nElem
							do iNode = 3, nNodes
								listConnec(iElem,iNode) = nElem+(nNodes-2)*(iElem-1)+iNode-1
							end do
						end do
						! Coordinates
						aux = hElem / real(pOrder,8)
						do iElem = 1, nElem
							do iNode = 3, nNodes
								iPoint = listConnec(iElem,iNode)
								x0 = xyz(listConnec(iElem,1))
								xyz(iPoint) = x0 + (iNode-2)*aux
							end do
						end do
					end if
			end subroutine generateMesh
end module mesh

module interpolation
	implicit none
		public :: lagrangeGrid, evalLagrangePoly, evalLagrangePolyDeriv
		contains
			pure subroutine lagrangeGrid(pOrder,xi)
				implicit none
					integer(4), intent(in)  :: pOrder
					integer(4)              :: i
					real(4)   , intent(out) :: xi(pOrder+1)
					real(4)                 :: h
					xi(1) = -1.0
					xi(2) =  1.0
					h = 2.0/real(pOrder,8)
					if (pOrder > 1) then
						do i = 3, pOrder+1
							xi(i) = xi(1) + real(i-2,8)*h
						end do
					end if
			end subroutine lagrangeGrid
			pure subroutine evalLagrangePoly(pOrder,x,lpXi)
				implicit none
					integer(4), intent(in)  :: pOrder
					real(4)   , intent(in)  :: x
					real(4)   , intent(out) :: lpXi(pOrder+1)
					real(4)                 :: xi(pOrder+1)
					integer(4)              :: i, j
					call lagrangeGrid(pOrder,xi)
					do i = 1, pOrder+1
						lpXi(i) = 1.0
						do j = 1, pOrder+1
							if (j /= i) then
								lpXi(i) = lpXi(i) * (x-xi(j))/(xi(i)-xi(j))
							end if
						end do
					end do
			end subroutine evalLagrangePoly
			pure subroutine evalLagrangePolyDeriv(pOrder,x,dlpXi)
				implicit none
					integer(4), intent(in)  :: pOrder
					real(4)   , intent(in)  :: x
					real(4)   , intent(out) :: dlpXi(pOrder+1)
					real(4)                 :: xi(pOrder+1)
					integer(4)              :: i, j, k
					real(4)				    :: aux
					call lagrangeGrid(pOrder,xi)
					do i = 1, pOrder+1
						dlpXi(i) = 0.0
						do j = 1, pOrder+1
							if (j /= i) then
								aux = 1.0
								do k = 1, pOrder+1
									if (k /= i .and. k /= j) then
										aux = aux * (x-xi(k))/(xi(i)-xi(k))
									end if
								end do
								dlpXi(i) = dlpXi(i) + aux/(xi(i)-xi(j))
							end if
						end do
					end do
			end subroutine evalLagrangePolyDeriv
end module interpolation

module FEM
	use interpolation
	implicit none
		private
		public :: quadratureTable, evalElemInfo
		contains
			subroutine quadratureTable(pOrder,xgp,wgp)
				implicit none
					integer(4), intent(in)  :: pOrder
					real(4)   , intent(out) :: xgp(pOrder+1), wgp(pOrder+1)
					! Gauss-Legendre quadrature
					if (pOrder == 1) then
						xgp(1) = -0.577350269189626
						xgp(2) =  0.577350269189626
						wgp(1) =  1.0
						wgp(2) =  1.0
					else if (pOrder == 2) then
						xgp(1) = -0.774596669241483
						xgp(2) =  0.0
						xgp(3) =  0.774596669241483
						wgp(1) =  0.555555555555556
						wgp(2) =  0.888888888888889
						wgp(3) =  0.555555555555556
					else if (pOrder == 3) then
						xgp(1) = -0.861136311594053
						xgp(2) = -0.339981043584856
						xgp(3) =  0.339981043584856
						xgp(4) =  0.861136311594053
						wgp(1) =  0.347854845137454
						wgp(2) =  0.652145154862546
						wgp(3) =  0.652145154862546
						wgp(4) =  0.347854845137454
					else if (pOrder == 4) then
						xgp(1) = -0.906179845938664
						xgp(2) = -0.538469310105683
						xgp(3) =  0.0
						xgp(4) =  0.538469310105683
						xgp(5) =  0.906179845938664
						wgp(1) =  0.236926885056189
						wgp(2) =  0.478628670499366
						wgp(3) =  0.568888888888889
						wgp(4) =  0.478628670499366
						wgp(5) =  0.236926885056
					else if (pOrder == 5) then
						xgp(1) = -0.932469514203152
						xgp(2) = -0.661209386466265
						xgp(3) = -0.238619186083197
						xgp(4) =  0.238619186083197
						xgp(5) =  0.661209386466265
						xgp(6) =  0.932469514203152
						wgp(1) =  0.171324492379170
						wgp(2) =  0.360761573048139
						wgp(3) =  0.467913934572691
						wgp(4) =  0.467913934572691
						wgp(5) =  0.360761573048139
						wgp(6) =  0.171324492379170
					else if (pOrder == 6) then
						xgp(1) = -0.949107912342759
						xgp(2) = -0.741531185599394
						xgp(3) = -0.405845151377397
						xgp(4) =  0.0
						xgp(5) =  0.405845151377397
						xgp(6) =  0.741531185599394
						xgp(7) =  0.949107912342759
						wgp(1) =  0.129484966168870
						wgp(2) =  0.279705391489277
						wgp(3) =  0.381830050505119
						wgp(4) =  0.417959183673469
						wgp(5) =  0.381830050505119
						wgp(6) =  0.279705391489277
						wgp(7) =  0.129484966168870
					else if (pOrder == 7) then
						xgp(1) = -0.960289856497536
						xgp(2) = -0.796666477413627
						xgp(3) = -0.525532409916329
						xgp(4) = -0.183434642495650
						xgp(5) =  0.183434642495650
						xgp(6) =  0.525532409916329
						xgp(7) =  0.796666477413627
						xgp(8) =  0.960289856497536
						wgp(1) =  0.101228536290376
						wgp(2) =  0.222381034453374
						wgp(3) =  0.313706645877887
						wgp(4) =  0.362683783378362
						wgp(5) =  0.362683783378362
						wgp(6) =  0.313706645877887
						wgp(7) =  0.222381034453374
						wgp(8) =  0.101228536290376
					else if (pOrder == 8) then
						xgp(1) = -0.968160239507626
						xgp(2) = -0.836031107326636
						xgp(3) = -0.613371432700590
						xgp(4) = -0.324253423403809
						xgp(5) =  0.0
						xgp(6) =  0.324253423403809
						xgp(7) =  0.613371432700590
						xgp(8) =  0.836031107326636
						xgp(9) =  0.968160239507626
						wgp(1) =  0.081274388361574
						wgp(2) =  0.180648160694857
						wgp(3) =  0.260610696402935
						wgp(4) =  0.312347077040003
						wgp(5) =  0.330239355001260
						wgp(6) =  0.312347077040003
						wgp(7) =  0.260610696402935
						wgp(8) =  0.180648160694857
						wgp(9) =  0.081274388361574
					else if (pOrder == 9) then
						xgp(1) = -0.973906528517172
						xgp(2) = -0.865063366688985
						xgp(3) = -0.679409568299024
						xgp(4) = -0.433395394129247
						xgp(5) = -0.148874338981631
						xgp(6) =  0.148874338981631
						xgp(7) =  0.433395394129247
						xgp(8) =  0.679409568299024
						xgp(9) =  0.865063366688985
						xgp(10) =  0.973906528517172
						wgp(1) =  0.066671344308688
						wgp(2) =  0.149451349150581
						wgp(3) =  0.219086362515982
						wgp(4) =  0.269266719309996
						wgp(5) =  0.295524224714753
						wgp(6) =  0.295524224714753
						wgp(7) =  0.269266719309996
						wgp(8) =  0.219086362515982
						wgp(9) =  0.149451349150581
						wgp(10) =  0.066671344308688
					else if (pOrder == 10) then
						xgp(1) = -0.978228658146057
						xgp(2) = -0.887062599768095
						xgp(3) = -0.730152005574049
						xgp(4) = -0.519096129206812
						xgp(5) = -0.269543155952345
						xgp(6) =  0.0
						xgp(7) =  0.269543155952345
						xgp(8) =  0.519096129206812
						xgp(9) =  0.730152005574049
						xgp(10) =  0.887062599768095
						xgp(11) =  0.978228658146057
						wgp(1) =  0.055668567116174
						wgp(2) =  0.125580369464905
						wgp(3) =  0.186290210927734
						wgp(4) =  0.233193764591990
						wgp(5) =  0.262804544510247
						wgp(6) =  0.272925086777901
						wgp(7) =  0.262804544510247
						wgp(8) =  0.233193764591990
						wgp(9) =  0.186290210927734
						wgp(10) =  0.125580369464905
						wgp(11) =  0.055668567116174
					else if (pOrder > 10) then
						write(*,*) "Error: pOrder > 6 not implemented"
						stop 1
					end if
			end subroutine quadratureTable
			subroutine evalElemInfo(pOrder,wgp,Ngp,dNgp)
				implicit none
					integer(4), intent(in)  :: pOrder
					integer(4)              :: iGauss
					real(4)   , intent(out) :: wgp(porder+1), Ngp(porder+1,pOrder+1), dNgp(porder+1,pOrder+1)
					real(4)                 :: xgp(pOrder+1)
					! Generate Gauss-Legendre quadrature table
					call quadratureTable(pOrder,xgp,wgp)
					! Generate shape functions and derivatives at all Gauss points
					do iGauss = 1,pOrder+1
						call evalLagrangePoly(pOrder,xgp(iGauss),Ngp(iGauss,:))
						call evalLagrangePolyDeriv(pOrder,xgp(iGauss),dNgp(iGauss,:))
					end do
			end subroutine evalElemInfo
end module FEM

module elemOps
	implicit none
		public :: convec_omp, convec_acc
		contains
			subroutine convec_omp(nElem,nNodes,nGauss,nPoints,listConnec,wgp,Ngp,dNgp,Je,He,phi,Rconvec)
				implicit none
					integer(4), intent(in)  :: nNodes, nGauss
					integer(8), intent(in)  :: nElem, nPoints
					integer(8), intent(in)  :: listConnec(nElem,nNodes)
					real(4)   , intent(in)  :: wgp(nGauss), Ngp(nGauss,nNodes), dNgp(nGauss,nNodes)
					real(4)   , intent(in)  :: Je, He, phi(nPoints)
					real(4)   , intent(out) :: Rconvec(nPoints)
					integer(4)              :: iElem, iNode, iGauss, iPoint
					real(4)                 :: phiElem(nNodes), aux
					! Initialize residual
					Rconvec = 0.0
					! Loop over elements
					!$omp parallel do default(none) private(iElem,iNode,iGauss,iPoint,phiElem,aux) shared(nElem,nNodes,nGauss,nPoints,listConnec,wgp,Ngp,dNgp,Je,He,phi,Rconvec)
					do iElem = 1,nElem
						! Extract nodal values
						do iNode = 1,nNodes
							iPoint = listConnec(iElem,iNode)
							phiElem(iNode) = phi(iPoint)
						end do
						! Loop over Gauss points
						do iGauss = 1,nGauss
							! Compute derivatives
							aux = dot_product(dNgp(iGauss,:),phiElem(:))
							! Compute residual
							do iNode = 1,nNodes
								iPoint = listConnec(iElem,iNode)
								Rconvec(iPoint) = Rconvec(iPoint) + (wgp(iGauss) * Ngp(iGauss,iNode) * aux)
							end do
						end do
					end do
					!$omp end parallel do
			end subroutine convec_omp

			subroutine convec_acc(nElem,nNodes,nGauss,nPoints,listConnec,wgp,Ngp,dNgp,Je,He,phi,Rconvec)
				implicit none
					integer(4), intent(in)  :: nNodes, nGauss
					integer(8), intent(in)  :: nElem, nPoints
					integer(8), intent(in)  :: listConnec(nElem,nNodes)
					real(4)   , intent(in)  :: wgp(nGauss), Ngp(nGauss,nNodes), dNgp(nGauss,nNodes)
					real(4)   , intent(in)  :: Je, He, phi(nPoints)
					real(4)   , intent(out) :: Rconvec(nPoints)
					integer(4)              :: iElem, iNode, iGauss, iPoint
					real(4)                 :: aux
					! Initialize residual
					!$acc kernels
					Rconvec(:) = 0.0
					!$acc end kernels
					!$acc parallel loop gang
					do iElem = 1,nElem
						!$acc loop worker private(aux)
						do iGauss = 1,nGauss
							aux = 0.0
							!$acc loop seq
							do iNode = 1,nNodes
								aux = aux + (dNgp(iGauss,iNode) * phi(listConnec(iElem,iNode)))
							end do
							!$acc loop vector
							do iNode = 1,nNodes
								!$acc atomic update
								Rconvec(listConnec(iElem,iNode)) = Rconvec(listConnec(iElem,iNode)) + (wgp(iGauss) * Ngp(iGauss,iNode) * aux)
								!$acc end atomic
							end do
						end do
					end do
					!$acc end parallel loop
			end subroutine convec_acc
end module elemOps

program fem1d
	use mesh
	use FEM
	use elemOps
	use omp_lib
	use openacc
	implicit none
		integer(4), parameter   :: pOrder=10
		integer(4), parameter   :: nNodes=pOrder+1
		integer(4), parameter   :: nGauss=pOrder+1
		integer(4), parameter   :: nTimes=10
		integer(8), parameter   :: nElem=1e6*64
		integer(8), parameter   :: nPoints=(nElem+1) + (nElem*(pOrder-1))
		integer(8), allocatable :: listConnec(:,:)
		integer(4)              :: iNode, iGauss, iTime
		integer(8)              :: iElem, iPoint
		real(8)   , parameter   :: pi=3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986d0
		real(4)   , allocatable :: xyz(:), wgp(:), phi(:), Rconvec(:), Ngp(:,:), dNgp(:,:)
		real(4)                 :: hElem, Je, He
		real(4)                 :: t0, t1, time, avgTime
		! Display number of threads and number of nvidia devices
		write(*,'(a)') "!----------------------------------------!"
		write(*,'(a,i4)') "Number of threads: ", omp_get_max_threads()
		write(*,'(a,i4)') "Number of NVIDIA devices: ", acc_get_num_devices(acc_device_nvidia)
		! Generate the mesh
		allocate(listConnec(nElem,nNodes))
		allocate(xyz(nPoints))
		call generateMesh(nElem,pOrder,nNodes,nPoints,listConnec,xyz,hElem)
		write(*,'(a)') "!----------------------------------------!"
		write(*,'(a,i0)') "Number of elements: ", nElem
		write(*,'(a,i0)') "Number of points: ", nPoints
		! Generate FEM info
		allocate(wgp(nGauss))
		allocate(Ngp(nGauss,nNodes))
		allocate(dNgp(nGauss,nNodes))
		! Element Jacobian and inverse (special scenario for 1D eelements)
		Je = hElem/2.0
		He = 1.0/Je
		call evalElemInfo(pOrder,wgp,Ngp,dNgp)
		! Generate initial dummy data
		allocate(phi(nPoints))
		do iPoint = 1,nPoints
			phi(iPoint) = 100.0*sin(2*10*pi*xyz(iPoint))
		end do
		! Call and time the omp convective kernel n times
		allocate(Rconvec(nPoints))
		avgTime = 0.0
		do iTime = 1,nTimes
			call CPU_TIME(t0)
			call convec_omp(nElem,nNodes,nGauss,nPoints,listConnec,wgp,Ngp,dNgp,Je,He,phi,Rconvec)
			call CPU_TIME(t1)
			write(*,*) "Completed in ", t1-t0, "s"
			time = t1-t0
			avgTime = avgTime + time
		end do
		avgTime = avgTime / real(nTimes,8)
		write(*,*) "Avergae time: ", avgTime, "s"
		print *, sum(Rconvec)
		! Call and time the acc convective kernel n times
		allocate(Rconvec(nPoints))
		avgTime = 0.0
		!$acc enter data copyin(listConnec,wgp,Ngp,dNgp,Je,He,phi)
		do iTime = 1,nTimes
			call CPU_TIME(t0)
			call convec_acc(nElem,nNodes,nGauss,nPoints,listConnec,wgp,Ngp,dNgp,Je,He,phi,Rconvec)
			call CPU_TIME(t1)
			write(*,*) "Completed in ", t1-t0, "s"
			time = t1-t0
			avgTime = avgTime + time
		end do
		!$acc exit data delete(listConnec,wgp,Ngp,dNgp,Je,He,phi) copyout(Rconvec)
		avgTime = avgTime / real(nTimes,8)
		write(*,*) "Avergae time: ", avgTime, "s"
		print *, sum(Rconvec)
end program fem1d