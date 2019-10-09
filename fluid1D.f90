program fluid1D

	implicit none
	real, parameter :: x0=0.0, xL=512.0, dx=0.25, tFinal=200, dt = 0.01
	real, parameter :: D=0.1, Eb=-1.0, xb=128.0
	integer :: i,N,iter
	double precision, allocatable, dimension(:) :: x, ne, np, E, E_CF, neNew, npNew, ENew
	double precision, allocatable, dimension(:) :: af, df, s, snew, afnew, dfnew
	real :: time
	
	
	
	N = int((xL - x0)/dx)
	allocate(x(N), ne(N), np(N), E(N), E_CF(N+1), neNew(N), npNew(N), ENew(N))
	allocate(af(N), afnew(N), df(N), dfnew(N), s(N), snew(N))
	
	do i=1,N
		x(i) = dx + (i-1)*dx
	end do
	
	!Initial conditions
	ne = 0.01*exp(-(x - xb)**2)
	np = 0.01*exp(-(x - xb)**2)
	!print *, ne
	E_CF(N+1) = Eb
	!print *,E_CF
	call calc_electricField(N, dx, ne, np, E, E_CF)
	!print *, E
	
	!Time integration---------------------------------------------------------------
	time = 0.0
	iter = 0
	print *, 'Starting integration..'
	do while (time < tFinal)
		
		
		call calc_advectionFlux(N, ne, E, dx, af)
		call calc_source(N, ne, E, s)	
		call calc_diffusionFlux(N, ne, E, dx, D, df)
		neNew = ne + dt*(af + df + s)
		npNew = np + dt*(s)
		call calc_electricField(N, dx, neNew, npNew, ENew, E_CF)
		
		call calc_advectionFlux(N, neNew, ENew, dx, afnew)
		call calc_source(N, neNew, ENew, snew)
		call calc_diffusionFlux(N, neNew, ENew, dx, D, dfnew)
		ne = ne + 0.5*dt*(afnew + dfnew + snew + af + df + s)
		np = np + 0.5*dt*(snew + s)
		call calc_electricField(N, dx, ne, np, E, E_CF)
		!print *, E
		
		if (sum(ne) .gt. 1e10) stop 'Solution Diverging!'
		print *, time + dt
		!Writing data
		if (mod(iter,50) .le. 1e-15) then
			call writeData(iter, x, ne, np, E, N)
		end if
		time = time + dt
		iter = iter + 1
	end do
	print *, "Integration done!"
	
	deallocate(x, ne, np, E, E_CF, neNew, npNew, ENew)
	deallocate(af, df, s, snew, afnew, dfnew)

!Subroutines used-------------------------------------------------------------------------------------------------------------------------

	contains

	subroutine calc_advectionFlux(N, ne, E, dx, advectionFlux)
		implicit none
		integer, intent(in) :: N
		double precision, dimension(N), intent(in) :: ne, E
		real, intent(in) :: dx
		double precision, dimension(N), intent(out) :: advectionFlux
		integer :: i
		
		do i=2,N
			advectionFlux(i) = (ne(i)*E(i) - ne(i-1)*E(i-1))/dx
		end do
		advectionFlux(1) = (ne(1)*E(1))/dx !ne and E at -1 node are taken to be zero
	end subroutine calc_advectionFlux
	
	
	subroutine calc_source(N, ne, E, S)
		implicit none
		integer, intent(in) :: N
		double precision, dimension(N), intent(in) :: ne, E
		double precision, dimension(N), intent(out) :: S
		integer :: i
		do i=2,N
			S(i) = ne(i)*abs(E(i))*exp(-1.0/abs(E(i)))
		end do
		S(1) = ne(1)*abs(E(1))*exp(-1.0/abs(E(1)))
	end subroutine calc_source
	
	subroutine calc_diffusionFlux(N, ne, E, dx, D, diffusionFlux)
		implicit none
		integer, intent(in) :: N
		double precision, dimension(N), intent(in) :: ne, E
		real, intent(in) :: dx, D
		double precision, dimension(N), intent(out) :: diffusionFlux
		integer :: i
		do i=2,N-1
			diffusionFlux(i) = (D/dx**2)*(ne(i+1) - 2.0*ne(i) + ne(i-1))
		end do
		diffusionFlux(1) = (D/dx**2)*(ne(2) - 2.0*ne(1))
		diffusionFlux(N) = (D/dx**2)*(ne(N-1) - ne(N))
		
	end subroutine calc_diffusionFlux
	
	subroutine calc_electricField(N, dx, ne, np, ECC, ECF)
		implicit none
		integer, intent(in) :: N
		double precision, dimension(N), intent(in) :: ne, np
		double precision, dimension(N), intent(out) :: ECC
		double precision, dimension(N+1), intent(out) :: ECF
		real, intent(in) :: dx
		integer :: i
		do i=1,N
			ECF(N-i+1) = ECF(N-i+2) + dx*(ne(N-i+1) - np(N-i+1))
		end do
		do i=1,N
			ECC(i) = 0.5*(ECF(i) + ECF(i+1)) 
		end do
	end subroutine calc_electricField
	
	subroutine writeData(timestep, x, ne, np, ECC, N)
		implicit none
		integer, intent(in) :: N, timestep
		double precision, dimension(N), intent(in) :: ne, np, x, ECC
		integer :: i
		character(21) :: filename
 	    	write(filename, '(a,i0.6,a)') 'fortran_op/', timestep, '.txt'
		open(1, file= filename, status='new')
		write(1, *) 'x ','sigma ', 'rho ', 'E '
		do i=1,N
			write(1, *) x(i), e(i), np(i), E(i)
		end do
		close(1)		
	end subroutine writeData
	

end program fluid1D





