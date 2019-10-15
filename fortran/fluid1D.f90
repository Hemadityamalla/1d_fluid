!Description:Code to propagate a negative streamer with a constant background electric field.
!            The equations and parameters used are all in their non-dimensional form
!	     We use a finite volume discretization method as follows:
!	     Advection term: Upwind discretization
!	     Diffusion term: 2nd order central dscretization
!            The time integration is performed using the trapezoidal rule (two stage)
!Author: Hemaditya


program fluid1D
	implicit none
	!Initializing the parameters used for the simulation
	real, parameter :: x0=0.0, xL=512.0, dx=0.05, tFinal=200.0, dt = 0.01
	real, parameter :: D=0.1, Eb=-1.0, xb=31.0
	integer :: i,N,iter !Initializing parameters for arrays and loops

	!Initializing the arrays for the various fields
	double precision, allocatable, dimension(:) :: x, ne, np, E, E_CF, neNew, npNew, ENew
	double precision, allocatable, dimension(:) :: af, df, s, snew, afnew, dfnew
	real :: time, tempVar1, tempVar2, maxCFL !Variables used to compute the cfl number and to keep track of time
	real, dimension(20000) :: tStep, fPos !Array to collect the position of the front
	
	
	
	N = int((xL - x0)/dx) !Evaluating the number of cells

	!Allocating the arrays
	allocate(x(N), ne(N), np(N), E(N), E_CF(N+1), neNew(N), npNew(N), ENew(N))
	allocate(af(N), afnew(N), df(N), dfnew(N), s(N), snew(N))
	
	!Computing the cell center coordinates
	do i=1,N
		x(i) = dx + (i-1)*dx
	end do
	
	!Initial conditions
	ne = 0.01*exp(-(x - xb)**2)
	np = 0.01*exp(-(x - xb)**2)
	E_CF(N+1) = Eb
	call calc_electricField(N, dx, ne, np, E, E_CF) !Calculating the electric field at the cell centers
	
	!Time integration---------------------------------------------------------------
	time = 0.0
	iter = 0
        tempVar1 = 0.0
        tempVar2 = 0.0
        maxCFL = 0.0
	print *, 'Starting integration..'
	do while (time < tFinal)
		
		!Calculating the fluxes and source terms for t_n/2
		call calc_advectionFlux(N, ne, E, dx, af)
		call calc_source(N, ne, E, s)	
		call calc_diffusionFlux(N, ne, E, dx, D, df)

		!Upate solutions for t_n/2
		neNew = ne + dt*(af + df + s)
		npNew = np + dt*(s)
		call calc_electricField(N, dx, neNew, npNew, ENew, E_CF)
		
		!Calculating the fluxes and source terms for t_n+1
		call calc_advectionFlux(N, neNew, ENew, dx, afnew)
		call calc_source(N, neNew, ENew, snew)
		call calc_diffusionFlux(N, neNew, ENew, dx, D, dfnew)

		!Update solution for t_n+1
		ne = ne + 0.5*dt*(afnew + dfnew + snew + af + df + s)
		np = np + 0.5*dt*(snew + s)
		call calc_electricField(N, dx, ne, np, E, E_CF)
		
		!Check and abort if solution is blowing up
		if (sum(ne) .gt. 1e10) stop 'Solution Diverging!'
                do i=1,N
                        if (isnan(ne(i))) stop 'ne is NaN'
                end do
		!Printing the time step and the maximum CFL for each timestep
		print *, time + dt, maxCFL

		!Writing data every 50 iterations
		if (mod(iter,50) .le. 1e-15) then
			call writeData(iter, x, ne, np, E, N)
		end if
               
		!Update time
		time = time + dt
		iter = iter + 1

                !Computing the front velocity using the maximum value of the electron density
                if (mod(iter,100) .le. 1e-15) then
                tempVar1 = tempVar2
                tempVar2 = x(maxloc(ne, dim=1))
                end if
                maxCFL = max(maxCFL, maxval(abs(E)*(dx/dt))) !Computing the max CFL

		!Collecting the front positiions
		if (iter .le. 20000) then
		tStep(iter) = time
		fPos(iter) =  x(maxloc(ne, dim=1))
		end if	
        end do
	print *, "Integration done!"
        print *, (tempVar2 - tempVar1)/(100*dt), iter !This computes the final front velocity

	!Writing the times and front positions to a file
	open(8, file= 'frontPos.dat', status='new')
	write(8, *) 'time ', 'position '
	do i=1,20000
		write(8, *) tStep(i), fPos(i)
	end do
	close(8)
	
	!Deallocating memory for the arrays
	deallocate(x, ne, np, E, E_CF, neNew, npNew, ENew)
	deallocate(af, df, s, snew, afnew, dfnew)

!Subroutines used-------------------------------------------------------------------------------------------------------------------------

	contains
	!To calculate the advection flux using the Upwind scheme
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
	
	!To calculate the simple source term
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
	
	!To  calculate the diffusion flux using the 2nd order central discretization
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
	
	!To calculate the electric fields at the cell faces and the cell centers
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
	
	!To write the data for a given timestep
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
			write(1, *) x(i), ne(i), np(i), E(i)
		end do
		close(1)		
	end subroutine writeData
end program fluid1D





