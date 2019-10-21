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
	real, parameter :: x0=0.0, xL=512.0, dx=0.1, tFinal=200.0, dt = 0.01
	real, parameter :: D=0.1, Eb=-1.0, xb=31.0
	integer :: i,N,iter, nGhost !Initializing parameters for arrays and loops

	!Initializing the arrays for the various fields
	double precision, allocatable, dimension(:) :: x, ne, np, E, E_CF, neNew, npNew, ENew
	double precision, allocatable, dimension(:) :: af, df, s, snew, afnew, dfnew
	real :: time, tempVar1, tempVar2, maxCFL !Variables used to compute the cfl number and to keep track of time
	real, dimension(20000) :: tStep, fPos !Array to collect the position of the front
	
	
	
	N = int((xL - x0)/dx) !Evaluating the number of cells
	nGhost = 4 !Number of ghost cells 2 for each boundary in this case

	!Allocating the arrays and initializing everything with zero
	allocate(x(N), ne(N+nGhost), np(N+nGhost), E(N+nGhost), E_CF(N+1), neNew(N+nGhost), npNew(N+nGhost), ENew(N+nGhost))
	allocate(af(N), afnew(N), df(N), dfnew(N), s(N), snew(N))
	ne = 0.0
	np = 0.0
	E_CF = 0.0
	E = 0.0
	neNew = 0.0
	npNew = 0.0
	ENew = 0.0
	
	
	!Computing the cell center coordinates
	do i=1,N
		x(i) = dx + (i-1)*dx
	end do
	
	!Initial conditions
	ne(3:N+2) = 0.01*exp(-(x - xb)**2)
	np(3:N+2) = 0.01*exp(-(x - xb)**2)
	E_CF(N+1) = Eb
	call calc_electricField(N,nGhost, dx, ne, np, E, E_CF) !Calculating the electric field at the cell centers
	!Ghost cells
	ne(1) = -ne(4)
	ne(2) = -ne(3)
	np(1) = -np(4)
	np(2) = -np(3)
	ne(N+3) = ne(N+2)
	ne(N+4) = ne(N+1)
	np(N+3) = ne(N+2)
	np(N+4) = ne(N+1)
	
	!Time integration---------------------------------------------------------------
	time = 0.0
	iter = 0
        tempVar1 = 0.0
        tempVar2 = 0.0
        maxCFL = 0.0
	print *, 'Starting integration..'
	do while (time < tFinal)
		
		!Calculating the fluxes and source terms for t_n/2
		call calc_advectionFlux(N,nGhost, ne, E_CF, dx, af)
		call calc_source(N,nGhost, ne, E, s)	
		call calc_diffusionFlux(N,nGhost, ne, E, dx, D, df)

		!Upate solutions for t_n/2
		neNew(3:N+2) = ne(3:N+2) + dt*(af + df + s)
		npNew(3:N+2) = np(3:N+2) + dt*(s)
		call calc_electricField(N,nGhost, dx, neNew, npNew, ENew, E_CF)
		!Ghost cells
		neNew(1) = -neNew(4)
		neNew(2) = -neNew(3)
		npNew(1) = -npNew(4)
		npNew(2) = -npNew(3)
		neNew(N+3) = neNew(N+2)
		neNew(N+4) = nenew(N+1)
		npNew(N+3) = npNew(N+2)
		npNew(N+4) = npNew(N+1)
		
		!Calculating the fluxes and source terms for t_n+1
		call calc_advectionFlux(N,nGhost, neNew, E_CF, dx, afnew)
		call calc_source(N,nGhost, neNew, ENew, snew)
		call calc_diffusionFlux(N,nGhost, neNew, ENew, dx, D, dfnew)

		!Update solution for t_n+1
		ne(3:N+2) = ne(3:N+2) + 0.5*dt*(afnew + dfnew + snew + af + df + s)
		np(3:N+2) = np(3:N+2) + 0.5*dt*(snew + s)
		call calc_electricField(N,nGhost, dx, ne, np, E, E_CF)
		!Ghost cells
		ne(1) = -ne(4)
		ne(2) = -ne(3)
		np(1) = -np(4)
		np(2) = -np(3)
		ne(N+3) = ne(N+2)
		ne(N+4) = ne(N+1)
		np(N+3) = ne(N+2)
		np(N+4) = ne(N+1)
		
		!Check and abort if solution is blowing up
		if (sum(ne) .gt. 1e10) stop 'Solution Diverging!'
                do i=1,N
                        if (isnan(ne(i))) stop 'ne is NaN'
                end do
		!Printing the time step and the maximum CFL for each timestep
		print *, time + dt, maxCFL

		!Writing data every 50 iterations
		if (mod(iter,50) .le. 1e-15) then
			call writeData(iter, x, ne(3:N+2), np(3:N+2), E(3:N+2), N)
		end if
               
		!Update time
		time = time + dt
		iter = iter + 1

                !Computing the front velocity using the maximum value of the electron density
                if (mod(iter,100) .le. 1e-15) then
                tempVar1 = tempVar2
                tempVar2 = x(maxloc(ne(3:N+2), dim=1))
                end if
                maxCFL = max(maxCFL, maxval(abs(E(3:N+2))*(dx/dt))) !Computing the max CFL

		!Collecting the front positiions
		if (iter .le. 20000) then
		tStep(iter) = time
		fPos(iter) =  x(maxloc(ne(3:N+2), dim=1))
		end if
		!return
        end do
	print *, "Integration done!"
        print *, (tempVar2 - tempVar1)/(100*dt), iter !This computes the final front velocity

	!Writing the times and front positions to a file
! 	open(8, file= 'frontPos.dat', status='new')
! 	write(8, *) 'time ', 'position '
! 	do i=1,20000
! 		write(8, *) tStep(i), fPos(i)
! 	end do
! 	close(8)
	
	!Deallocating memory for the arrays
	deallocate(x, ne, np, E, E_CF, neNew, npNew, ENew)
	deallocate(af, df, s, snew, afnew, dfnew)


!Subroutines and functions used-------------------------------------------------------------------------------------------------------------------------

	contains
	
	!To calculate the limiter
	double precision function phi(r)
		
		implicit none
		double precision, intent(in) :: r
		
		!phi = 0.0 !Upwind 
		!phi = 1.0 !Lax-Wendroff
		!phi = r !Beam-Warming
		!phi = 0.5*(1.0 + r) !Fromm
		
		!phi = (r + abs(r))/(1 + abs(r)) !van Leer
		!phi = maxval((/0.0, min(1.0, 2*r), min(2.0, r)/)) !supberbee
		!phi = int(r .gt. 0.0)*min(r, 1.0) + int(r .le. 0.0)*0.0 !minmod
		phi = max(0.0, min(1.0/3.0 + r/6.0, r)) !koren
		
	end function phi
	
	!To calculate the advection flux using the Upwind scheme
	subroutine calc_advectionFlux(N,Ng, ne, ECF, dx, advectionFlux)
		implicit none
		integer, intent(in) :: N, Ng
		double precision, dimension(N+Ng), intent(in) :: ne
		double precision, dimension(N+1), intent(in) :: ECF
		real, intent(in) :: dx
		double precision, dimension(N), intent(out) :: advectionFlux
		integer :: i
		double precision :: thetar, thetal, rf, lf
		real, parameter :: eps = 1e-10
		double precision, dimension(N+1) :: F
		
		!do i=2,N
		!	advectionFlux(i) = (ne(i)*E(i) - ne(i-1)*E(i-1))/dx
		!end do
		!advectionFlux(1) = (ne(1)*E(1))/dx !ne and E at -1 node are taken to be zero
		do i=1,N+1
			thetar = (ne((i - 1) + 3) - ne((i-1)+2)+eps)/(ne((i-1) + 2) - ne(i)+eps)
			rf = min(ECF(i), 0.0)*(ne((i-1) + 2) + phi(thetar)*(ne((i-1) + 2) - ne(i)))
			thetal = (ne((i-1) + 2) - ne((i-1) + 3)+eps)/(ne((i-1) + 3) - ne((i-1) + 4)+eps)
			lf = max(ECF(i), 0.0)*(ne((i-1) + 3) + phi(thetal)*(ne((i-1) + 3) - ne((i-1) + 4)))
			!if (i .eq. N+1) print *, phi(thetar)
			F(i) = rf + lf
		end do
		advectionFlux = (F(2:N+1) - F(1:N))/dx
		
	end subroutine calc_advectionFlux
	
	!To calculate the simple source term
	subroutine calc_source(N,Ng, ne, E, S)
		implicit none
		integer, intent(in) :: N, Ng
		double precision, dimension(N+Ng), intent(in) :: ne, E
		double precision, dimension(N), intent(out) :: S
		integer :: i
		!do i=2,N
		!	S(i) = ne(i)*abs(E(i))*exp(-1.0/abs(E(i)))
		!end do
		!S(1) = ne(1)*abs(E(1))*exp(-1.0/abs(E(1)))
		
		S = ne(3:N+2)*abs(E(3:N+2))*exp(-1.0/abs(E(3:N+2)))
	end subroutine calc_source
	
	!To  calculate the diffusion flux using the 2nd order central discretization
	subroutine calc_diffusionFlux(N,Ng, ne, E, dx, D, diffusionFlux)
		implicit none
		integer, intent(in) :: N, Ng
		double precision, dimension(N+Ng), intent(in) :: ne, E
		real, intent(in) :: dx, D
		double precision, dimension(N), intent(out) :: diffusionFlux
		integer :: i
		!do i=2,N-1
		!	diffusionFlux(i) = (D/dx**2)*(ne(i+1) - 2.0*ne(i) + ne(i-1))
		!end do
		!diffusionFlux(1) = (D/dx**2)*(ne(2) - 2.0*ne(1))
		!diffusionFlux(N) = (D/dx**2)*(ne(N-1) - ne(N))
		
		diffusionFlux = (D/dx**2)*(ne(4:N+3) - 2.0*ne(3:N+2) + ne(2:N+1)) !Second order central discretization
		
	end subroutine calc_diffusionFlux
	
	!To calculate the electric fields at the cell faces and the cell centers
	subroutine calc_electricField(N,Ng, dx, ne, np, ECC, ECF)
		implicit none
		integer, intent(in) :: N, Ng
		double precision, dimension(N+Ng), intent(in) :: ne, np
		double precision, dimension(N+Ng), intent(out) :: ECC
		double precision, dimension(N+1), intent(out) :: ECF
		real, intent(in) :: dx
		integer :: i
		do i=1,N
			ECF(N-i+1) = ECF(N-i+2) + dx*(ne(N-i+3) - np(N-i+3))
		end do
		!ECF(1:N) = ECF(2:N+1) + dx*(ne(2:N+1) - np(2:N+1))-- this is not giving correct values
		!do i=1,N
		!	ECC(i) = 0.5*(ECF(i) + ECF(i+1)) 
		!end do
		ECC(3:N+2) = 0.5*(ECF(1:N) + ECF(2:N+1))

	end subroutine calc_electricField
	
	!To write the data for a given timestep
	subroutine writeData(timestep, x, ne, np, ECC, N)
		implicit none
		integer, intent(in) :: N, timestep
		double precision, dimension(N), intent(in) :: x, ne, np, ECC
		integer :: i
		character(21) :: filename
 	    	write(filename, '(a,i0.6,a)') 'fortran_op/', timestep, '.txt'
		open(1, file= filename, status='new')
		write(1, *) 'x ','sigma ', 'rho ', 'E '
		do i=1,N
			write(1, *) x(i), ne(i), np(i), ECC(i)
		end do
		close(1)		
	end subroutine writeData
end program fluid1D





