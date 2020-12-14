module hessian_update
use definitions, only : dp
!use diagonalization, only : symmetric_upper

implicit none
private
public:: update_bfgs, update_sr1, update_psb, update_bofill

contains

! Hessian are stored in upper diagonal form (see diagonalization.F90)
!
!	Mat (symmetric,NxN) * Vec (Nx1) => Out (Nx1)
!	Out(1:N) = .0
!	do i = 1, N
!		t = .0
!		do j = 1, N
!			if( j < i )
!				k = j + i * ( i - 1 ) / 2
!			else
!				k = i + j * ( j - 1 ) / 2
!			t = t + Mat(k) * Vec(j)
!		end do
!		Out(i) = t
!	end do
!
!	Vec (Nx1) * Vec (1xN) => Out (symmetric,NxN)
!	Out(1:N*(N+1)/2) = .0
!	k = 0
!	do i = 1, N
!		do j = 1, i
!			k = k + 1
!			Out(k) = Vec(i) * Vec(j)
!		end do
!	end do
!

subroutine update_bfgs( n, dx, dg, hes )
!\begin{equation}
!\begin{split}
!s_{k}&=x_{k+1}-x_{k}\\
!y_{k}&=\nabla f \left( x_{k+1} \right)-\nabla f \left( x_{k} \right)\\
!H_{k+1}&=H_{k}-\frac{H_{k}s_{k}\left(H_{k}s_{k}\right)^{T}}{s_{k}^{T}H_{k}s_{k}}+
!         \frac{y_{k}y_{k}^{T}}{y_{k}^{T}s_{k}}
!\end{split}
!\end{equation}

	implicit none
	integer, intent(in) :: n
	real( kind=dp ), dimension(1:n), intent(in) :: dx, dg
	real( kind=dp ), dimension(1:n*(n+1)/2), intent(inout) :: hes

	real( kind=dp ), dimension(:), allocatable :: vec
	real( kind=dp ) :: t0, tys, tsvs
	integer :: i, j, k

	allocate( vec(1:n) )
	tys  = .0_dp
	tsvs = .0_dp
	do i = 1, n
		t0 = .0_dp
		do j = 1, n
			if( j < i ) then
				k = j + i * ( i - 1 ) / 2
			else
				k = i + j * ( j - 1 ) / 2
			end if
			t0 = t0 + dx(j) * hes(k)
		end do
		vec(i) = t0
		tsvs = tsvs + dx(i) * t0
		tys  = tys  + dg(i) * dx(i)
	end do
	k = 0
	do i = 1, n
		do j = 1, i
			k = k + 1
			hes(k) = hes(k) - vec(i) * vec(j) / tsvs + dg(i) * dg(j) / tys
		end do
	end do
	deallocate( vec )
	write( 6, "(a)" ) "Hessian updated by means of Broyden-Fletcher-Goldfarb-Shanno (BFGS)"
end subroutine


subroutine update_sr1( n, dx, dg, hes )
!\begin{equation}
!\begin{split}
!s_{k}&=x_{k+1}-x_{k}\\
!y_{k}&=\nabla f \left( x_{k+1} \right)-\nabla f \left( x_{k} \right)\\
!H_{k+1}&=H_{k}-\frac{\left(H_{k}s_{k}-y_{k}\right)\left(H_{k}s_{k}-y_{k}\right)^{T}}
!         {\left(H_{k}s_{k}-y_{k}\right)^{T}s_{k}}
!\end{split}
!\end{equation}

	implicit none
	integer, intent(in) :: n
	real( kind=dp ), dimension(1:n), intent(in) :: dx, dg
	real( kind=dp ), dimension(1:n*(n+1)/2), intent(inout) :: hes

	real( kind=dp ), dimension(:), allocatable :: vec
	real( kind=dp ) :: tvs, t0
	integer :: i, j, k

	allocate( vec(1:n) )
	tvs = .0_dp
	do i = 1, n
		t0 = .0_dp
		do j = 1, n
			if( j < i ) then
				k = j + i * ( i - 1 ) / 2
			else
				k = i + j * ( j - 1 ) / 2
			end if
			t0 = t0 + dx(j) * hes(k)
		end do
		vec(i) = t0 - dg(i)
		tvs = tvs + dx(i) * vec(i)
	end do
	k = 0
	do i = 1, n
		do j = 1, i
			k = k + 1
			hes(k) = hes(k) - vec(i) * vec(j) / tvs
		end do
	end do
	deallocate( vec )
	write( 6, "(a)" ) "Hessian updated by means of Symmetric Rank-one (SR1)"
end subroutine


subroutine update_psb( n, dx, dg, hes )
!\begin{equation}
!\begin{split}
!s_{k}&=x_{k+1}-x_{k}\\
!y_{k}&=\nabla f \left( x_{k+1} \right)-\nabla f \left( x_{k} \right)\\
!H_{k+1}&=H_{k}+\left(H_{k}s_{k}-y_{k}\right)^{T}s_{k}
!               \frac{s_{k}s_{k}^{T}}{\left(s_{k}^{T}s_{k}\right)^{2}}-
!               \frac{\left(H_{k}s_{k}-y_{k}\right)s_{k}^{T}+
!                s_{k}\left(H_{k}s_{k}-y_{k}\right)^{T}}
!               {s_{k}^{T}s_{k}}\\
!\end{split}
!\end{equation}

	implicit none
	integer, intent(in) :: n
	real( kind=dp ), dimension(1:n), intent(in) :: dx, dg
	real( kind=dp ), dimension(1:n*(n+1)/2), intent(inout) :: hes

	real( kind=dp ), dimension(:), allocatable :: vec
	real( kind=dp ) :: t0, tss, tvs
	integer :: i, j, k

	allocate( vec(1:n) )
	tss = .0_dp
	tvs = .0_dp
	do i = 1, n
		t0 = .0_dp
		do j = 1, n
			if( j < i ) then
				k = j + i * ( i - 1 ) / 2
			else
				k = i + j * ( j - 1 ) / 2
			end if
			t0 = t0 + dx(j) * hes(k)
		end do
		vec(i) = t0 - dg(i)
		tss = tss + dx(i)  * dx(i)
		tvs = tvs + vec(i) * dx(i)
	end do
	t0 = tss * tss
	k = 0
	do i = 1, n
		do j = 1, i
			k = k + 1
			hes(k) = hes(k) + dx(i) * dx(j) * tvs / t0 - ( vec(i) * dx(j) + vec(j) * dx(i) ) / tss
		end do
	end do
	deallocate( vec )
	write( 6, "(a)" ) "Hessian updated by means of Powell-symmetric-Broyden (PsB)"
end subroutine


subroutine update_bofill( n, dx, dg, hes )
!\begin{equation}
!\begin{split}
!s_{k}&=x_{k+1}-x_{k}\\
!y_{k}&=\nabla f \left( x_{k+1} \right)-\nabla f \left( x_{k} \right)\\
!\phi&=\frac{\left(\left(H_{k}s_{k}-y_{k}\right)^{T}s_{k}\right)^{2}}
!          {\left(H_{k}s_{k}-y_{k}\right)^{T}\left(H_{k}s_{k}-y_{k}\right)s_{k}^{T}s_{k}}\\
!H_{k}&=\phi H_{k}^{SR1}+(1-\phi)H_{k}^{PsB}
!\end{split}
!\end{equation}

	implicit none
	integer, intent(in) :: n
	real( kind=dp ), dimension(1:n), intent(in) :: dx, dg
	real( kind=dp ), dimension(1:n*(n+1)/2), intent(inout) :: hes

	real( kind=dp ), dimension(:), allocatable :: vec
	real( kind=dp ) :: phi, tvv, tvs, tss, t0
	integer :: i, j, k

	allocate( vec(1:n) )
	tvv = .0_dp
	tvs = .0_dp
	tss = .0_dp
	do i = 1, n
		t0 = .0_dp
		do j = 1, n
			if( j < i ) then
				k = j + i * ( i - 1 ) / 2
			else
				k = i + j * ( j - 1 ) / 2
			end if
			t0 = t0 + dx(j) * hes(k)
		end do
		vec(i) = t0 - dg(i)
		tvv = tvv + vec(i) * vec(i)
		tss = tss + dx(i)  * dx(i)
		tvs = tvs + vec(i) * dx(i)
	end do
	phi = tvs * tvs / ( tvv * tss )
	t0  = tss * tss
	k = 0
	do i = 1, n
		do j = 1, i
			k = k + 1
			hes(k) = hes(k) - phi * ( vec(i) * vec(j) / tvs ) &
							+ ( 1._dp - phi ) * ( dx(i) * dx(j) * tvs / t0 - &
										( vec(i) * dx(j) + vec(j) * dx(i) ) / tss )
		end do
	end do
	deallocate( vec )
	write( 6, "(a,f12.6)" ) "Hessian updated by means of Combined Bofill (Phi):", phi
end subroutine


subroutine update_dfp( n, dx, dg, hes )
!\begin{equation}
!\begin{split}
!s_{k}&=x_{k+1}-x_{k}\\
!y_{k}&=\nabla f \left( x_{k+1} \right)-\nabla f \left( x_{k} \right)\\
!H_{k+1}&=H_{k}+\left(H_{k}s_{k}-y_{k}\right)^{T}s_{k}
!               \frac{y_{k}y_{k}^{T}}{\left(y_{k}^{T}s_{k}\right)^{2}}-
!               \frac{\left(H_{k}s_{k}-y_{k}\right)y_{k}^{T}+
!                y_{k}\left(H_{k}s_{k}-y_{k}\right)^{T}}
!               {y_{k}^{T}s_{k}}\\
!\end{split}
!\end{equation}

	implicit none
	integer, intent(in) :: n
	real( kind=dp ), dimension(1:n), intent(in) :: dx, dg
	real( kind=dp ), dimension(1:n*(n+1)/2), intent(inout) :: hes

	real( kind=dp ), dimension(:), allocatable :: vec
	real( kind=dp ) :: t0, tgs, tvs, tvg
	integer :: i, j, k

	allocate( vec(1:n) )
	tgs = .0_dp
	tvs = .0_dp
	tvg = .0_dp
	do i = 1, n
		t0 = .0_dp
		do j = 1, n
			if( j < i ) then
				k = j + i * ( i - 1 ) / 2
			else
				k = i + j * ( j - 1 ) / 2
			end if
			t0 = t0 + dx(j) * hes(k)
		end do
		vec(i) = t0 - dg(i)
		tgs = tgs + dg(i)  * dx(i)
		tvs = tvs + vec(i) * dx(i)
		tvg = tvg + vec(i) * dg(i)
	end do
	t0 = tgs * tgs
	k = 0
	do i = 1, n
		do j = 1, i
			k = k + 1
			hes(k) = hes(k) + dg(i) * dg(j) * tvs / t0 - ( vec(i) * dg(j) + vec(j) * dg(i) ) / tgs
		end do
	end do
	deallocate( vec )
	write( 6, "(a)" ) "Hessian updated by means of Davidon-Fletcher-Powell (DFP)"
end subroutine

end module
