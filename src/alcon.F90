program alcon
implicit none

! declaration of variables
#include "alcon_dec.h90"

! start timing
call date_and_time(cdate(1), cdate(2))

! initialization
call SlepcInitialize(PETSC_NULL_CHARACTER, ierr)
CHKERRQ(ierr)
call MPI_Comm_rank(MPI_COMM_WORLD, mype, ierr)
CHKERRQ(ierr)

! load data into matrixes
nrad = 0
nfft = 0
if (mype == 0) then
  open (facd, file = 'alcon.dat', status = 'old')
  read (facd, *) nrad, nfft, nprofiledata, nfftdata
  write (*, *) "nrad = ", nrad
  write (*, *) "nfft = ", nfft
  write (*, *) "noffdiag = ", noffdiag
  if (nfft < noffdiag) then
    call PetscPrintf(MPI_COMM_WORLD, "Warning: nfft < noffdiag, not enough data to calculate specified coupling.\n", ierr)
    CHKERRQ(ierr)
  endif
  allocate (acdprofile(nprofiledata, nrad), acdfft(0 : nfft, nfftdata, nrad))
  read (facd, *) acdprofile
  read (facd, *) acdfft
  close (facd)
  write (*, *) "Monitor profile data on 1st radial grid:"
  write (*, *) "sqrt(tor/torw) = ", acdprofile(1, 1)
  write (*, *) "q = ", acdprofile(2, 1)
  write (*, *) "gq + I = ", acdprofile(3, 1)
  write (*, *) "beta = ", acdprofile(4, 1)
  write (*, *) "rho_M = ", acdprofile(5, 1)
  write (*, *) "Monitor FFT quantities:"
  write (*, *) "1st radial grid:"
  write (*, *) "(hcoef)_0 = ", acdfft(0, 1, 1)
  write (*, *) "(hcoef)_1 = ", acdfft(1, 1, 1)
  write (*, *) "(jcoef)_0 = ", acdfft(0, 2, 1)
  write (*, *) "(jcoef)_1 = ", acdfft(1, 2, 1)
  write (*, *) "(kcoef)_0 = ", acdfft(0, 3, 1)
  write (*, *) "(kcoef)_1 = ", acdfft(1, 3, 1)
  write (*, *) "(lcoef)_0 = ", acdfft(0, 4, 1)
  write (*, *) "(lcoef)_1 = ", acdfft(1, 4, 1)
  write (*, *) "(ncoef)_0 = ", acdfft(0, 5, 1)
  write (*, *) "(ncoef)_1 = ", acdfft(1, 5, 1)
  write (*, *) "Radial grid #:", iraddiag
  write (*, *) "(hcoef)_0 = ", acdfft(0, 1, iraddiag)
  write (*, *) "(hcoef)_1 = ", acdfft(1, 1, iraddiag)
  write (*, *) "(jcoef)_0 = ", acdfft(0, 2, iraddiag)
  write (*, *) "(jcoef)_1 = ", acdfft(1, 2, iraddiag)
  write (*, *) "(kcoef)_0 = ", acdfft(0, 3, iraddiag)
  write (*, *) "(kcoef)_1 = ", acdfft(1, 3, iraddiag)
  write (*, *) "(lcoef)_0 = ", acdfft(0, 4, iraddiag)
  write (*, *) "(lcoef)_1 = ", acdfft(1, 4, iraddiag)
  write (*, *) "(ncoef)_0 = ", acdfft(0, 5, iraddiag)
  write (*, *) "(ncoef)_1 = ", acdfft(1, 5, iraddiag)
endif

intbuf(1) = nrad
intbuf(2) = nfft
call MPI_Bcast(intbuf, 2, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
CHKERRQ(ierr)
nrad = intbuf(1)
nfft = intbuf(2)

! initialization of PETSc and SLEPc variables
#include "alcon_init.h90"

fgpdat = 0

! main loop
do irad = 1, nrad
!do irad = iraddiag - 1, iraddiag
  write (msg, *) "Solving radial grid: irad = ", irad, "\n"
  call PetscPrintf(MPI_COMM_WORLD, msg, ierr)
  CHKERRQ(ierr)
  if (mype == 0) then
    ! set values for matA and matB
    do i = 0, m2 - m1
      if (finitebeta == 2) then
        ! beta * G^dagger * M *G part
        num_ind = 0
        indices(num_ind) = i
        values(num_ind) = acdprofile(4, irad) * (dble(n) * acdprofile(2, irad) - dble(m1 + i))**2 / acdprofile(3, irad)
        num_ind = num_ind + 1
        call MatSetValues(matA, 1, i, num_ind, indices, values, INSERT_VALUES, ierr)

        num_ind = 0
        do j = max(0, i - min(nfft, noffdiag)), min(m2 - m1, i + min(nfft, noffdiag))
          ! L part
          indices(num_ind) = j
          if (j < i) then
            values(num_ind) = dconjg(acdfft(i - j, 4, irad))
          else
            values(num_ind) = acdfft(j - i, 4, irad)
          endif
          num_ind = num_ind + 1
        enddo
        call MatSetValues(matB, 1, i, num_ind, indices, values, INSERT_VALUES, ierr)
      else
        num_ind = 0
        do j = max(0, i - min(nfft, noffdiag)), min(m2 - m1, i + min(nfft, noffdiag))
          ! G^dagger * H * G part
          indices(num_ind) = j
          if (j < i) then
            values(num_ind) = dconjg(acdfft(i - j, 1, irad))
          else
            values(num_ind) = acdfft(j - i, 1, irad)
          endif
          values(num_ind) = values(num_ind) &
            * (dble(n) * acdprofile(2, irad) - dble(m1 + i)) &
            * (dble(n) * acdprofile(2, irad) - dble(m1 + j))
          ! slow sound approximation (low beta), A = G^dagger * H * G part + N
          if (finitebeta == 1) then
            if (j < i) then
              values(num_ind) = values(num_ind) + dconjg(acdfft(i - j, 5, irad))
            else
              values(num_ind) = values(num_ind) + acdfft(j - i, 5, irad)
            endif
          endif
          num_ind = num_ind + 1
          ! -beta * K part
          if (finitebeta > 2) then
            indices(num_ind) = m2 - m1 + 1 + j
            if (j < i) then
              values(num_ind) = -acdprofile(4, irad) * dconjg(acdfft(i - j, 3, irad))
            else
              values(num_ind) = -acdprofile(4, irad) * acdfft(j - i, 3, irad)
            endif
            num_ind = num_ind + 1
          endif
        enddo
        call MatSetValues(matA, 1, i, num_ind, indices, values, INSERT_VALUES, ierr)
        CHKERRQ(ierr)
        ! J part
        num_ind = 0
        do j = max(0, i - min(nfft, noffdiag)), min(m2 - m1, i + min(nfft, noffdiag))
          indices(num_ind) = j
          if (j < i) then
            values(num_ind) = dconjg(acdfft(i - j, 2, irad))
          else
            values(num_ind) = acdfft(j - i, 2, irad)
          endif
          num_ind = num_ind + 1
        enddo
        call MatSetValues(matB, 1, i, num_ind, indices, values, INSERT_VALUES, ierr)
        CHKERRQ(ierr)
      endif ! else of (finitebeta == 2)
      if (finitebeta > 2) then
        ! beta * G^dagger * M *G part
        num_ind = 0
        indices(num_ind) = m2 - m1 + 1 + i
        values(num_ind) = acdprofile(4, irad) * (dble(n) * acdprofile(2, irad) - dble(m1 + i))**2 / acdprofile(3, irad)
        num_ind = num_ind + 1
        call MatSetValues(matA, 1, m2 - m1 + 1 + i, num_ind, indices, values, INSERT_VALUES, ierr)

        num_ind = 0
        do j = max(0, i - min(nfft, noffdiag)), min(m2 - m1, i + min(nfft, noffdiag))
          ! K part
          indices(num_ind) = j
          if (j < i) then
            values(num_ind) = dconjg(acdfft(i - j, 3, irad))
          else
            values(num_ind) = acdfft(j - i, 3, irad)
          endif
          num_ind = num_ind + 1
          ! L part
          indices(num_ind) = m2 - m1 + 1 + j
          if (j < i) then
            values(num_ind) = dconjg(acdfft(i - j, 4, irad))
          else
            values(num_ind) = acdfft(j - i, 4, irad)
          endif
          num_ind = num_ind + 1
        enddo
        call MatSetValues(matB, 1, m2 - m1 + 1 + i, num_ind, indices, values, INSERT_VALUES, ierr)
      endif ! (finitebeta > 2)
    enddo ! i = 0, m2 - m1
  endif ! (mype == 0)
  call MatAssemblyBegin(matA, MAT_FINAL_ASSEMBLY, ierr)
  CHKERRQ(ierr)
  call MatAssemblyBegin(matB, MAT_FINAL_ASSEMBLY, ierr)
  CHKERRQ(ierr)
  call MatAssemblyEnd(matA, MAT_FINAL_ASSEMBLY, ierr)
  CHKERRQ(ierr)
  call MatAssemblyEnd(matB, MAT_FINAL_ASSEMBLY, ierr)
  CHKERRQ(ierr)

  ! print out matA, matB on irad == 1 and irad == iraddiag
  if (irad == 1 .or. irad == iraddiag) then
    call PetscPrintf(MPI_COMM_WORLD, "matA:\n", ierr)
    CHKERRQ(ierr)
    call MatView(matA, PETSC_VIEWER_STDOUT_WORLD, ierr)
    CHKERRQ(ierr)
    call PetscPrintf(MPI_COMM_WORLD, "matB:\n", ierr)
    CHKERRQ(ierr)
    call MatView(matB, PETSC_VIEWER_STDOUT_WORLD, ierr)
    CHKERRQ(ierr)
  endif

  ! initialize acsolver
  call EPSCreate(MPI_COMM_WORLD, acsolver, ierr)
  CHKERRQ(ierr)
  !call EPSSetProblemType(acsolver, EPS_PGNHEP, ierr)
  !CHKERRQ(ierr)
  call EPSSetWhichEigenpairs(acsolver, EPS_SMALLEST_MAGNITUDE, ierr)
  CHKERRQ(ierr)
  if (finitebeta > 2) then
    call EPSSetDimensions(acsolver, (m2 - m1 + 1) * 2, PETSC_DECIDE, PETSC_DECIDE, ierr)
    CHKERRQ(ierr)
  elseif (finitebeta == 2) then
    call EPSSetDimensions(acsolver, m2 - m1 + 1, PETSC_DECIDE, PETSC_DECIDE, ierr)
    CHKERRQ(ierr)
  else
    call EPSSetDimensions(acsolver, nev, PETSC_DECIDE, PETSC_DECIDE, ierr)
    CHKERRQ(ierr)
  endif
  call EPSSetOperators(acsolver, matA, matB, ierr)
  CHKERRQ(ierr)
  call EPSSetType(acsolver, EPSLAPACK, ierr)
  CHKERRQ(ierr)
  ! solve the eigenvalue problem
  call EPSSolve(acsolver, ierr)
  CHKERRQ(ierr)
  ! write out information about the solver on irad == 1
  if (irad == 1) then
    call EPSGetIterationNumber(acsolver, its, ierr)
    CHKERRQ(ierr)
    write (msg, *) "# of iterations: ", its, "\n"
    call PetscPrintf(MPI_COMM_WORLD, msg, ierr)
    CHKERRQ(ierr)
    call EPSGetType(acsolver, tname, ierr)
    CHKERRQ(ierr)
    write (msg, *) "Solution method: ", tname, "\n"
    call PetscPrintf(MPI_COMM_WORLD, msg, ierr)
    CHKERRQ(ierr)
    call EPSGetDimensions(acsolver, nev_out, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, ierr)
    CHKERRQ(ierr)
    write (msg, *) "# of requested eigenvalues: ", nev_out, "\n"
    call PetscPrintf(MPI_COMM_WORLD, msg, ierr)
    CHKERRQ(ierr)
    call EPSGetTolerances(acsolver, tol, maxit, ierr)
    CHKERRQ(ierr)
    write (msg, *) "Stopping condition: tol = ", tol, ", maxit = ", maxit, "\n"
    call PetscPrintf(MPI_COMM_WORLD, msg, ierr)
    CHKERRQ(ierr)
  endif 

  call EPSGetConverged(acsolver, nconv, ierr)
  CHKERRQ(ierr)
  write (msg, *) "# of converged solutions: nconv = ", nconv, "\n"
  call PetscPrintf(MPI_COMM_WORLD, msg, ierr)
  CHKERRQ(ierr)
  ! analyze converged solutions
  do i = 0, nconv - 1
    call EPSGetEigenpair(acsolver, i, evalre, evalim, evecre, evecim, ierr)
    CHKERRQ(ierr)
!    if (i == 0 .and. irad == 1) then
!      call VecView(evecre, PETSC_VIEWER_STDOUT_WORLD, ierr)
!      CHKERRQ(ierr)
!      call VecView(evecim, PETSC_VIEWER_STDOUT_WORLD, ierr)
!      CHKERRQ(ierr)
!    endif
    ! if the absolute value of omega's imaginary part is small enough, then output
    if (abs(dimag(sqrt(evalre))) < dble(sqrt(evalre)) * imagrealratiocutoff) then
      call VecAbs(evecre, ierr)
      CHKERRQ(ierr)
      call VecMax(evecre, j, evecreabsmax, ierr)
      CHKERRQ(ierr)
      ! write useful data to files
!      if (j <= m2 - m1 .and. mype == 0) then
      if (mype == 0) then
!        if (j > m2 - m1) j = j - (m2 - m1 + 1)
        omegarescaled = dble(sqrt(evalre / acdprofile(5, irad))) * omegascale
        if (omegarescaled < omegacutoff) then
          if (fgpdat(j) == 0) then
            fgpdat(j) = 1000 + j
            if (finitebeta > 2) then
              if (j < m2 - m1 + 1) then
                if (m1 + j < 0) then
                  write (msg, '("a_n", i3.3, "m", i3.2, ".gp.dat")') n, m1 + j
                else
                  write (msg, '("a_n", i3.3, "m", i3.3, ".gp.dat")') n, m1 + j
                endif
              else
                if (m1 + j - (m2 - m1 + 1) < 0) then
                  write (msg, '("s_n", i3.3, "m", i3.2, ".gp.dat")') n, m1 + j - (m2 - m1 + 1)
                else
                  write (msg, '("s_n", i3.3, "m", i3.3, ".gp.dat")') n, m1 + j - (m2 - m1 + 1)
                endif
              endif
            else
              if (m1 + j < 0) then
                write (msg, '("n", i3.3, "m", i3.2, ".gp.dat")') n, m1 + j
              else
                write (msg, '("n", i3.3, "m", i3.3, ".gp.dat")') n, m1 + j
              endif
            endif
            open (fgpdat(j), file = msg, status = 'replace')
          endif
          write (fgpdat(j), *) acdprofile(1, irad), omegarescaled
          if (i == 0) then
            write (*, *) "Lowest omega: ", omegarescaled
            write (*, *) "Corresponding imaginary part: ", dimag(sqrt(evalre / acdprofile(5, irad))) * omegascale
            if (j < m2 - m1 + 1) then
              write (*, *) "Corresponding m: ", m1 + j
            else
              write (*, *) "Corresponding m: ", m1 + j - (m2 - m1 + 1)
            endif
          endif
        endif ! (omegarescaled < omegacutoff)
      endif ! (mype == 0)
    endif ! (dble(evalre) > 0)
  enddo ! do i = 0, nconv - 1
  ! destroy acsolver
  call EPSDestroy(acsolver, ierr)
  CHKERRQ(ierr)
enddo

! finalization of PETSc and SLEPc variables
#include "alcon_final.h90"

call SlepcFinalize(ierr)
CHKERRQ(ierr)
end program alcon

