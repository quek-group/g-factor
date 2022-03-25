!
! Copyright (C) 2001-2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
SUBROUTINE c_bands( iter )
  !----------------------------------------------------------------------------
  !
  ! ... Driver routine for Hamiltonian diagonalization routines
  ! ... It reads the Hamiltonian and an initial guess of the wavefunctions
  ! ... from a file and computes initialization quantities for the
  ! ... diagonalization routines.
  !
  USE kinds,                ONLY : DP
  USE io_global,            ONLY : stdout
  USE io_files,             ONLY : iunhub, iunwfc, nwordwfc, nwordwfcU
  USE buffers,              ONLY : get_buffer, save_buffer, close_buffer
  USE klist,                ONLY : nkstot, nks, xk, ngk, igk_k
  USE uspp,                 ONLY : vkb, nkb
  USE gvect,                ONLY : g
  USE wvfct,                ONLY : et, nbnd, npwx, current_k
  USE control_flags,        ONLY : ethr, isolve, restart
  USE ldaU,                 ONLY : lda_plus_u, U_projection, wfcU
  USE lsda_mod,             ONLY : current_spin, lsda, isk
  USE wavefunctions_module, ONLY : evc
  USE bp,                   ONLY : lelfield
  USE mp_pools,             ONLY : npool, kunit, inter_pool_comm
  USE mp,                   ONLY : mp_sum
  USE check_stop,           ONLY : check_stop_now
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT (in) :: iter
  !
  ! ... local variables
  !
  REAL(DP) :: avg_iter
  ! average number of H*psi products
  INTEGER :: ik_, ik, nkdum, ios
  ! ik : counter on k points
  ! ik_: k-point already done in a previous run
  LOGICAL :: exst
  !
  CALL start_clock( 'c_bands' )
  !
  ik_ = 0
  avg_iter = 0.D0
  IF ( restart ) CALL restart_in_cbands(ik_, ethr, avg_iter, et )
  !
  ! ... If restarting, calculated wavefunctions have to be read from file
  ! ... (not needed for a single k-point: this is done in wfcinit, 
  ! ...  directly from file, in order to avoid wasting memory)
  !
  DO ik = 1, ik_
     IF ( nks > 1 .OR. lelfield ) &
        CALL get_buffer ( evc, nwordwfc, iunwfc, ik )
  END DO
  !
  IF ( isolve == 0 ) THEN
     WRITE( stdout, '(5X,"Davidson diagonalization with overlap")' )
  ELSE IF ( isolve == 1 ) THEN
     WRITE( stdout, '(5X,"CG style diagonalization")')
  ELSE
     CALL errore ( 'c_bands', 'invalid type of diagonalization', isolve)
  END IF
  !
  ! ... For each k point diagonalizes the hamiltonian
  !
  k_loop: DO ik = ik_+1, nks
     !
     ! ... Set k-point, spin, kinetic energy, needed by Hpsi
     !
     current_k = ik
     IF ( lsda ) current_spin = isk(ik)
     call g2_kin( ik )
     !
     ! ... More stuff needed by the hamiltonian: nonlocal projectors
     !
     IF ( nkb > 0 ) CALL init_us_2( ngk(ik), igk_k(1,ik), xk(1,ik), vkb )
     !
     ! ... read in wavefunctions from the previous iteration
     !
     IF ( nks > 1 .OR. lelfield ) &
          CALL get_buffer ( evc, nwordwfc, iunwfc, ik )
     !
     ! ... Needed for LDA+U
     !
     IF ( nks > 1 .AND. lda_plus_u .AND. (U_projection .NE. 'pseudo') ) &
          CALL get_buffer ( wfcU, nwordwfcU, iunhub, ik )
     !
     ! ... diagonalization of bands for k-point ik
     !
     call diag_bands ( iter, ik, avg_iter )
     !
     ! ... save wave-functions to be used as input for the
     ! ... iterative diagonalization of the next scf iteration
     ! ... and for rho calculation
     !
     IF ( nks > 1 .OR. lelfield ) &
          CALL save_buffer ( evc, nwordwfc, iunwfc, ik )
     !
     ! ... beware: with pools, if the number of k-points on different
     ! ... pools differs, make sure that all processors are still in
     ! ... the loop on k-points before checking for stop condition
     !
     nkdum  = kunit * ( nkstot / kunit / npool )
     !
     IF (ik .le. nkdum) THEN
        IF (check_stop_now()) THEN
           CALL save_in_cbands(ik, ethr, avg_iter, et )
           RETURN
        END IF
     ENDIF
     !
  END DO k_loop
  !
  CALL mp_sum( avg_iter, inter_pool_comm )
  avg_iter = avg_iter / nkstot
  !
  WRITE( stdout, &
       '( 5X,"ethr = ",1PE9.2,",  avg # of iterations =",0PF5.1 )' ) &
       ethr, avg_iter
  !
  CALL stop_clock( 'c_bands' )
  !
  RETURN
  !
END SUBROUTINE c_bands
!
!----------------------------------------------------------------------------
SUBROUTINE diag_bands( iter, ik, avg_iter )
  !----------------------------------------------------------------------------
  !
  ! ... Driver routine for diagonalization at each k-point
  ! ... Two types of iterative diagonalizations are currently used:
  ! ... a) Davidson algorithm (all-band)
  ! ... b) Conjugate Gradient (band-by-band)
  ! ...
  ! ... internal procedures :
  !
  ! ... diag_bands_gamma(): optimized algorithms for gamma sampling of the BZ
  ! ...                    (real Hamiltonian)
  ! ... diag_bands_k()    : general algorithm for arbitrary BZ sampling
  ! ...                     (complex Hamiltonian)
  ! ... test_exit_cond()  : the test on the iterative diagonalization
  !
  !
  USE kinds,                ONLY : DP
  USE buffers,              ONLY : get_buffer
  USE io_global,            ONLY : stdout
  USE io_files,             ONLY : nwordwfc, iunefieldp, iunefieldm
  USE uspp,                 ONLY : vkb, nkb, okvan
  USE gvect,                ONLY : gstart
  USE wvfct,                ONLY : g2kin, nbndx, et, nbnd, npwx, btype
  USE control_flags,        ONLY : ethr, lscf, max_cg_iter, isolve, &
                                   gamma_only, use_para_diag
  USE noncollin_module,     ONLY : noncolin, npol
  USE wavefunctions_module, ONLY : evc
  USE g_psi_mod,            ONLY : h_diag, s_diag
  USE scf,                  ONLY : v_of_0
  USE bp,                   ONLY : lelfield, evcel, evcelp, evcelm, bec_evcel,&
                                   gdir, l3dstring, efield, efield_cry
  USE becmod,               ONLY : bec_type, becp, calbec, &
                                   allocate_bec_type, deallocate_bec_type
  USE klist,                ONLY : nks, ngk
  USE mp_bands,             ONLY : nproc_bgrp, intra_bgrp_comm, inter_bgrp_comm, &
                                   set_bgrp_indices, my_bgrp_id, nbgrp
  USE mp,                   ONLY : mp_sum, mp_bcast
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN)  :: iter, ik
  !
  REAL (KIND=DP), INTENT(INOUT) :: avg_iter
  !
  REAL (KIND=DP) :: cg_iter
  ! (weighted) number of iterations in Conjugate-Gradient
  INTEGER :: npw, ig, dav_iter, ntry, notconv
  ! number of iterations in Davidson
  ! number or repeated call to diagonalization in case of non convergence
  ! number of notconverged elements
  INTEGER :: ierr, ipw, ibnd, ibnd_start, ibnd_end
  !
  LOGICAL :: lrot
  ! .TRUE. if the wfc have already be rotated
  !
  ALLOCATE( h_diag( npwx, npol ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( ' diag_bands ', ' cannot allocate h_diag ', ABS(ierr) )
  !
  ALLOCATE( s_diag( npwx, npol ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( ' diag_bands ', ' cannot allocate s_diag ', ABS(ierr) )
  !
  ipw=npwx
  CALL mp_sum(ipw, intra_bgrp_comm)
  IF ( nbndx > ipw ) &
     CALL errore ( 'diag_bands', 'too many bands, or too few plane waves',1)
  !
  ! ... allocate space for <beta_i|psi_j> - used in h_psi and s_psi
  !
  CALL allocate_bec_type ( nkb, nbnd, becp, intra_bgrp_comm )
  !
  npw = ngk(ik)
  IF ( gamma_only ) THEN
     !
     CALL diag_bands_gamma()
     !
  ELSE
     !
     CALL diag_bands_k()
     !
  END IF
  !
  ! ... deallocate work space
  !
  CALL deallocate_bec_type ( becp )
  DEALLOCATE( s_diag )
  DEALLOCATE( h_diag )
  !
  IF ( notconv > MAX( 5, nbnd / 4 ) ) THEN
     !
     CALL errore( 'c_bands', &
          & 'too many bands are not converged', 1 )
     !
  ELSE IF ( notconv > 0 ) THEN
     !
     WRITE( stdout, '(5X,"c_bands: ",I2, &
               &   " eigenvalues not converged")' ) notconv
     !
  END IF
  !
  RETURN
  !
CONTAINS
  !
  ! ... internal procedures
  !
  !-----------------------------------------------------------------------
  SUBROUTINE diag_bands_gamma()
    !-----------------------------------------------------------------------
    !
    ! ... Diagonalization of a real Hamiltonian
    !
    IMPLICIT NONE
    !
    IF ( isolve == 1 ) THEN
       !
       ! ... Conjugate-Gradient diagonalization
       !
       ! ... h_diag is the precondition matrix
       !
       FORALL( ig = 1 : npw )
          !
          h_diag(ig,1) = 1.D0 + g2kin(ig) + SQRT( 1.D0 + ( g2kin(ig) - 1.D0 )**2 )
          !
       END FORALL
       !
       ntry = 0
       !
       CG_loop : DO
          !
          lrot = ( iter == 1 .AND. ntry == 0 )
          !
          IF ( .NOT. lrot ) THEN
             !
             CALL rotate_wfc ( npwx, npw, nbnd, gstart, nbnd, evc, npol, okvan, evc, et(1,ik) )
             !
             avg_iter = avg_iter + 1.D0
             !
          END IF
          !
          CALL rcgdiagg( npwx, npw, nbnd, evc, et(1,ik), btype(1,ik), &
               h_diag, ethr, max_cg_iter, .NOT. lscf, notconv, cg_iter )
          !
          avg_iter = avg_iter + cg_iter
          !
          ntry = ntry + 1
          !
          ! ... exit condition
          !
          IF ( test_exit_cond() ) EXIT  CG_loop
          !
       END DO CG_loop
       !
    ELSE
       !
       ! ... Davidson diagonalization
       !
       ! ... h_diag are the diagonal matrix elements of the
       ! ... hamiltonian used in g_psi to evaluate the correction
       ! ... to the trial eigenvectors
       !
       h_diag(1:npw, 1) = g2kin(1:npw) + v_of_0
       !
       CALL usnldiag( npw, h_diag, s_diag )
       !
       ntry = 0
       !
       david_loop: DO
          !
          lrot = ( iter == 1 )
          !
          IF ( use_para_diag ) then
             !
!             ! make sure that all processors have the same wfc
             CALL pregterg( npw, npwx, nbnd, nbndx, evc, ethr, &
                         okvan, gstart, et(1,ik), btype(1,ik), &
                         notconv, lrot, dav_iter )
             !
          ELSE
             !
             CALL regterg ( npw, npwx, nbnd, nbndx, evc, ethr, &
                         okvan, gstart, et(1,ik), btype(1,ik), &
                         notconv, lrot, dav_iter )
          END IF
          !
          avg_iter = avg_iter + dav_iter
          !
          ntry = ntry + 1
          !
          ! ... exit condition
          !
          IF ( test_exit_cond() ) EXIT  david_loop
          !
       END DO david_loop
       !
    END IF
    !
    RETURN
    !
  END SUBROUTINE diag_bands_gamma
  !
  !-----------------------------------------------------------------------
  SUBROUTINE diag_bands_k()
    !-----------------------------------------------------------------------
    !
    ! ... Complex Hamiltonian diagonalization
    !
    IMPLICIT NONE
    !
    ! ... here the local variables
    !
    INTEGER :: ipol
    REAL(dp) :: eps
    !  --- Define a small number ---
    eps=0.000001d0
    !
    IF ( lelfield ) THEN
       !
       ! ... save wave functions from previous iteration for electric field
       !
       evcel = evc
       !
       !... read projectors from disk
       !
       if(.not.l3dstring .and. ABS(efield)>eps ) then
          CALL get_buffer (evcelm(:,:,gdir), nwordwfc, iunefieldm, ik+(gdir-1)*nks)
          CALL get_buffer (evcelp(:,:,gdir), nwordwfc, iunefieldp, ik+(gdir-1)*nks)
       else
          do ipol=1,3
             if(ABS(efield_cry(ipol))>eps) then
                CALL get_buffer (evcelm(:,:,ipol), nwordwfc, iunefieldm, ik+(ipol-1)*nks)
                CALL get_buffer (evcelp(:,:,ipol), nwordwfc, iunefieldp, ik+(ipol-1)*nks)
             endif
          enddo
       endif
       !
       IF ( okvan ) THEN
          !
          call allocate_bec_type(nkb,nbnd,bec_evcel)
          
          !
          CALL calbec(npw, vkb, evcel, bec_evcel)
          !
       ENDIF
       !
    END IF
    !
    IF ( isolve == 1 ) THEN
       !
       ! ... Conjugate-Gradient diagonalization
       !
       ! ... h_diag is the precondition matrix
       !
       h_diag = 1.D0
       !
       FORALL( ig = 1 : npwx )
          !
          h_diag(ig,:) = 1.D0 + g2kin(ig) + SQRT( 1.D0 + ( g2kin(ig) - 1.D0 )**2 )
          !
       END FORALL
       !
       ntry = 0
       !
       CG_loop : DO
          !
          lrot = ( iter == 1 .AND. ntry == 0 )
          !
          IF ( .NOT. lrot ) THEN
             !
             CALL rotate_wfc ( npwx, npw, nbnd, gstart, nbnd, evc, npol, okvan, evc, et(1,ik) )
             !
             avg_iter = avg_iter + 1.D0
             !
          END IF
          !
          CALL ccgdiagg( npwx, npw, nbnd, npol, evc, et(1,ik), btype(1,ik), &
               h_diag, ethr, max_cg_iter, .NOT. lscf, notconv, cg_iter )
          !
          avg_iter = avg_iter + cg_iter
          !
          ntry = ntry + 1
          !
          ! ... exit condition
          !
          IF ( test_exit_cond() ) EXIT  CG_loop
          !
       END DO CG_loop
       !
    ELSE
       !
       ! ... Davidson diagonalization
       !
       ! ... h_diag are the diagonal matrix elements of the
       ! ... hamiltonian used in g_psi to evaluate the correction
       ! ... to the trial eigenvectors
       !
       DO ipol = 1, npol
          !
          h_diag(1:npw, ipol) = g2kin(1:npw) + v_of_0
          !
       END DO
       !
       CALL usnldiag( npw, h_diag, s_diag )
       !
       ntry = 0
       !
       david_loop: DO
          !
          lrot = ( iter == 1 )
          !
          IF ( use_para_diag ) then
             !
             CALL pcegterg( npw, npwx, nbnd, nbndx, npol, evc, ethr, &
                         okvan, et(1,ik), btype(1,ik), &
                         notconv, lrot, dav_iter )
             !
          ELSE
             !
             CALL cegterg ( npw, npwx, nbnd, nbndx, npol, evc, ethr, &
                         okvan, et(1,ik), btype(1,ik), &
                         notconv, lrot, dav_iter )
          END IF
          !
          avg_iter = avg_iter + dav_iter
          !
          ! ... save wave-functions to be used as input for the
          ! ... iterative diagonalization of the next scf iteration
          ! ... and for rho calculation
          !
          ntry = ntry + 1
          !
          ! ... exit condition
          !
          IF ( test_exit_cond() ) EXIT david_loop
          !
       END DO david_loop
       !
    END IF
    !
    IF ( lelfield .AND. okvan ) call deallocate_bec_type( bec_evcel) 
    !
    RETURN
    !
  END SUBROUTINE diag_bands_k
  !
  !-----------------------------------------------------------------------
  FUNCTION test_exit_cond()
    !-----------------------------------------------------------------------
    !
    ! ... this logical function is .TRUE. when iterative diagonalization
    ! ... is converged
    !
    IMPLICIT NONE
    !
    LOGICAL :: test_exit_cond
    !
    !
    test_exit_cond = .NOT. ( ( ntry <= 5 ) .AND. &
         ( ( .NOT. lscf .AND. ( notconv > 0 ) ) .OR. &
         (       lscf .AND. ( notconv > 5 ) ) ) )
    !
  END FUNCTION test_exit_cond
  !
END SUBROUTINE diag_bands
!
!----------------------------------------------------------------------------
SUBROUTINE c_bands_efield ( iter )
  !----------------------------------------------------------------------------
  !
  ! ... Driver routine for Hamiltonian diagonalization under an electric field
  !
  USE noncollin_module,     ONLY : noncolin, npol
  USE kinds,                ONLY : DP
  USE bp,                   ONLY : nberrycyc, fact_hepsi, &
                                   evcel, evcelp, evcelm, gdir, l3dstring,&
                                   efield, efield_cry
  USE klist,                ONLY : nks
  USE wvfct,                ONLY : nbnd, npwx
  USE io_global,            ONLY : stdout
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT (in) :: iter
  !
  INTEGER :: inberry, ipol, ierr
  !
  !
  ALLOCATE( evcel ( npol*npwx, nbnd ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( ' c_bands_efield ', ' cannot allocate evcel ', ABS( ierr ) )
  ALLOCATE( evcelm( npol*npwx, nbnd, 3  ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( ' c_bands_efield ', ' cannot allocate evcelm ', ABS( ierr ) )
  ALLOCATE( evcelp( npol*npwx, nbnd, 3 ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( ' c_bands_efield ', ' cannot allocate evcelp ', ABS( ierr ) )
  ALLOCATE( fact_hepsi(nks, 3), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( ' c_bands_efield ', ' cannot allocate fact_hepsi ', ABS( ierr ) )
  !
  DO inberry = 1, nberrycyc
     !
     !...set up electric field hermitean operator
     !
     FLUSH(stdout)
     if(.not.l3dstring) then
        CALL h_epsi_her_set (gdir, efield)
     else
        do ipol=1,3
           CALL h_epsi_her_set(ipol, efield_cry(ipol))
        enddo
     endif
     FLUSH(stdout)
     !
     CALL c_bands( iter )
     !
  END DO
  !
  DEALLOCATE( fact_hepsi )
  DEALLOCATE( evcelp )
  DEALLOCATE( evcelm )
  DEALLOCATE( evcel  )
  !
  RETURN
  !
END SUBROUTINE c_bands_efield
!
SUBROUTINE c_bands_nscf( )
  !----------------------------------------------------------------------------
  !
  ! ... Driver routine for Hamiltonian diagonalization routines
  ! ... specialized to non-self-consistent calculations (no electric field)
  !
  USE kinds,                ONLY : DP
  USE io_global,            ONLY : stdout
  USE io_files,             ONLY : iunhub, iunwfc, nwordwfc, nwordwfcU
  USE buffers,              ONLY : get_buffer, save_buffer, close_buffer
  USE basis,                ONLY : starting_wfc
  USE klist,                ONLY : nkstot, nks, xk, ngk, igk_k, nelec
  USE uspp,                 ONLY : vkb, nkb
  USE gvect,                ONLY : g, ig_l2g,ngm,ngm_g
  USE wvfct,                ONLY : et, nbnd, npwx, current_k,nbndx, g2kin
  USE becmod,   ONLY : bec_type, becp, calbec, allocate_bec_type, deallocate_bec_type
  USE constants,            ONLY : tpi
  USE cell_base,            ONLY : tpiba,bg
  USE electrons_base,       ONLY : nelt,nel
  USE g_psi_mod,            ONLY : h_diag, s_diag
  USE control_flags,        ONLY : ethr, restart, isolve, io_level, iverbosity
  USE ldaU,                 ONLY : lda_plus_u, U_projection, wfcU
  USE lsda_mod,             ONLY : current_spin, lsda, isk
  USE wavefunctions_module, ONLY : evc
  USE mp_pools,             ONLY : npool, kunit, inter_pool_comm
  USE mp,                   ONLY : mp_sum, mp_bcast
  USE check_stop,           ONLY : check_stop_now
  USE mp_bands,      ONLY : intra_bgrp_comm, inter_bgrp_comm, root_bgrp_id, nbgrp, my_bgrp_id
  !
  IMPLICIT NONE
  !
  REAL(DP) :: avg_iter, ethr_
  ! average number of H*psi products
  INTEGER :: ik_, ik, nkdum, ios, ierr
  ! ik_: k-point already done in a previous run
  ! ik : counter on k points
  LOGICAL :: exst
  !
  REAL(DP), EXTERNAL :: get_clock
  !
  COMPLEX(DP) :: vectx(npwx),vecty(npwx),pupkx(npwx),pupky(npwx),overlapx,overlapy
  COMPLEX (DP), allocatable :: xpsi(:,:,:), xhpsi(:,:,:)
  REAL(DP) :: phaser,phasei,phasery,phaseiy,berrycurv,qx,qy,berrym
  INTEGER :: iiii,iiij,iiix,iiiy,iiiv 
  !
  CALL start_clock( 'c_bands' )
  !
  ik_ = 0
  avg_iter = 0.D0
  IF ( restart ) CALL restart_in_cbands(ik_, ethr, avg_iter, et )
  !
  ! ... If restarting, calculated wavefunctions have to be read from file
  !
  DO ik = 1, ik_
     CALL get_buffer ( evc, nwordwfc, iunwfc, ik )
  END DO
  !
  IF ( isolve == 0 ) THEN
     WRITE( stdout, '(5X,"Davidson diagonalization with overlap")' )
  ELSE IF ( isolve == 1 ) THEN
     WRITE( stdout, '(5X,"CG style diagonalization")')
  ELSE
     CALL errore ( 'c_bands', 'invalid type of diagonalization', isolve)
  END IF
  !
  ! ... For each k point (except those already calculated if restarting)
  ! ... diagonalizes the hamiltonian
  !
  k_loop: DO ik = ik_+1, nks
     !
     ! ... Set k-point, spin, kinetic energy, needed by Hpsi
     !
     current_k = ik
     IF ( lsda ) current_spin = isk(ik)
     WRITE(stdout,*) ik,current_spin
     call g2_kin( ik )
     ! 
     ! ... More stuff needed by the hamiltonian: nonlocal projectors
     !
     IF ( nkb > 0 ) CALL init_us_2( ngk(ik), igk_k(1,ik), xk(1,ik), vkb )
     !
     ! ... Needed for LDA+U
     !
     IF ( nks > 1 .AND. lda_plus_u .AND. (U_projection .NE. 'pseudo') ) &
          CALL get_buffer ( wfcU, nwordwfcU, iunhub, ik )
     !
     ! ... calculate starting  wavefunctions
     !
     IF ( iverbosity > 0 ) WRITE( stdout, 9001 ) ik
     !
     IF ( TRIM(starting_wfc) == 'file' ) THEN
        !
        CALL get_buffer ( evc, nwordwfc, iunwfc, ik )
        !
     ELSE
        !
        CALL init_wfc ( ik )
        !
     END IF
     !
     ! ... diagonalization of bands for k-point ik
     !
     call diag_bands ( 1, ik, avg_iter )
     !
     ! ... save wave-functions (unless disabled in input)
     !
!     IF ( io_level > -1 ) CALL save_buffer ( evc, nwordwfc, iunwfc, ik )
     !
! Xuan Fengyuan modify the code for g factor calculation starts here
     WRITE(stdout,*) 'nelec',nelec
     iiiv=NINT(nelec/2)+0
          WRITE(stdout,*) 'iiiv',iiiv
     IF ( ik .eq. 1 ) THEN
     WRITE( stdout, * ) ' ngk',ngk(ik)
     phaser = 0.0d0
     allocate(xpsi(npwx,1,nbndx), STAT=ierr)
     allocate(xhpsi(npwx,1,nbndx), STAT=ierr)
  xhpsi = ( 0.D0, 0.D0 )
  xpsi  = ( 0.D0, 0.D0 )
  xpsi(:,1,1:nbnd) = evc(:,1:nbnd)
     do iiii = 1,npwx
      vectx(iiii) = evc(iiii,iiiv)
     end do
     WRITE( stdout, * ) ngm,ngm_g,nbnd,nbndx,npwx
     WRITE( stdout, * ) 'XUAN,ngm,ngm_g,nbnd,nbndx,npwx'
     FLUSH(stdout)
     WRITE( stdout, * ) 'xpsi',xpsi(100,1,iiiv),xpsi(1005,1,iiiv)
     FLUSH(stdout)
  CALL allocate_bec_type ( nkb, nbnd, becp, intra_bgrp_comm )
      CALL h_psi( npwx, npwx, nbnd, xpsi, xhpsi )
  CALL deallocate_bec_type ( becp )
     WRITE( stdout, * ) 'check h_psi, xhpsi/xpsi, compare with eigenvalue',xhpsi(100,1,iiiv)/xpsi(100,1,iiiv)
     deallocate (xpsi)
     deallocate (xhpsi)
     END IF
     IF ( ik .eq. 2 ) THEN
     WRITE( stdout, * ) ' ngk',ngk(ik)
     phaser = 0.0d0
     do iiii = 1,npwx
     vecty(iiii) = evc(iiii,iiiv)
     end do
     END IF
     IF ( ik .eq. 3 ) THEN
     WRITE( stdout, * ) ' ngk',ngk(ik)
     phaser = 0.0d0
       berrycurv = 0.0d0
       berrym=0.0d0
!       qx=0.001*sqrt(3.d0)*2*3.14159/(4.651*0.52918)
!       qy=0.001*2*3.14159/(4.651*0.52918)
!BN
       qx=0.001*sqrt(3.d0)*2*3.14159/(4.651)
       qy=0.001*2*3.14159/(4.651)
!MoS2
!       qx=0.0001*1.1434*sqrt(3.d0)*2*3.14159/(5.9723)
!       qy=0.0001*1.1434*2*3.14159/(5.9723)
!WSe2
!       qx=0.0001*1.145218*sqrt(3.d0)*2*3.14159/(6.2172)
!       qy=0.0001*1.145218*2*3.14159/(6.2172)
!general
       qx=0.00001*sqrt(bg(1,2)**2+bg(2,2)**2)*sqrt(3.d0)*tpiba
       qy=0.00001*sqrt(bg(1,2)**2+bg(2,2)**2)*tpiba
!       qy= 0.001*bg(3,3)*tpiba       
!       qx=0.0002*tpiba
!       qy=qx
!       qy=0.0002*sqrt(bg(1,2)**2+bg(2,2)**2+bg(3,2)**2)*tpiba
!       WRITE(stdout,*) qx
!       WRITE(stdout,*) qy
       phaser = 0.0d0
       phasei = 0.0d0
       phasery = 0.0d0
       phaseiy = 0.0d0
     do iiii = 1,npwx
      iiix = 0
      iiiy = 0
      do iiij = 1,npwx
       if (igk_k(iiij,1) .eq. igk_k(iiii,3)) iiix = iiij
       if (igk_k(iiij,2) .eq. igk_k(iiii,3)) iiiy = iiij 
      end do
      if ((iiix .eq. 0) .or. (iiiy .eq. 0) )      WRITE( stdout, * ) 'WRONG'
      FLUSH(stdout)
      phaser = phaser + REAL(evc(iiii,iiiv))*REAL(vectx(iiix))+AIMAG(evc(iiii,iiiv))*AIMAG(vectx(iiix))
      phasei = phasei + REAL(evc(iiii,iiiv))*AIMAG(vectx(iiix))-AIMAG(evc(iiii,iiiv))*REAL(vectx(iiix))
      phasery = phasery + REAL(evc(iiii,iiiv))*REAL(vecty(iiiy))+AIMAG(evc(iiii,iiiv))*AIMAG(vecty(iiiy))
      phaseiy = phaseiy + REAL(evc(iiii,iiiv))*AIMAG(vecty(iiiy))-AIMAG(evc(iiii,iiiv))*REAL(vecty(iiiy))
     end do
     overlapx=CMPLX(phaser/sqrt(phaser**2+phasei**2),-phasei/sqrt(phaser**2+phasei**2))
     overlapy=CMPLX(phasery/sqrt(phasery**2+phaseiy**2),-phaseiy/sqrt(phasery**2+phaseiy**2))
     do iiii = 1,npwx
      iiix = 0
      iiiy = 0
      do iiij = 1,npwx
       if (igk_k(iiij,1) .eq. igk_k(iiii,3)) iiix = iiij
       if (igk_k(iiij,2) .eq. igk_k(iiii,3)) iiiy = iiij
      end do
      if ((iiix .eq. 0) .or. (iiiy .eq. 0) )      WRITE( stdout, * ) 'WRONG'
      FLUSH(stdout)
      pupkx(iiii)=(overlapx*vectx(iiix)-evc(iiii,iiiv))/qx
      pupky(iiii)=(overlapy*vecty(iiiy)-evc(iiii,iiiv))/qy
     end do
     allocate(xpsi(npwx,1,nbndx), STAT=ierr)
     allocate(xhpsi(npwx,1,nbndx), STAT=ierr)
  xhpsi = ( 0.D0, 0.D0 )
  xpsi  = ( 0.D0, 0.D0 )
  xpsi(:,1,1:nbnd) = evc(:,1:nbnd)
  xpsi(:,1,iiiv) = pupky(:)
     WRITE( stdout, * ) 'xpsi',xpsi(100,1,iiiv),xpsi(1005,1,iiiv)
  CALL allocate_bec_type ( nkb, nbnd, becp, intra_bgrp_comm )
      CALL h_psi( npwx, npwx, nbnd, xpsi, xhpsi )
  CALL deallocate_bec_type ( becp )
     WRITE( stdout, * ) 'xhpsi',xhpsi(100,1,iiiv),xhpsi(1005,1,iiiv)
     do iiii = 1,npwx
      berrycurv=berrycurv-2*(REAL(pupkx(iiii))*AIMAG(pupky(iiii))-AIMAG(pupkx(iiii))*REAL(pupky(iiii)))
      berrym=berrym-2*(REAL(pupkx(iiii))*AIMAG(xhpsi(iiii,1,iiiv))-AIMAG(pupkx(iiii))*REAL(xhpsi(iiii,1,iiiv)))
     end do
      berrym=berrym-et(iiiv,ik)*berrycurv
!      berrym=berrym+0.0762*berrycurv
!      berrym=berrym-0.1213*berrycurv
!      berrym=berrym-0.1114*berrycurv
!      berrym=berrym+0.0575*berrycurv
!      berrym=berrym+0.0768*berrycurv
      berrym=-berrym/2.0d0
     WRITE( stdout, * ) 'overlapx',phaser,phasei,phaser**2+phasei**2,overlapx
     WRITE( stdout, * ) 'overlapy',phasery,phaseiy,phasery**2+phaseiy**2,overlapy
     WRITE( stdout, * ) 'eigen value in Ryd',et(iiiv,ik)
     WRITE( stdout, * ) 'berry curvature in bohr**2, in A**2',berrycurv,berrycurv*0.52918**2
     WRITE( stdout, * ) 'orbtial magnetization gfactor unitless',berrym
     phaser = 0.0d0
     do iiii = 1,npwx
        phaser = phaser + g2kin(iiii)*(REAL(evc(iiii,iiiv))**2.0+AIMAG(evc(iiii,iiiv))**2.0)
     end do
     WRITE( stdout, * ) 'kinetic energy of VBM', phaser
      overlapx = 0.0d0
     do iiii = 1,npwx
      overlapx = overlapx+  REAL(pupkx(iiii))*REAL(pupkx(iiii))+AIMAG(pupkx(iiii))*AIMAG(pupkx(iiii))
     end do
     do iiii = 1,npwx
      evc(iiii,iiiv) = pupkx(iiii)/overlapx
     end do
     WRITE( stdout, * ) 'normalization of dudk', overlapx
     deallocate (xpsi)
     deallocate (xhpsi)
     END IF
! end of the modification for g factor calculation
     !
     IF ( io_level > -1 ) CALL save_buffer ( evc, nwordwfc, iunwfc, ik )
     ! ... beware: with pools, if the number of k-points on different
     ! ... pools differs, make sure that all processors are still in
     ! ... the loop on k-points before checking for stop condition
     !
     nkdum  = kunit * ( nkstot / kunit / npool )
     IF (ik .le. nkdum) THEN
        !
        ! ... stop requested by user: save restart information,
        ! ... save wavefunctions to file
        !
        IF (check_stop_now()) THEN
           CALL save_in_cbands(ik, ethr, avg_iter, et )
           RETURN
        END IF
     ENDIF
     !
     ! report about timing
     !
     IF ( iverbosity > 0 ) THEN
        WRITE( stdout, 9000 ) get_clock( 'PWSCF' )
        FLUSH( stdout )
     ENDIF
     !
  END DO k_loop
  !
  CALL mp_sum( avg_iter, inter_pool_comm )
  avg_iter = avg_iter / nkstot
  !
  WRITE( stdout, '(/,5X,"ethr = ",1PE9.2,",  avg # of iterations =",0PF5.1)' ) &
       ethr, avg_iter
  !
  CALL stop_clock( 'c_bands' )
  !
  RETURN
  !
  ! formats
  !
9001 FORMAT(/'     Computing kpt #: ',I5 )
9000 FORMAT( '     total cpu time spent up to now is ',F10.1,' secs' )
  !
END SUBROUTINE c_bands_nscf
