!----------------------------------------------------------------------
program main
  !-----------------------------------------------------------------------
  use tb_module

  implicit none
  integer :: is,n,m,ipol,na,ik

  call tb_read()
  call tb_read_hr()
  call tb_init()

  IF (lk_uniform)  THEN
     call tb_gen_kp_uniform()
  ELSEIF (lk_line) THEN 
     call tb_gen_kp()
  ENDIF

  call tb_eigen()
  IF (lpdos) THEN
    call tb_dos()
  ENDIF
  IF (lBC) THEN
    call tb_BC()
  ENDIF

! call tb_Floquet()
  IF (lw) THEN
    call tb_shift_TRB()
    IF (lspin) THEN
      call tb_spin_shift_TRB()
    ENDIF
   !call tb_shift()
 !ELSEIF (lk) THEN 
 !  IF (myrank .eq. 0) write(6,*) 'start tb_shift_k'
 !  call tb_shift_k()
 !ELSEIF (lvc) THEN 
 !  IF (myrank .eq. 0) write(6,*) 'start tb_shift_vc'
 !  call tb_shift_vc()
  ELSEIF (lkvc) THEN 
    call tb_shift_kvc()
  ENDIF

  STOP
end program main
