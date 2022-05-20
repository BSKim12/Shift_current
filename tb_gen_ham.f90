         SUBROUTINE tb_gen_ham

         USE tb_module , ONLY : num_wann,nrpts,deg_rpts, &
                                wsgp,Hmn,nk,kp,Hmnk,nwat, &
                                nat,nw_s,nw_e,at,cell,pi, & 
                                et,myrank,nprocs,ierr,pi2zi,soc
         IMPLICIT NONE
         INTEGER ::ir,irpts,iat,jat,jk,t1s,t1e,t2s,t2e,iwann,jwann,ik,nwann
         REAL :: kdotr,r(3),tau1(3),tau2(3),deltaH
         COMPLEX(kind=8) :: ratio

         INCLUDE 'mpif.h'

         Hmnk = CMPLX(0.0d0,0.0d0)

         IF (soc .GT. 0) THEN
           nwann = num_wann / 2
           DO ik = 1+myrank, nk, nprocs
             DO iat = 1, nat
               tau1(1:3) = at(1:3,iat) 
               t1s = nw_s(iat)
               t1e = nw_e(iat)
               DO jat = 1, nat
                 tau2(1:3) = at(1:3,jat)
                 t2s = nw_s(jat)
                 t2e = nw_e(jat)
                 DO irpts = 1, nrpts 
                   r(1) = dble(wsgp(irpts,1)) + tau2(1) - tau1(1)
                   r(2) = dble(wsgp(irpts,2)) + tau2(2) - tau1(2)
                   r(3) = dble(wsgp(irpts,3)) + tau2(3) - tau1(3)     
                   kdotr = (kp(1,ik)*r(1)+ kp(2,ik)*r(2)+kp(3,ik)*r(3))
                   ratio = exp(pi2zi*kdotr)/dble(deg_rpts(irpts))
                   Hmnk(t1s:t1e,t2s:t2e,ik) = Hmnk(t1s:t1e,t2s:t2e,ik)+ &
                                              Hmn(t1s:t1e,t2s:t2e,irpts)*ratio
                   Hmnk(t1s:t1e,t2s+nwann:t2e+nwann,ik) = Hmnk(t1s:t1e,t2s+nwann:t2e+nwann,ik)+ &
                                              Hmn(t1s:t1e,t2s+nwann:t2e+nwann,irpts)*ratio
                   Hmnk(t1s+nwann:t1e+nwann,t2s:t2e,ik) = Hmnk(t1s+nwann:t1e+nwann,t2s:t2e,ik)+ &
                                              Hmn(t1s+nwann:t1e+nwann,t2s:t2e,irpts)*ratio
                   Hmnk(t1s+nwann:t1e+nwann,t2s+nwann:t2e+nwann,ik) = Hmnk(t1s+nwann:t1e+nwann,t2s+nwann:t2e+nwann,ik)+ &
                                              Hmn(t1s+nwann:t1e+nwann,t2s+nwann:t2e+nwann,irpts)*ratio
                 ENDDO
               ENDDO
             ENDDO
           ENDDO

         ELSE
           DO ik = 1+myrank, nk, nprocs
             DO iat = 1, nat
               tau1(1:3) = at(1:3,iat) 
               t1s = nw_s(iat)
               t1e = nw_e(iat)
               DO jat = 1, nat
                 tau2(1:3) = at(1:3,jat)
                 t2s = nw_s(jat)
                 t2e = nw_e(jat)
                 DO irpts = 1, nrpts 
                   r(1) = dble(wsgp(irpts,1)) + tau2(1) - tau1(1)
                   r(2) = dble(wsgp(irpts,2)) + tau2(2) - tau1(2)
                   r(3) = dble(wsgp(irpts,3)) + tau2(3) - tau1(3)     
                   kdotr = (kp(1,ik)*r(1)+ kp(2,ik)*r(2)+kp(3,ik)*r(3))
                   ratio = exp(pi2zi*kdotr)/dble(deg_rpts(irpts))
                   Hmnk(t1s:t1e,t2s:t2e,ik) = Hmnk(t1s:t1e,t2s:t2e,ik)+ &
                                              Hmn(t1s:t1e,t2s:t2e,irpts)*ratio
                 ENDDO
               ENDDO
             ENDDO
           ENDDO
         ENDIF !> SOC



         CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)
         DO iwann = 1,num_wann 
           DO jwann = 1,num_wann 
             DO ik = 1+myrank, nprocs, nk
               deltaH = abs(Hmnk(iwann,jwann,ik)-conjg(Hmnk(jwann,iwann,ik)))
               IF( deltaH.GT.1.0d-4 ) WRITE(6,*) 'non-Hermit', iwann,jwann,ik,deltaH
             ENDDO
           ENDDO
         ENDDO

         CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)



         END SUBROUTINE tb_gen_ham
