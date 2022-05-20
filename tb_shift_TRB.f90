         SUBROUTINE tb_shift_TRB

         USE tb_module , ONLY : nk,kp,Hmnk,et,num_wann,alat, &
                                wsgp,deg_rpts,cell,nrpts,Hmnk, &
                                Hmn,pi,bohr2ang,BC,fermi,dk,&
                                shift_de,shift_a,shift_grid,volume,&
                                ev2ry,ry2ev,ang2bohr,ion,nat,at,myrank,&
                                nprocs,ierr,tb_eta,nw_s,nw_e,bk,kp_cart,soc,&
                                pi2zi,zi

         IMPLICIT NONE
         COMPLEX(kind=8) :: U(num_wann,num_wann),Ud(num_wann,num_wann)
         COMPLEX(kind=8) :: vx(num_wann,num_wann)
         COMPLEX(kind=8) :: vy(num_wann,num_wann)
         COMPLEX(kind=8) :: vz(num_wann,num_wann)
         COMPLEX(kind=8) :: mat_tmp(num_wann,num_wann)
         COMPLEX(kind=8) :: vmat(num_wann,num_wann,3),vmat_t(num_wann,num_wann,nk,3)
         COMPLEX(kind=8), ALLOCATABLE :: Bmat(:,:)
         COMPLEX(kind=8), ALLOCATABLE :: nsc(:,:,:,:), msc(:,:,:,:), shift_copy(:,:,:,:)
         COMPLEX(kind=8), ALLOCATABLE :: ee(:,:,:), ee_copy(:,:,:)
         COMPLEX(kind=8), ALLOCATABLE :: nic(:,:,:,:), mic(:,:,:,:), inject_copy(:,:,:,:)
         COMPLEX(kind=8) :: polar(3)
         COMPLEX(kind=8) :: ee_nm(3,3),nsc_nm(3,3,3),msc_nm(3,3,3)
         COMPLEX(kind=8) :: nic_nm(3,3,3), mic_nm(3,3,3)
         COMPLEX(kind=8) :: vmat_t2(num_wann,num_wann,nk,3)
         COMPLEX(kind=8) :: wmat(num_wann,num_wann,3,3)
         COMPLEX(kind=8) :: vdmat(num_wann,num_wann,3)
         COMPLEX(kind=8) :: amat(num_wann,num_wann,3)
         COMPLEX(kind=8) :: rmat(num_wann,num_wann,3,3)
         COMPLEX(kind=8) :: ratio
         REAL(kind=8) :: delta_e(num_wann,num_wann),del_e1,del_e2,del_e
         REAL(kind=8) :: ion_gauge(3)
         REAL(kind=8) :: ion_pol(3),tau1(3),tau2(3)
         REAL(kind=8) :: r_cart(3),r(3),kdotr,wg(num_wann)
         INTEGER :: iwann,irpts,jwann,igrid,kwann
         INTEGER :: ik,ipol,jpol,kpol,nwann
         INTEGER :: t1s,t1e,t2s,t2e,iat,jat
         INTEGER :: ivpt(3), info

         INCLUDE 'mpif.h'

         ALLOCATE( nsc(3,3,3,shift_grid), msc(3,3,3,shift_grid), shift_copy(3,3,3,shift_grid) )
         ALLOCATE( ee(3,3,shift_grid), ee_copy(3,3,shift_grid) )
         ALLOCATE( nic(3,3,3,shift_grid), mic(3,3,3,shift_grid), inject_copy(3,3,3,shift_grid) )
         ALLOCATE( Bmat(num_wann, num_wann) )
         et =0.0d0
         vmat_t = CMPLX(0.0d0,0.0d0)
         vmat_t2 = CMPLX(0.0d0,0.0d0)
         polar = CMPLX(0.0d0,0.0d0)
         ion_pol = 0.0d0
         ion_gauge = CMPLX(0.0d0,0.0d0)
         ee = CMPLX(0.0d0,0.0d0)
         nsc = CMPLX(0.0d0,0.0d0)
         msc = CMPLX(0.0d0,0.0d0)
         nic = CMPLX(0.0d0,0.0d0)
         mic = CMPLX(0.0d0,0.0d0)
         Hmnk = CMPLX(0.0d0,0.0d0) 
         Bmat = CMPLX(0.0d0,0.0d0)
         Hmn=Hmn*ev2ry
         CALL tb_gen_ham

     IF( myrank .EQ. 0 ) THEN
      WRITE(6,*) 
      WRITE(6,*) '===== shift current start ====='
      WRITE(6,1300) 'SMEARING FACTOR :', shift_a
      WRITE(6,1300) 'Delta E :', shift_de
      WRITE(6,1301) '# of GRIDS :', shift_grid
      WRITE(6,*) 'The unit of shift current is microA/V^2'
      WRITE(6,*) '========================='
      WRITE(6,*) 
     ENDIF

         DO ik = 1+myrank, nk, nprocs
           write(6,100) 'SHIFT-CURRENT K(',ik,') = ',kp(1:3,ik) 
           vmat = CMPLX(0.0d0,0.0d0) 
           wmat = dCMPLX(0.0d0,0.0d0)
           vdmat = CMPLX(0.0d0,0.0d0) 
           amat = CMPLX(0.0d0,0.0d0) 
           rmat = CMPLX(0.0d0,0.0d0) 
           delta_e = 0.0d0
           wg = 0.0d0
           vx = CMPLX(0.0d0,0.0d0) 
           vy = CMPLX(0.0d0,0.0d0) 
           vz = CMPLX(0.0d0,0.0d0) 
           U  = CMPLX(0.0d0,0.0d0) 
           Ud = CMPLX(0.0d0,0.0d0) 
           ee_nm = CMPLX(0.0d0,0.0d0)

           CALL tb_diag('V','U',num_wann,Hmnk(:,:,ik),et(:,ik))
  
           U(:,:) = Hmnk(:,:,ik)
           Ud(:,:) = CONJG(TRANSPOSE(U(:,:)))
 
           IF (soc .GT. 0) THEN
             nwann = num_wann / 2
             DO iat = 1, nat
               tau1(1:3) = at(1:3,iat) 
               t1s = nw_s(iat)
               t1e = nw_e(iat)
               DO jat = 1, nat
                 tau2(1:3) = at(1:3,jat)
                 t2s = nw_s(jat)
                 t2e = nw_e(jat)
                 DO irpts = 1, nrpts 
                     r(1) = dble(wsgp(irpts,1) + tau2(1) - tau1(1))
                     r(2) = dble(wsgp(irpts,2) + tau2(2) - tau1(2))
                     r(3) = dble(wsgp(irpts,3) + tau2(3) - tau1(3))     
                     kdotr = (kp(1,ik)*r(1)+ kp(2,ik)*r(2)+ kp(3,ik)*r(3))
                     ratio = exp(pi2zi*kdotr)/dble(deg_rpts(irpts))
                     r_cart(1) = cell(1,1)*r(1)+ cell(1,2)*r(2)+cell(1,3)*r(3)
                     r_cart(2) = cell(2,1)*r(1)+ cell(2,2)*r(2)+cell(2,3)*r(3)
                     r_cart(3) = cell(3,1)*r(1)+ cell(3,2)*r(2)+cell(3,3)*r(3)
                     r_cart = r_cart * alat 

                     vx(t1s:t1e,t2s:t2e)=vx(t1s:t1e,t2s:t2e)+ zi*r_cart(1) &
                                    *Hmn(t1s:t1e,t2s:t2e,irpts)*ratio
                     vx(t1s:t1e,t2s+nwann:t2e+nwann)=vx(t1s:t1e,t2s+nwann:t2e+nwann)+ zi*r_cart(1) &
                                    *Hmn(t1s:t1e,t2s+nwann:t2e+nwann,irpts)*ratio
                     vx(t1s+nwann:t1e+nwann,t2s:t2e)=vx(t1s+nwann:t1e+nwann,t2s:t2e)+ zi*r_cart(1) &
                                    *Hmn(t1s+nwann:t1e+nwann,t2s:t2e,irpts)*ratio
                     vx(t1s+nwann:t1e+nwann,t2s+nwann:t2e+nwann)=vx(t1s+nwann:t1e+nwann,t2s+nwann:t2e+nwann)+ zi*r_cart(1) &
                                    *Hmn(t1s+nwann:t1e+nwann,t2s+nwann:t2e+nwann,irpts)*ratio

                     vy(t1s:t1e,t2s:t2e)=vy(t1s:t1e,t2s:t2e)+ zi*r_cart(2) &
                                    *Hmn(t1s:t1e,t2s:t2e,irpts)*ratio
                     vy(t1s:t1e,t2s+nwann:t2e+nwann)=vy(t1s:t1e,t2s+nwann:t2e+nwann)+ zi*r_cart(2) &
                                    *Hmn(t1s:t1e,t2s+nwann:t2e+nwann,irpts)*ratio
                     vy(t1s+nwann:t1e+nwann,t2s:t2e)=vy(t1s+nwann:t1e+nwann,t2s:t2e)+ zi*r_cart(2) &
                                    *Hmn(t1s+nwann:t1e+nwann,t2s:t2e,irpts)*ratio
                     vy(t1s+nwann:t1e+nwann,t2s+nwann:t2e+nwann)=vy(t1s+nwann:t1e+nwann,t2s+nwann:t2e+nwann)+ zi*r_cart(2) &
                                    *Hmn(t1s+nwann:t1e+nwann,t2s+nwann:t2e+nwann,irpts)*ratio

                     vz(t1s:t1e,t2s:t2e)=vz(t1s:t1e,t2s:t2e)+ zi*r_cart(3) &
                                    *Hmn(t1s:t1e,t2s:t2e,irpts)*ratio
                     vz(t1s:t1e,t2s+nwann:t2e+nwann)=vz(t1s:t1e,t2s+nwann:t2e+nwann)+ zi*r_cart(3) &
                                    *Hmn(t1s:t1e,t2s+nwann:t2e+nwann,irpts)*ratio
                     vz(t1s+nwann:t1e+nwann,t2s:t2e)=vz(t1s+nwann:t1e+nwann,t2s:t2e)+ zi*r_cart(3) &
                                    *Hmn(t1s+nwann:t1e+nwann,t2s:t2e,irpts)*ratio
                     vz(t1s+nwann:t1e+nwann,t2s+nwann:t2e+nwann)=vz(t1s+nwann:t1e+nwann,t2s+nwann:t2e+nwann)+ zi*r_cart(3) &
                                    *Hmn(t1s+nwann:t1e+nwann,t2s+nwann:t2e+nwann,irpts)*ratio
                 ENDDO
               ENDDO
             ENDDO
           ELSE
             DO iat = 1, nat
               tau1(1:3) = at(1:3,iat) 
               t1s = nw_s(iat)
               t1e = nw_e(iat)
               DO jat = 1, nat
                 tau2(1:3) = at(1:3,jat)
                 t2s = nw_s(jat)
                 t2e = nw_e(jat)
                 DO irpts = 1, nrpts 
                     r(1) = dble(wsgp(irpts,1) + tau2(1) - tau1(1))
                     r(2) = dble(wsgp(irpts,2) + tau2(2) - tau1(2))
                     r(3) = dble(wsgp(irpts,3) + tau2(3) - tau1(3))     
                     kdotr = (kp(1,ik)*r(1)+ kp(2,ik)*r(2)+ kp(3,ik)*r(3))
                     ratio = exp(pi2zi*kdotr)/dble(deg_rpts(irpts))
                     r_cart(1) = cell(1,1)*r(1)+ cell(1,2)*r(2)+cell(1,3)*r(3)
                     r_cart(2) = cell(2,1)*r(1)+ cell(2,2)*r(2)+cell(2,3)*r(3)
                     r_cart(3) = cell(3,1)*r(1)+ cell(3,2)*r(2)+cell(3,3)*r(3)
                     r_cart = r_cart * alat 

                     vx(t1s:t1e,t2s:t2e)=vx(t1s:t1e,t2s:t2e)+ zi*r_cart(1) &
                                    *Hmn(t1s:t1e,t2s:t2e,irpts)*ratio

                     vy(t1s:t1e,t2s:t2e)=vy(t1s:t1e,t2s:t2e)+ zi*r_cart(2) &
                                    *Hmn(t1s:t1e,t2s:t2e,irpts)*ratio

                     vz(t1s:t1e,t2s:t2e)=vz(t1s:t1e,t2s:t2e)+ zi*r_cart(3) &
                                    *Hmn(t1s:t1e,t2s:t2e,irpts)*ratio
                 ENDDO
               ENDDO
             ENDDO
           ENDIF


           call mat_mul(num_wann, vx, U, Bmat)
           call mat_mul(num_wann, Ud, Bmat, vx)
           call mat_mul(num_wann, vy, U, Bmat)
           call mat_mul(num_wann, Ud, Bmat, vy)
           call mat_mul(num_wann, vz, U, Bmat)
           call mat_mul(num_wann, Ud, Bmat, vz)


           vmat(:,:,1) = (vx(:,:))
           vmat(:,:,2) = (vy(:,:))
           vmat(:,:,3) = (vz(:,:))
           vmat_t(:,:,ik,1) = (vx(:,:))
           vmat_t(:,:,ik,2) = (vy(:,:))
           vmat_t(:,:,ik,3) = (vz(:,:))



!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!wmat

           IF (soc .GT. 0) THEN
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
                   kdotr = (kp(1,ik)*r(1)+ kp(2,ik)*r(2)+ kp(3,ik)*r(3))
                   ratio = exp(pi2zi*kdotr)/dble(deg_rpts(irpts))

                   r_cart(1) = cell(1,1)*r(1)+ cell(1,2)*r(2)+cell(1,3)*r(3)
                   r_cart(2) = cell(2,1)*r(1)+ cell(2,2)*r(2)+cell(2,3)*r(3)
                   r_cart(3) = cell(3,1)*r(1)+ cell(3,2)*r(2)+cell(3,3)*r(3)
                   r_cart = r_cart * alat
                   
                   DO ipol = 1,3
                     DO jpol = 1,3
                       wmat(t1s:t1e,t2s:t2e,ipol,jpol)=wmat(t1s:t1e,t2s:t2e,ipol,jpol)- &
                                                              r_cart(ipol)*r_cart(jpol) &
                                                              *Hmn(t1s:t1e,t2s:t2e,irpts)*ratio
                       wmat(t1s:t1e,t2s+nwann:t2e+nwann,ipol,jpol)=wmat(t1s:t1e,t2s+nwann:t2e+nwann,ipol,jpol)- &
                                                              r_cart(ipol)*r_cart(jpol) &
                                                              *Hmn(t1s:t1e,t2s+nwann:t2e+nwann,irpts)*ratio
                       wmat(t1s+nwann:t1e+nwann,t2s:t2e,ipol,jpol)=wmat(t1s+nwann:t1e+nwann,t2s:t2e,ipol,jpol)- &
                                                              r_cart(ipol)*r_cart(jpol) &
                                                              *Hmn(t1s+nwann:t1e+nwann,t2s:t2e,irpts)*ratio
                       wmat(t1s+nwann:t1e+nwann,t2s+nwann:t2e+nwann,ipol,jpol)=wmat(t1s+nwann:t1e+nwann,t2s+nwann:t2e+nwann,ipol,jpol)- &
                                                              r_cart(ipol)*r_cart(jpol) &
                                                              *Hmn(t1s+nwann:t1e+nwann,t2s+nwann:t2e+nwann,irpts)*ratio
                     ENDDO
                   ENDDO
                 ENDDO
               ENDDO
             ENDDO
           ELSE
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
                   kdotr = (kp(1,ik)*r(1)+ kp(2,ik)*r(2)+ kp(3,ik)*r(3))
                   ratio = exp(pi2zi*kdotr)/dble(deg_rpts(irpts))

                   r_cart(1) = cell(1,1)*r(1)+ cell(1,2)*r(2)+cell(1,3)*r(3)
                   r_cart(2) = cell(2,1)*r(1)+ cell(2,2)*r(2)+cell(2,3)*r(3)
                   r_cart(3) = cell(3,1)*r(1)+ cell(3,2)*r(2)+cell(3,3)*r(3)
                   r_cart = r_cart * alat
                   
                   DO ipol = 1,3
                     DO jpol = 1,3
                       wmat(t1s:t1e,t2s:t2e,ipol,jpol)=wmat(t1s:t1e,t2s:t2e,ipol,jpol)- &
                                                              r_cart(ipol)*r_cart(jpol) &
                                                              *Hmn(t1s:t1e,t2s:t2e,irpts)*ratio
                     ENDDO
                   ENDDO
                 ENDDO
               ENDDO
             ENDDO
           ENDIF


           DO ipol = 1,3
             DO jpol = 1,3
               call mat_mul(num_wann, wmat(:,:,ipol,jpol), U, Bmat)
               call mat_mul(num_wann, Ud, Bmat, wmat(:,:,ipol,jpol))
             ENDDO
           ENDDO


        END PROGRAM tb_shift_TRB
