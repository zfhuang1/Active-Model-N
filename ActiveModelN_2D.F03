! Active Model N, incorporating nonreciprocal force interactions and torque
! 2D particle density field psi (i.e., rho)
! local polarization vector field (Px, Py)

! use a numerical algorithm with exponential propagation procedure, for conserved psi=psi0+dpsi
!     (exact solution for linear terms in Fourier space and approximation for nonlinearities, with initial transient time)
! pseudospectral method (using FFTW3)
  
! using OpenMP (To use single-core only: comment out "define open_MP")
! double precision

! Here B1=A1=>B1, B7=A7=>B2, B12=A12=>B3, B14=A14=>B4, B20=A20=>B5 of the paper

# define P0_nonzero ! initial nonzero P0, with Px0^2 + Py0^2 = P0_s^2
!# define P0_zero ! initial Px0=Py0=0

!# define save_conf ! save conf files

# define open_MP ! use openMP for multi-threading


module global_variables
! by default all variables in the module are public (if not being set explicitly as private)
  use, intrinsic :: iso_c_binding
  implicit none
  include 'fftw3.f03'
  
  save ! save all the variables in the module

  integer, parameter :: ly=512, lx=512

  real(C_DOUBLE), pointer :: psi(:,:),Px(:,:),Py(:,:)
  complex(C_DOUBLE_COMPLEX), pointer :: psiq(:,:),Pxq(:,:),Pyq(:,:)
  real(8) :: psi0,pi,B1,B7_2,B12,B14,B20_4,B1_2,B1B20,B1B20_4,Dr2_B20,B1_2_4Dr
  real(8) :: qx(lx),qy(ly/2+1),qx2(lx),qy2(ly/2+1)
  double precision, allocatable :: q2(:,:),qxy(:,:)
  double complex, allocatable :: lambda(:,:,:,:),chi(:,:,:,:)
  double complex, allocatable :: alpha_ij(:,:,:,:) ! coefficients of 3 PDEs
  double complex, allocatable :: sigma(:,:,:)
  integer, allocatable :: isigma(:,:)

! for FFTW parameters
  type(C_PTR) :: plan2d_f,plan2d_b
  real(8) :: scale2d_b

  character(:), allocatable :: file ! to be automatically allocated

end module global_variables


Program active_model_N_2D
  use global_variables
  implicit none

  type(C_PTR) :: p_alloc_psi,p_alloc_psiq ! for allocating aligned memory in fftw
  double precision, allocatable :: P2(:,:)
  real(8) :: dt,dti,dtime
  real(8) :: dx,dy,Px0,Py0,P02,P0_s,P0
  real(8) :: vartheta,alpha,F0,F1_F0,b,c,Dr,B7,B20
  real(8) :: B1psi0,B1_B10q2
  real(8) :: time,tmax,time0,t_i
  real(8) :: noise,noise0,ran3,yy
  real(8) :: aa,bb,cc,qq,Delta_
  real(8) :: ave_Px,ave_Py,ave_P
  integer :: i,j,iter,idum
  integer :: nend,nout,nt,raout,nout_conf,n_i
  integer :: iter0,n_0
  
# ifdef open_MP
  ! use multi-threaded fftw with openMP
  integer iret,nthreads,omp_get_max_threads
# endif
  
! vartheta=7*pi/9 (140^o)
  file='ahd512al2F01r0001b003c2Dr01theta140psi15Px0Py0' ! 512x512,def P0_nonzero,Px0=Py0=P0/sqrt(2),dx=1,dt=0.05,tmax=20000
  !file='ahd512al2F01r0001b003c2Dr01theta140psi15Px0Py0_0' ! 512x512,def P0_zero,Px0=Py0=0,dx=1,dt=0.02,tmax=20000

! parameters entering the Eqs
!
  psi0=1.5d0 ! average density
  Px0=0
  Py0=0

  pi=DACOS(-1.0d0)
  vartheta=7*pi/9 ! =140^o
  
  Dr=0.1d0
  
! for h(r) = F0 exp(-r^2/a^2) - F1 exp(-r/alpha*a)
  alpha=2
  
  F0=0.1d0
  F1_F0=1.d-3 ! r=F1/F0; alpha=2
  
  b=0.03d0
  c=2
  
  B1=pi*F0*(1-2*alpha*alpha*F1_F0)*sin(vartheta)
  B7=0.5d0*pi*F0*(sqrt(pi)-8*F1_F0*alpha**3)*vartheta
  B12=(pi/8)*F0*(1-12*F1_F0*alpha**4)*sin(vartheta)
  B14=(3*pi/32)*F0*(sqrt(pi)-64*F1_F0*alpha**5)*vartheta
  B20=-2*pi*b*c*c*vartheta

  if(abs(vartheta-pi).lt.1.e-10) then ! if vartheta=pi
     B1=0
     B12=0
  endif
  
  P02=-(4/B20)*(psi0+1/(2*B20)) ! = Px0^2+Py0^2 = P0_s^2
  if(P02 < 0) then
     P0_s=0
     Px0=0
     Py0=0
  else
     P0_s=sqrt(P02)
# ifdef P0_nonzero
     Py0=P0_s/sqrt(2.d0)
     Px0=P02-Py0*Py0
     if(Px0 < 0) then
        Px0=0
        if(abs(abs(Py0)-P0_s) > 1.e-10) Py0=P0_s
     else
        Px0=sqrt(Px0)
     endif
# endif
  endif
  
# ifdef P0_zero
  ! initial Px0=Py0=0
  Px0=0
  Py0=0
# endif
  
  B1_2=2*B1
  B7_2=B7/2
  B1B20=B1*B20
  B1B20_4=B1B20/4 ! =B1*B20/4
  Dr2_B20=2*Dr*B20
  B20_4=B20/4
  B1_2_4Dr=B1*B1/(4*Dr)
  B1psi0=B1*psi0
  
  noise=0.1d0 ! initial random noise (percentage)
!  noise=0
  noise0=noise

  dx=1
  dy=dx

  time0=0
  iter0=1
  t_i=1 ! time of initial transient; with t_i/dt=integer
  dti=0.01d0 ! used for initial transient up to t=t_i

! discretization of time
!  dt=0.02d0
  dt=0.05d0

!  tmax=1000 ! tmax: maximum time
!  nt=50
  tmax=20000
  nt=500

  n_i=NINT(t_i/dti)
  n_0=NINT(t_i/dt)-n_i ! used for output
  nend=NINT((tmax-t_i)/dt)+n_i
  nout=NINT(nt/dt)
  raout=nout/10 ! for output of intermediate results
  if(raout.lt.1) raout=1
  nout_conf=nout*5 ! for saving conf files
  
! fft related items: qx and qy
  do j=1,ly/2+1
     qy(j)=2*pi*(j-1)/(ly*dy)
! use finite difference in Fourier space for discrete derivatives
     qy2(j)=(2/(dy*dy))*(1-cos(qy(j)*dy))
     !qy2(j)=qy(j)*qy(j)
     qy(j)=sin(qy(j)*dy)/dy
  enddo
  
  do i=1,lx
     if (i.le.lx/2+1) then
        qx(i)=2*pi*(i-1)/(lx*dx)
     else
        qx(i)=-2*pi*(lx-i+1)/(lx*dx)
     endif
! use finite difference in Fourier space for discrete derivatives
     qx2(i)=(2/(dx*dx))*(1-cos(qx(i)*dx))
     !qx2(i)=qx(i)*qx(i)
     qx(i)=sin(qx(i)*dx)/dx
  enddo

  idum=-796921 ! random number seed; <0 for ran3()

! write parameters in a separate file
!
  open(unit=2,file=file//'_para.dat')
  write(2,*) 'use a pseudospectral method for psi, Px, Py active nonreciprocal equations'
  write(2,*) 'use an exponential propagation procedure in Fourier space'
  write(2,*) '    exact solution for linear terms in Fourier space and approximation for nonlinearities'
  write(2,*) '    expanding nonlinearities up to 1st order of t'
  write(2,*)
  write(2,*) 'lx=',lx, '  ly=',ly
  write(2,*) 'psi0=',psi0, '  Dr=',Dr
# ifdef P0_zero
  write(2,*) 'def P0_zero: initial Px0=Py0=0'
# endif
# ifdef P0_nonzero
  write(2,*) 'def P0_nonzero: initial nonzero P0, with Px0^2 + Py0^2 = P0_s^2'
# endif
  write(2,*) 'Px0=',Px0, '  Py0=',Py0
  write(2,*) '|P|_s = ', P0_s ! =sqrt(-(4/B20)*(psi0+1/(2*B20))) or =0
  
  write(2,*) 'vartheta=',vartheta, '  vartheta(^o)=',vartheta*180/pi
  write(2,*) 'alpha=',alpha
  write(2,*) 'F0=',F0, '  F1/F0=',F1_F0
  write(2,*) 'b=',b, '  c=',c 
  write(2,*) 'B1=',B1, '  B7=',B7
  write(2,*) 'B14=',B14, '  B20=',B20
  write(2,*) 'B12=',B12
  write(2,*)
  write(2,*) 'dt=',dt, '  dti=',dti
  write(2,*) 'dx=',dx, '  dy=',dy
  write(2,*) 'tmax=',tmax, '  t_transient=',t_i
  write(2,*) 'nend=',nend, '  n_transient=',n_i, '  raout=',raout
  write(2,*) 'nout=',nout, '  nt=',nt
  write(2,*) '  data output at t=t_transient, and every ',raout*dt
# ifdef save_conf
  write(2,*) '  save conf files every ',nout_conf*dt
# endif
  write(2,*) 'time0=',time0, '  iter0=',iter0, '  idum=',idum
  write(2,*) 'noise(%)=',noise,' noise0=',noise0
  write(2,*)
  write(2,*) 'dynamics without noise (only initial noise)'
# ifdef open_MP 
  write(2,*) 'use openMP for multi-threading'
# else
  write(2,*) 'no openMP'
# endif
  close(2)

  p_alloc_psi=fftw_alloc_real(int(ly * lx, C_SIZE_T))
  call c_f_pointer(p_alloc_psi, psi, [ly,lx])
  p_alloc_psi=fftw_alloc_real(int(ly * lx, C_SIZE_T))
  call c_f_pointer(p_alloc_psi, Px, [ly,lx])
  p_alloc_psi=fftw_alloc_real(int(ly * lx, C_SIZE_T))
  call c_f_pointer(p_alloc_psi, Py, [ly,lx])
  p_alloc_psiq=fftw_alloc_complex(int((ly/2+1) * lx, C_SIZE_T))
  call c_f_pointer(p_alloc_psiq, psiq, [ly/2+1,lx])
  p_alloc_psiq=fftw_alloc_complex(int((ly/2+1) * lx, C_SIZE_T))
  call c_f_pointer(p_alloc_psiq, Pxq, [ly/2+1,lx])
  p_alloc_psiq=fftw_alloc_complex(int((ly/2+1) * lx, C_SIZE_T))
  call c_f_pointer(p_alloc_psiq, Pyq, [ly/2+1,lx])

  allocate(q2(ly/2+1,lx),qxy(ly/2+1,lx))
  allocate(lambda(3,3,ly/2+1,lx),chi(3,3,ly/2+1,lx))
  allocate(P2(ly,lx))

# ifdef open_MP
  ! use multi-threaded fftw with openMP
  nthreads=omp_get_max_threads()
  iret=fftw_init_threads()
  if(iret.eq.0) then
     write(*,*) 'iret=0: error during thread initialization'
     stop
  endif
  call fftw_plan_with_nthreads(nthreads)
  open(unit=2,file=file//'_para.dat',position='append')
  write(2,*)
  write(2,*) 'use multi-threaded fftw; nthreads=',nthreads
  close(2)
# endif
  
  yy=ran3(idum) ! initialize random numbers (after calling it idum=1 always, untill idum is reset as a negative value)
 
! for FFTW parameters
  scale2d_b=1.d0/(ly*lx) ! for backward FFT

! initialization for FFTW (out of place; create the plans before initializing the inputs)
  plan2d_f=fftw_plan_dft_r2c_2d(lx,ly,psi,psiq,FFTW_MEASURE) ! note the reverse order of dimensions
  plan2d_b=fftw_plan_dft_c2r_2d(lx,ly,psiq,psi,FFTW_MEASURE)

! set up the homogeneous initial conditions
     do i=1,lx
        do j=1,ly
           !psi(j,i)=noise0*(ran3(idum)-0.5d0) ! here psi=dpsi (full psi=psi0+dpsi)
           psi(j,i)=psi0*noise*(ran3(idum)-0.5d0) ! noise as percentage, here <psi>=<dpsi>=0
      
           if(abs(Px0).lt.1.d-10) then
              Px(j,i)=Px0+noise0*(ran3(idum)-0.5d0) ! if Px0 ~ 0
           else
              Px(j,i)=Px0*(1+noise*(ran3(idum)-0.5d0))
              !Px(j,i)=Px0+noise*(ran3(idum)-0.5d0)
           endif
           if(abs(Py0).lt.1.d-10) then
              Py(j,i)=Py0+noise0*(ran3(idum)-0.5d0) ! if Py0 ~ 0
           else
              Py(j,i)=Py0*(1+noise*(ran3(idum)-0.5d0))
              !Py(j,i)=Py0+noise*(ran3(idum)-0.5d0)
           endif
        enddo
     enddo

     call fftw_execute_dft_r2c(plan2d_f,psi,psiq) ! forward FFT
     psiq(1,1)=0 ! <dpsi>=0, <psi>=psi0
! backward FFT (note: the input array is overwritten for all c2r transforms)
     call fftw_execute_dft_c2r(plan2d_b,psiq,psi)
     psi=psi*scale2d_b
     call fftw_execute_dft_r2c(plan2d_f,Px,Pxq) ! forward FFT
     call fftw_execute_dft_r2c(plan2d_f,Py,Pyq) ! forward FFT
     Pxq(1,1)=Px0*lx*ly ! <Px>=Px0
     Pyq(1,1)=Py0*lx*ly ! <Py>=Py0
     call fftw_execute_dft_c2r(plan2d_b,Pxq,Px)
     Px=Px*scale2d_b
     call fftw_execute_dft_c2r(plan2d_b,Pyq,Py)
     Py=Py*scale2d_b
     
  call fftw_execute_dft_r2c(plan2d_f,psi,psiq) ! forward FFT
  call fftw_execute_dft_r2c(plan2d_f,Px,Pxq) ! forward FFT
  call fftw_execute_dft_r2c(plan2d_f,Py,Pyq) ! forward FFT

  time=time0

  if(iter0.le.n_i) then
     dtime=dti ! for initial transient
  else
     dtime=dt
  endif
  
  allocate(alpha_ij(3,3,ly/2+1,lx))
  allocate(sigma(3,ly/2+1,lx))
  allocate(isigma(ly/2+1,lx))

  isigma=0
  alpha_ij=0
  sigma=0
  aa=1
  do i=1,lx
     do j=1,ly/2+1
        q2(j,i)=qx2(i)+qy2(j)
        qxy(j,i)=qx(i)*qy(j)
        alpha_ij(1,1,j,i)=q2(j,i)*(B14*psi0*q2(j,i)-(1+B7*psi0)) ! here psi=dpsi (full psi=psi0+dpsi)
        alpha_ij(1,2,j,i)=cmplx(0,B1psi0*qx(i))
        alpha_ij(1,3,j,i)=cmplx(0,B1psi0*qy(j))
        alpha_ij(2,1,j,i)=4*alpha_ij(1,2,j,i)
        alpha_ij(3,1,j,i)=4*alpha_ij(1,3,j,i)
        alpha_ij(2,2,j,i)=-Dr*(1+2*B20*psi0)-(1+B1_2_4Dr*psi0*psi0)*q2(j,i) ! B1_2_4Dr=B1*B1/(4*Dr)
        alpha_ij(2,3,j,i)=0
        alpha_ij(3,2,j,i)=0
        alpha_ij(3,3,j,i)=alpha_ij(2,2,j,i)
        B1_B10q2=4*B1psi0*B1psi0*q2(j,i) ! =4*psi0^2*B1^2*q^2, i.e., =-(alpha12*alpha21+alpha13*alpha31)
! solutions of the characteristic cubic equation to get sigma1, sigma2, sigma3
        if(abs(B1_B10q2) < 1.e-20) then ! when B1_B10q2 = -(alpha12*alpha21+alpha13*alpha31) = 0
           sigma(1,j,i)=alpha_ij(1,1,j,i)
           sigma(2,j,i)=alpha_ij(2,2,j,i)
           sigma(3,j,i)=sigma(2,j,i)
           if(abs(sigma(1,j,i)-sigma(2,j,i)).lt.1.e-20) then
              isigma(j,i)=3 ! 3 equal roots: sigma1=sigma2=sigma3
           else
              isigma(j,i)=2 ! 2 equal roots: sigma1 \neq sigma2 = sigma3
           endif
        else
           sigma(1,j,i)=alpha_ij(2,2,j,i)
           Delta_=(alpha_ij(1,1,j,i)-alpha_ij(2,2,j,i))**2-4*B1_B10q2
           if(abs(Delta_) < 1.e-20) then ! Delta=0
              sigma(2,j,i)=(alpha_ij(1,1,j,i)+alpha_ij(2,2,j,i))/2
              sigma(3,j,i)=sigma(2,j,i)
              isigma(j,i)=2 ! 2 equal roots: sigma1 \neq sigma2 = sigma3
           else ! see Numerical Recipes p178 for solutions of aa*x^2+bb*x+cc=0
              isigma(j,i)=1 ! 3 distinct roots
              bb=-(alpha_ij(1,1,j,i)+alpha_ij(2,2,j,i)) ! (here aa=1)
              if(Delta_ > 0) then ! 3 distinct real roots
                 cc=alpha_ij(1,1,j,i)*alpha_ij(2,2,j,i)+B1_B10q2
                 if(bb < 0) then
                    qq=-(bb-sqrt(Delta_))/2
                    sigma(2,j,i)=qq/aa
                    sigma(3,j,i)=cc/qq
                 else ! if bb > 0
                    qq=-(bb+sqrt(Delta_))/2
                    sigma(3,j,i)=qq/aa
                    sigma(2,j,i)=cc/qq
                 endif
              else ! Delta < 0; 1 real and 2 conjugate complex roots
                 sigma(2,j,i)=(cmplx(-bb,sqrt(-Delta_)))/(2*aa)
                 sigma(3,j,i)=conjg(sigma(2,j,i))
              endif
           endif
        endif
     enddo
  enddo
  
  call coefficients_eqs(dtime) ! obtain the coefficients used in numerical solution

  if(iter0.gt.n_i+1) deallocate(alpha_ij,sigma,isigma)

  ave_Px=sum(Px)*scale2d_b
  ave_Py=sum(Py)*scale2d_b
  P2=sqrt(Px*Px+Py*Py)
  ave_P=sum(P2)*scale2d_b
  
  open(unit=2,file=file//'_out.dat',position='append')
  write(2,*)
  write(2,*) 't=',time
  write(2,*) 'dpsi=',maxval(psi(1:ly,1:lx)),minval(psi(1:ly,1:lx)) ! here psi=dpsi (full psi=psi0+dpsi)
  write(2,*) '  <dpsi>=',sum(psi)*scale2d_b
  write(2,*) 'Px=',maxval(Px(1:ly,1:lx)),minval(Px(1:ly,1:lx))
  write(2,*) '  <Px>=',ave_Px
  write(2,*) 'Py=',maxval(Py(1:ly,1:lx)),minval(Py(1:ly,1:lx))
  write(2,*) '  <Py>=',ave_Py
  write(2,*) '|P|=',maxval(P2(1:ly,1:lx)),minval(P2(1:ly,1:lx))
  write(2,*) '  <|P|>=',ave_P
  close(2)

  open(unit=5,file=file//'_aveP.dat')
  write(5,*) '# time,   ave_Px,   ave_Py,   ave_P'
  write(5,5) time,ave_Px,ave_Py,ave_P
5 format(f12.3,3(1x,es15.8))
  close(5)
  
! start running

  do iter=iter0,nend

    if(iter.le.n_i) then ! for initial transient only
       time=iter*dti
    else
       time=(iter-n_i)*dt+t_i
       if(iter.eq.n_i+1) then
          dtime=dt
          call coefficients_eqs(dtime) ! obtain the coefficients used in numerical solution
          deallocate(alpha_ij,sigma,isigma)
       endif
    endif
  
    call eqs(time)
    
    if((iter.eq.n_i).or.(iter.gt.n_i.and.mod(iter+n_0,raout).eq.0).or.(iter.eq.nend)) then 

      ave_Px=sum(Px)*scale2d_b
      ave_Py=sum(Py)*scale2d_b
      P2=sqrt(Px*Px+Py*Py)
      ave_P=sum(P2)*scale2d_b
      open(unit=2,file=file//'_out.dat',position='append')
      write(2,*)
      write(2,*) 't=',time
      write(2,*) 'dpsi=',maxval(psi(1:ly,1:lx)),minval(psi(1:ly,1:lx)) ! here psi=dpsi (full psi=psi0+dpsi)
      write(2,*) '  <dpsi>=',sum(psi)*scale2d_b
      write(2,*) 'Px=',maxval(Px(1:ly,1:lx)),minval(Px(1:ly,1:lx))
      write(2,*) '  <Px>=',ave_Px
      write(2,*) 'Py=',maxval(Py(1:ly,1:lx)),minval(Py(1:ly,1:lx))
      write(2,*) '  <Py>=',ave_Py
      write(2,*) '|P|=',maxval(P2(1:ly,1:lx)),minval(P2(1:ly,1:lx))
      write(2,*) '  <|P|>=',ave_P
      close(2)
    
      open(unit=5,file=file//'_aveP.dat',position='append')
      write(5,5) time,ave_Px,ave_Py,ave_P
      close(5)

# ifdef save_conf
    ! save conf files
    if((iter.gt.n_i.and.mod(iter+n_0,nout_conf).eq.0).or.(iter.eq.nend)) then
       open(unit=11,file=file//'_conf.dat',form='unformatted')
       write(11) psi,Px,Py
       close(11)
       open(unit=2,file=file//'_out.dat',position='append')
       write(2,*) 't=',time, '  iter=',iter,' for output of conf files'
       close(2)
    endif
# endif

   endif ! for raout

  enddo ! for iter=iter0,nend
  
  deallocate(q2,qxy,P2)
  deallocate(lambda,chi)
  call fftw_free(p_alloc_psi)
  call fftw_free(p_alloc_psiq)

  call fftw_destroy_plan(plan2d_f)
  call fftw_destroy_plan(plan2d_b)

end Program active_model_N_2D


! obtain the coefficients used in numerical solution (this subroutine kept the same for different set of eqs)
!subroutine coefficients_eqs(dtime,sigma,alpha_ij,isigma)
subroutine coefficients_eqs(dtime)
  use global_variables
  implicit none

  real(8) :: dtime
  !complex(8) :: sigma(3,ly/2+1,lx),alpha_ij(3,3,ly/2+1,lx)
  !integer :: isigma(ly/2+1,lx)
  complex(8) :: LambdaABC(3,3),sigma123,sigma_kl
  complex(8) :: sig_dt,exp_sig(3),cf_1(3),cf2_1(3),cf_2,cf2_2,cf_3,cf2_3
  integer :: i,j,k,l,ii,jj,kk

  lambda=0
  chi=0
  do i=1,lx
     do j=1,ly/2+1
        do ii=1,3
           sig_dt=sigma(ii,j,i)*dtime
           exp_sig(ii)=exp(sig_dt)
           if(abs(sigma(ii,j,i)).gt.1.d-5) then
              cf_1(ii)=(exp_sig(ii)-1)/sigma(ii,j,i)
              cf2_1(ii)=(exp_sig(ii)-(1+sig_dt))/(sigma(ii,j,i)*sig_dt)
           else ! when sigma --> 0
              cf_1(ii)=dtime*(1.d0+sig_dt*(1.d0+sig_dt/3.d0)/2.d0)
              cf2_1(ii)=0.5d0*dtime*(1.d0+sig_dt*(1.d0+sig_dt/4.d0)/3.d0)
           endif
        enddo
        
        if(isigma(j,i).eq.1) then ! 3 distinct roots: sigma1 \neq sigma2 \neq sigma3
           sigma123=(sigma(1,j,i)-sigma(2,j,i))*(sigma(1,j,i)-sigma(3,j,i))*(sigma(2,j,i)-sigma(3,j,i))
           do ii=1,3 ! corresponding to sigma1,sigma2,sigma3
              k=ii+1
              if(k.gt.3) k=k-3
              l=ii+2
              if(l.gt.3) l=l-3
              sigma_kl=sigma(k,j,i)-sigma(l,j,i)

              ! coefficients to get nAq(t+dt)
              LambdaABC(1,ii)=((sigma(ii,j,i)-alpha_ij(2,2,j,i))*(sigma(ii,j,i)-alpha_ij(3,3,j,i)) &
                              -alpha_ij(2,3,j,i)*alpha_ij(3,2,j,i))*sigma_kl
              LambdaABC(2,ii)=(alpha_ij(1,2,j,i)*(sigma(ii,j,i)-alpha_ij(3,3,j,i)) &
                              +alpha_ij(1,3,j,i)*alpha_ij(3,2,j,i))*sigma_kl
              LambdaABC(3,ii)=(alpha_ij(1,3,j,i)*(sigma(ii,j,i)-alpha_ij(2,2,j,i)) &
                              +alpha_ij(1,2,j,i)*alpha_ij(2,3,j,i))*sigma_kl
              do jj=1,3
                 lambda(1,jj,j,i)=lambda(1,jj,j,i)+LambdaABC(jj,ii)*exp_sig(ii)
                 chi(1,jj,j,i)=chi(1,jj,j,i)+LambdaABC(jj,ii)*cf_1(ii)
              enddo

              ! coefficients to get nBq(t+dt)
              LambdaABC(1,ii)=(alpha_ij(2,1,j,i)*(sigma(ii,j,i)-alpha_ij(3,3,j,i)) &
                              +alpha_ij(2,3,j,i)*alpha_ij(3,1,j,i))*sigma_kl
              LambdaABC(2,ii)=((sigma(ii,j,i)-alpha_ij(1,1,j,i))*(sigma(ii,j,i)-alpha_ij(3,3,j,i)) &
                              -alpha_ij(1,3,j,i)*alpha_ij(3,1,j,i))*sigma_kl
              LambdaABC(3,ii)=(alpha_ij(2,3,j,i)*(sigma(ii,j,i)-alpha_ij(1,1,j,i)) &
                              +alpha_ij(1,3,j,i)*alpha_ij(2,1,j,i))*sigma_kl
              do jj=1,3
                 lambda(2,jj,j,i)=lambda(2,jj,j,i)+LambdaABC(jj,ii)*exp_sig(ii)
                 chi(2,jj,j,i)=chi(2,jj,j,i)+LambdaABC(jj,ii)*cf_1(ii)
              enddo
              
              ! coefficients to get nCq(t+dt)
              LambdaABC(1,ii)=(alpha_ij(3,1,j,i)*(sigma(ii,j,i)-alpha_ij(2,2,j,i)) &
                              +alpha_ij(2,1,j,i)*alpha_ij(3,2,j,i))*sigma_kl
              LambdaABC(2,ii)=(alpha_ij(3,2,j,i)*(sigma(ii,j,i)-alpha_ij(1,1,j,i)) &
                              +alpha_ij(1,2,j,i)*alpha_ij(3,1,j,i))*sigma_kl
              LambdaABC(3,ii)=((sigma(ii,j,i)-alpha_ij(1,1,j,i))*(sigma(ii,j,i)-alpha_ij(2,2,j,i)) &
                              -alpha_ij(1,2,j,i)*alpha_ij(2,1,j,i))*sigma_kl
              do jj=1,3
                 lambda(3,jj,j,i)=lambda(3,jj,j,i)+LambdaABC(jj,ii)*exp_sig(ii)
                 chi(3,jj,j,i)=chi(3,jj,j,i)+LambdaABC(jj,ii)*cf_1(ii)
              enddo
           enddo
           
           do kk=1,3
              do jj=1,3
                 lambda(kk,jj,j,i)=lambda(kk,jj,j,i)/sigma123
                 chi(kk,jj,j,i)=chi(kk,jj,j,i)/sigma123
              enddo
           enddo
        endif
        
        if(isigma(j,i).eq.2) then ! 2 equal roots: sigma1 \neq sigma2 = sigma3
           sig_dt=sigma(2,j,i)*dtime
           if(abs(sigma(2,j,i)).gt.1.d-5) then
              cf_2=(1+(sig_dt-1)*exp_sig(2))/(sigma(2,j,i)*sigma(2,j,i))
              cf2_2=(2+sig_dt+(sig_dt-2)*exp_sig(2))/(sigma(2,j,i)*sigma(2,j,i)*sig_dt)
           else ! when sigma --> 0
              cf_2=dtime*dtime*(0.5d0+sig_dt/3+sig_dt*sig_dt/8)
              cf2_2=dtime*dtime*(1.d0/6.d0+sig_dt/12+sig_dt*sig_dt/40)
           endif

           sigma123=sigma(1,j,i)-sigma(2,j,i)
           
           do kk=1,3 ! coefficients to get nAq(t+dt),nBq(t+dt),nCq(t+dt)
           do ii=1,2 ! corresponding to sigma1,sigma2=sigma3
              if(kk.eq.1) then ! for nAq(t+dt)
                 LambdaABC(1,ii)=(sigma(ii,j,i)-alpha_ij(2,2,j,i))*(sigma(ii,j,i)-alpha_ij(3,3,j,i)) &
                                 -alpha_ij(2,3,j,i)*alpha_ij(3,2,j,i)
                 LambdaABC(2,ii)=alpha_ij(1,2,j,i)*(sigma(ii,j,i)-alpha_ij(3,3,j,i)) &
                                 +alpha_ij(1,3,j,i)*alpha_ij(3,2,j,i)
                 LambdaABC(3,ii)=alpha_ij(1,3,j,i)*(sigma(ii,j,i)-alpha_ij(2,2,j,i)) &
                                 +alpha_ij(1,2,j,i)*alpha_ij(2,3,j,i)
              endif
              if(kk.eq.2) then ! for nBq(t+dt)
                 LambdaABC(1,ii)=alpha_ij(2,1,j,i)*(sigma(ii,j,i)-alpha_ij(3,3,j,i)) &
                                 +alpha_ij(2,3,j,i)*alpha_ij(3,1,j,i)
                 LambdaABC(2,ii)=(sigma(ii,j,i)-alpha_ij(1,1,j,i))*(sigma(ii,j,i)-alpha_ij(3,3,j,i)) &
                                 -alpha_ij(1,3,j,i)*alpha_ij(3,1,j,i)
                 LambdaABC(3,ii)=alpha_ij(2,3,j,i)*(sigma(ii,j,i)-alpha_ij(1,1,j,i)) &
                                 +alpha_ij(1,3,j,i)*alpha_ij(2,1,j,i)
              endif
              if(kk.eq.3) then ! for nCq(t+dt)
                 LambdaABC(1,ii)=alpha_ij(3,1,j,i)*(sigma(ii,j,i)-alpha_ij(2,2,j,i)) &
                                 +alpha_ij(2,1,j,i)*alpha_ij(3,2,j,i)
                 LambdaABC(2,ii)=alpha_ij(3,2,j,i)*(sigma(ii,j,i)-alpha_ij(1,1,j,i)) &
                                 +alpha_ij(1,2,j,i)*alpha_ij(3,1,j,i)
                 LambdaABC(3,ii)=(sigma(ii,j,i)-alpha_ij(1,1,j,i))*(sigma(ii,j,i)-alpha_ij(2,2,j,i)) &
                                 -alpha_ij(1,2,j,i)*alpha_ij(2,1,j,i)
              endif
           enddo
           do jj=1,3
              if(jj.eq.kk) then
                 sigma_kl=-(sigma(1,j,i)-alpha_ij(kk,kk,j,i))
              else
                 sigma_kl=alpha_ij(kk,jj,j,i)
              endif
              lambda(kk,jj,j,i)=LambdaABC(jj,1)*exp_sig(1)-(sigma123*(sigma_kl+LambdaABC(jj,2)*dtime) &
                                                           +LambdaABC(jj,2))*exp_sig(2)
              chi(kk,jj,j,i)=LambdaABC(jj,1)*cf_1(1)-(LambdaABC(jj,2)+sigma_kl*sigma123)*cf_1(2) &
                            -LambdaABC(jj,2)*sigma123*cf_2
           enddo
           enddo ! do kk=1,3

           sigma123=sigma123**2
           do kk=1,3
              do jj=1,3
                 lambda(kk,jj,j,i)=lambda(kk,jj,j,i)/sigma123
                 chi(kk,jj,j,i)=chi(kk,jj,j,i)/sigma123
              enddo
           enddo           
        endif
        
        if(isigma(j,i).eq.3) then ! 3 equal roots: sigma1 = sigma2 = sigma3
           sig_dt=sigma(1,j,i)*dtime
           if(abs(sigma(1,j,i)).gt.1.d-5) then
              cf_2=(1+(sig_dt-1)*exp_sig(1))/(sigma(1,j,i)*sigma(1,j,i))
              cf2_2=(2+sig_dt+(sig_dt-2)*exp_sig(1))/(sigma(1,j,i)*sigma(1,j,i)*sig_dt)
              cf_3=(((sig_dt-1)**2+1)*exp_sig(1)-2)/(sigma(1,j,i)**3)
              cf2_3=(((sig_dt-2)**2+2)*exp_sig(1)-2*(sig_dt+3))/(sigma(1,j,i)**3*sig_dt)
           else ! when sigma --> 0
              cf_2=dtime*dtime*(0.5d0+sig_dt/3+sig_dt*sig_dt/8)
              cf2_2=dtime*dtime*(1.d0/6.d0+sig_dt/12+sig_dt*sig_dt/40)
              cf_3=dtime*dtime*dtime*(1.d0/3+sig_dt/4+sig_dt*sig_dt/10)
              cf2_3=dtime*dtime*dtime*(1.d0/12+sig_dt/20+sig_dt*sig_dt/60)
           endif

           ii=1
           do kk=1,3 ! coefficients to get nAq(t+dt),nBq(t+dt),nCq(t+dt)
              if(kk.eq.1) then ! for nAq(t+dt)
                 LambdaABC(1,ii)=(sigma(ii,j,i)-alpha_ij(2,2,j,i))*(sigma(ii,j,i)-alpha_ij(3,3,j,i)) &
                                 -alpha_ij(2,3,j,i)*alpha_ij(3,2,j,i)
                 LambdaABC(2,ii)=alpha_ij(1,2,j,i)*(sigma(ii,j,i)-alpha_ij(3,3,j,i)) &
                                 +alpha_ij(1,3,j,i)*alpha_ij(3,2,j,i)
                 LambdaABC(3,ii)=alpha_ij(1,3,j,i)*(sigma(ii,j,i)-alpha_ij(2,2,j,i)) &
                                 +alpha_ij(1,2,j,i)*alpha_ij(2,3,j,i)
              endif
              if(kk.eq.2) then ! for nBq(t+dt)
                 LambdaABC(1,ii)=alpha_ij(2,1,j,i)*(sigma(ii,j,i)-alpha_ij(3,3,j,i)) &
                                 +alpha_ij(2,3,j,i)*alpha_ij(3,1,j,i)
                 LambdaABC(2,ii)=(sigma(ii,j,i)-alpha_ij(1,1,j,i))*(sigma(ii,j,i)-alpha_ij(3,3,j,i)) &
                                 -alpha_ij(1,3,j,i)*alpha_ij(3,1,j,i)
                 LambdaABC(3,ii)=alpha_ij(2,3,j,i)*(sigma(ii,j,i)-alpha_ij(1,1,j,i)) &
                                 +alpha_ij(1,3,j,i)*alpha_ij(2,1,j,i)
              endif
              if(kk.eq.3) then ! for nCq(t+dt)
                 LambdaABC(1,ii)=alpha_ij(3,1,j,i)*(sigma(ii,j,i)-alpha_ij(2,2,j,i)) &
                                 +alpha_ij(2,1,j,i)*alpha_ij(3,2,j,i)
                 LambdaABC(2,ii)=alpha_ij(3,2,j,i)*(sigma(ii,j,i)-alpha_ij(1,1,j,i)) &
                                 +alpha_ij(1,2,j,i)*alpha_ij(3,1,j,i)
                 LambdaABC(3,ii)=(sigma(ii,j,i)-alpha_ij(1,1,j,i))*(sigma(ii,j,i)-alpha_ij(2,2,j,i)) &
                                 -alpha_ij(1,2,j,i)*alpha_ij(2,1,j,i)
              endif
              do jj=1,3
                 if(jj.eq.kk) then
                    sigma_kl=-(sigma(1,j,i)-alpha_ij(kk,kk,j,i))
                    sig_dt=1
                 else
                    sigma_kl=alpha_ij(kk,jj,j,i)
                    sig_dt=0
                 endif
                 lambda(kk,jj,j,i)=(sig_dt+(sigma_kl+LambdaABC(jj,1)*dtime/2)*dtime)*exp_sig(1)
                 chi(kk,jj,j,i)=sig_dt*cf_1(1)+sigma_kl*cf_2+LambdaABC(jj,1)*cf_3/2
              enddo
           enddo ! do kk=1,3
        endif
     enddo
  enddo
  
  open(unit=2,file=file//'_out.dat',position='append')
  write(2,*)
  write(2,*) 'dtime=',dtime, '  from lambda, chi formulae'
  write(2,*) 'lambda11,12,13 at q=0:',lambda(1,1,1,1),lambda(1,2,1,1),lambda(1,3,1,1)
  write(2,*) 'lambda21,22,23 at q=0:',lambda(2,1,1,1),lambda(2,2,1,1),lambda(2,3,1,1)
  write(2,*) 'lambda31,32,33 at q=0:',lambda(3,1,1,1),lambda(3,2,1,1),lambda(3,3,1,1)
  write(2,*) '|lambda|11,12,13 max:'
  write(2,*) maxval(abs(lambda(1,1,:,:))),maxval(abs(lambda(1,2,:,:))),maxval(abs(lambda(1,3,:,:)))
  write(2,*) '|chi|11,12,13 max:'
  write(2,*) maxval(abs(chi(1,1,:,:))),maxval(abs(chi(1,2,:,:))),maxval(abs(chi(1,3,:,:)))
  write(2,*) '|lambda|21,22,23 max:'
  write(2,*) maxval(abs(lambda(2,1,:,:))),maxval(abs(lambda(2,2,:,:))),maxval(abs(lambda(2,3,:,:)))
  write(2,*) '|chi|21,22,23 max:'
  write(2,*) maxval(abs(chi(2,1,:,:))),maxval(abs(chi(2,2,:,:))),maxval(abs(chi(2,3,:,:)))
  write(2,*) '|lambda|31,32,33 max:'
  write(2,*) maxval(abs(lambda(3,1,:,:))),maxval(abs(lambda(3,2,:,:))),maxval(abs(lambda(3,3,:,:)))
  write(2,*) '|chi|31,32,33 max:'
  write(2,*) maxval(abs(chi(3,1,:,:))),maxval(abs(chi(3,2,:,:))),maxval(abs(chi(3,3,:,:)))
  close(2)

  return
end subroutine coefficients_eqs


! main routine that solves the mean field eqs (this subroutine kept the same for different set of eqs)
subroutine eqs(time)
  use global_variables
  implicit none

  real(8) :: time
  integer :: i,j,ly_2

  complex(C_DOUBLE_COMPLEX), pointer :: psi0q(:,:),Px0q(:,:),Py0q(:,:),psiq_tmp(:,:)
  complex(C_DOUBLE_COMPLEX), pointer :: Nonlinpsi_q(:,:),NonlinPx_q(:,:),NonlinPy_q(:,:)
! for allocating aligned memory in fftw (need different pointers for different arrays
! so that the memory can be released at the end of each call (otherwise the used memory will accumulate))
  type(C_PTR) :: p_alloc_psi0q,p_alloc_Nonlinpsi_q,p_alloc_psiq_tmp
  type(C_PTR) :: p_alloc_Px0q,p_alloc_NonlinPx_q
  type(C_PTR) :: p_alloc_Py0q,p_alloc_NonlinPy_q
    
  ly_2=ly/2

  p_alloc_psi0q=fftw_alloc_complex(int((ly_2+1) * lx, C_SIZE_T))
  call c_f_pointer(p_alloc_psi0q, psi0q, [ly_2+1,lx])
  p_alloc_Px0q=fftw_alloc_complex(int((ly_2+1) * lx, C_SIZE_T))
  call c_f_pointer(p_alloc_Px0q, Px0q, [ly_2+1,lx])
  p_alloc_Py0q=fftw_alloc_complex(int((ly_2+1) * lx, C_SIZE_T))
  call c_f_pointer(p_alloc_Py0q, Py0q, [ly_2+1,lx])
  p_alloc_psiq_tmp=fftw_alloc_complex(int((ly_2+1) * lx, C_SIZE_T))
  call c_f_pointer(p_alloc_psiq_tmp, psiq_tmp, [ly_2+1,lx])
  p_alloc_Nonlinpsi_q=fftw_alloc_complex(int((ly_2+1) * lx, C_SIZE_T))
  call c_f_pointer(p_alloc_Nonlinpsi_q, Nonlinpsi_q, [ly_2+1,lx])
  p_alloc_NonlinPx_q=fftw_alloc_complex(int((ly_2+1) * lx, C_SIZE_T))
  call c_f_pointer(p_alloc_NonlinPx_q, NonlinPx_q, [ly_2+1,lx])
  p_alloc_NonlinPy_q=fftw_alloc_complex(int((ly_2+1) * lx, C_SIZE_T))
  call c_f_pointer(p_alloc_NonlinPy_q, NonlinPy_q, [ly_2+1,lx])
  
!  calculate the nonlinear terms
  call nonlin_q(Nonlinpsi_q,NonlinPx_q,NonlinPy_q)

! First step: predictor
! predictor values of psiq, Pxq, Pyq

!$OMP PARALLEL PRIVATE(i,j)
!$OMP DO
!(Note: DO loop index variable is by default PRIVATE, all other variables are by default SHARED unless declared)
  do i=1,lx
     do j=1,ly_2+1
        psi0q(j,i)=lambda(1,1,j,i)*psiq(j,i)+lambda(1,2,j,i)*Pxq(j,i)+lambda(1,3,j,i)*Pyq(j,i) &
                 +chi(1,1,j,i)*Nonlinpsi_q(j,i)+chi(1,2,j,i)*NonlinPx_q(j,i)+chi(1,3,j,i)*NonlinPy_q(j,i)
        Px0q(j,i)=lambda(2,1,j,i)*psiq(j,i)+lambda(2,2,j,i)*Pxq(j,i)+lambda(2,3,j,i)*Pyq(j,i) &
                 +chi(2,1,j,i)*Nonlinpsi_q(j,i)+chi(2,2,j,i)*NonlinPx_q(j,i)+chi(2,3,j,i)*NonlinPy_q(j,i)
        Py0q(j,i)=lambda(3,1,j,i)*psiq(j,i)+lambda(3,2,j,i)*Pxq(j,i)+lambda(3,3,j,i)*Pyq(j,i) &
                 +chi(3,1,j,i)*Nonlinpsi_q(j,i)+chi(3,2,j,i)*NonlinPx_q(j,i)+chi(3,3,j,i)*NonlinPy_q(j,i)
     enddo
  enddo
!$OMP END DO
!$OMP END PARALLEL

!$OMP PARALLEL DO
  do i=lx/2+2,lx
     psi0q(1,i)=conjg(psi0q(1,lx-i+2)) ! real data treatment
     Px0q(1,i)=conjg(Px0q(1,lx-i+2))
     Py0q(1,i)=conjg(Py0q(1,lx-i+2))
  enddo
!$OMP END PARALLEL DO

  psiq_tmp=psi0q
! backward FFT (note: the input array is overwritten)
  call fftw_execute_dft_c2r(plan2d_b,psiq_tmp,psi)
  psi=psi*scale2d_b
  psiq=psi0q
  
  psiq_tmp=Px0q
! backward FFT (note: the input array is overwritten)
  call fftw_execute_dft_c2r(plan2d_b,psiq_tmp,Px)
  Px=Px*scale2d_b
  Pxq=Px0q

  psiq_tmp=Py0q
! backward FFT (note: the input array is overwritten)
  call fftw_execute_dft_c2r(plan2d_b,psiq_tmp,Py)
  Py=Py*scale2d_b
  Pyq=Py0q
  
  if((maxval(psi(1:ly,1:lx)).gt.1.d5).or.(maxval(Px(1:ly,1:lx)).gt.1.d5).or.(maxval(Py(1:ly,1:lx)).gt.1.d5)) then
     open(unit=2,file=file//'_out.dat',position='append')
     write(2,*) 'psi,Px,Py diverge at t=',time
     write(2,*) 'dpsi=',maxval(psi(1:ly,1:lx)),minval(psi(1:ly,1:lx))
     write(2,*) 'Px=',maxval(Px(1:ly,1:lx)),minval(Px(1:ly,1:lx))
     write(2,*) 'Py=',maxval(Py(1:ly,1:lx)),minval(Py(1:ly,1:lx))
     close(2)
     stop
  endif

  call fftw_free(p_alloc_psi0q)
  call fftw_free(p_alloc_psiq_tmp)
  call fftw_free(p_alloc_Nonlinpsi_q)
  call fftw_free(p_alloc_Px0q)
  call fftw_free(p_alloc_NonlinPx_q)
  call fftw_free(p_alloc_Py0q)
  call fftw_free(p_alloc_NonlinPy_q)
  return
end subroutine eqs


!  calculate the nonlinear terms
subroutine nonlin_q(Nonlin_psiq,Nonlin_Pxq,Nonlin_Pyq)
  use global_variables
  implicit none

  integer :: i,j
  complex(C_DOUBLE_COMPLEX) :: Nonlin_psiq(ly/2+1,lx),Nonlin_Pxq(ly/2+1,lx),Nonlin_Pyq(ly/2+1,lx)
  real(C_DOUBLE), pointer :: Nonlin_1(:,:),Nonlin_2(:,:),Nonlin_3(:,:),dx2psi(:,:),dy2psi(:,:)
  real(C_DOUBLE), pointer :: dxPx_dyPy(:,:),dxPydyPx(:,:),Nonlin_psi0_1(:,:),Nonlin_psi0_2(:,:)
  complex(C_DOUBLE_COMPLEX), pointer :: psiPq(:,:),psiP2q(:,:),psiP_q(:,:)
  type(C_PTR) :: p_alloc_Nonlin_1,p_alloc_Nonlin_2,p_alloc_Nonlin_3 ! for allocating aligned memory in fftw
  type(C_PTR) :: p_alloc_dx2psi,p_alloc_dy2psi,p_alloc_psiPq,p_alloc_psiP2q,p_alloc_psiP_q
  type(C_PTR) :: p_alloc_dxPx_dyPy,p_alloc_dxPydyPx,p_alloc_Nonlin_psi0_1,p_alloc_Nonlin_psi0_2
  double precision, allocatable :: PxPy(:,:),Px2_Py2(:,:),Nonlin_Px_psi0(:,:),Nonlin_Py_psi0(:,:)
  double precision, allocatable :: Px2(:,:),Py2(:,:)

  p_alloc_Nonlin_1=fftw_alloc_real(int(ly * lx, C_SIZE_T))
  call c_f_pointer(p_alloc_Nonlin_1, Nonlin_1, [ly,lx])
  p_alloc_Nonlin_2=fftw_alloc_real(int(ly * lx, C_SIZE_T))
  call c_f_pointer(p_alloc_Nonlin_2, Nonlin_2, [ly,lx])
  p_alloc_Nonlin_3=fftw_alloc_real(int(ly * lx, C_SIZE_T))
  call c_f_pointer(p_alloc_Nonlin_3, Nonlin_3, [ly,lx])
  p_alloc_dx2psi=fftw_alloc_real(int(ly * lx, C_SIZE_T))
  call c_f_pointer(p_alloc_dx2psi, dx2psi, [ly,lx])
  p_alloc_dy2psi=fftw_alloc_real(int(ly * lx, C_SIZE_T))
  call c_f_pointer(p_alloc_dy2psi, dy2psi, [ly,lx])
  p_alloc_dxPx_dyPy=fftw_alloc_real(int(ly * lx, C_SIZE_T))
  call c_f_pointer(p_alloc_dxPx_dyPy, dxPx_dyPy, [ly,lx])
  p_alloc_dxPydyPx=fftw_alloc_real(int(ly * lx, C_SIZE_T))
  call c_f_pointer(p_alloc_dxPydyPx, dxPydyPx, [ly,lx])
  p_alloc_psiPq=fftw_alloc_complex(int((ly/2+1) * lx, C_SIZE_T))
  call c_f_pointer(p_alloc_psiPq, psiPq, [ly/2+1,lx])
  p_alloc_psiP2q=fftw_alloc_complex(int((ly/2+1) * lx, C_SIZE_T))
  call c_f_pointer(p_alloc_psiP2q, psiP2q, [ly/2+1,lx])
  p_alloc_psiP_q=fftw_alloc_complex(int((ly/2+1) * lx, C_SIZE_T))
  call c_f_pointer(p_alloc_psiP_q, psiP_q, [ly/2+1,lx])

  p_alloc_Nonlin_psi0_1=fftw_alloc_real(int(ly * lx, C_SIZE_T))
  call c_f_pointer(p_alloc_Nonlin_psi0_1, Nonlin_psi0_1, [ly,lx])
  p_alloc_Nonlin_psi0_2=fftw_alloc_real(int(ly * lx, C_SIZE_T))
  call c_f_pointer(p_alloc_Nonlin_psi0_2, Nonlin_psi0_2, [ly,lx])

  allocate(PxPy(ly,lx),Px2_Py2(ly,lx),Nonlin_Px_psi0(ly,lx),Nonlin_Py_psi0(ly,lx))
  allocate(Px2(ly,lx),Py2(ly,lx))

!$OMP PARALLEL DO
  do i=1,lx
     do j=1,ly/2+1
        Nonlin_Pxq(j,i)=CMPLX(0.d0,-qx(i)*q2(j,i))*psiq(j,i) ! \partial_x \nabla^2 psi
        Nonlin_Pyq(j,i)=CMPLX(0.d0,-qy(j)*q2(j,i))*psiq(j,i) ! \partial_y \nabla^2 psi
     enddo
  enddo
!$OMP END PARALLEL DO
! backward FFT (note: the input array is overwritten)
  call fftw_execute_dft_c2r(plan2d_b,Nonlin_Pxq,Nonlin_1) ! = \partial_x \nabla^2 psi
  call fftw_execute_dft_c2r(plan2d_b,Nonlin_Pyq,Nonlin_2) ! = \partial_y \nabla^2 psi
  Nonlin_1=Nonlin_1*scale2d_b
  Nonlin_2=Nonlin_2*scale2d_b
  
  Nonlin_1=psi*(B14*Nonlin_1+B1*Px)
  Nonlin_2=psi*(B14*Nonlin_2+B1*Py)
  
!$OMP PARALLEL DO
  do i=1,lx
     do j=1,ly/2+1
        Nonlin_Pxq(j,i)=-qx2(i)*psiq(j,i) ! qx2=qx*qx
        Nonlin_Pyq(j,i)=-qy2(j)*psiq(j,i) ! qy2=qy*qy
        Nonlin_psiq(j,i)=-qxy(j,i)*psiq(j,i) ! qxy=qx*qy
     enddo
  enddo
!$OMP END PARALLEL DO
! backward FFT (note: the input array is overwritten)
  call fftw_execute_dft_c2r(plan2d_b,Nonlin_Pxq,dx2psi) ! = \partial_x^2 psi
  dx2psi=dx2psi*scale2d_b
  call fftw_execute_dft_c2r(plan2d_b,Nonlin_Pyq,dy2psi) ! = \partial_y^2 psi
  dy2psi=dy2psi*scale2d_b
  call fftw_execute_dft_c2r(plan2d_b,Nonlin_psiq,Nonlin_3) ! = \partial_x\partial_y psi
  Nonlin_3=Nonlin_3*scale2d_b
  Nonlin_1=Nonlin_1+B12*(Px*(3*dx2psi+dy2psi)+2*Py*Nonlin_3)
  Nonlin_2=Nonlin_2+B12*(Py*(dx2psi+3*dy2psi)+2*Px*Nonlin_3)

!$OMP PARALLEL DO
  do i=1,lx
     do j=1,ly/2+1
        Nonlin_Pxq(j,i)=CMPLX(0.d0,qx(i))*psiq(j,i) ! \partial_x psi
        Nonlin_Pyq(j,i)=CMPLX(0.d0,qy(j))*psiq(j,i) ! \partial_y psi
     enddo
  enddo
!$OMP END PARALLEL DO
! backward FFT (note: the input array is overwritten)
  call fftw_execute_dft_c2r(plan2d_b,Nonlin_Pxq,dx2psi) ! = \partial_x psi
  dx2psi=dx2psi*scale2d_b
  call fftw_execute_dft_c2r(plan2d_b,Nonlin_Pyq,dy2psi) ! = \partial_y psi
  dy2psi=dy2psi*scale2d_b
  
!$OMP PARALLEL DO
  do i=1,lx
     do j=1,ly/2+1
        Nonlin_Pxq(j,i)=CMPLX(0.d0,qx(i))*Pxq(j,i) ! \partial_x Px
        Nonlin_Pyq(j,i)=CMPLX(0.d0,qy(j))*Pyq(j,i) ! \partial_y Py
     enddo
  enddo
!$OMP END PARALLEL DO
! backward FFT (note: the input array is overwritten)
  call fftw_execute_dft_c2r(plan2d_b,Nonlin_Pxq,dxPx_dyPy) ! = \partial_x Px
  dxPx_dyPy=dxPx_dyPy*scale2d_b
  call fftw_execute_dft_c2r(plan2d_b,Nonlin_Pyq,Nonlin_3) ! = \partial_y Py
  Nonlin_3=Nonlin_3*scale2d_b
  Nonlin_Px_psi0=-Px*(3*dxPx_dyPy+5*Nonlin_3) ! = -Px*(3\partial_x Px + 5\partial_y Py)
  Nonlin_Py_psi0=-Py*(5*dxPx_dyPy+3*Nonlin_3) ! = -Py*(5\partial_x Px + 3\partial_y Py)
  
  dxPx_dyPy=dxPx_dyPy-Nonlin_3 ! = \partial_x Px - \partial_y Py
  
!$OMP PARALLEL DO
  do i=1,lx
     do j=1,ly/2+1
        Nonlin_Pxq(j,i)=CMPLX(0.d0,qx(i))*Pyq(j,i) ! \partial_x Py
        Nonlin_Pyq(j,i)=CMPLX(0.d0,qy(j))*Pxq(j,i) ! \partial_y Px
     enddo
  enddo
!$OMP END PARALLEL DO
! backward FFT (note: the input array is overwritten)
  call fftw_execute_dft_c2r(plan2d_b,Nonlin_Pxq,dxPydyPx) ! = \partial_x Py
  dxPydyPx=dxPydyPx*scale2d_b
  call fftw_execute_dft_c2r(plan2d_b,Nonlin_Pyq,Nonlin_3) ! = \partial_y Px
  Nonlin_3=Nonlin_3*scale2d_b
  Nonlin_Px_psi0=Nonlin_Px_psi0+Py*(5*dxPydyPx-3*Nonlin_3) ! to be used in the psi0 term of Nonlin_Pxq
  Nonlin_Py_psi0=Nonlin_Py_psi0+Px*(5*Nonlin_3-3*dxPydyPx) ! to be used in the psi0 term of Nonlin_Pyq

  dxPydyPx=dxPydyPx+Nonlin_3 ! = \partial_x Py + \partial_y Px

  Nonlin_psi0_1=dx2psi*dxPx_dyPy+dy2psi*dxPydyPx ! used in the psi0 term inside B1_2_4Dr*psi0 of Nonlin_Pxq
  Nonlin_psi0_2=dx2psi*dxPydyPx-dy2psi*dxPx_dyPy ! used in the psi0 term inside B1_2_4Dr*psi0 of Nonlin_Pyq

  dxPx_dyPy=psi*dxPx_dyPy ! = psi*(\partial_x Px - \partial_y Py)
  dxPydyPx=psi*dxPydyPx ! = psi*(\partial_x Py + \partial_y Px)

  PxPy=Px*Py
  Px2=Px*Px
  Py2=Py*Py
  Px2_Py2=Px2-Py2 ! =Px*Px-Py*Py
  
  Nonlin_3=dx2psi*(dxPx_dyPy+Px*dx2psi)+dy2psi*(dxPydyPx+Px*dy2psi) ! used in B1_2_4Dr term of Nonlin_Pxq
  call fftw_execute_dft_r2c(plan2d_f,Nonlin_1,Nonlin_Pxq) ! forward FFT (used for Nonlin_psiq below)
  Nonlin_1=dx2psi*(dxPydyPx+Py*dx2psi)+dy2psi*(Py*dy2psi-dxPx_dyPy) ! used in B1_2_4Dr term of Nonlin_Pyq
  call fftw_execute_dft_r2c(plan2d_f,Nonlin_2,Nonlin_Pyq) ! forward FFT (used for Nonlin_psiq below)
  
!$OMP PARALLEL DO
  do i=1,lx
     do j=1,ly/2+1
        Nonlin_psiq(j,i)=CMPLX(0.d0,qx(i))*Nonlin_Pxq(j,i)+CMPLX(0.d0,qy(j))*Nonlin_Pyq(j,i)
     enddo
  enddo
!$OMP END PARALLEL DO
  Nonlin_2=psi*psi
  call fftw_execute_dft_r2c(plan2d_f,Nonlin_2,psiPq) ! forward FFT
  Nonlin_psiq=Nonlin_psiq-B7_2*q2*psiPq ! Nonlin_psiq calculated

  Nonlin_2=B1_2*psi*(psi-B20_4*Px2_Py2) ! B1_2=2*B1, B20_4=B20/4, Px2_Py2=Px^2-Py^2
  call fftw_execute_dft_r2c(plan2d_f,Nonlin_2,Nonlin_Pxq) ! forward FFT
  Nonlin_2=B1_2*psi*(psi+B20_4*Px2_Py2)
  call fftw_execute_dft_r2c(plan2d_f,Nonlin_2,Nonlin_Pyq) ! forward FFT
  Nonlin_2=B1B20*psi*PxPy ! B1B20=B1*B20
  call fftw_execute_dft_r2c(plan2d_f,Nonlin_2,psiPq) ! forward FFT; psiPq=B1*B20*(psi*Px*Py)_q
  
  PxPy=psi*Px
  call fftw_execute_dft_r2c(plan2d_f,PxPy,psiP_q) ! forward FFT
  psiP_q=-q2*psiP_q ! =-q^2*(psi*Px)_q
  psiP2q=psiP_q
  ! backward FFT (note: the input array is overwritten)
  call fftw_execute_dft_c2r(plan2d_b,psiP2q,PxPy)
  PxPy=PxPy*scale2d_b  ! = \nabla^2(psi*Px)
  
  Px2_Py2=Px2+Py2 ! = Px^2+Py^2, Px2=Px*Px, Py2=Py*Py calculated above
  Nonlin_2=Dr2_B20*(psi+B20_4*Px2_Py2) ! =2*Dr*B20*(psi+B20/4*P^2)
  
  PxPy=B1B20_4*(Px2_Py2*dx2psi+Py*dxPydyPx+Px*dxPx_dyPy)-Px*Nonlin_2+B1_2_4Dr*(Nonlin_3+PxPy*psi+psi0*Nonlin_psi0_1)

  psiP2q=-q2*Pxq
  ! backward FFT (note: the input array is overwritten)
  call fftw_execute_dft_c2r(plan2d_b,psiP2q,Nonlin_psi0_1) ! = \nabla^2 Px
  Nonlin_psi0_1=Nonlin_psi0_1*scale2d_b

  PxPy=PxPy+psi0*((B1B20_4)*Nonlin_Px_psi0+B1_2_4Dr*psi*Nonlin_psi0_1)
      ! Nonlin_Px_psi0=Py*(5\partial_x Py - 3\partial_y Px) - Px*(3\partial_x Px + 5\partial_y Py)
  call fftw_execute_dft_r2c(plan2d_f,PxPy,psiP2q) ! forward FFT
!$OMP PARALLEL DO
  do i=1,lx
     do j=1,ly/2+1
        Nonlin_Pxq(j,i)=CMPLX(0.d0,qx(i))*Nonlin_Pxq(j,i)-CMPLX(0.d0,qy(j))*psiPq(j,i)+psiP2q(j,i)+B1_2_4Dr*psi0*psiP_q(j,i)
                       ! B1_2_4Dr=B1*B1/(4*Dr), psiP_q=-q^2*(psi*Px)_q
     enddo
  enddo
!$OMP END PARALLEL DO

  PxPy=psi*Py
  call fftw_execute_dft_r2c(plan2d_f,PxPy,psiP_q) ! forward FFT
  psiP_q=-q2*psiP_q ! =-q^2*(psi*Py)_q
  psiP2q=psiP_q
  ! backward FFT (note: the input array is overwritten)
  call fftw_execute_dft_c2r(plan2d_b,psiP2q,PxPy)
  PxPy=PxPy*scale2d_b  ! = \nabla^2(psi*Py)
  
  PxPy=B1B20_4*(Px2_Py2*dy2psi+Px*dxPydyPx-Py*dxPx_dyPy)-Py*Nonlin_2+B1_2_4Dr*(Nonlin_1+PxPy*psi+psi0*Nonlin_psi0_2)
  psiP2q=-q2*Pyq
  ! backward FFT (note: the input array is overwritten)
  call fftw_execute_dft_c2r(plan2d_b,psiP2q,Nonlin_psi0_2) ! = \nabla^2 Py
  Nonlin_psi0_2=Nonlin_psi0_2*scale2d_b

  PxPy=PxPy+psi0*((B1B20_4)*Nonlin_Py_psi0+B1_2_4Dr*psi*Nonlin_psi0_2)
      ! Nonlin_Py_psi0=Px*(5\partial_y Px - 3\partial_x Py) - Py*(5\partial_x Px + 3\partial_y Py)
  call fftw_execute_dft_r2c(plan2d_f,PxPy,psiP2q) ! forward FFT
!$OMP PARALLEL DO
  do i=1,lx
     do j=1,ly/2+1
        Nonlin_Pyq(j,i)=CMPLX(0.d0,qy(j))*Nonlin_Pyq(j,i)-CMPLX(0.d0,qx(i))*psiPq(j,i)+psiP2q(j,i)+B1_2_4Dr*psi0*psiP_q(j,i)
                       ! B1_2_4Dr=B1*B1/(4*Dr), psiP_q=-q^2*(psi*Py)_q
     enddo
  enddo
!$OMP END PARALLEL DO

  call fftw_free(p_alloc_Nonlin_1)
  call fftw_free(p_alloc_Nonlin_2)
  call fftw_free(p_alloc_Nonlin_3)
  call fftw_free(p_alloc_dx2psi)
  call fftw_free(p_alloc_dy2psi)
  call fftw_free(p_alloc_dxPx_dyPy)
  call fftw_free(p_alloc_dxPydyPx)
  call fftw_free(p_alloc_psiPq)
  call fftw_free(p_alloc_psiP2q)
  call fftw_free(p_alloc_psiP_q)
  call fftw_free(p_alloc_Nonlin_psi0_1)
  call fftw_free(p_alloc_Nonlin_psi0_2)
  deallocate(PxPy,Px2_Py2,Nonlin_Px_psi0,Nonlin_Py_psi0)
  deallocate(Px2,Py2)

  return
end subroutine nonlin_q
