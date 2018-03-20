module konstants
   implicit none

   ! Allgemeine Konstanten
   real(kind=8),parameter :: Ru = 8.314459848 ! J/mol/K
end module

module dry_air
   use konstants
   implicit none

   ! Stoffeigenschaften: trockene Luft
   real(kind=8),parameter :: M=0.028949 ! kg/mol
   real(kind=8),parameter :: kappa=1.4 ! Isentropenexponent
   real(kind=8),parameter :: R = Ru/M ! J/kg/K
   real(kind=8),parameter :: p_krit = 3.766d6 ! Pa
   real(kind=8),parameter :: t_krit = 132.52 ! K
   real(kind=8),parameter :: rho_krit = 313 ! kg/m³
   real(kind=8),parameter :: omega = 0.03405 ! azentrischer Faktor
   real(kind=8),parameter :: kappa_SRK = 0.48508 + 1.5517*omega - 0.15613*omega*omega
   real(kind=8),parameter :: lambda = 0.0262 ! Leitfähigkeit in W/m/K
   real(kind=8),parameter :: cp = kappa*R/(kappa-1.0) ! J/kg/K
   real(kind=8),parameter :: eta = 17.1d-06 ! kg/m/s
end module

module ZGKonstants
   implicit none

   integer,parameter :: calcTempIG = 1
   integer,parameter :: calcPressIG = 2
   integer,parameter :: calcDensIg = 3
end module

module limKonstants
   implicit none

   integer,parameter :: noLim = 0
   integer,parameter :: minmod = 1
   integer,parameter :: superbee = 2
   integer,parameter :: vanLeer = 3
   integer,parameter :: charm = 4
end module

module fluxKonstants
   implicit none

   integer,parameter :: LLF = 1
   integer,parameter :: AUSM = 2
   integer,parameter :: AUSMplus = 3
end module

module preConKonstants
   implicit none

   integer,parameter :: yang = 1
   integer,parameter :: beam = 2
end module

module solver_vars
   implicit none

   ! Solvervariablen
   integer:: nx ! Gitterpunkte
   integer :: inner, outer ! Zeitschritt
   integer :: no ! Iterationsschritte (outer)
   integer :: ni ! Iterationsschritte (inner)
   real(kind=8) :: cfp
   real(kind=8) :: cfl
   real(kind=8) :: dx
   real(kind=8) :: dt
   real(kind=8) :: dtau
   real(kind=8),allocatable,dimension(:) :: x
   real(kind=8),allocatable,dimension(:,:) :: Q,Q_old
   real(kind=8),allocatable,dimension(:,:,:) :: Q_recon
   real(kind=8),allocatable,dimension(:) :: Q_imp
   real(kind=8),allocatable,dimension(:,:) :: Z
   real(kind=8),allocatable,dimension(:) :: Z_imp
   real(kind=8),allocatable,dimension(:,:) :: LHS_imp
   real(kind=8),allocatable,dimension(:) :: RHS_imp
   real(kind=8),allocatable,dimension(:,:) :: F
   real(kind=8),dimension(3,3) :: LHS
   real(kind=8),dimension(3) :: RHS
   real(kind=8),dimension(3,3) :: Trans ! Transformationsmatrix
   real(kind=8),dimension(3,3) :: PreCon ! PreCon-Matrix (Gamma)
   real(kind=8),dimension(3,3) :: A_w,A_e ! Jacobimatrix west/east
   real(kind=8),dimension(3) :: S
   real(kind=8),dimension(3) :: error
   real(kind=8) :: errorP,targetP
   real(kind=8) :: a
   real(kind=8) :: eps

end module

module routines

   contains

   subroutine doGit(x0,xn,nx,x,dx)
      implicit none

      real(kind=8),dimension(nx),intent(out) :: x
      real(kind=8),intent(in) :: x0,xn
      real(kind=8),intent(out) :: dx
      integer,intent(in) :: nx

      ! Lokale Variablen
      integer :: i

      x(1) = x0
      x(nx) = xn
      dx = (xn-x0)/(nx-1)
      do i=2,nx-1
         x(i) = x(i-1) + dx
      end do

   end subroutine

   subroutine startSod(Q,Z,kappa,R,nx)
      implicit none

      real(kind=8),dimension(nx,3),intent(out) :: Q,Z
      real(kind=8),intent(in) :: kappa,R
      integer,intent(in) :: nx

      ! Lokale Variablen
      integer :: i

   do i=1,nx/2
      Z(i,1) = 1.0d05
      Z(i,2) = 0.0
      Z(i,3) = 300.0
      call calcConsVar(Q(i,:),Z(i,:),kappa,R)
   end do
   do i=nx/2+1,nx
      Z(i,1) = 1.0d04
      Z(i,2) = 0.0
      Z(i,3) = 300.0
      call calcConsVar(Q(i,:),Z(i,:),kappa,R)
   end do
   end subroutine

   subroutine idealGas(p,rho,t,r,n)
      use ZGKonstants
      implicit none

      real(kind=8) :: p,rho,t
      real(kind=8),intent(in) :: r
      integer,intent(in) :: n

      if (n==calcPressIG) then 
         p = rho*t*r
      else if (n==calcDensIG) then
         rho=p/r/t
      else if (n==calcTempIG) then
         t=p/rho/r
      else
         write(*,*) "Error in idealGas: falscher Übergabewert für n"
         STOP 1
      end if
   end subroutine

   subroutine calcPrimVar(Q,Z,kappa,R)
      use ZGKonstants
      implicit none

      real(kind=8),dimension(3),intent(in) :: Q
      real(kind=8),dimension(3),intent(out) :: Z
      real(kind=8),intent(in) :: kappa
      real(kind=8),intent(in) :: R

      Z(2) = Q(2)/Q(1)
      Z(1) = (kappa-1.0)*(Q(3) - 0.5*Q(1)*Z(2)*Z(2))
      call idealGas(Z(1),Q(1),Z(3),R,calcTempIG)
   end subroutine

   subroutine calcConsVar(Q,Z,kappa,R)
      use ZGKonstants
      implicit none

      real(kind=8),dimension(3),intent(out) :: Q
      real(kind=8),dimension(3),intent(in) :: Z
      real(kind=8),intent(in) :: kappa
      real(kind=8),intent(in) :: R

      real(kind=8) :: cp

      cp = kappa*R/(kappa-1.0)

      call idealGas(Z(1),Q(1),Z(3),R,calcDensIG)
      Q(2) = Q(1)*Z(2)
      Q(3) = cp*Q(1)*Z(3) - Z(1) + 0.5*Q(1)*Z(2)**2
   end subroutine

   subroutine calcFLUX(Q,nx,kappa,R,F,typ)
      use fluxKonstants
      implicit none

      real(kind=8),dimension(2,nx-1,3),intent(in) :: Q ! 1=L,2=R;nx-1;3
      integer,intent(in) :: nx
      real(kind=8),intent(in) :: kappa
      real(kind=8),intent(in) :: R
      real(kind=8),dimension(nx-1,3),intent(out) :: F ! 1=west 2=east
      integer,intent(in) :: typ

      ! Lokale Variablen
      integer :: i

      if (typ == LLF) then
         do i=1,nx-1
            call calcllfflux(Q(1,i,:),Q(2,i,:),kappa,R,F(i,:))
         end do
      else if (typ == AUSM) then
         do i=1,nx-1
            call calcAusmflux(Q(1,i,:),Q(2,i,:),kappa,R,F(i,:))
         end do
      else if (typ == AUSMplus) then
         do i=1,nx-1
            call calcAusmplusflux(Q(1,i,:),Q(2,i,:),kappa,R,F(i,:))
         end do
      else
         write(*,*) "Unbekannter Fluss, Programm wird beendet"
         stop 1
      end if

   end subroutine

   subroutine reconQ(Q,Q_recon,nx,ord,lim)
      implicit none

      real(kind=8),dimension(nx,3),intent(in) :: Q
      real(kind=8),dimension(2,nx-1,3),intent(out) :: Q_recon ! 1=L,2=R;nx-1;3
      integer,intent(in) :: nx
      integer,intent(in) :: ord
      integer,intent(in) :: lim

      ! Lokale Variablen
      integer :: i

      Q_recon(1,1,:) = Q(1,:)
      Q_recon(2,nx-1,:) = Q(nx,:)

      ! Linke Werte
      do i=2,nx-1
         call reconQL(Q(i-1:i+1,:),Q_recon(1,i,:),ord,lim)
      end do

      ! Rechte Werte
      do i=1,nx-2
         call reconQR(Q(i:i+2,:),Q_recon(2,i,:),ord,lim)
      end do

   end subroutine

   subroutine reconQR(Q,QR,ord_in,lim) !MUSCL?
      use limKonstants
      implicit none

      real(kind=8),dimension(3,3),intent(in) :: Q
      real(kind=8),dimension(3),intent(out) :: QR
      integer,intent(in) :: ord_in
      integer,intent(in) :: lim

      ! Lokale Variablen
      integer :: i,ord
      real(kind=8),dimension(3) :: phi
      real(kind=8),dimension(3) :: r
      real(kind=8),dimension(3) :: beta

      if (ord_in == 1) then
         QR = Q(2,:)
      else
         if (minval(abs(q(3,:)-q(2,:))) == 0.0) then
            r = 1d0
            ord = 2
         else
            r = (q(2,:)-q(1,:))/(q(3,:)-q(2,:))
            ord = ord_in
         end if
         if (ord == 2) then
            beta = 1d0
         else if (ord == 3) then
            beta = (1d0 +2d0*r)/3d0
         else
            write(*,*) "Nur Ordnungen 1-3 Implementiert"
            stop 1
         end if

         do i=1,3
            if (q(2,i)-q(1,i) < 1d-8) then
               r(i) = 0d0
            else
               r(i) = (q(3,i)-q(2,i))/(q(2,i)-q(1,i))
            end if
         end do
         call limiter(r,phi,lim)
         QR = Q(2,:) - 0.5d0*phi*beta*(Q(3,:) - Q(2,:))
      end if

   end subroutine

   subroutine reconQL(Q,QL,ord_in,lim) !MUSCL?
      use limKonstants
      implicit none

      real(kind=8),dimension(3,3),intent(in) :: Q
      real(kind=8),dimension(3),intent(out) :: QL
      integer,intent(in) :: ord_in
      integer,intent(in) :: lim

      ! Lokale Variablen
      integer :: i,ord
      real(kind=8),dimension(3) :: phi
      real(kind=8),dimension(3) :: r
      real(kind=8),dimension(3) :: beta

      if (ord_in == 1) then
         QL = Q(2,:)
      else
         if (minval(abs(q(2,:)-q(1,:))) == 0.0) then
            r = 1d0
            ord = 2
         else
            r = (q(3,:)-q(2,:))/(q(2,:)-q(1,:))
            ord = ord_in
         end if
         if (ord == 2) then
            beta = 1d0
         else if (ord == 3) then
            beta = (1d0 +2d0*r)/3d0
         else
            write(*,*) "Nur Ordnungen 1-3 Implementiert"
            stop 1
         end if

         do i=1,3
            if (q(3,i)-q(2,i) < 1d-8) then
               r(i) = 0d0
            else
               r(i) = (q(2,i)-q(1,i))/(q(3,i)-q(2,i))
            end if
         end do
         call limiter(r,phi,lim)
         QL = Q(2,:) + 0.5d0*beta*phi*(Q(2,:) - Q(1,:))
      end if

   end subroutine

   subroutine reconstructLeftValue (Q1,Q2,Q3,QL)
      ! Rekonstruktion des Zellwandwertes auf der rechnten Seite nach 3.Ordnung
      ! Für i+1/2 gilt: Q1 = Q(i-1); Q2 = Q(i); Q3 = Q(i+1)
      implicit none

      real(kind=8),dimension(3),intent(in) :: Q1
      real(kind=8),dimension(3),intent(in) :: Q2
      real(kind=8),dimension(3),intent(in) :: Q3
      real(kind=8),dimension(3),intent(out) :: QL

      QL = 0.375*Q3 + 0.75*Q2 - 0.125*Q1
   end subroutine

   subroutine reconstructRightValue (Q1,Q2,Q3,QR)
      ! Rekonstruktion des Zellwandwertes auf der linken Seite nach 3.Ordnung
      ! Für i+1/2 gilt: Q1 = Q(i); Q2 = Q(i+1); Q3 = Q(i+2)
      implicit none

      real(kind=8),dimension(3),intent(in) :: Q1
      real(kind=8),dimension(3),intent(in) :: Q2
      real(kind=8),dimension(3),intent(in) :: Q3
      real(kind=8),dimension(3),intent(out) :: QR

      QR = -0.125*Q3 + 0.75*Q2 + 0.375*Q1
   end subroutine

   subroutine calcStdFlux(Q,kappa,R,F)
      implicit none

      real(kind=8),dimension(3),intent(in) :: Q
      real(kind=8),intent(in) :: kappa
      real(kind=8),intent(in) :: R
      real(kind=8),dimension(3),intent(out) :: F

      !Lokale Variablen
      real(kind=8),dimension(3) :: Z

      call calcPrimVar(Q,Z,kappa,R)

      F(1) = Q(2)
      F(2) = Q(1)*Z(2)**2 + Z(1)
      F(3) = Z(2)*(Q(3) + Z(1))
   end subroutine
   
   !subroutine calcLLFFlux(Q1,Q2,kappa,R,lambda,dx,dt,F)
   !subroutine calcLLFFlux(Q1,Q2,kappa,R,dx,dt,F)
   subroutine calcLLFFlux(Q1,Q2,kappa,R,F)
      use solver_vars, only: dx,dt ! DONT CHANGE!
      implicit none

      real(kind=8),dimension(3),intent(in) :: Q1
      real(kind=8),dimension(3),intent(in) :: Q2
      real(kind=8),intent(in) :: kappa
      real(kind=8),intent(in) :: R
      !real(kind=8),intent(in) :: lambda
      !real(kind=8),intent(in) :: dx
      !real(kind=8),intent(in) :: dt
      real(kind=8),dimension(3),intent(out) :: F

      !Lokale Variablen
      real(kind=8),dimension(3) :: Z1,Z2
      real(Kind=8) :: a

      a = dx/dt
      call calcPrimVar(Q1,Z1,kappa,R)
      call calcPrimVar(Q2,Z2,kappa,R)

      F(1) = 0.5*(Q1(2)+Q2(2)) - 0.5*a*(Q2(1)-Q1(1))
      F(2) = 0.5*(Q1(1)*Z1(2)**2 + Q2(1)*Z2(2)**2 + Z1(1) + Z2(1)) &
            -0.5*a*(Q2(2) - Q1(2))
      F(3) = 0.5*(Z1(2)*(Q1(3)+Z1(1)) + Z2(2)*(Q2(3)+Z2(1))) &
            -0.5*a*(Q2(3) - Q1(3)) !+ lambda*(Z2(3)-Z1(3))/dx
   end subroutine
   
   subroutine calcFDFlux(Q1,Q2,kappa,R,F)
      implicit none

      real(kind=8),dimension(3),intent(in) :: Q1
      real(kind=8),dimension(3),intent(in) :: Q2
      real(kind=8),intent(in) :: kappa
      real(kind=8),intent(in) :: R
      !real(kind=8),intent(in) :: lambda
      real(kind=8),dimension(3),intent(out) :: F

      !Lokale Variablen
      real(kind=8),dimension(3) :: Z1,Z2
      real(kind=8) :: E1
      real(kind=8) :: E2

      call calcPrimVar(Q1,Z1,kappa,R)
      call calcPrimVar(Q2,Z2,kappa,R)
      E1 = Q1(3)/Q1(1)
      E2 = Q2(3)/Q2(1)

      F(1) = 2.0*Q2(1)*Z2(2) - Z2(2)*Q1(1) - Z1(2)*Q2(1)
      F(2) = 2.0*Q2(1)*Z2(2)*(Z2(2)-Z1(2)) + Z2(2)**2*(Q2(1)-Q1(1)) + Z2(1) - Z1(1)
      F(3) = ((Q2(3)+Z2(1))*(Z2(2)-Z1(2))  &
                     + Z2(2)*Q2(1)*(E2-E1) + Z2(2)*E2*(Q2(1)-Q1(1))  &
                     + Z2(2)*(Z2(1)-Z1(1))) ! - lambda*(Z2(3)-Z1(3))/dx
   end subroutine

   subroutine calcAUSMfLux(QL,QR,kappa,R,F)
      implicit none

      real(kind=8),dimension(3),intent(in) :: QL
      real(kind=8),dimension(3),intent(in) :: QR
      real(kind=8),intent(in) :: kappa
      real(kind=8),intent(in) :: R
      real(kind=8),dimension(3),intent(out) :: F

      !Lokale Variablen
      real(kind=8),dimension(3) :: ZL,ZR,Phi
      real(kind=8) :: aL,aR,ht
      real(kind=8) :: ML,MR,MpL,MmR,M
      real(kind=8) :: p,PpL,PmR

      call calcPrimVar(QL,ZL,kappa,R)
      call calcPrimVar(QR,ZR,kappa,R)

      aL = sqrt(kappa*R*ZL(3))
      aR = sqrt(kappa*R*ZR(3))
      ML = ZL(2)/aL
      MR = ZR(2)/aR

      if (abs(ML) > 1d0) then
         MpL = 0.5d0*(ML + abs(ML))
         PpL = 0.5d0*(1 + sign(1d0,ML))
      else
         MpL = 0.25d0*(ML + 1d0)**2
         PpL = 0.25d0*(ML + 1d0)**2 *(2d0 - ML)
      end if
      if (abs(MR) > 1d0) then
         MmR = 0.5d0*(MR - abs(MR))
         PmR = 0.5d0*(1 - sign(1d0,MR))
      else
         MmR = -0.25d0*(MR - 1d0)**2
         PmR = 0.25d0*(MR - 1d0)**2 *(2d0 + MR)
      end if

      M = MpL + MmR
      p = PpL*ZL(1) + PmR*ZR(1)

      if (M >= 0) then
         ht = aL**2/(kappa-1) + 0.5d0*ZL(2)**2
         Phi(1) = QL(1)*aL
         Phi(2) = Phi(1)*ZL(2)
         Phi(3) = Phi(1)*ht
      else 
         ht = aR**2/(kappa-1) + 0.5d0*ZR(2)**2
         Phi(1) = QR(1)*aR
         Phi(2) = Phi(1)*ZR(2)
         Phi(3) = Phi(1)*ht
      end if

      F = M*Phi
      F(2) = F(2) + p

   end subroutine

   subroutine calcAUSMPlusFLux(QL,QR,kappa,R,F)
      implicit none

      real(kind=8),dimension(3),intent(in) :: QL
      real(kind=8),dimension(3),intent(in) :: QR
      real(kind=8),intent(in) :: kappa
      real(kind=8),intent(in) :: R
      !real(kind=8),intent(in) :: lambda
      real(kind=8),dimension(3),intent(out) :: F

      !Lokale Variablen
      real(kind=8),dimension(3) :: ZL,ZR
      real(kind=8) :: a,a_crit_L,a_crit_R
      real(kind=8) :: HL,HR,H_t_L,H_t_R
      real(kind=8) :: ML,MR,MpL,MmR,Mpress,M2,M,Mp,Mm,M0
      real(kind=8) :: p,PpL,PmR,Pu
      real(kind=8) :: fa,alpha

      !Parameter
      real(kind=8),parameter :: Kp = 0.25d0
      real(kind=8),parameter :: Ku = 0.75d0
      real(kind=8),parameter :: sigma = 1d0
      real(kind=8),parameter :: beta = 0.125d0
      real(kind=8),parameter :: Mref = 0.1d0

      call calcPrimVar(QL,ZL,kappa,R)
      call calcPrimVar(QR,ZR,kappa,R)

      HL = (QL(3)+ZL(1))/QL(1)! + 0.5d0*ZL(2)**2
      HR = (QR(3)+ZR(1))/QR(1)! + 0.5d0*ZR(2)**2

      H_t_L = kappa/(kappa-1d0)*ZL(1)/QL(1) + 0.5d0*ZL(2)**2
      H_t_R = kappa/(kappa-1d0)*ZR(1)/QR(1) + 0.5d0*ZR(2)**2

      a_crit_L = sqrt(2d0*(kappa-1d0)/(kappa+1d0)*H_t_L)
      a_crit_R = sqrt(2d0*(kappa-1d0)/(kappa+1d0)*H_t_R)

      a = min(a_crit_L**2/max(a_crit_L,ZL(2)),a_crit_R**2/max(a_crit_R,-ZR(2)))

      ML = ZL(2)/a
      MR = ZR(2)/a
      M2 = 0.5*(ZL(2)**2 + ZR(2)**2)/a**2 !a**2 richtig ?!

      if (abs(ML)<1d0) then
         MpL = 0.25d0*(ML+1d0)**2 + beta*(ML**2-1d0)**2
      else
         MpL = 0.5d0*(ML + abs(ML))
      end if
      if (abs(MR)<1d0) then
         MmR = -0.25d0*(MR-1d0)**2 - beta*(MR**2-1d0)**2
      else
         MmR = 0.5d0*(MR - abs(MR))
      end if

      M0 = sqrt(min(1d0,max(M2,Mref**2)))
      fa = M0*(2d0 - M0)
      Mpress = -Kp/fa*max(1d0-sigma*M2,0d0)*(ZR(1)-ZL(1))/(0.5*(QL(1)+QR(1))*a**2)

      M = MpL + MmR + Mpress
      Mp = 0.5d0*(M + abs(M))
      Mm = 0.5d0*(M - abs(M))

      alpha = 3d0/16d0*(-4d0 +5d0*fa**2)
      if (abs(ML)<1d0) then
         PpL = 0.5d0*(ML+1d0)**2 -0.25d0*ML*(ML+1d0)**2 +alpha*ML*(ML**2-1d0)**2
      else
         PpL = 0.5d0/ML*(ML + abs(ML))
      end if
      if (abs(MR)<1d0) then
         PmR = 0.5d0*(MR-1d0)**2 +0.25d0*MR*(MR-1d0)**2 -alpha*MR*(MR**2-1d0)**2
      else
         PmR = 0.5d0/MR*(MR + abs(MR))
      end if

      Pu = -fa*Ku*PpL*PmR*(QL(1)+QR(1))*a*(ZR(2)-ZL(2))
      p = PpL*ZL(1) + PmR*ZR(1) + Pu

      F(1) = a*(Mp*QL(1) +Mm*QR(1))
      F(2) = a*(Mp*QL(1)*ZL(2) +Mm*QR(1)*ZR(2)) + p
      F(3) = a*(Mp*QL(1)*HL +Mm*QR(1)*HR)

   end subroutine
   
   subroutine calcQuell(Q,eta,S)
      implicit none

      real(kind=8),dimension(3),intent(in) :: Q
      real(kind=8),intent(in) :: eta
      real(kind=8),dimension(3),intent(out) :: S

      ! Lokale Variablen
      real(kind=8),parameter :: d_rohr = 0.1 !m
      real(kind=8) :: Re !Reynoldszahl

      Re = Q(2)*d_rohr/eta

      S(1) = 0.0
      S(2) = -100.0
      !S(2) = -64.0/Re*1000.0 ! Laminarströmun im Rohr
      S(3) = 0.0
   end subroutine

   subroutine calcJacobiKonv(Q,kappa,jacobi)
      use dry_Air , only: R
      implicit none

      real(kind=8),dimension(3),intent(in) :: Q
      real(kind=8),intent(in) :: kappa
      real(kind=8),dimension(3,3),intent(out) :: jacobi

      ! Lokale Variablen
      real(kind=8),dimension(3) :: Z

      call calcPrimVar(Q,Z,kappa,R)

      jacobi(1,1) = 0.0
      jacobi(1,2) = 1.0
      jacobi(1,3) = 0.0

      jacobi(2,1) = (kappa-3d0)/2d0*Z(2)**2
      jacobi(2,2) = (3d0-kappa)*Z(2)
      jacobi(2,3) = kappa-1d0

      jacobi(3,1) = -kappa*Z(2)*Q(3)/Q(1) + (kappa-1d0)*Z(2)**3
      jacobi(3,2) = kappa*Q(3)/Q(1) - 3d0/2d0*(kappa-1d0)*Z(2)**2
      jacobi(3,3) = kappa*Z(2)
   end subroutine

   subroutine calcJacobiKonvPrimitivIG(Z,R,kappa,jacobi)
      implicit none

      real(kind=8),dimension(3),intent(in) :: Z
      real(kind=8),intent(in) :: R
      real(kind=8),intent(in) :: kappa
      real(kind=8),dimension(3,3),intent(out) :: jacobi

      !Lokale Variablen
      real(kind=8),dimension(3) :: Q
      real(kind=8) :: cp
      real(kind=8) :: Htot
      real(kind=8) :: rho_P
      real(kind=8) :: rho_T

      call calcConsVar(Q,Z,kappa,R)
      cp    = kappa*R/(kappa-1.0)
      ! IG-Gleichungen:
      Htot  = (Q(3) + Z(1))/Q(1)
      rho_P = 1.0/R/Z(3)
      rho_T = -Q(1)/Z(3)

      jacobi(1,1) = Z(2)*rho_P
      jacobi(1,2) = Q(1)
      jacobi(1,3) = Z(2)*rho_T

      jacobi(2,1) = Z(2)*Z(2)*rho_P + 1.0
      jacobi(2,2) = 2.0*Q(1)*Z(2)
      jacobi(2,3) = Z(2)*Z(2)*rho_T

      jacobi(3,1) = Z(2)*(rho_P*Htot - 1.0)
      jacobi(3,2) = Q(3) + Z(1) + Q(1)*Z(2)**2
      jacobi(3,3) = Z(2)*(rho_T*Htot + Q(1)*cp)
   end subroutine

   subroutine calcJacobiDiff(Q,kappa,R,jacobi)
      implicit none

      real(kind=8),dimension(3),intent(in) :: Q
      real(kind=8),intent(in) :: kappa
      real(kind=8),intent(in) :: R
      real(kind=8),dimension(3,3),intent(out) :: jacobi

      jacobi(1,1) = 0.0
      jacobi(1,2) = 0.0
      jacobi(1,3) = 0.0

      jacobi(2,1) = 0.0
      jacobi(2,2) = 0.0
      jacobi(2,3) = 0.0

      jacobi(3,1) = (kappa-1.0)/R*Q(3)
      jacobi(3,2) = -(kappa-1.0)/R*Q(2)
      jacobi(3,3) = (kappa-1.0)/R*Q(1)
   end subroutine

   !subroutine calcJacobiQuell(Q,jacobi)
      !implicit none

      !real(kind=8),dimension(3),intent(in) :: Q
      !real(kind=8),dimension(3,3),intent(out) :: jacobi

      !jacobi(1,1) = 0.0
      !jacobi(1,2) = 0.0
      !jacobi(1,3) = 0.0

      !jacobi(2,1) = 0.0 !-konst/2.0/Q(2)/Q(2)
      !jacobi(2,2) = 0.0 !+konst*Q(1)/Q(2)
      !jacobi(2,3) = 0.0

      !jacobi(3,1) = 0.0
      !jacobi(3,2) = 0.0
      !jacobi(3,3) = 0.0
   !end subroutine

   subroutine calcTransMatIG(Q,R,kappa,T)
      implicit none

      real(kind=8),dimension(3),intent(in) :: Q
      real(kind=8),intent(in) :: R
      real(kind=8),intent(in) :: kappa
      real(kind=8),dimension(3,3),intent(out) :: T

      ! Lokale Variablen
      real(kind=8),dimension(3) :: Z
      real(kind=8) :: cp
      real(kind=8) :: H

      call calcPrimVar(Q,Z,kappa,R)
      cp = kappa*R/(kappa-1.0)
      H = cp*Z(3)+Z(2)**2/2.0

      T(1,1) = 1.0/R/Z(3)
      T(1,2) = 0.0
      T(1,3) = -Z(1)/R/Z(3)/Z(3)

      T(2,1) = T(1,1)*Z(2)
      T(2,2) = Q(1)
      T(2,3) = T(1,3)*Z(2)

      T(3,1) = T(1,1)*H-1.0
      T(3,2) = Q(2)
      T(3,3) = T(1,3)*H + Q(1)*cp
   end subroutine

   subroutine calcPreConMatIG(Q,R,kappa,eps,T,mode)
      use preConKonstants
      implicit none

      real(kind=8),dimension(3),intent(in) :: Q
      real(kind=8),intent(in) :: R
      real(kind=8),intent(in) :: kappa
      real(kind=8),intent(in) :: eps
      real(kind=8),dimension(3,3),intent(out) :: T
      integer,intent(in) :: mode

      ! Lokale Variablen
      real(kind=8) :: Theta
      real(kind=8) :: Ur
      real(kind=8),dimension(3) :: Z
      real(kind=8) :: cp
      real(kind=8) :: a,a2
      real(kind=8) :: H

      call calcPrimVar(Q,Z,kappa,R)
      cp = kappa*R/(kappa-1.0)
      H = cp*Z(3)+Z(2)**2/2.0
      a2 = kappa*R*Z(3) ! Quadrat der Schallgeschwindigkeit
      a = sqrt(a2)

      if (mode == yang) then
         ! aus Yang 2004 (gamma/beta)
         Theta = (1+(kappa-1)*eps)/(eps*a2)
      else if (mode == beam) then
         ! aus Beam 1995
         if (abs(Z(2)) < 1d-5*a2) then
            Ur = 1d-5*sqrt(a2)
         else if (abs(Z(2)) < a) then
            Ur = abs(Z(2))
         else
            Ur = a
         end if
         Theta = 1d0/(Ur**2) + 1d0/(cp*Z(3))
      else
         write(*,*) "Unbekannter mode in calcPreConMatIG"
         stop 1
      end if

      T(1,1) = Theta
      T(1,2) = 0.0
      T(1,3) = -Z(1)/R/Z(3)/Z(3)

      T(2,1) = Theta*Z(2)
      T(2,2) = Q(1)
      T(2,3) = T(1,3)*Z(2)

      T(3,1) = Theta*H - 1.0
      T(3,2) = Q(2)
      T(3,3) = T(1,3)*H + Q(1)*cp
   end subroutine

   subroutine calcEps(Z,l_x,dt,R,kappa,eps)
      implicit none

      real(kind=8),dimension(3),intent(in) :: Z
      real(kind=8),intent(in) :: l_x
      real(kind=8),intent(in) :: dt
      real(kind=8),intent(in) :: R
      real(kind=8),intent(in) :: kappa
      real(kind=8),intent(out) :: eps

      ! Lokale Variablen
      real(kind=8),parameter :: e = 1d-5
      real(kind=8) :: eps_inv
      real(kind=8) :: eps_dt
      real(kind=8) :: eps_vis
      real(kind=8) :: M
      real(kind=8) :: a

      a = sqrt(kappa*R*Z(3))
      M = Z(2)/a

      if(M<=e) then
         eps_inv = e*e
      else if (M<1) then
         eps_inv = 2.0*M*M
      else
         eps_inv = 1.0
      end if

      eps_dt = (l_x/3.1415926/dt)**2/a**2 + M*M
      eps_vis = 0.0 ! Keine Viskosität
      eps = min(1d0,eps_inv) !min(1d0,max(eps_inv,eps_dt,eps_vis))
   end subroutine

   subroutine limiter(r,phi,lim)
      use limKonstants
      implicit none

      real(kind=8),dimension(3),intent(in) :: r
      real(kind=8),dimension(3),intent(out) :: phi
      integer,intent(in) :: lim

      if (lim == noLim) then
         phi = 1d0
      else if (lim == minmod) then
         phi = max(0d0,min(1d0,r))
      else if (lim == superbee) then
         phi = max(0d0,min(2.0*r,1d0),min(r,2d0))
      else if (lim == vanLeer) then
         phi = (r + abs(r))/(1d0 + abs(r))
      else if (lim == charm) then
         if (minval(r)>0) then
            phi = (r*(3d0*r+1d0))/(r+1d0)**2
         else
            phi = 0
         end if
      else
         write(*,*) "Limiter nicht erkannt, STOP"
         stop 1
      end if

   end subroutine

   subroutine matinv(A,B)
      ! Invertierung einer 3x3 Matrix
      implicit none
      !! Performs a direct calculation of the inverse of a 3×3 matrix.
      real(kind=8), intent(in) :: A(3,3)   !! Matrix
      real(kind=8), intent(out) :: B(3,3)   !! Matrix
      real(kind=8)             :: detinv
   
      ! Calculate the inverse determinant of the matrix
      detinv = 1/(A(1,1)*A(2,2)*A(3,3) - A(1,1)*A(2,3)*A(3,2)&
                - A(1,2)*A(2,1)*A(3,3) + A(1,2)*A(2,3)*A(3,1)&
                + A(1,3)*A(2,1)*A(3,2) - A(1,3)*A(2,2)*A(3,1))
   
      ! Calculate the inverse of the matrix
      B(1,1) = +detinv * (A(2,2)*A(3,3) - A(2,3)*A(3,2))
      B(2,1) = -detinv * (A(2,1)*A(3,3) - A(2,3)*A(3,1))
      B(3,1) = +detinv * (A(2,1)*A(3,2) - A(2,2)*A(3,1))
      B(1,2) = -detinv * (A(1,2)*A(3,3) - A(1,3)*A(3,2))
      B(2,2) = +detinv * (A(1,1)*A(3,3) - A(1,3)*A(3,1))
      B(3,2) = -detinv * (A(1,1)*A(3,2) - A(1,2)*A(3,1))
      B(1,3) = +detinv * (A(1,2)*A(2,3) - A(1,3)*A(2,2))
      B(2,3) = -detinv * (A(1,1)*A(2,3) - A(1,3)*A(2,1))
      B(3,3) = +detinv * (A(1,1)*A(2,2) - A(1,2)*A(2,1))
   end subroutine

   subroutine solve (LHS,RHS,dimen) ! Gauß-Elimination
      implicit none
   
      integer,intent(in) :: dimen
      real(kind=8),dimension (dimen, dimen),intent (in) :: LHS
      real(kind=8),dimension (dimen),intent (inout) :: RHS
      
      ! Lokale Variablen
      
      integer :: i,j,k
      real(kind=8),dimension (dimen,dimen) :: A
      real(kind=8),dimension (dimen) :: b
      
      b = RHS
      A = LHS

      ! Obere Diagonalmatrix
      do i=1,dimen-1
         if (A(i,i) == 0) then
            write(*,*) "Matrix hat nicht vollen Rang! STOP"
            stop 1
         else
            b(i) = b(i)/A(i,i)
            do j=dimen,i,-1
               A(i,j) = A(i,j)/A(i,i)
            end do
         end if
         do j= i+1, dimen
            b(j) = A(i,i)*b(j) - A(j,i)*b(i)
            do k=dimen,i+1,-1
               A(j,k) = A(j,k)*A(i,i) - A(i,k)*A(j,i)
            end do
            A(j,i) = 0d0
         end do
      end do
      if (A(dimen,dimen) == 0) then
         write(*,*) "Matrix hat nicht vollen Rang! STOP"
         stop 1
      else
         b(dimen) = b(dimen)/A(dimen,dimen)
         A(dimen,dimen) = 1.0
      end if

      ! Rückwärts einsetzen
      do i=dimen-1,1,-1
         do j=i+1,dimen
            b(i) = b(i) - A(i,j)*b(j)
         end do
      end do

      RHS = b
   end subroutine

   subroutine solveGS(LHS,RHS,dimen)
      implicit none
   
      integer,intent(in) :: dimen
      real(kind=8),dimension (dimen, dimen),intent (in) :: LHS
      real(kind=8),dimension (dimen),intent (inout) :: RHS
      
      ! Lokale Variablen
      
      integer :: i,j
      real(kind=8),dimension (dimen,dimen) :: D,LD,UD

      ! LU Matrix System erstellen
      D = 0.0
      LD = 0.0
      UD = 0.0
      do i=1,dimen
         D(i,i) = LHS(i,i)
         LD(i,i) = LHS(i,i)
         UD(i,i) = LHS(i,i)
      end do
      do i=2,dimen
         do j=1,i-1
            LD(i,j) = LHS(i,j)
            UD(j,i) = LHS(j,i)
         end do
      end do

      call lSweep(LD,RHS,dimen)
      do i=1,dimen
         RHS(i) = D(i,i)*RHS(i)
      end do
      call uSweep(UD,RHS,dimen)

      ! Solve
   end subroutine

   subroutine lSweep(A,b,dimen)
      implicit none

      integer,intent(in) :: dimen
      real(kind=8),dimension (dimen, dimen),intent (in) :: A
      real(kind=8),dimension (dimen),intent (inout) :: b

      !Lokale Variablen
      integer :: i,j

      ! Vorwärts einsetzen
      do i=1,dimen
         if (A(i,i) == 0) then
            write(*,*) "Diagonalelement ",i," ist 0; STOP"
            stop 1
         end if
         do j=1,i-1
            b(i) = b(i) - A(i,j)*b(j)
         end do
         b(i) = b(i)/A(i,i)
      end do
   end subroutine

   subroutine uSweep(A,b,dimen)
      implicit none

      integer,intent(in) :: dimen
      real(kind=8),dimension (dimen, dimen),intent (in) :: A
      real(kind=8),dimension (dimen),intent (inout) :: b

      !Lokale Variablen
      integer :: i,j

      ! Rückwärts einsetzen
      do i=dimen,1,-1
         if (A(i,i) == 0) then
            write(*,*) "Diagonalelement ",i," ist 0; STOP"
            stop 1
         end if
         do j=i+1,dimen
            b(i) = b(i) - A(i,j)*b(j)
         end do
         b(i) = b(i)/A(i,i)
      end do
   end subroutine

   subroutine datOut(x,Q,Z,nx,nStep)
      implicit none

      real(kind=8),dimension(nx),intent(in) :: x
      real(kind=8),dimension(nx,3),intent(in) :: Q
      real(kind=8),dimension(nx,3),intent(in) :: Z
      integer,intent(in) :: nx
      integer,intent(in) :: nStep

      ! Lokale Variablen
      integer :: i

      do i=1,nx
         write(1000+nStep,*) x(i),Q(i,:),Z(i,:)
      end do
   end subroutine
end module
