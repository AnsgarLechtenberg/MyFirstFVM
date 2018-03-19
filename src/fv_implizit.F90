program fv_implizit
   use Konstants
   use ZGKonstants
   use fluxKonstants
   use limKonstants
   use dry_air
   use solver_vars
   use routines

   implicit none

   ! Iterationsvariablen
   integer :: i,info
   integer,allocatable,dimension(:) :: ipiv

   ! Konstanten setzen
   cfl = 0.5
   no = 160
   nx = 501

   allocate(ipiv(3*nx))
   allocate(x(nx),Q(nx,3),Q_recon(2,nx-1,3),Z(nx,3),Q_imp(3*nx),F(nx-1,3))
   allocate(LHS_imp(3*nx,3*nx),RHS_imp(3*nx))

   ! 1D-Gitter
   call doGit(0d0,1d1,nx,x,dx)

   ! Anfangsbedingungen
   call startSod(Q,Z,kappa,R,nx)
   call datOut(x,Q,Z,nx,0)

   do i=1,nx
      Q_imp((i-1)*3+1) = Q(i,1)
      Q_imp((i-1)*3+2) = Q(i,2)
      Q_imp((i-1)*3+3) = Q(i,3)
   end do

   ! Solve
   do outer=1,no
      write(*,*) "Zeitschritt: ",outer

      dt = dx*cfl/sqrt(kappa*R*minval(Z(:,3))) ! Zeitschritt

      ! RHS und LHS nullen
      RHS_imp = 0.0
      LHS_imp = 0.0

      ! Zero-Gradient Rand
      RHS_imp(1:3) = 0.0
      RHS_imp(nx-2:nx) = 0.0
      LHS_imp(1,1) = 1.0
      LHS_imp(2,2) = 1.0
      LHS_imp(3,3) = 1.0
      LHS_imp(1,4) = -1.0
      LHS_imp(2,5) = -1.0
      LHS_imp(3,6) = -1.0
      LHS_imp(3*nx,3*nx) = 1.0
      LHS_imp(3*nx-1,3*nx-1) = 1.0
      LHS_imp(3*nx-2,3*nx-2) = 1.0
      LHS_imp(3*nx,3*nx-3) = -1.0
      LHS_imp(3*nx-1,3*nx-4) = -1.0
      LHS_imp(3*nx-2,3*nx-5) = -1.0

      call reconQ(Q,Q_recon,nx,3)
      call calcFlux(Q_recon,nx,kappa,R,F,AUSM)
      
      do i=2,nx-1 ! Schleife über Zellen
         ! Berechnung RHS
         RHS_imp((i-1)*3+1) = +1.0/dx*(F(i-1,1)-F(i,1))
         RHS_imp((i-1)*3+2) = +1.0/dx*(F(i-1,2)-F(i,2))
         RHS_imp((i-1)*3+3) = +1.0/dx*(F(i-1,3)-F(i,3))

         ! Berechnung LHS
         LHS_imp((i-1)*3+1,(i-1)*3+1) = +1.0/dt
         LHS_imp((i-1)*3+2,(i-1)*3+2) = +1.0/dt
         LHS_imp((i-1)*3+3,(i-1)*3+3) = +1.0/dt

         call calcJacobiKonv(Q(i-1,:),kappa,A_w)
         call calcJacobiKonv(Q(i,:),kappa,A_e)

         LHS_imp((i-1)*3+1:(i-1)*3+3,(i-2)*3+1:(i-2)*3+3) = -A_w/dx/1d0
         LHS_imp((i-1)*3+1:(i-1)*3+3,(i)*3+1:(i)*3+3) = A_e/dx/1d0
      end do
      
      ! Lösen des Gleichungssystems
      call dgesv(3*nx,1,LHS_imp,3*nx,ipiv,RHS_imp,3*nx,info)

      Q_imp = Q_imp + RHS_imp

      ! Übernhame der alten Werte
      do i=1,nx
         Q(i,1) = Q_imp((i-1)*3+1)
         Q(i,2) = Q_imp((i-1)*3+2)
         Q(i,3) = Q_imp((i-1)*3+3)
         call calcPrimVar(Q(i,:),Z(i,:),kappa,R)
      end do

      call datOut(x,Q,Z,nx,outer)
   end do
end program
