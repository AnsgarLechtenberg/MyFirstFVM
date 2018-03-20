program fv_preCon_implizit
   use konstants
   use ZGKonstants
   use fluxKonstants
   use limKonstants
   use preConKonstants
   use dry_air
   use solver_vars
   use routines

   implicit none

   ! Iterationsvariablen
   integer :: i,info!,ieps,icfp
   integer :: totIter
   integer,allocatable,dimension(:) :: ipiv

   ! Solverparameter
   no = 20
   ni = 50
   nx = 201
   cfl = 2.0
   cfp = 0.7

   allocate(ipiv(3*nx))
   allocate(x(nx),Q(nx,3),Q_old(nx,3),Q_recon(2,nx-1,3),Z(nx,3),Z_imp(3*nx)  &
           ,LHS_imp(3*nx,3*nx),RHS_imp(3*nx),F(nx-1,3))

   ! 1D-Gitter
   call doGit(0d0,1d1,nx,x,dx)

   ! Anfangsbedingungen
   call startSod(Q,Z,kappa,R,nx)
   do i=1,nx
      Z_imp((i-1)*3+1) = Z(i,1)
      Z_imp((i-1)*3+2) = Z(i,2)
      Z_imp((i-1)*3+3) = Z(i,3)
   end do

   call datOut(x,Q,Z,nx,0)

   ! Solve

   totIter = 0
   do outer=1,no ! Äußere Iteraton

      Q_old = Q

      dt = dx*cfl/sqrt(kappa*R*minval(Z(:,3))) ! Zeitschritt
      write(*,*) "Zeitschritt:",outer,dt

      inner = 0
      errorP = 1
      targetP = 0
      do while ((errorP>targetP) .and. (inner<ni))
         inner = inner + 1
         totIter = totIter + 1
         error = 0.0

         ! RHS und LHS nullen
         RHS_imp = 0.0
         LHS_imp = 0.0

         ! Zero-Gradient Rand
         RHS_imp(1) = 0.0
         RHS_imp(2) = 0.0
         RHS_imp(3) = 0.0
         RHS_imp(nx-2) = 0.0
         RHS_imp(nx-1) = 0.0
         RHS_imp(nx) = 0.0
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

         call reconQ(Q,Q_recon,nx,3,minmod)
         call calcFlux(Q_recon,nx,kappa,R,F,llf)

         do i=2,nx-1 ! Schleife über Zellen
            RHS_imp((i-1)*3+1) = (Q_old(i,1) - Q(i,1))/dt + 1.0/dx*(F(i-1,1)-F(i,1))
            RHS_imp((i-1)*3+2) = (Q_old(i,2) - Q(i,2))/dt + 1.0/dx*(F(i-1,2)-F(i,2))
            RHS_imp((i-1)*3+3) = (Q_old(i,3) - Q(i,3))/dt + 1.0/dx*(F(i-1,3)-F(i,3))
   
            ! Berechnung LHS
            call calcTransMatIG(Q(i,:),R,kappa,Trans)
            !call calcEps(Z(i,:),x(nx)-x(1),dt,R,kappa,eps)
            eps = 1.0
            call calcPreConMatIG(Q(i,:),R,kappa,eps,PreCon,yang)
            a = sqrt(kappa*R*Z(i,3))
            dtau = 2.0*cfp*a*dt     &
                  /(abs(Z(i,2))*(1+eps)+sqrt(Z(i,2)**2*(1.0-eps)**2 + 4.0*eps*a*a))

            LHS_imp((i-1)*3+1:(i-1)*3+3,(i-1)*3+1:(i-1)*3+3) = PreCon/dtau + Trans/dt
   
            !call calcJacobiKonvPrimitivIG(0.5*(Z(i,:)+Z(i-1,:)),R,kappa,A_w)
            !call calcJacobiKonvPrimitivIG(0.5*(Z(i,:)+Z(i+1,:)),R,kappa,A_e)
            call calcJacobiKonvPrimitivIG(Z(i-1,:),R,kappa,A_w)
            call calcJacobiKonvPrimitivIG(Z(i+1,:),R,kappa,A_e)

            LHS_imp((i-1)*3+1:(i-1)*3+3,(i-2)*3+1:(i-2)*3+3) = -A_w/dx/2.0
            LHS_imp((i-1)*3+1:(i-1)*3+3,(i)*3+1:(i)*3+3) = A_e/dx/2.0
         end do ! Schleife über Zellen

         ! Lösen des Gleichungssystems
         call dgesv(3*nx,1,LHS_imp,3*nx,ipiv,RHS_imp,3*nx,info)
         !call solveGS(LHS_imp,RHS_imp,3*nx)

         ! Variablen update
         Z_imp = Z_imp + RHS_imp

         ! Übernhame der alten Werte
         do i=1,nx
            error(1) = error(1) + RHS_imp((i-1)*3+1)**2
            error(2) = error(2) + RHS_imp((i-1)*3+2)**2
            error(3) = error(3) + RHS_imp((i-1)*3+3)**2
            Z(i,1) = Z_imp((i-1)*3+1)
            Z(i,2) = Z_imp((i-1)*3+2)
            Z(i,3) = Z_imp((i-1)*3+3)
         end do

         error = sqrt(error)/nx
         errorP = error(1)
         if (inner==1) then
            targetP = errorP/1000.0
            write(*,*) inner,error
         end if

         do i=1,nx
            call calcConsVar(Q(i,:),Z(i,:),kappa,R)
         end do
      end do ! Innere Iteration

      write(*,*) inner,error
      call datOut(x,Q,Z,nx,outer)

   end do ! Äußere Iteration

   write(*,*) "Gesamt benötigte Iterationen: ", eps,cfp, totIter

!666 FORMAT(9(F10.2,1X))
end program
