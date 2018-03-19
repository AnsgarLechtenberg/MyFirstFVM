program fv_preCon_explizit
   use konstants
   use ZGKonstants
   use fluxKonstants
   use limKonstants
   use preConKonstants
   use dry_air
   use solver_vars
   use routines

   implicit none

   real(kind=8) :: test

   ! Iterationsvariablen
   integer :: i,info,iterTot
   integer,dimension(3) :: ipiv

   ! Solverparameter
   no = 20
   ni = 50
   nx = 201
   cfl = 1.5
   cfp = 0.5

   allocate(x(nx),Q(nx,3),Q_old(nx,3),Q_recon(2,nx-1,3),Z(nx,3),F(nx-1,3))

   ! 1D-Gitter
   call doGit(0d0,1d1,nx,x,dx)

   ! Anfangsbedingungen
   call startSod(Q,Z,kappa,R,nx)

   call datOut(x,Q,Z,nx,0)

   ! Solve
   test = 0d0
   outer = 0
   iterTot = 0
   !do while (test<0.005) ! Äußere Iteraton
   do outer=1,no ! Äußere Iteraton
      test = test +dt
      !outer = outer +1
      Q_old = Q

      dt = dx*cfl/sqrt(kappa*R*minval(Z(:,3))) ! Zeitschritt

      inner = 0
      errorP = 1
      targetP = 0
      do while ((errorP>targetP) .and. (inner<ni))
         iterTot = iterTot + 1
         inner = inner + 1
         error = 0.0

         call reconQ(Q,Q_recon,nx,3)
         call calcFlux(Q_recon,nx,kappa,R,F,llf)

         do i=2,nx-1 ! Schleife über Zellen
            RHS = (Q_old(i,:) - Q(i,:))/dt +1.0/dx*(F(i-1,:)-F(i,:))
   
            call calcTransMatIG(Q(i,:),R,kappa,Trans)
            !call calcEps(Z(i,:),x(nx)-x(1),dt,R,kappa,eps)
            eps = 1d0
            call calcPreConMatIG(Q(i,:),R,kappa,eps,PreCon,yang)
            a = sqrt(kappa*R*Z(i,3))
            dtau = 2.0*cfp*a*dt/     &
                  (abs(Z(i,2))*(1+eps)+sqrt(Z(i,2)**2*(1.0-eps)**2 + 4.0*eps*a*a))
            LHS = PreCon/dtau + Trans/dt
   
            ! Lösen des Gleichungssytems
            call dgesv(3,1,LHS,3,ipiv,RHS,3,info)
            !call solve(LHS,RHS,3)
   
            ! Variablen update
            error = error + RHS*RHS
            Z(i,:) = Z(i,:) + RHS
         end do ! Schleife über Zellen
         ! Zero-Gradient Rand
         Z(1,:) = Z(2,:)
         Z(nx,:) = Z(nx-1,:)
         do i=1,nx
            call calcConsVar(Q(i,:),Z(i,:),kappa,R)
         end do
         error = sqrt(error)/nx
         errorP = error(1)
         if (inner==1) then
            targetP = errorP/1000.0
            write(*,*) inner,error
         end if
      end do ! Innere Iteration
         write(*,*) inner,error

      ! Prim in Conserv Variablen / Ausgabe
      do i=1,nx
         call calcConsVar(Q(i,:),Z(i,:),kappa,R)
      end do
      call datOut(x,Q,Z,nx,outer)
   end do ! Äußere Iteration
   write(*,*) cfp,test, iterTot
   !end do

!666 FORMAT((F5.3),(1x),6(F15.5,1X))
end program
