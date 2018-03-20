program fv_explizit
   use konstants
   use ZGKonstants
   use fluxKonstants
   use limKonstants
   use dry_air
   use solver_vars
   use routines

   implicit none
   integer :: i

   ! Solverparameter
   cfl = 0.4
   no = 200
   nx = 501
   allocate(x(nx),Q(nx,3),Q_old(nx,3),Q_recon(2,nx-1,3),Z(nx,3),F(nx-1,3))

   ! 1D-Gitter
   call doGit(0d0,1d1,nx,x,dx)

   ! Anfangsbedingungen
   call startSod(Q,Z,kappa,R,nx)

   call datOut(x,Q,Z,nx,0)

   ! Solve
   do outer=1,no
      write(*,*) "Zeitschritt: ",outer
      ! Übernhame der alten Werte
      Q_old = Q
      dt = dx*cfl/sqrt(kappa*R*minval(Z(:,3))) ! Zeitschritt

      call reconQ(Q,Q_recon,nx,3,charm)
      call calcFlux(Q_recon,nx,kappa,R,F,AUSM)

      do i=2,nx-1
         RHS = dt/dx*(F(i-1,:)-F(i,:))
         Q(i,:) = Q_old(i,:) + RHS
         call calcPrimVar(Q(i,:),Z(i,:),kappa,R)
      end do
      ! Zero-Gradient Rand
      Q(1,:) = Q(2,:)
      Q(nx,:) = Q(nx-1,:)
      call calcPrimVar(Q(1,:),Z(1,:),kappa,R)
      call calcPrimVar(Q(nx,:),Z(nx,:),kappa,R)

      call datOut(x,Q,Z,nx,outer)
   end do
end program
