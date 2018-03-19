program test
   use routines
   implicit none

   integer,parameter :: dimen = 3
   integer :: i
   real(kind=8),dimension(dimen,dimen) :: LHS
   real(kind=8),dimension(dimen) :: RHS
   real(kind=8) :: start,ende

   LHS = 0.0
   do i = 1,dimen
      LHS(i,i) = 1.5
      RHS(i) = real(i)
   end do
   LHS(1,2) = 1.0
   LHS(2,3) = 1.0
   LHS(3,1) = 2.0

   write(*,*) LHS
   write(*,*)
   write(*,*) RHS

   call cpu_time(start)
   call solveGS(LHS,RHS,dimen)
   call cpu_time(ende)
   write(*,*) "Time = ",ende-start,"seconds."

   write(*,*) rhs
end program
