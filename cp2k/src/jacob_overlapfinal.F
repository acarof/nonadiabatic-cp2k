!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations
!   !
!   Copyright (C) 2000 - 2015  CP2K developers group
!   !
!-----------------------------------------------------------------------------!

!!!!!!!!!! CHANGE_AC !!!!!!!!!!!!!

! *****************************************************************************
!> \brief Interface for the adiabatic force calculations
!> \par History
!>      carof, 2015.12.11: creation, put only ffsh_aom subroutines
! *****************************************************************************


module jacob_overlapfinal
  USE kinds,                           ONLY: default_string_length,&
                                             dp
#include "./base/base_uses.f90"

      implicit none

      PUBLIC

      interface project_onto_unitvectors
        module procedure project_monomer_a
        module procedure project_monomer_b
      end interface

contains
      
      subroutine center_of_mass(xyz, mass, com)

      REAL(KIND=dp), intent(in) :: xyz(:,:)
      REAL(KIND=dp), intent(in) :: mass(:)
      REAL(KIND=dp), dimension(3), intent(out) :: com

      integer :: natm

      natm = size(xyz,2)

      call DGEMV('n', 3, natm, 1.0d0, xyz, 3, mass, 1, 0.0d0, com, 1)

      com = com / sum(mass)
      return
      end subroutine


      subroutine calc_rvecs2(xyz, com, rvecs)
!     calculates unnormalized connection vectors between
!     atomic coordinates and COM

      REAL(KIND=dp), dimension(3), intent(in) :: com
      REAL(KIND=dp), intent(in) :: xyz(:,:)

      REAL(KIND=dp), intent(out) :: rvecs(:,:)

      integer :: i

      do i=1,3
        rvecs(i,:) = xyz(i,:) - com(i)
      end do

      return
      end subroutine

      subroutine connect_list(xyz, atomlist, connectlist)

!     Calculates connectivity list for atoms identifying their
!     nearest neighbours for the sake of calculating local conjugation plane

      REAL(KIND=dp), intent(in) :: xyz(:,:) 
      integer, intent(in) :: atomlist(:)
      integer, intent(out) :: connectlist(:,:)
      integer :: i,j,k
      REAL(KIND=dp) :: distvec(3), DNRM2

      do i=1,size(atomlist)

         j=1

         connectlist(1,i)=atomlist(i)

         do k=1,size(xyz,2)

            distvec=xyz(:,atomlist(i)) - xyz(:,k)

!           Current nearest neighbour cutoff is 3.00 angs
            

            if (atomlist(i) .ne. k .and. DNRM2(3,distvec,1) < 3.00D0) then

                               j = j + 1
                connectlist(j,i) = k
            end if

         end do
   
         !sp2 carbons, N in pyrrole
         if(j .eq. 4) then

         !N in pyridine O in ethers
         elseif(j .eq. 3) then
             connectlist(j,i) = atomlist(i)
         !O in ketones
         elseif(j .eq. 2) then

            do k=1,size(xyz,2)

               distvec=xyz(:,connectlist(2,i)) - xyz(:,k)

!           Current nearest neighbour cutoff is 3.00 a. u.

            if (DNRM2(3,distvec,1) > 1.0D0 .and. atomlist(i) .ne. k &
&               .and. DNRM2(3,distvec,1) < 3.00D0) then

                               j = j + 1
                connectlist(j,i) = k
            end if
            end do

            connectlist(j,i) =  atomlist(i)

         else

             CPABORT("No such functional group ")
             !write(*,*) 'No such functional group' 
!             stop
            connectlist(j,i) =  atomlist(i)

         end if

      !print*, (connectlist(k,i), k=1,j)

      end do

      end subroutine

      subroutine connect_list2(xyz, atomlist, connectlist1, connectlist2)

!     Very similar to connect_list but also creates a connectlist2 array 
!     where the neighbours are at their true position rather than the one 
!     for which atomlist is needed. All those atoms that are not listed in
!     atomlist thus do not take part in the conjugation are denoted by 0 in the
!     first column but all connections to relevant atoms are listed. 
!     

      REAL(KIND=dp), intent(in) :: xyz(:,:)
      integer, intent(in) :: atomlist(:)
      integer, intent(out) :: connectlist1(:,:), connectlist2(:,:)
      integer :: i,j,k,l
      REAL(KIND=dp) :: distvec(3), DNRM2

      connectlist1=0
      connectlist2=0
!print*, atomlist
      do i=1,size(atomlist)

         j=1

         l=atomlist(i)

         connectlist1(1,l)=atomlist(i)
         connectlist2(1,i)=atomlist(i)

        ! print*,l

         do k=1,size(xyz,2)

            distvec=xyz(:,atomlist(i)) - xyz(:,k)

!print*, xyz(:,atomlist(i))
!print*, xyz(:,k)
            !print*, atomlist(i), k, DNRM2(3,distvec,1)


!           Current nearest neighbour cutoff is 3.00 a. u.


            if (atomlist(i) .ne. k .and. DNRM2(3,distvec,1) < 3.5D0) then

                               j = j + 1
            
                connectlist1(j,l) = k
                connectlist2(j,i) = k


            end if

         end do

         !print*, connectlist2(:,i), j

         !sp2 carbons, N in pyrrole
         if(j .eq. 4) then
!print*, "j = 4"

         !N in pyridine O in ethers
         elseif(j .eq. 3) then
!print*, "j = 3"
             connectlist1(j+1,l) = atomlist(i)
             connectlist2(j+1,i) = atomlist(i)
         !O in ketones
         elseif(j .eq. 2) then
!print*, "j = 2"
            do k=1,size(xyz,2)

               distvec=xyz(:,connectlist1(2,l)) - xyz(:,k)

!              Current nearest neighbour cutoff is 3.00 a. u.

               if (DNRM2(3,distvec,1) > 1.0D0 .and. atomlist(i) .ne. k &
&                  .and. DNRM2(3,distvec,1) < 3.5D0) then

                               j = j + 1
                connectlist1(j,l) = k
                connectlist2(j,i) = k
            end if
            end do

            connectlist1(j+1,l) =  atomlist(i)
            connectlist2(j+1,i) =  atomlist(i)
         else

             CPABORT("No such functional group")
             !write(*,*) 'No such functional group'
             !print*, "No such functional group"
             !print *, j
            !connectlist1(j+1,l) =  atomlist(i)
            !connectlist2(j+1,i) =  atomlist(i)
            ! stop

         end if

      end do

      do i=1,size(connectlist2, 2)

         do j=2,4

            if(connectlist1(1,connectlist2(j,i)) .eq. 0) then

               k=2

               do while (connectlist1(k,connectlist2(j,i)) > 0)

                  k=k+1

               end do

               connectlist1(k,connectlist2(j,i))=connectlist2(1,i)         

            end if        

         end do

         !print*, connectlist2(:,i)

      end do

      end subroutine


      subroutine calc_rvecs(xyz, rvecs, connectlist)
!     Calculates the p_pi direction in the general case for conjugated atoms
!     using the plane of the 3 connecting atoms

      REAL(KIND=dp), intent(in) :: xyz(:,:)
      integer, intent(in) :: connectlist(:,:)
      REAL(KIND=dp), intent(out) :: rvecs(:,:)
      integer :: i, k
      REAL(KIND=dp) :: neighbourlist(3,3)
      REAL(KIND=dp) :: veca(3), vecb(3), DNRM2
      
!print*, "  "
!print*, "xyz", xyz(:,1)
!print*, "connectlist", connectlist

      rvecs=0

      do i=1,size(connectlist,2)

         do k = 2,size(connectlist,1)
!print*, i,k,connectlist(k,i)
!print*, xyz(:,connectlist(k,i))
              neighbourlist(:,k-1) = xyz(:,connectlist(k,i))

         end do

!print*, neighbourlist(:,1)
!print*, neighbourlist(:,2)
!print*, neighbourlist(:,3)

         veca=neighbourlist(:,1)-neighbourlist(:,3)
         vecb=neighbourlist(:,2)-neighbourlist(:,3)

         rvecs(1,connectlist(1,i)) = veca(2)*vecb(3)-veca(3)*vecb(2)
         rvecs(2,connectlist(1,i)) = veca(3)*vecb(1)-veca(1)*vecb(3)
         rvecs(3,connectlist(1,i)) = veca(1)*vecb(2)-veca(2)*vecb(1)
!print*, "veca", veca
!print*, "vecb", vecb
         

!print*, 290,  "calc_recvs", i, rvecs(:,connectlist(1,i))
         rvecs(:,connectlist(1,i)) = rvecs(:,connectlist(1,i))/&
&        DNRM2(3 , rvecs(:,connectlist(1,i)) , 1)
!print*, "DNRM2", i, DNRM2(3 , rvecs(:,connectlist(1,i)) , 1) 
!print*, "calc_recvs", i, rvecs(:,connectlist(1,i))

      end do


      end subroutine 
      

      subroutine calc_unit_vecs(xyz, xyzb, dist, unitvecs)
!      calculates unit vectors with ez along the connection 
!      line between two atoms of monomer A and B for 
!      a single atom in A

      REAL(KIND=dp), intent(in) :: xyz(:)
      REAL(KIND=dp), intent(in) :: xyzb(:,:)

      REAL(KIND=dp), intent(out) :: dist(:)
      REAL(KIND=dp), intent(out) :: unitvecs(:,:,:)

      REAL(KIND=dp) :: DNRM2

      integer :: i,j

      integer :: midx(1)

      !calculate ez
      do i=1,size(xyzb,2)

        unitvecs(:,i,3) = xyzb(:,i) - xyz
        dist(i) = DNRM2(3, unitvecs(1,i,3), 1)

        if (dist(i) .ge. epsilon(0.d0)) then
          unitvecs(:,i,3) = unitvecs(:,i,3) / dist(i)
        else
          unitvecs(:,i,3) = 0.d0
        end if
      end do

      !calculate ey

      do i=1,size(xyzb,2)
        if (dist(i) .gt. epsilon(0.d0)) then
          midx = maxloc(abs(unitvecs(:,i,3)))
          !check if ez conincides with ex, ey or ez of the fixed frame coordinate system 
          if ((abs(unitvecs((mod(midx(1),3)+1),i,3)) .le. epsilon(0.0d0)) &
&             .and. (abs(unitvecs((mod((midx(1)+1),3)+1),i,3)) .le. epsilon(0.d0))) then
            unitvecs(midx,i,2) = 0.0d0
            unitvecs((mod(midx,3)+1),i,2) = 0.0d0
            unitvecs((mod((midx+1),3)+1),i,2) = 1.0d0
          else
            do j=1,3
              if (midx(1) .ne. j) then
                unitvecs(j,i,2) = 1.0d0
              else
                unitvecs(j,i,2) = - (unitvecs((mod(j,3)+1),i,3) + unitvecs((mod((j+1),3)+1),i,3)) &
                                   & / unitvecs(j,i,3)
              end if
            end do
          end if
          unitvecs(:,i,2) = unitvecs(:,i,2) / DNRM2(3, unitvecs(1,i,2),1)
        else
          unitvecs(:,i,2) = 0.d0
        end if
      end do
      
      !calculate ex
      !cross product between ey and ez

      do i=1,size(xyzb,2)
        if (dist(i) .gt. epsilon(0.d0)) then
          unitvecs(1,i,1) = unitvecs(2,i,2) * unitvecs(3,i,3) &
                          & - unitvecs(3,i,2) * unitvecs(2,i,3)
          unitvecs(2,i,1) = unitvecs(3,i,2) * unitvecs(1,i,3) &
                          & - unitvecs(1,i,2) * unitvecs(3,i,3)
          unitvecs(3,i,1) = unitvecs(1,i,2) * unitvecs(2,i,3) &
                          & - unitvecs(2,i,2) * unitvecs(1,i,3)

          unitvecs(:,i,1) = unitvecs(:,i,1) / DNRM2(3, unitvecs(1,i,1),1)
        else
          unitvecs(:,i,1) = 0.d0
        end if

      end do

     
      return
      end subroutine


      subroutine project_monomer_b(unitvecs, rvecs, expcoeff)
!     Does the projection for each atom in monomer B

      REAL(KIND=dp), dimension(:,:), intent(in) :: rvecs
      REAL(KIND=dp), dimension(:,:,:), intent(in) :: unitvecs

      REAL(KIND=dp), dimension(:,:), intent(out) :: expcoeff

      REAL(KIND=dp), dimension(:,:), allocatable :: prvecs
      !real(8) :: DNRM2

      integer :: stat_err

      integer :: i,j,k

      allocate(prvecs(size(unitvecs,2),3), stat=stat_err)
      if (stat_err .ne. 0) CPABORT("Problems in allocation")

      prvecs = 0.0d0

      do k=1,3
        do j=1,size(unitvecs,2)
          do i=1,3
            prvecs(j,k) = prvecs(j,k) + unitvecs(i,j,k) * rvecs(i,j)
          end do
        end do
      end do

      expcoeff = transpose(prvecs)

!     Only needed if calc_rvecs is used

      !Normalize
      !do i=1,size(unitvecs,2)
      !  print*, expcoeff(1,i), expcoeff(2,i), expcoeff(3,i)
      !  expcoeff(:,i) = expcoeff(:,i) / DNRM2(3,expcoeff(1,i),1)
      !  print*, expcoeff(1,i), expcoeff(2,i), expcoeff(3,i)
      !end do

      deallocate(prvecs)
      return
      end subroutine project_monomer_b


      subroutine project_monomer_a(unitvecs, rvec, expcoeff)
!     Does the projection for one atom in monomer A

      REAL(KIND=dp), dimension(:), intent(in) :: rvec
      REAL(KIND=dp), dimension(:,:,:), intent(in) :: unitvecs

      REAL(KIND=dp), dimension(:,:), intent(out) :: expcoeff

 
      REAL(KIND=dp), dimension(:,:), allocatable :: prvecs
      !real(8) :: DNRM2, norm

      integer :: stat_err

      integer :: i,j,k
!print*, "   "
!print*, 434
!print*, unitvecs
!print*, " "
!print*, rvec
!print*, " " 

      allocate(prvecs(size(unitvecs,2),3), stat=stat_err)
      !if (stat_err .ne. 0) stop
      if (stat_err .ne. 0) CPABORT(" ")

      prvecs = 0.d0


      do k=1,3
        do j=1,size(unitvecs,2)
          do i=1,3
            prvecs(j,k) = prvecs(j,k) + unitvecs(i,j,k) * rvec(i)
          end do
        end do
      end do

      expcoeff = transpose(prvecs)

!      Only needed if calc rvecs 2 is used

      !Normalize
      !do i=1,size(unitvecs,2)
      !  norm = DNRM2(3,expcoeff(1,i),1)
      !  if (norm .gt. epsilon(0.d0)) then
      !    expcoeff(:,i) = expcoeff(:,i) / norm
      !  else
      !    expcoeff(:,i) = 0.d0 
      !  end if
      !end do

      deallocate(prvecs)
      return
      end subroutine project_monomer_a

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                                   !  
!     The following set of subroutines are to calculate the overlap between atomic  !
!     orbitals according to Mulliken. Currently, the possibility to calculate the   !
!     the overlap including the s orbitals as well is implemented but disabled.     !
!                                                                                   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



      subroutine calc_overlap_p(expcoeffa, expcoeffb, mua, musb, dist, sarray)

!     Subroutine calculating the overlap between p orbitals only currently
!     currently this is the only one called.
 
      REAL(KIND=dp), intent(in) :: expcoeffa(:,:)
      REAL(KIND=dp), intent(in) :: expcoeffb(:,:)
      REAL(KIND=dp), intent(in) :: mua, musb(:)
      REAL(KIND=dp), intent(in) :: dist(:)
 
      REAL(KIND=dp), intent(out) :: sarray(:)
 
      REAL(KIND=dp), allocatable :: c(:,:)
      REAL(KIND=dp) :: p, t
 
      integer :: stat_err
      integer :: i
      allocate(c(3,size(expcoeffa,2)), stat=stat_err)
      !if (stat_err .ne. 0) stop
      if (stat_err .ne. 0) CPABORT(" ")
!print*, 502, expcoeffa
!print*, "   " 
!print*, 502, expcoeffb 
      c = expcoeffa * expcoeffb
      do i=1,size(expcoeffa,2)
        if (dist(i) .gt. abs(epsilon(0.d0))) then
          p = (mua + musb(i))/2.0D0*dist(i)
          t = (mua - musb(i))/(mua + musb(i)) !abs???


          !sarray(i) = (c(1,i) + c(2,i)) * S_pi(dist(i)) + c(3,i) * S_sigma(dist(i))
!print*, "    "
!print*, 512, "sarray, i", i, sarray(i)
          sarray(i) = (c(1,i) + c(2,i)) * S_ppi(p,t) + c(3,i) * S_psigma(p,t)
!print*, 512, "c, i", i, c(1,i), c(2,i), c(3,i)
!print*, 512, "c, p,t", p,t, S_ppi(p,t), S_psigma(p,t)

!print*, 512, "sarray, i", i, sarray(i)
!print*, "    "
        else
         sarray(i) = 1.d0
        end if
      end do
      deallocate(c)
      return
      end subroutine calc_overlap_p


      subroutine calc_overlap_s(mua, musb, dist, sarray)

!     Calculates overlap between s orbitals currently not called

      REAL(KIND=dp), intent(in) :: mua, musb(:)
      REAL(KIND=dp), intent(in) :: dist(:)

      REAL(KIND=dp), intent(out) :: sarray(:)


      REAL(KIND=dp) :: p,t
      integer :: i

      do i=1,size(musb)
        if (dist(i) .gt. abs(epsilon(0.d0))) then
           p = (mua + musb(i))/2.0D0*dist(i)
           t = (mua - musb(i))/(mua + musb(i)) !abs???
           sarray(i) = S_ss(p,t) 
        else
           sarray(i) = 1.d0 
        end if
      end do
      return
      end subroutine calc_overlap_s

      subroutine calc_overlap_sp(expcoeff, mua, musb, dist, sarray)

!     Calculates overlap between s and p orbitals currently not called

      REAL(KIND=dp), intent(in) :: mua, musb(:)
      REAL(KIND=dp), intent(in) :: dist(:), expcoeff(:)

      REAL(KIND=dp), intent(out) :: sarray(:)

      REAL(KIND=dp) :: p,t
      integer :: i

      do i=1,size(musb)
        if (dist(i) .gt. abs(epsilon(0.d0))) then
           p = (abs(mua) + abs(musb(i)))/2.0D0*dist(i)
           t = (mua - musb(i))/(abs(mua) + abs(musb(i)))

           sarray(i) = expcoeff(i)*S_sp(p,t) 
        else
           sarray(i) = 0.d0
        end if
      end do


      return
      end subroutine calc_overlap_sp

!     Series of subroutines to calculate individual atomic overlaps
!     the if statement is used to distinguish between homo- and heteroatomic 
!     overlaps. The former one needing much simpler formulae is much faster.

      function S_ss(p,t)
!     Caculates s overlap integral
      REAL(KIND=dp) :: S_ss
      REAL(KIND=dp), intent(in) :: p, t

      if (abs(t) .lt. epsilon(0.D0)) then

          S_ss = dexp(-p)/9.0D0 *(9.0D0 + 9.0D0*p + 4.0d0*p*p + 0.20D0*p*p*p)

      else

          S_ss = 1/48.0D0 * p**5.0D0 *(1.0D0 - t**2)**(5.0D0/2.0D0)*&
&                (A(4,p)*B(0,p,t)-2.0D0*A(2,p)*B(2,p,t)-A(0,p)*B(4,p,t))

      end if

      end function

      function S_sp(p,t)
!     Caculates s-p(sigma) overlap integral
      REAL(KIND=dp) :: S_sp

      REAL(KIND=dp), intent(in) :: p, t

      REAL(KIND=dp) :: m

      m = 1.0

      if (abs(t) .lt. epsilon(0.D0)) then
      
          S_sp = sqrt(3.0D0)*dexp(-p)*(15.0D0*p + 15.0D0*p*p + 7.0D0*p*p*p + 2.0D0*p*p*p*p)/90.0D0 

      else

          S_sp =  1/60.0D0 * sqrt(3.0D0)*p*p*p*p*p*(1 - t*t)**(2.50D0)*&
&                 (A(3,p)*(B(0,p,t)-B(2,p,t)) + A(2,p)*(B(4,p,t)-B(2,p,t))&
&                 + B(1,p,t)*(A(2,p)-A(4,p)) + B(3,p,t)*(A(2,p)-A(0,p)))

      end if

      end function


      function S_psigma(p,t)
!     Calculates p_sigma-p_sigma overlap integral 

      REAL(KIND=dp) :: S_psigma
      REAL(KIND=dp), intent(in) :: p, t
 
      if (abs(t) .lt. epsilon(0.D0)) then

          S_psigma = - dexp(-p)/15.0D0 * (p*p*p*p + 2.0D0*p*p*p - 3.0D0*p*p - 15.0D0*p - 15.0D0)    

      else
 
          S_psigma =  1/16.0D0 * p**5 *(1 - t**2)**(5.0D0/2.0D0)*&
&                 (B(2,p,t)*(A(0,p)+A(4,p)) - A(2,p)*(B(0,p,t)+B(4,p,t)))

      end if
      end function

      function S_ppi(p,t)
!     Calculates p_pi-p_pi overlap integral

      REAL(KIND=dp) :: S_ppi
      REAL(KIND=dp), intent(in) :: p, t
 

      if (abs(t) .lt. epsilon(0.D0)) then

         S_ppi = dexp(-p) / 15.0D0 * (p*p*p + 6.0D0 * p*p + 15.0D0 * p + 15.0D0)

      else

         S_ppi =  1/32.0D0 * p*p*p*p*p*(1 - (t*t))**(5.0D0/2.0D0)*&
&                 (A(4,p)*(B(0,p,t)-B(2,p,t))&
&                 +A(2,p)*(B(4,p,t)-B(0,p,t))&
&                 +A(0,p)*(B(2,p,t)-B(4,p,t)))

      end if

      end function

!     Subroutines factorial, A and B are only needed for hetero atompair calculations or
!     s-p overlaps.

      function factorial(k)

      integer :: factorial, i 
      integer, intent(in) :: k

      factorial=1

      if (k .ge. 1) then

         do i=1,k

            factorial=factorial*i

         end do

      elseif (k < 0) then
          
!        print*, 'Negative number argumented to factorial!'
       ! stop
         CPABORT("Negative number argumented to factorial!")

      end if

      end function  

      function A(k,p)

      REAL(KIND=dp) :: A
      REAL(KIND=dp), intent(in) :: p
      integer, intent(in) :: k 
      integer :: mu

      A = 0

      do mu=1,(k+1)

          A = A+dexp(-1*p)*factorial(k)/factorial(k-mu+1)/(p**mu)

      end do

      end function  

      function B(k,p,t)

      REAL(KIND=dp) :: B
      REAL(KIND=dp), intent(in) :: p, t
      integer, intent(in) :: k
      integer :: mu

      B=0 

      do mu=1,(k+1)

         B = B -dexp(-1*p*t)*factorial(k)/factorial(k-mu+1)/((p*t)**mu)&
&               -dexp(p*t)*factorial(k)/factorial(k-mu+1)/((p*t)**mu)*(-1)**(k-mu)

      end do


      end function

!     Original Felix atomic overlap progrrams currently not used

      function S_pi(d)
      !Calculates p-pi overlap integral
      REAL(KIND=dp) :: S_pi
      REAL(KIND=dp), intent(in) :: d

      REAL(KIND=dp) :: mu
      REAL(KIND=dp) :: p,b

      !include 'slater_decayconstants.h'

      !optimised pi decay constant in a.u.
      mu = 1.00D0 !533 !3125 !0.89 !0.90636d0 !0.856016d0 
      ! optimised scaling factor
      b =  1.0D0 !/1.304 !1 !0.136225d0*3.231  !change of factor by Fruzsi fitted to C60

      p = mu * d
      S_pi = b * dexp(-p) / 15.0D0 * (p*p*p + 6.0D0 * p*p + 15.0D0 * p + 15.0D0) 
      end function

      
      function S_sigma(d)
      !Calculates p-sigma overlap integral
      REAL(KIND=dp) :: S_sigma
      REAL(KIND=dp), intent(in) :: d

      REAL(KIND=dp) :: mu
      REAL(KIND=dp) :: p,b

!      !include 'slater_decayconstants.h'

!      !optimised sigma decay constant in a.u.
      mu = 1.00D0 !533 !3125 !0.89 !0.90636d0
!      ! optimised scaling factor
      b =  1.0D0 !/1.304 !0.136225d0*3.231 !0.23615d0*3.231  !change of factor by Fruzsi fitted to C60

      p = mu * d
      S_sigma = - b/15.0D0 * exp(-p) * (p*p*p*p + 2.0D0*p*p*p - 3.0D0*p*p - 15.0D0*p - 15.0D0) 
      end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine calc_sab(prexyza, connlista, precoeffap, &
&                         sab, prexyzb, connlistb, precoeffbp)

!     Original S_ab calculation program can calculate an pi_conjugated molecule
!
!     - Slater coefficients for carbon can be read in from external data file
!       by enabling lines 779-786
!     - Inclusion of s orbitals can be achieved by enabling lines 876-883,916- 920, 
!       930-944, 949, 974-981, 1011-1017 and 1036-1049


      REAL(KIND=dp), dimension(:,:),           intent(in) ::    prexyza
      REAL(KIND=dp), dimension(:),             intent(inout) :: precoeffap
      integer, dimension(:,:),           intent(in) ::    connlista
      REAL(KIND=dp), dimension(:,:), optional, intent(in) ::    prexyzb
      REAL(KIND=dp), dimension(:),   optional, intent(in) ::    precoeffbp
      integer, dimension(:,:), optional, intent(in) ::    connlistb
      REAL(KIND=dp),                 optional, intent(out) ::   sab

      REAL(KIND=dp), dimension(:), allocatable :: stc
      REAL(KIND=dp), dimension(:), allocatable :: dist
      REAL(KIND=dp), dimension(:,:), allocatable :: rvecsa, rvecsb
      REAL(KIND=dp), dimension(:,:), allocatable :: prervecsa, prervecsb 
      REAL(KIND=dp), dimension(:), allocatable :: coeffap, coeffbp
      REAL(KIND=dp), dimension(:,:), allocatable :: sarrayp 
      REAL(KIND=dp), dimension(:,:), allocatable :: expcoeffa, expcoeffb
      REAL(KIND=dp), dimension(:,:), allocatable :: xyza, xyzb
      REAL(KIND=dp), dimension(:), allocatable ::  musbp
      REAL(KIND=dp) :: muap 
      REAL(KIND=dp), dimension(:,:,:), allocatable :: unitvecs
      REAL(KIND=dp), dimension(8,2) :: slaters
 
      integer :: stat_err !, ioerr, ioread
      !integer :: i,j,k,z
      integer :: i,j
 
!      REAL(KIND=dp) :: com(3)
      REAL(KIND=dp) :: c
      REAL(KIND=dp) :: DDOT

      slaters(:,1)=(/0.0D0, 0.0D0, 3.0D0, 0.9560D0, 1.2881D0, 1.100D0, 1.9237D0, 2.2458D0/)
      slaters(:,2)=(/0.0D0, 0.0D0, 0.0D0, 0.0000D0, 1.2107D0, 1.000D0, 1.6000D0, 2.2266D0/)

      !open(unit=69, file='slats', status='old', action='read', iostat=ioerr)
      !if (ioerr .ne. 0) then
      !  print*, 'Could not open file slats'
      !  stop
      !else
      !   read(69,*, iostat=ioread) slaters(6,1), slaters(6,2)
      !end if
      !close(69)

      allocate(prervecsa(3,size(prexyza,2)),stat=stat_err)
      !if (stat_err .ne. 0) stop
      if (stat_err .ne. 0) CPABORT(" ")

      allocate(xyza(4,size(connlista,2)),stat=stat_err)
      !if (stat_err .ne. 0) stop
      if (stat_err .ne. 0) CPABORT(" ")

      allocate(coeffap(size(connlista,2)),stat=stat_err)
      !if (stat_err .ne. 0) stop
      if (stat_err .ne. 0) CPABORT(" ")

      !allocate(coeffas(size(connlista,2)),stat=stat_err)
      !if (stat_err .ne. 0) stop

      allocate(rvecsa(3,size(connlista,2)),stat=stat_err)
      !if (stat_err .ne. 0) stop
      if (stat_err .ne. 0) CPABORT(" ")

      allocate(stc(size(xyza,2)),stat=stat_err)
      !if (stat_err .ne. 0) stop
      if (stat_err .ne. 0) CPABORT(" ")


      call calc_rvecs(prexyza(2:4,:), prervecsa, connlista)

      do i=1,size(connlista,2)

           xyza(:,i)=prexyza(:,connlista(1,i))
           coeffap(i)=precoeffap(connlista(1,i))
           !coeffas(i)=precoeffas(connlista(1,i))
           rvecsa(:,i)=prervecsa(:,connlista(1,i))

      end do
!print*, "calc_sab", coeffap

      deallocate(prervecsa)
 
      ! These variables are recaculated for each atom in A

      ! Inherent for B atom

      if (present(prexyzb) .and. present(precoeffbp) .and. &
&         present(sab) .and. present(connlistb)) then

         allocate(unitvecs(3,size(connlistb,2),3),stat=stat_err)
         !if (stat_err .ne. 0) stop
         if (stat_err .ne. 0) CPABORT(" ")

         allocate(dist(size(connlistb,2)),stat=stat_err)
         !if (stat_err .ne. 0) stop
         if (stat_err .ne. 0) CPABORT(" ")

         allocate(prervecsb(3,size(prexyzb,2)),stat=stat_err)
         !if (stat_err .ne. 0) stop
         if (stat_err .ne. 0) CPABORT(" ")
 
         allocate(xyzb(4, size(connlistb,2)),stat=stat_err)
         !if (stat_err .ne. 0) stop
         if (stat_err .ne. 0) CPABORT(" ")

         allocate(expcoeffa(3,size(xyzb,2)),stat=stat_err)
         !if (stat_err .ne. 0) stop
         if (stat_err .ne. 0) CPABORT(" ")

         allocate(coeffbp(size(xyzb,2)),stat=stat_err)
         !if (stat_err .ne. 0) stop
         if (stat_err .ne. 0) CPABORT(" ")
  
         !allocate(coeffbs(size(xyzb,2)),stat=stat_err)
         !if (stat_err .ne. 0) stop

         allocate(expcoeffb(3,size(xyzb,2)),stat=stat_err)
         !if (stat_err .ne. 0) stop
         if (stat_err .ne. 0) CPABORT(" ")

         allocate(rvecsb(3,size(connlistb,2)),stat=stat_err)
         !if (stat_err .ne. 0) stop
         if (stat_err .ne. 0) CPABORT(" ")

         call calc_rvecs(prexyzb(2:4,:), prervecsb, connlistb)

         do i=1,size(connlistb,2)

             xyzb(:,i)=prexyzb(:,connlistb(1,i))
             coeffbp(i)=precoeffbp(connlistb(1,i))
             !coeffbs(i)=precoeffbs(connlistb(1,i))
             rvecsb(:,i)=prervecsb(:,connlistb(1,i))

        end do  

  
        deallocate(prervecsb)
       
        allocate(musbp(size(xyzb,2)),stat=stat_err)
        !if (stat_err .ne. 0) stop
        if (stat_err .ne. 0) CPABORT(" ")

        allocate(sarrayp(size(xyzb,2),size(xyza,2)), stat=stat_err)
        !if (stat_err .ne. 0) stop
        if (stat_err .ne. 0) CPABORT(" ")

        !allocate(sarrays(size(xyza,2),size(xyzb,2)), stat=stat_err)
        !if (stat_err .ne. 0) stop

        !allocate(sarraysp(size(xyza,2),size(xyzb,2)), stat=stat_err)
        !if (stat_err .ne. 0) stop

        !allocate(sarrayps(size(xyza,2),size(xyzb,2)), stat=stat_err)
        !if (stat_err .ne. 0) stop

        do i=1,size(xyzb,2)

         !musbs(i)=slaters(int(xyzb(1,i)),1)
         musbp(i)=slaters(int(xyzb(1,i)),2)

        end do


        do i=1,size(xyza,2)

          muap=slaters(int(xyza(1,i)),2)
         
          !calculate unitvectors for each atom in monomer A
          !and distance array
          call calc_unit_vecs(xyza(2:4,i), xyzb(2:4,:), dist, unitvecs)

          !Do the expansion for monomer A
          call project_onto_unitvectors(unitvecs, rvecsa(:,i), expcoeffa)

          !Do the expansion for monomer B

          call project_onto_unitvectors(unitvecs, rvecsb, expcoeffb)

          !Calculate the overlaps

          call calc_overlap_p(expcoeffa, expcoeffb, muap, musbp, dist, sarrayp(:,i))

          !call calc_overlap_s(muas, musbs, dist, sarrays(:,i))

          !call calc_overlap_sp(expcoeffb(3,:), muas, musbp, dist, sarraysp(:,i))

          !call calc_overlap_sp(expcoeffa(3,:), muap, musbs, dist, sarrayps(:,i))

        end do
        

        call DGEMV('t', size(sarrayp,1), size(sarrayp,2), 1.0d0, sarrayp,&
                   & size(sarrayp,1), coeffbp, 1, 0.d0, stc, 1)

        sab = DDOT(size(coeffap,1),coeffap,1,stc,1)
!print*, "calc_sab", coeffap

        !print*, size(sarrayp,1), size(sarrayp,2)

        !call DGEMV('t', size(sarrays,2), size(sarrays,2), 1.0d0, sarrays,&
        !            & size(sarrays,2), coeffbs, 1, 0.d0, stc, 1)

        !sab = sab+DDOT(size(coeffas,1),coeffas,1,stc,1)


        !call DGEMV('t', size(sarraysp,1), size(sarraysp,2), 1.0d0, sarraysp,&
        !           & size(sarrayp,2), coeffbp, 1, 0.d0, stc, 1)

        !sab = sab + DDOT(size(coeffas,1),coeffas,1,stc,1)

        !call DGEMV('t', size(sarrayps,2), size(sarrayps,2), 1.0d0, sarrayps,&
        !             & size(sarrayps,2), coeffbs, 1, 0.d0, stc, 1)

        !sab = sab+DDOT(size(coeffap,1),coeffap,1,stc,1)


        deallocate(rvecsb,  expcoeffa, expcoeffb, sarrayp, dist, unitvecs, coeffbp, musbp)
        !deallocate(sarrays, sarrayps, sarraysp, coeffbs, musbs)


      elseif ((present(prexyzb) .and. present(precoeffbp) .and. &
&          present(sab) .and. present(connlistb)) .eqv. .false.) then

        !Normalize

        allocate(unitvecs(3,size(xyza,2),3),stat=stat_err)
        !if (stat_err .ne. 0) stop
        if (stat_err .ne. 0) CPABORT(" ")

        allocate(dist(size(xyza,2)),stat=stat_err)
        !if (stat_err .ne. 0) stop
        if (stat_err .ne. 0) CPABORT(" ")

        allocate(musbp(size(xyza,2)),stat=stat_err)
        !if (stat_err .ne. 0) stop
        if (stat_err .ne. 0) CPABORT(" ")

        allocate(sarrayp(size(xyza,2),size(xyza,2)), stat=stat_err)
        !if (stat_err .ne. 0) stop
        if (stat_err .ne. 0) CPABORT(" ")

        allocate(expcoeffa(3,size(xyza,2)),stat=stat_err)
        !if (stat_err .ne. 0) stop
        if (stat_err .ne. 0) CPABORT(" ")

        allocate(expcoeffb(3,size(xyza,2)),stat=stat_err)
        !if (stat_err .ne. 0) stop
        if (stat_err .ne. 0) CPABORT(" ")


        sarrayp=0.0D0

        do i=1,size(xyza,2)

         musbp(i)=slaters(int(xyza(1,i)),2)

        end do

          do i=2,size(xyza,2)


             muap=slaters(int(xyza(1,i)),2)

!            Only calculates lower triangle of the overlap matrix in the normalisation step
!            then mirrors it (all real) over the diagonal 
 
             j=i-1        

             dist=0.0D0

             unitvecs=0.0D0

             expcoeffa=0.0D0

             expcoeffb=0.0D0

             call calc_unit_vecs(xyza(2:4,i), xyza(2:4,1:j), dist(1:j), unitvecs(:,1:j,:))    

             !print*, 'unitvecs', dist(2), (unitvecs(k, 2, 3), k=1,3)

             call project_onto_unitvectors(unitvecs(:,1:j,:), rvecsa(:,i), expcoeffa(:,1:j))

             call project_onto_unitvectors(unitvecs(:,1:j,:), rvecsa(:,1:j), expcoeffb(:,1:j))

             !print*, 'expcoeffa' 

             !do z=1,size(expcoeffa, 2)

             !      print*, (expcoeffa(k,z), k=1,3)
           
             !end do

             !print*, 'expcoeffa'


             !do z=1,size(expcoeffb, 2)

             !      print*, (expcoeffb(k,z), k=1,3)

             !end do
             
             !print*, 'dist'

             !do z=1,size(dist)

             !      print*, dist(z)

             !end do


             !do z=1,size(sarrayp,2)

             !    print*, (sarrayp(k,z), k=1,size(sarrayp, 1))

             !end do

             !print*, 'mua', muap

             !print*, 'mub', musbp


             !print*, 'Some info', dist(25),  (expcoeffb(k,25), k=1,3)

             !print*,'inp', (expcoeffb(k,2), k=1,3), (expcoeffa(k,2), k=1,3), dist(2)

             call calc_overlap_p(expcoeffa(:,1:j), expcoeffb(:,1:j), muap, musbp(1:j), dist(1:j), sarrayp(1:j,i))
          
             !print*, 'dist, sarray', dist(6),  sarrayp(6,i)

             !print*, 'dist, sarray', i, dist(6),  sarrayp(6,i)

             !call calc_overlap_p(expcoeffa, expcoeffb, muap, musbp, dist, sarrayp(:,i))
             
             !call calc_overlap_s(muas, musbs, dist, sarrays(:,i))

             !call calc_overlap_sp(expcoeffb(3,:), muas, musbp, dist, sarraysp(:,i))

             !call calc_overlap_sp(expcoeffa(3,:), muap, musbs, dist, sarrayps(:,i))
   

          end do

          !print*, lement 6, 7', sarrayp(2,3)

          sarrayp=sarrayp+transpose(sarrayp)         

          do i=1,size(xyza,2)

              sarrayp(i,i)=1.0D0

          end do

          !do i=1,size(sarrayp,2)

          ! print*, (sarrayp(k,i), k=1,size(sarrayp, 1))

          !end do

          !print*, 'line 8', (sarrayp(k,8), k=1,size(sarrayp, 1))

          call DGEMV('t', size(xyza,2), size(xyza,2), 1.0d0, sarrayp,&
                     & size(xyza,2), coeffap, 1, 0.d0, stc, 1)
!print*, "calc_sab", coeffap

          c = DDOT(size(xyza,2),coeffap,1,stc,1)
!print*, "calc_sab", coeffap

          !call DGEMV('t', size(xyza,2), size(xyza,2), 1.0d0, sarrays,&
          !           & size(xyza,2), coeffas, 1, 0.d0, stc, 1)

          !c = c+DDOT(size(xyza,2),coeffas,1,stc,1)

          !call DGEMV('t', size(sarraysp,1), size(sarraysp,2), 1.0d0, sarraysp,&
          !         & size(sarraysp,2), coeffap, 1, 0.d0, stc, 1)

          !c = c+DDOT(size(coeffas,1),coeffas,1,stc,1)

          !call DGEMV('t', size(sarrayps,1), size(sarrayps,2), 1.0d0, sarrayps,&
          !           & size(sarrayps,2), coeffas, 1, 0.d0, stc, 1)

          !c = c+DDOT(size(coeffap,1),coeffap,1,stc,1)

          !print*, 'Sarray element 8, 9', sarrayp(2,3)           

          !write(*,*) 'This is c ',sqrt(c)

          coeffap = coeffap / sqrt(c) 
!print*, "calc_sab", coeffap
          !coeffas = coeffas / sqrt(c)

          do i=1,size(connlista,2)

              precoeffap(connlista(1,i)) = coeffap(i)
 
          end do
         
         deallocate(sarrayp, expcoeffa, expcoeffb, dist, unitvecs, musbp) 
         !deallocate(sarrayps, sarraysp, sarrays)


      else

!         print*, 'Wrong argument passing to calc_sab' 

      end if

     deallocate(rvecsa, xyza, coeffap, stc)


      return
      end subroutine calc_sab


      subroutine calc_sab_turbo(prexyza, atomlista, connlistal, precoeffap, prervecsa, sarray, linecol,&
&                         sab, prexyzb, atomlistb, connlistbl, precoeffbp, prervecsb)

!     An alternative version for calc_sab. This version requires a variable called linecol
!     which tells the sytem in the normalisation case which columns and rows are to be 
!     recalculated. In the overlap case it only changes given rows.
!
!     The options to call external 'slats' file and the inclusion of s orbitals are omitted here.
!     The former feature is needless since by this point one should have finalised mu value. The 
!     latter one was omitted as we found the s contribution insignficant and computatiionally demanding !
!
!     This code only works with the longer version of the connectivity list where there is an index
!     to index correspondence between connectivity list and nuclear coordinates.


      REAL(KIND=dp), dimension(:,:),           intent(in) ::    prexyza
      REAL(KIND=dp), dimension(:),             intent(inout) :: precoeffap
      integer, dimension(:),             intent(in) ::    atomlista 
      integer, dimension(:,:),           intent(in) ::    connlistal
      REAL(KIND=dp), dimension(:,:),           intent(inout) :: prervecsa
      REAL(KIND=dp), dimension(:,:),           intent(inout) ::    sarray
      integer, dimension(:),             intent(in) ::    linecol
      
      REAL(KIND=dp), dimension(:,:), optional, intent(in) ::    prexyzb
      REAL(KIND=dp), dimension(:),   optional, intent(in) ::    precoeffbp
      REAL(KIND=dp), dimension(:,:), optional, intent(inout) :: prervecsb
      integer, dimension(:), optional, intent(in) ::    atomlistb
      integer, dimension(:,:), optional, intent(in) ::    connlistbl
      REAL(KIND=dp),                 optional, intent(out) ::   sab
      

      REAL(KIND=dp), dimension(:), allocatable :: stc
      REAL(KIND=dp), dimension(:), allocatable :: dist
      REAL(KIND=dp), dimension(:,:), allocatable :: xyza
      REAL(KIND=dp), dimension(:,:), allocatable :: rvecsb, rvecsa
      REAL(KIND=dp), dimension(:), allocatable :: coeffap, coeffbp 
      REAL(KIND=dp), dimension(:,:), allocatable :: expcoeffa, expcoeffb
      REAL(KIND=dp), dimension(:,:), allocatable :: minipre
      integer, dimension(:,:), allocatable :: minicon
      integer, dimension(:,:), allocatable :: grandconb
      REAL(KIND=dp), dimension(:), allocatable :: musbp
      REAL(KIND=dp) :: muap
      REAL(KIND=dp), dimension(:,:,:), allocatable :: unitvecs
      REAL(KIND=dp), dimension(8,2) :: slaters

      !integer :: stat_err, ioerr, ioread
      integer :: stat_err
      integer :: i,j,k

!      REAL(KIND=dp) :: rvec(3,1)

      !REAL(KIND=dp) :: com(3)
      REAL(KIND=dp) :: c
      REAL(KIND=dp) :: DDOT

!print*, 1221, "sarray", sarray
      slaters(:,1)=(/0.0D0, 0.0D0, 3.0D0, 0.9560D0, 1.2881D0, 1.100D0, 1.9237D0, 2.2458D0/)
      slaters(:,2)=(/0.0D0, 0.0D0, 0.0D0, 0.0000D0, 1.2107D0, 1.000D0, 1.6000D0, 2.2266D0/)

      allocate(coeffap(size(connlistal,2)),stat=stat_err)
      !if (stat_err .ne. 0) stop
      if (stat_err .ne. 0) CPABORT(" ")
     
      allocate(xyza(4, size(connlistal,2)),stat=stat_err)
      !if (stat_err .ne. 0) stop
      if (stat_err .ne. 0) CPABORT(" ")

      allocate(rvecsa(3, size(prexyza,2)),stat=stat_err)
      !   if (stat_err .ne. 0) stop
      if (stat_err .ne. 0) CPABORT(" ")

      allocate(stc(size(prexyza,2)),stat=stat_err) !????? or b
      !if (stat_err .ne. 0) stop
      if (stat_err .ne. 0) CPABORT(" ")

      allocate(minicon(4,size(linecol)),stat=stat_err)
      !if (stat_err .ne. 0) stop
      if (stat_err .ne. 0) CPABORT(" ")

      allocate(minipre(3,size(prervecsa,2)),stat=stat_err)
      !if (stat_err .ne. 0) stop
      if (stat_err .ne. 0) CPABORT(" ")

      !print*, 'One size fits all?', size(connlistal,2), size(prexyza,2), size(linecol),size(prervecsa,2), &
!&             size(sarray,1),  size(sarray,2) 

      !Identifying p_pi directions


      do i=1,size(linecol)
   
         j=linecol(i)

         minicon(:,i)=connlistal(:,j)

         !print*, (minicon(k,i), k=1,size(minicon,1))

      end do
!print*, "  "
!print*, 1276
!print*, "  "
         call calc_rvecs(prexyza(2:4,:), minipre , minicon)          


      !Updating rvecs for pass

      do i=1,size(linecol)
print*, "minipre", minipre(:,linecol(i))
           prervecsa(:,linecol(i))=minipre(:,linecol(i))

           !print*, (prervecsa(k,linecol(i)), k=1,size(prervecsa,1))         

      end do

print*, minipre
print*, "hey"
print*, minicon
     ! SIZE(minipre)
      deallocate(minipre)
      deallocate(minicon)
      !Selecting only relevant atoms

      xyza=prexyza
      coeffap=precoeffap 
!print*, "green", coeffap
      rvecsa=prervecsa
!print*, "rvecsa", rvecsa(:,1)

      ! These variables are recaculated for each atom in A

      ! Inherent for B atom

      if (present(prexyzb) .and. present(precoeffbp) .and. &
&          present(sab) .and. present(connlistbl)) then


         !print*, 'waaaargh', size(connlistbl,2) bloody cyclopentadiene fails in next line

         allocate(unitvecs(3,size(connlistbl,2),3),stat=stat_err)
         !if (stat_err .ne. 0) stop
         if (stat_err .ne. 0) CPABORT(" ")

         !print*, 'waaaargh'

         allocate(dist(size(connlistbl,2)),stat=stat_err)
         !if (stat_err .ne. 0) stop
         if (stat_err .ne. 0) CPABORT(" ")

         allocate(rvecsb(3, size(prexyzb,2)),stat=stat_err)
         !if (stat_err .ne. 0) stop
         if (stat_err .ne. 0) CPABORT(" ")

         allocate(expcoeffa(3,size(prexyzb,2)),stat=stat_err)
         !if (stat_err .ne. 0) stop
         if (stat_err .ne. 0) CPABORT(" ")

         allocate(coeffbp(size(prexyzb,2)),stat=stat_err)
         !if (stat_err .ne. 0) stop
         if (stat_err .ne. 0) CPABORT(" ")

         allocate(expcoeffb(3,size(prexyzb,2)),stat=stat_err)
         !if (stat_err .ne. 0) stop
         if (stat_err .ne. 0) CPABORT(" ")

         if(size(linecol) .eq. size(atomlista)) then

!        In the very first step rvecsb needs initialisation 

             allocate(grandconb(4,size(atomlistb)),stat=stat_err)
             !if (stat_err .ne. 0) stop
             if (stat_err .ne. 0) CPABORT(" ")

             do i=1,size(atomlistb)

                 grandconb(:,i)=connlistbl(:,atomlistb(i))

             end do

             call calc_rvecs(prexyzb(2:4,:), prervecsb, grandconb)

             deallocate(grandconb)
             

         end if

         coeffbp=precoeffbp
         rvecsb=prervecsb


         allocate(musbp(size(prexyzb,2)),stat=stat_err)
         !if (stat_err .ne. 0) stop
         if (stat_err .ne. 0) CPABORT(" ")
         
         do i=1,size(prexyzb,2)

            musbp(i)=slaters(int(prexyzb(1,i)),2)
 
         end do

        do j=1,size(linecol)

            i=linecol(j)


          muap=slaters(int(prexyza(1,i)),2)

          dist=0

          call calc_unit_vecs(prexyza(2:4,i), prexyzb(2:4,:), dist, unitvecs)



          !Do the expansion for monomer A
          call project_onto_unitvectors(unitvecs, rvecsa(:,i), expcoeffa)


          !Do the expansion for monomer B
          call project_onto_unitvectors(unitvecs, rvecsb, expcoeffb)

          !Calculate the overlaps

          call calc_overlap_p(expcoeffa, expcoeffb, muap, musbp, dist, sarray(:,i))


        end do

!print*, 1371
!print*, "stc", stc
!print*, "sarray", sarray
!print*, coeffbp
        call DGEMV('t', size(sarray,1), size(sarray,2), 1.0d0, sarray,&
                   & size(sarray,1), coeffbp, 1, 0.d0, stc, 1)
!print*, "stc", stc

!print*, "green", coeffap
!print*, "stc", stc
        sab = DDOT(size(coeffap,1),coeffap,1,stc,1)
!print*, "green", coeffap

        deallocate(expcoeffb, expcoeffa, unitvecs, dist, musbp, rvecsb, coeffbp)

   
       elseif ((present(prexyzb) .and. present(precoeffbp) .and. &
&         present(sab) .and. present(connlistbl)) .eqv. .false.) then

        !Normalize

        allocate(unitvecs(3,size(prexyza,2),3),stat=stat_err)
        !if (stat_err .ne. 0) stop
        if (stat_err .ne. 0) CPABORT(" ")

        allocate(dist(size(prexyza,2)),stat=stat_err)
        !if (stat_err .ne. 0) stop
        if (stat_err .ne. 0) CPABORT(" ")

        allocate(musbp(size(prexyza,2)),stat=stat_err)
        !if (stat_err .ne. 0) stop
        if (stat_err .ne. 0) CPABORT(" ")

        allocate(expcoeffa(3,size(prexyza,2)),stat=stat_err)
        !if (stat_err .ne. 0) stop
        if (stat_err .ne. 0) CPABORT(" ")

        allocate(expcoeffb(3,size(prexyza,2)),stat=stat_err)
        !if (stat_err .ne. 0) stop
        if (stat_err .ne. 0) CPABORT(" ")

        do i=1,size(prexyza,2)

         musbp(i)=slaters(int(prexyza(1,i)),2)

        end do

          do j=1,size(linecol)

             i=linecol(j)

             muap=slaters(int(prexyza(1,i)),2)

             call calc_unit_vecs(prexyza(2:4,i), prexyza(2:4,:), dist, unitvecs)
   
             !print*, 'unitvecs', dist(6), (unitvecs(k, 6, 3), k=1,3)

             !print*, 'Some info', prexyza(2:4,i), dist(33), (unitvecs(k, 33, 3), k=1,3)
!print*, "rvecsa", rvecsa(:,i)
             call project_onto_unitvectors(unitvecs, rvecsa(:,i), expcoeffa)
!print*, 1435, "expcoeffa", expcoeffa
             call project_onto_unitvectors(unitvecs, rvecsa, expcoeffb)
!print*, 1435, "expcoeffb", expcoeffb

             !print*, 'exp', (expcoeffb(k,10), k=1,3), (expcoeffa(k,10), k=1,3), dist(10)
!print*, 1430, sarray
            
             call calc_overlap_p(expcoeffa, expcoeffb, muap, musbp, dist, sarray(:,i))
!print*, 1433, sarray

             !print*, 'dist, sarray', dist(7),  sarray(7,i)

             !print*, 'dist, sarray', i,  dist(7),  sarray(7,i)

             !print*, (sarray(k,i), k=1,size(sarray, 1))

             do k=1,size(sarray,2)

                sarray(i,k)=sarray(k,i)
             
             end do
 
             

          end do


          !do i=1,size(sarray,2)

          ! print*, (sarray(k,i), k=1,size(sarray, 1))

          !end do

          !print*, 'line 10', (sarray(k,10), k=1,size(sarray, 1))
          !print*, 'Sarray element 8, 9', sarray(8, 9)

!print*, "green", coeffap
!print*, "stc", stc
!print*, "sarray", sarray
          call DGEMV('t', size(prexyza,2), size(prexyza,2), 1.0d0, sarray,&
                     & size(prexyza,2), coeffap, 1, 0.d0, stc, 1)
!print*, "green", coeffap
!print*, "stc", stc

          c = DDOT(size(prexyza,2),coeffap,1,stc,1)
!print*, "green", coeffap
!print*, size(prexyza,2)
!print*, coeffap
!print*, stc
!print*, "green", "c", c
          !write(*,*) 'This is c ', sqrt(c)

          coeffap = coeffap / sqrt(c)
!print*, "green", coeffap

          precoeffap = coeffap

          deallocate(expcoeffb, expcoeffa, unitvecs, dist, musbp)

      else

!         print*, 'Wrong argument passing to calc_sab'

      end if

      !deallocate(stc)
      deallocate(rvecsa, xyza, coeffap)

      return

      end subroutine calc_sab_turbo




      subroutine calc_dRSab(xyza, atomlista, connlistal, coeffap, coeffas, &
&                       nacva, nacvb, xyzb, atomlistb, connlistbl, coeffbp, coeffbs, &
&                       aneighbour, bneighbour)

!     The aim of this subroutine is to calculate the non adiabatic coupling vector 
!     using a 2 point scheme of numerical derivation (f(R+dR)-f(R))/dR. 3-point scheme
!     did not make a huge difference in results but makes the calculation much slower.

      REAL(KIND=dp), dimension(:,:), intent(in) ::    xyza
      REAL(KIND=dp), dimension(:),   intent(inout) :: coeffap, coeffas
      !REAL(KIND=dp), dimension(:), allocatable  :: coeffapdr, coeffasdr
      REAL(KIND=dp), dimension(:), allocatable  :: coeffapdr
      integer, dimension(:), intent(in) ::    atomlista
      integer, dimension(:,:), intent(in) ::    connlistal 
      REAL(KIND=dp), dimension(:,:), intent(in) ::    xyzb
      REAL(KIND=dp), dimension(:),   intent(inout) ::    coeffbp, coeffbs
      !REAL(KIND=dp), dimension(:), allocatable  ::    coeffbpdr, coeffbsdr
      REAL(KIND=dp), dimension(:), allocatable  ::    coeffbpdr
      integer, dimension(:), intent(in) ::    atomlistb
      integer, dimension(:,:), intent(in) ::    connlistbl 
      REAL(KIND=dp), dimension(:,:), intent(out) ::   nacva, nacvb
      integer, dimension(:,:), intent(in) ::  aneighbour
      integer, dimension(:,:), intent(in) ::  bneighbour


      REAL(KIND=dp), dimension(:,:), allocatable ::   sarrayn, sarrayo
      REAL(KIND=dp), dimension(:,:), allocatable ::   sarrayndr, sarrayodr
      REAL(KIND=dp), dimension(:,:), allocatable ::   rvecsa, rvecsadr
      REAL(KIND=dp), dimension(:,:), allocatable ::   rvecsb, rvecsbdr
      integer, dimension(5) :: dummylinecol 
      integer, dimension(:), allocatable :: linecol
 
      REAL(KIND=dp) :: dR
      !REAL(KIND=dp) :: dR, start, finish
      REAL(KIND=dp), dimension(:,:), allocatable :: xyzadr, xyzbdr
      !REAL(KIND=dp) :: sabRDRR, sabRR, sabRRDR  
      REAL(KIND=dp) :: sabRDRR, sabRR  
 
      integer :: stat_err, natoma, natomb
      !integer :: i, j, k, l, z, f, r
      integer :: i, j, k, z, f

      natoma=size(xyza,2)
      natomb=size(xyzb,2)

      allocate(rvecsa(3,natoma),stat=stat_err)
      !if (stat_err .ne. 0) stop
      if (stat_err .ne. 0) CPABORT(" ")

      allocate(rvecsadr(3,natoma),stat=stat_err)
      !if (stat_err .ne. 0) stop
      if (stat_err .ne. 0) CPABORT(" ")

      allocate(rvecsb(3,natomb),stat=stat_err)
      !if (stat_err .ne. 0) stop
      if (stat_err .ne. 0) CPABORT(" ")

      allocate(rvecsbdr(3,natomb),stat=stat_err)
      !if (stat_err .ne. 0) stop
      if (stat_err .ne. 0) CPABORT(" ")

      allocate(xyzadr(4,natoma),stat=stat_err)
      !if (stat_err .ne. 0) stop
      if (stat_err .ne. 0) CPABORT(" ")

      allocate(xyzbdr(4,natomb),stat=stat_err)
      !if (stat_err .ne. 0) stop
      if (stat_err .ne. 0) CPABORT(" ")

      allocate(coeffapdr(natoma),stat=stat_err)
      !if (stat_err .ne. 0) stop
      if (stat_err .ne. 0) CPABORT(" ")

      !allocate(coeffasdr(natoma),stat=stat_err)
      !if (stat_err .ne. 0) stop

      allocate(coeffbpdr(natomb),stat=stat_err)
      !if (stat_err .ne. 0) stop
      if (stat_err .ne. 0) CPABORT(" ")

      !allocate(coeffbsdr(natomb),stat=stat_err)
      !if (stat_err .ne. 0) stop


      rvecsa=0.0D0
      rvecsb=0.0d0


      dR=0.0001D0*1.889725989D0 ! in  Bohr

!     Initialising the atomic overlap matrices

      allocate(sarrayn(natoma,natoma),stat=stat_err)
      !if (stat_err .ne. 0) stop
      if (stat_err .ne. 0) CPABORT(" ")

      allocate(sarrayndr(natoma,natoma),stat=stat_err)
      !if (stat_err .ne. 0) stop
      if (stat_err .ne. 0) CPABORT(" ")

      allocate(sarrayo(natomb,natoma),stat=stat_err)
      !if (stat_err .ne. 0) stop
      if (stat_err .ne. 0) CPABORT(" ")

      allocate(sarrayodr(natomb,natoma),stat=stat_err)
      !if (stat_err .ne. 0) stop
      if (stat_err .ne. 0) CPABORT(" ")

      sarrayo=0.0D0
      sarrayn=0.0D0

!     Normalising coefficients and initialising atomic overlap matrix for wavefunction 
!     normalisation 

!print*, "blue", coeffap
!print*, 1603, sarrayn
      call calc_sab_turbo(xyza, atomlista, connlistal, coeffap, rvecsa, sarrayn, atomlista)
!print*, 1605, sarrayn
!print*, "blue", coeffap


!     Calculating Sab(R,R) and initialising p_pi directions and atomic overlap matrix sarrayo

      call calc_sab_turbo2(xyza, atomlista, connlistal, coeffap, rvecsa,&
&                                 sarrayo, atomlista, aneighbour(:,2:(size(aneighbour,2))),&
&                     sabRR, xyzb, atomlistb, connlistbl, coeffbp, rvecsb)
!print*, "blue", coeffap

!     Looping all over the structure 

      do i=1,size(connlistal, 2)
          
          f=1

          dummylinecol=0

!     Only picking the elements that are relevant to the calculation using neighbourlists
          
          do k=1,4

             z=connlistal(k,i)+1
             !print*, z, aneighbour(1,z), f*(sign(1,(aneighbour(1,z)-1))+1)/2+1
             dummylinecol(f*(sign(1,(aneighbour(1,z)-1))+1)/2+1)=z-1
             f=f+(sign(1,(aneighbour(1,z)-1))+1)/2
          end do

          allocate(linecol(f-1),stat=stat_err)
          !if (stat_err .ne. 0) stop
          if (stat_err .ne. 0) CPABORT(" ")

          linecol=dummylinecol(2:f)

          do j=1,3


              xyzadr=xyza
              xyzadr(j+1,i)=xyza(j+1,i)+dR 
              coeffapdr=coeffap        

              sarrayndr=sarrayn
              sarrayodr=sarrayo
              rvecsadr=rvecsa

              !if(i .eq. 2 .and. j .eq. 2) then

              !   do r=1,size(xyzadr, 2)

              !      print*, coeffap(i), coeffbp(i)
              !      print*, (xyzadr(k,r), k=1,size(xyzadr,1))
              !      print*, (sarrayo(k,i), k=1,size(sarrayo,1))        

              !   end do

              !end if
  
              call calc_sab_turbo(xyzadr, atomlista, connlistal, coeffapdr, rvecsadr, sarrayndr, linecol)

              call calc_sab_turbo2(xyzadr, atomlista, connlistal, coeffapdr, rvecsadr,&
&                                 sarrayodr, linecol, aneighbour(:,2:(size(aneighbour,2))),&
&                                 sabRDRR, xyzb, atomlistb, connlistbl, coeffbp, rvecsb)

             !if(i .eq. 5 .and. j .eq. 1) then
         
             !   print*, 'sab(R+dR)', sabRDRR, 'sab(R)', sabRR, 'and', (sabRDRR-sabRR)/dR

             !end if


              nacva(j,i)=(abs(sabRDRR)-abs(sabRR))/dR

          end do

          deallocate(linecol)

              !print*, nacva(1,i),  nacva(2,i), nacva(3,i)
      end do

!     These parameters change if natoma .neq. natomb

      deallocate(sarrayn, sarrayndr, sarrayo, sarrayodr)

      allocate(sarrayn(natomb,natomb),stat=stat_err)
!if (stat_err .ne. 0) stop
      if (stat_err .ne. 0) CPABORT(" ")

      allocate(sarrayndr(natomb,natomb),stat=stat_err)
      !if (stat_err .ne. 0) stop
      if (stat_err .ne. 0) CPABORT(" ")

      allocate(sarrayo(natoma,natomb),stat=stat_err)
      !if (stat_err .ne. 0) stop
      if (stat_err .ne. 0) CPABORT(" ")

      allocate(sarrayodr(natoma,natomb),stat=stat_err)
      !if (stat_err .ne. 0) stop
      if (stat_err .ne. 0) CPABORT(" ")

      sarrayo=0.0D0
      sarrayn=0.0D0

!     Repeat all over for nacvb 

      call calc_sab_turbo(xyzb, atomlistb, connlistbl, coeffbp, rvecsb, sarrayn, atomlistb)

!print*, "blue", coeffap
      call calc_sab_turbo2(xyzb, atomlistb, connlistbl, coeffbp, rvecsb, sarrayo,&
&                          atomlistb, bneighbour(:,2:(size(bneighbour,2))),&
&                         sabRR, xyza, atomlista, connlistal, coeffap, rvecsa)

!print*, "blue", coeffap
      do i=1,size(connlistbl,2)

          f=1

          dummylinecol=0

          do k=1,4

             z=connlistbl(k,i)+1
             !print*, z, bneighbour(1,z), f*(sign(1,(bneighbour(1,z)-1))+1)/2+1
             dummylinecol(f*(sign(1,(bneighbour(1,z)-1))+1)/2+1)=z-1
             f=f+(sign(1,(bneighbour(1,z)-1))+1)/2
          end do

          allocate(linecol(f-1),stat=stat_err)
          if (stat_err .ne. 0) stop
      !if (stat_err .ne. 0) CPABORT(" ")

          linecol=dummylinecol(2:f)

          do j=1,3
              
              xyzbdr=xyzb
              xyzbdr(j+1,i)=xyzb(j+1,i)+dR
              coeffbpdr=coeffbp
              rvecsbdr=rvecsb

              sarrayndr=sarrayn
              sarrayodr=sarrayo

              call calc_sab_turbo(xyzbdr, atomlistb, connlistbl, coeffbpdr, rvecsbdr, sarrayndr, linecol)


!print*, "blue", coeffap
              call calc_sab_turbo2(xyzbdr, atomlistb, connlistbl, coeffbpdr,& 
&                                  rvecsbdr, sarrayodr, linecol, bneighbour(:,2:(size(bneighbour,2))),&
&                                  sabRDRR, xyza, atomlista, connlistal,coeffap, rvecsa)
!print*, "blue", coeffap
 
              nacvb(j,i)=(sabRDRR-sabRR)/dR
          
          end do

          deallocate(linecol)

          !print*, nacvb(1,i),  nacvb(2,i), nacvb(3,i)

      end do

      deallocate(sarrayo, sarrayodr, sarrayndr, sarrayn, coeffapdr, coeffbpdr)
      deallocate(xyzbdr, xyzadr, rvecsadr, rvecsbdr, rvecsa, rvecsb)

      end subroutine

      subroutine neighbouratoms(screenlist, neighbour)

!     Sums up screenlist creating a neighbourlist where the
!     column index correpsonds to the atomic index in molecule
!     a and the numbers are the atomic indices in molecule b 

      integer, dimension(:,:), intent(out) ::  neighbour
      integer, dimension(:,:), intent(in) ::   screenlist

      integer :: i, j, l

      do i=1,size(screenlist,2)

         l=0

         do j=1,size(screenlist,1)

              l=l+screenlist(j,i)

              neighbour(l*screenlist(j,i)+1, i)=j

         end do
     
!        First row of the neighbourlist contains the number of significant atomic
!        overlaps for each atom.
 
         neighbour(1,i)=l
        
      end do

      end subroutine neighbouratoms


      subroutine newconnlist(connectlist, neighbourlist, linecol, nonzeros)

!     Updates connectivity list according to neighbour list

      implicit none
      integer, dimension(:,:), intent(in) :: connectlist, neighbourlist
      integer, allocatable, dimension(:,:) :: connect, neighb

      integer, intent(out) :: nonzeros
      integer, intent(out) :: linecol
      integer :: i, j, stat_err     

      allocate(connect(size(connectlist,1),size(connectlist,2)),stat=stat_err)
      !if (stat_err .ne. 0) stop
      if (stat_err .ne. 0) CPABORT(" ")

      allocate(neighb(size(neighbourlist,1),size(neighbourlist,2)),stat=stat_err)
      !if (stat_err .ne. 0) stop
      if (stat_err .ne. 0) CPABORT(" ")

      j=0

      do i=1,size(connectlist, 2)

          if(neighbourlist(1,i) > 0) then
             j = j + 1
             connect(:,j)=connectlist(:,i)
             neighb(:,j)=neighbourlist(:,i)       

          end if

      end do

      nonzeros=j


      deallocate(connect, neighb)
      end subroutine newconnlist


      subroutine sab_screen_new(prexyza, atomlista, connlistal, precoeffap, linecol,&
&                         prexyzb, atomlistb, connlistbl, precoeffbp, aneighbour,&
&                         bneighbour, minAO)

!     Creates neighbourlists for monomer a and b based on individual atomic overlaps.
!     The subroutine is based on calc_sab_turbo

      REAL(KIND=dp), dimension(:,:),           intent(in) :: prexyza
      REAL(KIND=dp), dimension(:),             intent(in) :: precoeffap
      integer, dimension(:,:),           intent(in) :: connlistal 
      integer, dimension(:),           intent(in) :: atomlista
      integer, dimension(:,:),           intent(out) ::  aneighbour, bneighbour
      REAL(KIND=dp), dimension(:,:), allocatable ::    sarray
      integer, dimension(:),             intent(in) ::    linecol
      
      REAL(KIND=dp), dimension(:,:),           intent(in) ::    prexyzb
      REAL(KIND=dp), dimension(:),             intent(in) ::    precoeffbp
      integer, dimension(:,:),           intent(in) ::    connlistbl 
      integer, dimension(:),             intent(in) ::   atomlistb
     ! REAL(KIND=dp) ::   sab
      integer, dimension(:,:), allocatable :: preaneighbour, prebneighbour
      integer, dimension(:,:), allocatable :: screenlist1, screenlist2 

      REAL(KIND=dp), dimension(:), allocatable :: stc
      REAL(KIND=dp), dimension(:), allocatable :: dist
      REAL(KIND=dp), dimension(:,:), allocatable :: xyza
      REAL(KIND=dp), dimension(:,:), allocatable :: rvecsb, rvecsa
      REAL(KIND=dp), dimension(:), allocatable :: coeffap, coeffbp 
      REAL(KIND=dp), dimension(:,:), allocatable :: expcoeffa, expcoeffb
      REAL(KIND=dp), dimension(:,:), allocatable :: minipre
      integer, dimension(:,:), allocatable :: minicon
      integer, dimension(:,:), allocatable :: grandconb
      REAL(KIND=dp), dimension(:), allocatable :: musbp
      REAL(KIND=dp) :: muap
      REAL(KIND=dp), intent(in) :: minAO
      REAL(KIND=dp), dimension(:,:,:), allocatable :: unitvecs
      REAL(KIND=dp), dimension(8,2) :: slaters

      !integer :: stat_err, ioerr, ioread
      integer :: stat_err
      integer :: i,j,k

      !REAL(KIND=dp) :: rvec(3,1)

      !REAL(KIND=dp) :: com(3)
      !REAL(KIND=dp) :: c
      !REAL(KIND=dp) :: DDOT

!print*, "   "
!print*, "sab_screen_new"
!print*, prexyza(1,:)
!print*, prexyza(2,:)
!print*, prexyza(3,:)
!print*, prexyza(4,:)
!print*, linecol
      slaters(:,1)=(/0.0D0, 0.0D0, 3.0D0, 0.9560D0, 1.2881D0, 1.100D0, 1.9237D0, 2.2458D0/)
      slaters(:,2)=(/0.0D0, 0.0D0, 0.0D0, 0.0000D0, 1.2107D0, 1.000D0, 1.6000D0, 2.2266D0/)

      allocate(coeffap(size(connlistal,2)),stat=stat_err)
      !if (stat_err .ne. 0) stop
      if (stat_err .ne. 0) CPABORT(" ")
     
      allocate(xyza(4, size(connlistal,2)),stat=stat_err)
      !if (stat_err .ne. 0) stop
      if (stat_err .ne. 0) CPABORT(" ")

      allocate(expcoeffa(3,size(prexyzb,2)),stat=stat_err)
      !if (stat_err .ne. 0) stop
      if (stat_err .ne. 0) CPABORT(" ")

      allocate(rvecsa(3, size(prexyza,2)),stat=stat_err)
      !if (stat_err .ne. 0) stop
      if (stat_err .ne. 0) CPABORT(" ")

      allocate(stc(size(prexyza,2)),stat=stat_err) !????? or b
      !if (stat_err .ne. 0) stop
      if (stat_err .ne. 0) CPABORT(" ")

      allocate(minicon(4,size(linecol)),stat=stat_err)
      !if (stat_err .ne. 0) stop
      if (stat_err .ne. 0) CPABORT(" ")

      allocate(minipre(3,size(prexyza,2)),stat=stat_err)
      !if (stat_err .ne. 0) stop
      if (stat_err .ne. 0) CPABORT(" ")
 
      allocate(screenlist1(size(prexyzb,2),size(prexyza,2)),stat=stat_err)
      !if (stat_err .ne. 0) stop
      if (stat_err .ne. 0) CPABORT(" ")

      allocate(screenlist2(size(prexyza,2),size(prexyzb,2)),stat=stat_err)
      !if (stat_err .ne. 0) stop
      if (stat_err .ne. 0) CPABORT(" ")

      screenlist1=0
      screenlist2=0

!     Identifying p_pi directions

      do i=1,size(linecol)
   
         j=linecol(i)

         minicon(:,i)=connlistal(:,j)

      end do

!print*, "call calc_rvecs(prexyza(2:4,:), minipre , minicon)"          
      call calc_rvecs(prexyza(2:4,:), minipre , minicon)          

!     Updating rvecs for pass

      do i=1,size(linecol)

           rvecsa(:,linecol(i))=minipre(:,linecol(i))

      end do

      deallocate(minipre, minicon)

!     Selecting only relevant atoms

      xyza=prexyza
      coeffap=precoeffap 

!     These variables are recaculated for each atom in A

!     Inherent for B atom

      allocate(unitvecs(3,size(connlistbl,2),3),stat=stat_err)
      !if (stat_err .ne. 0) stop
      if (stat_err .ne. 0) CPABORT(" ")

      allocate(dist(size(connlistbl,2)),stat=stat_err)
      !if (stat_err .ne. 0) stop
      if (stat_err .ne. 0) CPABORT(" ")

      allocate(rvecsb(3, size(prexyzb,2)),stat=stat_err)
      !if (stat_err .ne. 0) stop
      if (stat_err .ne. 0) CPABORT(" ")

      allocate(coeffbp(size(prexyzb,2)),stat=stat_err)
      !if (stat_err .ne. 0) stop
      if (stat_err .ne. 0) CPABORT(" ")

      allocate(expcoeffb(3,size(prexyzb,2)),stat=stat_err)
      !if (stat_err .ne. 0) stop
      if (stat_err .ne. 0) CPABORT(" ")

      allocate(sarray(size(prexyzb,2),size(prexyza,2)),stat=stat_err)
      !if (stat_err .ne. 0) stop
      if (stat_err .ne. 0) CPABORT(" ")


      sarray=0.0D0

      allocate(grandconb(4,size(atomlistb)),stat=stat_err)
      !if (stat_err .ne. 0) stop
      if (stat_err .ne. 0) CPABORT(" ")

      do i=1,size(atomlistb)

         grandconb(:,i)=connlistbl(:,atomlistb(i))

      end do

      call calc_rvecs(prexyzb(2:4,:), rvecsb, grandconb)
 
      deallocate(grandconb)


      coeffbp=precoeffbp

      allocate(musbp(size(prexyzb,2)),stat=stat_err)
      !if (stat_err .ne. 0) stop
      if (stat_err .ne. 0) CPABORT(" ")
         
      do i=1,size(prexyzb,2)

         musbp(i)=slaters(int(prexyzb(1,i)),2)
 
      end do


      do j=1,size(linecol)

         i=linecol(j)

         muap=slaters(int(prexyza(1,i)),2)

         dist=0

         call calc_unit_vecs(prexyza(2:4,i), prexyzb(2:4,:), dist, unitvecs)

!        Do the expansion for monomer A

         call project_onto_unitvectors(unitvecs, rvecsa(:,i), expcoeffa)

!        Do the expansion for monomer B

         call project_onto_unitvectors(unitvecs, rvecsb, expcoeffb)

!        Calculate the overlaps

         call calc_overlap_p(expcoeffa, expcoeffb, muap, musbp, dist, sarray(:,i))


         !minAO=1.0D-17

!        If the atomic overlap between to atoms of either of the molecule is smaller than
!        a certain 'minAO' value the given (aindex, bindex) position in screenlist becomes
!        0 and otherwise it is 1. 

         do k=1,size(sarray,1)

            screenlist1(k,i)=(int(sign(1.0D0, (abs(sarray(k,i)*coeffbp(k)*coeffap(i))-minAO)))+1)/2

         end do

       end do

      allocate(preaneighbour(size(prexyzb,2)+1,size(prexyza,2)),stat=stat_err)
      !if (stat_err .ne. 0) stop
      if (stat_err .ne. 0) CPABORT(" ")

!     Creating aneighbour list

      call neighbouratoms(screenlist1, aneighbour)

      deallocate(preaneighbour)


!     Creating bneighbour list

      screenlist2=transpose(screenlist1)

      allocate(prebneighbour(size(xyza,2)+1,size(prexyzb,2)),stat=stat_err)
      !if (stat_err .ne. 0) stop
      if (stat_err .ne. 0) CPABORT(" ")

      call neighbouratoms(screenlist2, bneighbour)

      !do i=1,size(bneighbour,2)

      !   print*, (bneighbour(j,i), j=1,size(bneighbour,1))

      !end do

      deallocate(prebneighbour)

      deallocate(expcoeffb)
      deallocate(rvecsb)
      deallocate(stc,dist,expcoeffa,unitvecs)
      deallocate(rvecsa)

      return
      end subroutine sab_screen_new

      subroutine calc_sab_turbo2(prexyza, atomlista, connlistal, precoeffap, prervecsa,&
&                                sarray, linecol, aneighbour, sab, prexyzb, atomlistb, &
&                                connlistbl, precoeffbp, prervecsb)

!     Updated version of calc_sab_turbo using neighbourlist to reduce calculation time 
!     Sadly it cannot be used for norm calculation (atoms too close)


      REAL(KIND=dp), dimension(:,:),           intent(in) ::    prexyza
      REAL(KIND=dp), dimension(:),             intent(inout) :: precoeffap
      integer, dimension(:,:),           intent(in) ::    connlistal
      integer, dimension(:),             intent(in) ::  atomlista
      integer, dimension(:,:),           intent(in) ::  aneighbour
      REAL(KIND=dp), dimension(:,:),           intent(inout) :: prervecsa
      REAL(KIND=dp), dimension(:,:),           intent(inout) ::    sarray
      integer, dimension(:),             intent(in) ::    linecol
      
      REAL(KIND=dp), dimension(:,:), optional, intent(in) ::    prexyzb
      REAL(KIND=dp), dimension(:),   optional, intent(in) ::    precoeffbp
      REAL(KIND=dp), dimension(:,:), optional, intent(inout) :: prervecsb
      integer, dimension(:,:), optional, intent(in) ::    connlistbl
      integer, dimension(:),   optional, intent(in) ::    atomlistb
      REAL(KIND=dp),                 optional, intent(out) ::   sab
      

      REAL(KIND=dp), dimension(:), allocatable :: stc
      REAL(KIND=dp), dimension(:), allocatable :: dist
      REAL(KIND=dp), dimension(:,:), allocatable :: xyza, xyzb
      REAL(KIND=dp), dimension(:,:), allocatable :: rvecsb, rvecsa
      REAL(KIND=dp), dimension(:), allocatable :: coeffap, coeffbp 
      REAL(KIND=dp), dimension(:,:), allocatable :: expcoeffa, expcoeffb
      REAL(KIND=dp), dimension(:,:), allocatable :: minipre
      integer, dimension(:,:), allocatable :: minicon
      REAL(KIND=dp), dimension(:), allocatable :: musbp, sarrayline
      REAL(KIND=dp) :: muap
      REAL(KIND=dp), dimension(:,:,:), allocatable :: unitvecs
      integer, dimension(:,:), allocatable :: grandconb
      REAL(KIND=dp), dimension(8,2) :: slaters

      !integer :: stat_err, ioerr, ioread
      integer :: stat_err
      integer :: i,j,k

      !REAL(KIND=dp) :: rvec(3,1)

!      REAL(KIND=dp) :: com(3)
      REAL(KIND=dp) :: c
      REAL(KIND=dp) :: DDOT

      slaters(:,1)=(/0.0D0, 0.0D0, 3.0D0, 0.9560D0, 1.2881D0, 1.100D0, 1.9237D0, 2.2458D0/)
      slaters(:,2)=(/0.0D0, 0.0D0, 0.0D0, 0.0000D0, 1.2107D0, 1.000D0, 1.6000D0, 2.2266D0/)

      allocate(coeffap(size(connlistal,2)),stat=stat_err)
      !if (stat_err .ne. 0) stop
      if (stat_err .ne. 0) CPABORT(" ")
     
      allocate(xyza(4, size(connlistal,2)),stat=stat_err)
      !if (stat_err .ne. 0) stop
      if (stat_err .ne. 0) CPABORT(" ")

      allocate(rvecsa(3, size(prexyza,2)),stat=stat_err)
      !   if (stat_err .ne. 0) stop
      if (stat_err .ne. 0) CPABORT(" ")

      allocate(stc(size(prexyza,2)),stat=stat_err) !????? or b
      !if (stat_err .ne. 0) stop
      if (stat_err .ne. 0) CPABORT(" ")

      allocate(minicon(4,size(linecol)),stat=stat_err)
      !if (stat_err .ne. 0) stop
      if (stat_err .ne. 0) CPABORT(" ")

      allocate(minipre(3,size(prervecsa,2)),stat=stat_err)
      !if (stat_err .ne. 0) stop
      if (stat_err .ne. 0) CPABORT(" ")

!     Identifying p_pi directions

      do i=1,size(linecol)
   
         j=linecol(i)

         minicon(:,i)=connlistal(:,j)

      end do

      call calc_rvecs(prexyza(2:4,:), minipre , minicon)          

!     Updating rvecs for pass

      do i=1,size(linecol)

           prervecsa(:,linecol(i))=minipre(:,linecol(i))

      end do

      deallocate(minipre, minicon)

!     Selecting only relevant atoms

      xyza=prexyza
      coeffap=precoeffap 
      rvecsa=prervecsa

!     These variables are recaculated for each atom in A

!     Inherent for B atom

      if (present(prexyzb) .and. present(precoeffbp) .and. &
&         present(sab) .and. present(connlistbl)) then

!         When the subroutine is called with the full atomlist it calculates
!         rvecsb that does not happen when rvecsb is alredy calculated e. g. 
!         in nacv calculation subroutine  

          if(size(linecol) .eq. size(atomlista)) then

             allocate(grandconb(4,size(atomlistb)),stat=stat_err)
             !if (stat_err .ne. 0) stop
             if (stat_err .ne. 0) CPABORT(" ")

             do i=1,size(atomlistb)

                 grandconb(:,i)=connlistbl(:,atomlistb(i))

             end do

             call calc_rvecs(prexyzb(2:4,:), prervecsb, grandconb)

             deallocate(grandconb)

          end if

          allocate(coeffbp(size(prexyzb,2)),stat=stat_err)
          !if (stat_err .ne. 0) stop
          if (stat_err .ne. 0) CPABORT(" ")

          coeffbp=precoeffbp

          do j=1,size(linecol)
 
             i=linecol(j)


!            Calculating only those elelments that are present in the relevant
!            line of aneighbour.

             allocate(xyzb(4,aneighbour(1,i)),stat=stat_err)
             !if (stat_err .ne. 0) stop
             if (stat_err .ne. 0) CPABORT(" ")

             do k=1,aneighbour(1,i)

                 xyzb(:,k)=prexyzb(:,aneighbour(k+1,i))

             end do

             allocate(expcoeffa(3,size(xyzb,2)),stat=stat_err)
             !if (stat_err .ne. 0) stop
             if (stat_err .ne. 0) CPABORT(" ")

             allocate(unitvecs(3,size(xyzb,2),3),stat=stat_err)
             !if (stat_err .ne. 0) stop
             if (stat_err .ne. 0) CPABORT(" ")

             allocate(dist(size(xyzb,2)),stat=stat_err)
             !if (stat_err .ne. 0) stop
             if (stat_err .ne. 0) CPABORT(" ")

             allocate(rvecsb(3, size(xyzb,2)),stat=stat_err)
             !if (stat_err .ne. 0) stop
             if (stat_err .ne. 0) CPABORT(" ")

             allocate(expcoeffb(3,size(xyzb,2)),stat=stat_err)
             !if (stat_err .ne. 0) stop
             if (stat_err .ne. 0) CPABORT(" ")

             allocate(sarrayline(size(xyzb,2)),stat=stat_err)
             !if (stat_err .ne. 0) stop
             if (stat_err .ne. 0) CPABORT(" ")

             do k=1,aneighbour(1,i)

                 rvecsb(:,k)=prervecsb(:,aneighbour(k+1,i))

             end do

             allocate(musbp(size(prexyzb,2)),stat=stat_err)
             !if (stat_err .ne. 0) stop
             if (stat_err .ne. 0) CPABORT(" ")
         
             do k=1,size(xyzb,2)

                 musbp(k)=slaters(int(xyzb(1,k)),2)
 
             end do

             muap=slaters(int(prexyza(1,i)),2)

             dist=0


             call calc_unit_vecs(prexyza(2:4,i), xyzb(2:4,:), dist, unitvecs)

!            Do the expansion for monomer A
             call project_onto_unitvectors(unitvecs, rvecsa(:,i), expcoeffa)

!            Do the expansion for monomer B
             call project_onto_unitvectors(unitvecs, rvecsb, expcoeffb)

!            Calculate the overlaps

             call calc_overlap_p(expcoeffa, expcoeffb, muap, musbp, dist, sarrayline)

!            Changing relelvant elements of atomic overlap matrix

             do k=1,aneighbour(1,i)

                sarray(aneighbour(k+1,i),i)=sarrayline(k)

             end do

             deallocate(expcoeffa, expcoeffb, rvecsb, sarrayline)
             deallocate(dist, unitvecs, musbp, xyzb)

        end do

        call DGEMV('t', size(sarray,1), size(sarray,2), 1.0d0, sarray,&
                   & size(sarray,1), coeffbp, 1, 0.d0, stc, 1)

        sab = DDOT(size(coeffap,1),coeffap,1,stc,1)

   
       elseif ((present(prexyzb) .and. present(precoeffbp) .and. &
&         present(sab) .and. present(connlistbl)) .eqv. .false.) then

        !Normalize

        allocate(unitvecs(3,size(prexyza,2),3),stat=stat_err)
        !if (stat_err .ne. 0) stop
        if (stat_err .ne. 0) CPABORT(" ")

        allocate(dist(size(prexyza,2)),stat=stat_err)
        !if (stat_err .ne. 0) stop
        if (stat_err .ne. 0) CPABORT(" ")

        allocate(musbp(size(prexyza,2)),stat=stat_err)
        !if (stat_err .ne. 0) stop
        if (stat_err .ne. 0) CPABORT(" ")

        allocate(expcoeffa(3,size(prexyza,2)),stat=stat_err)
        !if (stat_err .ne. 0) stop
        if (stat_err .ne. 0) CPABORT(" ")

        allocate(expcoeffb(3,size(prexyza,2)),stat=stat_err)
        !if (stat_err .ne. 0) stop
        if (stat_err .ne. 0) CPABORT(" ")

        do i=1,size(prexyza,2)

            musbp(i)=slaters(int(prexyza(1,i)),2)

        end do

        do j=1,size(linecol)

             i=linecol(j)

             muap=slaters(int(prexyza(1,i)),2)

             call calc_unit_vecs(prexyza(2:4,i), prexyza(2:4,:), dist, unitvecs)

             call project_onto_unitvectors(unitvecs, rvecsa(:,i), expcoeffa)

             call project_onto_unitvectors(unitvecs, rvecsa, expcoeffb)

             call calc_overlap_p(expcoeffa, expcoeffb, muap, musbp, dist, sarray(:,i))

             do k=1,size(sarray,2)

                sarray(i,k)=sarray(k,i)
             
             end do

        end do

        call DGEMV('t', size(prexyza,2), size(prexyza,2), 1.0d0, sarray,&
                     & size(prexyza,2), coeffap, 1, 0.d0, stc, 1)

        c = DDOT(size(prexyza,2),coeffap,1,stc,1)

        coeffap = coeffap / sqrt(c)

        precoeffap = coeffap

        deallocate(expcoeffb)

      else

!         print*, 'Wrong argument passing to calc_sab'

      end if

      deallocate(stc)
      deallocate(rvecsa)
      return
      end subroutine calc_sab_turbo2

      subroutine calc_sab_with_s(prexyza, connlista, precoeffap, precoeffas, &
&                         sab, prexyzb, connlistb, precoeffbp, precoeffbs)

!     Original S_ab calculation program can calculate an pi_conjugated molecule
!
!     - Slater coefficients for carbon can be read in from external data file
!       by enabling lines 779-786
!     - Inclusion of s orbitals can be achieved by enabling lines 876-883,916- 920, 
!       930-944, 949, 974-981, 1011-1017 and 1036-1049


      REAL(KIND=dp), dimension(:,:),           intent(in) ::    prexyza
      REAL(KIND=dp), dimension(:),             intent(inout) :: precoeffap, precoeffas
      integer, dimension(:,:),           intent(in) ::    connlista
      REAL(KIND=dp), dimension(:,:), optional, intent(in) ::    prexyzb
      REAL(KIND=dp), dimension(:),   optional, intent(in) ::    precoeffbp, precoeffbs
      integer, dimension(:,:), optional, intent(in) ::    connlistb
      REAL(KIND=dp),                 optional, intent(out) ::   sab

      REAL(KIND=dp), dimension(:), allocatable :: stc
      REAL(KIND=dp), dimension(:), allocatable :: dist
      REAL(KIND=dp), dimension(:,:), allocatable :: rvecsa, rvecsb
      REAL(KIND=dp), dimension(:,:), allocatable :: prervecsa, prervecsb 
      REAL(KIND=dp), dimension(:), allocatable :: coeffap, coeffas, coeffbp, coeffbs
      REAL(KIND=dp), dimension(:,:), allocatable :: sarrayp, sarrays, sarraysp, sarrayps
      REAL(KIND=dp), dimension(:,:), allocatable :: expcoeffa, expcoeffb
      REAL(KIND=dp), dimension(:,:), allocatable :: xyza, xyzb
      REAL(KIND=dp), dimension(:), allocatable :: musbs, musbp
      REAL(KIND=dp) :: muap, muas
      REAL(KIND=dp), dimension(:,:,:), allocatable :: unitvecs
      REAL(KIND=dp), dimension(8,2) :: slaters
 
      !integer :: stat_err, ioerr, ioread
      integer :: stat_err 
      !integer :: i,j,k,z
      integer :: i,j
 
      !REAL(KIND=dp) :: rvec(3,1)
 
      !REAL(KIND=dp) :: com(3)
      REAL(KIND=dp) :: c
      REAL(KIND=dp) :: DDOT

      slaters(:,1)=(/0.0D0, 0.0D0, 3.0D0, 0.9560D0, 1.2881D0, 1.100D0, 1.9237D0, 2.2458D0/)
      slaters(:,2)=(/0.0D0, 0.0D0, 0.0D0, 0.0000D0, 1.2107D0, 1.000D0, 1.6000D0, 2.2266D0/)

      !open(unit=69, file='slats', status='old', action='read', iostat=ioerr)
      !if (ioerr .ne. 0) then
      !  print*, 'Could not open file slats'
      !  stop
      !else
      !   read(69,*, iostat=ioread) slaters(6,1), slaters(6,2)
      !end if
      !close(69)

      allocate(prervecsa(3,size(prexyza,2)),stat=stat_err)
      !if (stat_err .ne. 0) stop
        if (stat_err .ne. 0) CPABORT(" ")

      allocate(xyza(4,size(connlista,2)),stat=stat_err)
      !if (stat_err .ne. 0) stop
        if (stat_err .ne. 0) CPABORT(" ")

      allocate(coeffap(size(connlista,2)),stat=stat_err)
      !if (stat_err .ne. 0) stop
        if (stat_err .ne. 0) CPABORT(" ")

      allocate(coeffas(size(connlista,2)),stat=stat_err)
      !if (stat_err .ne. 0) stop
        if (stat_err .ne. 0) CPABORT(" ")

      allocate(rvecsa(3,size(connlista,2)),stat=stat_err)
      !if (stat_err .ne. 0) stop
        if (stat_err .ne. 0) CPABORT(" ")

      allocate(stc(size(xyza,2)),stat=stat_err)
      !if (stat_err .ne. 0) stop
        if (stat_err .ne. 0) CPABORT(" ")


      call calc_rvecs(prexyza(2:4,:), prervecsa, connlista)

      do i=1,size(connlista,2)

           xyza(:,i)=prexyza(:,connlista(1,i))
           coeffap(i)=precoeffap(connlista(1,i))
           coeffas(i)=precoeffas(connlista(1,i))
           rvecsa(:,i)=prervecsa(:,connlista(1,i))

      end do

      deallocate(prervecsa)
 
      ! These variables are recaculated for each atom in A

      ! Inherent for B atom

      if (present(prexyzb) .and. present(precoeffbp) .and. &
&         present(precoeffbs) .and. present(sab) .and. present(connlistb)) then

         allocate(unitvecs(3,size(connlistb,2),3),stat=stat_err)
         !if (stat_err .ne. 0) stop
         if (stat_err .ne. 0) CPABORT(" ")

         allocate(dist(size(connlistb,2)),stat=stat_err)
         !if (stat_err .ne. 0) stop
         if (stat_err .ne. 0) CPABORT(" ")

         allocate(prervecsb(3,size(prexyzb,2)),stat=stat_err)
         !if (stat_err .ne. 0) stop
         if (stat_err .ne. 0) CPABORT(" ")
 
         allocate(xyzb(4, size(connlistb,2)),stat=stat_err)
         !if (stat_err .ne. 0) stop
         if (stat_err .ne. 0) CPABORT(" ")

         allocate(expcoeffa(3,size(xyzb,2)),stat=stat_err)
         !if (stat_err .ne. 0) stop
         if (stat_err .ne. 0) CPABORT(" ")

         allocate(coeffbp(size(xyzb,2)),stat=stat_err)
         !if (stat_err .ne. 0) stop
         if (stat_err .ne. 0) CPABORT(" ")
  
         allocate(coeffbs(size(xyzb,2)),stat=stat_err)
         !if (stat_err .ne. 0) stop
         if (stat_err .ne. 0) CPABORT(" ")

         allocate(expcoeffb(3,size(xyzb,2)),stat=stat_err)
         !if (stat_err .ne. 0) stop
         if (stat_err .ne. 0) CPABORT(" ")

         allocate(rvecsb(3,size(connlistb,2)),stat=stat_err)
         !if (stat_err .ne. 0) stop
         if (stat_err .ne. 0) CPABORT(" ")

         call calc_rvecs(prexyzb(2:4,:), prervecsb, connlistb)

         do i=1,size(connlistb,2)

             xyzb(:,i)=prexyzb(:,connlistb(1,i))
             coeffbp(i)=precoeffbp(connlistb(1,i))
             coeffbs(i)=precoeffbs(connlistb(1,i))
             rvecsb(:,i)=prervecsb(:,connlistb(1,i))

        end do  

  
        deallocate(prervecsb)
       
        allocate(musbs(size(xyzb,2)),stat=stat_err)
        !if (stat_err .ne. 0) stop
        if (stat_err .ne. 0) CPABORT(" ")

        allocate(musbp(size(xyzb,2)),stat=stat_err)
        !if (stat_err .ne. 0) stop
        if (stat_err .ne. 0) CPABORT(" ")

        allocate(sarrayp(size(xyzb,2),size(xyza,2)), stat=stat_err)
        !if (stat_err .ne. 0) stop
        if (stat_err .ne. 0) CPABORT(" ")

        allocate(sarrays(size(xyzb,2),size(xyza,2)), stat=stat_err)
        !if (stat_err .ne. 0) stop
        if (stat_err .ne. 0) CPABORT(" ")

        allocate(sarraysp(size(xyzb,2),size(xyza,2)), stat=stat_err)
        !if (stat_err .ne. 0) stop
        if (stat_err .ne. 0) CPABORT(" ")

        allocate(sarrayps(size(xyzb,2),size(xyza,2)), stat=stat_err)
        !if (stat_err .ne. 0) stop
        if (stat_err .ne. 0) CPABORT(" ")

        do i=1,size(xyzb,2)

         musbs(i)=slaters(int(xyzb(1,i)),1)
         musbp(i)=slaters(int(xyzb(1,i)),2)

        end do


        do i=1,size(xyza,2)

          muas=slaters(int(xyza(1,i)),1)
          muap=slaters(int(xyza(1,i)),2)
         
          !calculate unitvectors for each atom in monomer A
          !and distance array
          call calc_unit_vecs(xyza(2:4,i), xyzb(2:4,:), dist, unitvecs)

          !Do the expansion for monomer A
          call project_onto_unitvectors(unitvecs, rvecsa(:,i), expcoeffa)

          !Do the expansion for monomer B

          call project_onto_unitvectors(unitvecs, rvecsb, expcoeffb)

          !Calculate the overlaps

          call calc_overlap_p(expcoeffa, expcoeffb, muap, musbp, dist, sarrayp(:,i))

          call calc_overlap_s(muas, musbs, dist, sarrays(:,i))

          call calc_overlap_sp(expcoeffb(3,:), muas, musbp, dist, sarraysp(:,i))

          call calc_overlap_sp(expcoeffa(3,:), muap, musbs, dist, sarrayps(:,i))

        end do

        call DGEMV('t', size(sarrayp,1), size(sarrayp,2), 1.0d0, sarrayp,&
                   & size(sarrayp,1), coeffbp, 1, 0.d0, stc, 1)

        sab = DDOT(size(coeffap,1),coeffap,1,stc,1)

        !print*, size(sarrayp,1), size(sarrayp,2)

        call DGEMV('t', size(sarrays,1), size(sarrays,2), 1.0d0, sarrays,&
                    & size(sarrays,1), coeffbs, 1, 0.d0, stc, 1)

        sab = sab+DDOT(size(coeffas,1),coeffas,1,stc,1)


        call DGEMV('t', size(sarraysp,1), size(sarraysp,2), 1.0d0, sarraysp,&
                   & size(sarrayp,1), coeffbp, 1, 0.d0, stc, 1)

        sab = sab + DDOT(size(coeffas,1),coeffas,1,stc,1)

        call DGEMV('t', size(sarrayps,1), size(sarrayps,2), 1.0d0, sarrayps,&
                     & size(sarrayps,1), coeffbs, 1, 0.d0, stc, 1)

        sab = sab+DDOT(size(coeffap,1),coeffap,1,stc,1)


        deallocate(rvecsb,  expcoeffa, expcoeffb, sarrayp, dist, unitvecs, coeffbp, musbp)
        deallocate(sarrays, sarrayps, sarraysp, coeffbs, musbs)


      elseif ((present(prexyzb) .and. present(precoeffbp) .and. &
&         present(precoeffbs) .and. present(sab) .and. present(connlistb)) .eqv. .false.) then

        !Normalize

        allocate(unitvecs(3,size(xyza,2),3),stat=stat_err)
        !if (stat_err .ne. 0) stop
        if (stat_err .ne. 0) CPABORT(" ")

        allocate(dist(size(xyza,2)),stat=stat_err)
        !if (stat_err .ne. 0) stop
        if (stat_err .ne. 0) CPABORT(" ")

        allocate(musbs(size(xyza,2)),stat=stat_err)
        !if (stat_err .ne. 0) stop
        if (stat_err .ne. 0) CPABORT(" ")

        allocate(musbp(size(xyza,2)),stat=stat_err)
        !if (stat_err .ne. 0) stop
        if (stat_err .ne. 0) CPABORT(" ")

        allocate(sarrayp(size(xyza,2),size(xyza,2)), stat=stat_err)
        !if (stat_err .ne. 0) stop
        if (stat_err .ne. 0) CPABORT(" ")

        allocate(sarrays(size(xyza,2),size(xyza,2)), stat=stat_err)
        !if (stat_err .ne. 0) stop
        if (stat_err .ne. 0) CPABORT(" ")

        allocate(sarraysp(size(xyza,2),size(xyza,2)), stat=stat_err)
        !if (stat_err .ne. 0) stop
        if (stat_err .ne. 0) CPABORT(" ")

        allocate(sarrayps(size(xyza,2),size(xyza,2)), stat=stat_err)
        !if (stat_err .ne. 0) stop
        if (stat_err .ne. 0) CPABORT(" ")

        allocate(expcoeffa(3,size(xyza,2)),stat=stat_err)
        !if (stat_err .ne. 0) stop
        if (stat_err .ne. 0) CPABORT(" ")

        allocate(expcoeffb(3,size(xyza,2)),stat=stat_err)
        !if (stat_err .ne. 0) stop
        if (stat_err .ne. 0) CPABORT(" ")


        sarrayp=0

        do i=1,size(xyza,2)

         musbs(i)=slaters(int(xyza(1,i)),1)
         musbp(i)=slaters(int(xyza(1,i)),2)

        end do

          do i=2,size(xyza,2)

           dist=0.0D0

           unitvecs=0.0D0

           expcoeffa=0.0D0

           expcoeffb=0.0D0



             muas=slaters(int(xyza(1,i)),1)
             muap=slaters(int(xyza(1,i)),2)

!            Only calculates lower triangle of the overlap matrix in the normalisation step
!            then mirrors it (all real) over the diagonal 
 
             j=i-1        

             call calc_unit_vecs(xyza(2:4,i), xyza(2:4,1:j), dist(1:j), unitvecs(:,1:j,:))    
           
             !print*, 'unitvecs', dist(2), (unitvecs(k, 2, 3), k=1,3)

             call project_onto_unitvectors(unitvecs(:,1:j,:), rvecsa(:,i), expcoeffa(:,1:j))

             call project_onto_unitvectors(unitvecs(:,1:j,:), rvecsa(:,1:j), expcoeffb(:,1:j))

             !print*, 'expcoeffa'

             !do z=1,size(expcoeffa, 2)

             !      print*, (expcoeffa(k,z), k=1,3)

             !end do

             !print*, 'expcoeffb'


             !do z=1,size(expcoeffb, 2)

             !      print*, (expcoeffb(k,z), k=1,3)

             !end do

             !print*, 'dist'

             !do z=1,size(dist)

             !      print*, dist(z)

             !end do

             !do z=1,size(sarrayp,2)
             
             !     print*, (sarrayp(k,z), k=1,size(sarrayp, 1))

             !end do

             !print*, 'mua', muap

             !print*, 'mub', musbp



             call calc_overlap_p(expcoeffa(:,1:j), expcoeffb(:,1:j), muap, musbp(1:j), dist(1:j), sarrayp(1:j,i))
            !print*, 'exp', (expcoeffb(k,1), k=1,3), (expcoeffa(k,1), k=1,3), dist(1)

             !call calc_overlap_p(expcoeffa, expcoeffb, muap, musbp, dist, sarrayp(:,i))
             
             call calc_overlap_s(muas, musbs(1:j), dist(1:j), sarrays(1:j,i))

             call calc_overlap_sp(expcoeffb(3,1:j), muas, musbp(1:j), dist(1:j), sarraysp(1:j,i))

             call calc_overlap_sp(expcoeffa(3,1:j), muap, musbs(1:j), dist(1:j), sarrayps(1:j,i))
   

          end do

          sarrayp=sarrayp+transpose(sarrayp)         
          sarraysp = sarraysp + transpose(sarraysp)
          sarrayps = sarrayps + transpose(sarrayps)
          sarrays = sarrayp + transpose(sarrays)


          do i=1,size(xyza,2)

              sarrayp(i,i)=1.0D0
              sarraysp(i,i)=0.0D0
              sarrayps(i,i)=0.0D0
              sarrays(i,i)=1.0D0

          end do


          !do i=1,size(sarrayp,2)

          ! print*, (sarrayp(k,i), k=1,size(sarrayp, 1))

          !end do




          call DGEMV('t', size(xyza,2), size(xyza,2), 1.0d0, sarrayp,&
                     & size(xyza,2), coeffap, 1, 0.d0, stc, 1)

          c = DDOT(size(xyza,2),coeffap,1,stc,1)

          call DGEMV('t', size(xyza,2), size(xyza,2), 1.0d0, sarrays,&
                     & size(xyza,2), coeffas, 1, 0.d0, stc, 1)

          c = c+DDOT(size(xyza,2),coeffas,1,stc,1)

          call DGEMV('t', size(sarraysp,1), size(sarraysp,2), 1.0d0, sarraysp,&
                   & size(sarraysp,2), coeffap, 1, 0.d0, stc, 1)

          c = c+DDOT(size(coeffas,1),coeffas,1,stc,1)

          call DGEMV('t', size(sarrayps,1), size(sarrayps,2), 1.0d0, sarrayps,&
                     & size(sarrayps,2), coeffas, 1, 0.d0, stc, 1)

          c = c+DDOT(size(coeffap,1),coeffap,1,stc,1)

           

          !write(*,*) 'This is c ',sqrt(c)

          coeffap = coeffap / sqrt(c) 
          coeffas = coeffas / sqrt(c)

          do i=1,size(connlista,2)

              precoeffap(connlista(1,i)) = coeffap(i)
              precoeffas(connlista(1,i)) = coeffas(i)

          end do
         
         deallocate(sarrayp, expcoeffa, expcoeffb, dist, unitvecs, musbp) 
         deallocate(sarrayps, sarraysp, sarrays)


      else

         print*, 'Wrong argument passing to calc_sab' 

      end if

      deallocate(rvecsa, xyza, coeffap, stc)


      return
      end subroutine calc_sab_with_s

      end module jacob_overlapfinal

