module hexamod

      implicit none

      private
      public subeffbdy
      public convexe, concave1, concave2, concave3, concave4, orienttocurve
      public singutest, orienttobrick, cornerPA, findPA, sizempm, numbermpm
      public loopsizempm, loopnbmpm
      public integerPA

contains

      subroutine subeffbdy(perim, orientable, newperim, orientableOUT, skippedneigh)
            implicit none

            integer, intent(inout) :: perim
            integer, intent(in) :: orientable(0:perim-1)
            integer, intent(inout) :: skippedneigh(0:perim-1)
            integer, allocatable, intent(out) :: orientableOUT(:)
            integer, intent(out) :: newperim

            integer :: n, angle, chooseedge, pos, concorn, position
            real :: rnum


            do n=0, perim-1

                  angle = orientable(n)-orientable(modulo(n+1,perim))
                  !print*, 'angle', angle

                  ! skippedneigh follows the bdy and counts the concave corners
                  if (angle == 1 .or. angle == -5) then

                        skippedneigh(n) = 1
                  else
                        skippedneigh(n) = 0

                  end if

            end do

            call random_seed()
            call random_number(rnum)
            chooseedge = floor(rnum*perim)

            do while (skippedneigh(chooseedge)>0)
                  call random_seed()
                  call random_number(rnum)
                  chooseedge = floor(rnum*perim)
            end do

            concorn = 0
            pos = 1

            do while (skippedneigh(modulo(chooseedge-pos, perim)) .ne. 0)
                  concorn=concorn+1
                  pos=pos+1
            end do

            if (concorn==0) then
                  call convexe(orientable, orientableOUT, perim, newperim, chooseedge)

            else if (concorn==1) then
                  call concave1(orientable, orientableOUT, perim, newperim, chooseedge)

            else if (concorn==2) then
                  call concave2(orientable, orientableOUT, perim, newperim, chooseedge)

            else if (concorn==3) then
                  call concave3(orientable, orientableOUT, perim, newperim, chooseedge)

            else
                  allocate(orientableOUT(0:perim-1))
                  orientableOUT(:) = orientable(:)
                  newperim=perim

            end if


      end subroutine subeffbdy

      subroutine convexe(orientable, orientableOUT, perim, newperim, chooseedge)

            implicit none
            integer, intent(inout) :: perim
            integer, intent(in) :: orientable(0:perim-1)
            integer, intent(in) :: chooseedge

            integer, intent(out), allocatable :: orientableOUT(:)
            integer, intent(out) :: newperim

            newperim=perim+4

            allocate(orientableOUT(0:newperim-1))

            if (chooseedge>0) then
                  orientableOUT(:chooseedge-1) = orientable(:chooseedge-1)
            end if

            if (chooseedge<perim-1) then
                  orientableOUT(chooseedge+5:) = orientable(chooseedge+1:)
            end if

            if (chooseedge==0) then
                  orientableOUT(5:newperim-1) = orientable(1:perim-1)
            end if

            if (chooseedge==perim-1) then
                  orientableOUT(0:newperim-6) = orientable(0:perim-2)
            end if


            orientableOUT(chooseedge)=modulo(orientable(chooseedge)+4,6)
            orientableOUT(modulo(chooseedge+1, newperim))=modulo(orientableOUT(chooseedge)+1,6)
            orientableOUT(modulo(chooseedge+2, newperim))=modulo(orientableOUT(chooseedge)+2,6)
            orientableOUT(modulo(chooseedge+3, newperim))=modulo(orientableOUT(chooseedge)+3,6)
            orientableOUT(modulo(chooseedge+4, newperim))=modulo(orientableOUT(chooseedge)+4,6)


      end subroutine convexe

      subroutine concave1(orientable, orientableOUT, perim, newperim, chooseedge)

            implicit none
            integer, intent(inout) :: perim
            integer, intent(in) :: orientable(0:perim-1)
            integer, intent(in) :: chooseedge

            integer, intent(out), allocatable :: orientableOUT(:)
            integer, intent(out) :: newperim

            newperim=perim+2

            allocate(orientableOUT(0:newperim-1))


            if (chooseedge > 1 .and. chooseedge<perim-1) then

                  orientableOUT(0:chooseedge-2) = orientable(0:chooseedge-2)
                  !orientableOUT(chooseedge-)=modulo(orientable(chooseedge)+1,6)
                  orientableOUT(chooseedge+3:newperim-1)=orientable(chooseedge+1:perim-1)

                  orientableOUT(modulo(chooseedge-1, newperim))=modulo(orientable(chooseedge)-1,6)
                  orientableOUT(chooseedge)=orientable(chooseedge)
                  orientableOUT(modulo(chooseedge+1, newperim))=modulo(orientableOUT(chooseedge)+1,6)
                  orientableOUT(modulo(chooseedge+2, newperim))=modulo(orientableOUT(chooseedge)+2,6)

            else

                  orientableOUT(0:newperim-5) = orientable(modulo(chooseedge+1,perim):modulo(perim-2+chooseedge,perim))

                  orientableOUT(newperim-4)=modulo(orientable(chooseedge)-1,6)
                  orientableOUT(newperim-3)=orientable(chooseedge)
                  orientableOUT(newperim-2)=modulo(orientable(chooseedge)+1,6)
                  orientableOUT(newperim-1)=modulo(orientable(chooseedge)+2,6)

            end if


      end subroutine concave1

      subroutine concave2(orientable, orientableOUT, perim, newperim, chooseedge)

            implicit none
            integer, intent(inout) :: perim
            integer, intent(in) :: orientable(0:perim-1)
            integer, intent(in) :: chooseedge

            integer, intent(out), allocatable :: orientableOUT(:)
            integer, intent(out) :: newperim

            newperim=perim

            allocate(orientableOUT(0:newperim-1))

            orientableOUT(:)=orientable(:)

            orientableOUT(modulo(chooseedge-2, newperim))=orientable(chooseedge)
            orientableOUT(modulo(chooseedge-1, newperim))=modulo(orientable(chooseedge)+1,6)
            orientableOUT(chooseedge)=modulo(orientable(chooseedge)+2,6)


      end subroutine concave2

      subroutine concave3(orientable, orientableOUT, perim, newperim, chooseedge)

            implicit none
            integer, intent(inout) :: perim
            integer, intent(in) :: orientable(0:perim-1)
            integer, intent(in) :: chooseedge

            integer, intent(out), allocatable :: orientableOUT(:)
            integer, intent(out) :: newperim

            newperim=perim-2

            allocate(orientableOUT(0:newperim-1))

            if (chooseedge > 3 .and. chooseedge<perim-1) then

                  orientableOUT(0:chooseedge-4) = orientable(0:chooseedge-4)
                  orientableOUT(chooseedge-3)=modulo(orientable(chooseedge)+1,6)
                  orientableOUT(chooseedge-2)=modulo(orientable(chooseedge)+2,6)
                  orientableOUT(chooseedge-1:newperim-1)=orientable(chooseedge+1:perim-1)

            else

                  orientableOUT(0:newperim-3) = orientable(modulo(chooseedge+1,perim):modulo(perim-4+chooseedge,perim))
                  orientableOUT(newperim-2)=modulo(orientable(chooseedge)+1,6)
                  orientableOUT(newperim-1)=modulo(orientable(chooseedge)+2,6)

            end if

            !if (chooseedge == 1) then
            !orientableOUT(modulo(chooseedge-3, newperim))=modulo(orientable(chooseedge)+1,6)
            !else if (chooseedge == 0) then


            !orientableOUT(modulo(chooseedge-2, newperim))=modulo(orientable(chooseedge)+2,6)



      end subroutine concave3

      subroutine concave4(orientable, orientableOUT, perim, newperim, chooseedge)

            implicit none
            integer, intent(inout) :: perim
            integer, intent(in) :: orientable(0:perim-1)
            integer, intent(in) :: chooseedge

            integer, intent(out), allocatable :: orientableOUT(:)
            integer, intent(out) :: newperim

            newperim=perim-4

            allocate(orientableOUT(0:newperim-1))

            if (chooseedge > 4 .and. chooseedge<perim-1) then

                  orientableOUT(0:chooseedge-5) = orientable(0:chooseedge-5)
                  orientableOUT(chooseedge-4)=modulo(orientable(chooseedge)+2,6)
                  orientableOUT(chooseedge-3:newperim-1)=orientable(chooseedge+1:perim-1)

            else

                  orientableOUT(0:perim-6) = orientable(modulo(chooseedge+1,perim):modulo(perim-5+chooseedge,perim))
                  orientableOUT(perim-5)=modulo(orientable(chooseedge)+2,6)

            end if

            !else if (chooseedge == 4) then
            !orientableOUT(0:perim-6) = orientable(5:perim-1)
            !else if (chooseedge == 3) then
            !orientableOUT(0:perim-6) = orientable(4:perim-2)
            !else if (chooseedge == 2) then
            !orientableOUT(0:perim-6) = orientable(3:perim-3)
            !else if (chooseedge == 1) then
            !orientableOUT(0:perim-6) = orientable(2:perim-4)
            !else if (chooseedge == 0) then
            !orientableOUT(0:perim-6) = orientable(1:perim-5)
            !else if (chooseedge == perim-1) then
            !orientableOUT(0:perim-6) = orientable(0:perim-6)


      end subroutine concave4

      subroutine orienttocurve(orientable, bdy, perim)
      implicit none
      integer, intent(in) :: perim
      integer, intent(in) :: orientable(0:perim-1)
      complex, intent(out), allocatable :: bdy(:)

      integer :: i
      complex :: dir
      real :: a, b, PI

      !perim = size(orientable)

      allocate(bdy(0:perim))

      bdy(0) = (0.,0.)
      PI=4.*ATAN(1.)
      a = cos(PI/6.)
      b = sin(PI/6.)

      do i=0,size(orientable)-1

            if (orientable(i)==0) then
                  dir = cmplx(a,b)
            else if (orientable(i)==1) then
                  dir = (0.,1.)
            else if (orientable(i)==2) then
                  dir = cmplx(-a,b)
            else if (orientable(i)==3) then
                  dir = cmplx(-a,-b)
            else if (orientable(i)==4) then
                  dir = (0.,-1.)
            else if (orientable(i)==5) then
                  dir = cmplx(a,-b)
            end if

            bdy(i+1) = bdy(i) + dir

      end do

      print*, bdy(perim)

      end subroutine orienttocurve

      subroutine singutest(orientable, perim, error)
      implicit none
      integer, intent(in) :: perim
      integer, intent(in) :: orientable(0:perim-1)

      logical, intent(inout) :: error

      integer :: test, i

      error = .false.

      outer : do i=0,perim-1

            test = orientable(i) - orientable(modulo(i+1, perim))

            if (test == 2 .or. test == -2) then

                  error = .true.
                  exit outer
                  print*, 'Ahah!'

            end if

      end do outer

      end subroutine singutest

      subroutine orienttobrick(orientable, bdy, orbrick, perim)
      implicit none
      integer, intent(in) :: perim
      integer, intent(in) :: orientable(0:perim-1)
      complex, intent(out), allocatable :: bdy(:)
      integer, intent(out), allocatable :: orbrick(:)

      integer :: i
      complex :: dir
      !real :: a, b, PI

      !perim = size(orientable)

      allocate(bdy(0:perim))
      allocate(orbrick(0:perim-1))


      bdy(0) = (0.,0.)

      do i=0,size(orientable)-1

            if (orientable(i)==0) then
                  dir = (1.,0.)
                  orbrick(i) = 0
            else if (orientable(i)==1) then
                  dir = (0.,1.)
                  orbrick(i) = 1
            else if (orientable(i)==2) then
                  dir = (-1.,0.)
                  orbrick(i) = 2
            else if (orientable(i)==3) then
                  dir = (-1.,0.)
                  orbrick(i) = 2
            else if (orientable(i)==4) then
                  dir = (0.,-1.)
                  orbrick(i) = 3
            else if (orientable(i)==5) then
                  dir = (1.,0.)
                  orbrick(i) = 0
            end if

            bdy(i+1) = bdy(i) + dir

      end do

      print*, bdy(perim)

      end subroutine orienttobrick


      subroutine cornerPA(orbrick, bdy, corners, poscorner, perim)
      implicit none
      integer, intent(in) :: perim
      integer, intent(in) :: orbrick(0:perim-1)
      complex, intent(in) :: bdy(0:perim)
      complex, intent(out), allocatable :: corners(:)
      integer, intent(out), allocatable :: poscorner(:)

      integer :: n, angle

      allocate(corners(0))
      allocate(poscorner(0))

      do n=0, perim-1

            angle = orbrick(n)-orbrick(modulo(n+1,perim))
            !print*, 'angle', angle

            ! corners follows the bdy and counts the concave corners
            if (angle == 1 .or. angle == -3) then

                  corners = [corners,bdy(n+1)]
                  poscorner = [poscorner,n]

            end if

      end do

      end subroutine cornerPA

      subroutine findPA(PA, orbrick, bdy, corners, poscorner, perim)
      implicit none
      integer, intent(in) :: perim
      integer, intent(in) :: orbrick(0:perim-1)
      complex, intent(in) :: bdy(0:perim)
      complex, intent(in) :: corners(:)
      integer, intent(in) :: poscorner(:)
      complex, allocatable, intent(out) :: PA(:)

      integer :: n,m, taille, i
      real :: height, start, end
      logical :: cond
      complex, allocatable :: interset(:)

      taille = size(corners)

      allocate(PA(0))

      do n=1, taille

      do m=n+1, taille

            cond = .false.

            !height=abs(imagpart(bdy(poscorner(n)+1)) - imagpart(bdy(poscorner(m)+1)))
            height=abs(imagpart(corners(n)) - imagpart(corners(m)))

            if (height<0.01) then

                  cond = .true.

                  start=min(realpart(corners(n)),realpart(corners(m)))
                  end=max(realpart(corners(n)),realpart(corners(m)))

                  if (start.ne.end) then

                        allocate(interset(0))

                        do i=1,floor(end-start)-1
                              interset = [interset,cmplx(start,imagpart(corners(n)))+i*(1.,0.)]
                        end do

                        if (any(interset == bdy(poscorner(n))) .or. &
                            any(interset == bdy(modulo(poscorner(n)+2, perim+1))) .or. &
                            any(interset == bdy(poscorner(m))) .or. &
                            any(interset == bdy(modulo(poscorner(m)+2, perim+1)))) then

                              cond = .false.

                        end if

                        !print*, 'interset', interset
                        !print*, 'n,m', n,m
                        !print*, 'bdy(poscorner(n))', bdy(poscorner(n))
                        !print*, 'bdy(modulo(poscorner(n)+2, perim+1))', bdy(modulo(poscorner(n)+2, perim+1))
                        !print*, 'bdy(poscorner(m))',bdy(poscorner(m))
                        !print*, 'bdy(modulo(poscorner(m)+2, perim+1)))', bdy(modulo(poscorner(m)+2, perim+1))


                        deallocate(interset)

                  end if

            end if

            if (cond .eqv. .true.) then

                  PA=[PA,cmplx(n,m)]

            end if


      end do
      end do

      end subroutine findPA

      recursive subroutine sizempm(PA, taille, sizeS, i, j, sizeStot)
            implicit none
            integer, intent(in) :: taille, i,j
            complex, intent(in) :: PA(:)
            integer, intent(out) :: sizeS
            integer, intent(inout), allocatable :: sizeStot(:,:)

            integer :: m, v, siv, svj, wiv
            logical :: cond


            if (i.ge.j .or. i>taille .or. j>taille .or. i<1 .or. j<1) then
                  m=0
            else
                  call sizempm(PA, taille, sizeS, i+1, j, sizeStot)
                  m = sizeS
            end if

            do v=i+1, j

                  if (any(PA == cmplx(i,v))) then

                        wiv=1

                        call sizempm(PA, taille, siv, i+1, v-1, sizeStot)
                        !print*, i+1,v-1,'siv', siv
                        sizeStot(i+1,v-1)=siv

                        call sizempm(PA, taille, svj, v+1, j, sizeStot)
                        sizeStot(v+1,j)=svj

                        m = max(m,wiv + siv + svj)

                   end if

            end do

            sizeS=m
            sizeStot(i,j)=m

      end subroutine sizempm

      recursive subroutine numbermpm(PA, taille, nijN, i, j, sizeStot)
            implicit none
            integer, intent(in) :: taille, i,j
            complex, intent(in) :: PA(:)
            integer, intent(out) :: nijN
            integer, intent(in) :: sizeStot(:,:)

            integer :: sij, sipj, v, sizeS, nN, siv, svj
            integer :: Niv, Nvj, wiv

            if (i.ge.j .or. i>taille .or. j>taille .or. i<1 .or. j<1) then
                  nijN=1
            else

                  !call sizempm(PA, taille, sij, i, j, sizeStot)
                  !call sizempm(PA, taille, sipj, i+1, j, sizeStot)

                  !if (sij==sipj) then
                  if (sizeStot(i,j)==sizeStot(i+1,j)) then
                        call numbermpm(PA, taille, nN, i+1, j, sizeStot)
                  else
                        nN = 0
                  end if

                  do v=i+1, j

                        if (any(PA == cmplx(i,v))) then

                              wiv=1
                              !call sizempm(PA, taille, siv, i+1, v-1, sizeStot)
                              !call sizempm(PA, taille, svj, v+1, j, sizeStot)

                              !if (sij == wiv+siv+svj) then
                              if (sizeStot(i,j)==wiv+sizeStot(i+1,v-1)+sizeStot(v+1,j)) then
                                    call numbermpm(PA, taille, Niv, i+1, v-1, sizeStot)
                                    call numbermpm(PA, taille, Nvj, v+1, j, sizeStot)
                                    nN = nN + Niv*Nvj
                              end if

                        end if

                 end do

                 nijN = nN
            end if

      end subroutine numbermpm

      subroutine loopsizempm(PA, taille, w, sizeStot)

            implicit none
            integer, intent(in) :: taille
            complex, intent(in) :: PA(:)
            !integer, intent(out) :: sizeS
            integer, intent(out), allocatable :: sizeStot(:,:), w(:,:)

            integer :: i, j, next, v, val          

            allocate(w(taille,taille))
            allocate(sizeStot(taille,taille))

            w(:,:)=0

            do i=1,taille
                  do j=i+1, taille
                        if (any(PA == cmplx(i,j))) then
                              w(i,j)=1
                        end if
                  end do 
            end do

            ! do i=1,taille
            !       print*, 'w', w(i,:)
            !   end do

            sizeStot(:,:)=0

            do i=1,taille-1

                  if (w(i,i+1)==1) then
                        sizeStot(i,i+1)=1
                  end if

            end do 

            do i=1,taille-2

                  sizeStot(i,i+2)=max(w(i,i+2),w(i+1,i+2),w(i,i+1))

            end do 

            do i=1,taille-3

                  sizeStot(i,i+3)=max(sizeStot(i+1,i+3),sizeStot(i,i+2),w(i,i+3)+sizeStot(i+1,i+2),w(i,i+1)+sizeStot(i+2,i+3))

            end do 

            do next=4,taille-1
                  do i=1, taille-next
                        sizeStot(i,i+next)=0
                        do v=i+1,i+next

                              if (w(i,v)==0) then 
                                    val=0
                              else 
                                    val=1+sizeStot(i+1,max(i+1,v-1))+sizeStot(min(i+next,v+1),i+next)
                              end if 
                              
                              sizeStot(i,i+next)=max(sizeStot(i,i+next),sizeStot(i+1,i+next),val)
                        end do
                  end do 
            end do

      end subroutine

      subroutine loopnbmpm(PA, taille, w, sizeStot, nbmpm)

            implicit none
            integer, intent(in) :: taille
            complex, intent(in) :: PA(:)
            integer, intent(in) :: sizeStot(:,:), w(:,:)
            integer, intent(out), allocatable :: nbmpm(:,:)

            integer :: i, j, next, v

            allocate(nbmpm(taille,taille))

            do i=1,taille
                  do j=i,taille
                        nbmpm(j,i)=1
                  end do 
            end do

            do i=1,taille-1
                  ! if (w(i,i+1)==1) then 
                  !       nbmpm(i,i+1)=1
                  ! else 
                  !       nbmpm(i,i+1)=0
                  ! end if 
                  nbmpm(i,i+1)=1
            end do 

            do next=2,taille-1
                  do i=1,taille-next

                        if (sizeStot(i,i+next)==sizeStot(i+1,i+next)) then 
                              nbmpm(i,i+next)=nbmpm(i+1,i+next)
                        else 
                              nbmpm(i,i+next)=0
                        end if

                        do v=i+1,i+next
                              if (w(i,v)==1.and.sizeStot(i,i+next)==1+sizeStot(i+1,max(i+1,v-1))+sizeStot(min(i+next,v+1),i+next)) &
                              then 
                              
                                    nbmpm(i,i+next)=nbmpm(i,i+next)+nbmpm(i+1,max(i+1,v-1))*nbmpm(min(v+1,i+next),i+next)
                              
                              end if 
                        end do
                  end do
            end do

      end subroutine

      subroutine integerPA(PA, taille, intPA)

            implicit none
            integer, intent(in) :: taille
            complex, intent(in) :: PA(:)
            integer, intent(out), allocatable :: intPA(:,:)

            integer :: i, sizePA

            sizePA = size(PA(:))

            !print*, sizePA

            allocate(intPA(sizePA,2))

            do i=1,sizePA
                  intPA(i,1)=int(realpart(PA(i)))
                  intPA(i,2)=int(imagpart(PA(i)))
            end do

      end subroutine


end module hexamod
