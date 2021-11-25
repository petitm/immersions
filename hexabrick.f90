program hexabrick

use hexamod

implicit none

integer :: nbox, perim, newperim
integer :: addbox
integer :: chooseedge, nb
real :: rnum, x, y
!complex, allocatable :: bdystart(:), bdyend(:), tempbdystart(:), tempbdyend(:), bdystartOUT(:), bdyendOUT(:)
integer, allocatable :: orientable(:), temporientable(:), orientableOUT(:), skippedneigh(:)
integer, allocatable :: orbrick(:)
!complex, allocatable :: effbdystart(:), effbdyend(:)
complex, allocatable :: bdy(:)
complex, allocatable :: corners(:)
integer, allocatable :: poscorner(:)
complex, allocatable :: PA(:)
integer :: indexi, indexj, sizeS, nijN, taille, k, l
integer, allocatable :: sizeStot(:,:)

logical :: error

integer :: iter, niter, i

nbox = 100
niter = 1
!iter = 1

!open (1, file = 'testhexa.dat', status = 'new')

outer: do iter = 1, niter

        ! Constructing the immersed disk

        perim = 6

        ! Initialization

        allocate(orientable(0:perim-1))

        orientable(0)=0
        orientable(1)=1
        orientable(2)=2
        orientable(3)=3
        orientable(4)=4
        orientable(5)=5

        !print*, orientable

      do addbox = 1,nbox-1

        error = .true.

        do while (error .eqv. .true.)

              allocate(skippedneigh(0:perim-1))

              call subeffbdy(perim, orientable, newperim, orientableOUT, skippedneigh)

              call singutest(orientableOUT, newperim, error)

              if (error .eqv. .false.) then
                    perim = newperim
              end if

              deallocate(skippedneigh)

        end do

        call move_alloc(orientableOUT,orientable)

        !print*, 'orientable', orientable

      end do

        !print*, 'orientable', orientable

        !call orienttocurve(orientable, bdy, perim)
        call orienttobrick(orientable, bdy, orbrick, perim)

        call cornerPA(orbrick, bdy, corners, poscorner, perim)

        call findPA(PA, orbrick, bdy, corners, poscorner, perim)

        print*, corners
        print*, 'taille:', size(corners)
        print*, 'PA:', PA

        taille = size(corners)

        indexi=1
        indexj=taille
        !indexj=3

        allocate(sizeStot(taille,taille))
        do k=1,taille
        do l=k+1,taille
              sizeStot(l,k)=0
        end do
        end do

        call sizempm(PA, taille, sizeS, indexi, indexj, sizeStot)
        print*, 's(1,n)', sizeS
        print*, 'sizeStot', sizeStot

        !call numbermpm(PA, taille, nijN, indexi, indexj)
        !print*, 'n(1,n)', nijN

        !print*, orientable

        open (2, file = 'orientablehexa.dat', status = 'new')
        do i=0, size(bdy)-1
              x = nint(realpart(bdy(i)) * 1000.0) * 1E-3
              y = nint(imagpart(bdy(i)) * 1000.0) * 1E-3
              write(2,*) '(',x,',',y,')'
        end do
        close (2)

        !write(1,*) orientable

        deallocate(orientable)

  end do outer

  !open (1, file = 'randomWalknew.dat', status = 'new')
  !      do i=0, size(bdy)-1
  !         write(1,*) bdy(i)
  !      end do
  !close (1)


  !close (1)


end program hexabrick
