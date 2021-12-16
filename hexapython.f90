program hexabrick

    use hexamod
    
    implicit none
    
    integer :: nbox, perim, newperim
    integer :: addbox
    integer :: chooseedge, nb
    real :: rnum
    integer :: x, y
    !complex, allocatable :: bdystart(:), bdyend(:), tempbdystart(:), tempbdyend(:), bdystartOUT(:), bdyendOUT(:)
    integer, allocatable :: orientable(:), temporientable(:), orientableOUT(:), skippedneigh(:)
    integer, allocatable :: orbrick(:)
    !complex, allocatable :: effbdystart(:), effbdyend(:)
    complex, allocatable :: bdy(:)
    complex, allocatable :: corners(:)
    integer, allocatable :: poscorner(:)
    complex, allocatable :: PA(:)
    integer :: indexi, indexj, sizeS, nijN, taille, k, l
    integer, allocatable :: sizeStot1(:,:),sizeStot2(:,:)
    integer, allocatable :: w(:,:), nbmpm(:,:)
    integer, allocatable :: intPA(:,:)
    
    logical :: error
    
    integer :: iter, niter, i, nrun, run
    real :: start, finish
    
    call cpu_time(start)
        ! put code to test here
    
    
    nrun = 1
    niter = 1
    !iter = 1
    
    !open (1, file = 'testhexa.dat', status = 'new')
    !open (30, file = 'areaVperim200.dat', status = 'new')
    
    
    !$omp parallel 
    do run = 1,nrun
      outer: do iter = 1, niter

      !nbox = 500+(run-1)*500
      nbox = 100000
      
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

            !write(30,*) nbox,perim
            print*, nbox,perim
      
                  !print*, 'orientable', orientable
      
                  !call orienttocurve(orientable, bdy, perim)

            !       Functions for counting immersions __________________________________
                  
                  call orienttobrick(orientable, bdy, orbrick, perim)
      
                  call cornerPA(orbrick, bdy, corners, poscorner, perim)
      
                  call findPA(PA, orbrick, bdy, corners, poscorner, perim)
      
                  taille = size(corners)

                  call integerPA(PA, taille, intPA)



            !       !print*, intPA
      
            !       ! indexi=1
            !       ! indexj=taille
            !       !indexj=3
      
            ! !   allocate(sizeStot1(taille,taille))
            ! !   do k=1,taille
            ! !   do l=k+1,taille
            ! !         sizeStot1(l,k)=0
            ! !   end do
            ! !   end do
      
            ! !   call sizempm(PA, taille, sizeS, indexi, indexj, sizeStot1)
            ! !   print*, 's1(1,n)', sizeStot1
            ! !   do i=1,taille
            ! !      print*, 'sizeStot1', sizeStot1(i,:)
            ! !   end do
      
                  call loopsizempm(PA, taille, w, sizeStot2)
                  !print*, 's2(1,n)', sizeStot2
                  !do i=1,taille
                  !   print*, 'sizeStot2', sizeStot2(i,:)
                  !end do
                  
            ! !   do i=1,taille
            ! !       print*, 'sizeStot', sizeStot1(i,:)-sizeStot2(i,:)
            ! !   end do
      
                  call loopnbmpm(PA, taille, w, sizeStot2, nbmpm)
                  !do i=1,taille
                  !   print*, 'nbmpm', nbmpm(i,:)
                  !end do
                  if (nbmpm(1,taille)>1) then 
                        print*, 'nbmpm', nbmpm(1,taille)
                  end if
      
            !       !call numbermpm(PA, taille, nijN, indexi, indexj, sizeStot)
            !       !print*, 'n(1,n)', nijN
      
            !       !print*, orientable

            !       if (nbmpm(1,taille)==0) then 

            !             print*, corners
            !             print*, 'taille:', size(corners)
            !             print*, 'PA:', PA

            !             open (21, file = 'intPA.dat', status = 'new')
            !             do i=1,size(PA)
            !                   x = intPA(i,1)
            !                   y = intPA(i,2)
            !                   write(21,*) x,y
            !             end do
            !             close (21)
            
            !             open (2, file = 'orientablehexa.dat', status = 'new')
            !             do i=0, size(bdy)-1
            !                   x = nint(realpart(bdy(i)) * 1000.0) * 1E-3
            !                   y = nint(imagpart(bdy(i)) * 1000.0) * 1E-3
            !                   write(2,*) '(',x,',',y,')'
            !             end do
            !             close (2)

            !             exit outer

            !       end if

            !       End of Functions for counting immersions __________________________________
            
                  !write(1,*) orientable
      
                  deallocate(orientable)
                  deallocate(orbrick)
                  deallocate(PA)
                  deallocate(intPA)
                  deallocate(sizeStot2)
                  deallocate(w)
                  deallocate(nbmpm)
                  deallocate(corners)
                  deallocate(poscorner)
                  deallocate(bdy)
      
            end do outer
      end do 
      !$omp end parallel 
      
      !close (30)
    
      !open (1, file = 'randomWalknew.dat', status = 'new')
      !      do i=0, size(bdy)-1
      !         write(1,*) bdy(i)
      !      end do
      !close (1)
    
    
      !close (1)
      call cpu_time(finish)
    print '("Time = ",f6.3," seconds.")',finish-start
    
    
    end program hexabrick
    