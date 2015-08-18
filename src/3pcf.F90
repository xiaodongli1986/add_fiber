program ThreePCF
  use kdtree2_module
  use kdtree2_precision_module
  implicit none
#ifdef MPI
	include 'mpif.h'
#endif

  real, dimension(:,:), allocatable :: my_array, bins(:,:),tbins(:,:)
  double precision, allocatable ::  v1(:),v2(:)
  type(kdtree2), pointer :: tree

  double precision, allocatable :: Zddd(:),Zddr(:),Zdrr(:),Zrrr(:),crap(:), wgt1(:)
  double precision, allocatable :: Zdd(:,:),Zdr(:,:),Zrr(:,:),crap2(:,:)
  double precision :: ndd,ndr,nrr,dist,r1min,r1max,r2min,r2max
  integer :: k, i, j, l,d,chunk,nbins,ind
  character :: file1*100, file2*100
  type(kdtree2_result),allocatable :: resultsb(:)
  integer   ::  Ndata,Nrand,nn1,nn2,nnpairs,mm1,mm2
  double precision  :: aux,rmin,rmax,dr,odr,theta,odt,mid1,mid2
  logical :: wgt,saveran,loadran
  CHARACTER(LEN=100) :: outfile,ranfile
  integer :: myid , ntasks , ierr
#ifdef MPI 
  integer , dimension( MPI_STATUS_SIZE ) :: status
#endif
  integer , parameter :: master=0,msgtag1=11,msgtag2=12,msgtag3=13,msgtag4=14

! ! I n i t i a l i z a t i o n of the MPI environment
myid=0
ntasks=1

#ifdef MPI 
 call MPI_INIT( ierr )
 call MPI_COMM_RANK( MPI_COMM_WORLD , myid, ierr )
 call MPI_COMM_SIZE( MPI_COMM_WORLD , ntasks , ierr )
#endif

!--- set some default parameters ---
  d=3
  wgt=.false.
  nbins=0
  rmin=0.0
  rmax=0.0
  outfile='result.txt'
  loadran=.false.
  saveran=.false.

  if(myid==master) print*,'reading options'
  call parseOptions()

 allocate(bins(nbins+2,2))
 allocate(tbins(nbins,2))

   bins(2,1)=r1min
   bins(2,2)=r1max
   bins(1,1)=r2min
   bins(1,2)=r2max

   if(myid==master)  print*,'allocated bins'

  if(bins(2,1)<bins(1,1)) then
    print*,'side 1 must be larger than side 2'
    stop "exiting"
  endif


  call count_files()

  Allocate(my_array(d,Ndata+Nrand))

  if ( wgt ) then
    Allocate(wgt1(Ndata+Nrand))
  endif

  allocate(v1(d))
  allocate(v2(d))

  call read_files()

print*,' building tree on node', myid
tree => kdtree2_create(my_array,sort=.true.,rearrange=.true.)     ! this is how you create a tree.
print*,' built tree on node', myid

  call allocate_arrays()

call make_bins()

print*,'beginning data loop on thread:',myid

#ifdef MPI
chunk=floor(float(Ndata)/ntasks)+1
do i=max(myid*chunk,1),min(((myid+1)*chunk)-1,Ndata),1
#else
do i=1,Ndata,1
#endif

   call kdtree2_r_nearest_around_point(tp=tree,idxin=i,correltime=-1,r2=bins(2,2)*bins(2,2),&
&       nfound=nn2,nalloc=(Nrand+Ndata),results=resultsb)
   
      nn1=kdtree2_r_count_around_point(tp=tree,idxin=i,correltime=-1,r2=bins(2,1)*bins(2,1))      
      mm1=kdtree2_r_count_around_point(tp=tree,idxin=i,correltime=-1,r2=bins(1,1)*bins(1,1))
      mm2=kdtree2_r_count_around_point(tp=tree,idxin=i,correltime=-1,r2=bins(1,2)*bins(1,2))

      do j=nn1+1,nn2,1
        v1=my_array(:,resultsb(j)%idx)
        v1=[v1(1)-my_array(1,i),v1(2)-my_array(2,i),v1(3)-my_array(3,i)] ! connecting vector
        do k=mm1+1,mm2,1
!          v1=my_array(:,resultsb(j)%idx)
          v2=my_array(:,resultsb(k)%idx)

        if(resultsb(j)%idx==resultsb(k)%idx) goto 33

          call triplet_bin( )
        
          if(resultsb(j)%idx <= Ndata) then 
          
            if(resultsb(k)%idx <= Ndata) then
            
              if ( wgt ) then
                Zddd(ind)=Zddd(ind)+(wgt1(resultsb(j)%idx)*wgt1(resultsb(k)%idx)*wgt1(i))
              else
                Zddd(ind)=Zddd(ind)+1.d0
              endif

            else

              if ( wgt ) then
                Zddr(ind)=Zddr(ind)+(wgt1(resultsb(j)%idx)*wgt1(resultsb(k)%idx)*wgt1(i))
              else
                Zddr(ind)=Zddr(ind)+1.d0
              endif
            endif

          else 

            if(resultsb(k)%idx <=Ndata) then
!              if ( wgt ) then
!                Zddr(ind)=Zddr(ind)+(wgt1(resultsb(j)%idx)*wgt1(resultsb(k)%idx)*wgt1(i))
!              else
!                Zddr(ind)=Zddr(ind)+1.d0
!              endif
                continue
                
            else

              if ( wgt ) then
                Zdrr(ind)=Zdrr(ind)+(wgt1(resultsb(j)%idx)*wgt1(resultsb(k)%idx)*wgt1(i))
              else
                Zdrr(ind)=Zdrr(ind)+1.d0
              endif

            endif
          endif
33 continue
      enddo
   enddo


   call kdtree2_r_nearest_around_point(tp=tree,idxin=i,correltime=-1,r2=bins(3,2)*bins(3,2),&
&       nfound=nn2,nalloc=(Nrand+Ndata),results=resultsb)
   do k=nbins+2,1,-1
      nn1=kdtree2_r_count_around_point(tp=tree,idxin=i,correltime=-1,r2=bins(k,1)*bins(k,1))
      nn2=kdtree2_r_count_around_point(tp=tree,idxin=i,correltime=-1,r2=bins(k,2)*bins(k,2))

      do j=nn1+1,nn2,1
         if(resultsb(j)%idx <=Ndata) then            
            if ( wgt ) then
               Zdd(k,1)=Zdd(k,1)+(wgt1(resultsb(j)%idx)*wgt1(i))
            else
               Zdd(k,1)=Zdd(k,1)+1.d0
            endif
         else
            if ( wgt ) then
               Zdr(k,1)=Zdr(k,1)+(wgt1(resultsb(j)%idx)*wgt1(i))
            else
               Zdr(k,1)=Zdr(k,1)+1.d0
            endif
         endif
      enddo
    enddo

enddo

print*,'finished data loop in thread:', myid

#ifdef MPI
call MPI_Barrier(MPI_COMM_WORLD,ierr)
#endif

if(loadran) then
  if(myid==master) print*,'Random counts exist, skipping some calculation'
  goto 55
endif

!RR loop
print*,'beginning random loop on thread:',myid

#ifdef MPI
chunk=floor(float(Nrand)/ntasks)+1
do i=Ndata+max(myid*chunk,1),Ndata+min(((myid+1)*chunk)-1,Nrand),1
#else
do i=Ndata+1,Ndata+Nrand,1
#endif

   call kdtree2_r_nearest_around_point(tp=tree,idxin=i,correltime=-1,r2=bins(2,2)*bins(2,2),&
&       nfound=nn2,nalloc=(Nrand+Ndata),results=resultsb)

      nn1=kdtree2_r_count_around_point(tp=tree,idxin=i,correltime=-1,r2=bins(2,1)*bins(2,1))      
      mm1=kdtree2_r_count_around_point(tp=tree,idxin=i,correltime=-1,r2=bins(1,1)*bins(1,1))
      mm2=kdtree2_r_count_around_point(tp=tree,idxin=i,correltime=-1,r2=bins(1,2)*bins(1,2))

      do j=nn1+1,nn2,1
        v1=my_array(:,resultsb(j)%idx)
        v1=[v1(1)-my_array(1,i),v1(2)-my_array(2,i),v1(3)-my_array(3,i)] ! connecting vector
        
        do k=mm1+1,mm2,1
!          v1=my_array(:,resultsb(j)%idx)
          v2=my_array(:,resultsb(k)%idx)

        if(resultsb(j)%idx==resultsb(k)%idx) goto 34

         call triplet_bin( )

        if(resultsb(j)%idx <=Ndata) then !RD
!          if(resultsb(k)%idx <=Ndata) then !RDD
!            if ( wgt ) then
!               Zddr(ind)=Zddr(ind)+(wgt1(resultsb(j)%idx)*wgt1(resultsb(k)%idx)*wgt1(i))
!            else
!               Zddr(ind)=Zddr(ind)+1.d0
!            endif
!         else !RRD
!            if ( wgt ) then
!               Zdrr(ind)=Zdrr(ind)+(wgt1(resultsb(j)%idx)*wgt1(resultsb(k)%idx)*wgt1(i))
!            else
!               Zdrr(ind)=Zdrr(ind)+1.d0
!            endif
!          endif
            continue
                                   
         else !RR

          if(resultsb(k)%idx <=Ndata) then !DRR
!           if ( wgt ) then
!               Zdrr(ind)=Zdrr(ind)+(wgt1(resultsb(j)%idx)*wgt1(resultsb(k)%idx)*wgt1(i))
!            else
!               Zdrr(ind)=Zdrr(ind)+1.d0
!            endif
            continue
         else !RRR
            if ( wgt ) then
               Zrrr(ind)=Zrrr(ind)+(wgt1(resultsb(j)%idx)*wgt1(resultsb(k)%idx)*wgt1(i))
             else
               Zrrr(ind)=Zrrr(ind)+1.d0
            endif
         endif
       endif
34 continue
      enddo
   enddo

   call kdtree2_r_nearest_around_point(tp=tree,idxin=i,correltime=-1,r2=bins(3,2)*bins(3,2),&
&       nfound=nn2,nalloc=(Nrand+Ndata),results=resultsb)
   do k=nbins+2,1,-1
      nn1=kdtree2_r_count_around_point(tp=tree,idxin=i,correltime=-1,r2=bins(k,1)*bins(k,1))
      nn2=kdtree2_r_count_around_point(tp=tree,idxin=i,correltime=-1,r2=bins(k,2)*bins(k,2))
      do j=nn1+1,nn2,1
         if(resultsb(j)%idx >Ndata) then
            if ( wgt ) then
               Zrr(k,1)=Zrr(k,1)+(wgt1(resultsb(j)%idx)*wgt1(i))
            else
               Zrr(k,1)=Zrr(k,1)+1.d0
            endif

         else
            continue
!            if ( wgt ) then
!               Zdr(k,1)=Zdr(k,1)+(wgt1(resultsb(j)%idx)*wgt1(i))
!            else
!               Zdr(k,1)=Zdr(k,1)+1.d0
!            endif

         endif
      enddo
   enddo

enddo
print*,'finished RR loop on thread:', myid

55 continue

#ifdef MPI
call MPI_Barrier(MPI_COMM_WORLD,ierr)
call mpi_collect()
#endif

!call normalise_counts()
!call load_save_randoms()

#ifdef MPI
call MPI_Barrier(MPI_COMM_WORLD,ierr)
#endif

if(myid==master) then

call load_save_randoms()

do i=1,nbins+2
!   print*,bins(i,1),bins(i,2),Zdd(i,1),Zdr(i,1),Zrr(i,1), &
!& Zdd(i,1)/Zrr(i,1) - 1,Zdd(i,1)/Zdr(i,1) - 1, (Zdd(i,1)-2.0*Zdr(i,1)+Zrr(i,1))/Zrr(i,1)

  crap2(i,1)=(Zdd(i,1)-2.0*Zdr(i,1)+Zrr(i,1))/Zrr(i,1)
enddo

open(11,file=outfile,status='unknown')
write(11,*)'# theta_min, theta_max, r3_min, r3_max, r1_min, r1_max, r2_min, r2_max, DDD, DDR, DRR, &
& RRR, DD, DR, RR, xi_1, xi_2, xi_3, zeta, Q'
do i=nbins,1,-1
    crap(i)=(Zddd(i)-3.0*Zddr(i)+3.0*Zdrr(i) - Zrrr(i))/Zrrr(i)

  print*,(bins(2,1)-bins(1,2))+(dr*i),bins(1,1),bins(2,1),Zddd(i),Zddr(i),Zdrr(i), Zrrr(i),&
  &Zdd(i+2,1),Zdr(i+2,1),Zrr(i+2,1),crap(i)/(crap2(1,1)*crap2(2,1) + crap2(1,1)*crap2(2+i,1) + crap2(2,1)*crap2(2+i,1))
!  write(11,'(15(e14.7,1x))')(bins(i+2,2)+bins(i+2,1))/2.,bins(1,1),bins(2,1),Zddd(i),Zddr(i),Zdrr(i), Zrrr(i),&
!  &Zdd(i+2,1),Zdr(i+2,1),Zrr(i+2,1),crap2(1,1),crap2(2,1),crap2(2+i,1),crap(i),crap(i)/(crap2(1,1)*crap2(2,1) &
!  &+ crap2(1,1)*crap2(2+i,1) + crap2(2,1)*crap2(2+i,1))

  write(11,'(20(e12.5,1x))')acos(tbins(i,1)),acos(tbins(i,2)),bins(i+2,1),bins(i+2,2),bins(1,1),bins(1,2),&
  &bins(2,1),bins(2,2),Zddd(i),Zddr(i),Zdrr(i), Zrrr(i),&
  &Zdd(i+2,1),Zdr(i+2,1),Zrr(i+2,1),crap2(1,1),crap2(2,1),crap2(2+i,1),crap(i),crap(i)/(crap2(1,1)*crap2(2,1) &
  &+ crap2(1,1)*crap2(2+i,1) + crap2(2,1)*crap2(2+i,1))

enddo
close(11)
endif

! deallocate the memory for the data.
call kdtree2_destroy(tree)  

if(wgt) then
  deallocate(wgt1)
endif

call deallocate_arrays()

#ifdef MPI
call MPI_Barrier(MPI_COMM_WORLD,ierr)
call MPI_FINALIZE( ierr )
#endif

if(myid==master) stop "Exit... stage left!"

contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine parseOptions ()
     USE extension, ONLY : getArgument, nArguments
      implicit none

      logical lexist
      integer(kind=4)    :: i,n
      character(len=6)   :: opt
      character(len=100) :: arg

      n = nArguments()

      if (n < 8 .and. myid==master) then
        print*,' '
        print*,'Not enough input parameters. Please read the following help info'
        print*,' '
        call print_help
        stop
      end if

      do i = 1,n,2
         call getArgument(i,opt)
         if (i == n) then
            print '("option ",a2," has no argument")', opt
!            call fatal_error
         end if
         call getArgument(i+1,arg)
         select case (opt)
         case ('-gal')
            file1 = trim(arg)
         case ('-ran')
            file2 = trim(arg)
         case ('-out')
            outfile = trim(arg)
         case ('-r1min')
            read (arg,*) r1min
         case ('-r1max')
            read (arg,*) r1max
         case ('-r2min')
            read (arg,*) r2min
         case ('-r2max')
            read (arg,*) r2max
         case ('-nbins')
            read (arg,*) nbins
         case ('-wgt')
            wgt=.true.     
         case('-RR')        
            ranfile = trim(arg)
            inquire (file=ranfile,exist=lexist)
            if(lexist) then
               write(*,*) " Random count file exists!"
               loadran=.true.
            else
               write(*,*) " Random count file does not exist. It will be created!"
               saveran=.true.
            end if     
         case ('-help')
            call print_help
      stop

         case default
            print '("unknown option ",a6," ignored")', opt
            stop
         end select

      end do

   end subroutine parseOptions
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine print_help
print*,'PURPOSE: Code for calculating the 3PCF of a 3D point set. '
print*,'     '
print*,'     '
print*,'CALLING SEQUENCE:'
print*,'      3pcf [-gal gal_file] [-ran ran_file] [-out out_file] [-wgt WGT]'
print*,'           [-r1min Rmin] [-r1max Rmax] [-r2min Rmin] [-r2max Rmax]'
print*,'           [-nbins Nbins] [-RR RR_Counts]'
print*,' '
print*,'      eg: 3pcf -gal dr9.gal -ran dr9.ran -r1min 10.0 -r1max 12.0 -r2min 5.0 -r2max 6.0 -nbins 10'
print*,' '
print*,'INPUTS:'
print*,'       gal_file/ranfile = strings containing the name of a 3/4 column point file'
print*,'                 lines 1-3: XYZ comoving, optional 4: weight  '
print*,'   '
print*,'       out_file = filename for output correlation func. '
print*,'   '
print*,'       Rmin = Min. radial seperation'
print*,'   '
print*,'       Rmax = Max. radial seperation'
print*,' '
print*,'OPTIONAL:'
print*,'     Nbins = Number of radial bins to calculate'
print*,' '
print*,'       WGT = logical value to define if weights are to be included '
print*,'              In this case, the input gal/ran files should be 4 columns '
print*,' '
print*,' RR_Counts = string containing the name of a file regarding RR data counts.  '
print*,'             If file exists, then the RR counting will be skipped (saveing CPU time)'
print*,'             If it does not exist, the file will be created for future use.'
print*,'             Consider using this option if you are doing Monte Carlo runs with same random catalogue.' 
print*,' '
   end subroutine print_help
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine count_files ()
    implicit none

  Ndata=0
  open(11,file=file1,status='unknown')
12 continue
  read(11,*,err=12,end=13)aux
  Ndata=Ndata+1
  goto 12
13 close(11)
  if(myid==0) print*,'Preparing to read ',Ndata, 'data points'

  Nrand=0
  open(11,file=file2,status='unknown')
21 continue
  read(11,*,err=21,end=22)aux
  Nrand=Nrand+1
  goto 21
22 close(11)
  if(myid==0) print*,'Preparing to read ',Nrand, 'random points'

end subroutine count_files

  subroutine read_files ()
    implicit none

  if(myid==0)  print*, 'opening ',file1
  open(11,file=file1,status='unknown')
  do i=1,Ndata
    if ( wgt ) then
     read(11,*,end=14)my_array(1:3,i),wgt1(i)
    else
     read(11,*,end=14)my_array(1:3,i)
   endif

  enddo
14 close(11)
  if(myid==0)  print*,'Finished reading data file 1'  


  if(myid==0)  print*, 'opening ',file2
  open(11,file=file2,status='unknown')
  do i=Ndata+1,Ndata+Nrand
    if ( wgt ) then
     read(11,*,end=23)my_array(1:3,i),wgt1(i)
   else
     read(11,*,end=23)my_array(1:3,i)
   endif

  enddo
23 close(11)
  if(myid==0)  print*,'Finished reading data file 2'  

if(wgt) then
    wgt1(1:Ndata)=wgt1(1:Ndata)/sum(wgt1(1:Ndata))
    wgt1(Ndata+1:Ndata+Nrand)=wgt1(Ndata+1:Ndata+Nrand)/sum(wgt1(Ndata+1:Ndata+Nrand))
endif

end subroutine read_files
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine allocate_arrays ()
  implicit none
  allocate(resultsb(Ndata+Nrand))

  allocate(Zddd(nbins))
  allocate(Zddr(nbins))
  allocate(Zdrr(nbins))
  allocate(Zrrr(nbins))
  allocate(crap(nbins))

  crap(:)=0.0d0
  Zddd(:)=0.0d0
  Zddr(:)=0.0d0
  Zdrr(:)=0.0d0
  Zrrr(:)=0.0d0

  allocate(Zdd(nbins+2,1))
  allocate(Zdr(nbins+2,1))
  allocate(Zrr(nbins+2,1))
  allocate(crap2(nbins+2,1))

  crap2(:,:)=0d0
  Zdd(:,:)=0d0
  Zdr(:,:)=0d0
  Zrr(:,:)=0d0


end subroutine allocate_arrays
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine mpi_collect()

  if(myid==master) crap(:)=0.d0
call MPI_Barrier(MPI_COMM_WORLD,ierr)
call MPI_REDUCE( Zddd, crap, nbins, MPI_DOUBLE_PRECISION,MPI_SUM, master, MPI_COMM_WORLD, ierr )
if ( myid == master ) then !in master thread
Zddd=crap
endif
call MPI_Barrier(MPI_COMM_WORLD,ierr)

if(myid==master) crap(:)=0.d0
call MPI_Barrier(MPI_COMM_WORLD,ierr)
call MPI_REDUCE( Zddr, crap, nbins, MPI_DOUBLE_PRECISION,MPI_SUM, master, MPI_COMM_WORLD, ierr )
if ( myid == master ) then !in master thread
Zddr=crap
endif
call MPI_Barrier(MPI_COMM_WORLD,ierr)

if(myid==master) crap(:)=0.d0
call MPI_Barrier(MPI_COMM_WORLD,ierr)
call MPI_REDUCE( Zdrr, crap, nbins, MPI_DOUBLE_PRECISION,MPI_SUM, master, MPI_COMM_WORLD, ierr )
if ( myid == master ) then !in master thread
Zdrr=crap
endif
call MPI_Barrier(MPI_COMM_WORLD,ierr)

if(myid==master) crap(:)=0.d0
call MPI_Barrier(MPI_COMM_WORLD,ierr)
call MPI_REDUCE( Zrrr, crap, nbins, MPI_DOUBLE_PRECISION,MPI_SUM, master, MPI_COMM_WORLD, ierr )
if ( myid == master ) then !in master thread
Zrrr=crap
endif
call MPI_Barrier(MPI_COMM_WORLD,ierr)

if(myid==master) crap2(:,:)=0.d0
call MPI_Barrier(MPI_COMM_WORLD,ierr)
call MPI_REDUCE( Zdd, crap2, nbins+2, MPI_DOUBLE_PRECISION,MPI_SUM, master, MPI_COMM_WORLD, ierr )
if ( myid == master ) then !in master thread
Zdd=crap2
endif
call MPI_Barrier(MPI_COMM_WORLD,ierr)

if(myid==master) crap2(:,:)=0.d0
call MPI_Barrier(MPI_COMM_WORLD,ierr)
call MPI_REDUCE( Zdr, crap2, nbins+2, MPI_DOUBLE_PRECISION,MPI_SUM, master, MPI_COMM_WORLD, ierr )
if ( myid == master ) then !in master thread
Zdr=crap2
endif
call MPI_Barrier(MPI_COMM_WORLD,ierr)

if(myid==master) crap2(:,:)=0.d0
call MPI_Barrier(MPI_COMM_WORLD,ierr)
call MPI_REDUCE( Zrr, crap2, nbins+2, MPI_DOUBLE_PRECISION,MPI_SUM, master, MPI_COMM_WORLD, ierr )
if ( myid == master ) then !in master thread
Zrr=crap2
endif
call MPI_Barrier(MPI_COMM_WORLD,ierr)

end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine normalise_counts()
  if ( wgt ) then
   ndd=sum(wgt1(1:Ndata))
   nrr=sum(wgt1(Ndata+1:Ndata+Nrand))
   Zddd=Zddd/(ndd*ndd*ndd)
   Zddr=Zddr/(2*ndd*ndd*nrr)
   Zdrr=Zdrr/(1*ndd*nrr*nrr)
   Zrrr=Zrrr/(nrr*nrr*nrr)
else
   ndd=float(Ndata)
   nrr=float(Nrand)
   Zddd=Zddd/(ndd*ndd*ndd)
   Zddr=Zddr/(2*ndd*ndd*nrr)
   Zdrr=Zdrr/(1*ndd*nrr*nrr)
   Zrrr=Zrrr/(nrr*nrr*nrr)
endif

if ( wgt ) then
   nrr=sum(wgt1(Ndata+1:Ndata+Nrand))*(sum(wgt1(Ndata+1:Ndata+Nrand))-sum(wgt1(Ndata+1:Ndata+Nrand))/float(Nrand))
   ndd=sum(wgt1(1:Ndata))*(sum(wgt1(1:Ndata))-sum(wgt1(1:Ndata))/float(Ndata))
   ndr=sum(wgt1(1:Ndata))*sum(wgt1(Ndata+1:Ndata+Nrand))

   Zdd=Zdd/ndd
   Zrr=Zrr/nrr
   Zdr=Zdr/ndr
else
  Zdd=Zdd/(float(Ndata)*Ndata)
  Zrr=Zrr/(float(Nrand)*Nrand)
  Zdr=Zdr/(1.0*float(Ndata)*Nrand)
endif

end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine load_save_randoms()
if(saveran) then
  open(11,file=ranfile,status='unknown')
    write(11,*)Zrr(1,1)
    write(11,*)Zrr(2,1)
do i=1,nbins
    write(11,*)Zrr(i+2,1), Zrrr(i)
  enddo
  close(11)
endif

if(loadran) then
  open(11,file=ranfile,status='unknown')
    read(11,*)Zrr(1,1)
    read(11,*)Zrr(2,1)
  do i=1,nbins
    read(11,*)Zrr(i+2,1), Zrrr(i)
  enddo
  close(11)
endif

end subroutine load_save_randoms
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine triplet_bin( )

if (.false.) then !distance binning
    dist=sqrt((v1(1)-v2(1))**2 + (v1(2)-v2(2))**2 + (v1(3)-v2(3))**2)
    ind=floor((dist-bins(3,1))/dr)+1
endif

if (.true.) then !theta binning
!    v1=[v1(1)-my_array(1,i),v1(2)-my_array(2,i),v1(3)-my_array(3,i)] ! connecting vector
    v2=[v2(1)-my_array(1,i),v2(2)-my_array(2,i),v2(3)-my_array(3,i)] ! connecting vector

    theta = DOT_PRODUCT(v1,v2)/(SQRT(DOT_PRODUCT(v1,v1))*SQRT(DOT_PRODUCT(v2,v2)))

    ind=ceiling((theta+1.0)*odt)
!    print*,acos(theta)*180/3.14159,theta,ind

!l=1
!51  continue
!if(theta<tbins(l,1)) then
!    ind=l
!else
!    l=l+1
!    goto 51
!endif

endif

if(ind==0 .or. ind>nbins) print*,'tuple index out of bounds: ',ind, theta

end subroutine 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine make_bins()

  mid1=bins(1,2)!(bins(1,1)+bins(1,2))/2.
  mid2=bins(2,2)!(bins(2,1)+bins(2,2))/2.

  dr=(bins(2,2) + bins(1,2))
  dr=dr - (bins(1,1)-bins(2,2))

  dr=dr/float(nbins)
  odr=1./dr

  odt=1./(2.0/nbins)
  print*,dr

  do i=1,nbins
!!    bins(2+i,1)=abs(mid1-mid2)+((i-1)*dr)
!!    bins(2+i,2)=bins(2+i,1)+dr
!    bins(2+i,1)=bins(1,1)-bins(2,2)+((i-1)*dr)
!    bins(2+i,2)=bins(2+i,1)+dr
    
!    bins(2+i,1)=sqrt(mid1**2 + mid2**2 - 2.*mid1*mid2*cos(-1.0 + (i-1)/odt))
!    bins(2+i,2)=sqrt(mid1**2 + mid2**2 - 2.*mid1*mid2*cos(-1.0 + i/odt))
    bins(2+i,2)=sqrt(mid1**2 + mid2**2 - 2.*mid1*mid2*(-1+(i-1)/odt))
    bins(2+i,1)=sqrt(mid1**2 + mid2**2 - 2.*mid1*mid2*(-1+i/odt))
    tbins(i,1)=(mid1**2 + mid2**2 - bins(2+i,1)**2)/(2.*mid1*mid2)
    tbins(i,2)=(mid1**2 + mid2**2 - bins(2+i,2)**2)/(2.*mid1*mid2)
    
    if(myid==master) print*,bins(2+i,:),tbins(i,:), cos((i-1)/odt), cos(i/odt)
  enddo
  end subroutine make_bins
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine deallocate_arrays()
implicit none
if (allocated(my_array)) deallocate(my_array)
if (allocated(Zdd)) deallocate(Zdd)
if (allocated(Zdr)) deallocate(Zdr)
if (allocated(Zrr)) deallocate(Zrr)
if (allocated(crap)) deallocate(crap)
if (allocated(wgt1)) deallocate(wgt1)
if (allocated(Zddd)) deallocate(Zddd)
if (allocated(Zddr)) deallocate(Zddr)
if (allocated(Zdrr)) deallocate(Zdrr)
if (allocated(Zrrr)) deallocate(Zrrr)
if (allocated(bins)) deallocate(bins)
if (allocated(v1)) deallocate(v1)
if (allocated(v2)) deallocate(v2)
end subroutine deallocate_arrays
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end program  ThreePCF

