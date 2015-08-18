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
  type(kdtree2_result),allocatable :: resultsb(:), results2(:)
  integer   ::  Ndata,Nrand,nn1,nn2,nnpairs,mm1,mm2
  double precision  :: aux,rmin,rmax,dr,odr,theta,odt,mid1,mid2
  logical :: wgt
  CHARACTER(LEN=100) :: outfile
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

  if(myid==master) print*,'reading options'
  call parseOptions()

 allocate(bins(nbins+2,2))
 allocate(tbins(nbins,2))

!   bins(2,1)=r1min
!   bins(2,2)=r1max
!   bins(1,1)=r2min
!   bins(1,2)=r2max

if(logbins) then
  dr=(log10(rmax)-log10(rmin))/nbins
else
  dr=(rmax-rmin)/nbins
endif

  do i=1,nbins
    if(logbins) then
      bins(i,1)=10.D0**(log10(rmin) + dr*(float(i)-1.0))
      bins(i,2)=10.D0**(log10(rmin) + (dr*i))
    else
      bins(i,1)=rmin + dr*(float(i)-1.0)
      bins(i,2)=rmin + dr*i
    endif
    if(myid==0) print*,bins(i,:)
  enddo

   if(myid==master)  print*,'allocated bins'

!  if(bins(2,1)<bins(1,1)) then
!    print*,'side 1 must be larger than side 2'
!    stop "exiting"
!  endif


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

!  mid1=(bins(1,1)+bins(1,2))/2.
!  mid2=(bins(2,1)+bins(2,2))/2.

!  dr=(bins(2,2) + bins(1,2))
!  dr=dr - (bins(2,1)-bins(1,2))

!  dr=mid1+mid2
!  dr=dr-abs(mid1-mid2)

!  dr=dr/nbins
!  odr=1./dr
!  odt=1./(3.14159/nbins)
!  odt=1./(2.0/nbins)
!  print*,dr

!  do i=1,nbins
!    bins(2+i,1)=abs(mid1-mid2)+((i-1)*dr)
!    bins(2+i,2)=bins(2+i,1)+dr
!    bins(2+i,1)=sqrt(mid1**2 + mid2**2 - 2.*mid1*mid2*cos((i-1)/odt))
!    bins(2+i,2)=sqrt(mid1**2 + mid2**2 - 2.*mid1*mid2*cos(i/odt))
!    tbins(i,1)=(mid1**2 + mid2**2 - bins(2+i,1)**2)/(2.*mid1*mid2)
!    tbins(i,2)=(mid1**2 + mid2**2 - bins(2+i,2)**2)/(2.*mid1*mid2)
!    if(myid==master) print*,bins(2+i,:),tbins(i,:)
!  enddo


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
!      mm1=kdtree2_r_count_around_point(tp=tree,idxin=i,correltime=-1,r2=bins(1,1)*bins(1,1))
!      mm2=kdtree2_r_count_around_point(tp=tree,idxin=i,correltime=-1,r2=bins(1,2)*bins(1,2))

      do j=nn1+1,nn2,1
        v1=my_array(:,resultsb(j)%idx)
        v1=[v1(1)-my_array(1,i),v1(2)-my_array(2,i),v1(3)-my_array(3,i)] ! connecting vector

        do k=mm1+1,mm2,1
          v2=my_array(:,resultsb(k)%idx)
          v2=[v2(1)-my_array(1,i),v2(2)-my_array(2,i),v2(3)-my_array(3,i)] ! connecting vector
          theta=v1(1)*v2(1) + v1(2)*v2(2) + v1(3)*v2(3)
          theta=theta/((sqrt(v1(1)**2 + v1(2)**2 + v1(3)**2)) * (sqrt(v2(1)**2 +v2(2)**2 +v2(3)**2)))
!          theta=mod(theta/((sqrt(v1(1)**2 + v1(2)**2 + v1(3)**2)) * (sqrt(v2(1)**2 +v2(2)**2 +v2(3)**2))),1.0)

          l=1
51          continue
          if(theta>tbins(l,2)) then
            ind=l
          else
            l=l+1
            goto 51
          endif


!          ind=floor(acos(theta)*odt)+1
!          print*,theta,acos(theta),asin(theta)
!          dist=sqrt((v1(1)-v2(1))**2 + (v1(2)-v2(2))**2 + (v1(3)-v2(3))**2)
!          dist=dist-(bins(2,1)-bins(1,2))
!          ind=floor(dist*odr)+1

          if(resultsb(j)%idx <=Ndata) then 
          
            if(resultsb(k)%idx <=Ndata) then

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

              if ( wgt ) then
                Zddr(ind)=Zddr(ind)+(wgt1(resultsb(j)%idx)*wgt1(resultsb(k)%idx)*wgt1(i))
              else
                Zddr(ind)=Zddr(ind)+1.d0
              endif

            else

              if ( wgt ) then
                Zdrr(ind)=Zdrr(ind)+(wgt1(resultsb(j)%idx)*wgt1(resultsb(k)%idx)*wgt1(i))
              else
                Zdrr(ind)=Zdrr(ind)+1.d0
              endif

            endif
          endif

      enddo
   enddo


   call kdtree2_r_nearest_around_point(tp=tree,idxin=i,correltime=-1,r2=bins(nbins+2,2)*bins(nbins+2,2),&
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

        do k=mm1+1,mm2,1
          v2=my_array(:,resultsb(k)%idx)
          v2=[v2(1)-my_array(1,i),v2(2)-my_array(2,i),v2(3)-my_array(3,i)] ! connecting vector
          theta=v1(1)*v2(1) + v1(2)*v2(2) + v1(3)*v2(3)
          theta=theta/((sqrt(v1(1)**2 + v1(2)**2 + v1(3)**2)) * (sqrt(v2(1)**2 +v2(2)**2 +v2(3)**2)))
!          theta=mod(theta/((sqrt(v1(1)**2 + v1(2)**2 + v1(3)**2)) * (sqrt(v2(1)**2 +v2(2)**2 +v2(3)**2))),1.0)

          l=1
52        continue
          if(theta>tbins(l,2)) then
            ind=l
          else
            l=l+1
            goto 52
          endif

    !          ind=floor(acos(theta)*odt)+1
    !          dist=sqrt((v1(1)-v2(1))**2 + (v1(2)-v2(2))**2 + (v1(3)-v2(3))**2)
    !          dist=dist-(bins(2,1)-bins(1,2))
    !          ind=floor(dist*odr)+1

         if(resultsb(j)%idx <=Ndata) then 
          if(resultsb(k)%idx <=Ndata) then 
            if ( wgt ) then
               Zddr(ind)=Zddr(ind)+(wgt1(resultsb(j)%idx)*wgt1(resultsb(k)%idx)*wgt1(i))
            else
               Zddr(ind)=Zddr(ind)+1.d0
            endif
         else
            if ( wgt ) then
               Zdrr(ind)=Zdrr(ind)+(wgt1(resultsb(j)%idx)*wgt1(resultsb(k)%idx)*wgt1(i))
            else
               Zdrr(ind)=Zdrr(ind)+1.d0
            endif
          endif
            
         else 

          if(resultsb(k)%idx <=Ndata) then 
           if ( wgt ) then
               Zdrr(ind)=Zdrr(ind)+(wgt1(resultsb(j)%idx)*wgt1(resultsb(k)%idx)*wgt1(i))
            else
               Zdrr(ind)=Zdrr(ind)+1.d0
            endif
         else
            if ( wgt ) then
               Zrrr(ind)=Zrrr(ind)+(wgt1(resultsb(j)%idx)*wgt1(resultsb(k)%idx)*wgt1(i))
            else
               Zrrr(ind)=Zrrr(ind)+1.d0
            endif
         endif
       endif
      enddo
   enddo

   call kdtree2_r_nearest_around_point(tp=tree,idxin=i,correltime=-1,r2=bins(nbins+2,2)*bins(nbins+2,2),&
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

            if ( wgt ) then
               Zdr(k,1)=Zdr(k,1)+(wgt1(resultsb(j)%idx)*wgt1(i))
            else
               Zdr(k,1)=Zdr(k,1)+1.d0
            endif

         endif
      enddo
   enddo

enddo
print*,'finished RR loop on thread:', myid

#ifdef MPI
call mpi_collect()
#endif

if(myid==master) then
do i=1,nbins
   print*,(bins(2,1)-bins(1,2))+(dr*i), bins(1,1),bins(2,1),Zddd(i),Zddr(i),Zdrr(i), Zrrr(i)
enddo
endif

call normalise_counts()

if(myid==master) then
do i=1,nbins+2
!   print*,bins(i,1),bins(i,2),Zdd(i,1),Zdr(i,1),Zrr(i,1), &
!& Zdd(i,1)/Zrr(i,1) - 1,Zdd(i,1)/Zdr(i,1) - 1, (Zdd(i,1)-2.0*Zdr(i,1)+Zrr(i,1))/Zrr(i,1)

  crap2(i,1)=(Zdd(i,1)-2.0*Zdr(i,1)+Zrr(i,1))/Zrr(i,1)

enddo
endif


if(myid==master) then
open(11,file=outfile,status='unknown')
do i=1,nbins

  crap(i)=(Zddd(i)-3.0*Zddr(i)+3.0*Zdrr(i) - Zrrr(i))/Zrrr(i)
  print*,(bins(2,1)-bins(1,2))+(dr*i),bins(1,1),bins(2,1),Zddd(i),Zddr(i),Zdrr(i), Zrrr(i),&
  &Zdd(i+2,1),Zdr(i+2,1),Zrr(i+2,1),crap(i)/(crap2(1,1)*crap2(2,1) + crap2(1,1)*crap2(2+i,1) + crap2(2,1)*crap2(2+i,1))
  write(11,'(11(e14.7,1x))')(bins(2,1)-bins(1,2))+(dr*i),bins(1,1),bins(2,1),Zddd(i),Zddr(i),Zdrr(i), Zrrr(i),&
  &Zdd(i+2,1),Zdr(i+2,1),Zrr(i+2,1),crap(i)/(crap2(1,1)*crap2(2,1) + crap2(1,1)*crap2(2+i,1) + crap2(2,1)*crap2(2+i,1))

enddo
close(11)
endif

! deallocate the memory for the data.
call kdtree2_destroy(tree)  

if(wgt) then
  deallocate(wgt1)
endif

deallocate(my_array,Zddd,Zddr,Zdrr,Zrrr,crap,v1,v2)

#ifdef MPI
call MPI_Barrier(MPI_COMM_WORLD,ierr)
call MPI_FINALIZE( ierr )
#endif

if(myid==master) stop "Exit... stage left!"

contains

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine parseOptions ()
     USE extension, ONLY : getArgument, nArguments
      implicit none

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
         case ('-help')
            call print_help
      stop

         case default
            print '("unknown option ",a6," ignored")', opt
            stop
         end select

      end do

   end subroutine parseOptions
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine print_help
print*,'PURPOSE: Code for calculating the 3PCF of a 3D point set. '
print*,'     '
print*,'     '
print*,'CALLING SEQUENCE:'
print*,'      3pcf [-gal gal_file] [-ran ran_file] [-out out_file] [-wgt WGT]'
print*,'           [-r1min Rmin] [-r1max Rmax] [-r2min Rmin] [-r2max Rmax]'
print*,'           [-nbins Nbins] [-bins bin_file] '
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
print*,'       Nbins = Number of radial bins to calculate'
print*,' '
print*,'       WGT = logical value to define if weights are to be included '
print*,'              In this case, the input gal/ran files should be 4 columns '
print*,' '
print*,' '
!print*,'       bin_file = string containing the name of the user defined radial bins.'
!print*,'                  If the file exist, any supplied Rmin/Rmax/Nbins values will be ignored '
!print*,'                  The format should be 2 column of Lower and upper radial values'
!print*,' '
!print*,' '
   end subroutine print_help
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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

end subroutine read_files

subroutine allocate_arrays ()
  implicit none
  allocate(resultsb(Ndata+Nrand))
  allocate(results2(Ndata+Nrand))

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

subroutine normalise_counts()
  if ( wgt ) then
   ndd=sum(wgt1(1:Ndata))
   nrr=sum(wgt1(Ndata+1:Ndata+Nrand))
   Zddd=Zddd/(ndd*ndd*ndd)
   Zddr=Zddr/(3*ndd*ndd*nrr)
   Zdrr=Zdrr/(3*ndd*nrr*nrr)
   Zrrr=Zrrr/(nrr*nrr*nrr)
else
   ndd=float(Ndata)
   nrr=float(Nrand)
   Zddd=Zddd/(ndd*ndd*ndd)
   Zddr=Zddr/(3*ndd*ndd*nrr)
   Zdrr=Zdrr/(3*ndd*nrr*nrr)
   Zrrr=Zrrr/(nrr*nrr*nrr)
endif

if ( wgt ) then
   ndd=sum(wgt1(1:Ndata))**2
   nrr=sum(wgt1(Ndata+1:Ndata+Nrand))**2
   ndr=2.0*sqrt(ndd)*sqrt(nrr)

   Zdd=Zdd/ndd
   Zrr=Zrr/nrr
   Zdr=Zdr/ndr
else
  Zdd=Zdd/(float(Ndata)*Ndata)
  Zrr=Zrr/(float(Nrand)*Nrand)
  Zdr=Zdr/(2.0*float(Ndata)*Nrand)
endif

end subroutine

end program  ThreePCF
