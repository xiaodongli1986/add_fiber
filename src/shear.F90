program shear
  use kdtree2_module
  use kdtree2_precision_module
  implicit none
#ifdef MPI
	include 'mpif.h'
#endif

  double precision, dimension(:,:), allocatable :: my_array
  double precision, allocatable :: bins(:,:)

  type(kdtree2), pointer :: tree
  ! this is how you declare a tree in your main program
  double precision, allocatable :: Zdd(:,:),Zdr(:,:),Zrr(:,:),crap(:,:), wgt1(:),gamma(:,:)
  double precision, allocatable :: dd(:,:,:),drn(:,:,:),rr(:,:,:),xi(:,:,:), v1(:),v2(:)
  integer, allocatable :: boot(:)
  double precision :: ndd,ndr,nrr,sig,pi,mind,dist1,dist2
  integer :: k, i, j, d,chunk,nbins,ntbin,Nresample
  character :: file1*100, file2*100, outfile*100, binfile*100
  type(kdtree2_result),allocatable :: resultsb(:), results2(:)
  integer   ::  Ndata,Nrand,nn1,nn2,nnpairs,sigbin,pibin,tbin
  double precision :: aux,rmin,rmax,dr,odr,theta
  logical :: wgt,logbins,userbin,proj,resample,dosigpi,dotpcf,doap
  CHARACTER(LEN=100) :: arg,resample_method,ISOFLAG,DECOMP
  integer :: myid , ntasks , ierr , islave
  double precision, parameter :: pival=3.14159265358979323846

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
  call default_params()

!--- acquire input parameters ---
  call parseOptions()

  call count_files()

  Allocate(my_array(d,Ndata+Nrand))
  allocate(gamma(Ndata,2))
  if ( wgt ) then
    Allocate(wgt1(Ndata+Nrand))
  endif
  allocate(bins(nbins,2))

  call read_files()

print*,' building tree on node', myid
tree => kdtree2_create(my_array,sort=.true.,rearrange=.true.)     ! this is how you create a tree.
print*,' built tree on node', myid


if(ISOFLAG=='ISO') dotpcf=.true.
if(ISOFLAG=='ANISO') then
  if(DECOMP=='SIGPI') then
    dosigpi=.true.
  elseif(DECOMP=='SMU') then
      doap=.true.
  else
      stop "No valid decomposition method, I'm done."
  endif
endif


if(dotpcf) call tpcf_iso()


stop "Ciao for now!"

contains

!------------------------------------------------------------
!------------------------------------------------------------

subroutine tpcf_iso ()
if(proj) then
  print*,'Calculating the isotropic 2D correlation function'
else
  print*,'Calculating the isotropic 3D correlation function'
endif

call allocate_arrays_2pcf()

if(resample) call bootstrapper()

if(proj) then
  rmin=sqrt(2.d0*(1.d0 - cos(rmin*pival/180.d0)))
  rmax=sqrt(2.d0*(1.d0 - cos(rmax*pival/180.d0)))
endif

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


print*,'beginning DD loop on thread:',myid

#ifdef MPI
chunk=floor(float(Ndata)/ntasks)+1
do i=max(myid*chunk,1),min(((myid+1)*chunk)-1,Ndata),1
#else
do i=1,Ndata,1
#endif

   call kdtree2_r_nearest_around_point(tp=tree,idxin=i,correltime=-1,r2=bins(nbins,2)*bins(nbins,2),&
&       nfound=nn2,nalloc=(Nrand+Ndata),results=resultsb)
   
   do k=nbins,1,-1
      nn1=kdtree2_r_count_around_point(tp=tree,idxin=i,correltime=-1,r2=bins(k,1)*bins(k,1))
      do j=nn1+1,nn2,1
         if(resultsb(j)%idx <=Ndata) then            
            if ( wgt ) then
               if(resample) then
                  Zdd(k,boot(resultsb(j)%idx))=Zdd(k,boot(resultsb(j)%idx))+(wgt1(resultsb(j)%idx)*wgt1(i))
                  Zdd(k,boot(i))=Zdd(k,boot(i))+(wgt1(resultsb(j)%idx)*wgt1(i))
               else
                  Zdd(k,1)=Zdd(k,1)+(gamma(resultsb(j)%idx,1)*gamma(i,1))
                  Zdr(k,1)=Zdr(k,1)+(gamma(resultsb(j)%idx,2)*gamma(i,2))
!                  Zdd(k,1)=Zdd(k,1)+(gamma(resultsb(j)%idx,1)*gamma(i,1))+(gamma(resultsb(j)%idx,2)*gamma(i,2))!)/(1.d0*(nn2-nn1))
!                  Zdr(k,1)=Zdr(k,1)+(gamma(resultsb(j)%idx,1)*gamma(i,1))-(gamma(resultsb(j)%idx,2)*gamma(i,2))!)/(1.d0*(nn2-nn1))
                  Zrr(k,1)=Zrr(k,1)+1.d0*(nn2-nn1-1)
                  !(gamma(resultsb(j)%idx,1)*gamma(i,2))!/1.d0*(nn2-nn1)
               endif
!             else
!                if(resample) then
!                   Zdd(k,boot(resultsb(j)%idx))=Zdd(k,boot(resultsb(j)%idx))+1.d0
!                   Zdd(k,boot(i))=Zdd(k,boot(i))+1.d0
!                 else
!                   Zdd(k,1)=Zdd(k,1)+1.d0
!                endif
            endif
!          else
!             if ( wgt ) then
!                if(resample) then
!                   Zdr(k,boot(resultsb(j)%idx))=Zdr(k,boot(resultsb(j)%idx))+(wgt1(resultsb(j)%idx)*wgt1(i))
!                   Zdr(k,boot(i))=Zdr(k,boot(i))+(wgt1(resultsb(j)%idx)*wgt1(i))
!                 else
!                   Zdr(k,1)=Zdr(k,1)+(wgt1(resultsb(j)%idx)*wgt1(i))
!                endif
!             else
!                if(resample) then
!                   Zdr(k,boot(resultsb(j)%idx))=Zdr(k,boot(resultsb(j)%idx))+1.d0
!                   Zdr(k,boot(i))=Zdr(k,boot(i))+1.d0
!                 else
!                   Zdr(k,1)=Zdr(k,1)+1.d0
!                endif
!             endif
          endif
      enddo
      nn2=nn1

   enddo
enddo
print*,'finished DD loop in thread:', myid


! !RR loop
! print*,'beginning RR loop on thread:',myid
! #ifdef MPI
! chunk=floor(float(Nrand)/ntasks)+1
! do i=Ndata+max(myid*chunk,1),Ndata+min(((myid+1)*chunk)-1,Nrand),1
! #else
! do i=Ndata+1,Ndata+Nrand,1
! #endif
!    call kdtree2_r_nearest_around_point(tp=tree,idxin=i,correltime=-1,r2=bins(nbins,2)*bins(nbins,2),&
! &       nfound=nn2,nalloc=(Nrand+Ndata),results=resultsb)

!    do k=nbins,1,-1

!       nn1=kdtree2_r_count_around_point(tp=tree,idxin=i,correltime=-1,r2=bins(k,1)*bins(k,1))

!       do j=nn1+1,nn2,1
!          if(resultsb(j)%idx >Ndata) then
!             if ( wgt ) then
!                if(resample) then
!                   Zrr(k,boot(resultsb(j)%idx))=Zrr(k,boot(resultsb(j)%idx))+(wgt1(resultsb(j)%idx)*wgt1(i))
!                   Zrr(k,boot(i))=Zrr(k,boot(i))+(wgt1(resultsb(j)%idx)*wgt1(i))
!                 else
!                   Zrr(k,1)=Zrr(k,1)+(wgt1(resultsb(j)%idx)*wgt1(i))
!                endif
!             else
!                if(resample) then
!                   Zrr(k,boot(resultsb(j)%idx))=Zrr(k,boot(resultsb(j)%idx))+1.d0
!                   Zrr(k,boot(i))=Zrr(k,boot(i))+1.d0
!                else
!                   Zrr(k,1)=Zrr(k,1)+1.d0
!             endif
!           endif

!          else

!             if ( wgt ) then
!                 if(resample) then
!                   Zdr(k,boot(resultsb(j)%idx))=Zdr(k,boot(resultsb(j)%idx))+(wgt1(resultsb(j)%idx)*wgt1(i))
!                   Zdr(k,boot(i))=Zdr(k,boot(i))+(wgt1(resultsb(j)%idx)*wgt1(i))
!                 else
!                   Zdr(k,1)=Zdr(k,1)+(wgt1(resultsb(j)%idx)*wgt1(i))
!                 endif
!             else
!                if(resample) then
!                   Zdr(k,boot(resultsb(j)%idx))=Zdr(k,boot(resultsb(j)%idx))+1.d0
!                   Zdr(k,boot(i))=Zdr(k,boot(i))+1.d0
!                else
!                Zdr(k,1)=Zdr(k,1)+1.d0
!             endif

!          endif
!        endif
!       enddo
!       nn2=nn1

!    enddo
! enddo
! print*,'finished RR loop on thread:', myid

#ifdef MPI
call mpi_collect()
#endif

! if ( wgt ) then
!    ndd=sum(wgt1(1:Ndata))**2
!    nrr=sum(wgt1(Ndata+1:Ndata+Nrand))**2
!    ndr=2.0*sqrt(ndd)*sqrt(nrr)

!    Zdd=Zdd/ndd
!    Zrr=Zrr/ndd
!    Zdr=Zdr/ndd
! else

  Zdd=(Zdd+Zdr)/Zrr!(float(Ndata)*Ndata)
!  Zrr=Zrr/(float(Ndata)*Ndata)
  Zdr=(Zdd-Zdr)/Zrr!(float(Ndata)*Ndata)

! endif


!if(resample) then
! allocate(dd(Nresample,nbins,2))
! allocate(drn(Nresample,nbins,2))
! allocate(rr(Nresample,nbins,2))
! allocate(xi(Nresample,nbins,3))
! dd=0.d0
! drn=0.d0
! rr=0.d0

!  do i=1,Nresample

! ndd=0.d0
! ndr=0.d0
! nrr=0.d0

!   do j=1,Ndata
!    if(boot(j).ne.i) then
!      ndd=ndd+wgt1(j)
!    endif
!   enddo 
!   do j=Ndata+1,Ndata+Nrand
!    if(boot(j).ne.i) then
!      nrr=nrr+wgt1(j)
!    endif
!   enddo 

!   ndd=ndd**2
!   nrr=nrr**2
!   ndr=2.0*sqrt(ndd)*sqrt(nrr)

!   do j=1,Nresample
!     if(j.ne.i) then
!       do k=1,nbins
!         dd(i,k,1)=dd(i,k,1)+Zdd(k,j)
!         drn(i,k,1)=drn(i,k,1)+Zdr(k,j)
!         rr(i,k,1)=rr(i,k,1)+Zrr(k,j)
!       enddo
!     endif
!   enddo

!   do j=1,nbins
!     dd(i,j,2)=dd(i,j,1)/ndd
!     drn(i,j,2)=drn(i,j,1)/ndr
!     rr(i,j,2)=rr(i,j,1)/nrr

!     xi(i,j,1)=dd(i,j,2)/rr(i,j,2) -1.d0
!     xi(i,j,2)=dd(i,j,2)/drn(i,j,2) -1.d0
!     xi(i,j,3)=(dd(i,j,2)-2.*drn(i,j,2)+rr(i,j,2))/rr(i,j,2)
!   enddo
! enddo

! if(myid==master) then
! do i=1,nbins
!   print*,xi(1,i,1),xi(1,i,2),xi(1,i,3)
! enddo
! endif
! !endif



if(myid==master) then
open(11,file=outfile,status='unknown')
write(11,*)'R_min, R_max, DD, DR, RR, XI_natural, XI_Davis, XI_LS (most accurate)'
do i=1,nbins

!write(11,*) bins(i,1),bins(i,2), dd(1,i,2), drn(1,i,2),rr(1,i,2),xi(1,i,1),xi(1,i,2),xi(1,i,3)

   write(11,*) bins(i,1),bins(i,2),Zdd(i,1),Zdr(i,1),Zrr(i,1), &
& Zdd(i,1)/Zrr(i,1) - 1,Zdd(i,1)/Zdr(i,1) - 1, (Zdd(i,1)-2.0*Zdr(i,1)+Zrr(i,1))/Zrr(i,1)
   print*,bins(i,1),bins(i,2),Zdd(i,1),Zdr(i,1),Zrr(i,1), &
& Zdd(i,1)/Zrr(i,1) - 1,Zdd(i,1)/Zdr(i,1) - 1, (Zdd(i,1)-2.0*Zdr(i,1)+Zrr(i,1))/Zrr(i,1)

enddo
close(11)
endif

! deallocate the memory for the data.
call kdtree2_destroy(tree)  

if(wgt) then
  deallocate(wgt1)
endif

deallocate(my_array,Zdd,Zdr,Zrr,crap,gamma)
if(resample) deallocate(boot)

#ifdef MPI
call MPI_Barrier(MPI_COMM_WORLD,ierr)
call MPI_FINALIZE( ierr )
#endif

if(myid==master) then
  print*, "Exit... stage left!"
  return
else
  return
endif

end subroutine tpcf_iso

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine parseOptions ()
     USE extension, ONLY : getArgument, nArguments
      implicit none

      integer(kind=4)       :: i,n
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
            stop "Good day to you!"
         end if
         call getArgument(i+1,arg)
         select case (opt)
         case ('-gal')
            file1 = trim(arg)
         case ('-ran')
            file2 = trim(arg)
         case ('-out')
            outfile = trim(arg)
         case ('-rmin')
            read (arg,*) rmin
         case('-rmax')
            read (arg,*) rmax
         case('-nbins')
            read (arg,*) nbins
         case('-tbins')
            read (arg,*) ntbin
         case ('-wgt')
            wgt=.true.
!            wgt=trim(arg)
         case('-log')
            logbins=.true.
!            logbins=trim(arg)
         case ('-proj')
            proj=.true.
         case('-iso')
            ISOFLAG = trim(arg)
         case('-decp')        
            decomp = trim(arg)
         case ('-help')
            call print_help
      stop
         case default
            print '("unknown option ",a6," ignored")', opt
            stop
         end select
      end do

      call check_params()

   end subroutine parseOptions
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine print_help

print*,'CALLING SEQUENCE:'
print*,'      2pcf [-gal gal_file] [-ran ran_file] [-out out_file] [-wgt WGT]'
print*,'           [-iso ISO] [-decp DECOMP][-rmin Rmin] [-rmax Rmax] [-nbins Nbins]'
print*,'           [-log LOG] [-proj PROJ] '
print*,' '
print*,'      eg: 2pcf -gal dr9.gal -ran dr9.ran -rmin 1.0 -rmax 100.0 -nbins 20 -wgt .true.'
print*,' '
print*,'INPUTS:'
print*,'       gal_file/ranfile = strings containing the name of a 3/4 column point file'
print*,'                 lines 1-3: XYZ comoving, optional 4: weight  '
print*,'             OR  lines 1-2: RA,DEC(degs), optional 3: weight  '
print*,'   '
print*,'       out_file = filename for output correlation func. '
print*,'   '
print*,'       Rmin = Min. radial seperation (in units of input catalogue)'
print*,'              Rmin is fixed = 0 for anisotropic SIGMA-PI decomposition. May change in future release.'
print*,' '
print*,'       Rmax = Max. radial seperation (in units of input catalogue)'
print*,' '
print*,'OPTIONAL:'
print*,'       Nbins  = Number of radial bins to calculate. Default is 1.'
print*,'                In the anisotripic SIGMA-PI case, this will be the number of bins in both axes'
print*,' '
print*,'       NTbins = Number of theta/angualr bins to calculate, only if usnig -decmp SMU in anisotropic case'
print*,' '
print*,'       PROJ   = Logical, if true then do 2D angular correlations. Data must be RA,DEC in degs.'
print*,'                Possible values: .true.  OR  .false.'
print*,' '
print*,'       WGT    = logical value to define if weights are to be included'
print*,'                Possible values: .true.  OR  .false.'
print*,'                In this case, the input gal/ran files should be 4 columns'
print*,' '
print*,'       LOG    = Logical to define if bins should be logarithmically spaced'
print*,'                Possible values: .true.  OR  .false.'
print*,''
print*,'       ISO    = Flag for either istoropic (ISO) or anisotropic (ANISO) clustering. If ANISO, choose'
print*,'                decomposition method using -decp flag shown below.'
print*,''
print*,'       DECOMP = Decomposition of the anisotropic pair binning. SIGPI for the usual sigma-pi'
print*,'                or SMU for s-mu seperation, where is the pair distance and mu is the cosine'
print*,'                of the angle between the vector to gal1 and the connecting vector between '
print*,'                the galaxies pair.'
print*,''
   end subroutine print_help

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine check_params()

    if(rmin<0) stop "Rmin < 0, seems odd. I will stop now ;)"
    if(rmin>=rmax) stop "Rmin >= Rmax, seems odd. I will stop now ;)"
    if(nbins<=0) stop "Nbins < 1 , seems odd. I will stop now ;)"
    if(ntbin<=0) stop "Ntbins < 1, seems odd. I will stop now ;)"
    if(ISOFLAG=='ANISO'.and.(decomp.ne.'SIGPI'.and.decomp.ne.'SMU')) stop "Unknown Decomposition, seems odd. I will stop now ;)"

  end subroutine check_params

subroutine default_params()
  d=3
  wgt=.false.
  logbins=.false.
  proj=.false.
  resample=.false.
  resample_method='boot' !'jack' for jakknife or 'boot' for bootstrap
  outfile='result.txt'
  nbins=1
  ntbin=1
  rmin=0
  rmax=0
  Nresample=1
  ISOFLAG='ISO'

  dotpcf=.false.
  dosigpi=.false.
  doap=.false.

end subroutine default_params

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

!   Nrand=0
!   open(11,file=file2,status='unknown')
! 21 continue
!   read(11,*,err=21,end=22)aux
!   Nrand=Nrand+1
!   goto 21
! 22 close(11)
!   if(myid==0) print*,'Preparing to read ',Nrand, 'random points'

end subroutine count_files

   subroutine read_files ()
      implicit none
      double precision aux1, aux2

  if(myid==0)  print*, 'opening ',file1
  open(11,file=file1,status='unknown')

31  do i=1,Ndata

       read(11,*,err=31,end=14)aux1,aux2,gamma(i,1),gamma(i,2)
       call sph2cart(1.d0,aux1,aux2,my_array(1,i),my_array(2,i),my_array(3,i))

  enddo
14 close(11)
  if(myid==0)  print*,'Finished reading data file 1'  

!   if(myid==0)  print*, 'opening ',file2
!   open(11,file=file2,status='unknown')
! 32  do i=Ndata+1,Ndata+Nrand

!        read(11,*,err=32,end=23)aux1,aux2,gamma(i,1),gamma(i,2)
!        call sph2cart(1.d0,aux1,aux2,my_array(1,i),my_array(2,i),my_array(3,i))

!   enddo
! 23 close(11)
!   if(myid==0)  print*,'Finished reading data file 2'  

end subroutine read_files

  subroutine allocate_arrays_2pcf ()
  implicit none

!  allocate(gamma(Ndata,2))
  allocate(resultsb(Ndata+Nrand))
  allocate(results2(Ndata+Nrand))

  allocate(Zdd(nbins,Nresample))
  allocate(Zdr(nbins,Nresample))
  allocate(Zrr(nbins,Nresample))
  allocate(crap(nbins,Nresample))

  crap(:,:)=0d0
  Zdd(:,:)=0d0
  Zdr(:,:)=0d0
  Zrr(:,:)=0d0

  if(resample) allocate(boot(Ndata+Nrand))

  end subroutine allocate_arrays_2pcf


 subroutine sph2cart (radius,theta,phi,x,y,z)
      double precision radius,theta,phi,x,y,z
      
      z=radius*sin(phi)

      y=radius*cos(phi)*sin(theta)

      x=radius*cos(phi)*cos(theta)      
      return
    end subroutine sph2cart

    SUBROUTINE cart2sph(x,y,z,r,t,p)
     double precision t,p,r,x,y,z
     
     r=sqrt(x*x + y*y + z*z)
     t=atan2(y,x)
     if(t < 0.) then
        t=2.*3.14159+t
     endif
     p=asin(z/r)
     
     return
   end SUBROUTINE cart2sph

   subroutine bootstrapper()
    real ran

    do i=1,Ndata+Nrand
      call random_number(ran)
      boot(i)=floor(ran*Nresample)+1
    enddo

    return
   end subroutine bootstrapper

   subroutine mpi_collect()
    if(myid==master) crap(:,:)=0.d0
call MPI_Barrier(MPI_COMM_WORLD,ierr)
call MPI_REDUCE( Zdd, crap, nbins*Nresample, MPI_DOUBLE_PRECISION,MPI_SUM, master, MPI_COMM_WORLD, ierr )
if ( myid == master ) then !in master thread
Zdd=crap
endif
call MPI_Barrier(MPI_COMM_WORLD,ierr)

if(myid==master) crap(:,:)=0.d0
call MPI_Barrier(MPI_COMM_WORLD,ierr)
call MPI_REDUCE( Zdr, crap, nbins*Nresample, MPI_DOUBLE_PRECISION,MPI_SUM, master, MPI_COMM_WORLD, ierr )
if ( myid == master ) then !in master thread
Zdr=crap
endif
call MPI_Barrier(MPI_COMM_WORLD,ierr)

if(myid==master) crap(:,:)=0.d0
call MPI_Barrier(MPI_COMM_WORLD,ierr)
call MPI_REDUCE( Zrr, crap, nbins*Nresample, MPI_DOUBLE_PRECISION,MPI_SUM, master, MPI_COMM_WORLD, ierr )

if ( myid == master ) then !in master thread
Zrr=crap
endif
call MPI_Barrier(MPI_COMM_WORLD,ierr)
    end subroutine mpi_collect

end program  shear
