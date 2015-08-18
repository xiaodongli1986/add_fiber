program TwoPCF
  use kdtree2_module
  use kdtree2_precision_module
  implicit none
#ifdef MPI
    include 'mpif.h'
#endif

  real, dimension(:,:), allocatable :: my_array,vel

  type(kdtree2), pointer :: tree
  ! this is how you declare a tree in your main program
  double precision, allocatable :: Zdd(:,:),Zdr(:,:),Zrr(:,:),crap(:,:)
  real, allocatable :: wgt1(:),bins(:,:), mkd(:)
  real rmin,rmax,aux1,aux2,aux3
  double precision, allocatable :: dd(:,:,:),drn(:,:,:),rr(:,:,:),xi(:,:,:), v1(:),v2(:)
  integer, allocatable :: boot(:)
  double precision :: ndd,ndr,nrr,sig,pi,mind,dist1,dist2,nddn,ndrn,nrrn
  integer :: k, i, j, d,chunk,nbins,ntbin,Nresample
  character :: file1*600, file2*600, outfile*600
  type(kdtree2_result),allocatable :: resultsb(:)
  integer   ::  Ndata,Nrand,nn1,nn2,sigbin,pibin,tbin,bin,wgal,wran
  double precision :: aux,dr,odr,theta,var
  logical :: wgt,logbins,proj,resample,dosigpi,dotpcf,doap,saveran,loadran,readjk,marked
  CHARACTER(LEN=600) :: resample_method,ISOFLAG,DECOMP,ranfile*100
  integer :: myid , ntasks , ierr
  double precision, parameter :: pival=3.14159265358979323846
  real, allocatable :: av(:),sigma(:),delta(:,:),cov(:,:),Qmatrix(:,:)

    real :: den,radius,nex,v,totv,vd,dist
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

!--- analyse input files ---
  call count_files()

!--- allocate some initial memory ---
  Allocate(my_array(d,Ndata))
  allocate(wgt1(Ndata))
  allocate(mkd(Ndata))
  mkd(:)=1.0
  allocate(bins(nbins,2))
  if (resample) allocate(boot(Ndata+Nrand))

!--- load the input files ---
  call read_files()

!--- build the KD tree ---
if(myid==master) print*,' building tree on node', myid
tree => kdtree2_create(my_array,sort=.true.,rearrange=.true.)     ! this is how you create a tree.
if(myid==master) print*,' built tree on node', myid

!--- select the statistical analysis ---
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

!--- Do the statistical analysis ---

    call tpcf_iso()

if(myid==master) stop "Ciao for now!"

contains

!------------------------------------------------------------
!------------------------------------------------------------
!------------------------------------------------------------

subroutine tpcf_iso ()

if(myid==master) print*,'Beginning fiber collision simulation'

call allocate_arrays_2pcf()

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


do i=1,Ndata,1
wgt1(i)=1.0
enddo

print*,'starting gal loop'

#ifdef MPI
chunk=floor(float(Ndata)/ntasks)+1
do i=max(myid*chunk,1),min(((myid+1)*chunk)-1,Ndata),1
#else
do i=1,Ndata,1
#endif

   if(wgt1(i)>=1.0) then 
   
       call kdtree2_r_nearest_around_point(tp=tree,idxin=i,correltime=-1,r2=bins(nbins,2)*bins(nbins,2),&
&      nfound=nn2,nalloc=(Ndata),results=resultsb)
   
       do k=1,nn2
            if(resultsb(k)%idx==i) goto 33
            wgt1(i)=wgt1(i)+wgt1(resultsb(k)%idx)
            wgt1(resultsb(k)%idx)=0.0
33 continue            
        enddo  
   
   endif   
!   stop
enddo 
  
print*,'finished DD loop in thread:', myid

#ifdef MPI
call MPI_Barrier(MPI_COMM_WORLD,ierr)
#endif

if(myid==master) then

open(11,file=outfile,status='unknown')

do i=1,Ndata,1

   call cart2sph(my_array(1,i),my_array(2,i),my_array(3,i),aux1,aux2,aux3)
   write(11,'(3(e14.7,1x))') aux2,aux3,wgt1(i)

enddo
 close(11)
endif


! deallocate the memory for the data.
    call kdtree2_destroy(tree)  

    call deallocate_arrays()
    
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
      logical lexist
      integer(kind=4)       :: i,n
      character(len=6)   :: opt
      character(len=600) :: arg

      n = nArguments()

      if (n < 6 .and. myid==master) then
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
            if(trim(arg)=='.true.') then
              wgt=.true.
              if(myid==master) print*,'Using weighted points.'
            elseif(trim(arg)=='.false.') then
              wgt=.false.
            else
                stop "Error using -wgt flag. Either .true. or .false."
            endif
         case ('-mkd')
            if(trim(arg)=='.true.') then
              marked=.true.
              if(myid==master) print*,'Using marked statistics'
            elseif(trim(arg)=='.false.') then
              marked=.false.
            else
                stop "Error using -mkd flag. Either .true. or .false."
            endif
         case('-log')
            if(trim(arg)=='.true.') then
              logbins=.true.
              if(myid==master)  print*,'Using logarithmic binning scheme'
            elseif(trim(arg)=='.false.') then
              logbins=.false.
            else
                stop "Error using -log flag. Either .true. or .false."
            endif
         case ('-proj')
            proj=.true.
         case('-iso')
            ISOFLAG = trim(arg)
         case('-decp') 
            decomp = trim(arg)
         case('-err') 
             if (trim(arg)=='boot' .or. trim(arg)=='jack' .or. trim(arg)=='read') then
                continue
             else 
                print*,'-err flag requires either "boot" for bootstrap, "jack" for jackknife,'
                print*,' or "read" for reading error subsample from data file, expected last column.'
                stop
            endif
            resample_method=trim(arg)
            resample=.true.
            if(myid==master) print*," I will calculate the statistical covariance"
            if(myid==master) print*," I will ",trim(resample_method)," the data for covariance"
            if (trim(arg)=='read') readjk=.true.
         case('-nerr') 
            read (arg,*) Nresample
            if(myid==master) print*," Using ", Nresample," resamples for the statistical covariance."
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

      call check_params()

   end subroutine parseOptions
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine print_help

print*,'CALLING SEQUENCE:'
print*,'      add_fiber [-gal gal_file] [-out out_file] [-wgt WGT]'
print*,'                [-rmax Rmax]'
print*,' '
print*,'      eg: 2pcf -gal dr9.gal -ran dr9.ran -rmin 1.0 -rmax 100.0 -nbins 20 -wgt .true.'
print*,' '
print*,'INPUTS:'
print*,'       gal_file = string containing the name of a 2 column point file'
print*,'                 columns 1-2:  RA, DEC (in degrees)'
print*,'   '
print*,'       out_file = filename for output correlation func. '
print*,'   '
print*,'       Rmax = Fiber collision radius (in units of degs)'
print*,' '
print*,'OPTIONAL:'
print*,'       WGT    = logical value to define if weights are to be included'
print*,'                Possible values: .true.  OR  .false.'
print*,'                In this case, the input gal/ran files should be 4 columns'
print*,' '
   end subroutine print_help

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine check_params()

    if(rmin<0) stop "Rmin < 0, seems odd. I will stop now ;)"
    if(rmin>=rmax) stop "Rmin >= Rmax, seems odd. I will stop now ;)"
    if(nbins<=0) stop "Nbins < 1 , seems odd. I will stop now ;)"
    if(ntbin<=0) stop "Ntbins < 1, seems odd. I will stop now ;)"
    if(ISOFLAG=='ANISO'.and.(decomp.ne.'SIGPI'.and.decomp.ne.'SMU')) stop "Unknown Decomposition, seems odd. I will stop now ;)"
    if(resample_method=='jack') stop "Jackknife not yet implemented, you will have to JK sample yourself and use the read option."
    if(resample.and.(loadran.or.saveran)) stop "Option conflict: resampling and loading/saving random cannot work together."
    if(resample .and. Nresample.le.1) stop "Need more than one resample. Change Nerr on input. I will stop now."
    
  end subroutine check_params

subroutine default_params()
  d=3
  wgt=.false.
  logbins=.false.
  proj=.true.
  resample=.false.
  resample_method='boot' !'jack' for jakknife or 'boot' for bootstrap or 'read' for reading
  outfile='result.txt'
  nbins=1
  ntbin=1
  rmin=0
  rmax=0
  Nresample=1
  ISOFLAG='ISO'
  readjk=.false.
  marked=.false.

  dotpcf=.false.
  dosigpi=.false.
  doap=.false.

  loadran=.false.
  saveran=.false.

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

  Nrand=0


end subroutine count_files

   subroutine read_files ()
      implicit none
      real aux1, aux2

  if(myid==0)  print*, 'opening ',file1
  open(11,file=file1,status='unknown')

31  do i=1,Ndata

    read(11,*,err=31,end=14)aux1,aux2
    call sph2cart(1.0,aux1*3.14159/180.0,aux2*3.14159/180.0,my_array(1,i),my_array(2,i),my_array(3,i))

  enddo
14 close(11)
  if(myid==0)  print*,'Finished reading data file 1'  

end subroutine read_files

  subroutine allocate_arrays_2pcf ()
  implicit none
  allocate(resultsb(Ndata))

  allocate(Zdd(nbins,Nresample))
  allocate(Zdr(nbins,Nresample))
  allocate(Zrr(nbins,Nresample))
  allocate(crap(nbins,Nresample))

  crap(:,:)=0d0
  Zdd(:,:)=0d0
  Zdr(:,:)=0d0
  Zrr(:,:)=0d0


  end subroutine allocate_arrays_2pcf


subroutine deallocate_arrays()
implicit none
if (allocated(my_array)) deallocate(my_array)
if (allocated(Zdd)) deallocate(Zdd)
if (allocated(Zdr)) deallocate(Zdr)
if (allocated(Zrr)) deallocate(Zrr)
if (allocated(crap)) deallocate(crap)
if (allocated(av)) deallocate(av)
if (allocated(sigma)) deallocate(sigma)
if (allocated(delta)) deallocate(delta)
if (allocated(cov)) deallocate(cov)
if (allocated(wgt1)) deallocate(wgt1)
if (allocated(boot)) deallocate(boot)
if (allocated(bins)) deallocate(bins)
if (allocated(v1)) deallocate(v1)
if (allocated(v2)) deallocate(v2)

end subroutine deallocate_arrays

subroutine sph2cart (radius,theta,phi,x,y,z)
implicit none
real radius,theta,phi,x,y,z
      
z=radius*sin(phi)
y=radius*cos(phi)*sin(theta)
x=radius*cos(phi)*cos(theta)      
      
return
end subroutine sph2cart

SUBROUTINE cart2sph(x,y,z,r,t,p)
implicit none
real t,p,r,x,y,z
     
r=sqrt(x*x + y*y + z*z)
t=atan2(y,x)

if(t < 0.) then
    t=2.*3.14159+t
endif
p=asin(z/r)
return
end SUBROUTINE cart2sph

end program  TwoPCF
