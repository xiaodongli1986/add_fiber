program density
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
  real rmin,rmax,aux1,aux2,aux3,dh,pos(3),minx,maxx,miny,maxy,minz,maxz
  double precision, allocatable :: dd(:,:,:),drn(:,:,:),rr(:,:,:),xi(:,:,:), v1(:),v2(:)
  integer, allocatable :: boot(:)
  double precision :: ndd,ndr,nrr,sig,pi,mind,dist1,dist2,nddn,ndrn,nrrn
  integer :: k, i, j, d,chunk,nbins,ntbin,Nresample,nx,ny,nz
  character :: file1*100, file2*100, outfile*100
  type(kdtree2_result),allocatable :: resultsb(:)
  integer   ::  Ndata,Nrand,nn1,nn2,sigbin,pibin,tbin,bin,wgal,wran
  double precision :: aux,dr,odr,theta,var
  logical :: wgt,logbins,proj,resample,dosigpi,dotpcf,doap,saveran,loadran,readjk,marked
  CHARACTER(LEN=100) :: resample_method,ISOFLAG,DECOMP,ranfile*100
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
  Allocate(my_array(d,Ndata+Nrand))


!--- load the input files ---
  call read_files()

!--- build the KD tree ---
if(myid==master) print*,' building tree on node', myid
tree => kdtree2_create(my_array,sort=.true.,rearrange=.true.)     ! this is how you create a tree.
if(myid==master) print*,' built tree on node', myid

!--- select the statistical analysis ---
!--- Do the statistical analysis ---

    call tpcf_iso()

if(myid==master) stop "Ciao for now!"

contains

!------------------------------------------------------------
!------------------------------------------------------------
!------------------------------------------------------------

subroutine tpcf_iso ()

allocate(resultsb(Ndata+Nrand))

dh=15.0
minx=minval(my_array(1,:))
maxx=maxval(my_array(1,:))
miny=minval(my_array(2,:))
maxy=maxval(my_array(2,:))
minz=minval(my_array(3,:))
maxz=maxval(my_array(3,:))

nx=int((maxx-minx)/dh)
ny=int((maxy-miny)/dh)
nz=int((maxz-minz)/dh)

print*,'xmin/xmax:',minx,maxx
print*,'ymin/ymax:',miny,maxy
print*,'zmin/zmax:',minz,maxz

print*,'nx,ny,nz:',nx,ny,nz

open(11,file=outfile,status='unknown')


do i=1,nx
    do j=1,ny
        do k=1,nz
        
        pos(1)=minx+ (i-1)*dh + dh/2.
        pos(2)=miny+ (j-1)*dh + dh/2.
        pos(3)=minz+ (k-1)*dh + dh/2.

        call kdtree2_n_nearest(tp=tree,qv=pos(:),nn=21,results=resultsb)

        call spline()

        write(11,'(4(e14.7,1x))') pos(1),pos(2),pos(3),den               
  !      print*,i,j,k,pos(1),pos(2),pos(3),den 

        enddo
    enddo
enddo

 close(11)


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
      character(len=100) :: arg

      n = nArguments()

      if (n < 2 .and. myid==master) then
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

         case ('-out')
            outfile = trim(arg)

         case ('-help')
            call print_help
      stop
         case default
            print '("unknown option ",a6," ignored")', opt
            stop
         end select
      end do

!      call check_params()

   end subroutine parseOptions
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine print_help

print*,'CALLING SEQUENCE:'
print*,'      2pcf [-gal gal_file] [-ran ran_file] [-out out_file] [-wgt WGT]'
print*,'           [-iso ISO] [-decp DECOMP][-rmin Rmin] [-rmax Rmax] [-nbins Nbins]'
print*,'           [-tbins NTbins] [-log LOG] [-proj PROJ] [-RR RR_Counts] [-err Err]'
print*,'           [-nerr NERR][-mkd MKD]'
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
print*,'    RR_Counts = string containing the name of a file regarding RR data counts.  '
print*,'                If file exists, then the RR counting will be skipped (saveing CPU time)'
print*,'                If it does not exist, the file will be created for future use.'
print*,'                Consider using this option if you are doing Monte Carlo runs with same random catalogue.' 
print*,''
print*,'        ERR  =  The error treatment, "boot" will create bootstrap samples and use them to construct '
print*,'                a covariance estimation. "jack" will create jackknife samples (NOT IMPLEMENTED YET)'
print*,'                "read" will read the error subsample from the data, expected last column (integer!)'
print*,''
print*,'       NERR  =  Number of jackknife or bootstrap samples to create. Covariance and correlation '
print*,'                coefficients contained within the file covariance.out'
print*,''
print*,'       MKD    = Logical to define if marked statistics should be used.'
print*,'                Possible values: .true.  OR  .false.'
print*,''
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

    read(11,*,err=31,end=14)my_array(1,i),my_array(2,i),my_array(3,i)
!    call sph2cart(1.0,aux1*3.14159/180.0,aux2*3.14159/180.0,my_array(1,i),my_array(2,i),my_array(3,i))

  enddo
14 close(11)
  if(myid==0)  print*,'Finished reading data file 1'  

end subroutine read_files


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

subroutine spline()
implicit none
integer ns,jj
real h,h3,q,q2,q3

        ns=20
        den=0.
        h=resultsb(ns+1)%dis *0.5
        h3=h*h*h
        do jj=1,ns
          q=resultsb(jj)%dis/h
          q2=q*q
          q3=q*q*q
          if (q.le.1) then
            den = den+1.-1.5*q2+0.75*q3
          else if (q.le.2) then
            den = den+0.25*(2-q)*(2-q)*(2-q)
          endif
          !print*,j,den(i)
        enddo
        den=den/h3/pival
end subroutine spline

end program  density
