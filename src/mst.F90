program TwoPCF
  use kdtree2_module
  use kdtree2_precision_module
  implicit none
#ifdef MPI
    include 'mpif.h'
#endif

  real, dimension(:,:), allocatable :: my_array

  type(kdtree2), pointer :: tree
  ! this is how you declare a tree in your main program
  double precision, allocatable :: Zdd(:,:),Zdr(:,:),Zrr(:,:),crap(:,:)
  real, allocatable :: wgt1(:),bins(:,:)
  real rmin,rmax
  double precision, allocatable :: dd(:,:,:),drn(:,:,:),rr(:,:,:),xi(:,:,:), v1(:),v2(:)
  integer, allocatable :: boot(:)
  double precision :: ndd,ndr,nrr,sig,pi,mind,dist1,dist2,nddn,ndrn,nrrn
  integer :: k, i, j, d,chunk,nbins,ntbin,Nresample
  character :: file1*100, file2*100, outfile*100
  type(kdtree2_result),allocatable :: resultsb(:)
  integer   ::  Ndata,Nrand,nn1,nn2,sigbin,pibin,tbin,bin
  double precision :: aux,dr,odr,theta
  logical :: wgt,logbins,proj,resample,dosigpi,dotpcf,doap,saveran,loadran,readjk
  CHARACTER(LEN=100) :: resample_method,ISOFLAG,DECOMP,ranfile*100
  integer :: myid , ntasks , ierr
  double precision, parameter :: pival=3.14159265358979323846
  real, allocatable :: av(:),sigma(:),delta(:,:),cov(:,:),Qmatrix(:,:)

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
  if ( wgt )  allocate(wgt1(Ndata+Nrand))
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
if(dotpcf) call tpcf_iso()

if(dosigpi) call sigma_pi()

if(doap) call ap_test()

if(myid==master) stop "Ciao for now!"

contains

!------------------------------------------------------------
subroutine mst ()
if(myid==master) print*,'Calculating the minimum spanning tree'


do i=1,Ndata,1
    chunk=1
    tree=tree+1
    wgt1(i)=1      

    continue
    
    do j=1,chunk
    
    call kdtree2_n_nearest_around_point(tp=tree,idxin=i,correltime=-1,nn=chunk,results=resultsb)

    if(wgt1(resultsb(1)%idx) .ne. 0 .and. wgt1(resultsb(1)%idx) .ne. tree) then
    print*, 'bugger!'
    endif
    
    if(resultsb(1)%dis<=rmin) then
        wgt1(resultsb(1)%idx)=tree
    endif


   call kdtree2_r_nearest_around_point(tp=tree,idxin=i,correltime=-1,r2=2*rmax*rmax,&
&       nfound=nn2,nalloc=(Nrand+Ndata),results=resultsb)


end subroutine mst

subroutine sigma_pi ()
if(myid==master) print*,'Calculating the anisotropic correlation function decomposed into sigma-pi'

  ntbin=nbins
  rmin=0.0         ! will fix in the near future!!!!!!!!

  call allocate_arrays_sig_pi()

if(resample) then 
    if(readjk) then
        continue   
    else
        call bootstrapper()
    endif
endif

  dr=(rmax-rmin)/nbins
  odr=1./dr

  do i=1,nbins
    bins(i,1)=rmin + dr*(float(i)-1.0)
    bins(i,2)=rmin + dr*i
    if(myid==0) print*,bins(i,:)
  enddo

print*,'beginning DD loop on thread:',myid

#ifdef MPI
chunk=floor(float(Ndata)/ntasks)+1
do i=max(myid*chunk,1),min(((myid+1)*chunk)-1,Ndata),1
#else
do i=1,Ndata,1
#endif

   call kdtree2_r_nearest_around_point(tp=tree,idxin=i,correltime=-1,r2=2*rmax*rmax,&
&       nfound=nn2,nalloc=(Nrand+Ndata),results=resultsb)

   nn1=kdtree2_r_count_around_point(tp=tree,idxin=i,correltime=-1,r2=2*rmin*rmin)

    v1=my_array(1:3,i)
    dist1=sqrt(v1(1)**2+v1(2)**2+v1(3)**2)

     do j=nn1+1,nn2

        v2=my_array(1:3,resultsb(j)%idx)
        theta=v1(1)*v2(1) + v1(2)*v2(2) + v1(3)*v2(3)
        theta=(mod(abs(theta/((sqrt(v1(1)**2 + v1(2)**2 + v1(3)**2)) * (sqrt(v2(1)**2 +v2(2)**2 +v2(3)**2)))),1.0))

        dist2=sqrt(v2(1)**2 +v2(2)**2 +v2(3)**2)
        mind=min(dist1,dist2)
        sig=abs(dist1-dist2)
        pi=sqrt(2.d0*(mind*mind - mind*mind*theta))

        sigbin=floor(sig*odr)+1
        pibin=floor(pi*odr)+1
        if(sigbin .gt. nbins .or. sigbin.le.0 .or. pibin.gt.nbins .or. pibin.le.0 .or. sig==0 .or. pi==0) goto 41
        bin=nbins*(pibin-1) + sigbin

         if(resultsb(j)%idx <=Ndata) then

            if ( wgt ) then

               if(resample) then
                  Zdd(bin,boot(resultsb(j)%idx))=Zdd(bin,boot(resultsb(j)%idx))+(wgt1(resultsb(j)%idx)*wgt1(i))
                  Zdd(bin,boot(i))=Zdd(bin,boot(i))+(wgt1(resultsb(j)%idx)*wgt1(i))
               else
                  Zdd(bin,1)=Zdd(bin,1)+(wgt1(resultsb(j)%idx)*wgt1(i))
               endif

            else
               if(resample) then
                  Zdd(bin,boot(resultsb(j)%idx))=Zdd(bin,boot(resultsb(j)%idx))+1.d0
                  Zdd(bin,boot(i))=Zdd(bin,boot(i))+1.d0
               else
                  Zdd(bin,1)=Zdd(bin,1)+1.d0
                endif
            endif

         else

            if ( wgt ) then
                if(resample) then
                  Zdr(bin,boot(resultsb(j)%idx))=Zdr(bin,boot(resultsb(j)%idx))+(wgt1(resultsb(j)%idx)*wgt1(i))
                  Zdr(bin,boot(i))=Zdr(bin,boot(i))+(wgt1(resultsb(j)%idx)*wgt1(i))
               else
                  Zdr(bin,1)=Zdr(bin,1)+(wgt1(resultsb(j)%idx)*wgt1(i))
                endif
            else
              if(resample) then
                  Zdr(bin,boot(resultsb(j)%idx))=Zdr(bin,boot(resultsb(j)%idx))+1.d0
                  Zdr(bin,boot(i))=Zdr(bin,boot(i))+1.d0
               else
                 Zdr(bin,1)=Zdr(bin,1)+1.d0
               endif
            endif
            
         endif
41       continue
   enddo
enddo
print*,'finished DD loop in thread:', myid

if(loadran) then
  if(myid==master) print*,'Random counts exist, skipping some calculation'
  goto 57
endif

!RR loop
print*,'beginning RR loop on thread:',myid
#ifdef MPI
chunk=floor(float(Nrand)/ntasks)+1
do i=Ndata+max(myid*chunk,1),Ndata+min(((myid+1)*chunk)-1,Nrand),1
#else
do i=Ndata+1,Ndata+Nrand,1
#endif

   call kdtree2_r_nearest_around_point(tp=tree,idxin=i,correltime=-1,r2=2*rmax*rmax,&
& nfound=nn2,nalloc=(Nrand+Ndata),results=resultsb)

   nn1=kdtree2_r_count_around_point(tp=tree,idxin=i,correltime=-1,r2=2*rmin*rmin)

    v1=my_array(1:3,i) 
    dist1=sqrt(v1(1)*v1(1)+v1(2)*v1(2)+v1(3)*v1(3))   

      do j=nn1+1,nn2

        v2=my_array(1:3,resultsb(j)%idx)
        theta=v1(1)*v2(1) + v1(2)*v2(2) + v1(3)*v2(3)
        theta=(mod(abs(theta/((sqrt(v1(1)**2 + v1(2)**2 + v1(3)**2)) * (sqrt(v2(1)**2 +v2(2)**2 +v2(3)**2)))),1.0))

        dist2=sqrt(v2(1)*v2(1)+v2(2)*v2(2)+v2(3)*v2(3))

        mind=min(dist1,dist2)
        sig=abs(dist1-dist2)
        pi=sqrt(2.d0*(mind*mind - mind*mind*theta))

        sigbin=floor(sig*odr)+1
        pibin=floor(pi*odr)+1

        bin=nbins*(pibin-1) + sigbin

        if(sigbin .gt. nbins .or. sigbin.le.0 .or. pibin.gt.nbins .or. pibin.le.0 .or. sig==0 .or. pi==0) goto 42

        if(resultsb(j)%idx > Ndata) then
            if ( wgt ) then
               if(resample) then
                  Zrr(bin,boot(resultsb(j)%idx))=Zrr(bin,boot(resultsb(j)%idx))+(wgt1(resultsb(j)%idx)*wgt1(i))
                  Zrr(bin,boot(i))=Zrr(bin,boot(i))+(wgt1(resultsb(j)%idx)*wgt1(i))
               else
                  Zrr(bin,1)=Zrr(bin,1)+(wgt1(resultsb(j)%idx)*wgt1(i))
               endif
            else
              if(resample) then
                  Zrr(bin,boot(resultsb(j)%idx))=Zrr(bin,boot(resultsb(j)%idx))+1.d0
                  Zrr(bin,boot(i))=Zrr(bin,boot(i))+1.d0
               else
                 Zrr(bin,1)=Zrr(bin,1)+1.d0
               endif
            endif

         endif

42        continue

   enddo
enddo
print*,'finished RR loop on thread:', myid

#ifdef MPI
call MPI_Barrier(MPI_COMM_WORLD,ierr)
#endif

57 continue

#ifdef MPI
if(myid==master) print*,'MPI: collecting results'
call mpi_collect_aniso()
#endif

#ifdef MPI
call MPI_Barrier(MPI_COMM_WORLD,ierr)
call MPI_FINALIZE( ierr )
#endif

if(myid==master) then
    print*,'Analysing resamples'
    if(resample) then
     call analyse_samples()
    endif
endif

if(myid==master) then
    print*,'Normalising counts'
call normalise_counts()

open(11,file=outfile,status='unknown')
if(resample) then
    write(11,'(A)')'Sig_min, Sig_max, Pi_min, Pi_max, XI_LS, error'
else
    write(11,'(A)')'Sig_min, Sig_max, Pi_min, Pi_max, DD, DR, RR, XI_natural, XI_Davis, XI_LS (most accurate)'
endif

do i=1,nbins
  do j=1,ntbin

     bin=nbins*(i-1) + j

if(resample) then

   write(11,'(6(e14.7,1x))') (i-1)*dr,i*dr,(j-1)*dr,j*dr,(Zdd(bin,1)-2.0*Zdr(bin,1)+Zrr(bin,1))/Zrr(bin,1),sigma(bin)
   print*,(i-1)*dr,i*dr,(j-1)*dr,j*dr,av(bin),sigma(bin)

else

   write(11,'(10(e14.7,1x))')(i-1)*dr,i*dr,(j-1)*dr,j*dr,Zdd(bin,1),Zdr(bin,1),Zrr(bin,1), &
& Zdd(bin,1)/Zrr(bin,1) - 1,Zdd(bin,1)/Zdr(bin,1) - 1, (Zdd(bin,1)-2.0*Zdr(bin,1)+Zrr(bin,1))/Zrr(bin,1)
   print*,(i-1)*dr,i*dr,(j-1)*dr,j*dr,Zdd(bin,1),Zdr(bin,1),Zrr(bin,1), &
& Zdd(bin,1)/Zrr(bin,1) - 1,Zdd(bin,1)/Zdr(bin,1) - 1, (Zdd(bin,1)-2.0*Zdr(bin,1)+Zrr(bin,1))/Zrr(bin,1)
endif

  enddo
enddo
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

end subroutine sigma_pi

!------------------------------------------------------------

subroutine ap_test ()
if(myid==master) print*,'Calculating the anisotropic correlation function decomposed into dist - mu coordinates'

allocate(v1(d))
allocate(v2(d))

call allocate_arrays_ap_test()
  
if(resample) then 
    if(readjk) then
        continue   
    else
        call bootstrapper()
    endif
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

    v1=my_array(1:3,i)
   
   do k=nbins,1,-1
      nn1=kdtree2_r_count_around_point(tp=tree,idxin=i,correltime=-1,r2=bins(k,1)*bins(k,1))
      do j=nn1+1,nn2,1

        v2=my_array(1:3,resultsb(j)%idx)         ! position of 2nd gal
        v2=[v2(1)-v1(1),v2(2)-v1(2),v2(3)-v1(3)] ! connecting vector
        theta=v1(1)*v2(1) + v1(2)*v2(2) + v1(3)*v2(3)
        theta=(mod(abs(theta/((sqrt(v1(1)**2 + v1(2)**2 + v1(3)**2)) * (sqrt(v2(1)**2 +v2(2)**2 +v2(3)**2)))),1.0))
        tbin=floor(theta*ntbin)+1

        bin=(k-1)*ntbin + tbin 

         if(resultsb(j)%idx <=Ndata) then
            if ( wgt ) then
            
               if(resample) then
                  Zdd(bin,boot(resultsb(j)%idx))=Zdd(bin,boot(resultsb(j)%idx))+(wgt1(resultsb(j)%idx)*wgt1(i))
                  Zdd(bin,boot(i))=Zdd(bin,boot(i))+(wgt1(resultsb(j)%idx)*wgt1(i))
               else
                  Zdd(bin,1)=Zdd(bin,1)+(wgt1(resultsb(j)%idx)*wgt1(i))
               endif
            
            else
               if(resample) then
                  Zdd(bin,boot(resultsb(j)%idx))=Zdd(bin,boot(resultsb(j)%idx))+1.d0
                  Zdd(bin,boot(i))=Zdd(bin,boot(i))+1.d0
               else
                  Zdd(bin,1)=Zdd(bin,1)+1.d0
               endif
            
            endif
            
         else

            if ( wgt ) then
                if(resample) then
                  Zdr(bin,boot(resultsb(j)%idx))=Zdr(bin,boot(resultsb(j)%idx))+(wgt1(resultsb(j)%idx)*wgt1(i))
                  Zdr(bin,boot(i))=Zdr(bin,boot(i))+(wgt1(resultsb(j)%idx)*wgt1(i))
               else
                  Zdr(bin,1)=Zdr(bin,1)+(wgt1(resultsb(j)%idx)*wgt1(i))
                endif
            else
              if(resample) then
                  Zdr(bin,boot(resultsb(j)%idx))=Zdr(bin,boot(resultsb(j)%idx))+1.d0
                  Zdr(bin,boot(i))=Zdr(bin,boot(i))+1.d0
               else
                 Zdr(bin,1)=Zdr(bin,1)+1.d0
               endif

            endif
            
         endif
      enddo
      nn2=nn1

   enddo
enddo
print*,'finished DD loop in thread:', myid


if(loadran) then
  if(myid==master) print*,'Random counts exist, skipping some calculation'
  goto 56
endif

!RR loop
print*,'beginning RR loop on thread:',myid

#ifdef MPI
chunk=floor(float(Nrand)/ntasks)+1
do i=Ndata+max(myid*chunk,1),Ndata+min(((myid+1)*chunk)-1,Nrand),1
#else
do i=Ndata+1,Ndata+Nrand,1
#endif

   call kdtree2_r_nearest_around_point(tp=tree,idxin=i,correltime=-1,r2=bins(nbins,2)*bins(nbins,2),&
& nfound=nn2,nalloc=(Nrand+Ndata),results=resultsb)

    v1=my_array(1:3,i) 

   do k=nbins,1,-1
      nn1=kdtree2_r_count_around_point(tp=tree,idxin=i,correltime=-1,r2=bins(k,1)*bins(k,1))
      do j=nn1+1,nn2,1

        v2=my_array(1:3,resultsb(j)%idx)         ! position of 2nd gal
        v2=[v2(1)-v1(1),v2(2)-v1(2),v2(3)-v1(3)] ! connecting vector
        theta=v1(1)*v2(1) + v1(2)*v2(2) + v1(3)*v2(3)
        theta=(mod(abs(theta/((sqrt(v1(1)**2 + v1(2)**2 + v1(3)**2)) * (sqrt(v2(1)**2 +v2(2)**2 +v2(3)**2)))),1.0))
        tbin=floor(theta*ntbin)+1

        bin=(k-1)*ntbin + tbin 

         if(resultsb(j)%idx >Ndata) then
            if ( wgt ) then
               if(resample) then
                  Zrr(bin,boot(resultsb(j)%idx))=Zrr(bin,boot(resultsb(j)%idx))+(wgt1(resultsb(j)%idx)*wgt1(i))
                  Zrr(bin,boot(i))=Zrr(bin,boot(i))+(wgt1(resultsb(j)%idx)*wgt1(i))
               else
                  Zrr(bin,1)=Zrr(bin,1)+(wgt1(resultsb(j)%idx)*wgt1(i))
               endif
            else
              if(resample) then
                  Zrr(bin,boot(resultsb(j)%idx))=Zrr(bin,boot(resultsb(j)%idx))+1.d0
                  Zrr(bin,boot(i))=Zrr(bin,boot(i))+1.d0
               else
                 Zrr(bin,1)=Zrr(bin,1)+1.d0
               endif

            endif

         endif
      enddo
      nn2=nn1

   enddo
enddo
print*,'finished RR loop on thread:', myid

56 continue

#ifdef MPI
call mpi_collect_aniso()
#endif

#ifdef MPI
call MPI_Barrier(MPI_COMM_WORLD,ierr)
call MPI_FINALIZE( ierr )
#endif

if(myid==master) then
    print*,'Analysing resamples'
    if(resample) then
     call analyse_samples()
    endif
endif

if(myid==master) then
  call normalise_counts()

open(11,file=outfile,status='unknown')
if (resample) then
write(11,*)'#  cos(theta)_min, cos(theta)_max, Rmin, Rmax,XI_LS, error'
else
write(11,*)'#  cos(theta)_min,        cos(theta)_max,          Rmin,                   Rmax,                &
& DD (normalised),        DR (normalised),        RR (normalised),         XI_natural,            XI_Davis, &
&              XI_LS (most accurate)'
endif
print*,'#  cos(theta)_min,        cos(theta)_max,          Rmin,                   Rmax,                  &
& DD (normalised),        DR (normalised),        RR (normalised),         XI_natural,            XI_Davis, &
&              XI_LS (most accurate)'
do i=1,nbins
  do j=ntbin,1,-1

        bin=(i-1)*ntbin + j

if (resample) then
   write(11,'(7(e14.7,1x))')1.d0-(float(j)/ntbin),1.d0-(float(j-1)/ntbin),bins(i,1),bins(i,2), &
&   (Zdd(bin,1)-2.0*Zdr(bin,1)+Zrr(bin,1))/Zrr(bin,1),sigma(bin)
else
   write(11,'(10(e14.7,1x))')1.d0-(float(j)/ntbin),1.d0-(float(j-1)/ntbin),bins(i,1),bins(i,2),Zdd(bin,1),Zdr(bin,1),Zrr(bin,1), &
& Zdd(bin,1)/Zrr(bin,1) - 1,Zdd(bin,1)/Zdr(bin,1) - 1, (Zdd(bin,1)-2.0*Zdr(bin,1)+Zrr(bin,1))/Zrr(bin,1)
endif
   print*,1.d0-(float(j)/ntbin),1.d0-(float(j-1)/ntbin),bins(i,1),bins(i,2),Zdd(bin,1),Zdr(bin,1),Zrr(bin,1), &
& Zdd(bin,1)/Zrr(bin,1) - 1,Zdd(bin,1)/Zdr(bin,1) - 1, (Zdd(bin,1)-2.0*Zdr(bin,1)+Zrr(bin,1))/Zrr(bin,1)

  enddo
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

end subroutine ap_test

!------------------------------------------------------------

subroutine tpcf_iso ()

if(proj) then
  if(myid==master) print*,'Calculating the isotropic 2D correlation function'
else
  if(myid==master) print*,'Calculating the isotropic 3D correlation function'
endif

call allocate_arrays_2pcf()

if(resample) then
    if(readjk) then
        if(maxval(boot).ne.Nresample) stop "why is boot samples not equal to number of resamples?!"
        continue   
    else
        call bootstrapper()
    endif
endif

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
                  Zdd(k,1)=Zdd(k,1)+(wgt1(resultsb(j)%idx)*wgt1(i))
               endif
            else
               if(resample) then
                  Zdd(k,boot(resultsb(j)%idx))=Zdd(k,boot(resultsb(j)%idx))+1.d0
                  Zdd(k,boot(i))=Zdd(k,boot(i))+1.d0
                else
                  Zdd(k,1)=Zdd(k,1)+1.d0
               endif
            endif
         else
            if ( wgt ) then
               if(resample) then
                  Zdr(k,boot(resultsb(j)%idx))=Zdr(k,boot(resultsb(j)%idx))+(wgt1(resultsb(j)%idx)*wgt1(i))
                  Zdr(k,boot(i))=Zdr(k,boot(i))+(wgt1(resultsb(j)%idx)*wgt1(i))
                else
                  Zdr(k,1)=Zdr(k,1)+(wgt1(resultsb(j)%idx)*wgt1(i))
               endif
            else
               if(resample) then
                  Zdr(k,boot(resultsb(j)%idx))=Zdr(k,boot(resultsb(j)%idx))+1.d0
                  Zdr(k,boot(i))=Zdr(k,boot(i))+1.d0
                else
                  Zdr(k,1)=Zdr(k,1)+1.d0
               endif
            endif
         endif
      enddo
      nn2=nn1

   enddo
enddo
print*,'finished DD loop in thread:', myid

if(loadran) then
  if(myid==master) print*,'Random counts exist, skipping some calculation'
  goto 55
endif

#ifdef MPI
call MPI_Barrier(MPI_COMM_WORLD,ierr)
#endif

!RR loop
print*,'beginning RR loop on thread:',myid
#ifdef MPI
chunk=floor(float(Nrand)/ntasks)+1
do i=Ndata+max(myid*chunk,1),Ndata+min(((myid+1)*chunk)-1,Nrand),1
#else
do i=Ndata+1,Ndata+Nrand,1
#endif
   call kdtree2_r_nearest_around_point(tp=tree,idxin=i,correltime=-1,r2=bins(nbins,2)*bins(nbins,2),&
&       nfound=nn2,nalloc=(Nrand+Ndata),results=resultsb)

   do k=nbins,1,-1

      nn1=kdtree2_r_count_around_point(tp=tree,idxin=i,correltime=-1,r2=bins(k,1)*bins(k,1))

      do j=nn1+1,nn2,1
         if(resultsb(j)%idx >Ndata) then
            if ( wgt ) then
               if(resample) then
                  Zrr(k,boot(resultsb(j)%idx))=Zrr(k,boot(resultsb(j)%idx))+(wgt1(resultsb(j)%idx)*wgt1(i))
                  Zrr(k,boot(i))=Zrr(k,boot(i))+(wgt1(resultsb(j)%idx)*wgt1(i))
                else
                  Zrr(k,1)=Zrr(k,1)+(wgt1(resultsb(j)%idx)*wgt1(i))
               endif
            else
               if(resample) then
                  Zrr(k,boot(resultsb(j)%idx))=Zrr(k,boot(resultsb(j)%idx))+1.d0
                  Zrr(k,boot(i))=Zrr(k,boot(i))+1.d0
               else
                  Zrr(k,1)=Zrr(k,1)+1.d0
            endif
          endif

       endif
      enddo
      nn2=nn1

   enddo

enddo
print*,'finished RR loop on thread:', myid

#ifdef MPI
call MPI_Barrier(MPI_COMM_WORLD,ierr)
#endif

55 continue

#ifdef MPI
call mpi_collect()
#endif

#ifdef MPI
call MPI_Barrier(MPI_COMM_WORLD,ierr)
call MPI_FINALIZE( ierr )
#endif

if(myid==master) then
    if(resample) then
     call analyse_samples()
    endif
endif


if(myid==master) then

call normalise_counts()
   
open(11,file=outfile,status='unknown')
if(resample) then
write(11,'(A)')'R_min, R_max, XI_LS (most accurate), Xi_error'
else
write(11,'(A)')'R_min, R_max, DD, DR, RR, XI_natural, XI_Davis, XI_LS (most accurate)'
endif

do i=1,nbins
if(resample) then

   write(11,'(5(e14.7,1x))') bins(i,1),bins(i,2),av(i),sigma(i)
   print*,bins(i,1),bins(i,2),(Zdd(i,1)-2.0*Zdr(i,1)+Zrr(i,1))/Zrr(i,1),sigma(i)

else

   write(11,'(8(e14.7,1x))') bins(i,1),bins(i,2),Zdd(i,1),Zdr(i,1),Zrr(i,1), &
& Zdd(i,1)/Zrr(i,1) - 1,Zdd(i,1)/Zdr(i,1) - 1, (Zdd(i,1)-2.0*Zdr(i,1)+Zrr(i,1))/Zrr(i,1)
   print*,bins(i,1),bins(i,2),Zdd(i,1),Zdr(i,1),Zrr(i,1), &
& Zdd(i,1)/Zrr(i,1) - 1,Zdd(i,1)/Zdr(i,1) - 1, (Zdd(i,1)-2.0*Zdr(i,1)+Zrr(i,1))/Zrr(i,1)
endif
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
            if(trim(arg)=='.true.') then
              wgt=.true.
              if(myid==master) print*,'Using weighted points.'
            elseif(trim(arg)=='.false.') then
              wgt=.false.
            else
                stop "Error using -wgt flag. Either .true. or .false."
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
print*,'      2pcf [-gal gal_file] [-ran ran_file] [-out out_file] [-wgt WGT]'
print*,'           [-iso ISO] [-decp DECOMP][-rmin Rmin] [-rmax Rmax] [-nbins Nbins]'
print*,'           [-tbins NTbins] [-log LOG] [-proj PROJ] [-RR RR_Counts] [-err Err]'
print*,'           [-nerr NERR]'
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
    
  end subroutine check_params

subroutine default_params()
  d=3
  wgt=.false.
  logbins=.false.
  proj=.false.
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
      real aux1, aux2

  if(myid==0)  print*, 'opening ',file1
  open(11,file=file1,status='unknown')

31  do i=1,Ndata

    if(proj) then
      if ( wgt ) then
       read(11,*,err=31,end=14)aux1,aux2,wgt1(i)
      else
       read(11,*,err=31,end=14)aux1,aux2
      endif
      call sph2cart(1.0,aux1,aux2,my_array(1,i),my_array(2,i),my_array(3,i))
    else
      if ( wgt ) then
        if (readjk) then
            read(11,*,err=31,end=14)my_array(1:3,i),wgt1(i),boot(i)
        else
            read(11,*,err=31,end=14)my_array(1:3,i),wgt1(i)
        endif
      else
        if (readjk) then
           read(11,*,err=31,end=14)my_array(1:3,i),boot(i)
        else
           read(11,*,err=31,end=14)my_array(1:3,i)
      endif
    endif 
    endif
  enddo
14 close(11)
  if(myid==0)  print*,'Finished reading data file 1'  

  if(myid==0)  print*, 'opening ',file2
  open(11,file=file2,status='unknown')
32  do i=Ndata+1,Ndata+Nrand

    if(proj) then
      if ( wgt ) then
       read(11,*,err=32,end=23)aux1,aux2,wgt1(i)
      else
       read(11,*,err=32,end=23)aux1,aux2
      endif
      call sph2cart(1.0,aux1,aux2,my_array(1,i),my_array(2,i),my_array(3,i))
    else
      if ( wgt ) then
        if (readjk) then
           read(11,*,err=32,end=23)my_array(1:3,i),wgt1(i),boot(i)
        else
           read(11,*,err=32,end=23)my_array(1:3,i),wgt1(i)
        endif
      else
        if (readjk) then
           read(11,*,err=32,end=23)my_array(1:3,i),boot(i)
        else
           read(11,*,err=32,end=23)my_array(1:3,i)
        endif
     endif
   endif

  enddo
23 close(11)
  if(myid==0)  print*,'Finished reading data file 2'  

end subroutine read_files

  subroutine allocate_arrays_2pcf ()
  implicit none
  allocate(resultsb(Ndata+Nrand))

  allocate(Zdd(nbins,Nresample))
  allocate(Zdr(nbins,Nresample))
  allocate(Zrr(nbins,Nresample))
  allocate(crap(nbins,Nresample))

  crap(:,:)=0d0
  Zdd(:,:)=0d0
  Zdr(:,:)=0d0
  Zrr(:,:)=0d0

  if(resample) then
    allocate(av(Nbins))
    allocate(sigma(Nbins))
    allocate(delta(Nresample,Nbins))
    allocate(cov(Nbins,Nbins))
    allocate(Qmatrix(Nresample,nbins))
  endif 

  end subroutine allocate_arrays_2pcf

  subroutine allocate_arrays_sig_pi ()
  implicit none
  allocate(v1(d))
  allocate(v2(d))

  allocate(resultsb(Ndata+Nrand))

  allocate(Zdd(nbins*ntbin,Nresample))
  allocate(Zdr(nbins*ntbin,Nresample))
  allocate(Zrr(nbins*ntbin,Nresample))
  allocate(crap(nbins*ntbin,Nresample))

  crap(:,:)=0.d0
  Zdd(:,:)=0.d0
  Zdr(:,:)=0.d0
  Zrr(:,:)=0.d0
  
  if(resample) then
    allocate(av(Nbins*ntbin))
    allocate(sigma(Nbins*ntbin))
    allocate(delta(Nresample,Nbins*ntbin))
    allocate(cov(Nbins*ntbin,Nbins*ntbin))
    allocate(Qmatrix(Nresample,nbins*ntbin))
  endif 
  
end subroutine allocate_arrays_sig_pi

subroutine allocate_arrays_ap_test ()
  implicit none
  allocate(resultsb(Ndata+Nrand))

  allocate(Zdd(nbins*ntbin,Nresample))
  allocate(Zdr(nbins*ntbin,Nresample))
  allocate(Zrr(nbins*ntbin,Nresample))
  allocate(crap(nbins*ntbin,Nresample))

  crap(:,:)=0d0
  Zdd(:,:)=0d0
  Zdr(:,:)=0d0
  Zrr(:,:)=0d0

  if(resample) then
    allocate(av(Nbins*ntbin))
    allocate(sigma(Nbins*ntbin))
    allocate(delta(Nresample,Nbins*ntbin))
    allocate(cov(Nbins*ntbin,Nbins*ntbin))
    allocate(Qmatrix(Nresample,nbins*ntbin))
  endif 
  
end subroutine allocate_arrays_ap_test

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
implicit none
real ran

    do i=1,Ndata+Nrand
      call random_number(ran)
      boot(i)=floor(ran*Nresample)+1
    enddo

return
end subroutine bootstrapper

subroutine mpi_collect()
implicit none

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

subroutine mpi_collect_aniso()
implicit none

    if(myid==master) crap(:,:)=0.d0
call MPI_Barrier(MPI_COMM_WORLD,ierr)
call MPI_REDUCE( Zdd, crap, nbins*ntbin*Nresample, MPI_DOUBLE_PRECISION,MPI_SUM, master, MPI_COMM_WORLD, ierr )
if ( myid == master ) then !in master thread
Zdd=crap
endif
call MPI_Barrier(MPI_COMM_WORLD,ierr)

if(myid==master) crap(:,:)=0.d0
call MPI_Barrier(MPI_COMM_WORLD,ierr)
call MPI_REDUCE( Zdr, crap, nbins*ntbin*Nresample, MPI_DOUBLE_PRECISION,MPI_SUM, master, MPI_COMM_WORLD, ierr )
if ( myid == master ) then !in master thread
Zdr=crap
endif
call MPI_Barrier(MPI_COMM_WORLD,ierr)

if(myid==master) crap(:,:)=0.d0
call MPI_Barrier(MPI_COMM_WORLD,ierr)
call MPI_REDUCE( Zrr, crap, nbins*ntbin*Nresample, MPI_DOUBLE_PRECISION,MPI_SUM, master, MPI_COMM_WORLD, ierr )

if ( myid == master ) then !in master thread
Zrr=crap
endif
call MPI_Barrier(MPI_COMM_WORLD,ierr)
    end subroutine mpi_collect_aniso

subroutine normalise_counts()
implicit none

if ( wgt ) then
   ndd=sum(wgt1(1:Ndata))**2
   nrr=sum(wgt1(Ndata+1:Ndata+Nrand))**2
   ndr=sqrt(ndd)*sqrt(nrr)

   Zdd=Zdd/ndd
   Zrr=Zrr/nrr
   Zdr=Zdr/ndr
else
  Zdd=Zdd/(float(Ndata)*(Ndata-1))
  Zrr=Zrr/(float(Nrand)*(Nrand-1))
  Zdr=Zdr/(float(Ndata)*(Nrand-1))
endif

if(saveran) then
  open(11,file=ranfile,status='unknown')
do i=1,nbins
  do j=ntbin,1,-1
        bin=(i-1)*ntbin + j
    write(11,*)Zrr(bin,1)
    enddo 
  enddo
  close(11)
endif

if(loadran) then
  open(11,file=ranfile,status='unknown')
  do i=1,nbins
    do j=ntbin,1,-1
        bin=(i-1)*ntbin + j
    read(11,*)Zrr(bin,1)
    enddo
  enddo
  close(11)
endif


end subroutine normalise_counts

subroutine covariance()
implicit none
integer count

delta=0.0

count=nbins*ntbin

print*,'count= ',count

do i=1,Nresample
    do j=1,count
      Qmatrix(i,j)=xi(i,j,3)
    enddo
enddo

print*,' calculate mean'
av=0.0
do i=1,count
   av(i)=sum(Qmatrix(1:Nresample,i))/float(Nresample)
enddo

print*,' calculate variance'
sigma=0.0   
do i=1,count
   do j=1,Nresample
      sigma(i)=sigma(i)+(Qmatrix(j,i)-av(i))**2.0
   enddo
   sigma(i)=sqrt(sigma(i)*(float(Nresample)-1.)/float(Nresample))
enddo

print*,' calculate delta'
do j=1,Nresample
   do i=1,count
      delta(j,i)=(Qmatrix(j,i)-av(i))/(sigma(i))
   enddo
enddo

cov=0.0
print*,' writing covariance'
open(18,file='covariance.out',status='unknown')
do i=1,count
   do j=1,count
      do k=1,Nresample
         cov(i,j)=cov(i,j)+((delta(k,i)*delta(k,j))*(float(Nresample)-1.)/float(Nresample))
      enddo
!      if(ISOFLAG=='ANISO') then
!       if (DECOMP=='SMU') then
!          write(18,'(4F8.2,2E11.3)')1.d0-(float(j)/ntbin),1.d0-(float(j-1)/ntbin),bins(i,1),&
!          &bins(i,2),sigma(i)*sigma(j)*cov(i,j),cov(i,j)
!        else
!          write(18,'(4F8.2,2E11.3)')(i-1)*dr,i*dr,(j-1)*dr,j*dr,sigma(i)*sigma(j)*cov(i,j),cov(i,j)
!        endif
!      else
!          write(18,'(4F8.2,2E11.3)')bins(i,1:2),bins(j,1:2),sigma(i)*sigma(j)*cov(i,j),cov(i,j)
!      endif
   enddo
   write(18,*)(sigma(i)*sigma(j)*cov(i,j), j=1,count)
!   write(12,*)theta(i,1:2),sig(i)
!   write(13,*)theta(i,1:2),av(i)
enddo
close(18)

end subroutine covariance

subroutine analyse_samples()
implicit none

    print*,'allocating arrays'
    allocate(dd(Nresample,nbins*ntbin,2))
    allocate(drn(Nresample,nbins*ntbin,2))
    allocate(rr(Nresample,nbins*ntbin,2))
    allocate(xi(Nresample,nbins*ntbin,3))
    dd=0.d0
    drn=0.d0
    rr=0.d0

    if(.not.wgt) then
        allocate(wgt1(Ndata+Nrand))
        wgt1=1.0
    endif

    print*,'gathering resamples'
    do i=1,Nresample

        ndd=0.d0
        ndr=0.d0
        nrr=0.d0

        do j=1,Ndata
            if(boot(j).ne.i) then
                 ndd=ndd+wgt1(j)
            endif
        enddo 
        do j=Ndata+1,Ndata+Nrand
           if(boot(j).ne.i) then
                nrr=nrr+wgt1(j)
           endif
        enddo 
    
        nddn=(sum(wgt1(1:Ndata))-0)*sum(wgt1(1:Ndata)-0)
        nrrn=(sum(wgt1(Ndata+1:Ndata+Nrand))-0)*(sum(wgt1(Ndata+1:Ndata+Nrand))-0)
        ndrn=(sum(wgt1(1:Ndata))-0)*(sum(wgt1(Ndata+1:Ndata+Nrand))-0)

        do j=1,Nresample
            if(j.ne.i) then
              do k=1,nbins*ntbin
                dd(i,k,1)=dd(i,k,1)+Zdd(k,j)
                drn(i,k,1)=drn(i,k,1)+Zdr(k,j)
                rr(i,k,1)=rr(i,k,1)+Zrr(k,j)
              enddo
            endif
        enddo

        do j=1,nbins*ntbin
            dd(i,j,2)=dd(i,j,1)/nddn
            drn(i,j,2)=drn(i,j,1)/ndrn
            rr(i,j,2)=rr(i,j,1)/nrrn

            xi(i,j,1)=dd(i,j,2)/rr(i,j,2) -1.d0
            xi(i,j,2)=dd(i,j,2)/drn(i,j,2) -1.d0
            xi(i,j,3)=(dd(i,j,2)-2.*drn(i,j,2)+rr(i,j,2))/rr(i,j,2)
!            print*,i,j,dd(i,j,2),drn(i,j,2),rr(i,j,2),xi(i,j,1:3)
        enddo
    enddo

    print*,'doing covariance'
       call covariance()
       print*,'done with covariance'
       do j=2,Nresample
        do i=1,nbins*ntbin
            Zdd(i,1)=Zdd(i,1)+Zdd(i,j)
            Zdr(i,1)=Zdr(i,1)+Zdr(i,j)
            Zrr(i,1)=Zrr(i,1)+Zrr(i,j)
        enddo
       enddo   

        print*,'writing resample files'
         open(17,file='jk.txt',status='unknown')
         do j=1,nbins*ntbin
            write(17,*) (xi(i,j,3), i=1,nresample)
         enddo
         close(17)
    
end subroutine analyse_samples


end program  TwoPCF
