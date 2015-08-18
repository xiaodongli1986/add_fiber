module utils
implicit none

contains

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

    subroutine getArgument(index, argument)
      use global_vars
      ! ===========================================================
      integer(kind=4), intent(in) :: index
      character(len=*), intent(out) :: argument
      integer(kind=4) :: i1
      ! ===========================================================
      i1 = index + arg_shift
      call getarg(i1, argument)

      return
    end subroutine getArgument

    ! ===========================================================
    function nArguments() result(narg)
      use global_vars
      ! ===========================================================
      integer(kind=4) :: narg
      ! ===========================================================

      narg = iargc() - arg_shift
!VF      narg = nargs() - arg_shift

      return
    end function nArguments