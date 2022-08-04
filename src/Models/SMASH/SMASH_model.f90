module SMASH_model

!~**********************************************************************
!~* Purpose: Interface to SMASH model
!~**********************************************************************
!~* Programmer: Ben Renard, INRAE Aix
!~**********************************************************************
!~* Last modified: 02/08/2022
!~**********************************************************************
!~* Comments:
!~**********************************************************************
!~* References:
!~**********************************************************************
!~* 2Do List:
!~**********************************************************************
!~* Quick description of public procedures:
!~*     1. SMASH_GetParNumber, number of inferred parameters
!~*     2. SMASH_Apply, apply model
!~*     3. SMASH_XtraRead, read Xtra model information
!~**********************************************************************

use kinds_dmsl_kit ! numeric kind definitions from DMSL
use Geodesy_tools,only:rasterGridType,interpolateOnGrid
use DataRW_tools,only:WriteRasterASC,ReadRasterASC,DatRead

implicit none
Private
public :: SMASH_GetParNumber, SMASH_Run, SMASH_Load, SMASH_XtraRead

! Module variables
type(rasterGridType)::SMASHgrid
real(mrk), allocatable::SMASHanchors(:,:)
character(250),allocatable::SMASHparNames(:)

Contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine SMASH_GetParNumber(npar,err,mess)
!^**********************************************************************
!^* Purpose: number of inferred parameters
!^*          = number of anchor points + 1 smoothing parameter per SMASH par
!^**********************************************************************
!^* Programmer: Ben Renard, INRAE Aix
!^**********************************************************************
!^* Last modified: 02/08/2022
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*     1.anchors, definition of n anchor points. Array with n rows and
!^*                3 columns: x, y, index of SMASH par the row refers to
!^* OUT
!^*     1. npar, par. number
!^*     2.err, error code; <0:Warning, ==0:OK, >0: Error
!^*     3.mess, error message
!^**********************************************************************
integer(mik), intent(out)::npar,err
character(*),intent(out)::mess
! locals
character(250),parameter::procname='SMASH_GetParNumber'

err=0;mess=''
npar=size(SMASHanchors,dim=1)+nint(maxval(SMASHanchors(:,3)),mik)

end subroutine SMASH_GetParNumber

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine SMASH_Run(runScript,projectDir,thetaDir,QSIMfile,theta,Y,feas,err,mess)
!^**********************************************************************
!^* Purpose: Run SMASH
!^**********************************************************************
!^* Programmer: Ben Renard, INRAE Aix
!^**********************************************************************
!^* Last modified:
!^* - 02/08/2022
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*     1. XXX
!^* OUT
!^*     1. XXX
!^*     3. feas, feasible?
!^*     4.err, error code; <0:Warning, ==0:OK, >0: Error
!^*     5.mess, error message
!^**********************************************************************
use utilities_dmsl_kit,only:getSpareUnit
character(*), intent(in)::runScript,projectDir,thetaDir,QSIMfile
real(mrk), intent(in)::theta(:)
real(mrk), intent(out)::Y(:,:)
logical, intent(out)::feas
integer(mik), intent(out)::err
character(*),intent(out)::mess
!locals
character(250),parameter::procname='SMASH_Run'
character(len_uLongStr)::cmdString,QSIMfullPath
integer(mik)::k,npar,i,nmask,runStatus,unt
real(mrk)::smooth,foo
logical::mask(size(SMASHanchors,dim=1))
real(mrk),allocatable::anc(:,:), M(:,:)

! Init
err=0;mess='';feas=.true.;Y=undefRN

k=1
! Make sense of parameter vector and write spatialized parameters
npar=nint(maxval(SMASHanchors(:,3)),mik)
do i=1,npar
    ! Smoothing parameter used in IDW interpolation
    smooth=theta(k)
    ! Perform interpolation from anchors to grid
    mask=(SMASHanchors(:,3)==i);nMask=count(mask)
    if(allocated(anc)) deallocate(anc); allocate(anc(nMask,3))
    anc(:,1)=pack(SMASHanchors(:,1),mask)
    anc(:,2)=pack(SMASHanchors(:,2),mask)
    anc(:,3)=theta( (k+1) : (k+nMask) )
    k=k+nMask+1
    call interpolateOnGrid(pts=anc,grid=SMASHgrid,method='IDW',par=(/smooth/),M=M,err=err,mess=mess)
    if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
    ! Write spatialized parameter to file in thetaDir
    call WriteRasterASC(grid=SMASHgrid,M=M,file=trim(projectDir)//trim(thetaDir)//trim(SMASHparNames(i))//'.asc',err=err,mess=mess)
    if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
enddo

! Run SMASH through Python script
cmdString='python3'//' '//trim(runScript)//' '//trim(projectDir)//' '//trim(projectDir)//trim(thetaDir)
call execute_command_line (trim(cmdString),exitStat=runStatus,cmdMsg=mess)
if(runStatus/=0) then;feas=.false.;return;endif

! Read simulated streamflow
QSIMfullPath=trim(projectDir)//trim(QSIMfile)
call getSpareUnit(unt,err,mess)
if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
open(unit=unt,file=QSIMfullPath,status='old', iostat=err)
if(err>0) then;mess=trim(procname)//':problem opening QSIM file: '//trim(QSIMfullPath);return;endif
! headers
do i=1,2
    read(unt,*,iostat=err)
    if(err>0) then;mess=trim(procname)//':problem reading headers of QSIM file: '//trim(QSIMfullPath);return;endif
enddo
do i=1,size(Y,dim=1)
    read(unt,*,iostat=err) foo, Y(i,:)
    if(err>0) then;mess=trim(procname)//':problem reading QSIM file: '//trim(QSIMfullPath);return;endif
enddo
end subroutine SMASH_Run

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine SMASH_Load(loadScript,projectDir,precipDir,petDir,err,mess)
!^**********************************************************************
!^* Purpose: Load SMASH
!^**********************************************************************
!^* Programmer: Ben Renard, INRAE Aix
!^**********************************************************************
!^* Last modified:
!^* - 02/08/2022
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*     1. XXX
!^* OUT
!^*     1. XXX
!^*     3. feas, feasible?
!^*     4.err, error code; <0:Warning, ==0:OK, >0: Error
!^*     5.mess, error message
!^**********************************************************************
character(*), intent(in)::loadScript,projectDir,precipDir,petDir
integer(mik), intent(out)::err
character(*),intent(out)::mess
!locals
character(250),parameter::procname='SMASH_Load'
character(len_uLongStr)::cmdString

! Init
err=0;mess=''

! Run Python script that loads SMASH and saves the related model object in projectDir
cmdString='python3'//' '//trim(loadScript)//' '//trim(projectDir)//' '//trim(precipDir)//' '//trim(petDir)
call execute_command_line (trim(cmdString),exitStat=err,cmdMsg=mess)
if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif

end subroutine SMASH_Load

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine SMASH_XtraRead(file,xtra,err,mess)
!^**********************************************************************
!^* Purpose: Read Xtra information
!^**********************************************************************
!^* Programmer: Ben Renard, INRAE Aix
!^**********************************************************************
!^* Last modified: 02/08/2022
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*     1. file, Xtra file
!^* OUT
!^*     1. xtra, xtra information
!^*     2. err, error code; <0:Warning, ==0:OK, >0: Error
!^*     3. mess, error message
!^**********************************************************************
use types_dmsl_kit, only:data_ricz_type
use utilities_dmsl_kit,only:getSpareUnit
character(*), intent(in)::file
type(data_ricz_type),intent(out):: xtra
integer(mik), intent(out)::err
character(*),intent(out)::mess
! locals
character(250),parameter::procname='SMASH_XtraRead'
integer(mik)::unt,k
character(250):: head(3)
real(mrk),pointer::foo(:,:)

err=0;mess=''

! String information: python scripts and directories
if(associated(xtra%cp1)) nullify(xtra%cp1); allocate(xtra%cp1(7))
call getSpareUnit(unt,err,mess)
if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
open(unit=unt,file=trim(file), status='old', iostat=err)
if(err>0) then;mess=trim(procname)//':problem opening file '//trim(file);return;endif
k=0
k=k+1
read(unt,*,iostat=err) xtra%cp1(k) ! loadScript
if(err/=0) then;mess=trim(procname)//':problem reading file '//trim(file);return;endif
k=k+1
read(unt,*,iostat=err) xtra%cp1(k) ! runScript
if(err/=0) then;mess=trim(procname)//':problem reading file '//trim(file);return;endif
k=k+1
read(unt,*,iostat=err) xtra%cp1(k) ! projectDir
if(err/=0) then;mess=trim(procname)//':problem reading file '//trim(file);return;endif
k=k+1
read(unt,*,iostat=err) xtra%cp1(k) ! precipDir
if(err/=0) then;mess=trim(procname)//':problem reading file '//trim(file);return;endif
k=k+1
read(unt,*,iostat=err) xtra%cp1(k) ! petDir
if(err/=0) then;mess=trim(procname)//':problem reading file '//trim(file);return;endif
k=k+1
read(unt,*,iostat=err) xtra%cp1(k) ! BaMstuffDir
if(err/=0) then;mess=trim(procname)//':problem reading file '//trim(file);return;endif
k=k+1
read(unt,*,iostat=err) xtra%cp1(k) ! QSIM file
if(err/=0) then;mess=trim(procname)//':problem reading file '//trim(file);return;endif
close(unt)
! SMASH parameter names (hard-coded for the moment)
if(allocated(SMASHparNames)) deallocate(SMASHparNames); allocate(SMASHparNames(4))
SMASHparNames=(/'cp ','ctr','cr ','ml '/)

! Read grid in BaMstuff directory
call ReadRasterASC(file=trim(xtra%cp1(3))//trim(xtra%cp1(6))//'grid.asc',gridOnly=.true.,grid=SMASHgrid,err=err,mess=mess)
if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif

! Read anchors in BaMstuff directory
call DatRead(file=trim(xtra%cp1(3))//trim(xtra%cp1(6))//'anchors.txt',ncol=3,y=foo,headers=head,err=err,mess=mess)
if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
if(allocated(SMASHanchors)) deallocate(SMASHanchors);allocate(SMASHanchors(size(foo,dim=1),size(foo,dim=2)))
SMASHanchors=foo

end subroutine SMASH_XtraRead
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!==============!
! Private subs !
!==============!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module SMASH_model