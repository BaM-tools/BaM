module MAGE_model

!~**********************************************************************
!~* Purpose: Interface to MAGE model
!~**********************************************************************
!~* Programmer: Ben Renard and Felipe Mendez, INRAE Aix & INRAE Lyon
!~**********************************************************************
!~* Last modified: 29/08/2022
!~**********************************************************************
!~* Comments:
!~**********************************************************************
!~* References:
!~**********************************************************************
!~* 2Do List:
!~**********************************************************************
!~* Quick description of public procedures:
!~*     1. MAGE_GetParNumber, number of inferred parameters
!~*     2. MAGE_Run, run MAGE
!~*     3. MAGE_XtraRead, read Xtra model information
!~*     4. MAGE_Setup, complete model setup after reading Xtra information
!~**********************************************************************

use kinds_dmsl_kit ! numeric kind definitions from DMSL

implicit none
Private
public :: MAGE_GetParNumber, MAGE_Run, MAGE_XtraRead, MAGE_SetUp

! Module variables
character(len_uLongStr)::RUGfile='' ! .RUG file specifying rugosities
real(mrk), allocatable::RUGx(:) ! x's at which rugosities are specified
integer(mik), allocatable::RUGb(:) ! "Bief" (reach?) index
real(mrk), allocatable::ZxKmin(:,:) ! covariates Z's of the regression Kmin(x)=a1Z1(x)+...+apZp(x). Dimension size(RUGb)*p.
real(mrk), allocatable::ZxKmoy(:,:) ! covariates Z's of the regression Kmoy(x)=a1Z1(x)+...+apZp(x). Dimension size(RUGb)*p.
logical::applyExpKmin,applyExpKmoy ! apply exponential transformation to computed K's to ensure positivity?
character(len_uLongStr),allocatable::outFiles(:) ! list of files where MAGE outputs (Qsim) are read

Contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine MAGE_GetParNumber(npar,err,mess)
!^**********************************************************************
!^* Purpose: number of inferred parameters
!^*          = number of "troncons" * 2 (main chanel and floodway)
!^**********************************************************************
!^* Programmer: Ben Renard and Felipe Mendez, INRAE Aix & INRAE Lyon
!^**********************************************************************
!^* Last modified:29/08/2022
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*     nothing
!^* OUT
!^*     1. npar, par. number
!^*     2.err, error code; <0:Warning, ==0:OK, >0: Error
!^*     3.mess, error message
!^**********************************************************************
integer(mik), intent(out)::npar,err
character(*),intent(out)::mess
! locals
character(250),parameter::procname='MAGE_GetParNumber'

err=0;mess=''
npar=size(ZxKmin,dim=2)+size(ZxKmoy,dim=2)

end subroutine MAGE_GetParNumber

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine MAGE_Run(exeFile,projectDir,REPfile,theta,Y,feas,err,mess)
!^**********************************************************************
!^* Purpose: Run MAGE
!^**********************************************************************
!^* Programmer: Ben Renard and Felipe Mendez, INRAE Aix & INRAE Lyon
!^**********************************************************************
!^* Last modified: 29/08/2022
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*     1. exeFile, MAGE executable
!^*     2. projectDir, MAGE project directory
!^*     3. REPfile, MAGE .REP file
!^*     4. theta, parameter vector
!^* OUT
!^*     1. Y, discharge computed at requested points
!^*     3. feas, feasible?
!^*     4. err, error code; <0:Warning, ==0:OK, >0: Error
!^*     5. mess, error message
!^**********************************************************************
use utilities_dmsl_kit,only:getSpareUnit
use DataRW_tools,only:DatRead
character(*), intent(in)::exeFile,projectDir,REPfile
real(mrk), intent(in)::theta(:)
real(mrk), intent(out)::Y(:,:)
logical, intent(out)::feas
integer(mik), intent(out)::err
character(*),intent(out)::mess
!locals
character(250),parameter::procname='MAGE_Run'
integer(mik),parameter::outFileNCOL=2
integer(mik)::i,n,unt,k
character(250)::forma,head(outFileNCOL)
character(len_uLongStr)::cmdString
real(mrk), pointer::foo(:,:)
real(mrk)::Kmin(size(RUGb)),Kmoy(size(RUGb))

! Init
err=0;mess='';feas=.true.;Y=undefRN
n=size(RUGb)
forma='(A1,I3,6x,2F10.0,2F10.2)'

! Size check
k=size(theta)
if(k /= (size(ZxKmin,dim=2)+size(ZxKmoy,dim=2))) then
    mess=trim(procname)//': length of theta should be equal to the total number of columns in Zfiles'
    err=1;return
endif

! Compute Kmin and Kmoy from regression
Kmin=matmul(ZxKmin, theta( 1:size(ZxKmin,dim=2) ) )
Kmoy=matmul(ZxKmoy, theta( (size(ZxKmin,dim=2)+1):k ) )
if(applyExpKmin) Kmin=exp(Kmin)
if(applyExpKmoy) Kmoy=exp(Kmoy)
! Write Kmin and Kmoy in .RUG file
call getSpareUnit(unt,err,mess)
if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
open(unit=unt,file=trim(projectDir)//trim(RUGfile),status='replace',iostat=err)
if(err/=0) then;mess=trim(procname)//':problem creating RUG file';return;endif
write(unt,'(A)') '**************** STRICKLERS **********************'
write(unt,'(A)') '*Bief        x deb     x fin     K min     K moy  '
write(unt,'(A)') '*---======----------==========----------=========='
do i=1,n
    write(unt,forma) 'K',RUGb(i),RUGx(i),RUGx(i+1),Kmin(i),Kmoy(i)
enddo
close(unt)

! Call mage
cmdString='cd '//trim(projectDir)//'&&'//trim(exeFile)//' '//trim(REPfile)
call execute_command_line (trim(cmdString),wait=.true.,exitStat=err,cmdMsg=mess)
if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif

! Read outputs into Y
do i=1,size(outFiles)
    call DatRead(file=outFiles(i),ncol=outFileNCOL,y=foo,headers=head,err=err,mess=mess)
    if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
    if( size(foo,dim=1) /= size(Y,dim=1) ) then ! MAGE didn't run properly
        feas=.false.;return
    endif
    Y(:,i)=foo(:,outFileNCOL)
enddo

end subroutine MAGE_Run

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine MAGE_XtraRead(file,xtra,err,mess)
!^**********************************************************************
!^* Purpose: Read Xtra information
!^**********************************************************************
!^* Programmer: Ben Renard and Felipe Mendez, INRAE Aix & INRAE Lyon
!^**********************************************************************
!^* Last modified: 29/08/2022
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
character(250),parameter::procname='MAGE_XtraRead'
integer(mik)::unt

err=0;mess=''

call getSpareUnit(unt,err,mess)
if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
open(unit=unt,file=trim(file), status='old', iostat=err)
if(err>0) then;mess=trim(procname)//':problem opening file '//trim(file);return;endif
read(unt,*,iostat=err) xtra%cs1 ! 1. exe file
if(err/=0) then;mess=trim(procname)//':problem reading file '//trim(file);return;endif
read(unt,*,iostat=err) xtra%cs2 ! 2. Project Directory
if(err/=0) then;mess=trim(procname)//':problem reading file '//trim(file);return;endif
read(unt,*,iostat=err) xtra%cs3 ! 3. REP file
if(err/=0) then;mess=trim(procname)//':problem reading file '//trim(file);return;endif
! K-Z regressions
if(associated(xtra%cp1)) nullify(xtra%cp1);allocate(xtra%cp1(2))
read(unt,*,iostat=err) xtra%cp1(1) ! 4. Z file (covariates for Kmin-Z regression)
if(err/=0) then;mess=trim(procname)//':problem reading file '//trim(file);return;endif
read(unt,*,iostat=err) xtra%is1    ! 5. Number of columns in Z file
if(err/=0) then;mess=trim(procname)//':problem reading file '//trim(file);return;endif
read(unt,*,iostat=err) xtra%ls1    ! 6. apply exponential transformation to computed  Kmin's?
if(err/=0) then;mess=trim(procname)//':problem reading file '//trim(file);return;endif
read(unt,*,iostat=err) xtra%cp1(2) ! 7. Z file (covariates for Kmoy-Z regression)
if(err/=0) then;mess=trim(procname)//':problem reading file '//trim(file);return;endif
read(unt,*,iostat=err) xtra%is2    ! 8. Number of columns in Z file
if(err/=0) then;mess=trim(procname)//':problem reading file '//trim(file);return;endif
read(unt,*,iostat=err) xtra%ls2    ! 9. apply exponential transformation to computed  Kmoy's?
if(err/=0) then;mess=trim(procname)//':problem reading file '//trim(file);return;endif
close(unt)

end subroutine MAGE_XtraRead

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine MAGE_setUp(projectDir,REPfile,Zfile,ZnCol,doExp,err,mess)
!^**********************************************************************
!^* Purpose: Finalize MAGE setup after having read Xtra information. In
!^*          particular this sub defines module variables (RUG stuff and outFiles)
!^**********************************************************************
!^* Programmer: Ben Renard and Felipe Mendez, INRAE Aix & INRAE Lyon
!^**********************************************************************
!^* Last modified: 24/02/2023
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*     1. projectDir, MAGE project directory
!^*     2. REPfile, MAGE .REP file
!^*     3. Zfile, containing ovariates for K-Z regression
!^*     4. doExp, apply exponential transformation to computed  K's?
!^* OUT
!^*     1.err, error code; <0:Warning, ==0:OK, >0: Error
!^*     2.mess, error message
!^**********************************************************************
use utilities_dmsl_kit,only:getNumItemsInFile,getSpareUnit
character(*), intent(in)::projectDir,REPfile,Zfile(2)
integer(mik), intent(in)::ZnCol(2)
logical, intent(in)::doExp(2)
integer(mik), intent(out)::err
character(*),intent(out)::mess
!locals
character(250),parameter::procname='MAGE_setUp'
integer(mik)::unt,i,n,k,j
character(250)::foo,prefix,var,ib,pk,mode,dt0,forma
character(len_uLongStr)::line,fname
character(len_uLongStr),allocatable::fnames(:)
real(mrk)::c1,c2
real(mrk),allocatable::Zx(:,:)

! Init
err=0;mess=''

! Open REP file
call getSpareUnit(unt,err,mess)
if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
open(unit=unt,file=trim(projectDir)//trim(REPfile), status='old', iostat=err)
if(err>0) then;mess=trim(procname)//':problem opening file '//trim(projectDir)//trim(REPfile);return;endif
call getNumItemsInFile(unt=unt,preRwnd=.true.,nskip=0,nitems=n,postPos=0,jchar=foo,err=err,message=mess)
if(err>0) then;mess=trim(procname)//trim(mess);return;endif

! Interpret REP file
k=0
if(allocated(fnames)) deallocate(fnames);allocate(fnames(n))
do i=1,n
    read(unt,'(A)',iostat=err) line
    if(err>0) then;mess=trim(procname)//':problem reading file '//trim(REPfile);return;endif
    if(line(1:3)=='RUG') RUGfile=trim(line(5:))
    if(line(1:3)=='CSV') then
        read(line,*,iostat=err) foo,prefix,var,ib,pk,mode,dt0
        if(err>0) then;mess=trim(procname)//':problem reading file '//trim(REPfile);return;endif
        k=k+1
        fname=trim(projectDir)//trim(prefix)//'_'//trim(var)//'_'//trim(ib)//'_'//trim(pk)//'.csv'
        fnames(k)=fname
    endif
enddo
if(k==0) then;err=1;mess=trim(procname)//'FATAL: no CSV file found.';return;endif
if(allocated(outFiles)) deallocate(outFiles);allocate(outFiles(k))
outFiles=fnames(1:k)
close(unt)

! Read RUG file
if(RUGfile=='') then;err=1;mess=trim(procname)//'FATAL: no .RUG file found';return;endif
call getSpareUnit(unt,err,mess)
if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
open(unit=unt,file=trim(projectDir)//trim(RUGfile), status='old', iostat=err)
if(err>0) then;mess=trim(procname)//':problem opening file '//trim(projectDir)//trim(RUGfile);return;endif
call getNumItemsInFile(unt=unt,preRwnd=.true.,nskip=3,nitems=n,postPos=0,jchar=foo,err=err,message=mess)
if(err>0) then;mess=trim(procname)//trim(mess);return;endif
if(allocated(RUGx)) deallocate(RUGx)
if(allocated(RUGb)) deallocate(RUGb)
allocate(RUGx(n+1),RUGb(n))
rewind(unt)
do i=1,3;read(unt,*);enddo
forma='(1x,I3,6x,4F10.0)'
do i=1,n
    read(unt,forma,iostat=err) RUGb(i),RUGx(i),RUGx(i+1),c1,c2
    if(err/=0) then;mess=trim(procname)//':problem reading file '//trim(projectDir)//trim(RUGfile);return;endif
enddo
close(unt)

! Read Z file
n=size(RUGb)
do j=1,2
    if(Zfile(j)=='') then
        if(allocated(Zx)) deallocate(Zx);allocate(Zx(n,n))
        Zx(:,:)=0._mrk
        forall(i=1:n)Zx(i,i)=1._mrk
    else
        call getSpareUnit(unt,err,mess)
        if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
        open(unit=unt,file=trim(Zfile(j)), status='old', iostat=err)
        if(err>0) then;mess=trim(procname)//':problem opening file '//trim(Zfile(j));return;endif
        call getNumItemsInFile(unt=unt,preRwnd=.true.,nskip=0,nitems=k,postPos=0,jchar=foo,err=err,message=mess)
        if(err>0) then;mess=trim(procname)//trim(mess);return;endif
        if(k /= (n+1)) then
            mess=trim(procname)//':Zfile should have n+1 rows (n rows of the RUG table + one header line). Zfile='//trim(Zfile(j))
            err=1;return
        endif
        if(allocated(Zx)) deallocate(Zx);allocate(Zx(n,ZnCol(j)))
        rewind(unt)
        read(unt,*) ! header
        do i=1,n
            read(unt,*,iostat=err) Zx(i,:)
            if(err/=0) then;mess=trim(procname)//':problem reading file '//trim(Zfile(j));return;endif
        enddo
        close(unt)
    endif
    if(j==1) then ! Applies to Kmin
        if(allocated(ZxKmin)) deallocate(ZxKmin);allocate(ZxKmin(size(Zx,dim=1),size(Zx,dim=2)))
        ZxKmin=Zx
    else ! Applies to Kmoy
        if(allocated(ZxKmoy)) deallocate(ZxKmoy);allocate(ZxKmoy(size(Zx,dim=1),size(Zx,dim=2)))
        ZxKmoy=Zx
    endif
enddo
if(allocated(Zx)) deallocate(Zx)

! exp. transformation
applyExpKmin=doExp(1)
applyExpKmoy=doExp(2)

end subroutine MAGE_setUp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!==============!
! Private subs !
!==============!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module MAGE_model
