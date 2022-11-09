module BaRatin_model

!~**********************************************************************
!~* Purpose: BaRatin rating curve (Bonnifait Formalism)
!~**********************************************************************
!~* Programmer: Ben Renard, Cemagref Lyon
!~**********************************************************************
!~* Last modified: 21/08/2017
!~**********************************************************************
!~* Comments:
!~**********************************************************************
!~* References:
!~**********************************************************************
!~* 2Do List:
!~**********************************************************************
!~* Quick description of public procedures:
!~*     1. BaRatin_GetParNumber, number of parameters of the RC
!~*     2. BaRatin_Apply, compute Q=f(H|theta)
!~*     3. BaRatin_XtraRead, read Bonnifait matrix
!~**********************************************************************

use kinds_dmsl_kit ! numeric kind definitions from DMSL

implicit none
Private
public :: BaRatin_GetParNumber,BaRatin_Apply,BaRatin_XtraRead,&
          BaRatin_computeQ,CheckControlMatrix

Contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
pure subroutine BaRatin_GetParNumber(ControlMatrix,npar,err,mess)

!^**********************************************************************
!^* Purpose: number of parameters of the RC
!^**********************************************************************
!^* Programmer: Ben Renard, Cemagref Lyon
!^**********************************************************************
!^* Last modified: 21/08/2017
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*     1. ControlMatrix, control Matrix
!^* OUT
!^*     1. npar, par. number
!^*     2.err, error code; <0:Warning, ==0:OK, >0: Error
!^*     3.mess, error message
!^**********************************************************************

integer(mik), intent(in)::ControlMatrix(:,:)
integer(mik), intent(out)::npar,err
character(*),intent(out)::mess
! locals
integer(mik)::ncontrol

err=0;mess='';npar=undefIN

ncontrol=size(ControlMatrix,1) ! it is assumed that the size of ControlMatrix has been checked outside of this sub
npar=3*ncontrol

end subroutine BaRatin_GetParNumber
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine BaRatin_Apply(H,theta,ControlMatrix,Q,b,feas,err,mess)

!^**********************************************************************
!^* Purpose:
!^**********************************************************************
!^* Programmer: Ben Renard, Cemagref Lyon
!^**********************************************************************
!^* Last modified:
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*     1. H, water stage
!^*     2. theta, parameters
!^*     3. ControlMatrix, control Matrix in Laurent's framework
!^* OUT
!^*     1. Q, RC-computed runoff
!^*     2. b, offsets
!^*     3. feas, feasible?
!^*     4.err, error code; <0:Warning, ==0:OK, >0: Error
!^*     5.mess, error message
!^**********************************************************************

real(mrk), intent(in)::H(:), theta(:)
integer(mik), intent(in)::ControlMatrix(:,:)
real(mrk), intent(out)::Q(:),b(:)
logical, intent(out)::feas(:)
integer(mik), intent(out)::err
character(*),intent(out)::mess
!locals
character(250),parameter::procname='BaRatin_Apply'
integer(mik)::nH,ncontrol,i,p,range,t
real(mrk), allocatable::a(:),c(:),k(:)
logical::ok

err=0;mess='';Q=undefRN;feas=.true.

! Get number of H values
nH=size(H)
! Get number of controls (=number of ranges)
ncontrol=size(ControlMatrix,1) ! it is assumed that the size of ControlMatrix has been checked outside of this sub

! check size of theta
if(size(theta)/=3*ncontrol) then
    err=1;mess=trim(procname)//':Fatal:Incorrect size [theta]';feas=.false.;return
endif

! check size of b
if(size(b)/=ncontrol) then
    err=1;mess=trim(procname)//':Fatal:Incorrect size [b]';feas=.false.;return
endif

! Allocate and assign values to k,a,c
allocate(a(ncontrol),c(ncontrol),k(ncontrol))
do i=1,ncontrol
    p=3*(i-1)
    k(i)=theta(p+1);a(i)=theta(p+2);c(i)=theta(p+3)
enddo

! Apply Continuity Constraint & check feasability
call ApplyContinuity(a,c,k,ControlMatrix,b,ok,err,mess)
if(.not.ok) then;feas=.false.;return;endif

! Everything checked OK, can finally apply RC equation!
call BaRatin_computeQ(H,b,a,c,k,ControlMatrix,Q,feas)

end subroutine BaRatin_Apply
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine BaRatin_XtraRead(file,xtra,err,mess)

!^**********************************************************************
!^* Purpose: Read Xtra information for RC: Bonnifait matrix
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified:21/08/2017
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
use utilities_dmsl_kit,only:getSpareUnit,getNumItemsInFile
use types_dmsl_kit, only:data_ricz_type
character(*), intent(in)::file
type(data_ricz_type),intent(out):: xtra
integer(mik), intent(out)::err
character(*),intent(out)::mess
! locals
character(250),parameter::procname='BaRatin_XtraRead'
character(250)::foo
integer(mik)::unt,i,n
logical::feas

err=0;mess=''
call getSpareUnit(unt,err,mess)
if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
open(unit=unt,file=trim(file), status='old', iostat=err)
if(err>0) then;mess=trim(procname)//':problem opening file '//trim(file);return;endif
call getNumItemsInFile(unt=unt,preRwnd=.false.,nskip=0,nitems=n,postPos=0,&
                       jchar=foo,err=err,message=mess)
if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
if(associated(xtra%rpm1)) nullify(xtra%rpm1);allocate(xtra%rpm1(n,n))
do i=1,n
    read(unt,*,iostat=err) xtra%rpm1(i,:)
    if(err/=0) then;mess=trim(procname)//':problem reading file '//trim(file);return;endif
enddo
close(unt)
call CheckControlMatrix(nint(xtra%rpm1),feas,err,mess)
if(err/=0 .or. (.not.feas) ) then
    mess=trim(procname)//':invalid control matrix';
    err=1;return
endif
end subroutine BaRatin_XtraRead
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine BaRatin_computeQ(H,b,a,c,k,ControlMatrix,Q,feas)

!^**********************************************************************
!^* Purpose: Apply RC equation. All checks on parameters are assumed
!^*          already done, and continuity equation has already been applied
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified: 10/01/2018
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*     1. H, water stages
!^*     2. b, offset parameters
!^*     3. a, coefficient parameters
!^*     4. c, exponent parameters
!^*     5. k, activation stages
!^*     6. ControlMatrix, control Matrix in Laurent's framework
!^* OUT
!^*     1. Q, RC-computed runoff
!^*     2. feas, feasible?
!^**********************************************************************

real(mrk), intent(in)::H(:),b(:),a(:),c(:),k(:)
integer(mik), intent(in)::ControlMatrix(:,:)
real(mrk), intent(out)::Q(:)
logical, intent(out)::feas(:)
!locals
integer(mik)::t,i,range

Q=undefRN;feas=.true.

do t=1,size(H)
    ! Retrieve range corresponding to stage H
    range=GetHRange(H(t),k)
    ! Apply RC equation
    Q(t)=0._mrk
    range_loop: do i=1,range ! if range == 0 (h smaller than first activation), loop is inactive and Q=0
        if(ControlMatrix(range,i)==1) then
            if((H(t)-b(i))<0._mrk) then
                Q(t)=0._mrk;exit range_loop
            endif
            Q(t)=Q(t)+a(i)*(H(t)-b(i))**c(i)
        else
            if((H(t)-b(i))<0._mrk) then
                write(*,*) 'WARNING: You should NEVER EVER arrive there!!!'
                write(*,*) 'If you see this message, please contact developpers !!!'
                feas(t)=.false.;exit range_loop
            endif
        endif
    enddo range_loop
enddo

end subroutine BaRatin_computeQ
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
pure subroutine CheckControlMatrix(ControlMatrix,feas,err,mess)

!^**********************************************************************
!^* Purpose: Check Control Matrix is OK
!^**********************************************************************
!^* Programmer: Laurent Bonnifait & Ben Renard, Cemagref Lyon
!^**********************************************************************
!^* Last modified:
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*     1. ControlMatrix
!^* OUT
!^*     1.feas
!^*     2.err, error code; <0:Warning, ==0:OK, >0: Error
!^*     3.mess, error message
!^**********************************************************************

integer(mik), intent(in)::ControlMatrix(:,:)
integer(mik), intent(out)::err
character(*),intent(out)::mess
logical(mlk), intent(out)::feas
! locals
character(250),parameter::procname='CheckControlMatrix'
integer(mik)::ncontrol,nrange,i
logical(mlk)::mask(size(ControlMatrix,dim=2))

err=0;mess='';feas=.true.

ncontrol=size(ControlMatrix,dim=2)
nrange=size(ControlMatrix,dim=1)

if(ncontrol/=nrange) then
    err=1;mess=trim(procname)//':FATAL:ncontrol/=nrange';feas=.false.;return
endif

if (.not. all( ControlMatrix==0 .or. ControlMatrix==1) ) then
    err=2;mess=trim(procname)//':FATAL:Some elements of ControlMatrix differ from 0/1';feas=.false.;return
endif

if(nrange>1) then ! check matrix is lower-diagonal
    do i=1,nrange-1
        if( any(ControlMatrix(i,i+1:)>0) ) then
            err=3;mess=trim(procname)//':FATAL:ControlMatrix should be lower-diagonal';feas=.false.;return
        endif
    enddo
endif

if(nrange>1) then ! check each change of range activates exactly one control
    do i=2,nrange
        mask=ControlMatrix(i,:) - ControlMatrix(i-1,:)>0
        if (count(mask)/=1) then
            err=4;mess=trim(procname)//':FATAL: Exactly one control should be activated when changing range';feas=.false.;return
        endif
    enddo
endif

do i=1,nrange ! check that change to the kth range activates the kth control
    if(ControlMatrix(i,i)/=1) then
        err=5;mess=trim(procname)//':FATAL: change to the kth range should activate the kth control';feas=.false.;return
    endif
enddo

end subroutine CheckControlMatrix
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!==============!
! Private subs !
!==============!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
pure subroutine ApplyContinuity(a,c,k,ControlMatrix,b,feas,err,mess)

!^**********************************************************************
!^* Purpose: Apply Continuity constraints & evaluate feasability
!^**********************************************************************
!^* Programmer: Laurent Bonnifait & Ben Renard, Cemagref Lyon
!^**********************************************************************
!^* Last modified:
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*     1.a
!^*     2.c
!^*     3.k
!^*     4.ControlMatrix
!^* OUT
!^*     1.feas, feasability
!^*     2.err, error code; <0:Warning, ==0:OK, >0: Error
!^*     3.mess, error message
!^* INOUT
!^*     1.b
!^**********************************************************************

real(mrk), intent(in)::a(:),c(:),k(:)
integer(mik),intent(in):: ControlMatrix(:,:)
real(mrk), intent(inout)::b(:)
logical(mlk), intent(out)::feas
integer(mik), intent(out)::err
character(*),intent(out)::mess
! locals
integer(mik)::ncontrol,i,j
real(mrk)::som

err=0;mess='';feas=.true.

ncontrol=size(a)

if(a(1)<=0._mrk .or. c(1)==0._mrk) then
    feas=.false.;return
endif
b(1)=k(1)
if(ncontrol<=1) return ! nothing more to do

do j=2,ncontrol
    ! feasability checks
    if(k(j)<=k(j-1)) then
        feas=.false.;return
    endif
    if( any( k(j)<b(1:(j-1)) ) .or. a(j)<=0._mrk .or. c(j)==0._mrk) then
        feas=.false.;return
    endif
    som=0._mrk
    do i=1,j-1
        som=som+real(ControlMatrix(j-1,i)-ControlMatrix(j,i),mrk)*a(i)*(k(j)-b(i))**c(i)
    enddo
    if(som<0._mrk) then
        feas=.false.;return
    endif
    b(j)=k(j)-((1._mrk/a(j))*som)**(1._mrk/c(j))
enddo

end subroutine ApplyContinuity
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function GetHRange(H,k)

!#**********************************************************************
!#* Purpose: find the range where H belongs
!#**********************************************************************
!#* Programmer: Laurent Bonnifait & Ben Renard, Cemagref Lyon
!#**********************************************************************
!#* Last modified:
!#**********************************************************************
!#* Comments:
!#**********************************************************************
!#* References:
!#**********************************************************************
!#* 2Do List:
!#**********************************************************************
!#* IN
!#*     1. H
!#*     2. k
!#* OUT
!#*     1.RC_General_GetRange
!#**********************************************************************

real(mrk), intent(in)::H,k(:)
integer(mik)::GetHRange
!locals
integer(mik)::nk,i

nk=size(k)
do i=1,nk
    If(H<k(i)) then
        GetHRange=i-1;return
    else
        cycle
    endif
enddo

GetHRange=nk

end function GetHRange
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module BaRatin_model
