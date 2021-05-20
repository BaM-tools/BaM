module SWOT_model

!~**********************************************************************
!~* Purpose: Rating curve from space using SWOT data
!~**********************************************************************
!~* Programmer: Ben Renard, Irstea Lyon
!~**********************************************************************
!~* Last modified: 29/03/2018
!~**********************************************************************
!~* Comments:
!~**********************************************************************
!~* References:
!~**********************************************************************
!~* 2Do List: 
!~**********************************************************************
!~* Quick description of public procedures:
!~*		1. SWOT_GetParNumber, number of parameters of the RC
!~*		2. SWOT_Apply, compute Q=f(h,S,B?|theta)
!~*		3. SWOT_XtraRead, read Xtra information
!~**********************************************************************

use kinds_dmsl_kit ! numeric kind definitions from DMSL

implicit none
Private
public :: SWOT_GetParNumber,SWOT_Apply,SWOT_XtraRead

Character(100), parameter, PUBLIC:: & ! Catalogue of available rating curves
                    SWOT_hS="SWOT_hS", & ! use only h and S from the satellite as inputs
                    SWOT_hSB="SWOT_hSB" ! use h, S and B from the satellite as inputs

Contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
pure subroutine SWOT_GetParNumber(ID,npar,err,mess)

!^**********************************************************************
!^* Purpose: number of parameters
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified: 29/03/2018
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*		1. ID, ID of the SWOT function
!^* OUT
!^*		1. npar, par. number
!^*		2.err, error code; <0:Warning, ==0:OK, >0: Error
!^*		3.mess, error message
!^**********************************************************************

character(*), intent(in)::ID
integer(mik), intent(out)::npar,err
character(*),intent(out)::mess
! locals
character(250),parameter::procname='SWOT_GetParNumber'

err=0;mess='';npar=undefIN
select case(trim(ID))
case(SWOT_hS,SWOT_hSB)
    npar=5
case default
    err=1;mess=trim(procname)//': Fatal: Unavailable [ID]'
end select

end subroutine SWOT_GetParNumber
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine SWOT_Apply(ID,IN,theta,OUT,feas,err,mess)

!^**********************************************************************
!^* Purpose: apply SWOT model
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified:29/03/2018
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*		1. ID, ID of SWOT formula
!^*		2. IN, input vector (depends on ID)
!^*		3. theta, parameter vector
!^* OUT
!^*		1. OUT, discharge
!^*		2. feas, feasible?
!^*		3. err, error code; <0:Warning, ==0:OK, >0: Error
!^*		4. mess, error message
!^**********************************************************************
character(*), intent(in)::ID
real(mrk), intent(in)::IN(:,:),theta(:)
real(mrk), intent(out)::OUT(:)
logical, intent(out)::feas(:)
integer(mik), intent(out)::err
character(*),intent(out)::mess
!locals
character(250),parameter::procname='SWOT_Apply'
real(mrk)::K,B0,S0,h0,M
real(mrk)::ht(size(IN,dim=1)),St(size(IN,dim=1)),Bt(size(IN,dim=1))
integer(mik)::t,nobs

err=0;mess='';feas=.true.;OUT=undefRN
nobs=size(IN,dim=1)

! Make sense of inputs and theta
ht=IN(:,1);St=IN(:,2)
select case(trim(ID))
case(SWOT_hS)
    Bt=theta(2);B0=0._mrk
case(SWOT_hSB)
    Bt=IN(:,3);B0=theta(2)
case default
    err=1;mess=trim(procname)//': Fatal: Unavailable [ID]'
end select
K=theta(1)
S0=theta(3)
h0=theta(4)
M=theta(5)

! Check feasability
if(K<=0._mrk .or. M<=0._mrk) then
    feas=.false.;return
endif

! run model
do t=1,nobs
    ! feasability checks
    if (Bt(t)-B0<0._mrk .or. St(t)-S0<0._mrk) then;feas(t)=.false.;cycle;endif
    if (ht(t)-h0<=0._mrk) then;OUT(t)=0._mrk;cycle;endif
    OUT(t)=K*(Bt(t)-B0)*sqrt(St(t)-S0)*(ht(t)-h0)**M
enddo
end subroutine SWOT_Apply
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine SWOT_XtraRead(file,xtra,err,mess)

!^**********************************************************************
!^* Purpose: Read Xtra information for SWOT model: model ID 
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified:29/03/2018
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*		1. file, Xtra file
!^* OUT
!^*		1. xtra, xtra information
!^*		2. err, error code; <0:Warning, ==0:OK, >0: Error
!^*		3. mess, error message
!^**********************************************************************
use types_dmsl_kit, only:data_ricz_type
use utilities_dmsl_kit,only:getSpareUnit

character(*), intent(in)::file
type(data_ricz_type),intent(out):: xtra
integer(mik), intent(out)::err
character(*),intent(out)::mess
! locals
character(250),parameter::procname='SWOT_XtraRead'
integer(mik), parameter::nNewtonPar=5
integer(mik)::unt

err=0;mess=''
call getSpareUnit(unt,err,mess)
if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
open(unit=unt,file=trim(file), status='old', iostat=err)
if(err>0) then;mess=trim(procname)//':problem opening file '//trim(file);return;endif
read(unt,*,iostat=err) xtra%cs1 ! model ID
if(err/=0) then;mess=trim(procname)//':problem reading file '//trim(file);return;endif

end subroutine SWOT_XtraRead
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module SWOT_model