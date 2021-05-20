module StageFallDischarge_model

!~**********************************************************************
!~* Purpose: Catalogue of stage-fall-discharge rating curve (variable backwater)
!~**********************************************************************
!~* Programmer: Ben Renard, Irstea Lyon
!~**********************************************************************
!~* Last modified: 08/10/2015
!~**********************************************************************
!~* Comments:
!~**********************************************************************
!~* References:
!~**********************************************************************
!~* 2Do List: 
!~**********************************************************************
!~* Quick description of public procedures:
!~*		1. SFD_GetParNumber, number of parameters of the RC
!~*		2. SFD_Apply, compute Q=f(H1,H2|theta)
!~*		3.
!~**********************************************************************

use kinds_dmsl_kit ! numeric kind definitions from DMSL

implicit none
Private
public :: SFD_GetParNumber,SFD_Apply,SFD_XtraRead

Character(100), parameter, PUBLIC:: & ! Catalogue of available rating curves
                    SFD_Val_General="SFD_Val_General" ! General formalism derived by Val(entin Mansanarez): variable slope formula then standard channel control
                    ! TO DO:
                    !SFD_Val_SameChannel="SFD_Val_SameChannel",& ! restriction of the General formalism by assuming the width and roughness of the channel is identical in both regimes
                    !SFD_POR="SFD_POR" ! Petersen-Overleir and Reitan (2009) formalism

type,public:: SFD_NewtonOptionType
    real(mrk)::upperbound=undefRN,xscale=undefRN,fscale=undefRN,tolX=undefRN,tolF=undefRN
    integer(mik)::itmax=undefIN
end type SFD_NewtonOptionType

Contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
pure subroutine SFD_GetParNumber(ID,npar,err,mess)

!^**********************************************************************
!^* Purpose: number of parameters
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified: 08/10/2015
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*		1. ID, ID of the SFD function
!^* OUT
!^*		1. npar, par. number
!^*		2.err, error code; <0:Warning, ==0:OK, >0: Error
!^*		3.mess, error message
!^**********************************************************************

character(*), intent(in)::ID
integer(mik), intent(out)::npar,err
character(*),intent(out)::mess
! locals
character(250),parameter::procname='SFD_GetParNumber'

err=0;mess='';npar=undefIN
select case(trim(ID))
case(SFD_Val_General)
    npar=8
case default
    err=1;mess=trim(procname)//': Fatal: Unavailable [ID]'
end select

end subroutine SFD_GetParNumber
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine SFD_Apply(ID,IN,theta,NewtonOption,OUT,kappa,feas,err,mess)

!^**********************************************************************
!^* Purpose: apply SFD
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified:08/10/2015
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*		1. ID, ID of SFD formula
!^*		2. IN, input vector (h1,h2)
!^*		3. theta, parameter vector
!^*		4. [NewtonOption], options for Newton numerical resolution
!^* OUT
!^*		1. OUT, discharge
!^*		2. kappa, transition stage
!^*		3. feas, feasible?
!^*		4. err, error code; <0:Warning, ==0:OK, >0: Error
!^*		5. mess, error message
!^**********************************************************************
use types_dmsl_kit,only:data_ricz_type
use utilities_dmsl_kit, only:cleanPointers
use numerix_dmsl_kit, only:znewton

character(*), intent(in)::ID
real(mrk), intent(in)::IN(:,:),theta(:)
type(SFD_NewtonOptionType), intent(in), optional::NewtonOption
real(mrk), intent(out)::OUT(:),kappa(:)
logical, intent(out)::feas(:)
integer(mik), intent(out)::err
character(*),intent(out)::mess
!locals
character(250),parameter::procname='SFD_Apply'
real(mrk),parameter::eps=0.001_mrk
real(mrk)::KB,KBS0,h0,h0prime,M,Mprime,L,delta,fkappa,&
           h1,h2,x1,x2,foo(1)
type(data_ricz_type)::NewtonInfo
integer(mik)::ntheta,fcalls,t,nobs

err=0;mess='';feas=.true.;OUT=undefRN;kappa=undefRN;
ntheta=size(theta);nobs=size(IN,dim=1)

select case(trim(ID))
case(SFD_Val_General)
    if(.not.present(NewtonOption)) then
        err=1;mess=trim(procname)//': Fatal: [NewtonOption] required when [ID]=='//trim(ID)
    endif
    ! make sense of theta
    KB=theta(1)
    h0=theta(2)
    M=theta(3)
    L=theta(4)
    h0prime=theta(5)
    KBS0=theta(6)
    delta=theta(7)
    Mprime=theta(8)
    ! Check feasability
    if(KB<=0._mrk .or. M<=0._mrk .or. L<=0._mrk .or. KBS0<=0._mrk .or. Mprime<=0._mrk) then
        feas=.false.;return
    endif
    do t=1,nobs
        h1=IN(t,1);h2=IN(t,2)
        ! basic checks
        if(h1-h2-delta<=0._mrk) then;feas(t)=.false.;cycle;endif
        if ( (h1-h0<=0._mrk) .or. (h2-h0prime<=0._mrk) ) then;OUT(t)=0._mrk;cycle;endif
        ! Pack parameters and h2 into RICZ structure (for numerical resolution of transition stage k)
        if(associated(NewtonInfo%rp0)) nullify(NewtonInfo%rp0);allocate(NewtonInfo%rp0(ntheta+1))
        NewtonInfo%rp0(1:ntheta)=theta;NewtonInfo%rp0(ntheta+1)=h2
        ! get search interval 
        foo=maxval((/h0,h0prime,h2+delta/))
        x1=foo(1)+eps
        x2=NewtonOption%upperbound !+x1
        ! Apply numerical resolution
        call znewton(evalFunc=fNewton_Val_General,dataIN=NewtonInfo,x1=x1,x2=x2,&
                     tolX=NewtonOption%tolX,tolF=NewtonOption%tolF,xscale=NewtonOption%xscale,&
                     fscale=NewtonOption%fscale,itmax=NewtonOption%itmax,&
                     xroot=kappa(t),froot=fkappa,fcalls=fcalls,err=err,message=mess)
        if(err>0) then
            feas(t)=.false.
            call cleanPointers(NewtonInfo,what=1,err=err,message=mess)            
            cycle
        endif
        ! Does not find root, cycle
        if(fcalls>NewtonOption%itmax) then
            feas(t)=.false.
            call cleanPointers(NewtonInfo,what=1,err=err,message=mess)
            cycle
        endif
        ! Compute discharge
        if(h1<=kappa(t)) then ! variable-slope model
            OUT(t)=(KB*(h1-h0)**M) * sqrt((h1-h2-delta)/L)
        else !standard channel model
            OUT(t)=KBS0*(h1-h0prime)**Mprime 
        endif
        ! All good, but clean RICZ structure to avoid memory leak
        call cleanPointers(NewtonInfo,what=1,err=err,message=mess)
        if(err>0) then;mess=trim(procname)//':'//trim(mess);return;endif
    enddo
case default
    err=1;mess=trim(procname)//': Fatal: Unavailable [ID]'
end select

end subroutine SFD_Apply
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine SFD_XtraRead(file,xtra,err,mess)

!^**********************************************************************
!^* Purpose: Read Xtra information for SFD model: model type 
!^*          & Newton-Raphson parameters
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
character(250),parameter::procname='SFD_XtraRead'
integer(mik), parameter::nNewtonPar=5
integer(mik)::unt

err=0;mess=''
call getSpareUnit(unt,err,mess)
if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
open(unit=unt,file=trim(file), status='old', iostat=err)
if(err>0) then;mess=trim(procname)//':problem opening file '//trim(file);return;endif
read(unt,*,iostat=err) xtra%cs1 ! model ID
if(err/=0) then;mess=trim(procname)//':problem reading file '//trim(file);return;endif
if(associated(xtra%rp1)) nullify(xtra%rp1);allocate(xtra%rp1(nNewtonPar))
read(unt,*,iostat=err) xtra%rp1(1) ! Newton - upper bound
if(err/=0) then;mess=trim(procname)//':problem reading file '//trim(file);return;endif
read(unt,*,iostat=err) xtra%rp1(2) ! Newton - xscale
if(err/=0) then;mess=trim(procname)//':problem reading file '//trim(file);return;endif
read(unt,*,iostat=err) xtra%rp1(3) ! Newton - fscale
if(err/=0) then;mess=trim(procname)//':problem reading file '//trim(file);return;endif
read(unt,*,iostat=err) xtra%rp1(4) ! Newton - xtol
if(err/=0) then;mess=trim(procname)//':problem reading file '//trim(file);return;endif
read(unt,*,iostat=err) xtra%rp1(5) ! Newton - ftol
if(err/=0) then;mess=trim(procname)//':problem reading file '//trim(file);return;endif
read(unt,*,iostat=err) xtra%is1    ! Newton - maxiter
if(err/=0) then;mess=trim(procname)//':problem reading file '//trim(file);return;endif

end subroutine SFD_XtraRead
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!==============!
! Private subs !
!==============!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine fNewton_Val_General(dataIN,dataOUT,x,feas,fx,dfdxV,err,message)
! function to nullify by Newton method for SFD_Val_General
use kinds_dmsl_kit
use types_dmsl_kit,only:data_ricz_type
implicit none
type(data_ricz_type),intent(in),optional::dataIN
type(data_ricz_type),intent(inout),optional::dataOUT
real(mrk),intent(in)::x
logical(mlk),intent(out)::feas
real(mrk),intent(out),optional::fx,dfdxV(:)
integer(mik),intent(out)::err
character(*),intent(out)::message
    
feas=.true.;err=0;message='';
if(present(fx)) fx=undefRN
if(present(dfdxV)) dfdxV=undefRN 
if(present(fx))then
    fx=dataIN%rp0(1)*((x-dataIN%rp0(2))**dataIN%rp0(3))&
       *sqrt((x-dataIN%rp0(9)-dataIN%rp0(7))/dataIN%rp0(4))&
       -dataIN%rp0(6)*((x-dataIN%rp0(5))**dataIN%rp0(8))
endif
if(present(dfdxV))then
    dfdxV(1)=dataIN%rp0(1)*((x-dataIN%rp0(2))**(dataIN%rp0(3)-1))&
             /(2*dataIN%rp0(4)*sqrt((x-dataIN%rp0(9)-dataIN%rp0(7))&
             /dataIN%rp0(4)))*((2*dataIN%rp0(3)+1)*x&
             -(2*dataIN%rp0(3)*(dataIN%rp0(9)+dataIN%rp0(7))+dataIN%rp0(2)))&
             -dataIN%rp0(6)*dataIN%rp0(8)*((x-dataIN%rp0(5))**(dataIN%rp0(8)-1))
endif

end subroutine fNewton_Val_General
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module StageFallDischarge_model
