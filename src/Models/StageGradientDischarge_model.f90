module StageGradientDischarge_model

!~**********************************************************************
!~* Purpose: SGD model with Jones' formula for wide rectangular channel
!~**********************************************************************
!~* Programmer: Ben Renard, Irstea Lyon
!~**********************************************************************
!~* Last modified: 27/06/2018
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
public :: SGD_GetParNumber,SGD_Apply

Contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
pure subroutine SGD_GetParNumber(npar,err,mess)

!^**********************************************************************
!^* Purpose: number of parameters
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified: 27/06/2018
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* OUT
!^*		1. npar, par. number
!^*		2.err, error code; <0:Warning, ==0:OK, >0: Error
!^*		3.mess, error message
!^**********************************************************************

integer(mik), intent(out)::npar,err
character(*),intent(out)::mess
! locals
character(250),parameter::procname='SGD_GetParNumber'

err=0;mess='';npar=5

end subroutine SGD_GetParNumber
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine SGD_Apply(IN,theta,OUT,feas,err,mess)

!^**********************************************************************
!^* Purpose: apply SGD
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified:27/06/2018
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*		1. IN, input matrix (h,dhdt)
!^*		2. theta, parameter vector
!^* OUT
!^*		1. OUT, discharge
!^*		2. feas, feasible?
!^*		3. err, error code; <0:Warning, ==0:OK, >0: Error
!^*		4. mess, error message
!^**********************************************************************

real(mrk), intent(in)::IN(:,:),theta(:)
real(mrk), intent(out)::OUT(:)
logical, intent(out)::feas(:)
integer(mik), intent(out)::err
character(*),intent(out)::mess
!locals
character(250),parameter::procname='SGD_Apply'
real(mrk)::K,B,S0,h0,M,Q0,c,h,dhdt,corr
integer(mik)::ntheta,t,nobs

err=0;mess='';feas=.true.;OUT=undefRN
ntheta=size(theta);nobs=size(IN,dim=1)

! make sense of theta
K=theta(1)
B=theta(2)
S0=theta(3)
h0=theta(4)
M=theta(5)

! Check feasability
if(K<=0._mrk .or. M<=0._mrk .or. B<=0._mrk .or. S0<=0._mrk) then
    feas=.false.;return
endif
! compute discharge
do t=1,nobs
    h=IN(t,1);dhdt=IN(t,2)
    ! feasability check
    if(h-h0<=0._mrk) then;feas(t)=.false.;cycle;endif
    ! compute
    Q0=K*B*sqrt(S0)*(h-h0)**M
    c=M*K*sqrt(S0)*(h-h0)**(M-1._mrk)
    corr=1._mrk+(1._mrk/(c*S0))*dhdt
    if(corr<=0._mrk) then;feas(t)=.false.;cycle;endif
    OUT(t)=Q0*sqrt(corr)
enddo

end subroutine SGD_Apply
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!==============!
! Private subs !
!==============!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module StageGradientDischarge_model
