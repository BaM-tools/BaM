module SFDTidal2_model

!~**********************************************************************
!~* Purpose: Rating curve for twin gauges in tidal rivers with correction
!~* of water slope based on the estimation of the water surface temporal
!~* variations at a fictive station (2') closer to the primary station (1)
!~**********************************************************************
!~* Programmer: Ben Renard, Jerome Le Coz, Emeline Perret, Irstea Lyon
!~**********************************************************************
!~* Last modified: 01/08/2019
!~**********************************************************************
!~* Comments: single influenced channel control,
!~*           h1t and h2t have the same, constant time steps
!~*           time lag (Dt) expressed as number of time steps and Dt>0
!~**********************************************************************
!~* References:
!~**********************************************************************
!~* 2Do List: add slope h1t-h2t as derived parameter
!~*           transition to non-influenced channel control
!~*           extension to variable time steps with Dt real
!~*           correction of slope for highly transient cases (Seine)
!~**********************************************************************
!~* Quick description of public procedures:
!~*		1. SFDTidal2_GetParNumber, number of parameters of the RC
!~*		2. SFDTidal2_Apply, compute Q=f(h1t(t),h2t(t)|theta)
!~**********************************************************************

use kinds_dmsl_kit ! numeric kind definitions from DMSL

implicit none
Private
public :: SFDTidal2_GetParNumber,SFDTidal2_Apply

Contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
pure subroutine SFDTidal2_GetParNumber(npar,err,mess)

!^**********************************************************************
!^* Purpose: number of parameters of the SFDTidal2 model
!^**********************************************************************
!^* Programmer: Ben Renard, Jerome Le Coz, Emeline Perret, Irstea Lyon
!^**********************************************************************
!^* Last modified: 01/08/2019
!^**********************************************************************
!^* Comments: +1 parameter compared to the previous model SFDTidal_model
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^* OUT
!^*		1. npar, par. number
!^*		2. err, error code; <0:Warning, ==0:OK, >0: Error
!^*		3. mess, error message
!^**********************************************************************

integer(mik), intent(out)::npar,err
character(*),intent(out)::mess
! locals

err=0;mess='';npar=8

end subroutine SFDTidal2_GetParNumber
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine SFDTidal2_Apply(IN,theta,OUT,feas,err,mess)

!^**********************************************************************
!^* Purpose: apply SFDTidal2 model
!^**********************************************************************
!^* Programmer: Ben Renard, Jerome Le Coz, Emeline Perret, Irstea Lyon
!^**********************************************************************
!^* Last modified: 01/08/2019
!^**********************************************************************
!^* Comments: correction for the water slope calculation
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*		1. IN, input matrix (h1,h2)
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
character(250),parameter::procname='SFDTidal2_Apply'
real(mrk)::K,B,h0,M,L,alpha,delta
real(mrk)::h1t(size(IN,dim=1)),h2t(size(IN,dim=1))
integer(mik)::t,nobs,Dt,tlag

err=0;mess='';feas=.true.;OUT=undefRN
nobs=size(IN,dim=1)

! Make sense of inputs and theta
h1t=IN(:,1);h2t=IN(:,2)
Dt=theta(1)
K=theta(2)
B=theta(3)
h0=theta(4)
M=theta(5)
L=theta(6)
delta=theta(7)
alpha=theta(8)

! Check feasability
if(K<=0._mrk .or. B<=0._mrk .or. M<=0._mrk .or. L<=0._mrk .or. alpha<=0._mrk .or. alpha>=1._mrk) then
    feas=.false.;return
endif

! run model
do t=1,nobs
    ! feasability checks
    if (t-Dt<1) then
        tlag=1
    else if (t-Dt>nobs) then
        tlag=nobs
    else
        tlag=t-Dt
    endif
    if (h1t(t)-h0<=0._mrk) then
        OUT(t)=0._mrk !feas(t)=.false.
        cycle
    endif
    OUT(t)=sign(1._mrk,h1t(t)-alpha*h2t(tlag)-delta)*K*B*sqrt(abs(h1t(t)-alpha*h2t(tlag)-delta)/L)*(h1t(t)-h0)**M
enddo
end subroutine SFDTidal2_Apply
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module SFDTidal2_model
