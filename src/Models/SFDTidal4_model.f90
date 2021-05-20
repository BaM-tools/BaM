module SFDTidal4_model

!~**********************************************************************
!~* Purpose: Rating curve for twin gauges in tidal rivers with smoothing
!~* of the water slope
!~**********************************************************************
!~* Programmer: Ben Renard, Jerome Le Coz, Emeline Perret, Irstea Lyon
!~**********************************************************************
!~* Last modified: 07/08/2019
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
!~*		1. SFDTidal4_GetParNumber, number of parameters of the RC
!~*		2. SFDTidal4_Apply, compute Q=f(h1t(t),h2t(t)|theta)
!~**********************************************************************

use kinds_dmsl_kit ! numeric kind definitions from DMSL

implicit none
Private
public :: SFDTidal4_GetParNumber,SFDTidal4_Apply

Contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
pure subroutine SFDTidal4_GetParNumber(npar,err,mess)

!^**********************************************************************
!^* Purpose: number of parameters of the SFDTidal4 model
!^**********************************************************************
!^* Programmer: Ben Renard, Jerome Le Coz, Emeline Perret, Irstea Lyon
!^**********************************************************************
!^* Last modified: 07/08/2019
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

end subroutine SFDTidal4_GetParNumber
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine SFDTidal4_Apply(IN,theta,OUT,feas,err,mess)

!^**********************************************************************
!^* Purpose: apply SFDTidal4 model
!^**********************************************************************
!^* Programmer: Ben Renard, Jerome Le Coz, Emeline Perret, Irstea Lyon
!^**********************************************************************
!^* Last modified: 07/08/2019
!^**********************************************************************
!^* Comments: smoothing of the water slope calculation
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
!^*     5. Sf(t-delta_t), water slope
!^**********************************************************************
real(mrk), intent(in)::IN(:,:),theta(:)
real(mrk), intent(out)::OUT(:)
logical, intent(out)::feas(:)
integer(mik), intent(out)::err
character(*),intent(out)::mess
!locals
character(250),parameter::procname='SFDTidal4_Apply'
real(mrk)::K,B,h0,M,L,delta
real(mrk)::h1t(size(IN,dim=1)),h2t(size(IN,dim=1)),Sf(size(OUT,dim=1)),Smoy(size(OUT,dim=1))
integer(mik)::t,nobs,Dt,tlag,t1,t2,P

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
P=theta(8) !half-period for averaging water slope

! Check feasability
if(K<=0._mrk .or. B<=0._mrk .or. M<=0._mrk .or. L<=0._mrk .or. P<0._mrk) then
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
    Sf(t)=(h1t(tlag)-h2t(tlag)-delta)/L
enddo
do t=1,nobs
   if (t-P<1 .or. t-P>nobs .or. t+P<1 .or. t+P>nobs) then
       t1=t
       t2=t
   else
       t1=t-P
       t2=t+P
   endif
   Smoy(t)=sum(Sf(t1:t2))/size(Sf(t1:t2))
    if (h1t(t)-h0<=0._mrk) then
        OUT(t)=0._mrk !feas(t)=.false.
        cycle
    endif
    OUT(t)=sign(1._mrk,Smoy(t))*K*B*sqrt(abs(Smoy(t)))*(h1t(t)-h0)**M
enddo
end subroutine SFDTidal4_Apply
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module SFDTidal4_model
