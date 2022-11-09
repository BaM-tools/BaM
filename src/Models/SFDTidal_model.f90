module SFDTidal_model

!~**********************************************************************
!~* Purpose: Rating curve for twin gauges in tidal rivers
!~**********************************************************************
!~* Programmer: Ben Renard, Jerome Le Coz, Emeline Perret, Irstea Lyon
!~**********************************************************************
!~* Last modified: 08/06/2020
!~**********************************************************************
!~* Comments: single influenced channel control,
!~*           h1t and h2t have the same, constant time steps
!~*           time lag (Dt) expressed as number of time steps
!~**********************************************************************
!~* References:
!~**********************************************************************
!~* 2Do List: add slope h1t-h2t as derived parameter
!~*           transition to non-influenced channel control
!~*           extension to variable time steps with Dt real
!~*           correction of slope for highly transient cases (Seine)
!~**********************************************************************
!~* Quick description of public procedures:
!~*     1. SFDTidal_GetParNumber, number of parameters of the RC
!~*     2. SFDTidal_Apply, compute Q=f(h1t(t),h2t(t)|theta)
!~*     3. SFDTidal_XtraRead, read optional information - 2 options :
!^*          - (i) keep the water depth at the main station for the calculation
!^*          - (ii) keep the average water depth for the calculation (more common)
!~**********************************************************************

use kinds_dmsl_kit ! numeric kind definitions from DMSL

implicit none
Private
public :: SFDTidal_GetParNumber,SFDTidal_Apply

Contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
pure subroutine SFDTidal_GetParNumber(npar,err,mess)

!^**********************************************************************
!^* Purpose: number of parameters of the SFDTidal model
!^**********************************************************************
!^* Programmer: Ben Renard, Jerôme Le Coz, Irstea Lyon
!^**********************************************************************
!^* Last modified: 08/06/2020
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^* OUT
!^*     1. npar, par. number
!^*     2. err, error code; <0:Warning, ==0:OK, >0: Error
!^*     3. mess, error message
!^**********************************************************************

integer(mik), intent(out)::npar,err
character(*),intent(out)::mess
! locals

err=0;mess='';npar=8 !Dt,K,B1,B2,h0,M,L,deltaz

end subroutine SFDTidal_GetParNumber
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine SFDTidal_Apply(IN,theta,OUT,ComputationOption,feas,err,mess)

!^**********************************************************************
!^* Purpose: apply SFDTidal model
!^**********************************************************************
!^* Programmer: Ben Renard, Jerôme Le Coz, Emeline Perret, Irstea Lyon
!^**********************************************************************
!^* Last modified: 08/06/2020
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*     1. IN, input matrix (h1,h2)
!^*     2. theta, parameter vector
!^* OUT
!^*     1. OUT, discharge
!^*     2. feas, feasible?
!^*     3. err, error code; <0:Warning, ==0:OK, >0: Error
!^*     4. mess, error message
!^**********************************************************************
real(mrk), intent(in)::IN(:,:),theta(:)
real(mrk), intent(out)::OUT(:)
character(*), intent(in)::ComputationOption
logical, intent(out)::feas(:)
integer(mik), intent(out)::err
character(*),intent(out)::mess
!locals
character(250),parameter::procname='SFDTidal_Apply'
real(mrk)::K,B1,B2,h0,M,L,deltaz
real(mrk)::h1t(size(IN,dim=1)),h2t(size(IN,dim=1))
integer(mik)::t,nobs,Dt,tlag

err=0;mess='';feas=.true.;OUT=undefRN
nobs=size(IN,dim=1)

! Make sense of inputs and theta
h1t=IN(:,1);h2t=IN(:,2)
Dt=theta(1)
K=theta(2)
B1=theta(3)
B2=theta(4)
h0=theta(5)
M=theta(6)
L=theta(7)
deltaz=theta(8)

! Check feasability
if(K<=0._mrk .or. B1<=0._mrk .or. B2<=0._mrk .or. M<=0._mrk .or. L<=0._mrk) then
    feas=.false.;return
endif

! run model
do t=1,nobs
   if (t-Dt<1) then
        tlag=1
    else if (t-Dt>nobs) then
        tlag= nobs
    else
        tlag=t-Dt
    endif
   !computation between station 1 and station 2 and addition of a lag to obtain the discharge in station 1
   OUT(t)=sign(1._mrk,h1t(tlag)-h2t(tlag)-deltaz)*K*(B2+B1)/2* &
   sqrt(abs(h1t(tlag)-h2t(tlag)-deltaz)/L)*((h1t(tlag)+h2t(tlag))/2-h0)**M
   !here h0=h0' the mean bed level of the section
    if (h1t(t)-h0<=0._mrk) then
        OUT(t)=0._mrk !feas(t)=.false.
        cycle
    endif
enddo

end subroutine SFDTidal_Apply
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module SFDTidal_model
