module SFDTidal_Sw_correction_model

!~**********************************************************************
!~* Purpose: Rating curve for twin gauges in tidal rivers
!~*          Model based on SFDTidal but with a water slope correction
!~*          for highly transient cases (ex: Seine).
!~*          During the time when the tidal wave front moves across the
!~*          channel control, the slope computation is stopped and
!~*          the discharge is replaced by a linear interpolation
!~**********************************************************************
!~* Programmer: Emeline Perret, INRAE Lyon
!~**********************************************************************
!~* Last modified: 18/06/2021
!~**********************************************************************
!~* Comments: single influenced channel control,
!~*           h1t and h2t have the same, constant time steps
!~*           time lag (Dt) expressed as number of time steps
!~**********************************************************************
!~* References:
!~**********************************************************************
!~* 2Do List:
!~*           transition to non-influenced channel control
!~*           extension to variable time steps with Dt real
!~**********************************************************************
!~* Quick description of public procedures:
!~*     1. SFDTidal_Sw_correction_GetParNumber, number of parameters of the RC
!~*     2. SFDTidal_Sw_correction_Apply, compute Q=f(h1t(t),h2t(t)|theta) with the correction
!~**********************************************************************

use kinds_dmsl_kit ! numeric kind definitions from DMSL

implicit none
Private
public :: SFDTidal_Sw_correction_GetParNumber,SFDTidal_Sw_correction_Apply

Contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
pure subroutine SFDTidal_Sw_correction_GetParNumber(npar,err,mess)

!^**********************************************************************
!^* Purpose: number of parameters of the SFDTidal_Sw_correction model
!^**********************************************************************
!^* Programmer: Emeline Perret, INRAE Lyon
!^**********************************************************************
!^* Last modified: 18/06/2021
!^**********************************************************************
!^* Comments: addition of two parameters compared to SFDTidal_model
!^*           corresponding to two thresholds for stopping slope computation
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

err=0;mess='';npar=11 !Dt,K,B1,B2,h0,M,L,deltaz,time_step,th1,th2

end subroutine SFDTidal_Sw_correction_GetParNumber
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine SFDTidal_Sw_correction_Apply(IN,theta,Q,interp,feas,err,mess)

!^**********************************************************************
!^* Purpose: apply SFDTidal model
!^**********************************************************************
!^* Programmer: Emeline Perret, INRAE Lyon
!^**********************************************************************
!^* Last modified: 23/06/2021
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
!^*     1. Q, discharge
!^*     2. interp, variable indicating where the discharge linear interp is made
!^*     2. feas, feasible?
!^*     3. err, error code; <0:Warning, ==0:OK, >0: Error
!^*     4. mess, error message
!^**********************************************************************
real(mrk), intent(in)::IN(:,:),theta(:)
real(mrk), intent(out)::Q(:),interp(:)
logical, intent(out)::feas(:)
integer(mik), intent(out)::err
character(*),intent(out)::mess
!locals
character(250),parameter::procname='SFDTidal_Apply'
real(mrk)::K,B1,B2,h0,M,L,deltaz,time_step,th1,th2,slopechange(size(IN,dim=1)),Qmodel(size(IN,dim=1))
real(mrk)::h1t(size(IN,dim=1)),h2t(size(IN,dim=1)),dh1dt(size(IN,dim=1)),dh2dt(size(IN,dim=1))
integer(mik)::t,nobs,Dt,tlag,i,j,t1,t2,start_interp(size(IN,dim=1)),end_interp(size(IN,dim=1))

!Init
err=0;mess='';feas=.true.;Q=undefRN;Qmodel=undefRN;interp=undefRN;dh1dt=0;dh2dt=0;slopechange=0;i=1
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
time_step=theta(9)
th1=theta(10)
th2=theta(11)

! Check feasability
if(K<=0._mrk .or. B1<=0._mrk .or. B2<=0._mrk .or. M<=0._mrk .or. L<=0._mrk .or. th1<=0._mrk .or. th2<=0._mrk) then
    feas=.false.;return
endif

! 1 - run model >> Qmodel without correction
do t=1,nobs
   if (t-Dt<1) then
        tlag=1
    else if (t-Dt>nobs) then
        tlag= nobs
    else
        tlag=t-Dt
    endif
   !computation between station 1 and station 2 and addition of a lag to obtain the discharge in station 1
   Qmodel(t)=sign(1._mrk,h1t(tlag)-h2t(tlag)-deltaz)*K*(B2+B1)/2* &
   sqrt(abs(h1t(tlag)-h2t(tlag)-deltaz)/L)*((h1t(tlag)+h2t(tlag))/2-h0)**M
   !here h0=h0' the mean bed level of the section
    if (h1t(t)-h0<=0._mrk) then
        Qmodel(t)=0._mrk !feas(t)=.false.
        cycle
    endif
enddo

!2 - find the times where we need to stop the computation because the water slope computation is not appropriate
!compute dhdt
dh1dt(1)=0
dh2dt(1)=0
do t=2,nobs
  dh1dt(t)=(h1t(t)-h1t(t-1))/time_step
  dh2dt(t)=(h2t(t)-h2t(t-1))/time_step
enddo
!interp=1 : stop the computation and interpolate
do t=1,nobs
   if (t>1) then
        slopechange(t)=sign(1._mrk,(h1t(t)-h2t(t)-deltaz))-sign(1._mrk,(h1t(t-1)-h2t(t-1)-deltaz))
        if (slopechange(t)<0._mrk .and. sign(1._mrk,(h1t(t)-h2t(t)-deltaz))/=0) then !.and. dh2dt(t)>(th2*abs(dh1dt(t)))
             interp(t)=1
             start_interp(i)=t-1
             i=i+1
        else if (interp(t-1)==1 .and. dh2dt(t)>(th2*abs(dh1dt(t)))) then
             interp(t)=1
             end_interp(i-1)=t+1
        else if (interp(t-1)==1 .and. dh1dt(t)>(th1*abs(dh2dt(t)))) then
             interp(t)=1
             end_interp(i-1)=t+1
        else if (interp(t-1)==1 .and. (abs(dh1dt(t))-abs(dh2dt(t)))>=-1.479936e-17 &
            .and. (abs(dh1dt(t))-abs(dh2dt(t)))<=1.479936e-17) then
             interp(t)=1
             end_interp(i-1)=t+1
        else
             interp(t)=undefRN
        endif
   else
       interp(t)=undefRN
   endif
enddo


! 3 - Linear interpolation for the discharge during the stopping time : Qinterpolated
! & - Calculation of the output (a mix of Qmodel and Qinterpolated)
do j=1,size(start_interp,dim=1)
t1=start_interp(j)+Dt ! le lag a appliquer ici?
t2=end_interp(j)+Dt
do t=1,nobs
    if (t>t1 .and. t<t2) then
        Q(t)=Qmodel(t1)+((Qmodel(t2)-Qmodel(t1))*(t-t1)/(t2-t1)) !Qmodel(t1)+(a*(Qmodel(t2)-Qmodel(t1)))  a=(t-t1)/(t2-t1)
    else
        Q(t)=Qmodel(t)
    endif
enddo
Qmodel=Q
enddo

end subroutine SFDTidal_Sw_correction_Apply
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module SFDTidal_Sw_correction_model
