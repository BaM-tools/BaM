module TidalRemenieras_model

!~**********************************************************************
!~* Purpose: Rating curve for twin gauges in tidal rivers
!~           Discharge computation in unsteady and non-uniform flows
!~           Method similar to Remenieras (1949)
!~           Resolution of the polynom Q(i)-Q(i-1)/Dt + A_NP Q(i)^2 + B_NP Q(i) + C_NP = 0
!~**********************************************************************
!~* Programmer: Emeline Perret, INRAE Lyon
!~**********************************************************************
!~* Last modified: 08/06/2020
!~**********************************************************************
!~* Comments: single influenced channel control,
!~*           h1t and h2t have the same, constant time steps
!~*           time lag (Dt) expressed as number of time steps
!~**********************************************************************
!~* References:
!~**********************************************************************
!~* 2Do List:
!~**********************************************************************
!~* Quick description of public procedures:
!~*     1. Tidal2_GetParNumber, number of parameters of the RC
!~*     2. Tidal2_Apply, compute Q=f(h1t(t),h2t(t)|theta)
!~*     3. Tidal2_XtraRead, read optional information - 2 options :
!^*          - (i) keep the water depth and the wetted section at the main station for the calculation
!^*          - (ii) keep the average water depth and wetted section for the calculation (more common)
!~**********************************************************************

use kinds_dmsl_kit ! numeric kind definitions from DMSL

implicit none
Private
public :: TidalRemenieras_GetParNumber,TidalRemenieras_Apply

Contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine TidalRemenieras_GetParNumber(npar,err,mess)

!^**********************************************************************
!^* Purpose: number of parameters of the TidalRemenieras model
!^**********************************************************************
!^* Programmer: Emeline Perret, INRAE Lyon
!^**********************************************************************
!^* Last modified: 12/05/2020
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

err=0;mess='';npar=9 !K,B1,B2,zf1,zf2,deltaz,L,Dt,time_step

end subroutine TidalRemenieras_GetParNumber
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine TidalRemenieras_Apply(IN,theta,OUT,ComputationOption,feas,err,mess)!

!^**********************************************************************
!^* Purpose: apply TidalRemenieras model
!^**********************************************************************
!^* Programmer: Emeline Perret, INRAE Lyon
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
!^*     3. ComputationOption
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
character(250),parameter::procname='Tidal_Apply'
real(mrk),parameter::g=9.81_mrk
real(mrk)::K,B1,B2,zf1,zf2,deltaz,L,time_step
real(mrk)::z1t(size(IN,dim=1)),z2t(size(IN,dim=1)),h1(size(OUT,dim=1)),h2(size(OUT,dim=1)),h12(size(OUT,dim=1)),&
           dAdx(size(OUT,dim=1)),A1(size(OUT,dim=1)),A2(size(OUT,dim=1)),A12(size(OUT,dim=1)),&
           A_NP(size(OUT,dim=1)),B_NP(size(OUT,dim=1)),C_NP(size(OUT,dim=1)),delta_polynom(size(OUT,dim=1))
integer(mik)::t,nobs,Dt,tlag

err=0;mess='';feas=.true.;OUT=undefRN
nobs=size(IN,dim=1)

! Make sense of inputs and theta
z1t=IN(:,1);z2t=IN(:,2)
Dt=theta(1)
K=theta(2)
B1=theta(3)
B2=theta(4)
zf1=theta(5)
zf2=theta(6)
deltaz=theta(7)
L=theta(8)
time_step=theta(9)

! Check feasability
if(K<=0._mrk .or. B1<=0._mrk .or. B2<=0._mrk .or. L<=0._mrk .or. Dt<0._mrk) then
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
    if (z1t(t)-zf1<=0._mrk) then
        OUT(t)=0._mrk !feas(t)=.false.
        cycle
    endif

    h1(tlag)=z1t(tlag)-zf1 !water depth main station 1
    h2(tlag)=z2t(tlag)-zf2 !water depth secondary station 2
    h12(tlag)=(h1(tlag)+h2(tlag))/2 !water depth averaged over the section between station 1 and station 2
    A1(tlag)= B1*h1(tlag) !wetted section main station 1
    A2(tlag)=B2*h2(tlag) !wetted section secondary station 2
    A12(tlag)=(A1(tlag)+A2(tlag))/2 ! wetted section average over the section between station 1 and station 2
    dAdx(tlag)=(A2(tlag)-A1(tlag))/L
    !solve the following polynom : Q(i)-Q(i-1)/Dt + A_NP Q(i)^2 + B_NP Q(i) + C_NP = 0
    A_NP(t)=-1/A12(tlag)**2*(dAdx(tlag))+(g/(K**2*A12(tlag)*(h12(tlag))**(4/3)))
    B_NP(t)=(-2*(A12(tlag)-A12(tlag-1))/A12(tlag)/time_step)
    C_NP(t)=-g*A12(tlag)*(abs(z1t(tlag)-z2t(tlag)-deltaz)/L)
    if (tlag==1) then
        OUT(t)=sqrt(-C_NP(t)/A_NP(t))*sign(1._mrk,z1t(tlag)-z2t(tlag)-deltaz)!permanent
    else
        !calculate discriminent of the polynom
        delta_polynom(t)=(B_NP(t)+(1/time_step))**2-(4*A_NP(t)*(C_NP(t)-abs(OUT(t-1)/time_step)))!
        !Case discriminant<0: not feasable. two solutions (complex)
        if (delta_polynom(t)<0) then
         feas=.false.;return
        endif
        if (delta_polynom(t)==0) then !unique solution (racine double)
        OUT(t)=-(B_NP(t)+(1/time_step))/(2*A_NP(t))
        endif
        if (delta_polynom(t)>0) then
        !Case discriminant>0: two solutions (real)
            if ((z1t(tlag)-z2t(tlag)-deltaz)/L>0) then
                OUT(t)=(-(B_NP(t)+(1/time_step))+sqrt(delta_polynom(t)))/(2*A_NP(t))!solution 1
            else
                OUT(t)=-(-(B_NP(t)+(1/time_step))+sqrt(delta_polynom(t)))/(2*A_NP(t))  !(- solution 1) solution 2 gives nonphysical results
            endif
        endif
    endif
enddo


end subroutine TidalRemenieras_Apply
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module TidalRemenieras_model
