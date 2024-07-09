module SFDTidal_Qmec2_model

!~**********************************************************************
!~* Purpose: Rating curve for twin gauges in tidal rivers with the Qmec2 model
!~**********************************************************************
!~* Programmer: Felipe Mendez-Rios,Ben Renard, INRAE
!~**********************************************************************
!~* Last modified: 22/01/2024
!~**********************************************************************
!~* Comments: single influenced rectangular channel control,
!~*           h1t and h2t have the same, constant time steps
!~*           This is a modified version of the original Qmec model
!~**********************************************************************
!~* References:
!~**********************************************************************
!~* 2Do List:
!~**********************************************************************
!~* Quick description of public procedures:
!~*     1. SFDTidal_Qmec_GetParNumber, number of parameters of the RC
!~*     2. SFDTidal_Qmec_Apply, compute Q=f(h1t(t),h2t(t)|theta)
!~**********************************************************************

use kinds_dmsl_kit ! numeric kind definitions from DMSL

implicit none
Private
public :: SFDTidal_Qmec2_GetParNumber,SFDTidal_Qmec2_Apply

Contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
pure subroutine SFDTidal_Qmec2_GetParNumber(npar,err,mess)

!^**********************************************************************
!^* Purpose: number of parameters of the SFDTidal Qmec2 model
!^**********************************************************************
!^* Programmer: Felipe Mendez-Rios, Ben Renard, INRAE
!^**********************************************************************
!^* Last modified: 08/07/2024
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

err=0;mess='';npar=11 !Be, bevs,delta, ne, d1, d2, c, g, Q0, dx, dt

end subroutine SFDTidal_Qmec2_GetParNumber
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine SFDTidal_Qmec2_Apply(h1,h2,theta,Q,pressure_gradient,bottom_friction,advection,feas,err,mess)

!^**********************************************************************
!^* Purpose: apply SFDTidal Qmec2 model
!^**********************************************************************
!^* Programmer: Felipe Mendez-Rios, Ben Renard, INRAE
!^**********************************************************************
!^* Last modified: 08/07/2024
!^**********************************************************************
!^* Comments:
!^* Qmec original is calibrated using stages (h) from chart datum (d).
!^* Qmec BaM is calibrated using stages from sea level (y or Zw).
!^* The relation both variable is y = h + d
!^**********************************************************************
!^* References: 10.1029/2019JC015992
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*     1. h1, upstream stage time series
!^*     2. h2, downstream stage time series
!^*     3. theta, parameter vector
!^* OUT
!^*     1. Q, discharge
!^*     2. feas, feasible?
!^*     3. err, error code; <0:Warning, ==0:OK, >0: Error
!^*     4. mess, error message
!^**********************************************************************
real(mrk), intent(in)::h1(:),h2(:),theta(:)
real(mrk), intent(out)::Q(:),pressure_gradient(:), bottom_friction(:), advection(:)
logical, intent(out)::feas(:)
integer(mik), intent(out)::err
character(*),intent(out)::mess
!locals
character(250),parameter::procname='SFDTidal_Qmec2_Apply'
real(mrk)::Be, bevs, delta, ne, d1, d2, c, g,Q0, dx, dt, pg, bf, ad
real(mrk)::dydt1, dydt2
real(mrk)::y1(size(h1)),y2(size(h2)),dy(size(h1)),ym(size(h1)),Ae(size(h1)),Pe(size(h1)),Rhe(size(h1))
integer(mik)::nm, m

err=0;mess='';feas=.true.;Q=undefRN
pressure_gradient=undefRN;bottom_friction=undefRN;advection=undefRN
nm=size(h1)

! Check feasability
if(size(h2)/=nm .or. size(Q)/=nm) then
    err=1;mess=trim(procname)//':Fatal: Size mismatch [h1,h2]';feas=.false.;return
endif

! Make sense of inputs and theta
Be=theta(1)    !effective width B_e [m]
bevs=theta(2)  !effective riverbed elevation from sea level [m] (he_qmec = -bevs)
delta=theta(3) !riverbed error leveling [m] (dzeta_qmec = d2 - d1 + delta)
ne=theta(4)    !effective Manning's coefficient [s m^(-1/3)]
d1=theta(5)    !upstream chart datum [m]
d2=theta(6)    !downstream chart datum [m]
c=theta(7)     !exponent[]
g=theta(8)     !gravity[m/s2]
Q0=theta(9)    !initial discharge [m3/s]
dx=theta(10)   !distance between stations
dt=theta(11)   !time step

! Check feasability
if(Be<=0._mrk .or. ne<=0._mrk .or. c<=0._mrk .or. g<=0._mrk .or. dx<=0._mrk .or. dt<=0._mrk) then
    feas=.false.;return
endif

! Conversion of water level from chart datum level to absolute elevation (Zw)
y1 = h1+d1
y2 = h2+d2

! Intermediate calculations
dy = y2 - y1           ! Difference between the two stations   (dh_qmec + dzeta_qmec = dy + delta) (see pressure gradient calculation)
ym = (y1 + y2)/2._mrk  ! Mean between the two stations (hm_qmec = ym - (d1+d2)/2)

! Geometry parameters for a rectangular cross-section
Ae = Be*(ym - bevs)          ! Mean cross sectional area
Pe = Be + 2._mrk*(ym - bevs) ! Wetted perimeter
Rhe = Ae/Pe                  ! Hydraulic radius

! Check feasability
if(any(Ae<0._mrk) .or. any(Pe<0._mrk) .or. any(Rhe<0._mrk)) then
    feas=.false.;return  ! Consider setting Q=0 (BaRatin) ?
endif

! run model
Q(1)=Q0
pressure_gradient(1)=0
bottom_friction(1)=0
advection(1)=0

do m=2,nm
    if (m==2) then
        !First-order scheme for the first time-step with Q(1) = Q0
        dydt1 = (ym(m) - ym(m-1))/(dt)

        pg = -g*(Ae(m-1)*(dy(m-1)+delta)/dx) ! (dh_qmec + dzeta_qmec = dy + delta)

        bf = -g*((ne**2)/(Rhe(m-1)**(c)))*Q(m-1)*abs(Q(m-1))/Ae(m-1)

        ad = 2._mrk*(Q(m-1)/Ae(m-1))*Be*dydt1

    else
        dydt1 = (ym(m) - ym(m-2))/(2*dt)
        if (m == 3) then
            dydt2 = (ym(m-1) - ym(m-2))/dt
        else
            dydt2 = (ym(m-1) - ym(m-3))/(2*dt)
        endif
        !Second-order Adams-Bashforth scheme

        pg = (1.5_mrk)*(-g*Ae(m-1)*(dy(m-1)+delta)/dx) &
                           -(0.5_mrk)*(-g*Ae(m-2)*(dy(m-2)+delta)/dx)

        bf = (1.5_mrk)*(-g*((ne**2)/(Rhe(m-1)**(c)))*Q(m-1)*abs(Q(m-1))/Ae(m-1))&
                           -(0.5_mrk)*(-g*((ne**2)/(Rhe(m-2)**(c)))*Q(m-2)*abs(Q(m-2))/Ae(m-2))

        ad = (1.5_mrk)*(2._mrk*(Q(m-1)/Ae(m-1))*Be*dydt1) &
                           -(0.5_mrk)*(2._mrk*(Q(m-2)/Ae(m-2))*Be*dydt2)
    endif

    Q(m) = Q(m-1) + dt*(pg+bf+ad)
    pressure_gradient(m) = dt*pg
    bottom_friction(m) = dt*bf
    advection(m) = dt*ad

enddo

end subroutine SFDTidal_Qmec2_Apply
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module SFDTidal_Qmec2_model
