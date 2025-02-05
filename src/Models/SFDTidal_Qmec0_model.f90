module SFDTidal_Qmec0_model

!~**********************************************************************
!~* Purpose: Rating curve for twin gauges in tidal rivers with the
!~*          original Qmec model
!~**********************************************************************
!~* Programmer: Felipe Mendez-Rios,Ben Renard, INRAE
!~**********************************************************************
!~* Last modified: 05/02/2025
!~**********************************************************************
!~* Comments: single influenced rectangular channel control,
!~*           h1t and h2t have the same, constant time steps
!~**********************************************************************
!~* References:
!~**********************************************************************
!~* 2Do List:
!~**********************************************************************
!~* Quick description of public procedures:
!~*     1. SFDTidal_Qmec0_GetParNumber, number of parameters of the RC
!~*     2. SFDTidal_Qmec0_Apply, compute Q=f(h1t(t),h2t(t)|theta)
!~**********************************************************************

use kinds_dmsl_kit ! numeric kind definitions from DMSL

implicit none
Private
public :: SFDTidal_Qmec0_GetParNumber,SFDTidal_Qmec0_Apply

Contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
pure subroutine SFDTidal_Qmec0_GetParNumber(npar,err,mess)

!^**********************************************************************
!^* Purpose: number of parameters of the SFDTidal Qmec model
!^**********************************************************************
!^* Programmer: Felipe Mendez-Rios, Ben Renard, INRAE
!^**********************************************************************
!^* Last modified: 22/01/2024
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

err=0;mess='';npar=11 !Be, he, dzeta, ne, d1, d2, c, g, Q0, dx, dt

end subroutine SFDTidal_Qmec0_GetParNumber
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine SFDTidal_Qmec0_Apply(h1,h2,theta,Q,feas,err,mess)

!^**********************************************************************
!^* Purpose: apply SFDTidal Qmec model
!^**********************************************************************
!^* Programmer: Felipe Mendez-Rios, Ben Renard, INRAE
!^**********************************************************************
!^* Last modified: 22/01/2024
!^**********************************************************************
!^* Comments:
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
real(mrk), intent(out)::Q(:)
logical, intent(out)::feas(:)
integer(mik), intent(out)::err
character(*),intent(out)::mess
!locals
character(250),parameter::procname='SFDTidal_Qmec0_Apply'
real(mrk)::Be, he, dzeta, ne, d1, d2, c, g,Q0, dx, dt
real(mrk)::dhdt1, dhdt2, pressure_gradient, bottom_friction, advection
real(mrk)::y1(size(h1)),y2(size(h2)),dh(size(h1)),hm(size(h1)),Ae(size(h1)),Pe(size(h1)),Rhe(size(h1))
integer(mik)::nm, m

err=0;mess='';feas=.true.;Q=undefRN
nm=size(h1)

! Check feasability
if(size(h2)/=nm .or. size(Q)/=nm) then
    err=1;mess=trim(procname)//':Fatal: Size mismatch [h1,h2]';feas=.false.;return
endif

! Make sense of inputs and theta
Be=theta(1)    !effective width B_e [m]
he=theta(2)    !effective depth h_e [m]
dzeta=theta(3) !hydraulic slope dzeta [m]
ne=theta(4)    !effective Manning's coefficient [s m^(-1/3)]
d1=theta(5)    !upstream offset [m]
d2=theta(6)    !downstream offset [m]
c=theta(7)     !exponent[]
g=theta(8)     !gravity[m/s2]
Q0=theta(9)    !initial discharge [m3/s]
dx=theta(10)    !distance between stations
dt=theta(11)   !time step

! Check feasability
if(Be<=0._mrk .or. he<=0._mrk .or. dzeta>=0._mrk .or. ne<=0._mrk .or. c<=0._mrk .or. g<=0._mrk .or. dx<=0._mrk .or. dt<=0._mrk) then
    feas=.false.;return
endif

! water levels minus datums
y1 = h1-d1
y2 = h2-d2

! Intermediate calculations
dh = y2 - y1           ! Difference between the two stations
hm = (y1 + y2)/2._mrk  ! Mean between the two stations

! Geometry parameters for a rectangular cross-section
Ae = Be*(he + hm)     ! Mean cross sectional area
Pe = Be + 2._mrk*(he + hm) ! Wetted perimeter
Rhe = Ae/Pe           ! Hydraulic radius

! Check feasability
if(any(Ae<0._mrk) .or. any(Pe<0._mrk) .or. any(Rhe<0._mrk)) then
    feas=.false.;return
endif

! run model
Q(1)=Q0
do m=2,nm
    if (m==2) then
        !First-order scheme for the first time-step with Q(1) = Q0
        dhdt1 = (hm(m) - hm(m-1))/(dt)
        pressure_gradient = -g*(Ae(m-1)*(dh(m-1) + dzeta)/dx)

        bottom_friction   = -g*((ne**2)/(Rhe(m-1)**(c)))*Q(m-1)*abs(Q(m-1))/Ae(m-1)

        advection         = 2._mrk*(Q(m-1)/Ae(m-1))*Be*dhdt1

    else
        dhdt1 = (hm(m) - hm(m-2))/(2*dt)
        if (m == 3) then
            dhdt2 = (hm(m-1) - hm(m-2))/dt
        else
            dhdt2 = (hm(m-1) - hm(m-3))/(2*dt)
        endif
        !Second-order Adams-Bashforth scheme

        pressure_gradient = (1.5_mrk)*(-g*Ae(m-1)*(dh(m-1) + dzeta)/dx) &
                           -(0.5_mrk)*(-g*Ae(m-2)*(dh(m-2) + dzeta)/dx)

        bottom_friction   = (1.5_mrk)*(-g*((ne**2)/(Rhe(m-1)**(c)))*Q(m-1)*abs(Q(m-1))/Ae(m-1))&
                           -(0.5_mrk)*(-g*((ne**2)/(Rhe(m-2)**(c)))*Q(m-2)*abs(Q(m-2))/Ae(m-2))

        advection         = (1.5_mrk)*(2._mrk*(Q(m-1)/Ae(m-1))*Be*dhdt1) &
                           -(0.5_mrk)*(2._mrk*(Q(m-2)/Ae(m-2))*Be*dhdt2)
    endif

    Q(m) = Q(m-1) + dt*(pressure_gradient + bottom_friction + advection)
enddo

end subroutine SFDTidal_Qmec0_Apply
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module SFDTidal_Qmec0_model

