module SFDTidal_QmecMS_model

!~**********************************************************************
!~* Purpose: Rating curve for twin gauges in tidal rivers with the Qmec
!~*          model corresponding to the Manning-Strickler equation
!~**********************************************************************
!~* Programmer: Felipe Mendez-Rios,Ben Renard, INRAE
!~**********************************************************************
!~* Last modified: 08/04/2026
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
!~*     1. SFDTidal_QmecMS_GetParNumber, number of parameters of the RC
!~*     2. SFDTidal_QmecMS_Apply, compute Q=f(h1t(t),h2t(t)|theta)
!~**********************************************************************

use kinds_dmsl_kit ! numeric kind definitions from DMSL

implicit none
Private
public :: SFDTidal_QmecMS_GetParNumber,SFDTidal_QmecMS_Apply

Contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
pure subroutine SFDTidal_QmecMS_GetParNumber(npar,err,mess)

!^**********************************************************************
!^* Purpose: number of parameters of the SFDTidal QmecMS model
!^**********************************************************************
!^* Programmer: Felipe Mendez-Rios, Ben Renard, INRAE
!^**********************************************************************
!^* Last modified: 08/04/2026
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

err=0;mess='';npar=10 !y0, A0, be, delta, phi, d1, d2, c, g, dt

end subroutine SFDTidal_QmecMS_GetParNumber
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine SFDTidal_QmecMS_Apply(h1,h2,theta,Q,feas,err,mess)

!^**********************************************************************
!^* Purpose: apply SFDTidal QmecMS model
!^**********************************************************************
!^* Programmer: Felipe Mendez-Rios, Ben Renard, INRAE
!^**********************************************************************
!^* Last modified: 08/04/2026
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
real(mrk), intent(out)::Q(:)
logical, intent(out)::feas(:)
integer(mik), intent(out)::err
character(*),intent(out)::mess
!locals
character(250),parameter::procname='SFDTidal_QmecMS_Apply'
real(mrk)::y0, A0, be, delta, ne, phi, d1, d2, c, g,dx, width, R0
real(mrk)::y1(size(h1)),y2(size(h2)),dy(size(h1)),ym(size(h1)),Ae(size(h1)),Pe(size(h1)),Rhe(size(h1))
integer(mik)::nm, m,k

err=0;mess='';feas=.true.;Q=undefRN
nm=size(h1)

! Check feasability
if(size(h2)/=nm .or. size(Q)/=nm) then
    err=1;mess=trim(procname)//':Fatal: Size mismatch [h1,h2]';feas=.false.;return
endif

! Make sense of inputs and theta
k=0
k=k+1;y0=theta(k)     !typical water level at virtual effective section [m]. Should be fixed by user.
k=k+1;A0=theta(k)     !corresponding typical area A0 [m^2]
k=k+1;be=theta(k)     !effective riverbed elevation [m] (he_qmec = (d1+d2)/2 - be)
k=k+1;delta=theta(k)  !riverbed error leveling [m] (dzeta_qmec = d2 - d1 + delta)
k=k+1;phi=theta(k)    !f/(A0*R0) [2DO: compute unit]
k=k+1;d1=theta(k)     !upstream chart datum [m]
k=k+1;d2=theta(k)     !downstream chart datum [m]
k=k+1;c=theta(k)      !exponent[]. Should in general be fixed to 4/3
k=k+1;g=theta(k)      !gravity[m/s2]
k=k+1;dx=theta(k)     !distance between stations

! Check feasability
if(A0<=0._mrk .or. (y0-be)<=0._mrk .or. phi<=0._mrk .or. c<=0._mrk .or. g<=0._mrk .or. dx<=0._mrk) then
    feas=.false.;return
endif

! Conversion of water level from chart datum level to absolute elevation (Zw)
y1 = h1+d1
y2 = h2+d2

! Intermediate calculations
width=A0/(y0-be)
R0=A0/(width+2._mrk*(y0-be))
ne=sqrt( phi * (1._mrk / (8._mrk*g) ) * A0 * (R0**c) )
dy = y2 - y1           ! Difference between the two stations (dh_qmec + dzeta_qmec = dy + delta) (see pressure gradient calculation)
ym = (y1 + y2)/2._mrk  ! Mean between the two stations (hm_qmec = ym - (d1+d2)/2)

! Geometry parameters for a rectangular cross-section
Ae = width*(ym - be)         ! Mean cross sectional area   (hm_qmec + he_qmec) = ( ym - (d1+d2)/2 + (d1+d2)/2 - be )
Pe = width + 2._mrk*(ym - be) ! Wetted perimeter
Rhe = Ae/Pe                  ! Hydraulic radius

! Check feasability
if(any(Ae<0._mrk) .or. any(Pe<0._mrk) .or. any(Rhe<0._mrk)) then
    feas=.false.;return  ! Consider setting Q=0 (BaRatin) ?
endif

! run model
do m=1,nm
    Q(m) = (1._mrk/ne) * Rhe(m)**(0.5_mrk*c) * Ae(m) * sqrt(abs((dy(m)+delta)/dx))
    if((dy(m)+delta)>0._mrk) Q(m)=-1._mrk*Q(m)
enddo

end subroutine SFDTidal_QmecMS_Apply
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module SFDTidal_QmecMS_model
