module AlgaeBiomass_model

!~**********************************************************************
!~* Purpose: Growth/decline of biomass for algae
!~**********************************************************************
!~* Programmer: Ben Renard, Irstea Lyon - University of Adelaide
!~**********************************************************************
!~* Created: 01/08/2019
!~**********************************************************************
!~* Comments:
!~**********************************************************************
!~* References:
!~**********************************************************************
!~* 2Do List:
!~**********************************************************************
!~* Quick description of public procedures:
!~*     1. AlgaeBiomass_GetParNumber, number of parameters
!~*     2. AlgaeBiomass_Apply, apply biomass model
!~**********************************************************************

use kinds_dmsl_kit ! numeric kind definitions from DMSL

implicit none
Private
public :: AlgaeBiomass_GetParNumber,AlgaeBiomass_Apply

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
pure subroutine AlgaeBiomass_GetParNumber(npar,err,mess)
!^**********************************************************************
!^* Purpose: number of parameters of the model
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon - University of Adelaide
!^**********************************************************************
!^* Created: 10/09/2019
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* OUT
!^*     1. npar, par. number
!^*     2.err, error code; <0:Warning, ==0:OK, >0: Error
!^*     3.mess, error message
!^**********************************************************************
integer(mik), intent(out)::npar,err
character(*),intent(out)::mess
! locals
err=0;mess='';npar=18
end subroutine AlgaeBiomass_GetParNumber
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine AlgaeBiomass_Apply(temp,irradiance,depth,turbidity,removal,& ! inputs: water temperature, irradiance at the surface, water depth, turbidity, removal
                              theta,& ! parameters
                              biomass,biomass0,Fb,Ft,Fi,Fn,Rhyd,& ! outputs
                              feas,err,mess)
!^**********************************************************************
!^* Purpose: Biomass growth and decline
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon - University of Adelaide
!^**********************************************************************
!^* Last modified:
!^* - 10/09/2019 by E. Perret > modification of the limiting function Fi
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*     1. XXX
!^* OUT
!^*     1. XXX
!^*     3. feas, feasible?
!^*     4.err, error code; <0:Warning, ==0:OK, >0: Error
!^*     5.mess, error message
!^**********************************************************************
real(mrk), intent(in)::temp(:),irradiance(:),depth(:),turbidity(:),removal(:),theta(:)
real(mrk), intent(out)::biomass(:),biomass0(:),Fb(:),Ft(:),Fi(:),Fn(:),Rhyd(:)
logical, intent(out)::feas
integer(mik), intent(out)::err
character(*),intent(out)::mess
!locals
character(250),parameter::procname='AlgaeBiomass_Apply'
real(mrk)::dt,b0,bmin,bmax,mu0,Tmin,Topt,Tmax,Iopt,katt,kcw,Rsen,g,rho,J,tau0,bHyd,cHyd
integer(mik)::k,t,n
real(mrk)::bt,Texp,Iratio,tau(size(depth)),irradiancez(size(depth))

! Init
err=0;mess='';feas=.true.
biomass=undefRN;biomass0=undefRN;Fb=undefRN;Ft=undefRN;Fi=undefRN;Fn=undefRN;Rhyd=undefRN

!----------------------------------------
! Check inputs are ok
if(any(irradiance<0._mrk)) then
    mess=trim(procname)//':Irradiance<0 impossible'
    err=1;return
endif
if(any(depth<0._mrk)) then
    mess=trim(procname)//':depth<0 impossible'
    err=1;return
endif
if(any(turbidity<0._mrk)) then
    mess=trim(procname)//':turbidity<0 impossible'
    err=1;return
endif
if(any(removal<0._mrk)) then
    mess=trim(procname)//':removal<0 impossible'
    err=1;return
endif
if(any(removal>1._mrk)) then
    mess=trim(procname)//':removal>1 impossible'
    err=1;return
endif

!----------------------------------------
! make sense of parameter vector
k=0
! generic parameters
k=k+1;dt=theta(k)        ! time step
k=k+1;b0=theta(k)        ! initial biomass at t=0
k=k+1;bmin=theta(k)      ! minimum residual biomass
k=k+1;bmax=theta(k)      ! maximal biomass
! "growth" parameters
k=k+1;mu0=theta(k)       ! maximal growth rate (when all limiting factors are equal to 1)
k=k+1;Tmin=theta(k)      ! minimal temperature
k=k+1;Topt=theta(k)      ! optimal temperature
k=k+1;Tmax=theta(k)      ! maximal temperature
k=k+1;Iopt=theta(k)      ! optimal irradiance
k=k+1;katt=theta(k)      ! irradiance attenuation due to turbidity
k=k+1;kcw=theta(k)       ! irradiance attenuation due to clear water
! "decline" parameters
k=k+1;Rsen=theta(k)      ! senescence rate
k=k+1;g=theta(k)         ! gravity
k=k+1;rho=theta(k)       ! density of water in kg/m3
k=k+1;J=theta(k)         ! slope
k=k+1;tau0=theta(k)      ! critical shear stress in kg/(m*s2) or Pascal
k=k+1;bHyd=theta(k)      ! hydraulic exponent
k=k+1;cHyd=theta(k)      ! hydraulic removal coefficient
! derived parameters
Texp=(Topt-Tmin)/(Tmax-Topt)
!----------------------------------------
! parameter feasability
if(dt<=0._mrk .or. bmin<0._mrk .or. b0<bmin .or. b0>bmax .or. bmin>=bmax .or. &
   mu0<0._mrk .or. Tmin>=Tmax .or. Topt>=Tmax .or. Topt<=Tmin .or. &
   Iopt<=0._mrk .or. katt<0._mrk .or. kcw<0._mrk .or. Rsen<0._mrk .or. g<0._mrk .or. &
   rho<0._mrk .or. J<0._mrk .or. tau0<=0._mrk .or. bHyd<0._mrk .or. cHyd<0._mrk ) then
    feas=.false.;return
endif

!----------------------------------------
! loop through all time steps
n=size(temp)
bt=b0
tau=depth*g*J*rho
do t=1,n
    ! biomass limiting factor
    Fb(t)=1._mrk-bt/bmax
    ! temperature limiting factor
    if(temp(t)<=Tmin .or. temp(t)>=Tmax) then
        Ft(t)=0._mrk
    else
        Ft(t)=((Tmax-temp(t))/(Tmax-Topt)) * ((temp(t)-Tmin)/(Topt-Tmin))**Texp
        if(Ft(t)<=epsRe) Ft(t)=0._mrk ! avoid ridicously small numbers.
    endif
    ! irradiance limiting factor
    irradiancez(t)=irradiance(t)*exp(-(katt*turbidity(t)+kcw)*depth(t))
    Iratio=(irradiancez(t)/Iopt)
    Fi(t)=Iratio*exp(1-Iratio)
    ! nutrient limiting factor
    Fn(t)=1._mrk ! TODO: consider activating this?
    ! hydraulic removal
    if(tau(t)<=tau0) then
        Rhyd(t)=0._mrk
    else
        Rhyd(t)=cHyd * ((tau(t)-tau0)/tau0)**bHyd
    endif
    ! biomass at the end of the time step
    biomass(t)=bt &
               +(mu0*Fb(t)*Ft(t)*Fi(t)*Fn(t))*bt*dt &
               -(Rsen+removal(t)+Rhyd(t))*(bt-bmin)*dt
    ! make sure biomass is between bmin and bmax
    ! TODO: consider using a more elegant way of achieving this...
    biomass(t)=min(biomass(t),bmax)
    biomass(t)=max(biomass(t),bmin)
    ! update bt
    bt=biomass(t)
    biomass0(t)=biomass(t)/bmax
enddo

end subroutine AlgaeBiomass_Apply
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module AlgaeBiomass_model
