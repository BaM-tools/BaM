module DynamicVegetation_model

!~**********************************************************************
!~* Purpose: Dynamic model for rating curves affected by vegetation
!~**********************************************************************
!~* Programmer: Ben Renard & Emeline Perret, Irstea Lyon
!~**********************************************************************
!~* Last modified: 16/03/2020
!~**********************************************************************
!~* Comments:
!~**********************************************************************
!~* References:
!~**********************************************************************
!~* 2Do List:
!~**********************************************************************
!~* Quick description of public procedures:
!~*     1. DynamicVegetation_GetParNumber, number of parameters
!~*     2. DynamicVegetation_Apply, apply model
!~**********************************************************************

use kinds_dmsl_kit ! numeric kind definitions from DMSL

implicit none
Private
public :: DynamicVegetation_GetParNumber, DynamicVegetation_Apply, DynamicVegetation_XtraRead

! Define number of parameters for Q0 and bending submodels.
! NOTE: the slope S0 is NOT in the list of Q0 parameters because it's already
!       part of the parameters for the biomass submodel (used for hydraulic removal).
integer(mik),parameter:: npar_Q0=4,& ! Parameters for Q0: B, n1, b1, c1.
                         npar_bend=2 ! Parameters for bending: chi, U_chi
! Where (i.e. at which location) are gravity g and slope S0 in the list of biomass parameters?
integer(mik),parameter::iGravity=13,iS0=15

Contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine DynamicVegetation_GetParNumber(npar,err,mess)
!^**********************************************************************
!^* Purpose: number of parameters
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified: 16/03/2020
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*     Nothing
!^* OUT
!^*     1. npar, par. number
!^*     2.err, error code; <0:Warning, ==0:OK, >0: Error
!^*     3.mess, error message
!^**********************************************************************
use AlgaeBiomass_model,only: AlgaeBiomass_GetParNumber
integer(mik), intent(out)::npar,err
character(*),intent(out)::mess
! locals
character(250),parameter::procname='DynamicVegetation_GetParNumber'
integer(mik)::npar_biomass

err=0;mess='';npar=undefIN
! Get number of parameters for biomass submodel
call AlgaeBiomass_GetParNumber(npar_biomass,err,mess)
if(err>0) then;mess=trim(procname)//':'//trim(mess);return;endif
! Compute total number of parameters
npar=npar_Q0+npar_bend+npar_biomass
end subroutine DynamicVegetation_GetParNumber
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine DynamicVegetation_Apply(temp,irradiance,stage,turbidity,removal,& ! inputs: water temperature, irradiance at the surface, water stage, turbidity, removal
                                   theta,& ! parameters
                                   NewtonOption,& ! NewtonOption object - see Vegetation_model module.
                                   Q,biomass,biomass0,Fb,Ft,Fi,Fn,Rhyd,& ! outputs
                                   feas,err,mess)
!^**********************************************************************
!^* Purpose: Apply Dynamic Vegetation model, based on algae biomass submodel
!^**********************************************************************
!^* Programmer: Ben Renard & Emeline Perret, Irstea Lyon
!^**********************************************************************
!^* Last modified:
!^* - 16/03/2020
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
! modules needed for submodels
use AlgaeBiomass_model,only:AlgaeBiomass_Apply,AlgaeBiomass_GetParNumber
use Vegetation_model,only:Vegetation_fNewton,Vegetation_NewtonOptionType
! DMSL tools for Newton resolution
use types_dmsl_kit,only:data_ricz_type
use utilities_dmsl_kit, only:cleanPointers
use numerix_dmsl_kit, only:znewton
! Inputs/outputs
real(mrk), intent(in)::temp(:),irradiance(:),stage(:),turbidity(:),removal(:),theta(:)
type(Vegetation_NewtonOptionType), intent(in)::NewtonOption
real(mrk), intent(out)::Q(:),biomass(:),biomass0(:),Fb(:),Ft(:),Fi(:),Fn(:),Rhyd(:)
logical, intent(out)::feas(:)
integer(mik), intent(out)::err
character(*),intent(out)::mess
!locals
character(250),parameter::procname='DynamicVegetation_Apply'
integer(mik)::k,npar_biomass,t,nt,fcalls
real(mrk)::B,S0,n1,b1,c1,chi,U_chi,g,depth(size(stage)),y,Q0,gamma,gamma2,froot
logical::ok
type(data_ricz_type)::NewtonInfo

! Init
err=0;mess='';feas=.true.;Q=undefRN
biomass=undefRN;biomass0=undefRN;Fb=undefRN;Ft=undefRN;Fi=undefRN;Fn=undefRN;Rhyd=undefRN
! Get number of parameters for biomass submodel
call AlgaeBiomass_GetParNumber(npar_biomass,err,mess)
if(err>0) then;mess=trim(procname)//':'//trim(mess);return;endif
!----------------------------------------
! Start by making sense of the parameter vector theta
! The first npar_biomass elements are for the biomass submodel
k=npar_biomass
! Q0 parameters
k=k+1;B=theta(k)
k=k+1;n1=theta(k)
k=k+1;b1=theta(k)
k=k+1;c1=theta(k)
! bending parameters
k=k+1;chi=theta(k)
k=k+1;U_chi=theta(k)
! Get gravity g and slope S0 from list of biomass parameters (also needed for computing Q0 and Q)
S0=theta(iS0)
g=theta(iGravity)
!----------------------------------------
! call biomass submodel
depth=stage-b1
where(depth<=0._mrk) depth=0._mrk ! avoids triggering an error in biomass submodel.
call AlgaeBiomass_Apply(temp,irradiance,depth,turbidity,removal,& ! inputs
                        theta( 1:npar_biomass ),& ! parameters
                        biomass,biomass0,Fb,Ft,Fi,Fn,Rhyd,& ! outputs
                        ok,err,mess)
if(err>0) then;mess=trim(procname)//':'//trim(mess);return;endif
if(.not.ok) then;feas=.false.;return;endif
!----------------------------------------
! Compute discharge
nt=size(temp)
do t=1,nt
    y=depth(t) ! shorthand notation
    ! Handle zero-discharge case
    if(y<=0._mrk) then;Q(t)=0._mrk;cycle;endif
    ! Compute Q0
    Q0=((B*sqrt(S0))/n1)*(y**c1)
    ! Handle case with no biomass
    if(biomass(t)<=0._mrk) then;Q(t)=Q0;cycle;endif
    ! General case with non-zero biomass
    gamma=(biomass(t)*y**(1._mrk/3._mrk))/(8._mrk*g*(n1**2))
    gamma2=(B**chi)*(y**chi)*(U_chi**chi)
    ! Pack parameters Q0, gamma and chi into RICZ structure
    if(associated(NewtonInfo%rp0)) nullify(NewtonInfo%rp0);allocate(NewtonInfo%rp0(4))
    NewtonInfo%rp0(1)=Q0
    NewtonInfo%rp0(2)=gamma
    NewtonInfo%rp0(3)=gamma2
    NewtonInfo%rp0(4)=chi
    ! Apply numerical resolution
    call znewton(evalFunc=Vegetation_fNewton,dataIN=NewtonInfo,x1=0._mrk,x2=Q0,&
                 tolX=NewtonOption%tolX,tolF=NewtonOption%tolF,xscale=NewtonOption%xscale,&
                 fscale=NewtonOption%fscale,itmax=NewtonOption%itmax,&
                 xroot=Q(t),froot=froot,fcalls=fcalls,err=err,message=mess)
    if(err>0) then
        feas(t)=.false.
        call cleanPointers(NewtonInfo,what=1,err=err,message=mess)
        if(err>0) then;mess=trim(procname)//':'//trim(mess);return;endif
        cycle
    endif
    ! Does not find root, cycle
    if(fcalls>NewtonOption%itmax) then
        feas(t)=.false.
        call cleanPointers(NewtonInfo,what=1,err=err,message=mess)
        if(err>0) then;mess=trim(procname)//':'//trim(mess);return;endif
        cycle
    endif
    ! clean RICZ structure to avoid memory leak
    call cleanPointers(NewtonInfo,what=1,err=err,message=mess)
    if(err>0) then;mess=trim(procname)//':'//trim(mess);return;endif
enddo

end subroutine DynamicVegetation_Apply
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine DynamicVegetation_XtraRead(file,xtra,err,mess)
!^**********************************************************************
!^* Purpose: Read Xtra information (properties of Newton algorithm)
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified:16/03/2020
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*     1. file, Xtra file
!^* OUT
!^*     1. xtra, xtra information
!^*     2. err, error code; <0:Warning, ==0:OK, >0: Error
!^*     3. mess, error message
!^**********************************************************************
use types_dmsl_kit, only:data_ricz_type
use utilities_dmsl_kit,only:getSpareUnit
character(*), intent(in)::file
type(data_ricz_type),intent(out):: xtra
integer(mik), intent(out)::err
character(*),intent(out)::mess
! locals
character(250),parameter::procname='DynamicVegetation_XtraRead'
integer(mik), parameter::nNewtonPar=4
integer(mik)::unt

err=0;mess=''
call getSpareUnit(unt,err,mess)
if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
open(unit=unt,file=trim(file), status='old', iostat=err)
if(err>0) then;mess=trim(procname)//':problem opening file '//trim(file);return;endif
if(associated(xtra%rp1)) nullify(xtra%rp1);allocate(xtra%rp1(nNewtonPar))
read(unt,*,iostat=err) xtra%rp1(1) ! Newton - xscale
if(err/=0) then;mess=trim(procname)//':problem reading file '//trim(file);return;endif
read(unt,*,iostat=err) xtra%rp1(2) ! Newton - fscale
if(err/=0) then;mess=trim(procname)//':problem reading file '//trim(file);return;endif
read(unt,*,iostat=err) xtra%rp1(3) ! Newton - xtol
if(err/=0) then;mess=trim(procname)//':problem reading file '//trim(file);return;endif
read(unt,*,iostat=err) xtra%rp1(4) ! Newton - ftol
if(err/=0) then;mess=trim(procname)//':problem reading file '//trim(file);return;endif
read(unt,*,iostat=err) xtra%is2    ! Newton - maxiter
if(err/=0) then;mess=trim(procname)//':problem reading file '//trim(file);return;endif
close(unt)

end subroutine DynamicVegetation_XtraRead
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module DynamicVegetation_model
