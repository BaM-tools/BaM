module ModelLibrary_tools

!~**********************************************************************
!~* Purpose: Catalogue of models
!~**********************************************************************
!~* Programmer: Ben Renard, Irstea Lyon
!~**********************************************************************
!~* Last modified: 27/06/2018
!~**********************************************************************
!~* Comments: All the dirty interfacing with existing models is done
!~*           within this module. In particular models in modules
!~*           XXX_model have no idea of the objects manipulated here
!~*           (e.g. the MODEL object)
!~*
!~*           However, it is expected that any model to be used within
!~*           BAM provides the following three procedures:
!~*           1. XXX_GetParNumber: compute number of parameters of the model
!~*           2. XXX_Apply: apply (run) the model
!~*           3. XXX_XtraRead: load any information needed to run the model
!~**********************************************************************
!~* References:
!~**********************************************************************
!~* 2Do List:
!~**********************************************************************
!~* Quick description of public procedures:
!~*    1. GetModelParNumber, number of parameters of the model
!~*    2. ApplyModel, apply the model
!~*    3. XtraRead, read xtra information
!~*    4. XtraSetup, xtra setup (derived par and states in particular)
!~**********************************************************************

use kinds_dmsl_kit ! numeric kind definitions from DMSL
use types_dmsl_kit, only:data_ricz_type

! Available models
use Linear_model
use Mixture_model
use Orthorectification_model
use SedimentTransport_model
use StageFallDischarge_model
use StageGradientDischarge_model
use BaRatin_model
use BaRatinBAC_model
use GR4J_model
use Vegetation_model
use SWOT_model
use SuspendedLoad_model
use SFDTidal_model
use SFDTidal_Sw_correction_model
use SFDTidal2_model
use SFDTidalJones_model
use SFDTidal4_model
use Recession_model_h
use Segmentation_model
use TextFile_model
use AlgaeBiomass_model
use DynamicVegetation_model
use TidalODE_model
use TidalRemenieras_model
use SMASH_model

implicit none
Private
public :: GetModelParNumber,ApplyModel,XtraRead,XtraSetup,XtraCleanup

! Catalogue of available models
Character(100), parameter, PUBLIC:: &
                    MDL_Linear="MDL_Linear",& ! linear model OUT=IN*theta
                    MDL_Mixture="MDL_Mixture",& ! linear mixture (output is a weighted sum of input contributions)
                    MDL_BaRatin="MDL_BaRatin",& ! BaRatin (Bonnifait formalism), original k-a-c parameterization
                    MDL_BaRatinBAC="MDL_BaRatinBAC",& ! BaRatin (Bonnifait formalism), b-a-c parameterization
                    MDL_Ortho="MDL_Orthorectification",& ! Orthorectification
                    MDL_GR4J="MDL_GR4J",& ! GR4J (4 param + 2 initial state values)
                    MDL_Sediment="MDL_Sediment",& ! Sediment transport
                    MDL_SFD="MDL_SFD",& ! Stage-Fall-Discharge rating curve
                    MDL_SGD="MDL_SGD",& ! Stage-Gradient-Discharge rating curve
                    MDL_Vegetation="MDL_Vegetation",& ! Rating curve affected by vegetation
                    MDL_SWOT="MDL_SWOT",& ! Rating curve from space with SWOT
                    MDL_SuspendedLoad="MDL_SuspendedLoad",& ! Suspended Load
                    MDL_SFDTidal="MDL_SFDTidal",& ! Stage-Fall-Discharge rating curve for tidal rivers
                    MDL_SFDTidal_Sw_correction="MDL_SFDTidal_Sw_correction",& ! Stage-Fall-Discharge rating curve for tidal rivers with linear discharge interpolation during the tidal wave passage
                    MDL_SFDTidal2="MDL_SFDTidal2",& ! Stage-Fall-Discharge rating curve for tidal rivers with water slope correction
                    MDL_SFDTidalJones="MDL_SFDTidalJones",& ! Stage-Fall-Discharge rating curve for tidal rivers with water slope correction
                    MDL_SFDTidal4="MDL_SFDTidal4",& ! Stage-Fall-Discharge rating curve for tidal rivers with water slope correction
                    MDL_Recession_h="MDL_Recession_h", & ! Recession 2-exponential regression for h(t)
                    MDL_Segmentation="MDL_Segmentation",& ! Segmentation of a time series
                    MDL_TextFile="MDL_TextFile",& ! model written in a text file
                    MDL_Algae="MDL_AlgaeBiomass",& ! Biomass dynamics for algae
                    MDL_DynamicVegetation="MDL_DynamicVegetation",& ! Vegetation rating curve model, using Algae Biomass submodel
                    MDL_TidalODE="MDL_TidalODE",& ! Discharge for tidal rivers
                    MDL_TidalRemenieras="MDL_TidalRemenieras",& ! Discharge for tidal rivers
                    MDL_SMASH="MDL_SMASH" ! Interface to SMASH distributed hydrological model

! Model object
type, public:: ModelType ! the "model" object
    character(100)::ID='AintGotNoName' ! ID of the model - see above for the catalogue
    integer(mik)::nX=undefIN ! number of input variables of the model
    integer(mik)::nY=undefIN ! number of output variables of the model
    integer(mik)::ntheta=undefIN ! number of parameters for the model
    character(250),allocatable::parname(:) ! name of the parameters
    type(data_ricz_type)::xtra ! any extra info needed by the model (e.g. for BaRatin: control matrix, for orthorectif: distortion formula, etc.)
    real(mrk)::unfeasFlag=-666.666_mrk ! value denoting unfeasible simulations
    integer(mik)::nDpar=undefIN ! number of "derived" parameters, i.e. parameters that can be written as a function of theta (for instance: offsets in BaRatin)
    character(250),allocatable::DparName(:) ! name of derived parameters
    integer(mik)::nState=undefIN ! number of states, i.e. additional output variables that are not used for estimation but that can be of interest nonetheless (for instance: transition stage between backwter influenced/noninfluenced regimes in SFD; states of an hydrologic model)
    character(250),allocatable::StateName(:) ! name of states
end type ModelType

Contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine GetModelParNumber(model,npar,err,mess)
!^**********************************************************************
!^* Purpose: number of parameters of the model
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified: 10/07/2015
!^**********************************************************************
!^* Comments: Probably useless since the number of parameters is specified
!^*           in the Config_Model file.
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1. model, model object (at least model%ID should be specified)
!^* OUT
!^*    1. npar, par. number
!^*    2.err, error code; <0:Warning, ==0:OK, >0: Error
!^*    3.mess, error message
!^**********************************************************************
type (ModelType), intent(in)::model
integer(mik), intent(out)::npar,err
character(*),intent(out)::mess
! locals
character(250),parameter::procname='GetmodelParNumber'

err=0;mess='';npar=undefIN

select case(model%ID)
case(MDL_Linear)
    call LM_GetParNumber(model%nX,model%nY,npar,err,mess)
case(MDL_Mixture)
    call Mixture_GetParNumber(nS=model%xtra%is1,nEvt=model%xtra%is2,&
                              nElt=model%xtra%is3,missingSource=model%xtra%ls1,&
                              npar=npar,err=err,mess=mess)
case(MDL_BaRatin)
    call BaRatin_GetParNumber(nint(model%xtra%rpm1),npar,err,mess)
case(MDL_BaRatinBAC)
    call BaRatinBAC_GetParNumber(nint(model%xtra%rpm1),npar,err,mess)
case(MDL_GR4J)
    call GR4J_GetParNumber(npar,err,mess)
case(MDL_Ortho)
    call Ortho_GetParNumber(model%xtra%is1,npar,err,mess)
case(MDL_Sediment)
    call Sediment_GetParNumber(ID=model%xtra%cs1,CSS_ID=model%xtra%cs2,&
            PPO=model%xtra%cp1,CO=model%xtra%cs3,npar=npar,err=err,mess=mess)
case(MDL_SFD)
    call SFD_GetParNumber(ID=model%xtra%cs1,npar=npar,err=err,mess=mess)
case(MDL_SGD)
    call SGD_GetParNumber(npar=npar,err=err,mess=mess)
case(MDL_Vegetation)
    call Vegetation_GetParNumber(growthFormula=model%xtra%cs1,nYear=model%xtra%is1,&
                                 npar=npar,err=err,mess=mess)
case(MDL_SWOT)
    call SWOT_GetParNumber(ID=model%xtra%cs1,npar=npar,err=err,mess=mess)
case(MDL_SuspendedLoad)
    call SuspendedLoad_GetParNumber(npar=npar,err=err,mess=mess)
case(MDL_SFDTidal)
    call SFDTidal_GetParNumber(npar=npar,err=err,mess=mess)
case(MDL_SFDTidal_Sw_correction)
    call SFDTidal_Sw_correction_GetParNumber(npar=npar,err=err,mess=mess)
case(MDL_SFDTidal2)
    call SFDTidal2_GetParNumber(npar=npar,err=err,mess=mess)
case(MDL_SFDTidalJones)
    call SFDTidalJones_GetParNumber(npar=npar,err=err,mess=mess)
case(MDL_SFDTidal4)
    call SFDTidal4_GetParNumber(npar=npar,err=err,mess=mess)
case(MDL_Recession_h)
    call Recession_h_GetParNumber(npar=npar,err=err,mess=mess)
case(MDL_Segmentation)
    call Segmentation_GetParNumber(nS=model%xtra%is1,npar=npar,err=err,mess=mess)
case(MDL_TextFile)
    npar=model%ntheta
case(MDL_Algae)
    call AlgaeBiomass_GetParNumber(npar=npar,err=err,mess=mess)
case(MDL_DynamicVegetation)
    call DynamicVegetation_GetParNumber(npar=npar,err=err,mess=mess)
case default
    err=1;mess=trim(procname)//': Fatal: Unavailable [model%ID]';return
case(MDL_TidalODE)
    call TidalODE_GetParNumber(npar=npar,err=err,mess=mess)
case(MDL_TidalRemenieras)
    call TidalRemenieras_GetParNumber(npar=npar,err=err,mess=mess)
case(MDL_SMASH)
    call SMASH_GetParNumber(npar=npar,err=err,mess=mess)
end select
if(err>0) then;mess=trim(procname)//':'//trim(mess);return;endif
end subroutine GetmodelParNumber
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine ApplyModel(model,X,theta,Y,Dpar,state,feas,err,mess)
!^**********************************************************************
!^* Purpose: Apply the model
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified: 10/07/2015
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1. model, model object
!^*    2. X, model inputs
!^*    3. theta, model parameters
!^* OUT
!^*    1. Y, model output
!^*    2. Dpar, derived parameters
!^*    3. state, model states
!^*    4. feas, is computation feasible?
!^*    5.err, error code; <0:Warning, ==0:OK, >0: Error
!^*    6.mess, error message
!^**********************************************************************
type(ModelType), intent(in)::model
real(mrk), intent(in)::X(:,:), theta(:)
real(mrk), intent(out)::Y(:,:),Dpar(:),state(:,:)
logical, intent(out)::feas
integer(mik), intent(out)::err
character(*),intent(out)::mess
!locals
character(250),parameter::procname='ApplyModel'
integer(mik)::i
logical::vfeas(size(X,dim=1))

err=0;mess='';feas=.true.;vfeas=.true.
Y=model%unfeasFlag;Dpar=model%unfeasFlag;state=model%unfeasFlag

select case(model%ID)
case(MDL_Linear)
    call LM_Apply(X=X,theta=theta,Y=Y,err=err,mess=mess)
case(MDL_Mixture)
    call Mixture_Apply(nS=model%xtra%is1,nEvt=model%xtra%is2,&
                       nElt=model%xtra%is3,missingSource=model%xtra%ls1,&
                       X=X,theta=theta,Y=Y(:,1),theta_last=Dpar,&
                       feas=feas,err=err,mess=mess)
case(MDL_BaRatin)
    call BaRatin_Apply(H=X(:,1),theta=theta,ControlMatrix=nint(model%xtra%rpm1),&
                       Q=Y(:,1),b=Dpar,feas=vfeas,err=err,mess=mess)
case(MDL_BaRatinBAC)
    i=size(model%xtra%rpm1,dim=1)
    call BaRatinBAC_Apply(H=X(:,1),theta=theta,ControlMatrix=nint(model%xtra%rpm1),&
                       hmax=model%xtra%rs1,Q=Y(:,1),k=Dpar(1:i),ktype=Dpar((i+1):(2*i)),&
                       feas=vfeas,err=err,mess=mess)
case(MDL_Ortho)
    call Ortho_Apply(Disto_ID=model%xtra%cs1,Disto_npar=model%xtra%is1,IN=X,&
                     theta=theta,OUT=Y,Dpar=Dpar,feas=vfeas,err=err,mess=mess)
case(MDL_Sediment)
    call Sediment_Apply(ID=model%xtra%cs1,CSS_ID=model%xtra%cs2,&
                        PPO=model%xtra%cp1,CO=model%xtra%cs3,&
                        IN=X,theta=theta,pseudoval=model%xtra%rp1,&
                        OUT=Y(:,1),feas=vfeas,err=err,mess=mess)
case(MDL_SFD)
    call SFD_Apply(ID=model%xtra%cs1,IN=X,theta=theta,&
                   NewtonOption=SFD_NewtonOptionType(upperbound=model%xtra%rp1(1),&
                                             xscale=model%xtra%rp1(2),&
                                             fscale=model%xtra%rp1(3),&
                                             tolX=model%xtra%rp1(4),&
                                             tolF=model%xtra%rp1(5),&
                                             itmax=model%xtra%is1),&
                   OUT=Y(:,1),kappa=state(:,1),feas=vfeas,err=err,mess=mess)
case(MDL_SGD)
    call SGD_Apply(IN=X,theta=theta,OUT=Y(:,1),feas=vfeas,err=err,mess=mess)
case(MDL_GR4J)
    call GR4J_Apply(X=X,theta=theta,pcofile=model%xtra%cs1,Y=Y,state=state,feas=feas,err=err,mess=mess)
case(MDL_Vegetation)
    call Vegetation_Apply(growthFormula=model%xtra%cs1,nYear=model%xtra%is1,&
                          IN=X,theta=theta,&
                          NewtonOption=Vegetation_NewtonOptionType(&
                                             xscale=model%xtra%rp1(1),&
                                             fscale=model%xtra%rp1(2),&
                                             tolX=model%xtra%rp1(3),&
                                             tolF=model%xtra%rp1(4),&
                                             itmax=model%xtra%is2),&
                          OUT=Y,Dpar=Dpar,feas=vfeas,err=err,mess=mess)
case(MDL_SWOT)
    call SWOT_Apply(ID=model%xtra%cs1,IN=X,theta=theta,OUT=Y(:,1),&
                    feas=vfeas,err=err,mess=mess)
case(MDL_SuspendedLoad)
    call SuspendedLoad_Apply(h=X(:,1),theta=theta,Qs=Y(:,1),&
                             feas=vfeas,err=err,mess=mess)
case(MDL_SFDTidal)
    call SFDTidal_Apply(IN=X,theta=theta,OUT=Y(:,1),&
                        ComputationOption=model%xtra%cs1,feas=vfeas,err=err,mess=mess)
case(MDL_SFDTidal_Sw_correction)
    call SFDTidal_Sw_correction_Apply(IN=X,theta=theta,Q=Y(:,1),interp=state(:,1),&
                        feas=vfeas,err=err,mess=mess)
case(MDL_SFDTidal2)
    call SFDTidal2_Apply(IN=X,theta=theta,OUT=Y(:,1),&
                        feas=vfeas,err=err,mess=mess)
case(MDL_SFDTidalJones)
    call SFDTidalJones_Apply(IN=X,theta=theta,OUT=Y(:,1),&
                        feas=vfeas,err=err,mess=mess)
case(MDL_SFDTidal4)
    call SFDTidal4_Apply(IN=X,theta=theta,OUT=Y(:,1),&
                        feas=vfeas,err=err,mess=mess)
case(MDL_Recession_h)
    call Recession_h_Apply(time=X(:,1), theta=theta, h=Y(:,1),feas=feas,err=err,mess=mess)

case(MDL_Segmentation)
    call Segmentation_Apply(time=X(:,1),nS=model%xtra%is1, tmin=model%xtra%rs1, theta=theta, y= Y(:,1),feas=feas,err=err,mess=mess)
case(MDL_TextFile)
    call TxtMdl_Apply(IN=X,theta=theta,OUT=Y,feas=vfeas,err=err,mess=mess)
case(MDL_Algae)
    call AlgaeBiomass_Apply(temp=X(:,1),irradiance=X(:,2),depth=X(:,3),&
                            turbidity=X(:,4),removal=X(:,5),&
                            theta=theta,&
                            biomass=Y(:,1),biomass0=Y(:,2),&
                            Fb=state(:,1),Ft=state(:,2),Fi=state(:,3),&
                            Fn=state(:,4),Rhyd=state(:,5),&
                            feas=feas,err=err,mess=mess)
case(MDL_DynamicVegetation)
    call DynamicVegetation_Apply(temp=X(:,1),irradiance=X(:,2),stage=X(:,3),&
                                 turbidity=X(:,4),removal=X(:,5),&
                                 theta=theta,&
                                 NewtonOption=Vegetation_NewtonOptionType(&
                                             xscale=model%xtra%rp1(1),&
                                             fscale=model%xtra%rp1(2),&
                                             tolX=model%xtra%rp1(3),&
                                             tolF=model%xtra%rp1(4),&
                                             itmax=model%xtra%is2),&
                                 Q=Y(:,1),biomass=Y(:,2),biomass0=Y(:,3),&
                                 Fb=state(:,1),Ft=state(:,2),Fi=state(:,3),&
                                 Fn=state(:,4),Rhyd=state(:,5),&
                                 feas=vfeas,err=err,mess=mess)
case(MDL_TidalODE)
    call TidalODE_Apply(IN=X,theta=theta,OUT=Y(:,1),&
                        ComputationOption=model%xtra%cs1,feas=vfeas,err=err,mess=mess)!
case(MDL_TidalRemenieras)
    call TidalRemenieras_Apply(IN=X,theta=theta,OUT=Y(:,1),&
                        ComputationOption=model%xtra%cs1,feas=vfeas,err=err,mess=mess)!
case(MDL_SMASH)
    call SMASH_Run(projectDir=model%xtra%cp1(2),communicationDir=model%xtra%cp1(5),&
                   QSIMfile=model%xtra%cp1(6),theta=theta,Y=Y,feas=feas,err=err,mess=mess)
case default
    err=1;mess=trim(procname)//': Fatal: Unavailable [model%ID]'
end select

! handle errors & unfeasible runs
if(err>0) then;mess=trim(procname)//':'//trim(mess);return;endif
if(.not.feas) vfeas=.false.
feas=all(vfeas)
do i=1,model%nY
    where(.not.vfeas) Y(:,i)=model%unfeasFlag
enddo
do i=1,model%nState
    where(.not.vfeas) state(:,i)=model%unfeasFlag
enddo

end subroutine ApplyModel
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine XtraRead(file,ID,xtra,err,mess)
!^**********************************************************************
!^* Purpose: Read Xtra config files for a model
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified: 21/08/2017
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1. file, file where xtra information should be read
!^*    2. ID, model ID
!^* OUT
!^*    1.xtra, extra model info
!^*    2.err, error code; <0:Warning, ==0:OK, >0: Error
!^*    3.mess, error message
!^**********************************************************************
use types_dmsl_kit, only:data_ricz_type
character(*), intent(in)::file,ID
type(data_ricz_type),intent(out):: xtra
integer(mik), intent(out)::err
character(*),intent(out)::mess
! locals
character(250),parameter::procname='XtraRead'

err=0;mess=''

select case (trim(ID))
case(MDL_Linear)
    ! nothing to read
case(MDL_Mixture)
    call Mixture_XtraRead(file=file,xtra=xtra,err=err,mess=mess)
case(MDL_BaRatin)
    call BaRatin_XtraRead(file=file,xtra=xtra,err=err,mess=mess)
case(MDL_BaRatinBAC)
    call BaRatinBAC_XtraRead(file=file,xtra=xtra,err=err,mess=mess)
case(MDL_GR4J)
    call GR4J_XtraRead(file=file,xtra=xtra,err=err,mess=mess)
case(MDL_Ortho)
    call Ortho_XtraRead(file=file,xtra=xtra,err=err,mess=mess)
case(MDL_Sediment)
    call Sediment_XtraRead(file=file,xtra=xtra,err=err,mess=mess)
case(MDL_SFD)
    call SFD_XtraRead(file=file,xtra=xtra,err=err,mess=mess)
case(MDL_SGD)
    ! nothing to read
case(MDL_Vegetation)
    call Vegetation_XtraRead(file=file,xtra=xtra,err=err,mess=mess)
case(MDL_SWOT)
    call SWOT_XtraRead(file=file,xtra=xtra,err=err,mess=mess)
case(MDL_SuspendedLoad)
    ! nothing to read
case(MDL_SFDTidal)
    ! nothing to read
    !call SFDTidal_XtraRead(file=file,xtra=xtra,err=err,mess=mess)
case(MDL_SFDTidal_Sw_correction)
    ! nothing to read
case(MDL_SFDTidal2)
    ! nothing to read
case(MDL_SFDTidalJones)
    ! nothing to read
case(MDL_SFDTidal4)
    ! nothing to read
case(MDL_Recession_h)
    ! nothing to read
case(MDL_Segmentation)
    call Segmentation_XtraRead(file=file,xtra=xtra,err=err,mess=mess)
case(MDL_TextFile)
    call TxtMdl_load(file=file,err=err,mess=mess)
case(MDL_Algae)
    ! nothing to read
case(MDL_DynamicVegetation)
    call DynamicVegetation_XtraRead(file=file,xtra=xtra,err=err,mess=mess)
case(MDL_TidalODE)
    ! nothing to read
    !call TidalODE_XtraRead(file=file,xtra=xtra,err=err,mess=mess)
case(MDL_TidalRemenieras)
    ! nothing to read
    !call TidalRemenieras_XtraRead(file=file,xtra=xtra,err=err,mess=mess)
case(MDL_SMASH)
    call SMASH_XtraRead(file=file,xtra=xtra,err=err,mess=mess)
case default
    err=1;mess=trim(procname)//': Fatal: Unavailable [model%ID]'
end select

if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif

end subroutine XtraRead
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine XtraSetup(model,err,mess)
!^**********************************************************************
!^* Purpose: Xtra model setup (derived pars and states in particular)
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified: 22/08/2017
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* INOUT
!^*    1. model, model object
!^* OUT
!^*    1.err, error code; <0:Warning, ==0:OK, >0: Error
!^*    2.mess, error message
!^**********************************************************************
use utilities_dmsl_kit,only:number_string
type(ModelType), intent(inout)::model
integer(mik), intent(out)::err
character(*),intent(out)::mess
! locals
character(250),parameter::procname='XtraSetup'
integer(mik)::i

err=0;mess=''

select case (trim(model%ID))
case(MDL_Linear)
    model%nDpar=0;model%nState=0
    allocate(model%DparName(model%nDpar));allocate(model%StateName(model%nState))
case(MDL_Mixture)
    model%nDpar=model%xtra%is2;model%nState=0
    allocate(model%DparName(model%nDpar));allocate(model%StateName(model%nState))
    do i=1,model%nDpar
        model%DparName(i)='LastWeight'//trim(number_string(i))
    enddo
case(MDL_BaRatin)
    model%nDpar=model%ntheta/3 ! offsets
    model%nState=0
    allocate(model%DparName(model%nDpar));allocate(model%StateName(model%nState))
    do i=1,model%nDpar
        model%DparName(i)='b'//trim(number_string(i))
    enddo
case(MDL_BaRatinBAC)
    model%nDpar=(model%ntheta/3)*2 ! activation stages and their type
    model%nState=0
    allocate(model%DparName(model%nDpar));allocate(model%StateName(model%nState))
    do i=1,(model%ntheta/3)
        model%DparName(i)='k'//trim(number_string(i))
    enddo
    do i=1,(model%ntheta/3)
        model%DparName(model%ntheta/3+i)='ktype'//trim(number_string(i))
    enddo
case(MDL_GR4J)
    model%nDpar=0
    model%nState=2 ! production and routing stores
    allocate(model%DparName(model%nDpar));allocate(model%StateName(model%nState))
    model%StateName=(/'production','routing   '/)
case(MDL_Ortho)
    model%nDpar=20 ! rotation matrix + ortho-pars
    model%nState=0
    allocate(model%DparName(model%nDpar));allocate(model%StateName(model%nState))
    model%DparName=(/'r11','r12','r13','r21','r22','r23','r31','r32','r33',&
                     'a1 ','a2 ','a3 ','a4 ','b1 ','b2 ','b3 ','b4 ','c1 ','c2 ','c3 '/)
case(MDL_Sediment)
    model%nDpar=0;model%nState=0
    allocate(model%DparName(model%nDpar));allocate(model%StateName(model%nState))
case(MDL_SFD)
    model%nDpar=0
    model%nState=1 ! transition stage between backwter influenced/noninfluenced regimes
    allocate(model%DparName(model%nDpar));allocate(model%StateName(model%nState))
    model%StateName=(/'TransitionStage'/)
case(MDL_SGD)
    model%nDpar=0;model%nState=0
    allocate(model%DparName(model%nDpar));allocate(model%StateName(model%nState))
case(MDL_Vegetation)
    model%nState=0
    if(model%xtra%cs1==Growth_Yin1_temporal) then
        model%nDpar=1
    else
        model%nDpar=model%xtra%is1+1
    endif
    allocate(model%DparName(model%nDpar));allocate(model%StateName(model%nState))
    model%DparName(1)="a1"
    if(model%nDpar>1) then
        do i=2,model%nDpar
            model%DparName(i)='Topt_'//trim(number_string(i))
        enddo
    endif
case(MDL_SWOT)
    model%nDpar=0
    model%nState=0
    allocate(model%DparName(model%nDpar));allocate(model%StateName(model%nState))
case(MDL_SuspendedLoad)
    model%nDpar=0
    model%nState=0
    allocate(model%DparName(model%nDpar));allocate(model%StateName(model%nState))
case(MDL_SFDTidal)
    model%nDpar=0
    model%nState=0
    allocate(model%DparName(model%nDpar));allocate(model%StateName(model%nState))
case(MDL_SFDTidal_Sw_correction)
    model%nDpar=0
    model%nState=1
    allocate(model%DparName(model%nDpar));allocate(model%StateName(model%nState))
    model%StateName=(/'interp'/)
case(MDL_SFDTidal2)
    model%nDpar=0
    model%nState=0
    allocate(model%DparName(model%nDpar));allocate(model%StateName(model%nState))
case(MDL_SFDTidalJones)
    model%nDpar=0
    model%nState=0
    allocate(model%DparName(model%nDpar));allocate(model%StateName(model%nState))
case(MDL_SFDTidal4)
    model%nDpar=0
    model%nState=0
    allocate(model%DparName(model%nDpar));allocate(model%StateName(model%nState))
case(MDL_Recession_h)
    model%nDpar=0
    model%nState=0
    allocate(model%DparName(model%nDpar));allocate(model%StateName(model%nState))
case(MDL_TextFile)
    model%nDpar=0
    model%nState=0
    allocate(model%DparName(model%nDpar));allocate(model%StateName(model%nState))
case(MDL_Algae)
    model%nDpar=0
    model%nState=5
    allocate(model%DparName(model%nDpar));allocate(model%StateName(model%nState))
    model%StateName=(/'Fb   ','Ft   ','Fi   ','Fn   ','Rhyd '/)
case(MDL_DynamicVegetation)
    model%nDpar=0
    model%nState=5
    allocate(model%DparName(model%nDpar));allocate(model%StateName(model%nState))
    model%StateName=(/'Fb   ','Ft   ','Fi   ','Fn   ','Rhyd '/)
case(MDL_TidalODE)
    model%nDpar=0
    model%nState=0
    allocate(model%DparName(model%nDpar));allocate(model%StateName(model%nState))
case(MDL_TidalRemenieras)
    model%nDpar=0
    model%nState=0
    allocate(model%DparName(model%nDpar));allocate(model%StateName(model%nState))
case(MDL_SMASH)
    model%nDpar=0
    model%nState=0
    allocate(model%DparName(model%nDpar));allocate(model%StateName(model%nState))
    ! Load SMASH
    call SMASH_Load(loadScript=model%xtra%cp1(1),projectDir=model%xtra%cp1(2),&
                    precipDir=model%xtra%cp1(3),petDir=model%xtra%cp1(4),&
                    communicationDir=model%xtra%cp1(5),err=err,mess=mess)
    if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
case default
    err=1;mess=trim(procname)//': Fatal: Unavailable [model%ID]'
end select

end subroutine XtraSetup
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine XtraCleanup(model,err,mess)
!^**********************************************************************
!^* Purpose: Xtra model cleanup. Subroutine used for any cleanup action
!^*          that may be necessary after model use.
!^*          For many models this sub does nothing.
!^**********************************************************************
!^* Programmer: Ben Renard, INRAE Aix
!^**********************************************************************
!^* Last modified: 05/08/2022
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* INOUT
!^*    1. model, model object
!^* OUT
!^*    1.err, error code; <0:Warning, ==0:OK, >0: Error
!^*    2.mess, error message
!^**********************************************************************
use utilities_dmsl_kit,only:number_string
type(ModelType), intent(inout)::model
integer(mik), intent(out)::err
character(*),intent(out)::mess
! locals
character(250),parameter::procname='XtraCleanup'

err=0;mess=''

if(trim(model%ID)==MDL_SMASH) then ! Stops Python infinite loop
    call SMASH_Cleanup(projectDir=model%xtra%cp1(2),communicationDir=model%xtra%cp1(5),err=err,mess=mess)
    if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
endif

end subroutine XtraCleanup
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!==============!
! Private subs !
!==============!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module ModelLibrary_tools
