module BaM_tools

!~**********************************************************************
!~* Purpose: Bayesian Model fitting
!~**********************************************************************
!~* Programmer: Ben Renard, Irstea Lyon
!~**********************************************************************
!~* Last modified:10/07/2015
!~**********************************************************************
!~* Comments: A generalisation of BaRatin, for a model with possibly
!~*           several inputs/outputs
!~*           or alternatively, a simplified version of BATEAU_DK
!~**********************************************************************
!~* References:
!~**********************************************************************
!~* 2Do List:
!~**********************************************************************
!~* Quick description of public procedures:
!~*    1.
!~*    2.
!~*    3.
!~**********************************************************************

use kinds_dmsl_kit ! numeri c kind definitions from DMSL
use BayesianEstimation_tools, only:PriorListType
use ModelLibrary_tools, only:ModelType

implicit none
Private

public :: &! main subroutines to handle the probabilistic model behind BaM
          LoadBamObjects,BaM_Fit,BaM_ReadData,BaM_CookMCMC,BaM_SummarizeMCMC,BaM_computeDIC,&
          ! Read config files
          Config_Read,Config_Read_Model,Config_Read_Xtra,Config_Read_Data,&
          Config_Read_RemnantSigma,Config_Read_MCMC,Config_Read_RunOptions,&
          Config_Read_Residual,Config_Read_Cook,Config_Read_Summary,&
          Config_Read_Pred_Master,Config_Read_Pred,&
          Config_Finalize,&
          BaM_ConsoleMessage,BaM_PrintHelp,&
          ! Post-processing tools
          BaM_LoadMCMC,BaM_Residual,BaM_Prediction,BaM_ReadSpag,BaM_Cleanup

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! variables globally available to this module

! prior list
type,public::plist
    type(PriorListType), allocatable:: p(:) ! priors for remnant std - see BayesianEstimation_tools module for PriorListType definition
end type plist

! par list
type,public::parlist
    real(mrk), allocatable:: p(:) ! list of parameter values
end type parlist

! string list
type,public::slist
    character(250), allocatable:: s(:) ! list of strings
end type slist

! parameter (theta)
type,public::ParType
    character(250):: parname='' ! parameter name
    integer(mik):: parType=undefIN ! parameter type: fixed, standard, varying or stochastic
    integer(mik):: nval=undefIN ! number of distinct values (1 for STD and FIX)
    integer(mik),allocatable:: indx(:) ! indices associated with each observation (size nObs, should be between 1 and nval)
    type(PriorListType), allocatable:: prior(:) ! priors (size nval)
    real(mrk),allocatable:: init(:) ! initial values (size nval)
    character(250):: config='' ! config file (only for VAR and STOK parameters)
end type ParType

! Inference information
type:: InferenceType
    integer(mik)::nobs=undefIN ! number of observations used to fit the model
    integer(mik),allocatable::nremnant(:) ! number of parameters for remnant errors (size nY)
    logical,allocatable:: Xerror(:) ! is any input value affected by an error? if not, this is WLS-like fitting (size nX)
    integer(mik),allocatable::nXerror(:) ! number of input values to be estimated (size nX)
    integer(mik)::nInfer=undefIN ! total number of inferred parameters
    integer(mik)::nFit=undefIN ! Number of fitted parameters theta (sum of nvals)
    integer(mik)::nFix=undefIN ! Number of fixed parameters theta
    integer(mik)::nStd=undefIN ! Number of standard parameters theta
    integer(mik)::nVar=undefIN ! Number of varying parameters theta
    integer(mik)::nStok=undefIN ! Number of stochastic parameters theta
    type(ParType), allocatable:: theta(:) ! properties of theta parameters (size nTheta)
    integer(mik),allocatable::L(:,:) ! localisation matrix (size nobs * nTheta): L(t,i)=index where to find ith theta value at time t in a flatten vector theta.
    real(mrk), allocatable::X(:,:),Y(:,:) ! observed series of inputs/outputs (nobs*nX, nobs*nY)
    real(mrk), allocatable::Xu(:,:),Yu(:,:) ! series of stdevs for random observation errors (nobs*nX, nobs*nY)
    real(mrk), allocatable::Xb(:,:),Yb(:,:) ! series of stdevs for systematic observation errors (nobs*nX, nobs*nY)
    integer(mik), allocatable::Xbindx(:,:),Ybindx(:,:) ! series of indices for systematic observation errors (nobs*nX, nobs*nY)
    integer(mik),allocatable::nXb(:),nYb(:) ! number of input/output unknown biases (size nX/nY)
    type(slist), allocatable:: Parname_RemnantSigma(:) ! names for remnant std parameters (size nY)
    type(plist), allocatable:: Prior_RemnantSigma(:) ! priors for remnant std (size nY)
    type(PriorListType), allocatable:: Prior_Theta(:) ! priors for theta, flatten as a vector (size nFit)
    character(250),allocatable::RemnantSigma_funk(:) ! function for std of remnant errors (size nY)
    integer(mik)::nDpar ! number of derived parameters
    character(250),allocatable::DparName(:) ! names of derived parameters
    integer(mik),allocatable::DparIndx(:) ! indices corresponding to combinations of VAR par values, used to deduce Derived parameters
    real(mrk)::priorLogDet=undefRN ! log-determinant of the prior correlation matrix C
    real(mrk),allocatable::priorM(:,:) ! matrix (Cinv - I) used to compute the joint prior via a Gaussian copula
    real(mrk)::mv=-9999._mrk
end type InferenceType

! "Vector of Matrixes" - used for XspagType, see below
type::VoMtype
    real(mrk), allocatable:: m(:,:) ! spaghettis, Nobs*nspag
end type VoMtype

! Spaghettis of input variables - used for uncertainty propagation
type,public::XspagType
    integer(mik)::nobs=undefIN
    integer(mik),allocatable::nspag(:) ! size nX
    type(VoMtype), allocatable:: spag(:) ! size nX. The kth matrix is Nobs*nspag(k)
end type XspagType

integer(mik),parameter,public::messID_Read=-3,messID_Open=-2,messID_Write=-8
character(250),parameter::fmt_numeric='e15.6',fmt_string='A15'
integer(mik),parameter::FIX=-1,STD=0,VAR=1,STOK=2 ! flags for fixed, standard, varying and stochastic parameters
character(250),parameter::FIX_str='FIX',STD_str='STD',VAR_str='VAR',STOK_str='STOK'
real(mrk),parameter::bias_init=0._mrk,bias_initStd=0.1_mrk
type(InferenceType):: INFER
type(ModelType):: MODEL

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!=============================!
! Main BaM public subroutines !
!=============================!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine LoadBamObjects(X,Xu,Xb,Xbindx,& ! observed inputs and their uncertainties
                          Y,Yu,Yb,Ybindx,& ! observed outputs and their uncertainties
                          ID,&             ! Model ID
                          theta,&          ! theta parameter object
                          RemnantSigma_funk,& ! chosen f in { residual var = f(Qrc) }
                          Parname_RemnantSigma,& ! parameter names
                          Prior_RemnantSigma,& ! priors
                          priorCorrFile,& ! name of file containing prior correlations
                          infoFile,& ! name of file containing BaM information
                          xtra,& ! Xtra information
                          nstate,err,mess)! error handling
!^**********************************************************************
!^* Purpose: Load model object
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified:06/08/2018
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List: Additional size checks
!^**********************************************************************
!^* IN
!^*    1. X, observed inputs
!^*    2. Xu, standard deviations of random errors
!^*    3. Xb, standard deviations of systematic errors
!^*    4. Xbindx, index of systematic errors
!^*    5. Y, observed outputs
!^*    6. Yu, standard deviations of output errors
!^*    7. Yb, standard deviations of systematic errors
!^*    8. Ybindx, index of systematic errors
!^*    9. ID, ID of the model
!^*    10. theta, theta parameter object
!^*    11. RemnantSigma_funk, function for remnant sigma
!^*    12. Parname_RemnantSigma, parameter names for remnant sigma parameters
!^*    13. Prior_RemnantSigma, priors for remnant sigma parameters
!^*    14. priorCorrFile, name of file containing prior correlations
!^*    15. [xtra], optional xtra information needed by the model
!^* OUT
!^*    1.nstate, number of state variables
!^*    2.err, error code; <0:Warning, ==0:OK, >0: Error
!^*    3.mess, error message
!^**********************************************************************
use types_dmsl_kit, only:data_ricz_type
use ModelLibrary_tools, only:XtraSetup
use utilities_dmsl_kit,only:number_string,getSpareUnit
use linalg_dmsl_kit,only:choles_invrt

real(mrk), intent(in)::X(:,:),Xu(:,:),Xb(:,:)
real(mrk), intent(in)::Y(:,:),Yu(:,:),Yb(:,:)
integer(mik), intent(in)::Xbindx(:,:),Ybindx(:,:)
character(*), intent(in)::ID
type(ParType), intent(in):: theta(:)
type(slist), intent(in):: Parname_RemnantSigma(:)
type(plist), intent(in):: Prior_RemnantSigma(:)
character(*),intent(in)::RemnantSigma_funk(:),priorCorrFile,infoFile
type(data_ricz_type),intent(in),optional:: xtra
integer(mik), intent(out)::err,nstate
character(*),intent(out)::mess
! locals
character(250),parameter::procname='LoadBamObjects'
integer(mik)::i,j,k
integer(mik)::dummyIndx(size(X,dim=1)),dummyCombo(size(X,dim=1),size(theta))
logical::match,feas

err=0;mess='';nstate=undefIN
INFER%nobs=size(X,dim=1)
MODEL%nX=size(X,dim=2)
MODEL%nY=size(Y,dim=2)
! Size checks
if( size(Xu,dim=1)/=INFER%nobs .or. size(Yu,dim=1)/=INFER%nobs .or. size(Y,dim=1)/=INFER%nobs ) then
    err=1;mess=trim(procname)//':'//trim(BaM_Message(1));return
endif
MODEL%ID=ID
!----------------------
! Data
!----------------------
! Input/output data
if(allocated(INFER%X)) deallocate(INFER%X);allocate(INFER%X(INFER%nobs,MODEL%nX))
if(allocated(INFER%Y)) deallocate(INFER%Y);allocate(INFER%Y(INFER%nobs,MODEL%nY))
INFER%X=X
INFER%Y=Y
! Input/output uncertainties
if(allocated(INFER%Xu)) deallocate(INFER%Xu);allocate(INFER%Xu(INFER%nobs,MODEL%nX))
if(allocated(INFER%Yu)) deallocate(INFER%Yu);allocate(INFER%Yu(INFER%nobs,MODEL%nY))
INFER%Xu=Xu
INFER%Yu=Yu
if(allocated(INFER%Xb)) deallocate(INFER%Xb);allocate(INFER%Xb(INFER%nobs,MODEL%nX))
if(allocated(INFER%Yb)) deallocate(INFER%Yb);allocate(INFER%Yb(INFER%nobs,MODEL%nY))
INFER%Xb=Xb
INFER%Yb=Yb
if(allocated(INFER%Xbindx)) deallocate(INFER%Xbindx);allocate(INFER%Xbindx(INFER%nobs,MODEL%nX))
if(allocated(INFER%Ybindx)) deallocate(INFER%Ybindx);allocate(INFER%Ybindx(INFER%nobs,MODEL%nY))
INFER%Xbindx=Xbindx
INFER%Ybindx=Ybindx
! Determine whether some true inputs need to be estimated
if(allocated(INFER%Xerror)) deallocate(INFER%Xerror);allocate(INFER%Xerror(MODEL%nX))
if(allocated(INFER%nXerror)) deallocate(INFER%nXerror);allocate(INFER%nXerror(MODEL%nX))
do j=1,MODEL%nX
    INFER%Xerror(j)=any(Xu(:,j)>0._mrk)
    INFER%nXerror(j)=count(Xu(:,j)>0._mrk)
enddo
! Determine number of unknown biases
if(allocated(INFER%nXb)) deallocate(INFER%nXb);allocate(INFER%nXb(MODEL%nX))
if(allocated(INFER%nYb)) deallocate(INFER%nYb);allocate(INFER%nYb(MODEL%nY))
do j=1,MODEL%nX
    if(size(INFER%Xbindx(:,j))==0) then ! may happen for 0-data calibration
        INFER%nXb(j)=0
    else
        INFER%nXb(j)=maxval(INFER%Xbindx(:,j))
    endif
enddo
do j=1,MODEL%nY
   if(size(INFER%Ybindx(:,j))==0) then ! may happen for 0-data calibration
        INFER%nYb(j)=0
    else
        INFER%nYb(j)=maxval(INFER%Ybindx(:,j))
    endif
enddo
!----------------------
! Model
!----------------------
MODEL%ntheta=size(theta)
if(allocated(MODEL%parname)) deallocate(MODEL%parname);allocate(MODEL%parname(MODEL%ntheta))
MODEL%parname=theta(:)%parname
if(allocated(INFER%nremnant)) deallocate(INFER%nremnant);allocate(INFER%nremnant(MODEL%nY))
if(allocated(INFER%Parname_RemnantSigma)) deallocate(INFER%Parname_RemnantSigma)
allocate(INFER%Parname_RemnantSigma(MODEL%nY))
do j=1,MODEL%nY
    INFER%nremnant(j)=size(Prior_RemnantSigma(j)%p)
    if(allocated(INFER%Parname_RemnantSigma(j)%s)) deallocate(INFER%Parname_RemnantSigma(j)%s)
    allocate(INFER%Parname_RemnantSigma(j)%s(INFER%nremnant(j)))
    INFER%Parname_RemnantSigma(j)%s=Parname_RemnantSigma(j)%s
enddo
! Xtra info
if(present(xtra)) MODEL%xtra=xtra
! handle RemnantSigma
if(allocated(INFER%Prior_RemnantSigma)) deallocate(INFER%Prior_RemnantSigma)
allocate(INFER%Prior_RemnantSigma(MODEL%nY))
do j=1,MODEL%nY
    if(allocated(INFER%Prior_RemnantSigma(j)%p)) deallocate(INFER%Prior_RemnantSigma(j)%p)
    allocate(INFER%Prior_RemnantSigma(j)%p(INFER%nremnant(j)))
    INFER%Prior_RemnantSigma(j)%p=Prior_RemnantSigma(j)%p
enddo
if(allocated(INFER%RemnantSigma_funk)) deallocate(INFER%RemnantSigma_funk)
allocate(INFER%RemnantSigma_funk(MODEL%NY))
INFER%RemnantSigma_funk=RemnantSigma_funk
!----------------------
! Theta
!----------------------
if(allocated(INFER%theta)) deallocate(INFER%theta);allocate(INFER%theta(MODEL%ntheta))
INFER%nFit=0
do j=1,MODEL%nTheta
    INFER%theta(j)%parname=theta(j)%parname
    INFER%theta(j)%parType=theta(j)%parType
    INFER%theta(j)%nVal=theta(j)%nVal
    if(INFER%theta(j)%parType/=FIX) INFER%nFit=INFER%nFit+INFER%theta(j)%nVal
    if(allocated(INFER%theta(j)%indx)) deallocate(INFER%theta(j)%indx)
    allocate(INFER%theta(j)%indx(INFER%nobs))
    if(allocated(INFER%theta(j)%prior)) deallocate(INFER%theta(j)%prior)
    allocate(INFER%theta(j)%prior(INFER%theta(j)%nVal))
    if(allocated(INFER%theta(j)%init)) deallocate(INFER%theta(j)%init)
    allocate(INFER%theta(j)%init(INFER%theta(j)%nVal))
    INFER%theta(j)%indx=theta(j)%indx
    INFER%theta(j)%init=theta(j)%init
    do k=1,INFER%theta(j)%nVal
        INFER%theta(j)%prior(k)%dist=theta(j)%prior(k)%dist
        if(allocated(INFER%theta(j)%prior(k)%par)) deallocate(INFER%theta(j)%prior(k)%par)
        allocate(INFER%theta(j)%prior(k)%par(size(theta(j)%prior(k)%par)))
        INFER%theta(j)%prior(k)%par=theta(j)%prior(k)%par
    enddo
enddo
INFER%nFix=count(INFER%theta%partype==FIX)
INFER%nStd=count(INFER%theta%partype==STD)
INFER%nVar=count(INFER%theta%partype==VAR)
INFER%nStok=count(INFER%theta%partype==STOK)
! Flatten priors
if(allocated(INFER%Prior_Theta)) deallocate(INFER%Prior_Theta);allocate(INFER%Prior_Theta(INFER%nFit))
k=0
do j=1,MODEL%nTheta
    if(INFER%theta(j)%parType/=FIX) then
        do i=1,INFER%theta(j)%nVal
            k=k+1
            INFER%Prior_Theta(k)%dist=INFER%theta(j)%prior(i)%dist
            if(allocated(INFER%Prior_Theta(k)%par)) deallocate(INFER%Prior_Theta(k)%par)
            allocate(INFER%Prior_Theta(k)%par(size(INFER%theta(j)%prior(i)%par)))
            INFER%Prior_Theta(k)%par=INFER%theta(j)%prior(i)%par
        enddo
    endif
enddo
! Create localisation matrix
if(allocated(INFER%L)) deallocate(INFER%L);allocate(INFER%L(INFER%nobs,MODEL%ntheta))
INFER%L=0
k=0
do j=1,MODEL%nTheta
    if(INFER%theta(j)%parType/=FIX) then
        INFER%L(:,j)=k+INFER%theta(j)%indx
        k=k+INFER%theta(j)%nval
    endif
enddo
!----------------------
INFER%nInfer=INFER%nFit+sum(INFER%nremnant)+sum(INFER%nXerror)+sum(INFER%nXb)+sum(INFER%nYb)
! Xtra model setup (derived parameters, states)
call XtraSetup(MODEL,err,mess)
if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
nstate=MODEL%nstate
!----------------------
! Derived parameters
!----------------------
if(INFER%nVar+INFER%nSTOK==0) then
    INFER%nDpar=MODEL%nDpar
    if(allocated(INFER%DparName)) deallocate(INFER%DparName);allocate(INFER%DparName(INFER%nDpar))
    if(allocated(INFER%DparIndx)) deallocate(INFER%DparIndx);allocate(INFER%DparIndx(1))
    INFER%DparName=MODEL%DparName
    INFER%DparIndx=1
else
    ! enumarate all combinations of VAR par values
    dummyIndx(1)=1
    dummyCombo(1,:)=INFER%L(1,:)
    k=1
    do i=2,INFER%nobs
        match=.false.
        do j=1,k
            if(all(INFER%L(i,:)-dummyCombo(j,:)==0)) then; match=.true.;exit;endif
        enddo
        if(.not.match) then ! new combination, save it
            k=k+1
            dummyIndx(k)=i
            dummyCombo(k,:)=INFER%L(i,:)
        endif
    enddo
    if(allocated(INFER%DparIndx)) deallocate(INFER%DparIndx);allocate(INFER%DparIndx(k))
    INFER%DparIndx=dummyIndx(1:k)
    ! Dpar number and names
    INFER%nDpar=MODEL%nDpar*size(INFER%DparIndx)
    if(allocated(INFER%DparName)) deallocate(INFER%DparName);allocate(INFER%DparName(INFER%nDpar))
    if(INFER%nDpar>0) then
        k=0
        do i=1,size(INFER%DparIndx)
            do j=1,size(MODEL%DparName)
                k=k+1
                INFER%DparName(k)=trim(MODEL%DparName(j))//"_"//trim(number_string(i))
            enddo
        enddo
    endif
endif

!----------------------
! Optional prior correlation
!----------------------
! get dimension of matrix and allocate
k=INFER%nFit
do j=1,model%nY
    k=k+size(INFER%Prior_RemnantSigma(j)%p)
enddo
if(allocated(INFER%priorM)) deallocate(INFER%priorM);allocate(INFER%priorM(k,k))
call Config_Read_PriorCorr(file=trim(priorCorrFile),exist=match,&
                           C=INFER%priorM,err=err,mess=mess)
if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
if(.not.match) then ! user didn't specify anything,deallocate
    deallocate(INFER%priorM)
else
    ! Make sure all diagonal terms are equal to one
    do i=1,k
        if(INFER%priorM(i,i)/=1._mrk) then
            err=1
            mess=trim(procname)//':'//trim(BaM_message(14))
            return
        endif
    enddo
    ! invert correlation
    call choles_invrt(ainv=INFER%priorM,logDet=INFER%priorLogDet,&
                      posDefinite=feas,err=err,message=mess)
    if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
    if(.not.feas) then;err=1;mess=trim(procname)//':'//trim(BaM_message(15));return;endif
    ! compute Cinv-I
    do i=1,k
        INFER%priorM(i,i)=INFER%priorM(i,i)-1._mrk
    enddo
endif

! Conclude by writing INFO file
call writeBamInfo(infoFile,err,mess)
if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif

end subroutine LoadBamObjects

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine BaM_ReadData(file,nrow,ncol,nHeader,&
                        XCol,XuCol,XbCol,XbindxCol,&
                        YCol,YuCol,YbCol,YbindxCol,&
                        X,Y,Xu,Yu,Xb,Yb,Xbindx,Ybindx,&
                        err,mess)
!^**********************************************************************
!^* Purpose: Read data file
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified: 22/07/2015
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1.file, DataFile
!^*    2.nrow (WITHOUT header lines)
!^*    3.ncol
!^*    4.nHeader, number of skipped header lines
!^*    5 and beyond. Columns for X/Y and their uncertainties
!^* OUT
!^*    1.X,Y,Xu,Yu,Xb,Yb,Xbindx,Ybindx
!^*    2.err, error code; <0:Warning, ==0:OK, >0: Error
!^*    3.mess, error message
!^**********************************************************************
use utilities_dmsl_kit,only:getSpareUnit
character(*), intent(in)::file
integer(mik), intent(in)::nrow,ncol,nHeader,&
                          XCol(:),XuCol(:),XbCol(:),XbindxCol(:),&
                          YCol(:),YuCol(:),YbCol(:),YbindxCol(:)
real(mrk),intent(out)::X(:,:),Y(:,:),Xu(:,:),Yu(:,:),Xb(:,:),Yb(:,:)
integer(mik),intent(out)::Xbindx(:,:),Ybindx(:,:)
integer(mik), intent(out)::err
character(*),intent(out)::mess
!locals
character(250),parameter::procname='BaM_ReadData'
integer(mik)::i,unt
real(mrk)::W(nrow,ncol),&
           fooX(size(Xbindx,dim=1),size(Xbindx,dim=2)),&
           fooY(size(Ybindx,dim=1),size(Ybindx,dim=2))

err=0;mess='';

call getSpareUnit(unt,err,mess)
if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif

open(unit=unt,file=trim(file), status='old',iostat=err)
if(err>0) then;call BaM_ConsoleMessage(messID_Open,trim(file));endif

do i=1,nHeader
    read(unt,*,iostat=err)
    if(err/=0) then;call BaM_ConsoleMessage(messID_Read,trim(file));endif
enddo

do i=1,nrow
    read(unt,*,iostat=err) W(i,:)
    if(err/=0) then;call BaM_ConsoleMessage(messID_Read,trim(file));endif
enddo
close(unt)
X=W(:,XCol)
Y=W(:,YCol)
call Wextract(XuCol,W,Xu)
call Wextract(YuCol,W,Yu)
call Wextract(XbCol,W,Xb)
call Wextract(YbCol,W,Yb)
call Wextract(XbindxCol,W,fooX);Xbindx=nint(fooX)
call Wextract(YbindxCol,W,fooY);Ybindx=nint(fooY)

end subroutine BaM_ReadData

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine BaM_Fit(theta0,RemnantSigma0, &!initial values for teta and remnant std
                   theta_std0,RemnantSigma_std0,& ! initial values for the std of jump distribution
                   nAdapt,nCycles,&
                   MinMoveRate,MaxMoveRate,&
                   DownMult,UpMult,&
                   OutFile,MonitorFile, & ! Output and monitoring files (for MCMC samples)
                   err,mess)! error handling
!^**********************************************************************
!^* Purpose: Performs MCMC sampling
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified:10/07/2015
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1 and beyond: properties of the MCMC sampler
!^*    last. OutFile, Output file (for MCMC samples)
!^* OUT
!^*    1.err, error code; <0:Warning, ==0:OK, >0: Error
!^*    2.mess, error message
!^**********************************************************************
use MCMCStrategy_tools, only:Adaptive_Metro_OAAT
use utilities_dmsl_kit, only:number_string

! INPUTS
real(mrk), intent(in)::theta0(:),theta_std0(:)
type(parlist), intent(in)::RemnantSigma0(:),RemnantSigma_std0(:)
integer(mik), intent(in)::nAdapt,nCycles
real(mrk), intent(in)::MinMoveRate,MaxMoveRate,DownMult,UpMult
character(*), intent(in)::OutFile,MonitorFile
! OUTPUTS
integer(mik), intent(out)::err
character(*),intent(out)::mess
!locals
character(250),parameter::procname='BaM_Fit'
real(mrk)::fx
real(mrk), allocatable:: start(:), startStd(:),Dpar(:)
character(100), allocatable:: headers(:)

err=0;mess=''

! Size checks
if( size(theta0)/=size(theta_std0) .or. size(theta0)/=INFER%nFit .or. &
    size(theta_std0)/=INFER%nFit) then
    err=1;mess=trim(procname)//':'//trim(BaM_message(3));return
endif
if( size(RemnantSigma0)/=size(RemnantSigma_std0) .or. size(RemnantSigma0)/=size(INFER%Prior_RemnantSigma) .or. &
    size(RemnantSigma_std0)/=size(INFER%Prior_RemnantSigma)) then
    err=1;mess=trim(procname)//':'//trim(BaM_message(3));return
endif

! specify starting values & starting std's
if(allocated(start)) deallocate(start);allocate(start(INFER%nInfer))
if(allocated(startStd)) deallocate(startStd);allocate(startStd(INFER%nInfer))
call Packer(theta=theta0,remnant=RemnantSigma0,X=INFER%X,bias=bias_init,OUT=start)
call Packer(theta=theta_std0,remnant=RemnantSigma_std0,X=INFER%Xu,bias=bias_initStd,OUT=startStd)

! get headers
if(allocated(headers)) deallocate(headers);allocate(headers(INFER%nInfer+1+INFER%nDpar))
call GetHeaders(headers)
! allocate Dpars
if(allocated(Dpar)) deallocate(Dpar);allocate(Dpar(INFER%nDpar))

! GO!
call Adaptive_Metro_OAAT(f=Posterior_wrapper,x=start,&
                fx=fx,fAux=Dpar,std=startStd,&
                nAdapt=nAdapt,nCycles=nCycles,&
                MinMoveRate=MinMoveRate,MaxMoveRate=MaxMoveRate,&
                DownMult=DownMult,UpMult=UpMult,&
                OutFile=OutFile,MonitorFile=MonitorFile,headers=headers,&
                err=err,mess=mess)

end subroutine BaM_Fit

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine BaM_LoadMCMC(MCMCFile,burnFactor,Nslim,& ! Read properties
                        maxpost,mcmc,logpost,err,mess)!out
!^**********************************************************************
!^* Purpose: Load MCMC file
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified:03/08/2015
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1. MCMCFile, containing MCMC samples
!^*    2. burnFactor, burnin factor
!^*    3. Nslim, only one line every Nslim will be used
!^* OUT
!^*    1.maxpost, maxpost parameter vector
!^*    2.mcmc, mcmc parameter samples
!^*    3.LogPost, associated (unnormalized) log-posterior densities
!^*    4.err, error code; <0:Warning, ==0:OK, >0: Error
!^*    5.mess, error message
!^**********************************************************************
use DataRW_tools,only:DatRead
character(*), intent(in)::MCMCFile
integer(mik),intent(in)::Nslim
real(mrk),intent(in)::burnFactor
integer(mik), intent(out)::err
character(*),intent(out)::mess
real(mrk),pointer::maxpost(:),mcmc(:,:),logpost(:)
! locals
character(250),parameter::procname='BaM_LoadMCMC'
integer(mik)::N,Nburn,Nread,i,k(1)
character(250)::headers(INFER%nInfer+1+INFER%nDpar)
real(mrk), pointer::y(:,:)
integer(mik), allocatable::rows(:)

err=0;mess=''
! read whole MCMC file
call DatRead(file=trim(MCMCFile),ncol=INFER%nInfer+1+INFER%nDpar,y=y,headers=headers,err=err,mess=mess)
if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif

! Slim MCMC
N=size(y,dim=1)
Nburn=NINT(N*BurnFactor)
Nread=size((/(i,i= Nburn+1,N,Nslim)/))
if(allocated(rows)) deallocate(rows);allocate(rows(Nread))
rows=(/(i,i= Nburn+1,N,Nslim)/)
if(associated(mcmc)) nullify(mcmc);allocate(mcmc(Nread,INFER%nInfer+INFER%nDpar))
mcmc(:,1:INFER%nInfer)=y(rows,1:INFER%nInfer)
if(INFER%nDpar>0) mcmc(:,(INFER%nInfer+1):(INFER%nInfer+INFER%nDpar))=y(rows,(INFER%nInfer+1+1):(INFER%nInfer+1+INFER%nDpar))
if(associated(LogPost)) nullify(LogPost);allocate(LogPost(Nread))
LogPost=y(rows,INFER%nInfer+1)

! get maxpost parameter vector
if(associated(maxpost)) nullify(maxpost);allocate(maxpost(INFER%nInfer))
k=maxloc(logpost)
maxpost=mcmc(k(1),1:INFER%nInfer)

end subroutine BaM_LoadMCMC

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine BaM_CookMCMC(mcmc,logpost,OutFile,err,mess)
!^**********************************************************************
!^* Purpose: Cook MCMC samples
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified:25/01/2017
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1. logpost, (unnormalized) log-posterior density
!^*    2. OutFile, File where cooked MCMC samples are stored
!^* OUT
!^*    1.err, error code; <0:Warning, ==0:OK, >0: Error
!^*    2.mess, error message
!^* INOUT
!^*    1. mcmc, MCMC samples
!^**********************************************************************
use ModelLibrary_tools, only:MDL_Linear,MDL_BaRatin,MDL_GR4J,MDL_Ortho,&
                             MDL_Sediment,MDL_SFD
use utilities_dmsl_kit,only:number_string,getSpareUnit
real(mrk),intent(inout)::mcmc(:,:)
real(mrk),intent(in)::logpost(:)
character(*),intent(in)::OutFile
integer(mik), intent(out)::err
character(*),intent(out)::mess
! locals
character(250),parameter::procname='BaM_CookMCMC'
integer(mik)::i,n,p,ncol,unt
character(250)::fmt,sfmt
character(250)::head(INFER%nInfer+1+INFER%nDpar)

err=0;mess=''
n=size(mcmc,dim=1);p=size(mcmc,dim=2)
call GetHeaders(head)

! Write cooked MCMC to file
ncol=size(mcmc,dim=2)+1
fmt='('//trim(number_string(ncol))//trim(fmt_numeric)//')'
sfmt='('//trim(number_string(ncol))//trim(fmt_string)//')'
! open and write
call getSpareUnit(unt,err,mess)
if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
open(unit=unt,file=trim(OutFile), status='replace', iostat=err)
if(err>0) then;call BaM_ConsoleMessage(messID_Open,trim(OutFile));endif
write(unt,trim(sfmt)) head
if(err>0) then;call BaM_ConsoleMessage(messID_Write,trim(OutFile));endif
do i=1,n
    write(unt,trim(fmt)) mcmc(i,1:INFER%nInfer),logpost(i),mcmc(i,(INFER%nInfer+1):(INFER%nInfer+INFER%nDpar))
    if(err>0) then;call BaM_ConsoleMessage(messID_Write,trim(OutFile));endif
enddo
close(unt)
end subroutine BaM_CookMCMC

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine BaM_Residual(maxpost,OutFile,err,mess)
!^**********************************************************************
!^* Purpose: Runs the model on calibration data using maxpost par
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified: 04/08/2015
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1.maxpost, parameter vector
!^*    2.OutFile, result file
!^* OUT
!^*    1.err, error code; <0:Warning, ==0:OK, >0: Error
!^*    2.mess, error message
!^**********************************************************************
use utilities_dmsl_kit,only:getspareunit,number_string
real(mrk),intent(in)::maxpost(:)
character(*),intent(in)::OutFile
integer(mik), intent(out)::err
character(*),intent(out)::mess
!locals
character(250),parameter::procname='BaM_Residual'
real(mrk), allocatable::Xtrue(:,:),Ysim(:,:),Yunb(:,:),state(:,:),res(:,:),stdres(:,:),sfpar(:)
logical:: feas
integer(mik)::i,j,k,unt,ncol
real(mrk)::sig
character(250)::fmt,sfmt
character(250),allocatable::head(:)

err=0;mess=''
! run model
if(allocated(Xtrue)) deallocate(Xtrue)
if(allocated(Ysim)) deallocate(Ysim)
allocate(Xtrue(infer%nobs,model%nX),Ysim(infer%nobs,model%nY),state(infer%nobs,model%nState))
call GetSim_calibration(par=maxpost,Xtrue=Xtrue,Ysim=Ysim,state=state,feas=feas,err=err,mess=mess)
if(err>0) then;mess=trim(procname)//':'//trim(mess);return;endif
if(.not. feas) then
    err=1;mess=trim(procname)//':'//BaM_message(5);return
endif
! get unbiased outputs
if(allocated(Yunb)) deallocate(Yunb);allocate(Yunb(infer%nobs,model%nY))
Yunb=INFER%Y
k=INFER%nFit+sum(INFER%nremnant)+sum(INFER%nXerror)+sum(INFER%nXb)
do j=1,model%nY
    if(INFER%nYb(j)==0) cycle ! no bias to correct
    do i=1,infer%nobs
        if(INFER%Ybindx(i,j)>0 .and. INFER%Y(i,j)/=infer%mv) then ! correct obs
            Yunb(i,j)=INFER%Y(i,j)-INFER%Yb(i,j)*maxpost(k+INFER%Ybindx(i,j))
        endif
    enddo
    k=k+INFER%nYb(j)
enddo
! get residuals
if(allocated(res)) deallocate(res)
if(allocated(stdres)) deallocate(stdres)
allocate(res(infer%nobs,model%nY),stdres(infer%nobs,model%nY))
res=Yunb-Ysim
where(infer%Y==infer%mv) res=infer%mv
k=INFER%nFit
do j=1,model%nY
    if(allocated(sfpar)) deallocate(sfpar)
    allocate(sfpar(infer%nremnant(j)))
    sfpar=maxpost((k+1):(k+infer%nremnant(j)))
    k=k+infer%nremnant(j)
    do i=1,infer%nobs
        call Sigmafunk_Apply(infer%RemnantSigma_funk(j),sfpar, Ysim(i,j), sig, err, mess)
        if(err>0) then;mess=trim(procname)//':'//trim(mess);return;endif
        stdres(i,j)=res(i,j)/sqrt(sig**2+INFER%Yu(i,j)**2)
    enddo
enddo
where(infer%Y==infer%mv) stdres=infer%mv

! write to file
ncol=2*model%nX+5*model%nY+model%nState
fmt='('//trim(number_string(ncol))//trim(fmt_numeric)//')'
sfmt='('//trim(number_string(ncol))//trim(fmt_string)//')'
! assemble headers
if(allocated(head)) deallocate(head);allocate(head(ncol))
k=0
do i=1,model%nX
    head(k+1)="X"//trim(number_string(i))//"_obs";k=k+1
enddo
do i=1,model%nX
    head(k+1)="X"//trim(number_string(i))//"_true";k=k+1
enddo
do i=1,model%nY
    head(k+1)="Y"//trim(number_string(i))//"_obs";k=k+1
enddo
do i=1,model%nY
    head(k+1)="Y"//trim(number_string(i))//"_unbiased";k=k+1
enddo
do i=1,model%nY
    head(k+1)="Y"//trim(number_string(i))//"_sim";k=k+1
enddo
do i=1,model%nY
    head(k+1)="Y"//trim(number_string(i))//"_res";k=k+1
enddo
do i=1,model%nY
    head(k+1)="Y"//trim(number_string(i))//"_stdres";k=k+1
enddo
do i=1,model%nState
    head(k+1)=model%StateName(i);k=k+1
enddo

! open and write
call getSpareUnit(unt,err,mess)
if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
open(unit=unt,file=trim(OutFile), status='replace', iostat=err)
if(err>0) then;call BaM_ConsoleMessage(messID_Open,trim(OutFile));endif
write(unt,trim(sfmt)) head
if(err>0) then;call BaM_ConsoleMessage(messID_Write,trim(OutFile));endif
do i=1,infer%nobs
    write(unt,trim(fmt)) infer%X(i,:),Xtrue(i,:),infer%Y(i,:),Yunb(i,:),&
                         Ysim(i,:),res(i,:),stdres(i,:),state(i,:)
    if(err>0) then;call BaM_ConsoleMessage(messID_Write,trim(OutFile));endif
enddo
close(unt)

end subroutine BaM_Residual

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine BaM_Prediction(mcmc,maxpost,nsim,&
                          Xspag,DoParametric,DoRemnant,&
                          SpagFiles,DoTranspose,DoEnvelop,EnvelopFiles,&
                          DoState,SpagFiles_S,DoTranspose_S,DoEnvelop_S,EnvelopFiles_S,&
                          PrintCounter,MonitorFile,err,mess)
!^**********************************************************************
!^* Purpose: Runs the model on calibration data using maxpost par
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified: 09/10/2015
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1.[mcmc], MCMC simulations. If not present, will sample from prior
!^*    2.[maxpost], maxpost parameters. If not present, will compute the prior mode
!^*    3.[nsim], only used for prior sampling, not used is mcmc is provided
!^*    4.Xspag, spaghettis of inputs (size nX, each element comprising a matrix)
!^*    5.DoParametric, propagate parametric uncertainty? (0/1)
!^*    6.DoRemnant, size nY. Propagate remnant uncertainty? (0/1)
!^*    7.SpagFiles, size nY. Name of spag file for each output variable
!^*    8.DoTranspose, size nY, transpose spaghettis? (so that each column is one spag)
!^*    9.DoEnvelop, size nY, compute uncertainty envelops?
!^*    10.EnvelopFiles, size nY, Name of envelop file for each output variable
!^*    11.DoState, size model%nState, compute predicted states?
!^*    12.SpagFiles_S, size model%nState. Name of spag file for each state variable
!^*    13.DoTranspose_S, size model%nState, transpose state spaghettis? (so that each column is one spag)
!^*    14.DoEnvelop_S, size model%nState, compute state uncertainty envelops?
!^*    15.EnvelopFiles_S, size model%nState, Name of envelop file for each state variable
!^*    16.[PrintCounter], print a counter in console during computations?
!^*    17.[MonitorFile], file to monitor computations
!^* OUT
!^*    1.err, error code; <0:Warning, ==0:OK, >0: Error
!^*    2.mess, error message
!^**********************************************************************
use Distribution_tools, only: Generate,GAUSS
use utilities_dmsl_kit,only:getSpareUnit,number_string

real(mrk),intent(in),optional::mcmc(:,:),maxpost(:)
integer(mik),intent(in),optional::nsim
type(XspagType),intent(in)::Xspag
logical, intent(in)::DoParametric,DoRemnant(:),DoTranspose(:),DoEnvelop(:),&
                     DoState(:),DoTranspose_S(:),DoEnvelop_S(:)
logical, intent(in), optional::PrintCounter
character(*),intent(in)::SpagFiles(:),EnvelopFiles(:),&
                         SpagFiles_S(:),EnvelopFiles_S(:)
character(*),intent(in),optional::MonitorFile
integer(mik), intent(out)::err
character(*),intent(out)::mess
!locals
character(250),parameter::procname='BaM_Prediction'
logical, parameter::printC_def=.false.
integer(mik),parameter::Nrefresh=1 ! used for counter
logical::printC,IsMaxpost,feas
integer(mik)::Nrep,next,rep,i,t,indx,compt,unt(MODEL%nY),unt_S(MODEL%nState),f,untMonitor
real(mrk),allocatable::mode(:),MC(:,:),X(:,:),Ysim(:,:),Dpar(:),state(:,:)
real(mrk)::theta(INFER%nFit),sig,dev
character(250)::fmt

err=0;mess=''
! check feasability
if(any(INFER%theta%parType==VAR) .or. any(INFER%theta%parType==STOK)) then
    err=1;mess=trim(procname)//':'//trim(BaM_Message(13));return
endif
! handle optional counter
if(present(PrintCounter)) then;printC=PrintCounter;else;printC=printC_def;endif
! Number of Monte Carlo simulations
if(.not.DoParametric .and. all(.not.DoRemnant) .and. all(Xspag%nspag==1)) then ! This is a maxpost simulation!
    IsMaxpost=.true.
    Nrep=1
else
    IsMaxpost=.false.
    if(present(mcmc)) then
        Nrep=size(mcmc,dim=1)
    else
        if(.not.present(nsim)) then
            err=1;mess=trim(procname)//':'//trim(BaM_Message(7));return
        endif
        Nrep=Nsim
    endif
endif
allocate(mode(INFER%nInfer),MC(Nrep,INFER%nInfer))
! handle prior simulation
if(present(mcmc)) then ! use provided mcmc/maxpost
    mode=maxpost
    if(.not.IsMaxPost) MC=mcmc
else ! get prior mode and prior simulations
    if(IsMaxPost) then
        call BaM_GetPriorMC(mode=mode,err=err,mess=mess)
        if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
    else
        call BaM_GetPriorMC(MC=MC,mode=mode,err=err,mess=mess)
        if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
    endif
endif
! handle files
unt=undefIN;unt_S=undefIN
do f=1,MODEL%nY
    call getSpareUnit(unt(f),err,mess)
    if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
    open(unt(f),file=trim(SpagFiles(f)),status='replace')
enddo
fmt='('//trim(number_string(Xspag%nobs))//trim(fmt_numeric)//')'
if(any(DoState)) then
    do f=1,MODEL%nState
        if(DoState(f)) then
            call getSpareUnit(unt_S(f),err,mess)
            if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
            open(unt_S(f),file=trim(SpagFiles_S(f)),status='replace')
        endif
    enddo
endif
if(present(MonitorFile)) then
    call getSpareUnit(untMonitor,err,mess)
    if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
    open(unit=untMonitor,file=MonitorFile,status='replace')
endif

! GO!
allocate(X(Xspag%nobs,MODEL%nX),Ysim(Xspag%nobs,MODEL%nY))
allocate(Dpar(INFER%nDpar),state(Xspag%nobs,MODEL%nState))
do rep=1,Nrep
    ! reinitialize
    X=MODEL%unfeasFlag;Ysim=MODEL%unfeasFlag
    !print counter if requested
    if(printC .or. present(MonitorFile)) then
        if(rep==1) next=Nrefresh
        if(rep>=0.01*next*Nrep) then
            if(printC) then;write(*,'(I4,A)') next,'% DONE';endif
            if(present(MonitorFile)) call writeMonitor(untMonitor,rep,Nrep)
            next=next+Nrefresh
        endif
    endif
    ! Create input matrix
    do i=1,MODEL%nX
        ! Define spag index. If nspag<rep, cycle through X-spags
        indx=modulo(rep,Xspag%nspag(i));if(indx==0)indx=Xspag%nspag(i)
        ! pack spag into X matrix
        X(:,i)=Xspag%spag(i)%m(:,indx)
    enddo
    ! define theta vector and run model
    if(INFER%nFit>0) then
        if(DoParametric) then !from MCMC
            theta=MC(rep,1:INFER%nFit)
        else ! maxpost
            theta=mode(1:INFER%nFit)
        endif
    endif
    ! Run model
    call ApplyModel_BaM(model=MODEL,X=X,theta=theta,Y=Ysim,&
                    Dpar=Dpar,state=state,feas=feas,err=err,mess=mess)
    if(err>0) then
        mess=trim(procname)//':'//trim(mess)
        call CloseAllFiles(unt,unt_S,untMonitor)
        return
    endif

    ! write states to file
    do i=1,MODEL%nState
        if(DoState(i)) then
            write(unt_S(i),trim(fmt)) state(:,i)
        endif
    enddo

    ! add structural error
    compt=INFER%nFit
    do i=1,MODEL%nY
        if(DoRemnant(i)) then
            do t=1,Xspag%nobs
                if(Ysim(t,i)==MODEL%unfeasflag) cycle ! don't add an error to junk...
                ! compute sdev
                call Sigmafunk_Apply(funk=infer%RemnantSigma_funk(i),&
                                 par=MC(rep,(compt+1):(compt+INFER%nremnant(i))),&
                                 Y=Ysim(t,i),res=sig,&
                                 err=err,mess=mess)
                if(err>0) then
                    mess=trim(procname)//':'//trim(mess)
                    feas=.false.
                    call CloseAllFiles(unt,unt_S,untMonitor)
                    return
                endif
                ! generate error
                call Generate(DistId=GAUSS,par=(/0._mrk,sig/),&
                          gen=dev,feas=feas,err=err,mess=mess)
                if(err/=0) then
                    mess=trim(procname)//':'//trim(mess)
                    call CloseAllFiles(unt,unt_S,untMonitor)
                    return
                endif
                if(feas) Ysim(t,i)=Ysim(t,i)+dev
            enddo
        endif
        compt=compt+INFER%nremnant(i)
        write(unt(i),trim(fmt)) Ysim(:,i)
    enddo
enddo
deallocate(mode,MC,X,Ysim)
call CloseAllFiles(unt,unt_S)
! Post-process
if(present(MonitorFile)) call writeMonitor(untMonitor,Nrep,Nrep)
call BaM_ProcessSpag(Xspag%nobs,nrep,SpagFiles,DoTranspose,DoEnvelop,&
                     EnvelopFiles,PrintCounter,err,mess)
if(err/=0) then
    mess=trim(procname)//':'//trim(mess)
    if(present(MonitorFile)) close(untMonitor)
    return
endif
if(present(MonitorFile)) call writeMonitor(untMonitor,Nrep+1,Nrep)
if(any(DoState)) then
    if(present(MonitorFile)) call writeMonitor(untMonitor,Nrep+2,Nrep)
    call BaM_ProcessSpag(Xspag%nobs,nrep,pack(SpagFiles_S,DoState),&
                         pack(DoTranspose_S,DoState),pack(DoEnvelop_S,DoState),&
                         pack(EnvelopFiles_S,DoState),PrintCounter,err,mess)
    if(err/=0) then
        mess=trim(procname)//':'//trim(mess)
        if(present(MonitorFile)) close(untMonitor)
        return
    endif
    if(present(MonitorFile)) call writeMonitor(untMonitor,Nrep+3,Nrep)
endif
if(present(MonitorFile)) close(untMonitor)

end subroutine BaM_Prediction

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine BaM_ReadSpag(Xspag,Spag_Files,err,mess)
!^**********************************************************************
!^* Purpose: Read spaghettis of input variables
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified: 09/10/2015
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1.Spag_Files, Files containing spaghettis
!^* INOUT
!^*    1.Xspag, spaghetti object
!^* OUT
!^*    1.err, error code; <0:Warning, ==0:OK, >0: Error
!^*    2.mess, error message
!^**********************************************************************
use utilities_dmsl_kit,only:getSpareUnit
type(XspagType), intent(inout)::Xspag
character(*), intent(in)::Spag_Files(:)
integer(mik), intent(out)::err
character(*),intent(out)::mess
!locals
character(250),parameter::procname='BaM_ReadSpag'
integer(mik)::i,j,unt

err=0;mess=''
do i=1,model%nX
    if(allocated(Xspag%spag(i)%m)) deallocate(Xspag%spag(i)%m)
    allocate(Xspag%spag(i)%m(Xspag%nobs,Xspag%nspag(i)))
    call getSpareUnit(unt,err,mess)
    if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
    open(unit=unt,file=trim(Spag_Files(i)), status='old')
    do j=1,Xspag%nobs
        read(unt,*,iostat=err) Xspag%spag(i)%m(j,:)
        if(err>0) then;call BaM_ConsoleMessage(messID_Open,trim(Spag_Files(i)));endif
    enddo
    close(unt)
enddo

end subroutine BaM_ReadSpag

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine BaM_ProcessSpag(nobs,nrep,SpagFiles,DoTranspose,DoEnvelop,&
                           EnvelopFiles,PrintCounter,err,mess)
!^**********************************************************************
!^* Purpose: Processing spaghettis from a prediction experiment
!^*          (transposing & computing enveloppes)
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
!^*    1.nobs, number of observations (lines in spag file)
!^*    2.nrep, number of replications (columns in spag file)
!^*    3.SpagFiles, file containing prediction spaghettis
!^*    4.DoTranspose, transpose spaghettis?
!^*    5.DoEnvelop, compute uncertainty envelops?
!^*    6.EnvelopFiles, file where envelops are written
!^*    6.PrintCounter, print message during execution?
!^* OUT
!^*    1.err, error code; <0:Warning, ==0:OK, >0: Error
!^*    2.mess, error message
!^**********************************************************************
use EmpiricalStats_tools, only:GetEmpiricalStats,GetEmpiricalQuantile
use numerix_dmsl_kit, only:quicksort
use utilities_dmsl_kit,only:getSpareUnit,number_string

integer(mik),intent(in)::nobs,nrep
character(*),intent(in)::SpagFiles(:),EnvelopFiles(:)
logical,intent(in)::DoTranspose(:),DoEnvelop(:),PrintCounter
integer(mik), intent(out)::err
character(*),intent(out)::mess
! locals
character(250),parameter::procname='BaM_ProcessSpag'
integer(mik),parameter::nEnv=7
logical,parameter::removeMV=.true.
real(mrk),parameter::maxMVrate=0.75_mrk
character(250), dimension(nEnv),parameter::head_envelop=(/"Median       ",&
                                                       "q2.5         ",&
                                                       "q97.5        ",&
                                                       "q16          ",&
                                                       "q84          ",&
                                                       "Mean         ",&
                                                       "Stdev        "/)
integer(mik)::i,j,unt
real(mrk)::Env(nEnv)
character(250)::fmt,fmt2,fmt3
real(mrk),allocatable::Spag(:,:),arr(:)
logical::mask(nrep)

err=0;mess=''
fmt='('//trim(number_string(nrep))//trim(fmt_numeric)//')'
fmt2='('//trim(number_string(nEnv))//trim(fmt_numeric)//')'
fmt3='('//trim(number_string(nEnv))//trim(fmt_string)//')'

if(PrintCounter) call BaM_ConsoleMessage(15,'')
do i=1,size(SpagFiles)
    if( DoTranspose(i) .or. DoEnvelop(i) ) then ! need to allocate big Spag matrix if not done yet
        if(.not.allocated(Spag)) allocate(Spag(nrep,nobs))
    else ! nothing to do, cycle
        cycle
    endif
    ! read Yspag
    call getSpareUnit(unt,err,mess)
    if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
    open(unt,file=trim(SpagFiles(i)),status='old',iostat=err)
    if(err/=0) then;call BaM_ConsoleMessage(messID_Open,trim(SpagFiles(i)));endif
    do j=1,nrep
        read(unt,*,iostat=err) Spag(j,:)
        if(err/=0) then;call BaM_ConsoleMessage(messID_Read,trim(SpagFiles(i)));endif
    enddo
    close(unt)
    if(DoTranspose(i)) then
        ! rewrite Spag transposed
        call getSpareUnit(unt,err,mess)
        if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
        open(unt,file=trim(SpagFiles(i)),status='replace',iostat=err)
        if(err/=0) then;call BaM_ConsoleMessage(messID_Open,trim(SpagFiles(i)));endif
        do j=1,nobs
            write(unt,trim(fmt),iostat=err) Spag(:,j)
            if(err/=0) then;call BaM_ConsoleMessage(messID_Write,trim(SpagFiles(i)));endif
        enddo
        close(unt)
    endif
    if(DoEnvelop(i)) then
        if(nrep<2) then ! no point computing an envelop with a single rep!!!
            call BaM_ConsoleMessage(16,'')
        else
            ! Open envelop file
            call getSpareUnit(unt,err,mess)
            if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
            open(unt,file=trim(EnvelopFiles(i)),status='replace',iostat=err)
            if(err/=0) then;call BaM_ConsoleMessage(messID_Open,trim(EnvelopFiles(i)));endif
            write(unt,trim(fmt3),iostat=err) head_envelop
            if(err/=0) then;call BaM_ConsoleMessage(messID_Write,trim(EnvelopFiles(i)));endif
            ! Compute envelops
            do j=1,nobs
                mask=(Spag(:,j)==model%unfeasFlag) .or. (Spag(:,j)/=Spag(:,j)) ! the second condition aims at detecting NAs since NA/=NA is true
                if(removeMV .and. real(count(mask),mrk)/real(nrep,mrk) < maxMVrate) then ! remove unfeasible values and NAs
                    if(allocated(arr)) deallocate(arr);allocate(arr(count(.not.mask)))
                    arr=pack(Spag(:,j),.not.mask)
                else
                    if(allocated(arr)) deallocate(arr);allocate(arr(nrep))
                    arr=Spag(:,j)
                endif
                if(any(arr==model%unfeasFlag .or. arr/=arr)) then ! Do nothing
                    Env=model%unfeasFlag
                else
                    ! Sort
                    call quicksort(arr=arr,ascnd=.true.,err=err)
                    if(err>0) then;mess=trim(procname)//':'//trim(mess);return;endif
                    ! Get quantiles
                    call GetEmpiricalQuantile(p=0.5_mrk,x=arr,IsXSorted=.true.,q=Env(1),err=err,mess=mess)
                    if(err>0) then;mess=trim(procname)//':'//trim(mess);return;endif
                    call GetEmpiricalQuantile(p=0.025_mrk,x=arr,IsXSorted=.true.,q=Env(2),err=err,mess=mess)
                    if(err>0) then;mess=trim(procname)//':'//trim(mess);return;endif
                    call GetEmpiricalQuantile(p=0.975_mrk,x=arr,IsXSorted=.true.,q=Env(3),err=err,mess=mess)
                    if(err>0) then;mess=trim(procname)//':'//trim(mess);return;endif
                    call GetEmpiricalQuantile(p=0.16_mrk,x=arr,IsXSorted=.true.,q=Env(4),err=err,mess=mess)
                    if(err>0) then;mess=trim(procname)//':'//trim(mess);return;endif
                    call GetEmpiricalQuantile(p=0.84_mrk,x=arr,IsXSorted=.true.,q=Env(5),err=err,mess=mess)
                    if(err>0) then;mess=trim(procname)//':'//trim(mess);return;endif
                    ! Get mean and stdev
                    call GetEmpiricalStats(x=arr,mean=Env(6),std=Env(7),err=err,mess=mess)
                    if(err>0) then;mess=trim(procname)//':'//trim(mess);return;endif
                endif
                ! Write to file
                write(unt,trim(fmt2),iostat=err) Env
                if(err/=0) then;call BaM_ConsoleMessage(messID_Write,trim(EnvelopFiles(i)));endif
            enddo
            close(unt)
        endif
    endif
enddo

end subroutine BaM_ProcessSpag

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine BaM_SummarizeMCMC(mcmc,logpost,OutFile,err,mess)
!^**********************************************************************
!^* Purpose: summarize MCMC samples
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified:04/08/2015
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1. mcmc, MCMC samples
!^*    2. logpost, (unnormalized) log-posterior density
!^*    3. OutFile, File where summary is written
!^* OUT
!^*    1.err, error code; <0:Warning, ==0:OK, >0: Error
!^*    2.mess, error message
!^**********************************************************************
use utilities_dmsl_kit,only:number_string,getSpareUnit
use EmpiricalStats_tools,only:GetEmpiricalStats

real(mrk),intent(in)::mcmc(:,:),logpost(:)
character(*),intent(in)::OutFile
integer(mik), intent(out)::err
character(*),intent(out)::mess
! locals
character(250),parameter::procname='BaM_SummarizeMCMC'
character(250)::head(INFER%nInfer+1+INFER%nDpar)
character(250), dimension(16),parameter::lefters=(/"N              ",&
                                                   "Minimum        ",&
                                                   "Maximum        ",&
                                                   "Range          ",&
                                                   "Mean           ",&
                                                   "Median         ",&
                                                   "Q10%           ",&
                                                   "Q25%           ",&
                                                   "Q75%           ",&
                                                   "Q90%           ",&
                                                   "St.Dev.        ",&
                                                   "Variance       ",&
                                                   "CV             ",&
                                                   "Skewness       ",&
                                                   "Kurtosis       ",&
                                                   "MaxPost        "/)
integer(mik),parameter::nstat=16
integer(mik)::i,p,k(1),unt
real(mrk),allocatable::summary(:,:)
real(mrk)::maxpost(size(mcmc,dim=2))
character(250)::fmt,sfmt,h(size(head)-1)

err=0;mess=''
p=size(mcmc,dim=2)
allocate(summary(nstat,p))
call GetHeaders(head)
! get maxpost parameter vector
k=maxloc(logpost)
maxpost=mcmc(k(1),:)
Summary(nstat,:)=maxpost
do i=1,p
    call GetEmpiricalStats(x=mcmc(:,i), all15=Summary(1:(nstat-1),i),err=err,mess=mess)
    if(err>0) then; mess=trim(procname)//':'//trim(mess);return;endif
enddo
! write to file
fmt='('//trim(fmt_string)//','//trim(number_string(p))//trim(fmt_numeric)//')'
sfmt='('//trim(number_string(p+1))//trim(fmt_string)//')'
! open and write
call getSpareUnit(unt,err,mess)
if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
open(unit=unt,file=trim(OutFile), status='replace', iostat=err)
if(err>0) then;call BaM_ConsoleMessage(messID_Open,trim(OutFile));endif
if(INFER%nDpar==0) then
    h=head(1:INFER%nInfer)
else
    h=(/head(1:INFER%nInfer),head( (INFER%nInfer+1+1):(INFER%nInfer+1+INFER%nDpar) )/)
endif
write(unt,trim(sfmt)) " ",h
if(err>0) then;call BaM_ConsoleMessage(messID_Write,trim(OutFile));endif
do i=1,nstat
    write(unt,trim(fmt)) lefters(i),Summary(i,:)
    if(err>0) then;call BaM_ConsoleMessage(messID_Write,trim(OutFile));endif
enddo
close(unt)
end subroutine BaM_SummarizeMCMC

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine BaM_computeDIC(mcmc,logpost,DICfile,MCMCfile,err,mess)
!^**********************************************************************
!^* Purpose: DIC computations, DIC = D(thetaHat) + 2*(effective nPar)
!^**********************************************************************
!^* Programmer: Ben Renard, INRAE Aix
!^**********************************************************************
!^* Last modified:23/02/2023
!^**********************************************************************
!^* Comments: The maxpost is taken as thetaHat in this sub, because it
!^*           belongs to the MCMC sample and is hence necessarily feasible.
!^*           Note that it is more common in the statistical literature
!^*           to use the posterior mean, but as discussed by the DIC
!^*           creators (Spiegelhalter 2002) this is just a choice.
!^*
!^*       3 versions of the DIC are provided:
!^*           1. DIC1 is the original version of Spiegelhalter (2002),
!^*                   based on the posterior MEAN of deviance D(theta)
!^*              DIC1 = D(thetaHat) + 2*(Dmean-D(thetaHat))
!^*           2. DIC2 is the alternative definition mentionned in
!^*                   Gelman (2014) and Spiegelhalter (2014),
!^*                   based on the posterior VARIANCE of D(theta)
!^*              DIC2 = D(thetaHat) + 2*(0.5*Dvar)
!^*           3. DIC3 is a third option used in Pooley and Marion (2018),
!^*                   based on both the posterior MEAN and VARIANCE of D(theta).
!^*                   Unlike DIC1 and DIC2, it does not use the deviance
!^*                   D(thetaHat) computed at point-estimate thetaHat,
!^*                   and therefore avoids having to make a choice for
!^*                   thetaHat (posterior mode, mean, median etc.)
!^*              DIC3 = Dmean+0.5*Dvar = (DIC1+DIC2)/2
!^**********************************************************************
!^* References:
!^*  * Spiegelhalter, D.J., Best, N.G., Carlin, B.P., van der Linde, A.:
!^*    Bayesian measures of model complexity and fit (with discussion).
!^*    J. R. Stat. Soc. B (2002)
!^*  * Gelman, A., Hwang, J. & Vehtari, A. Understanding predictive
!^*    information criteria for Bayesian models. Stat Comput 24,
!^*    997–1016 (2014). https://doi.org/10.1007/s11222-013-9416-2
!^*  * Spiegelhalter, D.J., Best, N.G., Carlin, B.P. and van der Linde, A.
!^*    The deviance information criterion: 12 years on.
!^*    J. R. Stat. Soc. B (2014), https://doi.org/10.1111/rssb.12062
!^*  * Pooley, C. M., & Marion, G. Bayesian model evidence as a practical
!^*    alternative to deviance information criterion. Royal Society Open
!^*    Science (2018), https://doi.org/10.1098/rsos.171519
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1. mcmc, MCMC samples
!^*    2. logpost, (unnormalized) log-posterior density
!^*    3. DICFile, File where DIC and associated statistics are written
!^*    4. MCMCFile, File where extended MCMC is written. File not written if empty string
!^* OUT
!^*    1.err, error code; <0:Warning, ==0:OK, >0: Error
!^*    2.mess, error message
!^**********************************************************************
use utilities_dmsl_kit,only:getSpareUnit,number_string
use EmpiricalStats_tools,only:GetEmpiricalStats

real(mrk),intent(in)::mcmc(:,:),logpost(:)
character(*),intent(in)::DICfile, MCMCfile
integer(mik), intent(out)::err
character(*),intent(out)::mess
! locals
character(250),parameter::procname='BaM_computeDIC'
real(mrk)::theta(INFER%nFit),Xtrue(INFER%nobs,MODEL%nX)
type(parlist)::gamma(MODEL%nY),Xbias(MODEL%nX),Ybias(MODEL%nY)
real(mrk)::prior,hyper,lkh,dev(size(mcmc,dim=1)),Dmean,Dmed,Dvar,Ds,Dmaxpost
integer(mik)::i,k(1),unt0,unt1,n
character(250)::head(INFER%nInfer+1+INFER%nDpar)
character(250)::fmt,sfmt,h(INFER%nInfer+5)
logical::feas,isnull,rewriteMCMC
character(250), dimension(8),parameter::lefters=(/"DIC1           ",&
                                                  "DIC2           ",&
                                                  "DIC3           ",&
                                                  "D(maxpost)     ",&
                                                  "E[D]           ",&
                                                  "VAR[D]         ",&
                                                  "MEDIAN[D]      ",&
                                                  "SDEV[D]        "/)

err=0;mess=''
n=size(mcmc,dim=1)
!---------------
! Extended MCMC file
rewriteMCMC=.not.(MCMCfile=='')
if(rewriteMCMC) then
    call getSpareUnit(unt0,err,mess)
    if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
    open(unit=unt0,file=trim(MCMCfile),status='replace',iostat=err)
    if(err>0) then;call BaM_ConsoleMessage(messID_Open,trim(MCMCfile));endif
    fmt='('//trim(number_string(INFER%nInfer+5))//trim(fmt_numeric)//')'
    sfmt='('//trim(number_string(INFER%nInfer+5))//trim(fmt_string)//')'
    ! Write header
    call GetHeaders(head)
    h=(/head(1:INFER%nInfer),'LogPost','LogPrior','LogHyper','LogLkh','Deviance'/)
    write(unt0,trim(sfmt)) h
    if(err>0) then;call BaM_ConsoleMessage(messID_Write,trim(MCMCFile));endif
endif
!---------------
! compute deviances
do i=1,n
    call UnfoldParVector(mcmc(i,1:INFER%nInfer),theta,gamma,Xtrue,Xbias,Ybias,err,mess)
    if(err>0) then; mess=trim(procname)//':'//trim(mess);return;endif
    call GetPriorLogPdf(theta,gamma,Xtrue,Xbias,Ybias,MODEL,INFER,prior,feas,isnull,err,mess)
    if(err>0) then; mess=trim(procname)//':'//trim(mess);return;endif
    if ( (.not. feas) .or. (isnull) ) then
        mess=trim(procname)//':'//trim(BaM_message(16));return
    endif
    call GetHyperLogPdf(theta,gamma,Xtrue,Xbias,Ybias,MODEL,INFER,hyper,feas,isnull,err,mess)
    if(err>0) then; mess=trim(procname)//':'//trim(mess);return;endif
    if ( (.not. feas) .or. (isnull) ) then
        mess=trim(procname)//':'//trim(BaM_message(17));return
    endif
    lkh=logpost(i)-prior-hyper
    dev(i)=-2._mrk*lkh
    if(rewriteMCMC) then
        ! write to extended MCMC file
        write(unt0,trim(fmt)) mcmc(i,1:INFER%nInfer),logpost(i),prior,hyper,lkh,dev(i)
        if(err>0) then;call BaM_ConsoleMessage(messID_Write,trim(MCMCfile));endif
    endif
enddo
if(rewriteMCMC) close(unt0)
k=maxloc(logpost)
Dmaxpost=dev(k(1))
!---------------
! Compute DIC
call GetEmpiricalStats(x=dev,mean=Dmean,median=Dmed,std=Ds,var=Dvar,err=err,mess=mess)
if(err>0) then; mess=trim(procname)//':'//trim(mess);return;endif
!---------------
! write to file
fmt='('//trim(fmt_string)//','//trim(fmt_numeric)//')'
!sfmt='('//trim(number_string(p+1))//trim(fmt_string)//')'
call getSpareUnit(unt1,err,mess)
if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
open(unit=unt1,file=trim(DICfile),status='replace',iostat=err)
if(err>0) then;call BaM_ConsoleMessage(messID_Open,trim(DICfile));endif
! 3 DIC versions
write(unt1,trim(fmt)) lefters(1),2._mrk*Dmean-Dmaxpost
if(err>0) then;call BaM_ConsoleMessage(messID_Write,trim(DICfile));endif
write(unt1,trim(fmt)) lefters(2),Dmaxpost+Dvar
if(err>0) then;call BaM_ConsoleMessage(messID_Write,trim(DICfile));endif
write(unt1,trim(fmt)) lefters(3),Dmean+0.5_mrk*Dvar
if(err>0) then;call BaM_ConsoleMessage(messID_Write,trim(DICfile));endif
! DIC stats
write(unt1,trim(fmt)) lefters(4),Dmaxpost
if(err>0) then;call BaM_ConsoleMessage(messID_Write,trim(DICfile));endif
write(unt1,trim(fmt)) lefters(5),Dmean
if(err>0) then;call BaM_ConsoleMessage(messID_Write,trim(DICfile));endif
write(unt1,trim(fmt)) lefters(6),Dvar
if(err>0) then;call BaM_ConsoleMessage(messID_Write,trim(DICfile));endif
write(unt1,trim(fmt)) lefters(7),Dmed
if(err>0) then;call BaM_ConsoleMessage(messID_Write,trim(DICfile));endif
write(unt1,trim(fmt)) lefters(8),Ds
if(err>0) then;call BaM_ConsoleMessage(messID_Write,trim(DICfile));endif
close(unt1)

end subroutine BaM_computeDIC

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine BaM_Cleanup(err,mess)
!^**********************************************************************
!^* Purpose: performs final cleanup
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
!^* IN
!^* OUT
!^*    1.err, error code; <0:Warning, ==0:OK, >0: Error
!^*    2.mess, error message
!^**********************************************************************
use ModelLibrary_tools, only:XtraCleanup
integer(mik), intent(out)::err
character(*),intent(out)::mess
!locals
character(250),parameter::procname='BaM_Cleanup'

err=0;mess=''
call XtraCleanup(MODEL,err,mess)
end subroutine BaM_Cleanup

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine BaM_ConsoleMessage(id,mess)
!^**********************************************************************
!^* Purpose: console messages printed during BaM execution
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified:21/08/2017
!^**********************************************************************
!^* Comments: negative id's are for fatal error messages
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1. id, message id
!^*    2. mess, message that will be parsed into the printed message
!^*             - typically an error message issued by some subroutine
!^**********************************************************************

integer(mik),intent(in)::id
character(*), intent(in)::mess

select case(id)
case (1) ! Welcome!
    write(*,*) ''
    write(*,*) '*********************************'
    write(*,*) '*********************************'
    write(*,*) '**           BaM !!!           **'
    write(*,*) '*********************************'
    write(*,*) '* Bayesian Model fitting        *'
    write(*,*) '*********************************'
    write(*,*) '*********************************'
case(2)
    ! available
case(3)
    ! available
case(4)
    write(*,*) ''
    write(*,*) '*********************************'
    write(*,*) '**  MCMC sampling...           **'
    write(*,*) ''
case (5)
    write(*,*) 'Open MCMC file and monitor results...'
    write(*,*) 'MCMC File = ', trim(mess)
    write(*,*) ''
    write(*,*) '**  MCMC sampling : DONE!      **'
    write(*,*) '*********************************'
    write(*,*) ''
case(6)
    write(*,*) ''
    write(*,*) '*********************************'
    write(*,*) '**  Cooking MCMC samples...    **'
    write(*,*) ''
case(7)
    write(*,*) ''
    write(*,*) '** Cooking MCMC samples: DONE! **'
    write(*,*) '*********************************'
    write(*,*) ''
case(8)
    write(*,*) ''
    write(*,*) '*************************************'
    write(*,*) '** Summarizing MCMC samples...     **'
    write(*,*) ''
case(9)
    write(*,*) ''
    write(*,*) '** Summarizing MCMC samples: DONE! **'
    write(*,*) '*************************************'
    write(*,*) ''
case(10)
    write(*,*) ''
    write(*,*) '*********************************'
    write(*,*) '** Residual diagnostics...     **'
    write(*,*) ''
case(11)
    write(*,*) ''
    write(*,*) '** Residual diagnostics: DONE! **'
    write(*,*) '*********************************'
    write(*,*) ''
case(12)
    write(*,*) ''
    write(*,*) '***********************************'
    write(*,*) '** Prediction Experiments...     **'
    write(*,*) ''
case(13)
    write(*,*) ''
    write(*,*) '** Prediction Experiments: DONE! **'
    write(*,*) '***********************************'
    write(*,*) ''
case(14)
    write(*,*) '--------'
    write(*,*) 'Experiment number:'//trim(mess)
case(15)
    write(*,*) '~~~~~~~~'
    write(*,*) 'Processing spaghettis (transposing & computing envelops)'
case(16)
    write(*,*) 'WARNING: Computation of envelops.'
    write(*,*) 'A single spaghetti is available, no point computing envelop.'
    write(*,*) 'Computation is skipped'
case(17)
    write(*,*) 'WARNING: Prediction experiment has failed'
    write(*,*) 'Error message: '//trim(mess)
case(18)
    write(*,*) 'WARNING: Final cleanup has failed'
    write(*,*) 'Error message: '//trim(mess)
case(19)
    write(*,*) 'WARNING: No state prediction will be performed.'
case(100)
    write(*,*) ''
    write(*,*) '*********************************'
    write(*,*) '**  Reading config files..     **'
    write(*,*) ''
case(101)
    write(*,*) 'Model: '//trim(mess)
    write(*,*) '** Reading config files: DONE! **'
    write(*,*) '*********************************'
    write(*,*) ''
case(999) ! Ciao!
    write(*,*) ''
    write(*,*) '*********************************'
    write(*,*) '*********************************'
    write(*,*) '**           BaM !!!!          **'
    write(*,*) '*********************************'
    write(*,*) '*           All done!           *'
    write(*,*) '*********************************'
    write(*,*) '*********************************'
    write(*,*) ''
    ! write(*,*) 'Thanx for using me...'
    write(*,*) 'Press [enter] and I''ll go away'
case (-1) ! Fatal Error - general
    write(*,*) ""
    write(*,*) "BaM: a FATAL ERROR has occured:"
    write(*,*) trim(mess)
    write(*,*) "Execution will stop"
    write(*,*) ""
    call BaM_Fatal_Exit
case(-2) ! Fatal Error - Open File
    write(*,*) ""
    write(*,*) "BaM: a FATAL ERROR has occured"
    write(*,*) "while opening the following file."
    write(*,*) trim(mess)
    write(*,*) "Please check path to file."
    write(*,*) "Execution will stop"
    write(*,*) ""
    call BaM_Fatal_Exit
case(-3) ! Fatal Error - reading a Config file
    write(*,*) ""
    write(*,*) "BaM: a FATAL ERROR has occured"
    write(*,*) "while reading the config file:"
    write(*,*) trim(mess)
    write(*,*) "Execution will stop"
    write(*,*) ""
    call BaM_Fatal_Exit
case(-4) ! Fatal Error - Prior Model
    write(*,*) ""
    write(*,*) "BaM: a FATAL ERROR has occured"
    write(*,*) "while generating the prior model."
    write(*,*) "[error message] : "
    write(*,*) trim(mess)
    write(*,*) "Execution will stop"
    write(*,*) ""
    call BaM_Fatal_Exit
case(-5) ! Fatal Error - MCMC fit
    write(*,*) ""
    write(*,*) "BaM: a FATAL ERROR has occured"
    write(*,*) "while fitting the model."
    write(*,*) "[error message] : "
    write(*,*) trim(mess)
    write(*,*) "Execution will stop"
    write(*,*) ""
    call BaM_Fatal_Exit
case(-6) ! Fatal Error - Post-processing
    write(*,*) ""
    write(*,*) "BaM: a FATAL ERROR has occured"
    write(*,*) "while post-processing MCMC samples."
    write(*,*) "[error message] : "
    write(*,*) trim(mess)
    write(*,*) "Execution will stop"
    write(*,*) ""
    call BaM_Fatal_Exit
case(-7) ! Fatal Error - propagation
    write(*,*) ""
    write(*,*) "BaM: a FATAL ERROR has occured"
    write(*,*) "while propagating uncertainty."
    write(*,*) "[error message] : "
    write(*,*) trim(mess)
    write(*,*) "Execution will stop"
    write(*,*) ""
    call BaM_Fatal_Exit
case(-8) ! Fatal Error - writing to a file
    write(*,*) ""
    write(*,*) "BaM: a FATAL ERROR has occured"
    write(*,*) "while writting to the file:"
    write(*,*) trim(mess)
    write(*,*) "Execution will stop"
    write(*,*) ""
    call BaM_Fatal_Exit
case default
    write(*,*) ""
    write(*,*) "BaM: a FATAL ERROR has occured"
    write(*,*) "with unknown error ID..."
    write(*,*) "Execution will stop"
    write(*,*) ""
    call BaM_Fatal_Exit
end select

end subroutine BaM_ConsoleMessage


subroutine BaM_PrintHelp()
!^**********************************************************************
!^* Purpose: print help on command line arguments in console
!^**********************************************************************
!^* Programmer: Ben Renard, INRAE Aix
!^**********************************************************************
!^* Last modified:21/02/2023
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************

    write(*,'(a)') 'usage: BaM [OPTIONS]'
    write(*,'(a)') 'available options:'
    write(*,'(a)') '  -cf path, --config path:    set path to main config file'
    write(*,'(a)') '  -sd k, --seed k:            set seed to k (k should be an integer)'
    write(*,'(a)') '  -rd, --random:              randomize seed (=> non-reproducible MCMC runs)'
    write(*,'(a)') '  -v, --version:              print version information and exit'
    write(*,'(a)') '  -h, --help:                 print help and exit'
end subroutine BaM_PrintHelp


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!==============!
! Config Files !
!==============!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine Config_Read(file,&
                       workspace,&
                       Config_RunOptions,Config_Model,Config_Xtra,Config_Data,&
                       Config_RemnantSigma,Config_MCMC,&
                       Config_Cooking,Config_summary,&
                       Config_Residual,Config_Pred,&
                       err,mess)
!^**********************************************************************
!^* Purpose: Read main config file
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified: 24/07/2015
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1.file
!^* OUT
!^*    1. All parameters in Config_Model + prior lists
!^*    2.err, error code; <0:Warning, ==0:OK, >0: Error
!^*    3.mess, error message
!^**********************************************************************
use utilities_dmsl_kit,only:getSpareUnit,countSubstringInString
character(*), intent(in)::file
integer(mik), intent(out)::err
character(*),intent(out)::mess,workspace,Config_RunOptions,Config_Model,Config_Xtra,&
          Config_Data,Config_MCMC,Config_Cooking,Config_summary,&
          Config_Residual,Config_Pred
character(*),pointer::Config_RemnantSigma(:)
! locals
character(250),parameter::procname='Config_Read',sep=","
character(len_uLongStr)::foo
integer(mik)::n,unt

err=0;mess=''

call getSpareUnit(unt,err,mess)
if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
open(unit=unt,file=trim(file), status='old', iostat=err)
if(err>0) then;call BaM_ConsoleMessage(messID_Open,trim(file));endif
read(unt,*,iostat=err) workspace
if(err/=0) then;call BaM_ConsoleMessage(messID_Read,trim(file));endif
read(unt,*,iostat=err) Config_RunOptions
if(err/=0) then;call BaM_ConsoleMessage(messID_Read,trim(file));endif
read(unt,*,iostat=err) Config_Model
if(err/=0) then;call BaM_ConsoleMessage(messID_Read,trim(file));endif
read(unt,*,iostat=err) Config_Xtra
if(err/=0) then;call BaM_ConsoleMessage(messID_Read,trim(file));endif
read(unt,*,iostat=err) Config_Data
if(err/=0) then;call BaM_ConsoleMessage(messID_Read,trim(file));endif
! Remnant Sigma config files: First read to get number of files
read(unt,'(A)',iostat=err) foo
if(err/=0) then;call BaM_ConsoleMessage(messID_Read,trim(file));endif
n=countSubstringInString(foo,trim(sep))+1
if(associated(Config_RemnantSigma)) nullify(Config_RemnantSigma)
allocate(Config_RemnantSigma(n))
! Remnant Sigma config files: read
read(foo,*,iostat=err) Config_RemnantSigma
if(err/=0) then;call BaM_ConsoleMessage(messID_Read,trim(file));endif
! next config files
read(unt,*,iostat=err) Config_MCMC
if(err/=0) then;call BaM_ConsoleMessage(messID_Read,trim(file));endif
read(unt,*,iostat=err) Config_Cooking
if(err/=0) then;call BaM_ConsoleMessage(messID_Read,trim(file));endif
read(unt,*,iostat=err) Config_Summary
if(err/=0) then;call BaM_ConsoleMessage(messID_Read,trim(file));endif
read(unt,*,iostat=err) Config_Residual
if(err/=0) then;call BaM_ConsoleMessage(messID_Read,trim(file));endif
read(unt,*,iostat=err) Config_Pred
if(err/=0) then;call BaM_ConsoleMessage(messID_Read,trim(file));endif
close(unt)

end subroutine Config_Read

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine Config_Read_Model(file,&
                             ID,nX,nY,&
                             theta,&
                             err,mess)
!^**********************************************************************
!^* Purpose: Read config file for a model
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified: 06/08/2018
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1.file
!^* OUT
!^*    1. All parameters in Config_Model
!^*    2.err, error code; <0:Warning, ==0:OK, >0: Error
!^*    3.mess, error message
!^**********************************************************************
use utilities_dmsl_kit,only:getSpareUnit
character(*), intent(in)::file
integer(mik), intent(out)::err,nX,nY
Type(ParType), pointer::theta(:)
character(*),intent(out)::mess,ID
! locals
character(250),parameter::procname='Config_Read_Model'
integer(mik)::unt,i,j,ntheta

err=0;mess=''

call getSpareUnit(unt,err,mess)
if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
open(unit=unt,file=trim(file), status='old', iostat=err)
if(err>0) then;call BaM_ConsoleMessage(messID_Open,trim(file));endif
read(unt,*,iostat=err) ID
if(err/=0) then;call BaM_ConsoleMessage(messID_Read,trim(file));endif
read(unt,*,iostat=err) nX
if(err/=0) then;call BaM_ConsoleMessage(messID_Read,trim(file));endif
read(unt,*,iostat=err) nY
if(err/=0) then;call BaM_ConsoleMessage(messID_Read,trim(file));endif
read(unt,*,iostat=err) ntheta
if(err/=0) then;call BaM_ConsoleMessage(messID_Read,trim(file));endif
if(associated(theta)) nullify(theta);allocate(theta(nTheta))
do i=1,ntheta
    call Config_Read_Par(unt,theta(i),err,mess)
    if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
enddo
close(unt)
! check that all parameters have distinct names
if(ntheta>=2) then
    do i=2,ntheta
        do j=1,(i-1)
            if(trim(theta(i)%parname)==trim(theta(j)%parname)) then
                err=1;mess=trim(procname)//':'//trim(BaM_message(10));return
            endif
        enddo
    enddo
endif

end subroutine Config_Read_Model

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine Config_Read_Xtra(file,ID,xtra,err,mess)
!^**********************************************************************
!^* Purpose: Read Xtra config files for a model
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified: 21/08/2017
!^**********************************************************************
!^* Comments: Model-specific action => see ModelLibrary_tools + specific model modules
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1.file
!^*    2.ID, model ID
!^* OUT
!^*    1.xtra, extra model info
!^*    2.err, error code; <0:Warning, ==0:OK, >0: Error
!^*    3.mess, error message
!^**********************************************************************
use types_dmsl_kit, only:data_ricz_type
use ModelLibrary_tools, only:XtraRead
character(*), intent(in)::file
character(*),intent(in)::ID
type(data_ricz_type),intent(out):: xtra
integer(mik), intent(out)::err
character(*),intent(out)::mess
! locals
character(250),parameter::procname='Config_Read_Xtra'

err=0;mess=''
call XtraRead(file=file,ID=ID,xtra=xtra,err=err,mess=mess)
if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif

end subroutine Config_Read_Xtra

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine Config_Read_Data(file,&
                            DataFile,nHeader,nobs,ncol,&
                            XCol,XuCol,XbCol,XbindxCol,&
                            YCol,YuCol,YbCol,YbindxCol,&
                            err,mess)
!^**********************************************************************
!^* Purpose: Read Config_Data
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified: 22/07/2015
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1.file
!^* OUT
!^*    1. All parameters in Config_Data
!^*    2.err, error code; <0:Warning, ==0:OK, >0: Error
!^*    3.mess, error message
!^**********************************************************************
use utilities_dmsl_kit,only:getSpareUnit
character(*), intent(in)::file
integer(mik), intent(out)::err,nHeader,nobs,ncol,&
                           XCol(:),XuCol(:),XbCol(:),XbindxCol(:),&
                           YCol(:),YuCol(:),YbCol(:),YbindxCol(:)
character(*),intent(out)::mess,DataFile
! locals
character(250),parameter::procname='Config_Read_Data'
integer(mik)::unt

err=0;mess=''

call getSpareUnit(unt,err,mess)
if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
open(unit=unt,file=trim(file), status='old', iostat=err)
if(err/=0) then;call BaM_ConsoleMessage(messID_Open,trim(file));endif
read(unt,*,iostat=err) DataFile
if(err/=0) then;call BaM_ConsoleMessage(messID_Read,trim(file));endif
read(unt,*,iostat=err) nHeader
if(err/=0) then;call BaM_ConsoleMessage(messID_Read,trim(file));endif
read(unt,*,iostat=err) nobs
if(err/=0) then;call BaM_ConsoleMessage(messID_Read,trim(file));endif
read(unt,*,iostat=err) ncol
if(err/=0) then;call BaM_ConsoleMessage(messID_Read,trim(file));endif
read(unt,*,iostat=err) XCol
if(err/=0) then;call BaM_ConsoleMessage(messID_Read,trim(file));endif
read(unt,*,iostat=err) XuCol
if(err/=0) then;call BaM_ConsoleMessage(messID_Read,trim(file));endif
read(unt,*,iostat=err) XbCol
if(err/=0) then;call BaM_ConsoleMessage(messID_Read,trim(file));endif
read(unt,*,iostat=err) XbindxCol
if(err/=0) then;call BaM_ConsoleMessage(messID_Read,trim(file));endif
read(unt,*,iostat=err) YCol
if(err/=0) then;call BaM_ConsoleMessage(messID_Read,trim(file));endif
read(unt,*,iostat=err) YuCol
if(err/=0) then;call BaM_ConsoleMessage(messID_Read,trim(file));endif
read(unt,*,iostat=err) YbCol
if(err/=0) then;call BaM_ConsoleMessage(messID_Read,trim(file));endif
read(unt,*,iostat=err) YbindxCol
if(err/=0) then;call BaM_ConsoleMessage(messID_Read,trim(file));endif
close(unt)

end subroutine Config_Read_Data

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine Config_Read_RemnantSigma(files,&
                                    RemnantSigma_funk,&
                                    parname,RemnantSigma0,Prior_RemnantSigma,&
                                    err,mess)
!^**********************************************************************
!^* Purpose: Read Config_RemnantSigma
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified: 24/07/2015
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1.files (vector!)
!^* OUT
!^*    1.RemnantSigma_funk, functions for each output variable
!^*    2.parname, parameter names for each output variable
!^*    3.RemnantSigma0, starting parameter values for each output variable
!^*    4.Prior_RemnantSigma, priors for each output variable
!^*    5.err, error code; <0:Warning, ==0:OK, >0: Error
!^*    6.mess, error message
!^**********************************************************************
character(*), intent(in)::files(:)
Type(plist), intent(out):: Prior_RemnantSigma(:)
Type(parlist),intent(out):: RemnantSigma0(:)
Type(slist),intent(out):: parname(:)
integer(mik), intent(out)::err
character(*),intent(out)::mess,RemnantSigma_funk(:)
! locals
character(250),parameter::procname='Config_Read_RemnantSigma'
integer(mik)::i,n
Type(PriorListType), pointer:: Prior(:)
real(mrk),pointer:: theta0(:)
character(250), pointer::pname(:)

err=0;mess=''
n=size(files)
! size checks

do i=1,n
    call Config_Read_RemnantSigma_engine(files(i),&
                 RemnantSigma_funk(i),&
                 pname,theta0,prior,&
                 err,mess)
    if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
    if(allocated(parname(i)%s)) deallocate(parname(i)%s)
    allocate(parname(i)%s(size(theta0)))
    if(allocated(RemnantSigma0(i)%p)) deallocate(RemnantSigma0(i)%p)
    allocate(RemnantSigma0(i)%p(size(theta0)))
    if(allocated(Prior_RemnantSigma(i)%p)) deallocate(Prior_RemnantSigma(i)%p)
    allocate(Prior_RemnantSigma(i)%p(size(prior)))
    parname(i)%s=pname
    RemnantSigma0(i)%p=theta0
    Prior_RemnantSigma(i)%p=prior
enddo

end subroutine Config_Read_RemnantSigma

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine Config_Read_MCMC(file,theta0,RemnantSigma0,defaultStd,&
                            nAdapt,nCycles,InitStdMode,&
                            MinMoveRate,MaxMoveRate,DownMult,UpMult,&
                            theta_std0,RemnantSigma_std0,&
                            MCMCfile,&
                            err,mess)
!^**********************************************************************
!^* Purpose: Read Config_MCMC
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified: 20/07/2015
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1.file: config_mcmc file
!^*    2.theta0: starting point for model parameters
!^*    3.RemnantSigma0: starting point for remnant error parameters
!^*    4.defaultStd: when a starting point is zero, value given to the starting std
!^* OUT
!^*    1->15. All parameters in Config_MCMC + starting stds
!^*    16.err, error code; <0:Warning, ==0:OK, >0: Error
!^*    17.mess, error message
!^**********************************************************************
use utilities_dmsl_kit,only:getSpareUnit
character(*), intent(in)::file
real(mrk), intent(in)::theta0(:),defaultStd
type(parlist), intent(in)::RemnantSigma0(:)
integer(mik), intent(out)::err,nAdapt,nCycles,InitStdMode
real(mrk),intent(out)::MinMoveRate,MaxMoveRate,DownMult,UpMult
real(mrk),intent(out):: theta_std0(:)
type(parlist),intent(out)::RemnantSigma_std0(:) ! starting point/std for teta
character(*),intent(out)::mess,MCMCfile
! locals
character(250),parameter::procname='Config_Read_MCMC'
integer(mik)::unt,ntheta,j
real(mrk)::stdFactor

err=0;mess=''

call getSpareUnit(unt,err,mess)
if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
open(unit=unt,file=trim(file), status='old', iostat=err)
if(err>0) then;call BaM_ConsoleMessage(messID_Open,trim(file));endif
! result file
read(unt,*,iostat=err) MCMCfile
if(err/=0) then;call BaM_ConsoleMessage(messID_Read,trim(file));endif
! general properties of the sampler
read(unt,*,iostat=err) nAdapt
if(err/=0) then;call BaM_ConsoleMessage(messID_Read,trim(file));endif
read(unt,*,iostat=err) nCycles
if(err/=0) then;call BaM_ConsoleMessage(messID_Read,trim(file));endif
read(unt,*,iostat=err) MinMoveRate
if(err/=0) then;call BaM_ConsoleMessage(messID_Read,trim(file));endif
read(unt,*,iostat=err) MaxMoveRate
if(err/=0) then;call BaM_ConsoleMessage(messID_Read,trim(file));endif
read(unt,*,iostat=err) DownMult
if(err/=0) then;call BaM_ConsoleMessage(messID_Read,trim(file));endif
read(unt,*,iostat=err) UpMult
if(err/=0) then;call BaM_ConsoleMessage(messID_Read,trim(file));endif
read(unt,*,iostat=err) InitStdMode
if(err/=0) then;call BaM_ConsoleMessage(messID_Read,trim(file));endif
! Handling of initial jump stdev
ntheta=size(theta0)
if(InitStdMode==0) then ! initial stdevs from a multiplier
    read(unt,*,iostat=err) !cosmetics
    if(err/=0) then;call BaM_ConsoleMessage(messID_Read,trim(file));endif
    read(unt,*,iostat=err) stdFactor
    if(err/=0) then;call BaM_ConsoleMessage(messID_Read,trim(file));endif
    theta_std0=stdFactor*abs(theta0)
    where(theta_std0==0._mrk) theta_std0=defaultStd
    do j=1,size(RemnantSigma_std0)
        if(allocated(RemnantSigma_std0(j)%p)) deallocate(RemnantSigma_std0(j)%p)
        allocate(RemnantSigma_std0(j)%p(size(RemnantSigma0(j)%p)))
        RemnantSigma_std0(j)%p = stdFactor*abs(RemnantSigma0(j)%p)
        where(RemnantSigma_std0(j)%p==0._mrk) RemnantSigma_std0(j)%p=defaultStd
    enddo
else ! initial stdevs read from file as absolute values - botchy...
    read(unt,*,iostat=err) !cosmetics
    if(err/=0) then;call BaM_ConsoleMessage(messID_Read,trim(file));endif
    read(unt,*,iostat=err)
    if(err/=0) then;call BaM_ConsoleMessage(messID_Read,trim(file));endif
    read(unt,*,iostat=err) theta_std0
    if(err/=0) then;call BaM_ConsoleMessage(messID_Read,trim(file));endif
    do j=1,size(RemnantSigma_std0)
        if(allocated(RemnantSigma_std0(j)%p)) deallocate(RemnantSigma_std0(j)%p)
        allocate(RemnantSigma_std0(j)%p(size(RemnantSigma0(j)%p)))
        read(unt,*,iostat=err) RemnantSigma_std0(j)%p
        if(err/=0) then;call BaM_ConsoleMessage(messID_Read,trim(file));endif
    enddo
endif

close(unt)

end subroutine Config_Read_MCMC

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine Config_Read_RunOptions(file,&
                                  DoMCMC,DoSummary,&
                                  DoResidual,DoPred,&
                                  err,mess)
!^**********************************************************************
!^* Purpose: Read Config_RunOptions
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified: 04/08/2015
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1.file
!^* OUT
!^*    1. All parameters in Config_RunOptions
!^*    2.err, error code; <0:Warning, ==0:OK, >0: Error
!^*    3.mess, error message
!^**********************************************************************
use utilities_dmsl_kit,only:getSpareUnit
character(*), intent(in)::file
logical, intent(out)::DoMCMC,DoSummary,DoResidual,DoPred
integer(mik), intent(out)::err
character(*),intent(out)::mess
! locals
character(250),parameter::procname='Config_Read_RunOptions'
integer(mik)::unt
logical::exist

err=0;mess=''

call getSpareUnit(unt,err,mess)
if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
inquire(file=trim(file),exist=exist)
if(exist) then !
    open(unit=unt,file=trim(file), status='old', iostat=err)
    if(err>0) then;call BaM_ConsoleMessage(messID_Open,trim(file));endif
    read(unt,*,iostat=err) DoMCMC
    if(err/=0) then;call BaM_ConsoleMessage(messID_Read,trim(file));endif
    read(unt,*,iostat=err) DoSummary
    if(err/=0) then;call BaM_ConsoleMessage(messID_Read,trim(file));endif
    read(unt,*,iostat=err) DoResidual
    if(err/=0) then;call BaM_ConsoleMessage(messID_Read,trim(file));endif
    read(unt,*,iostat=err) DoPred
    if(err/=0) then;call BaM_ConsoleMessage(messID_Read,trim(file));endif
    close(unt)
else ! default settings - MCMC + MCMC summary + Maxpost_calib simulations
    DoMCMC=.true.
    DoSummary=.true.
    DoResidual=.true.
    DoPred=.false.
endif

end subroutine Config_Read_RunOptions

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine Config_Read_Residual(file,Residual_File,err,mess)
!^**********************************************************************
!^* Purpose: Read Config_Residual
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified: 04/08/2015
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1.file
!^* OUT
!^*    1. All parameters in Config_RunOptions
!^*    2.err, error code; <0:Warning, ==0:OK, >0: Error
!^*    3.mess, error message
!^**********************************************************************
use utilities_dmsl_kit,only:getSpareUnit
character(*), intent(in)::file
integer(mik), intent(out)::err
character(*),intent(out)::Residual_File,mess
! locals
character(250),parameter::procname='Config_Read_Residual'
integer(mik)::unt

err=0;mess=''

call getSpareUnit(unt,err,mess)
if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
open(unit=unt,file=trim(file), status='old', iostat=err)
if(err>0) then;call BaM_ConsoleMessage(messID_Open,trim(file));endif
read(unt,*,iostat=err) Residual_File
if(err/=0) then;call BaM_ConsoleMessage(messID_Read,trim(file));endif
close(unt)

end subroutine Config_Read_Residual

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine Config_Read_Cook(file,BurnFactor,nSlim,Cooking_File,err,mess)
!^**********************************************************************
!^* Purpose: Read Config_Cook
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified: 04/08/2015
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1.file
!^* OUT
!^*    1. All parameters in Config_RunOptions
!^*    2.err, error code; <0:Warning, ==0:OK, >0: Error
!^*    3.mess, error message
!^**********************************************************************
use utilities_dmsl_kit,only:getSpareUnit
character(*), intent(in)::file
real(mrk),intent(out)::BurnFactor
integer(mik), intent(out)::nSlim,err
character(*),intent(out)::Cooking_File,mess
! locals
character(250),parameter::procname='Config_Read_Cook'
integer(mik)::unt

err=0;mess=''

call getSpareUnit(unt,err,mess)
if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
open(unit=unt,file=trim(file), status='old', iostat=err)
if(err>0) then;call BaM_ConsoleMessage(messID_Open,trim(file));endif
read(unt,*,iostat=err) Cooking_File
if(err/=0) then;call BaM_ConsoleMessage(messID_Read,trim(file));endif
read(unt,*,iostat=err) BurnFactor
if(err/=0) then;call BaM_ConsoleMessage(messID_Read,trim(file));endif
read(unt,*,iostat=err) nSlim
if(err/=0) then;call BaM_ConsoleMessage(messID_Read,trim(file));endif
close(unt)

end subroutine Config_Read_Cook

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine Config_Read_Summary(file,Summary_File,DIC_file,xtendedMCMC_File,err,mess)
!^**********************************************************************
!^* Purpose: Read Config_Summary
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified: 04/08/2015
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1.file
!^* OUT
!^*    1. All parameters in Config_RunOptions
!^*    2.err, error code; <0:Warning, ==0:OK, >0: Error
!^*    3.mess, error message
!^**********************************************************************
use utilities_dmsl_kit,only:getSpareUnit
character(*), intent(in)::file
integer(mik), intent(out)::err
character(*),intent(out)::Summary_File,DIC_file,xtendedMCMC_File,mess
! locals
character(250),parameter::procname='Config_Read_Summary'
integer(mik)::unt

err=0;mess=''

call getSpareUnit(unt,err,mess)
if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
open(unit=unt,file=trim(file), status='old', iostat=err)
if(err>0) then;call BaM_ConsoleMessage(messID_Open,trim(file));endif
read(unt,*,iostat=err) Summary_File
if(err/=0) then;call BaM_ConsoleMessage(messID_Read,trim(file));endif
read(unt,*,iostat=err) DIC_File
if(err/=0) then;DIC_File='';err=0;endif ! will not do DIC computations. Ensures back-compatibility
read(unt,*,iostat=err) xtendedMCMC_File
if(err/=0) then;xtendedMCMC_File='';err=0;endif ! will not write extended MCMC. Ensures back-compatibility
close(unt)

end subroutine Config_Read_Summary

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine Config_Read_Pred_Master(file,npred,Pred_Files,err,mess)
!^**********************************************************************
!^* Purpose: Read Config_Pred_Master
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified: 09/10/2015
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1.file
!^* OUT
!^*    1. All parameters in Config_Pred_Master
!^*    2.err, error code; <0:Warning, ==0:OK, >0: Error
!^*    3.mess, error message
!^**********************************************************************
use utilities_dmsl_kit,only:getSpareUnit
character(*), intent(in)::file
integer(mik), intent(out)::npred,err
character(*), pointer::Pred_Files(:)
character(*),intent(out)::mess
! locals
character(250),parameter::procname='Config_Read_Pred_Master'
integer(mik)::i,unt

err=0;mess=''

call getSpareUnit(unt,err,mess)
if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
open(unit=unt,file=trim(file), status='old', iostat=err)
if(err>0) then;call BaM_ConsoleMessage(messID_Open,trim(file));endif
read(unt,*,iostat=err) npred
if(err/=0) then;call BaM_ConsoleMessage(messID_Read,trim(file));endif
if(associated(Pred_Files)) nullify(Pred_Files);allocate(Pred_Files(npred))
do i=1,npred
    read(unt,*,iostat=err) Pred_Files(i)
    if(err/=0) then;call BaM_ConsoleMessage(messID_Read,trim(file));endif
enddo
close(unt)

end subroutine Config_Read_Pred_Master

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine Config_Read_Pred(file,XSpag_Files,nobs,nspag,&
                            DoParametric,DoRemnant,nsim,&
                            YSpag_Files,&
                            DoTranspose,DoEnvelop,Envelop_Files,PrintCounter,&
                            DoState,SpagFiles_S,DoTranspose_S,DoEnvelop_S,EnvelopFiles_S,&
                            err,mess)
!^**********************************************************************
!^* Purpose: Read Config_Pred
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified: 09/10/2015
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1.file
!^* OUT
!^*    1. All parameters in Config_Pred_Master
!^*    2.err, error code; <0:Warning, ==0:OK, >0: Error
!^*    3.mess, error message
!^**********************************************************************
use utilities_dmsl_kit,only:getSpareUnit
character(*), intent(in)::file
integer(mik), intent(out)::nobs,nspag(:),nsim,err
logical,intent(out)::DoParametric,DoRemnant(:),DoTranspose(:),DoEnvelop(:),PrintCounter,&
                     DoState(:),DoTranspose_S(:),DoEnvelop_S(:)
character(*),intent(out)::XSpag_Files(:),YSpag_Files(:),Envelop_Files(:),mess,&
                          SpagFiles_S(:),EnvelopFiles_S(:)
! locals
character(250),parameter::procname='Config_Read_Pred'
integer(mik)::unt,nstate

err=0;mess=''

call getSpareUnit(unt,err,mess)
if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
open(unit=unt,file=trim(file), status='old', iostat=err)
if(err>0) then;call BaM_ConsoleMessage(messID_Open,trim(file));endif
read(unt,*,iostat=err) XSpag_Files
if(err/=0) then;call BaM_ConsoleMessage(messID_Read,trim(file));endif
read(unt,*,iostat=err) nobs
if(err/=0) then;call BaM_ConsoleMessage(messID_Read,trim(file));endif
read(unt,*,iostat=err) nspag
if(err/=0) then;call BaM_ConsoleMessage(messID_Read,trim(file));endif
read(unt,*,iostat=err) DoParametric
if(err/=0) then;call BaM_ConsoleMessage(messID_Read,trim(file));endif
read(unt,*,iostat=err) DoRemnant
if(err/=0) then;call BaM_ConsoleMessage(messID_Read,trim(file));endif
read(unt,*,iostat=err) Nsim
if(err/=0) then;call BaM_ConsoleMessage(messID_Read,trim(file));endif
read(unt,*,iostat=err) YSpag_Files
if(err/=0) then;call BaM_ConsoleMessage(messID_Read,trim(file));endif
read(unt,*,iostat=err) DoTranspose
if(err/=0) then;call BaM_ConsoleMessage(messID_Read,trim(file));endif
read(unt,*,iostat=err) DoEnvelop
if(err/=0) then;call BaM_ConsoleMessage(messID_Read,trim(file));endif
read(unt,*,iostat=err) Envelop_Files
if(err/=0) then;call BaM_ConsoleMessage(messID_Read,trim(file));endif
read(unt,*,iostat=err) PrintCounter
if(err/=0) then;call BaM_ConsoleMessage(messID_Read,trim(file));endif
nstate=size(DoState)
if(nstate>0) then ! read whether state prediction is required
    read(unt,*,iostat=err) DoState
    if(err/=0) then ! Do not fail, just skip state prediction with a warning
        call BaM_ConsoleMessage(19,'')
        DoState=.false.
        err=0
    endif
    if(any(Dostate)) then ! at least one state prediction, keep reading
        read(unt,*,iostat=err) SpagFiles_S
        if(err/=0) then;call BaM_ConsoleMessage(messID_Read,trim(file));endif
        read(unt,*,iostat=err) DoTranspose_S
        if(err/=0) then;call BaM_ConsoleMessage(messID_Read,trim(file));endif
        read(unt,*,iostat=err) DoEnvelop_S
        if(err/=0) then;call BaM_ConsoleMessage(messID_Read,trim(file));endif
        read(unt,*,iostat=err) EnvelopFiles_S
        if(err/=0) then;call BaM_ConsoleMessage(messID_Read,trim(file));endif
    endif
endif
close(unt)

end subroutine Config_Read_Pred

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine Config_Finalize(workspace,&
                           datafile,nrow,ncol,nHeader,&
                           theta,&
                           theta0,err,mess)
!^**********************************************************************
!^* Purpose: Finalize configuration handling VAR/STOK parameters and specifying theta0
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified: 06/08/2018
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1.worspace
!^*    2.datafile
!^*    3.nrow,ncol,nHeader: properties of datafile
!^* INOUT
!^*    1.theta
!^* OUT
!^*    1.theta0, flatten vector of initial values
!^*    2.err, error code; <0:Warning, ==0:OK, >0: Error
!^*    3.mess, error message
!^**********************************************************************
use utilities_dmsl_kit,only:getSpareUnit
character(*), intent(in)::workspace,datafile
integer(mik),intent(in)::nrow,ncol,nHeader
Type(ParType), intent(inout)::theta(:)
real(mrk), pointer::theta0(:)
integer(mik), intent(out)::err
character(*),intent(out)::mess
! locals
character(250),parameter::procname='Config_Finalize'
integer(mik)::i,j,unt,col,n,nFit,k
character(250)::file,foo
Type(ParType)::par
real(mrk)::W(nrow,ncol)

err=0;mess='';

! if at least one VAR/STOK parameter, need to read data matrix
if(count(theta(:)%parType==VAR)+count(theta(:)%parType==STOK)>0) then
    ! read data matrix
    call getSpareUnit(unt,err,mess)
    if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
    open(unit=unt,file=trim(datafile), status='old',iostat=err)
    if(err>0) then;call BaM_ConsoleMessage(messID_Open,trim(datafile));endif
    do i=1,nHeader
        read(unt,*,iostat=err)
        if(err/=0) then;call BaM_ConsoleMessage(messID_Read,trim(datafile));endif
    enddo
    do i=1,nrow
        read(unt,*,iostat=err) W(i,:)
        if(err/=0) then;call BaM_ConsoleMessage(messID_Read,trim(datafile));endif
    enddo
    close(unt)
endif

! complete parameter configuration
do i=1,size(theta)
    select case(theta(i)%partype)
    case(FIX,STD)
        if(allocated(theta(i)%indx)) deallocate(theta(i)%indx);allocate(theta(i)%indx(nrow))
        theta(i)%indx=1
    case(VAR)
        call getSpareUnit(unt,err,mess)
        if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
        file=theta(i)%config
        open(unit=unt,file=trim(workspace)//trim(file), status='old', iostat=err)
        if(err>0) then;call BaM_ConsoleMessage(messID_Open,trim(file));endif
        read(unt,*,iostat=err) theta(i)%nval
        if(err/=0) then;call BaM_ConsoleMessage(messID_Read,trim(file));endif
        read(unt,*,iostat=err) col
        if(err/=0) then;call BaM_ConsoleMessage(messID_Read,trim(file));endif
        ! read each parameter block
        if(allocated(theta(i)%init)) deallocate(theta(i)%init);allocate(theta(i)%init(theta(i)%nval))
        if(allocated(theta(i)%prior)) deallocate(theta(i)%prior);allocate(theta(i)%prior(theta(i)%nval))
        do j=1,theta(i)%nval
            call Config_Read_Par_STD(unt,foo,theta(i)%init(j),theta(i)%prior(j),err,mess)
            if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
        enddo
        ! read indices in data file
        if(allocated(theta(i)%indx)) deallocate(theta(i)%indx);allocate(theta(i)%indx(nrow))
        theta(i)%indx=W(:,col)
        close(unt)
    case(STOK)
        err=1;mess=trim(procname)//':STOK parameters not implemented yet';return
    case default
        err=1;mess=trim(procname)//":"//trim(BaM_message(11));return
    end select
    ! check consistency between nval and indx
    if(any(theta(i)%indx<=0) .or. maxval(theta(i)%indx)>theta(i)%nval) then
        err=1;mess=trim(procname)//":"//trim(BaM_message(12));return
    endif
enddo

! specify theta0
nFit=sum(theta%nval)-count(theta%partype==FIX)
if(associated(theta0)) nullify(theta0);allocate (theta0(nFit))
k=0
do i=1,size(theta)
   if(theta(i)%partype/=FIX) then
       theta0( (k+1):(k+theta(i)%nval) )=theta(i)%init
       k=k+theta(i)%nval
   endif
enddo

end subroutine Config_Finalize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine Config_Read_PriorCorr(file,exist,C,err,mess)
!^**********************************************************************
!^* Purpose: Read prior corr file, if present
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified: 09/08/2018
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1.file
!^* OUT
!^*    1.exist, does file exist?
!^*    2.C, correlation matrix
!^*    3.err, error code; <0:Warning, ==0:OK, >0: Error
!^*    4.mess, error message
!^**********************************************************************
use utilities_dmsl_kit,only:getSpareUnit
character(*), intent(in)::file
logical,intent(out)::exist
real(mrk),intent(out)::C(:,:)
integer(mik), intent(out)::err
character(*),intent(out)::mess
! locals
character(250),parameter::procname='Config_Read_PriorCorr'
integer(mik)::unt,j

err=0;mess='';C=undefRN

INQUIRE(file=trim(file),exist=exist)
if(.not.exist) then
    return
else
    ! read matrix
    call getSpareUnit(unt,err,mess)
    if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
    open(unit=unt,file=trim(file),status='old',iostat=err)
    if(err>0) then;call BaM_ConsoleMessage(messID_Open,trim(file));endif
    do j=1,size(C,dim=1)
        read(unt,*,iostat=err) C(j,:)
        if(err/=0) then;call BaM_ConsoleMessage(messID_Read,trim(file));endif
    enddo
    close(unt)
endif

end subroutine Config_Read_PriorCorr

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!==============!
! Private subs !
!==============!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine GetPriorLogPdf(theta,gamma,Xtrue,Xbias,Ybias,& ! inferred quantities
                          model,infer,& ! model and inference objects
                          prior,feas, isnull,err,mess) ! outputs
!^**********************************************************************
!^* Purpose: compute log(prior(theta,gamma,Xtrue,biases))
!^**********************************************************************
!^* Programmer: Ben Renard, INRAE Aix
!^**********************************************************************
!^* Last modified: 23/02/2023
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1. theta, model parameters
!^*    2. gamma, remnant errors parameters
!^*    3. Xtrue, estimated "true" inputs
!^*    4. Xbias, input biases
!^*    5. Ybias, output biases
!^*    6. model, model object
!^*    7. infer, infer object
!^* OUT
!^*    1.prior, prior log-pdf
!^*    2.feas, feasible?
!^*    3.is null, is (natural) posterior = zero?
!^*    4.err, error code; <0:Warning, ==0:OK, >0: Error
!^*    5.mess, error message
!^**********************************************************************
use BayesianEstimation_tools, only:GetLogPrior, PriorListType
use Distribution_tools, only: GetPdf, GAUSS

real(mrk), intent(in)::theta(:),Xtrue(:,:)
type(parlist)::gamma(:),Xbias(:),Ybias(:)
type(ModelType),intent(in):: model
type(InferenceType),intent(in):: infer
real(mrk), intent(out)::prior
logical, intent(out)::feas, isnull
integer(mik), intent(out)::err
character(*),intent(out)::mess
!locals
character(250),parameter::procname='GetPriorLogPdf'
real(mrk)::logp
real(mrk),allocatable::v(:)
integer(mik)::i,j,k

!Init
err=0;mess='';feas=.true.;isnull=.false.;prior=undefRN

!----------------------------------------------------------------
! Get Prior
prior=0._mrk
k=0
! remnant parameters
do j=1,model%nY
    call GetLogPrior(teta=gamma(j)%p,PriorList=infer%Prior_RemnantSigma(j)%p,&
                lp=logp, feas=feas, isnull=isnull,err=err,mess=mess)
    if(err>0) then;mess=trim(procname)//':'//trim(mess);feas=.false.;return;endif
    if ( (.not. feas) .or. (isnull) ) return
    prior=prior+logp
    k=k+size(gamma(j)%p)
enddo

! theta
call GetLogPrior(teta=theta,PriorList=infer%Prior_theta,&
                lp=logp, feas=feas, isnull=isnull,err=err,mess=mess)
if(err>0) then;mess=trim(procname)//':'//trim(mess);feas=.false.;return;endif
if ( (.not. feas) .or. (isnull) ) return
prior=prior+logp
k=k+size(theta)

! prior correlation
if(allocated(infer%priorM)) then ! used has specified a prior correlation => Gaussian copula computation
    ! compute gaussian-transformed values v's
    if(allocated(v)) deallocate(v);allocate(v(k))
    k=0
    do i=1,size(theta)
        k=k+1
        call Gcop_getV(x=theta(i),distID=infer%Prior_theta(i)%dist,&
                       par=infer%Prior_theta(i)%par,&
                       v=v(k),feas=feas,err=err,mess=mess)
        if(err>0) then;mess=trim(procname)//':'//trim(mess);feas=.false.;return;endif
        if (.not. feas) return
    enddo
    do j=1,model%nY
        do i=1,size(gamma(j)%p)
            k=k+1
            call Gcop_getV(x=gamma(j)%p(i),distID=infer%Prior_RemnantSigma(j)%p(i)%dist,&
                           par=infer%Prior_RemnantSigma(j)%p(i)%par,&
                           v=v(k),feas=feas,err=err,mess=mess)
            if(err>0) then;mess=trim(procname)//':'//trim(mess);feas=.false.;return;endif
            if (.not. feas) return
        enddo
    enddo

    ! compute Gaussian copula term
    ! NOTE: the correct computation encompasses an additional term -0.5_mrk*infer%priorLogDet
    ! This term is not here because there seems to be a bug in DMSL's computation of the
    ! log-determinant within sub choles_invrt
    ! (more accurately: in choles_dcmp_engine1, if uninitialized variable posDefinite is set to
    ! false by the compiler, which seems to be the case with IVF, then computation does not proceed)
    ! The log-determinant is just a constant wrt to infer parameters, so this omission is not a problem
    ! in the context of Bayesian - MCMC estimation. But still, need to check with Dmitri.
    prior=prior-0.5_mrk*dot_product(matmul(v,infer%priorM),v)
endif


! input biases
do j=1,model%nX
    do i=1,infer%nXb(j)
        call GetPdf(DistId=GAUSS,x=Xbias(j)%p(i),par=(/0._mrk,1._mrk/),loga=.true.,&
                    pdf=logp,feas=feas,isnull=isnull,err=err,mess=mess)
        if(err>0) then;mess=trim(procname)//':'//trim(mess);feas=.false.;return;endif
        if ( (.not. feas) .or. (isnull) ) return
        prior=prior+logp
    enddo
enddo

! output biases
do j=1,model%nY
    do i=1,infer%nYb(j)
        call GetPdf(DistId=GAUSS,x=Ybias(j)%p(i),par=(/0._mrk,1._mrk/),loga=.true.,&
                    pdf=logp,feas=feas,isnull=isnull,err=err,mess=mess)
        if(err>0) then;mess=trim(procname)//':'//trim(mess);feas=.false.;return;endif
        if ( (.not. feas) .or. (isnull) ) return
        prior=prior+logp
    enddo
enddo
!----------------------------------------------------------------

end subroutine GetPriorLogPdf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine GetLogLikelihood(theta,gamma,Xtrue,Xbias,Ybias,& ! inferred quantities
                            model,infer,& ! model and inference objects
                            lkh,Dpar,feas, isnull,err,mess) ! outputs
!^**********************************************************************
!^* Purpose: compute log(p(data|theta,gamma,Xtrue,biases))
!^**********************************************************************
!^* Programmer: Ben Renard, INRAE Aix
!^**********************************************************************
!^* Last modified: 23/02/2023
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1. theta, model parameters
!^*    2. gamma, remnant errors parameters
!^*    3. Xtrue, estimated "true" inputs
!^*    4. Xbias, input biases
!^*    5. Ybias, output biases
!^*    6. model, model object
!^*    7. infer, infer object
!^* OUT
!^*    1.lkh, log-likelihood
!^*    2.Dpar, derived parameters
!^*    3.feas, feasible?
!^*    4.is null, is (natural) posterior = zero?
!^*    5.err, error code; <0:Warning, ==0:OK, >0: Error
!^*    6.mess, error message
!^**********************************************************************
use Distribution_tools, only: GetPdf, GAUSS

real(mrk), intent(in)::theta(:),Xtrue(:,:)
type(parlist)::gamma(:),Xbias(:),Ybias(:)
type(ModelType),intent(in):: model
type(InferenceType),intent(in):: infer
real(mrk), intent(out)::lkh
real(mrk), intent(out)::Dpar(:)
logical, intent(out)::feas, isnull
integer(mik), intent(out)::err
character(*),intent(out)::mess
!locals
character(250),parameter::procname='GetLogLikelihood'
real(mrk)::logp,mu,sig
real(mrk)::yhat(model%nY),bigY(infer%nObs,model%nY),&
           state(infer%nObs,model%nstate)
integer(mik)::n,i,j

!Init
err=0;mess='';feas=.true.;isnull=.false.;lkh=undefRN
n=infer%nObs

!----------------------------------------------------------------
lkh=0._mrk
! Run model for all data
call ApplyModel_BaM(model=model,X=Xtrue,theta=theta,Y=bigY,&
                Dpar=Dpar,state=state,feas=feas,err=err,mess=mess)
if(err>0) then;mess=trim(procname)//':'//trim(mess);return;endif
if(.not. feas) return
! Add up likelihood contributions
do i=1,n
    Yhat=bigY(i,:)
    do j=1,model%nY
        if(infer%Y(i,j)==infer%mv) cycle
        call Sigmafunk_Apply(infer%RemnantSigma_funk(j), gamma(j)%p, Yhat(j), sig, err, mess)
        if(err>0) then;mess=trim(procname)//':'//trim(mess);feas=.false.;return;endif
        if(sig<=0._mrk) then;
            feas=.false.;return
        else
            sig=sqrt(infer%Yu(i,j)**2+sig**2)
        endif
        ! Yobs ~ N(mu = Yhat+b*Bias, sigma2 = sigma_Y2 + sigma_remnant2)
        if(infer%Ybindx(i,j)==0) then
            mu=Yhat(j)
        else
            mu=Yhat(j)+infer%Yb(i,j)*Ybias(j)%p(infer%Ybindx(i,j))
        endif
        call GetPdf(DistId=GAUSS,x=infer%Y(i,j),par=(/mu,sig/),loga=.true.,pdf=logp,&
                feas=feas,isnull=isnull,err=err,mess=mess)
        if(err>0) then;mess=trim(procname)//':'//trim(mess);return;endif
        if( (.not. feas) .or. isnull) return
        lkh=lkh+logp
    enddo
enddo

end subroutine GetLogLikelihood

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine GetHyperLogPdf(theta,gamma,Xtrue,Xbias,Ybias,& ! inferred quantities
                         model,infer,& ! model and inference objects
                         hyper,feas,isnull,err,mess) ! outputs
!^**********************************************************************
!^* Purpose: compute log(p(Xobs|Xtrue,Xbias))
!^**********************************************************************
!^* Programmer: Ben Renard, INRAE Aix
!^**********************************************************************
!^* Last modified: 23/02/2023
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1. theta, model parameters
!^*    2. gamma, remnant errors parameters
!^*    3. Xtrue, estimated "true" inputs
!^*    4. Xbias, input biases
!^*    5. Ybias, output biases
!^*    6. model, model object
!^*    7. infer, infer object
!^* OUT
!^*    1.hyper, hyperdistribution log-pdf
!^*    2.feas, feasible?
!^*    3.is null, is (natural) posterior = zero?
!^*    4.err, error code; <0:Warning, ==0:OK, >0: Error
!^*    5.mess, error message
!^**********************************************************************
use Distribution_tools, only: GetPdf, GAUSS

real(mrk), intent(in)::theta(:),Xtrue(:,:)
type(parlist)::gamma(:),Xbias(:),Ybias(:)
type(ModelType),intent(in):: model
type(InferenceType),intent(in):: infer
real(mrk), intent(out)::hyper
logical, intent(out)::feas, isnull
integer(mik), intent(out)::err
character(*),intent(out)::mess
!locals
character(250),parameter::procname='GetHyperLogPdf'
real(mrk)::logp,mu
integer(mik)::n,i,j

!Init
err=0;mess='';feas=.true.;isnull=.false.;hyper=undefRN
n=infer%nObs

!----------------------------------------------------------------
hyper=0._mrk
do j=1,model%nX
    if(.not.infer%Xerror(j)) cycle
    do i=1,n
        if(infer%Xu(i,j) == 0._mrk) cycle
        ! get hypermean
        if(infer%Xbindx(i,j)==0) then
            mu=Xtrue(i,j)
        else
            mu=Xtrue(i,j)+infer%Xb(i,j)*Xbias(j)%p(infer%Xbindx(i,j))
        endif
        ! get hyper-pdf
        call GetPdf(DistId=GAUSS,x=infer%X(i,j),&
                    par=(/mu,infer%Xu(i,j)/),&
                    loga=.true.,pdf=logp,feas=feas,isnull=isnull,err=err,mess=mess)
        if(err>0) then;mess=trim(procname)//':'//trim(mess);return;endif
        if( (.not. feas) .or. isnull) return
        hyper=hyper+logp
    enddo
enddo
!----------------------------------------------------------------

end subroutine GetHyperLogPdf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine GetPostLogPdf(theta,gamma,Xtrue,Xbias,Ybias,& ! inferred quantities
                         model,infer,& ! model and inference objects
                         lp,Dpar,feas, isnull,err,mess) ! outputs
!^**********************************************************************
!^* Purpose: compute log(p(theta,gamma,Xtrue,biases|data)) [unnormalized]
!^**********************************************************************
!^* Programmer: Ben Renard, INRAE Aix
!^**********************************************************************
!^* Last modified: 23/02/2023
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1. theta, model parameters
!^*    2. gamma, remnant errors parameters
!^*    3. Xtrue, estimated "true" inputs
!^*    4. Xbias, input biases
!^*    5. Ybias, output biases
!^*    6. model, model object
!^*    7. infer, infer object
!^* OUT
!^*    1.lp, unnormalized posterior log-pdf
!^*    2.Dpar, derived parameters
!^*    3.feas, feasible?
!^*    4.is null, is (natural) posterior = zero?
!^*    5.err, error code; <0:Warning, ==0:OK, >0: Error
!^*    6.mess, error message
!^**********************************************************************
real(mrk), intent(in)::theta(:),Xtrue(:,:)
type(parlist)::gamma(:),Xbias(:),Ybias(:)
type(ModelType),intent(in):: model
type(InferenceType),intent(in):: infer
real(mrk), intent(out)::lp
real(mrk), intent(out)::Dpar(:)
logical, intent(out)::feas, isnull
integer(mik), intent(out)::err
character(*),intent(out)::mess
!locals
character(250),parameter::procname='GetPostLogPdf'
real(mrk)::prior,lkh, hyper

!Init
err=0;mess='';feas=.true.;isnull=.false.;lp=undefRN
!----------------------------------------------------------------
! Get Prior
call GetPriorLogPdf(theta=theta,gamma=gamma,Xtrue=Xtrue,Xbias=Xbias,Ybias=Ybias,& ! inferred quantities
                    model=model,infer=infer,& ! model and inference objects
                    prior=prior,feas=feas,isnull=isnull,err=err,mess=mess) ! outputs
if(err>0) then;mess=trim(procname)//':'//trim(mess);feas=.false.;return;endif
if ( (.not. feas) .or. (isnull) ) return
!----------------------------------------------------------------
! Get Likelihood
call GetLogLikelihood(theta=theta,gamma=gamma,Xtrue=Xtrue,Xbias=Xbias,Ybias=Ybias,& ! inferred quantities
                      model=model,infer=infer,& ! model and inference objects
                      lkh=lkh,Dpar=Dpar,feas=feas,isnull=isnull,err=err,mess=mess) ! outputs
if(err>0) then;mess=trim(procname)//':'//trim(mess);feas=.false.;return;endif
if ( (.not. feas) .or. (isnull) ) return
!----------------------------------------------------------------
! Get Hyper
call GetHyperLogPdf(theta=theta,gamma=gamma,Xtrue=Xtrue,Xbias=Xbias,Ybias=Ybias,& ! inferred quantities
                    model=model,infer=infer,& ! model and inference objects
                    hyper=hyper,feas=feas,isnull=isnull,err=err,mess=mess) ! outputs
if(err>0) then;mess=trim(procname)//':'//trim(mess);feas=.false.;return;endif
if ( (.not. feas) .or. (isnull) ) return
!----------------------------------------------------------------
! Get Posterior
lp = prior + lkh + hyper

end subroutine GetPostLogPdf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine UnfoldParVector(v,theta,gamma,Xtrue,Xbias,Ybias,err,mess)
!^**********************************************************************
!^* Purpose: Interpret a parameter vector and split it accordingly
!^**********************************************************************
!^* Programmer: Ben Renard, INRAE
!^**********************************************************************
!^* Last modified:23/02/2023
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1.v, parameter vector
!^* OUT
!^*    1.theta, model parameters
!^*    2.gamma, remnant errors parameters
!^*    3.Xtrue, true inputs if inferred
!^*    4.Xbias, input biases if any
!^*    5.Ybias, output biases if any
!^*    6.err, error code; <0:Warning, ==0:OK, >0: Error
!^*    7.mess, error message
!^**********************************************************************
real(mrk),intent(in)::v(:)
real(mrk),intent(out)::theta(INFER%nFit),Xtrue(INFER%nobs,MODEL%nX)
type(parlist),intent(out)::gamma(MODEL%nY),Xbias(MODEL%nX),Ybias(MODEL%nY)
integer(mik),intent(out)::err
character(*),intent(out)::mess
!locals
character(250),parameter::procname='UnfoldParVector'
integer(mrk)::k,i,j

! theta
if(INFER%nFit>0) theta=v(1:(INFER%nFit))
k=INFER%nFit
! remnant parameters
do i=1,MODEL%nY
    if(allocated(gamma(i)%p)) deallocate(gamma(i)%p);allocate(gamma(i)%p(INFER%nremnant(i)))
    gamma(i)%p=v( (k+1):(k+INFER%nremnant(i)) )
    k=k+INFER%nremnant(i)
enddo
! inferred true inputs
Xtrue=INFER%X
do j=1,MODEL%nX
    do i=1,INFER%nobs
        if(INFER%Xu(i,j) > 0._mrk) then
            Xtrue(i,j)=v(k+1)
            k=k+1
        endif
    enddo
enddo
! input biases
do j=1,MODEL%nX
    ! allocate biases
    if(allocated(Xbias(j)%p)) deallocate(Xbias(j)%p);allocate(Xbias(j)%p(INFER%nXb(j)))
    if(INFER%nXb(j)==0) cycle ! no bias
    Xbias(j)%p=v( (k+1):(k+INFER%nXb(j)) )
    ! additional loop to determine whether Xtrue is estimated or not
    do i=1,INFER%nobs
        if(INFER%Xbindx(i,j)>0) then ! bias on this time step
            if(INFER%Xu(i,j) == 0._mrk) then ! Xtrue is  NOT estimated - need to "unbias" Xobs
                Xtrue(i,j)=INFER%X(i,j)- INFER%Xb(i,j)*Xbias(j)%p(infer%Xbindx(i,j))
            endif
        endif
    enddo
    k=k+INFER%nXb(j)
enddo
! output biases
do j=1,MODEL%nY
    if(allocated(Ybias(j)%p)) deallocate(Ybias(j)%p);allocate(Ybias(j)%p(INFER%nYb(j)))
    if(INFER%nYb(j)==0) cycle ! no bias
    Ybias(j)%p=v( (k+1):(k+INFER%nYb(j)) )
    k=k+INFER%nXb(j)
enddo

end subroutine UnfoldParVector

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine Posterior_wrapper(x,feas,isnull,fx,fAux,err,mess)
!^**********************************************************************
!^* Purpose: wrapper to GetPostLogPdf to comply with MCMC interface
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea
!^**********************************************************************
!^* Last modified:10/07/2015
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1.x, parameters
!^* OUT
!^*    1.feas, feasible?
!^*    2.isnull, is post=0?
!^*    3.fx, posterior(par|data)
!^*    4.[fAux], auxiliaries
!^*    5.err, error code; <0:Warning, ==0:OK, >0: Error
!^*    6.mess, error message
!^**********************************************************************

real(mrk),intent(in)::x(:)
logical,intent(out)::feas,isnull
real(mrk),intent(out)::fx
real(mrk),intent(out),optional::fAux(:)
integer(mik),intent(out)::err
character(*),intent(out)::mess
!locals
character(250),parameter::procname='Posterior_wrapper'
real(mrk):: theta(INFER%nFit),Xtrue(INFER%nobs,MODEL%nX),Dpar(INFER%nDpar)
type(parlist)::gamma(MODEL%nY),Xbias(MODEL%nX),Ybias(MODEL%nY)
integer(mrk)::k,i,j

!---------------------------
! unfold x
call UnfoldParVector(v=x,theta=theta,gamma=gamma,Xtrue=Xtrue,Xbias=Xbias,Ybias=Ybias,err=err,mess=mess)
if(err>0) then;mess=trim(procname)//':'//trim(mess);feas=.false.;return;endif

!---------------------------
! Compute log-posterior
call GetPostLogPdf(theta=theta,gamma=gamma,Xtrue=Xtrue,Xbias=Xbias,Ybias=Ybias,&
                   model=MODEL,infer=INFER,&
                   lp=fx,Dpar=DPar,feas=feas, isnull=isnull,err=err,mess=mess)
if(err>0) then;mess=trim(procname)//':'//trim(mess);feas=.false.;return;endif
if(present(fAux)) fAux=Dpar

end subroutine Posterior_wrapper

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

pure subroutine Sigmafunk_Apply(funk,par,Y,res,err,mess)
!^**********************************************************************
!^* Purpose: Apply the selected Sigmafunk
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Created: 29/04/2013, last modified: 05/08/2022, added 'Power'
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1.funk, which function?? (e.g., 'Constant','Linear')
!^*    2.par, parameters of funk
!^*    3.Y, covariate of funk
!^* OUT
!^*    1.res, result
!^*    2.err, error code; <0:Warning, ==0:OK, >0: Error
!^*    3.mess, error message
!^**********************************************************************
character(*), intent(in)::funk
real(mrk), intent(in)::par(:),Y
real(mrk), intent(out)::res
integer(mik), intent(out)::err
character(*),intent(out)::mess
! locals
character(250),parameter::procname='Sigmafunk_Apply'

err=0;mess='';res=undefRN

select case(trim(funk))
case('Constant')
    res=par(1)
case('Linear')
    res=par(1)+par(2)*abs(Y)
case('Proportional')
    res=par(1)*abs(Y)
case('Power')
    res=par(1)+par(2)*(abs(Y)**par(3))
case('Exponential')
    res=par(1)+( par(3)-par(1) )*( 1._mrk - exp(- (abs(Y)/par(2))**1 ) )
case('Gaussian')
    res=par(1)+( par(3)-par(1) )*( 1._mrk - exp(- (abs(Y)/par(2))**2 ) )
case default
    err=1;mess=trim(procname)//":"//trim(BaM_message(2))
end select

end subroutine SigmaFunk_Apply

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

pure subroutine Sigmafunk_GetParNumber(funk, npar, err, mess)
!^**********************************************************************
!^* Purpose: Get number of parameters of the selected Sigmafunk
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Created: 29/04/2013, last modified: 05/08/2022, added 'Power'
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1.funk, which function?? (e.g., 'Linear')
!^* OUT
!^*    1.npar
!^*    2.err, error code; <0:Warning, ==0:OK, >0: Error
!^*    3.mess, error message
!^**********************************************************************
character(*), intent(in)::funk
integer(mik), intent(out)::npar
integer(mik), intent(out)::err
character(*),intent(out)::mess
!locals
character(250),parameter::procname='Sigmafunk_GetParNumber'

err=0;mess='';npar=undefIN
select case(trim(funk))
case('Constant','Proportional')
    npar=1
case('Linear')
    npar=2
case('Exponential', 'Gaussian', 'Power')
    npar=3
case default
    err=1;mess=trim(procname)//':'//trim(BaM_message(2))
end select

end subroutine Sigmafunk_GetParNumber

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine GetSim_calibration(par,Xtrue,Ysim,state,feas,err,mess)
!^**********************************************************************
!^* Purpose: Runs the model on calibration data
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified: 03/08/2015
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1.par, parameter vector
!^* OUT
!^*    1.Xtrue, simulated true inputs
!^*    2.Ysim, simulated values
!^*    3.state, model states
!^*    4.feas, feasible?
!^*    5.err, error code; <0:Warning, ==0:OK, >0: Error
!^*    6.mess, error message
!^**********************************************************************
real(mrk),intent(in)::par(:)
real(mrk),intent(out)::Xtrue(:,:),Ysim(:,:),state(:,:)
logical, intent(out)::feas
integer(mik), intent(out)::err
character(*),intent(out)::mess
!locals
character(250),parameter::procname='GetSim_calibration'
integer(mik)::i,j,k
real(mrk),allocatable::Dpar(:)
real(mrk)::theta(INFER%nFit)

err=0;mess='';Ysim=undefRN
! Get true inputs
Xtrue=INFER%X
k=INFER%nFit+sum(INFER%nremnant)
do j=1,MODEL%nX
    do i=1,INFER%nobs
        if(INFER%Xu(i,j) > 0._mrk) then
            Xtrue(i,j)=par(k+1)
            k=k+1
        endif
    enddo
enddo
! Xtra loop for unbiasing those input that are not inferred
do j=1,MODEL%nX
    if(INFER%nXb(j)==0) cycle ! no bias
    do i=1,INFER%nobs
        if(INFER%Xbindx(i,j)>0) then ! bias on this time step
            if(INFER%Xu(i,j) == 0._mrk) then ! Xtrue is  NOT estimated - need to "unbias" Xobs
                Xtrue(i,j)=INFER%X(i,j)- INFER%Xb(i,j)*par(k+infer%Xbindx(i,j))
            endif
        endif
    enddo
    k=k+INFER%nXb(j)
enddo

! Run model
allocate(Dpar(INFER%nDpar))
if(INFER%nFit>0) theta=par(1:INFER%nFit)
call ApplyModel_BaM(model=model,X=Xtrue,theta=theta,Y=Ysim,&
                Dpar=Dpar,state=state,feas=feas,err=err,mess=mess)
if(err>0) then;mess=trim(procname)//':'//trim(mess);return;endif
if(.not. feas) return

end subroutine GetSim_calibration

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

pure function BaM_message(id)
! returns error/warning messages
integer(mik),intent(in)::id
character(250)::BaM_message

select case(id)
case(1)
    BaM_message='size mismatch in input or output data'
case(2)
    BaM_message='unknown SigmaFunk'
case(3)
    BaM_message='size mismatch in parameters'
case(4)
    BaM_Message='incorrect number of parameters for remnant sigma function'
case(5)
    BaM_Message='problem while performing maxpost runs (calibration)'
case(6)
    BaM_Message='unknown model%ID'
case(7)
    BaM_Message='[nsim] compulsory when [mcmc] is not provided'
case(8)
    BaM_Message='prior mode is not feasible'
case(9)
    BaM_Message='impossible to generate from prior distribution'
case(10)
    BaM_Message='several parameters have identical names'
case(11)
    BaM_Message='unknown parameter type'
case(12)
    BaM_Message='problem with par%indx (all values in par%indx should be between 1 and par%nval)'
case(13)
    BaM_Message='Prediction with VAR or STOK parameters is not allowed - prediction aborted'
case(14)
    BaM_Message='Prior correlation matrix: all diagonal terms should be equal to 1'
case(15)
    BaM_Message='Prior correlation matrix: not positive definite'
case(16)
    BaM_Message='Parameter vector leads to zero or unfeasible prior'
case(17)
    BaM_Message='Parameter vector leads to zero or unfeasible hyper-distribution'
case(18)
    BaM_Message='FIX, VAR or STOK parameters are not allowed for remnant error parameters'
case default
    BaM_message='unknown message ID'
end select

end function BaM_message

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine BaM_Fatal_Exit
! action taken on fatal error
!read(*,*)
STOP
end subroutine BaM_Fatal_Exit

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine Config_Read_Par(unt,par,err,mess)
!^**********************************************************************
!^* Purpose: Read a parameter bloc in a config file
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified: 23/07/2015
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1.unt: unit of the file (already opened, not close within this sub)
!^* OUT
!^*    1.parname: parameter name
!^*    2.theta0: initial parameter value
!^*    3.prior: prior distribution
!^*    4.err, error code; <0:Warning, ==0:OK, >0: Error
!^*    5.mess, error message
!^**********************************************************************
use Distribution_tools, only:GetParNumber
integer(mik),intent(in)::unt
type(ParType),intent(out)::par
integer(mik), intent(out)::err
character(*),intent(out)::mess
!locals
character(250),parameter::procname='Config_Read_Par'
integer(mik)::np
real(mrk)::init
character(250)::priordist

read(unt,*,iostat=err) par%parname
if(err/=0) then;mess=trim(procname)//':ReadError';return;endif
read(unt,*,iostat=err) init
if(err/=0) then;mess=trim(procname)//':ReadError';return;endif
read(unt,*,iostat=err) priordist
if(err/=0) then;mess=trim(procname)//':ReadError';return;endif
select case (trim(priordist))
case(FIX_str)
    par%partype=FIX
    par%nval=1
    if(allocated(par%init)) deallocate(par%init);allocate(par%init(par%nval))
    par%init=init
    read(unt,*,iostat=err) ! read nothing
    if(err/=0) then;mess=trim(procname)//':ReadError';return;endif
    if(allocated(par%prior)) deallocate(par%prior);allocate(par%prior(par%nval))
    par%prior(1)%dist=trim(priordist)
    if(allocated(par%prior(1)%par)) deallocate(par%prior(1)%par);allocate(par%prior(1)%par(0))
case(VAR_str)
    par%partype=VAR
    read(unt,*,iostat=err) par%config ! read config file
    if(err/=0) then;mess=trim(procname)//':ReadError';return;endif
case(STOK_str)
    err=1;mess=trim(procname)//':STOK parameters not implemented yet';return
case default
    par%partype=STD
    par%nval=1
    if(allocated(par%init)) deallocate(par%init);allocate(par%init(par%nval))
    par%init=init
    if(allocated(par%prior)) deallocate(par%prior);allocate(par%prior(par%nval))
    par%prior(1)%dist=trim(priordist)
    call GetParNumber(DistID=trim(par%prior(1)%dist), npar=np, err=err, mess=mess)
    if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
    if(allocated(par%prior(1)%par)) deallocate(par%prior(1)%par);allocate(par%prior(1)%par(np))
    read(unt,*,iostat=err) par%prior(1)%par(:) ! read prior pars
    if(err/=0) then;mess=trim(procname)//':ReadError';return;endif
endselect
end subroutine Config_Read_Par

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine Config_Read_Par_STD(unt,parname,theta0,prior,err,mess)
!^**********************************************************************
!^* Purpose: Read a standard parameter bloc in a config file
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified: 23/07/2015
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1.unt: unit of the file (already opened, not close within this sub)
!^* OUT
!^*    1.parname: parameter name
!^*    2.theta0: initial parameter value
!^*    3.prior: prior distribution
!^*    4.err, error code; <0:Warning, ==0:OK, >0: Error
!^*    5.mess, error message
!^**********************************************************************
use Distribution_tools, only:GetParNumber
integer(mik),intent(in)::unt
character(*), intent(out)::parname
real(mrk),intent(out)::theta0
type(PriorListType),intent(out)::prior
integer(mik), intent(out)::err
character(*),intent(out)::mess
!locals
character(250),parameter::procname='Config_Read_Par_STD'
integer(mik)::np

read(unt,*,iostat=err) parname
if(err/=0) then;mess=trim(procname)//':ReadError';return;endif
read(unt,*,iostat=err) theta0
if(err/=0) then;mess=trim(procname)//':ReadError';return;endif
read(unt,*,iostat=err) prior%dist
if(err/=0) then;mess=trim(procname)//':ReadError';return;endif
prior%dist=trim(prior%dist)
if(trim(prior%dist)==FIX_str .or. trim(prior%dist)==VAR_str .or. trim(prior%dist)==STOK_str) then
    err=1;mess=trim(procname)//':'//trim(BaM_Message(18));return
else
    call GetParNumber(DistID=prior%dist, npar=np, err=err, mess=mess)
    if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
    if(allocated(prior%par)) deallocate(prior%par)
    allocate(prior%par(np))
    read(unt,*,iostat=err) prior%par(:)
    if(err/=0) then;mess=trim(procname)//':ReadError';return;endif
endif
end subroutine Config_Read_Par_STD

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine Config_Read_RemnantSigma_engine(file,&
                                           RemnantSigma_funk,&
                                           parname,RemnantSigma0,Prior_RemnantSigma,&
                                           err,mess)
!^**********************************************************************
!^* Purpose: Read Config_RemnantSigma (for a single output variable)
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified: 24/07/2015
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1.file
!^* OUT
!^*    1.RemnantSigma_funk
!^*    2.parname, parameter names
!^*    3.RemnantSigma0, starting parameter values
!^*    4.Prior_RemnantSigma, priors
!^*    5.err, error code; <0:Warning, ==0:OK, >0: Error
!^*    6.mess, error message
!^**********************************************************************
use utilities_dmsl_kit,only:getSpareUnit
use Distribution_tools, only:GetParNumber
character(*), intent(in)::file
Type(PriorListType), pointer:: Prior_RemnantSigma(:)
character(*),pointer:: parname(:)
real(mrk),pointer:: RemnantSigma0(:)
integer(mik), intent(out)::err
character(*),intent(out)::mess,RemnantSigma_funk
! locals
character(250),parameter::procname='Config_Read_RemnantSigma_engine'
integer(mik)::unt,i,np,np0,n
type(ParType)::par

err=0;mess=''

call getSpareUnit(unt,err,mess)
if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
open(unit=unt,file=trim(file), status='old', iostat=err)
if(err>0) then;call BaM_ConsoleMessage(messID_Open,trim(file));endif
read(unt,*,iostat=err) RemnantSigma_funk
if(err>0) then;call BaM_ConsoleMessage(messID_Read,trim(file));endif
read(unt,*,iostat=err) np
if(err>0) then;call BaM_ConsoleMessage(messID_Read,trim(file));endif
! check number of parameters is correct
call Sigmafunk_GetParNumber(funk=RemnantSigma_funk, npar=np0, err=err, mess=mess)
if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
if(np0/=np) then
    err=1;mess=trim(procname)//':'//BaM_message(4);return
endif
if(associated(parname)) nullify(parname);allocate(parname(np))
if(associated(RemnantSigma0)) nullify(RemnantSigma0);allocate(RemnantSigma0(np))
if(associated(Prior_RemnantSigma)) nullify(Prior_RemnantSigma);allocate(Prior_RemnantSigma(np))
!read each parameter block
do i=1,np
    call Config_Read_Par_STD(unt,parname(i),RemnantSigma0(i),Prior_RemnantSigma(i),err,mess)
    if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
enddo
close(unt)

end subroutine Config_Read_RemnantSigma_engine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine GetHeaders(head)
! get headers from MODEL object
use utilities_dmsl_kit,only:number_string
character(*),intent(out)::head(:)
! locals
integer(mik)::i,j,k

k=0
! theta
do i=1,MODEL%ntheta
    if(INFER%theta(i)%partype==STD) then
        k=k+1;head(k)=MODEL%parname(i)
    else if(INFER%theta(i)%partype==VAR) then
        do j=1,INFER%theta(i)%nval
            k=k+1;head(k)=trim(MODEL%parname(i))//"_"//trim(number_string(j))
        enddo
    endif
enddo
! remnant parameters
do i=1,MODEL%nY
    do j=1,INFER%nremnant(i)
       head(k+j)="Y"//trim(number_string(i))//"_"//trim(INFER%Parname_RemnantSigma(i)%s(j))
    enddo
    k=k+INFER%nremnant(i)
enddo
! inferred true inputs
do j=1,MODEL%nX
    do i=1,INFER%nObs
        if(INFER%Xu(i,j) > 0._mrk) then
            k=k+1
            head(k)="X"//trim(number_string(j))//"_TrueVal"//trim(number_string(i))
        endif
    enddo
enddo
! Input biases
do j=1,MODEL%nX
    do i=1,INFER%nXb(j)
        k=k+1
        head(k)="X"//trim(number_string(j))//"_Bias"//trim(number_string(i))
    enddo
enddo
! Output biases
do j=1,MODEL%nY
    do i=1,INFER%nYb(j)
        k=k+1
        head(k)="Y"//trim(number_string(j))//"_Bias"//trim(number_string(i))
    enddo
enddo
! Log-post
k=k+1;head(k)="LogPost"
! Derived parameters
do i=1,INFER%nDpar
    k=k+1;head(k)=INFER%DparName(i)
enddo
end subroutine GetHeaders

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine BaM_getPriorMC(MC,mode,err,mess)
!^**********************************************************************
!^* Purpose: Get prior Monte Carlo simulations and prior mode
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified: 09/10/2015
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^* OUT
!^*    1.[MC], Monte Carlo simulations
!^*    2.mode, prior mode
!^*    3.err, error code; <0:Warning, ==0:OK, >0: Error
!^*    4.mess, error message
!^**********************************************************************
use Distribution_tools, only:GetMode,generateSample
real(mrk), intent(out),optional::MC(:,:)
real(mrk), intent(out)::mode(:)
integer(mik), intent(out)::err
character(*),intent(out)::mess
!locals
character(250),parameter::procname='BaM_getPriorMC'
integer(mik)::k,i,j,compt
logical::feas

mode=undefRN;err=0;mess='';if(present(MC)) MC=undefRN;
! model parameters
do k=1,INFER%nFit
    ! get mode
    call GetMode(DistId=INFER%Prior_theta(k)%dist,&
                 par=INFER%Prior_theta(k)%par,&
                 m=mode(k),feas=feas,err=err,mess=mess)
    if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
    if(.not.feas) then
        err=1;mess=trim(procname)//':'//trim(BaM_Message(8));return
    endif
    ! Do MC simulations
    if(present(MC)) then
        call generateSample(DistId=INFER%Prior_theta(k)%dist,&
                        par=INFER%Prior_theta(k)%par,&
                        gen=MC(:,k),feas=feas,err=err,mess=mess)
        if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
        if(.not.feas) then
            err=1;mess=trim(procname)//':'//trim(BaM_Message(9));return
        endif
    endif
enddo
compt=INFER%nFit

! remnant errors
do i=1,MODEL%nY
    do j=1,INFER%nremnant(i)
        compt=compt+1
        ! get mode
        call GetMode(DistId=INFER%Prior_RemnantSigma(i)%p(j)%dist,&
                 par=INFER%Prior_RemnantSigma(i)%p(j)%par,&
                 m=mode(compt),feas=feas,err=err,mess=mess)
        if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
        if(.not.feas) then
            err=1;mess=trim(procname)//':'//trim(BaM_Message(8));return
        endif
        ! do MC simulations
        if(present(MC)) then
            call generateSample(DistId=INFER%Prior_RemnantSigma(i)%p(j)%dist,&
                            par=INFER%Prior_RemnantSigma(i)%p(j)%par,&
                            gen=MC(:,compt),feas=feas,err=err,mess=mess)
            if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
            if(.not.feas) then
                err=1;mess=trim(procname)//':'//trim(BaM_Message(9));return
            endif
        endif
    enddo
enddo

end subroutine BaM_getPriorMC
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine CloseAllFiles(unt,unt_S,untMonitor)
! close all files opened during prediction experiments
integer(mik),intent(in)::unt(:),unt_S(:)
integer(mik),intent(in),optional::untMonitor
! locals
integer(mik)::i
logical::opened

do i=1,size(unt)
    inquire (UNIT=unt(i), OPENED=opened)
    if(opened) close(unt(i))
enddo
do i=1,size(unt_S)
    inquire (UNIT=unt_S(i), OPENED=opened)
    if(opened) close(unt_S(i))
enddo
if(present(untMonitor)) then
    inquire (UNIT=untMonitor,OPENED=opened)
    if(opened) close(untMonitor)
endif
end subroutine CloseAllFiles

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine ApplyModel_BaM(model,X,theta,Y,Dpar,state,feas,err,mess)
! Same as the ApplyModel of ModelLibrary_tools, but smart handling of theta:
! * fill in with values of FIX parameters
! * seek the current value of VAR/STOK parameters
use ModelLibrary_tools, only:ApplyModel
type(ModelType), intent(in)::model
real(mrk), intent(in)::X(:,:), theta(:)
real(mrk), intent(out)::Y(:,:),Dpar(:),state(:,:)
logical, intent(out)::feas
integer(mik), intent(out)::err
character(*),intent(out)::mess
! locals
character(250),parameter::procname='ApplyModel_BaM'
real(mrk)::fullTheta(MODEL%nTheta)
integer(mik)::i,k,t
real(mrk)::dummyX(1,model%nX),dummyY(1,model%nY),&
           dummyState(1,model%nstate),dummyDpar(model%nDpar)

if((INFER%nVar + INFER%nStok)==0) then ! no varying/ stok parameters, can run model in one go
    if(INFER%nFix==0) then ! no fixed parameters, easy
        fullTheta=theta
    else ! at least one fixed parameter
        k=0
        do i=1,MODEL%nTheta
            if(INFER%theta(i)%partype==FIX) then
                fullTheta(i)=INFER%theta(i)%init(1)
            else if(INFER%theta(i)%partype==STD) then
                k=k+1;fullTheta(i)=theta(k)
            endif
        enddo
    endif
    ! Apply model
    call ApplyModel(model,X,fullTheta,Y,Dpar,state,feas,err,mess)
    if(err>0) then;mess=trim(procname)//':'//trim(mess);return;endif
    if(.not. feas) return
else ! at least one varying / stok parameter, run model step by step
    k=1
    do t=1,INFER%nobs
        do i=1,MODEL%nTheta
            if(INFER%theta(i)%partype==FIX) then
                fullTheta(i)=INFER%theta(i)%init(1)
            else
                fullTheta(i)=theta(INFER%L(t,i))
            endif
        enddo
        dummyX(1,:)=X(t,:)
        call ApplyModel(model=model,X=dummyX,theta=fullTheta,&
                Y=dummyY,Dpar=dummyDpar,state=dummyState,&
                feas=feas,err=err,mess=mess)
        if(err>0) then;mess=trim(procname)//':'//trim(mess);return;endif
        if(.not. feas) return
        Y(t,:)=dummyY(1,:)
        state(t,:)=dummyState(1,:)
        if(k<=size(INFER%DparIndx)) then
            if(t==INFER%DparIndx(k)) then ! new DPar combination, save it
                Dpar( ((k-1)*MODEL%nDpar+1) : (k*MODEL%nDpar) )=dummyDpar
                k=k+1
            endif
        endif
    enddo
endif

end subroutine ApplyModel_BaM

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine Wextract(col,W,OUT)
! extract the relevant columns from a big W in ReadData sub
integer(mik),intent(in)::col(:)
real(mrk),intent(in)::W(:,:)
real(mrk),intent(out)::OUT(:,:)
!locals
integer(mik)::i

OUT=undefRN
do i=1,size(col)
    if(col(i)==0) then
        OUT(:,i)=0._mrk
    else
        OUT(:,i)=W(:,col(i))
    endif
enddo

end subroutine Wextract

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine Packer(theta,remnant,X,bias,OUT)
! pack the various types of parameters into a vector
real(mrk),intent(in)::theta(:),X(:,:),bias
type(parlist), intent(in)::remnant(:)
real(mrk),intent(out)::OUT(:)
! locals
integer(mik)::i,j,k

! theta
OUT(1:INFER%nFit)=theta
k=INFER%nFit
! Remnant parameters
do i=1,MODEL%nY
    OUT( (k+1):(k+INFER%nremnant(i)) )=remnant(i)%p
    k=k+INFER%nremnant(i)
enddo
! True inputs
do j=1,MODEL%nX
    do i=1,INFER%nObs
        if(INFER%Xu(i,j) > 0._mrk) then
            OUT(k+1)=X(i,j)
            k=k+1
        endif
    enddo
enddo
! Input biases
do j=1,MODEL%nX
    if(INFER%nXb(j) > 0) then
        OUT( (k+1):(k+INFER%nXb(j)) )=bias
        k=k+INFER%nXb(j)
    endif
enddo
! Output biases
do j=1,MODEL%nY
    if(INFER%nYb(j) > 0) then
        OUT( (k+1):(k+INFER%nYb(j)) )=bias
        k=k+INFER%nYb(j)
    endif
enddo

end subroutine Packer

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine Gcop_getV(x,distID,par,v,feas,err,mess)
! compute v=qnorm(cdf(x)), the gaussianized x value used in a Gaussian copula
use Distribution_tools, only: GetCdf,GetQuantile,GAUSS
character(*), intent(in)::DistID
real(mrk), intent(in)::x, par(:)
real(mrk), intent(out)::v
logical, intent(out)::feas
integer(mik), intent(out)::err
character(*),intent(out)::mess
! locals
character(250),parameter::procname='Gcop_getV'
real(mrk)::u

!Init
v=UndefRN;feas=.true.;err=0;mess=''
call GetCdf(DistId=DistId,x=x,par=par,cdf=u,feas=feas,err=err,mess=mess)
if(err>0) then;mess=trim(procname)//':'//trim(mess);feas=.false.;return;endif
if (.not. feas) return
call GetQuantile(DistId=GAUSS,p=u,par=(/0._mrk,1._mrk/),q=v,feas=feas,err=err,mess=mess)
if(err>0) then
    err=0;feas=.false. ! turn error into infeasible: corresponds to a p-value of 0/1
endif

end subroutine Gcop_getV

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine writeMonitor(unt,i,n)
! Write monitoring file
use utilities_dmsl_kit, only:number_string
integer(mik),intent(in)::unt,i,n
rewind(unt)
write (unt,'(A)') trim(number_string(i))//'/'//trim(number_string(n))
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function collapse_c(x)
! collapse a vector of characters into a long string using a comma separator
character(*), intent(in)::x(:)
character(len_uLongStr)::collapse_c
integer(mik)::i
collapse_c='['
do i=1,size(x)
    collapse_c=trim(collapse_c)//'"'//trim(x(i))//'"'
    if(i<size(x)) then;collapse_c=trim(collapse_c)//',';endif
enddo
collapse_c=trim(collapse_c)//']'
end function collapse_c

function collapse_i(x)
! collapse a vector of integers into a long string using a comma separator
use utilities_dmsl_kit,only:number_string
integer(mik), intent(in)::x(:)
character(len_uLongStr)::collapse_i
integer(mik)::i
collapse_i='['
do i=1,size(x)
    collapse_i=trim(collapse_i)//trim(number_string(x(i)))
    if(i<size(x)) then;collapse_i=trim(collapse_i)//',';endif
enddo
collapse_i=trim(collapse_i)//']'
end function collapse_i

function collapse_s(x)
! collapse a vector of slist objects into a long string using a comma separator
type(slist), intent(in)::x(:)
character(len_uLongStr)::collapse_s
integer(mik)::i
collapse_s='['
do i=1,size(x)
    collapse_s=trim(collapse_s)//trim(collapse_c(x(i)%s))
    if(i<size(x)) then;collapse_s=trim(collapse_s)//',';endif
enddo
collapse_s=trim(collapse_s)//']'
end function collapse_s

function collapse_p(x)
! collapse a vector of parameters objects into a long string using a comma separator
type(ParType), intent(in)::x(:)
character(len_uLongStr)::collapse_p
integer(mik)::i
character(5)::foo
collapse_p='['
do i=1,size(x)
    select case(x(i)%parType)
    case(-1)
        foo=FIX_str
    case(0)
        foo=STD_str
    case(1)
        foo=VAR_str
    case(2)
        foo=STOK_str
    end select
    collapse_p=trim(collapse_p)//'"'//trim(foo)//'"'
    if(i<size(x)) then;collapse_p=trim(collapse_p)//',';endif
enddo
collapse_p=trim(collapse_p)//']'
end function collapse_p

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine writeBamInfo(infoFile,err,mess)
! Write INFO file
use utilities_dmsl_kit,only:number_string,getSpareUnit
character(*),intent(in)::infoFile
integer(mik), intent(out)::err
character(*),intent(out)::mess
! locals
character(250),parameter::procname='writeBamInfo'
integer(mik)::unt
character(250),allocatable::head(:)

if(allocated(head)) deallocate(head);allocate(head(INFER%nInfer+1+INFER%nDpar))
call GetHeaders(head)
call getSpareUnit(unt,err,mess)
if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
open(unit=unt,file=trim(infoFile), status='replace',iostat=err)
if(err>0) then;call BaM_ConsoleMessage(messID_Open,trim(infoFile));endif
write(unt,*) '"ID":','"'//trim(MODEL%ID)//'"'
if(err>0) then;call BaM_ConsoleMessage(messID_Write,trim(infoFile));endif
write(unt,*) '"nX":',trim(number_string(MODEL%nX))
if(err>0) then;call BaM_ConsoleMessage(messID_Write,trim(infoFile));endif
write(unt,*) '"nY":',trim(number_string(MODEL%nY))
if(err>0) then;call BaM_ConsoleMessage(messID_Write,trim(infoFile));endif
write(unt,*) '"parName":',trim(collapse_c(MODEL%parname))
if(err>0) then;call BaM_ConsoleMessage(messID_Write,trim(infoFile));endif
write(unt,*) '"parType":',trim(collapse_p(INFER%theta))
if(err>0) then;call BaM_ConsoleMessage(messID_Write,trim(infoFile));endif
write(unt,*) '"dParName":',trim(collapse_c(MODEL%dparname))
if(err>0) then;call BaM_ConsoleMessage(messID_Write,trim(infoFile));endif
write(unt,*) '"stateName":',trim(collapse_c(MODEL%statename))
if(err>0) then;call BaM_ConsoleMessage(messID_Write,trim(infoFile));endif
write(unt,*) '"structuralParName":',trim(collapse_s(INFER%Parname_RemnantSigma))
if(err>0) then;call BaM_ConsoleMessage(messID_Write,trim(infoFile));endif
write(unt,*) '"estimatedX":',trim(collapse_i(INFER%nXerror))
if(err>0) then;call BaM_ConsoleMessage(messID_Write,trim(infoFile));endif
write(unt,*) '"estimatedXbias":',trim(collapse_i(INFER%nXb))
if(err>0) then;call BaM_ConsoleMessage(messID_Write,trim(infoFile));endif
write(unt,*) '"estimatedYbias":',trim(collapse_i(INFER%nYb))
if(err>0) then;call BaM_ConsoleMessage(messID_Write,trim(infoFile));endif
write(unt,*) '"MCMCheaders":',trim(collapse_c(head))
if(err>0) then;call BaM_ConsoleMessage(messID_Write,trim(infoFile));endif
close(unt)
end subroutine writeBamInfo

end module BaM_tools
