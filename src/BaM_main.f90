program BaM_main

use kinds_dmsl_kit
use utilities_dmsl_kit,only:number_string
use uniran1_dmsl_mod,only: seed_uniran
use ModelLibrary_tools,only:modelType,XtraSetup,XtraCleanUp,applyModel
use BayesianEstimation_tools, only:priorListType
use BaM_tools, only:plist,parlist,slist,parType,&
                    LoadBamObjects,BaM_ReadData,BaM_Fit,BaM_CookMCMC,BaM_SummarizeMCMC,BaM_computeDIC,&
                    Config_Read,Config_Read_oneRun,Config_Read_Model,Config_Read_Xtra,Config_Read_Data,&
                    Config_Read_RemnantSigma,Config_Read_MCMC,Config_Read_RunOptions,&
                    Config_Read_Residual,Config_Read_Cook,Config_Read_Summary,&
                    Config_Read_Pred_Master,Config_Read_Pred,Config_Finalize,Config_Read_Inputs,BaM_WriteOutputs,&
                    BaM_ConsoleMessage,BaM_PrintHelp,BaM_LoadMCMC,BaM_Residual,&
                    BaM_Prediction,XspagType,BaM_ReadSpag,BaM_ReadInputs,BaM_Cleanup,BaM_getDataFile
implicit none

!-----------------------
! Constants
character(len_stdStrD),parameter::Config_file_def="Config_BaM.txt"
character(len_stdStrD),parameter::Prior_file_def="PriorSimulations.txt"
character(len_stdStrD),parameter::priorCorrFile="PriorCorrelation.txt"
character(len_stdStrD),parameter::infoFile="INFO_BaM.txt"
character(len_stdStrD),parameter::MonitorExt=".monitor"
character(len_stdStrD),parameter::version="1.1.0 June 2025"
real(mrk),parameter::defaultstd=0.1_mrk
!-----------------------
! Config files
character(len_vLongStr)::workspace,Config_file,filePath,filePath2
character(len_longStr)::Config_RunOptions,Config_Model,Config_Xtra,Config_Data,&
                Config_MCMC,Config_Cooking,Config_summary,&
                Config_Residual,Config_Pred_Master,Config_Inputs
character(len_longStr),pointer::Config_RemnantSigma(:),Config_Pred(:)
character(len_vlongStr),allocatable::Config_RemnantList(:)
!-----------------------
! run options
logical::DoMCMC,DoSummary,DoResidual,DoPred
!-----------------------
! MCMC properties
integer(mik)::nAdapt,nCycles,nSlim,InitStdMode
real(mrk):: BurnFactor,MinMoveRate,MaxMoveRate,DownMult,UpMult
character(len_longStr)::MCMCFile
real(mrk), allocatable::theta_std0(:)
type(parlist),allocatable::RemnantSigma_std0(:)
!-----------------------
! Data
integer(mik),allocatable::XCol(:),XuCol(:),XbCol(:),XbindxCol(:)
integer(mik),allocatable::YCol(:),YuCol(:),YbCol(:),YbindxCol(:)
real(mrk), allocatable::X(:,:),Y(:,:),Xu(:,:),Yu(:,:),Xb(:,:),Yb(:,:),state(:,:),Dpar(:)
integer(mik), allocatable::Xbindx(:,:),Ybindx(:,:)
!-----------------------
! Inference
Type(ParType), pointer::theta(:)
real(mrk), pointer::theta0(:)
type(parlist),allocatable::RemnantSigma0(:)
type(plist),allocatable::Prior_RemnantSigma(:)
type(slist),allocatable::ParName_RemnantSigma(:)
character(len_longStr),allocatable::RemnantSigma_funk(:)
type(ModelType)::model
!-----------------------
! Post-processing
real(mrk), pointer::maxpost(:),mcmc(:,:),logpost(:)
character(len_longStr)::Residual_File,Cooking_File,Summary_File,DIC_File,xtendedMCMC_File,priorFile,YFile
character(len_vlongStr)::dummy_File
!-----------------------
! Prediction
integer(mik)::npred
logical::DoParametric
logical,allocatable::DoRemnant(:),DoTranspose(:),DoEnvelop(:),DoState(:),&
                     DoTranspose_S(:),DoEnvelop_S(:)
character(len_longStr),allocatable::YSpag_Files(:),Envelop_Files(:),SpagFiles_S(:),EnvelopFiles_S(:),head(:)
character(len_vlongStr),allocatable::XSpag_Files(:),fp1(:),fp2(:),fp3(:),fp4(:)
logical::PrintCounter
type(XspagType)::Xspag
real(mrk), allocatable::vTheta(:)
!-----------------------
! Misc.
integer(mik)::i,j,err,nobs,nc,nhead,nsim,narg,seed,ntheta
logical::IsMCMCLoaded,earlyStop,savePrior,oneRun,feas
character(len_vLongStr)::mess,datafile,arg
!-----------------------

!---------------------------------------------------------------------
!---------------------------------------------------------------------
! WELCOME !!!
!---------------------------------------------------------------------
!---------------------------------------------------------------------
call BaM_ConsoleMessage(1,'')
IsMCMCLoaded=.false.
earlyStop=.false.
savePrior=.false.
oneRun=.false.
priorFile=Prior_file_def

!---------------------------------------------------------------------
!---------------------------------------------------------------------
! INTERPRET COMMAND LINE ARGUMENTS
!---------------------------------------------------------------------
!---------------------------------------------------------------------
i=1;Config_file=''
narg=command_argument_count()
do while (i<=narg)
     call get_command_argument(i,arg)
     select case (arg)
     case ('-cf', '--config')
        i=i+1
        if(i<=narg) then
            call get_command_argument(i,arg)
            Config_file=trim(arg);i=i+1
        else
            call BaM_ConsoleMessage(-1,'-cf requires a path to a file')
        endif
     case ('-sp', '--saveprior')
        savePrior=.true.;i=i+1
        if(i<=narg) then
            call get_command_argument(i,arg)
            priorFile=trim(arg);i=i+1
        else
            call BaM_ConsoleMessage(-1,'--saveprior requires a file name')
        endif
     case ('-sd', '--seed')
        i=i+1
        if(i<=narg) then
            call get_command_argument(i,arg)
            read(arg,*,iostat=err) seed
            if(err==0) then
                call seed_uniran(put=seed)
                i=i+1
            else
                call BaM_ConsoleMessage(-1,'-sd requires the seed as an integer number')
            endif
        else
            call BaM_ConsoleMessage(-1,'-sd requires the seed as an integer number')
        endif
     case ('-rd', '--random')
        call seed_uniran(CPUtime=.true.)
        i=i+1
     case ('-dr', '--dontrun')
        earlyStop=.true.
        i=i+1
     case ('-or', '--onerun')
        oneRun=.true.;i=i+1
        if(i<=narg) then
            call get_command_argument(i,arg)
            YFile=trim(arg);i=i+1
        else
            call BaM_ConsoleMessage(-1,'--onerun requires a file name')
        endif
    case ('-v', '--version')
        write(*,*) 'version: ', trim(version)
        STOP
     case ('-h', '--help')
        call BaM_PrintHelp()
        STOP
     case default
        write(*,*) 'Unrecognized command-line option: ', trim(arg)
        call BaM_PrintHelp()
        call BaM_ConsoleMessage(-1,'')
     end select
  end do

! Use default config file if not provided
if (Config_file=='') then ! workspace not passed through command line, try reading it from wkFile
    Config_file=trim(Config_file_def)
endif

!---------------------------------------------------------------------
!---------------------------------------------------------------------
! READ CONFIG FILES
!---------------------------------------------------------------------
!---------------------------------------------------------------------
call BaM_ConsoleMessage(100,'')

! Read main config file
if(oneRun) then
    call Config_Read_oneRun(trim(Config_file),&
                 workspace,Config_Model,Config_Xtra,Config_Inputs,&
                 err,mess)
else
    call Config_Read(trim(Config_file),&
                     workspace,&
                     Config_RunOptions,Config_Model,Config_Xtra,Config_Data,&
                     Config_RemnantSigma,Config_MCMC,&
                     Config_Cooking,Config_Summary,&
                     Config_Residual,Config_Pred_Master,&
                     err,mess)
endif
if(err>0) then; call BaM_ConsoleMessage(-1,trim(mess));endif

if(.not.oneRun) then
    ! Read run options
    filePath=trim(workspace)//trim(Config_RunOptions)
    call Config_Read_RunOptions(filePath,&
                                DoMCMC,DoSummary,&
                                DoResidual,DoPred,&
                                err,mess)
    if(err>0) then; call BaM_ConsoleMessage(-1,trim(mess));endif
endif

! model setup
filePath=trim(workspace)//trim(Config_Model)
call Config_Read_Model(filePath,&
                       model%ID,model%nX,model%nY,&
                       theta,&
                       err,mess)
if(err>0) then; call BaM_ConsoleMessage(-1,trim(mess));endif
model%ID="MDL_"//trim(model%ID)
ntheta=size(theta)
model%ntheta=ntheta
filePath=trim(workspace)//trim(Config_Xtra)
call Config_Read_Xtra(file=filePath,&
                      ID=model%ID,xtra=model%xtra,err=err,mess=mess)
if(err>0) then; call BaM_ConsoleMessage(-1,trim(mess));endif

if(oneRun) then
    ! Finish model setup (Dpar, states)
    call XtraSetup(model,err,mess)
    if(err>0) then; call BaM_ConsoleMessage(-1,trim(mess));endif
    ! READ X
    call Config_Read_Inputs(file=trim(workspace)//trim(Config_Inputs),&
                      DataFile=datafile,nHeader=nhead,nobs=nobs,ncol=nc,&
                      err=err,mess=mess)
    if(err>0) then; call BaM_ConsoleMessage(-1,trim(mess));endif
    call BaM_getDataFile(datafile,workspace,err,mess)
    if(err>0) then; call BaM_ConsoleMessage(-1,trim(mess));endif
    allocate(X(nobs,model%nX),Y(nobs,model%nY),head(model%nY),state(nobs,model%nState),Dpar(model%nDpar),vtheta(ntheta))
    call BaM_ReadInputs(file=trim(datafile),nhead=nhead,X=X,err=err,mess=mess)
    if(err>0) then; call BaM_ConsoleMessage(-1,trim(mess));endif
    ! Apply model
    do i=1,ntheta
        vtheta(i)=theta(i)%init(1)
    enddo
    call ApplyModel(model=model,X=X,theta=vtheta,Y=Y,Dpar=Dpar,state=state,feas=feas,err=err,mess=mess)
    if(err>0) then; call BaM_ConsoleMessage(-1,trim(mess));endif
    ! WRITE Y
    do i=1,model%nY
        head(i)='Y'//number_string(i)
    enddo
    filePath=trim(workspace)//trim(Yfile)
    call BaM_WriteOutputs(file=filePath,head=head,Y=Y,err=err,mess=mess)
    if(err>0) then; call BaM_ConsoleMessage(-1,trim(mess));endif
    ! model cleanup
    call XtraCleanup(model,err,mess)
    if(err>0) then; call BaM_ConsoleMessage(-1,trim(mess));endif
    ! Get out
    STOP
endif

! read data
allocate(Xcol(model%nX),XuCol(model%nX),XbCol(model%nX),XbindxCol(model%nX))
allocate(Ycol(model%nY),YuCol(model%nY),YbCol(model%nY),YbindxCol(model%nY))
filePath=trim(workspace)//trim(Config_Data)
call Config_Read_Data(file=filePath,&
                      DataFile=datafile,nHeader=nhead,nobs=nobs,ncol=nc,&
                      XCol=XCol,XuCol=XuCol,XbCol=XbCol,XbindxCol=XbindxCol,&
                      YCol=YCol,YuCol=YuCol,YbCol=YbCol,YbindxCol=YbindxCol,&
                      err=err,mess=mess)
if(err>0) then; call BaM_ConsoleMessage(-1,trim(mess));endif

! Does datafile exists?
call BaM_getDataFile(datafile,workspace,err,mess)
if(err>0) then; call BaM_ConsoleMessage(-1,trim(mess));endif

allocate(X(nobs,model%nX),Xu(nobs,model%nX),Xb(nobs,model%nX),Xbindx(nobs,model%nX))
allocate(Y(nobs,model%nY),Yu(nobs,model%nY),Yb(nobs,model%nY),Ybindx(nobs,model%nY))
call BaM_ReadData(file=datafile,nrow=nobs,ncol=nc,nHeader=nhead,&
                  XCol=XCol,XuCol=XuCol,XbCol=XbCol,XbindxCol=XbindxCol,&
                  YCol=YCol,YuCol=YuCol,YbCol=YbCol,YbindxCol=YbindxCol,&
                  X=X,Y=Y,Xu=Xu,Yu=Yu,Xb=Xb,Yb=Yb,Xbindx=Xbindx,Ybindx=Ybindx,&
                  err=err,mess=mess)
if(err>0) then; call BaM_ConsoleMessage(-1,trim(mess));endif

! set remnant error model
allocate(RemnantSigma_funk(model%nY),RemnantSigma0(model%nY),Prior_RemnantSigma(model%nY))
allocate(Config_RemnantList(model%nY),Parname_RemnantSigma(model%nY))
do i=1,model%nY
    Config_RemnantList(i)=trim(workspace)//trim(Config_RemnantSigma(i))
enddo
call Config_Read_RemnantSigma(files=Config_RemnantList,&
                      RemnantSigma_funk=RemnantSigma_funk,&
                      Parname=Parname_RemnantSigma,&
                      RemnantSigma0=RemnantSigma0,&
                      Prior_RemnantSigma=Prior_RemnantSigma,&
                      err=err,mess=mess)
if(err>0) then; call BaM_ConsoleMessage(-1,trim(mess));endif

! Finalize configuration
call Config_Finalize(workspace=trim(workspace),&
                     datafile=datafile,nrow=nobs,ncol=nc,nHeader=nHead,&
                     theta=theta,&
                     theta0=theta0,err=err,mess=mess)
if(err>0) then; call BaM_ConsoleMessage(-1,trim(mess));endif

! Load all info into BAM objects
call LoadBamObjects(X=X,Xu=Xu,Xb=Xb,Xbindx=Xbindx,& ! observed inputs and their uncertainties
                    Y=Y,Yu=Yu,Yb=Yb,Ybindx=Ybindx,& ! observed outputs and their uncertainties
                    ID=model%ID,&               ! Model ID
                    theta=theta,&               ! parameters theta
                    RemnantSigma_funk=RemnantSigma_funk,& ! chosen f in { residual var = f(Qrc) }
                    Parname_RemnantSigma=Parname_RemnantSigma,& ! names
                    Prior_RemnantSigma=Prior_RemnantSigma,& ! priors
                    priorCorrFile=trim(workspace)//trim(priorCorrFile),&
                    infoFile=trim(workspace)//trim(infoFile),&
                    xtra=model%xtra,&
                    nstate=model%nState,err=err,mess=mess)! error handling
if(err>0) then; call BaM_ConsoleMessage(-1,trim(mess));endif

if(earlyStop) then ! user just wanted to write the INFO file
    call BaM_ConsoleMessage(20, '')
    call BaM_Cleanup(err=err,mess=mess)
    if(err>0) then; call BaM_ConsoleMessage(18,trim(mess));endif
    call BaM_ConsoleMessage(999, '')
    STOP
endif

! MCMC Config file
allocate(theta_std0(size(theta0)),RemnantSigma_std0(model%nY))
filePath=trim(workspace)//trim(Config_MCMC)
call Config_Read_MCMC(filePath,&
                  theta0,RemnantSigma0,defaultStd,&
                  nAdapt,nCycles,InitStdMode,&
                  MinMoveRate,MaxMoveRate,DownMult,UpMult,&
                  theta_std0,RemnantSigma_std0,&
                  MCMCfile,&
                  err,mess)
if(err>0) then; call BaM_ConsoleMessage(-1,trim(mess));endif
! Also read Cook-MCMC to get burn & nslim
filePath=trim(workspace)//trim(Config_Cooking)
call Config_Read_Cook(filePath,&
                      burnFactor,nSlim,Cooking_File,err,mess)
if(err>0) then; call BaM_ConsoleMessage(-1,trim(mess));endif

call BaM_ConsoleMessage(101,model%ID)

!---------------------------------------------------------------------
!---------------------------------------------------------------------
! MCMC SAMPLING
!---------------------------------------------------------------------
!---------------------------------------------------------------------
if(DoMCMC) then
    call BaM_ConsoleMessage(4,'')
    ! Go!
    filePath=trim(workspace)//trim(MCMCfile)
    filePath2=trim(workspace)//trim(Config_MCMC)//trim(MonitorExt)
    call BaM_Fit(theta0=theta0,RemnantSigma0=RemnantSigma0, &!initial values for teta and remnant std
                 theta_std0=theta_std0,RemnantSigma_std0=RemnantSigma_std0,& ! initial values for the std of jump distribution
                 nAdapt=nAdapt,nCycles=nCycles,&
                 MinMoveRate=MinMoveRate,MaxMoveRate=MaxMoveRate,&
                 DownMult=DownMult,UpMult=UpMult,&
                 OutFile=filePath, & ! Output file (for MCMC samples)
                 MonitorFile=filePath2, & ! monitoring file (for MCMC samples)
                 err=err,mess=mess)
    if(err>0) then; call BaM_ConsoleMessage(-1,trim(mess));endif
    call BaM_ConsoleMessage(5, trim(trim(workspace)//trim(MCMCfile)))
    ! load MCMC samples and write cooked samples
    call BaM_ConsoleMessage(6,'')
    filePath=trim(workspace)//trim(MCMCfile)
    call BaM_LoadMCMC(MCMCFile=filePath,&
                  burnFactor=burnFactor,Nslim=Nslim,& ! Read properties
                  maxpost=maxpost,mcmc=mcmc,logpost=logpost,&
                  err=err,mess=mess)!out
    if(err>0) then; call BaM_ConsoleMessage(-1,trim(mess));endif
    IsMCMCLoaded=.true.
    filePath=trim(workspace)//trim(Cooking_File)
    call BaM_CookMCMC(mcmc=mcmc,logpost=logpost,&
                      OutFile=filePath,&
                      err=err,mess=mess)
    if(err>0) then; call BaM_ConsoleMessage(-1,trim(mess));endif
    call BaM_ConsoleMessage(7,'')
endif

!---------------------------------------------------------------------
!---------------------------------------------------------------------
! LOAD cooked MCMC (used for subsequent analyses) if not done yet
!---------------------------------------------------------------------
!---------------------------------------------------------------------
if(.not.IsMCMCLoaded) then
    if(DoSummary .or. DoResidual) then
        filePath=trim(workspace)//trim(Cooking_File)
        call BaM_LoadMCMC(MCMCFile=filePath,&
                      burnFactor=0._mrk,Nslim=1,& ! Read properties
                      maxpost=maxpost,mcmc=mcmc,logpost=logpost,&
                      err=err,mess=mess)!out
        if(err>0) then; call BaM_ConsoleMessage(-1,trim(mess));endif
        IsMCMCLoaded=.true.
    endif
endif

!---------------------------------------------------------------------
!---------------------------------------------------------------------
! MCMC SUMMARY
!---------------------------------------------------------------------
!---------------------------------------------------------------------
if(DoSummary) then
    call BaM_ConsoleMessage(8,'')
    ! config file
    filePath=trim(workspace)//trim(Config_Summary)
    call Config_Read_Summary(filePath,&
                          Summary_File,DIC_File,xtendedMCMC_File,err,mess)
    if(err>0) then; call BaM_ConsoleMessage(-1,trim(mess));endif
    ! summarize
    filePath=trim(workspace)//trim(Summary_File)
    call BaM_SummarizeMCMC(mcmc=mcmc,logpost=logpost,&
                      OutFile=filePath,&
                      err=err,mess=mess)
    if(err>0) then; call BaM_ConsoleMessage(-1,trim(mess));endif
    ! Compute DICs if requested
    if(.not.(DIC_File=='')) then
        if(xtendedMCMC_File=='') then
            dummy_File='' ! Extended MCMC with prior, lkh, hyper, deviance is not written
        else
            dummy_File=trim(workspace)//trim(xtendedMCMC_File)
        endif
        filePath=trim(workspace)//trim(DIC_File)
        call BaM_computeDIC(mcmc=mcmc,logpost=logpost,&
                            DICfile=filePath,&
                            MCMCfile=dummy_File,err=err,mess=mess)
        if(err>0) then; call BaM_ConsoleMessage(-1,trim(mess));endif
    endif
    call BaM_ConsoleMessage(9,'')

endif

!---------------------------------------------------------------------
!---------------------------------------------------------------------
! Residual diagnostics
!---------------------------------------------------------------------
!---------------------------------------------------------------------
if(DoResidual) then
    call BaM_ConsoleMessage(10,'')
    ! read config file
    filePath=trim(workspace)//trim(Config_Residual)
    call Config_Read_Residual(filePath,Residual_File,err,mess)
    if(err>0) then; call BaM_ConsoleMessage(-1,trim(mess));endif
    ! Perform runs
    filePath=trim(workspace)//trim(Residual_File)
    call BaM_Residual(maxpost=maxpost,&
                      OutFile=filePath,&
                      err=err,mess=mess)
    if(err>0) then; call BaM_ConsoleMessage(-1,trim(mess));endif
    call BaM_ConsoleMessage(11,'')
endif

!---------------------------------------------------------------------
!---------------------------------------------------------------------
! Prediction experiments
!---------------------------------------------------------------------
!---------------------------------------------------------------------
if(DoPred) then
    call BaM_ConsoleMessage(12,'')
    ! read master config file
    filePath=trim(workspace)//trim(Config_Pred_Master)
    call Config_Read_Pred_Master(filePath,npred,Config_Pred,err,mess)
    if(err>0) then; call BaM_ConsoleMessage(-1,trim(mess));endif
    allocate(XSpag_Files(model%nX))
    allocate(DoRemnant(model%nY),YSpag_Files(model%nY),DoTranspose(model%nY),&
             DoEnvelop(model%nY),Envelop_Files(model%nY))
    allocate(Xspag%nspag(model%nX))
    allocate(Xspag%spag(model%nX))
    allocate(DoState(model%nState),SpagFiles_S(model%nState),DoTranspose_S(model%nState),&
             DoEnvelop_S(model%nState),EnvelopFiles_S(model%nState))
    allocate(fp1(model%nY),fp2(model%nY),fp3(model%nState),fp4(model%nState))
    do i=1,npred
        call BaM_ConsoleMessage(14,trim(number_string(i))//'/'//trim(number_string(npred)))
        ! read individual pred config files
        filePath=trim(workspace)//trim(Config_Pred(i))
        call Config_Read_Pred(filePath,&
                            XSpag_Files,Xspag%nobs,Xspag%nspag,&
                            DoParametric,DoRemnant,nsim,YSpag_Files,&
                            DoTranspose,DoEnvelop,Envelop_Files,PrintCounter,&
                            DoState,SpagFiles_S,DoTranspose_S,DoEnvelop_S,EnvelopFiles_S,&
                            err,mess)
        if(err>0) then; call BaM_ConsoleMessage(-1,trim(mess));endif
        ! Do spaghetti files exist?
        do j=1,size(XSpag_Files)
            call BaM_getDataFile(XSpag_Files(j),workspace,err,mess)
            if(err>0) then; call BaM_ConsoleMessage(-1,trim(mess));endif
        enddo
        ! read spaghettis
        call BaM_ReadSpag(Xspag,XSpag_Files,err,mess)
        if(err>0) then; call BaM_ConsoleMessage(-1,trim(mess));endif
        ! Do Prediction
        if(nsim>0) then ! prior sampling
            fp1=trim(workspace)//YSpag_Files
            fp2=trim(workspace)//Envelop_Files
            fp3=trim(workspace)//SpagFiles_S
            fp4=trim(workspace)//EnvelopFiles_S
            filePath=trim(workspace)//trim(Config_Pred(i))//trim(MonitorExt)
            call BaM_Prediction(nsim=nsim,&
                          Xspag=Xspag,DoParametric=DoParametric,DoRemnant=DoRemnant,&
                          SpagFiles=fp1,PrintCounter=PrintCounter,&
                          DoTranspose=DoTranspose,DoEnvelop=DoEnvelop,&
                          EnvelopFiles=fp2,&
                          DoState=DoState,DoTranspose_S=DoTranspose_S,DoEnvelop_S=DoEnvelop_S,&
                          SpagFiles_S=fp3,&
                          EnvelopFiles_S=fp4,&
                          MonitorFile=filePath,& ! monitoring file
                          savePrior=savePrior,priorFile=trim(workspace)//trim(priorFile),&
                          err=err,mess=mess)
            if(err>0) then; call BaM_ConsoleMessage(17,trim(mess));endif
        else ! posterior sampling
            if(.not.IsMCMCLoaded) then
                filePath=trim(workspace)//trim(Cooking_File)
                call BaM_LoadMCMC(MCMCFile=filePath,&
                      burnFactor=0._mrk,Nslim=1,& ! Read properties
                      maxpost=maxpost,mcmc=mcmc,logpost=logpost,&
                      err=err,mess=mess)!out
                if(err>0) then; call BaM_ConsoleMessage(-1,trim(mess));endif
                IsMCMCLoaded=.true.
            endif
            fp1=trim(workspace)//YSpag_Files
            fp2=trim(workspace)//Envelop_Files
            fp3=trim(workspace)//SpagFiles_S
            fp4=trim(workspace)//EnvelopFiles_S
            filePath=trim(workspace)//trim(Config_Pred(i))//trim(MonitorExt)
            call BaM_Prediction(mcmc=mcmc,maxpost=maxpost,&
                          Xspag=Xspag,DoParametric=DoParametric,DoRemnant=DoRemnant,&
                          SpagFiles=fp1,PrintCounter=PrintCounter,&
                          DoTranspose=DoTranspose,DoEnvelop=DoEnvelop,&
                          EnvelopFiles=fp2,&
                          DoState=DoState,DoTranspose_S=DoTranspose_S,DoEnvelop_S=DoEnvelop_S,&
                          SpagFiles_S=fp3,&
                          EnvelopFiles_S=fp4,&
                          MonitorFile=filePath,& ! monitoring file
                          err=err,mess=mess)
            if(err>0) then; call BaM_ConsoleMessage(17,trim(mess));endif
        endif
    enddo
    call BaM_ConsoleMessage(13,'')
endif

!---------------------------------------------------------------------
!---------------------------------------------------------------------
! ALL DONE !!!
!---------------------------------------------------------------------
!---------------------------------------------------------------------
call BaM_Cleanup(err=err,mess=mess)
if(err>0) then; call BaM_ConsoleMessage(18,trim(mess));endif

call BaM_ConsoleMessage(999, '')
!read(*,*)
end program BaM_main
