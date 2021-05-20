program TEST_BaM_RatingCurve

use kinds_dmsl_kit
use ModelLibrary_tools,only:ApplyModel,GetModelParNumber,modelType,MDL_BaRatin
use BayesianEstimation_tools, only:CreateFlatPriorList,priorListType
use BaM_tools, only:plist,parlist,LoadBamObjects,BaM_Fit,BaM_ReadData,&
                    Config_Read_Model,Config_Read_Data,Config_Read_MCMC,&
                    Config_Read_Xtra,Config_Read_RemnantSigma
use DataRW_tools,only:DatRead
implicit none

!-----------------------
character(250),parameter::workspace="Tests\BaM_Sauze\"
real(mrk),parameter::defaultstd=0.1_mrk
!-----------------------
! MCMC properties
integer(mik)::nAdapt,nCycles,Nburn,Nread,nSlim,InitStdMode
real(mrk):: BurnFactor,MinMoveRate,MaxMoveRate,DownMult,UpMult
character(250)::MCMCFile
logical::saveMCMC
!-----------------------
real(mrk), allocatable::X(:,:),Y(:,:),Xsigma(:,:),Ysigma(:,:),sigma(:),mu(:),theta_std0(:)
real(mrk), pointer::theta0(:)
type(priorListType),pointer::Prior_theta(:)
type(plist),allocatable::Prior_RemnantSigma(:)
type(parlist),allocatable::RemnantSigma0(:),RemnantSigma_std0(:)
integer(mik),allocatable::XCol(:),XsigmaCol(:),YCol(:),YsigmaCol(:)
character(250),allocatable::RemnantSigma_funk(:)
type(ModelType)::model
integer(mik)::i,j,err,nt,nc,nhead
character(250)::mess,datafile
logical::feas
!-----------------------

! model setup
call Config_Read_Model(trim(workspace)//"Config_Model.txt",&
                       model%ID,model%nX,model%nY,model%ntheta,&
                       theta0,&
                       Prior_theta,&
                       err,mess)
model%ID="MDL_"//model%ID                      
call Config_Read_Xtra(files=(/trim(workspace)//"Config_ControlMatrix.txt"/),&
                      ID=model%ID,xtra=model%xtra,err=err,mess=mess) 
                       
! read data
allocate(Xcol(model%nX),XsigmaCol(model%nX))
allocate(Ycol(model%nY),YsigmaCol(model%nY))
call Config_Read_Data(file=trim(workspace)//"Config_Data.txt",&
                      DataFile=datafile,nHeader=nhead,nobs=nt,ncol=nc,&
                      XCol=XCol,XsigmaCol=XsigmaCol,&
                      YCol=YCol,YsigmaCol=YsigmaCol,&
                      err=err,mess=mess)

allocate(X(nt,model%nX),Xsigma(nt,model%nX))
allocate(Y(nt,model%nY),Ysigma(nt,model%nY))                 
call BaM_ReadData(file=datafile,nrow=nt,ncol=nc,nHeader=nhead,&
                  XCol=XCol,XsigmaCol=XsigmaCol,YCol=YCol,YsigmaCol=YsigmaCol,&
                  X=X,Y=Y,Xsigma=Xsigma,Ysigma=Ysigma,&
                  err=err,mess=mess)

! set remnant error model
allocate(RemnantSigma_funk(model%nY),RemnantSigma0(model%nY),Prior_RemnantSigma(model%nY))
call Config_Read_RemnantSigma(files=(/trim(workspace)//"Config_RemnantSigma.txt"/),&
                      RemnantSigma_funk=RemnantSigma_funk,&
                      RemnantSigma0=RemnantSigma0,&
                      Prior_RemnantSigma=Prior_RemnantSigma,&
                      err=err,mess=mess)

! Load all info into BAM objects
call LoadBamObjects(X=X,Xsigma=Xsigma,Y=Y,Ysigma=Ysigma,& ! observations and uncertainties
                        ID=model%ID,&               ! Model ID 
                        RemnantSigma_funk=RemnantSigma_funk,& ! chosen f in { residual var = f(Qrc) }
				        Prior_theta=Prior_theta,Prior_RemnantSigma=Prior_RemnantSigma,& ! priors
				        xtra=model%xtra,&
				        err=err,mess=mess)! error handling


! MCMC properties
allocate(theta_std0(model%ntheta),RemnantSigma_std0(model%nY))
call Config_Read_MCMC(trim(workspace)//"Config_MCMC.txt",&
                      theta0,RemnantSigma0,defaultStd,&
                      nAdapt,nCycles,Nburn,Nread,nSlim,InitStdMode,&
                      BurnFactor,MinMoveRate,MaxMoveRate,DownMult,UpMult,&
                      theta_std0,RemnantSigma_std0,&
                      saveMCMC,MCMCfile,&
                      err,mess)
! Go!
call BaM_Fit(theta0=theta0,RemnantSigma0=RemnantSigma0, &!initial values for teta and remnant std
			 theta_std0=theta_std0,RemnantSigma_std0=RemnantSigma_std0,& ! initial values for the std of jump distribution
			 nAdapt=nAdapt,nCycles=nCycles,& 
			 MinMoveRate=MinMoveRate,MaxMoveRate=MaxMoveRate,&
			 DownMult=DownMult,UpMult=UpMult,&
			 OutFile=trim(workspace)//trim(MCMCfile), & ! Output file (for MCMC samples)
			 err=err,mess=mess)	

write(*,*) "Finito!"
read(*,*)			        
end program TEST_BaM_RatingCurve