program TEST_BaM

use kinds_dmsl_kit
use Distribution_tools,only:GenerateSample,Generate
use ModelLibrary_tools,only:ApplyModel,GetModelParNumber,modelType,MDL_Linear
use BayesianEstimation_tools, only:CreateFlatPriorList,priorListType
use BaM_tools, only:plist,parlist,LoadBamObjects,BaM_Fit
implicit none

!-----------------------
integer(mik),parameter:: nobs=100,nX=3,nY=2
character(250),parameter:: Xdist="Gaussian", modelID=MDL_Linear
!-----------------------
integer(mik),parameter::nAdapt=100,nCycles=100
real(mrk),parameter::MinMoveRate=0.1_mrk,MaxMoveRate=0.5_mrk,DownMult=0.9_mrk,UpMult=1.1_mrk
character(250), parameter::OutFile="mcmc.txt"
!-----------------------
real(mrk)::X(nobs,nX),Y(nobs,nY),theta(nX*nY),sigma(nY),mu(nY),eps
type(ModelType)::model
type(priorListType)::Prior_theta(nX*nY)
type(plist)::Prior_RemnantSigma(nY)
type(parlist)::RemnantSigma0(nY),RemnantSigma_std0(nY)
integer(mik)::i,j,err
logical::feas
character(250)::RemnantSigma_funk(nY)
character(250)::mess

model%ID=modelID
model%nX=nX
model%nY=nY
call GetModelParNumber(model=model,npar=model%ntheta,err=err,mess=mess)
! generate theta
call GenerateSample(DistId="Gaussian",par=(/0._mrk,1._mrk/),gen=theta,feas=feas,err=err,mess=mess)
! generate sigmas
call GenerateSample(DistId="Gaussian",par=(/0._mrk,1._mrk/),gen=sigma,feas=feas,err=err,mess=mess)
sigma=sigma**2
! generate covariates
do j=1,nX
    call GenerateSample(DistId=Xdist,par=(/0._mrk,1._mrk/),gen=X(:,j),feas=feas,err=err,mess=mess)
enddo
! generate outputs
do i=1,nobs
    call ApplyModel(model=model,X=X(i,:),theta=theta,Y=mu,feas=feas,err=err,mess=mess)
    do j=1,nY
        call Generate(DistId="Gaussian",par=(/0._mrk,sigma(j)/),gen=eps,feas=feas,err=err,mess=mess)
        Y(i,j)=mu(j)+eps
    enddo
enddo
do j=1,nY
    RemnantSigma_funk(j)="Constant"
enddo

call CreateFlatPriorList(npar=model%ntheta,PriorList=Prior_theta,err=err,mess=mess)
do j=1,nY
    allocate(Prior_RemnantSigma(j)%p(1))
    call CreateFlatPriorList(npar=1,PriorList=Prior_RemnantSigma(j)%p,err=err,mess=mess)
enddo

call LoadBamObjects(X=X,Xsigma=0._mrk*X,Y=Y,Ysigma=0._mrk*Y,& ! observations and uncertainties
                        ID=modelID,&               ! Model ID 
                        RemnantSigma_funk=RemnantSigma_funk,& ! chosen f in { residual var = f(Qrc) }
				        Prior_theta=Prior_theta,Prior_RemnantSigma=Prior_RemnantSigma,& ! priors
				        err=err,mess=mess)! error handling

do j=1,nY
    allocate(RemnantSigma0(j)%p(1))
    allocate(RemnantSigma_std0(j)%p(1))
    RemnantSigma0(j)%p=1._mrk
    RemnantSigma_std0(j)%p=0.1_mrk
enddo
call BaM_Fit(theta0=theta,RemnantSigma0=RemnantSigma0, &!initial values for teta and remnant std
			 theta_std0=0._mrk*theta+0.1_mrk,RemnantSigma_std0=RemnantSigma_std0,& ! initial values for the std of jump distribution
			 nAdapt=nAdapt,nCycles=nCycles,& 
			 MinMoveRate=MinMoveRate,MaxMoveRate=MaxMoveRate,&
			 DownMult=DownMult,UpMult=UpMult,&
			 OutFile=OutFile, & ! Output file (for MCMC samples)
			 err=err,mess=mess)	

write(*,*) "Finito!"
read(*,*)			        
end program TEST_BaM