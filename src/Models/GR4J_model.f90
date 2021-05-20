module GR4J_model

!~**********************************************************************
!~* Purpose: GR4J model
!~**********************************************************************
!~* Programmer: Ben Renard, Irstea Lyon
!~**********************************************************************
!~* Last modified: 21/08/2017
!~**********************************************************************
!~* Comments:
!~**********************************************************************
!~* References:
!~**********************************************************************
!~* 2Do List: 
!~**********************************************************************
!~* Quick description of public procedures:
!~*		1. GR4J_GetParNumber, number of parameters (6)
!~*		2. GR4J_Apply, run GR4J (through all time steps)
!~*		3. GR4J_XtraLoad, load xtra information (pco file)
!~**********************************************************************

use kinds_dmsl_kit ! numeric kind definitions from DMSL

implicit none
Private
public :: GR4J_GetParNumber,GR4J_Apply,GR4J_XtraRead

Contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
pure subroutine GR4J_GetParNumber(npar,err,mess)

!^**********************************************************************
!^* Purpose: number of parameters of GR4J
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified: 21/08/2017
!^**********************************************************************
!^* Comments: 6 parameters because the initial values for 2 stores are
!^*           considered as parameters
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* OUT
!^*		1. npar, par. number
!^*		2.err, error code; <0:Warning, ==0:OK, >0: Error
!^*		3.mess, error message
!^**********************************************************************

integer(mik), intent(out)::npar,err
character(*),intent(out)::mess
! locals
err=0;mess='';npar=6

end subroutine GR4J_GetParNumber
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine GR4J_Apply(X,theta,pcofile,Y,state,feas,err,mess)

!^**********************************************************************
!^* Purpose: Run GR4J over a time series
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified:21/08/2017
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*		1. matrix of inputs (precip & pet, Nt*2)
!^*		2. theta, parameters
!^*		3. pcofile, PCO configuration file from DMSL
!^* OUT
!^*		1. Y, output matrix (Q, Nt*1)
!^*		2. state, states (production & routing, Nt*2)
!^*		3. feas, feasability flag
!^*		4. err, error code; <0:Warning, ==0:OK, >0: Error
!^*		5. mess, error message
!^**********************************************************************
use dynamicModelLibrary,only:GR4J_DMDL,DMDL_controlModel,DMDL_runModel,&
                        DMDL_setModel,DMDL_getModelInfo
real(mrk), intent(in)::X(:,:), theta(:)
character(*),intent(in)::pcofile
real(mrk), intent(out)::Y(:,:),state(:,:)
logical, intent(out)::feas
integer(mik), intent(out)::err
character(*),intent(out)::mess
!locals
integer(mik), parameter::qindx=1,prodindx=3,routingindx=4
character(250),parameter::procname='GR4J_Apply'
integer(mik)::nstate,npar
character(250)::foo
real(mrk),allocatable::par(:),istate(:)
integer(mik)::t,Nt

err=0;mess='';feas=.true.;Y=undefRN;state=undefRN
Nt=size(X,dim=1)

! Set model
call DMDL_setModel(modelID=(/GR4J_DMDL/),setupCmd="",err=err,message=mess) 
if(err>0) then;mess=trim(procname)//':'//trim(mess);return;endif
! Read PCO file for general info
call DMDL_getModelInfo(modelID=(/GR4J_DMDL/),infoCmd=trim(pcofile),modelName=foo,&
        nstate=nstate,npar=npar,err=err,message=mess)
if(err>0) then;mess=trim(procname)//':'//trim(mess);return;endif
! allocate par / state
if(allocated(par)) deallocate(par);allocate(par(npar))
if(allocated(istate)) deallocate(istate);allocate(istate(nstate))
par=(/1.0_mrk,0.0_mrk,theta/) ! in DMSL:(rmult,radd,X1,X2,X3,X4,initprod,initrouting)
! Set par & initial values
call DMDL_controlModel(modelID=(/GR4J_DMDL/),parIn=par,setS0in=.true.,feas=feas,err=err,message=mess)
if(err>0) then;mess=trim(procname)//':'//trim(mess);return;endif
if(.not.feas) return
! run it for all time steps
do t=1,Nt
    ! Note: the "dataProps" below works for GR4J but may need to be changed for other 
    !       hydro-model in DMSL. I'm not sure what these variables are doing...
    call DMDL_runModel(modelID=(/GR4J_DMDL/),&
                       runitCmd='',iT=undefIN,dataProps=(/1._mrk/),& 
                       input=X(t,:),state=istate,feas=feas,err=err,message=mess)
    if(err>0) then;mess=trim(procname)//':'//trim(mess);return;endif
    if(.not.feas) return
    Y(t,:)=(/istate(qindx)/)
    state(t,:)=(/istate(prodindx),istate(routingindx)/)
enddo
end subroutine GR4J_Apply
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine GR4J_XtraRead(file,xtra,err,mess)

!^**********************************************************************
!^* Purpose: pass address of PCO file into xtra
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified:21/08/2017
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*		1. file, PCO setup file
!^* OUT
!^*		1. xtra, xtra information
!^*		2. err, error code; <0:Warning, ==0:OK, >0: Error
!^*		3. mess, error message
!^**********************************************************************
use types_dmsl_kit, only:data_ricz_type
character(*), intent(in)::file
type(data_ricz_type),intent(out):: xtra
integer(mik), intent(out)::err
character(*),intent(out)::mess

err=0;mess='';xtra%cs1=trim(file)

end subroutine GR4J_XtraRead
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module GR4J_model