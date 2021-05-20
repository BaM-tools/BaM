module Mixture_model

!~**********************************************************************
!~* Purpose: Mixture model (geochemestry)
!~**********************************************************************
!~* Programmer: Ben Renard, Irstea Lyon
!~**********************************************************************
!~* Last modified: 29/01/2018
!~**********************************************************************
!~* Comments: Mixture model computing Yhat as a sum of contributions,
!~*           Yhat=sum(theta_i*X_i). Each theta_i is between 0 and 1,
!~*           and theta_last=1-sum[i=1,last-1](theta_i).
!~*           Typical usage: downstream concentration is equal to the 
!~*           weighted sum of upstream concentrations (sum(weights)=1).
!~* ADDED 05/04/2018: multi-event and missing contribution
!~**********************************************************************
!~* References:
!~**********************************************************************
!~* 2Do List: 
!~**********************************************************************
!~* Quick description of public procedures:
!~*    1. Mixture_GetParNumber, number of parameters
!~*    2. Mixture_Apply, compute Yhat=sum(theta_i*X_i)
!~**********************************************************************

use kinds_dmsl_kit ! numeric kind definitions from DMSL

implicit none
Private
public :: Mixture_GetParNumber,Mixture_Apply,Mixture_XtraRead

Contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
pure subroutine Mixture_GetParNumber(nS,nEvt,nElt,missingSource,npar,err,mess)

!^**********************************************************************
!^* Purpose: number of parameters of the Mixture model
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified: 05/04/2018
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1. nS, number of sources
!^*    2. nEvt, number of events
!^*    3. nElt, number of elements
!^*    4. missingSource, allow missing source?
!^* OUT
!^*    1. npar, par. number
!^*    2.err, error code; <0:Warning, ==0:OK, >0: Error
!^*    3.mess, error message
!^**********************************************************************

integer(mik), intent(in)::nS,nEvt,nElt
logical,intent(in)::missingSource
integer(mik), intent(out)::npar,err
character(*),intent(out)::mess
! locals

err=0;mess=''
if(missingSource) then
    npar=nS*nEvt+nElt
else
    npar=(nS-1)*nEvt
endif

end subroutine Mixture_GetParNumber
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Mixture_Apply(nS,nEvt,nElt,missingSource,X,theta,Y,theta_last,feas,err,mess)

!^**********************************************************************
!^* Purpose: apply the Mixture model
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified:05/04/2018
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1. nS, number of sources
!^*    2. nEvt, number of events
!^*    3. nElt, number of elements
!^*    4. missingSource, allow missing source?
!^*    5. input matrix
!^*    6. theta, parameter vector
!^* OUT
!^*    1. Y, output vector
!^*    2. theta_last, last thetas (one minus sum of others for each event)
!^*    3.feas, feasibility (all theta between 0 and 1)
!^*    4.err, error code; <0:Warning, ==0:OK, >0: Error
!^*    5.mess, error message
!^**********************************************************************

integer(mik), intent(in)::nS,nEvt,nElt
logical,intent(in)::missingSource
real(mrk), intent(in)::X(:,:), theta(:)
real(mrk), intent(out)::Y(:),theta_last(nEvt)
logical, intent(out)::feas
integer(mik), intent(out)::err
character(*),intent(out)::mess
!locals
integer(mik)::i,j,nX,nY
real(mrk)::dummytheta(nS),dummylast,sources(nElt,nS)

err=0;mess='';Y=undefRN;theta_last=undefRN

do j=1,nEvt
    if(missingSource) then
        dummytheta=theta( (j-1)*nS+1 : (j*nS) )
        if(sum(dummytheta)>1._mrk) then;feas=.false.;return;endif
        dummylast=1._mrk-sum(dummytheta)
    else
        dummytheta(1:(nS-1))=theta( (j-1)*(nS-1)+1 : (j*(nS-1)) )
        dummytheta(nS)=1._mrk-sum(dummytheta(1:(nS-1)))
        dummylast=dummytheta(nS)
    endif
    if(any(dummytheta<0._mrk) .or. any(dummytheta>1._mrk)) then
        feas=.false.;return
    endif
    do i=1,nS
        sources(:,i)=pack(X(:,2+i),nint(X(:,1))==j)
    enddo
    Y( (j-1)*nElt+1 : (j*nElt) )=matmul(sources,dummytheta)
    if(missingSource) then
        Y( (j-1)*nElt+1 : (j*nElt) ) = Y( (j-1)*nElt+1 : (j*nElt) ) + dummylast*theta( (nEvt*nS+1):)
    endif
    theta_last(j)=dummylast
enddo

end subroutine Mixture_Apply
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Mixture_XtraRead(file,xtra,err,mess)
!^**********************************************************************
!^* Purpose: Read Xtra information (nS,nEvt,nElt,missingSource)
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified:05/04/2018
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*		1. file, Xtra file
!^* OUT
!^*		1. xtra, xtra information
!^*		2. err, error code; <0:Warning, ==0:OK, >0: Error
!^*		3. mess, error message
!^**********************************************************************
use types_dmsl_kit, only:data_ricz_type
use utilities_dmsl_kit,only:getSpareUnit
character(*), intent(in)::file
type(data_ricz_type),intent(out):: xtra
integer(mik), intent(out)::err
character(*),intent(out)::mess
! locals
character(250),parameter::procname='Mixture_XtraRead'
integer(mik)::unt

err=0;mess=''
call getSpareUnit(unt,err,mess)
if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
open(unit=unt,file=trim(file), status='old', iostat=err)
if(err>0) then;mess=trim(procname)//':problem opening file '//trim(file);return;endif
read(unt,*,iostat=err) xtra%is1 ! nS
if(err/=0) then;mess=trim(procname)//':problem reading file '//trim(file);return;endif
read(unt,*,iostat=err) xtra%is2 ! nEvt
if(err/=0) then;mess=trim(procname)//':problem reading file '//trim(file);return;endif
read(unt,*,iostat=err) xtra%is3 ! nElt
if(err/=0) then;mess=trim(procname)//':problem reading file '//trim(file);return;endif
read(unt,*,iostat=err) xtra%ls1 ! missingSource
if(err/=0) then;mess=trim(procname)//':problem reading file '//trim(file);return;endif

close(unt)

end subroutine Mixture_XtraRead
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module Mixture_model