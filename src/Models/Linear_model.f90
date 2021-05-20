module Linear_model

!~**********************************************************************
!~* Purpose: Linear model
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
!~*		1. LM_GetParNumber, number of parameters (nX*nY)
!~*		2. LM_Apply, compute Y_hat=X*theta (theta packaged as a nX*nY matrix)
!~**********************************************************************

use kinds_dmsl_kit ! numeric kind definitions from DMSL

implicit none
Private
public :: LM_GetParNumber,LM_Apply

Contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
pure subroutine LM_GetParNumber(nX,nY,npar,err,mess)

!^**********************************************************************
!^* Purpose: number of parameters of the Linear model
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
!^*		1. nX, number of input variables
!^*		1. nY, number of output variables
!^* OUT
!^*		1. npar, par. number
!^*		2.err, error code; <0:Warning, ==0:OK, >0: Error
!^*		3.mess, error message
!^**********************************************************************

integer(mik), intent(in)::nX,nY
integer(mik), intent(out)::npar,err
character(*),intent(out)::mess
! locals

err=0;mess='';npar=nX*nY

end subroutine LM_GetParNumber
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine LM_Apply(X,theta,Y,err,mess)

!^**********************************************************************
!^* Purpose: apply the linear model
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
!^*		1. input vector
!^*		2. theta, parameters
!^* OUT
!^*		1. Y, output vector
!^*		3.err, error code; <0:Warning, ==0:OK, >0: Error
!^*		4.mess, error message
!^**********************************************************************

real(mrk), intent(in)::X(:,:), theta(:)
real(mrk), intent(out)::Y(:,:)
integer(mik), intent(out)::err
character(*),intent(out)::mess
!locals
integer(mik)::nX,nY
real(mrk), allocatable::dummytheta(:,:)

err=0;mess='';Y=undefRN;

nX=size(X,dim=2);nY=size(Y,dim=2)
if(allocated(dummytheta)) deallocate(dummytheta);allocate(dummytheta(nX,nY))
dummytheta=reshape(theta,(/nX,nY/))
Y=matmul(X,dummytheta)

end subroutine LM_Apply
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module Linear_model