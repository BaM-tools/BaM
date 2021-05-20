module SuspendedLoad_model

!~**********************************************************************
!~* Purpose: Suspended Load model
!~**********************************************************************
!~* Programmer: Ben Renard, Emeline Perret, Jerôme Le Coz, Irstea Lyon
!~**********************************************************************
!~* Last modified: 30/01/2019
!~**********************************************************************
!~* Comments:
!~**********************************************************************
!~* References:
!~**********************************************************************
!~* 2Do List: 
!~**********************************************************************
!~* Quick description of public procedures:
!~*		1. SuspendedLoad_GetParNumber, number of parameters (5 + optional ones)
!~*		2. SuspendedLoad_Apply, compute suspended load capacity = f(stage, parameters)
!~**********************************************************************

use kinds_dmsl_kit ! numeric kind definitions from DMSL

implicit none
Private
public :: SuspendedLoad_GetParNumber,SuspendedLoad_Apply

Contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
pure subroutine SuspendedLoad_GetParNumber(npar,err,mess)

!^**********************************************************************
!^* Purpose: number of parameters of the suspended load model
!^**********************************************************************
!^* Programmer: Ben Renard, Emeline Perret, Jerôme Le Coz, Irstea Lyon
!^**********************************************************************
!^* Last modified: 30/01/2019
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^* OUT
!^*		1. npar, par. number
!^*		2. err, error code; <0:Warning, ==0:OK, >0: Error
!^*		3. mess, error message
!^**********************************************************************

integer(mik), intent(out)::npar,err
character(*),intent(out)::mess
! locals

err=0;mess='';npar=5

end subroutine SuspendedLoad_GetParNumber
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine SuspendedLoad_Apply(h,theta,Qs,feas,err,mess)

!^**********************************************************************
!^* Purpose: apply the suspended load model
!^**********************************************************************
!^* Programmer: Ben Renard, Emeline Perret, Jerôme Le Coz, Irstea Lyon
!^**********************************************************************
!^* Last modified: 30/01/2019
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*		1. h, stage
!^*		2. theta, parameters
!^* OUT
!^*		1. Qs, suspended load capacity 
!^*		2. feas, feasible? 
!^*		3. err, error code; <0:Warning, ==0:OK, >0: Error
!^*		4. mess, error message
!^**********************************************************************

real(mrk), intent(in)::h(:),theta(:)
real(mrk), intent(out)::Qs(:)
logical, intent(out)::feas(:)
integer(mik), intent(out)::err
character(*),intent(out)::mess
!locals
integer(mik)::i,n
real(mrk)::num,den

! Initialize
err=0;mess='';feas=.true.;Qs=undefRN;
n=size(h)

! Check parameter feasibility
if( theta(1)<=0._mrk .or. theta(3)<=0._mrk .or. theta(5)<=0._mrk) then
    feas=.false.
    return
endif

! Compute
do i=1,n
    if(h(i)<theta(2) .or. h(i)<theta(4)) then 
        feas(i)=.false.
        !Qs(i)=0._mrk
    else
        num=(h(i)-theta(2))**theta(3)
        den=(h(i)-theta(4))**theta(5)
        Qs(i)=theta(1) * (num/den)
    endif 
enddo

end subroutine SuspendedLoad_Apply
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


end module SuspendedLoad_model