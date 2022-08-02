module Recession_model_h

!~**********************************************************************
!~* Purpose: Retrospective analysis using recession of stage time series
!~**********************************************************************
!~* Programmer: Matteo Darienzo, Irstea Lyon
!~**********************************************************************
!~* Last modified: 03/04/2018
!~**********************************************************************
!~* Comments:
!~**********************************************************************
!~* References:
!~**********************************************************************
!~* 2Do List:
!~**********************************************************************
!~* Quick description of public procedures:
!~*     1.
!~*     2.
!~*     3.
!~**********************************************************************




use kinds_dmsl_kit ! numeric kind definitions from DMSL

implicit none
Private
public :: Recession_h_GetParNumber,Recession_h_Apply

Contains






!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
pure subroutine Recession_h_GetParNumber(npar,err,mess)

!^**********************************************************************
!^* Purpose: number of parameters
!^**********************************************************************
!^* Programmer: Matteo Darienzo, Irstea Lyon
!^**********************************************************************
!^* Last modified: 03/04/2018
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^* OUT
!^*     1. npar, par. number
!^*     2.err, error code; <0:Warning, ==0:OK, >0: Error
!^*     3.mess, error message
!^**********************************************************************
integer(mik), intent(out)::npar,err
character(*),intent(out)::mess
! locals
character(250),parameter::procname='Recession_h_GetParNumber'

!err=0;mess='';npar=5   ! model full with two exponentials + asymptote = 5 parameters
!err=0;mess='';npar=3    ! model with one exponential + asymptote = 3 parameters
err=0;mess='';npar=4    ! model with two exponentials + asymptote and normalised = 4 parameters
!err=0;mess='';npar=2    ! model with 1 exponential + asymptote and normalised = 2 parameters
end subroutine Recession_h_GetParNumber








!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Recession_h_Apply(time,theta,h,feas,err,mess)

!^**********************************************************************
!^* Purpose: apply SFD
!^**********************************************************************
!^* Programmer: Matteo Darienzo, Irstea Lyon
!^**********************************************************************
!^* Last modified: 03/04/2018
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*     1. h, input vector (h1,h2)
!^*     2. theta, parameter vector
!^*     3. [NewtonOption], options for Newton numerical resolution
!^* OUT
!^*     1. Q, discharge
!^*     2. kappa, transition stage
!^*     3. feas, feasible?
!^*     4. err, error code; <0:Warning, ==0:OK, >0: Error
!^*     5. mess, error message
!^**********************************************************************
use types_dmsl_kit,only:data_ricz_type
use utilities_dmsl_kit, only:cleanPointers

real(mrk), intent(in):: time(:)
!real(mrk), intent(in):: theta(5)  !model without normalisation
real(mrk), intent(in):: theta(4)   !model with normalisation and 2 exp
!real(mrk), intent(in):: theta(2)   !model with normalisation and only one exponential
real(mrk), intent(out):: h(:)

!locals
character(250),parameter::procname='Recession_h_Apply'
logical, intent(out)::feas
integer(mik), intent(out)::err
character(*),intent(out)::mess

integer(mik)::nobs,t
!real(mrk):: a(3),b(2)
real(mrk):: a(2),b(2)  !model with normalisation and 2 exp
!real(mrk):: a(1),b(1)  !model with normalisation and only one exponential

err=0;
mess='';
feas = .TRUE.

nobs=size(h)
a(1)=theta(1)
b(1)=theta(2)
a(2)=theta(3)
b(2)=theta(4)
!a(3)=theta(5)

do t=1,nobs
    !h(t)= a(1)*exp(-b(1)*time(t)) + a(2)*exp(-b(2)*time(t)) + a(3)
    h(t)= a(1)*exp(-b(1)*time(t))+ (1-a(1)-a(2))*exp(-b(2)*time(t)) + a(2)  !model with normalisation
    !h(t)= (1- a(1))*exp(-b(1)*time(t))+ a(1)
enddo

!if (b(1) <= b(2)) then
!   feas = .FALSE.
!   return
!end if

!if ((a(1) <= a(2)) .or. (a(2) <= a(3))) then
!   feas = .FALSE.
!   return
!end if

!if ((a(1) <= a(2))) then
!   feas = .FALSE.
!   return
!end if

end subroutine Recession_h_Apply




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module Recession_model_h

