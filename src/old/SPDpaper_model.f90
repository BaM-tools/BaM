module SPDpaper_model

!~**********************************************************************
!~* Purpose: 3-control RC model used in SPD paper
!~**********************************************************************
!~* Programmer: Ben Renard & Matteo Darienzo, Irstea Lyon
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
!~*		1. 
!~*		2. 
!~*		3.
!~**********************************************************************

use kinds_dmsl_kit ! numeric kind definitions from DMSL

implicit none
Private
public :: SPDpaper_GetParNumber,SPDpaper_Apply

Contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
pure subroutine SPDpaper_GetParNumber(npar,err,mess)

!^**********************************************************************
!^* Purpose: number of parameters
!^**********************************************************************
!^* Programmer: Ben Renard & Matteo Darienzo, Irstea Lyon
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
!^*		1. npar, par. number
!^*		2.err, error code; <0:Warning, ==0:OK, >0: Error
!^*		3.mess, error message
!^**********************************************************************
integer(mik), intent(out)::npar,err
character(*),intent(out)::mess
! locals
character(250),parameter::procname='SPDpaper_GetParNumber'

err=0;mess='';npar=9

end subroutine SPDpaper_GetParNumber
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine SPDpaper_Apply(h,theta,Q,kappa,newton,feas,err,mess)

!^**********************************************************************
!^* Purpose: apply SFD
!^**********************************************************************
!^* Programmer: Ben Renard & Matteo Darienzo, Irstea Lyon
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
!^*		1. h, input vector (h1,h2)
!^*		2. theta, parameter vector
!^*		3. [NewtonOption], options for Newton numerical resolution
!^* OUT
!^*		1. Q, discharge
!^*		2. kappa, transition stage
!^*		3. feas, feasible?
!^*		4. err, error code; <0:Warning, ==0:OK, >0: Error
!^*		5. mess, error message
!^**********************************************************************
use types_dmsl_kit,only:data_ricz_type
use utilities_dmsl_kit, only:cleanPointers
use numerix_dmsl_kit, only:znewton

real(mrk), intent(in)::h(:),theta(9)
real(mrk), intent(out)::Q(:),kappa(3),newton(3)
logical, intent(out)::feas(:)
integer(mik), intent(out)::err
character(*),intent(out)::mess
!locals
character(250),parameter::procname='SPDpaper_Apply'
real(mrk),parameter::smallH=0.001_mrk,bigH=10000._mrk
integer(mik),parameter::itmax=10000
integer(mik)::nobs,t,fcalls
real(mrk)::a(3),b(3),c(3),delta,y2,alpha,y0,fkappa
type(data_ricz_type)::NewtonInfo


err=0;mess='';feas=.true.;Q=undefRN;kappa=undefRN;newton=undefRN
nobs=size(h)
! make sense of theta
a(1)=theta(1)
b(1)=theta(2)
c(1)=theta(3)
a(2)=theta(4)
b(2)=theta(5)
c(2)=theta(6)
a(3)=theta(7)
b(3)=theta(8)
c(3)=theta(9)
kappa(1)=b(1)
kappa(3)=b(3)
! Solve continuity constraint
delta=b(1)-b(2)
alpha=a(1)/a(2)
y0=delta*(c(2)/(c(2)-c(1)))
if(associated(NewtonInfo%rp0)) nullify(NewtonInfo%rp0);allocate(NewtonInfo%rp0(4))
NewtonInfo%rp0(1)=c(2)
NewtonInfo%rp0(2)=alpha
NewtonInfo%rp0(3)=c(1)
NewtonInfo%rp0(4)=delta

if(c(2)>c(1) .and. delta<=0._mrk) then
    newton(1)=2
else if(c(2)<c(1) .and. delta>=0._mrk) then 
    newton(1)=3
else if(c(2)>c(1) .and. delta>0._mrk) then 
    newton(1)=1
else if(c(2)<c(1) .and. delta<0._mrk) then 
    newton(1)=4
endif

if(newton(1)==2 .or. newton(1)==3) then 
    call znewton(evalFunc=fNewton,dataIN=NewtonInfo,x1=0._mrk,x2=bigH,&
             tolX=0.00001_mrk,tolF=0.0000000001_mrk,xscale=5._mrk,&
             fscale=1._mrk,itmax=itmax,&
             xroot=kappa(2),froot=fkappa,fcalls=fcalls,err=err,message=mess)
    if(err>0) then
        feas=.false.
        call cleanPointers(NewtonInfo,what=1,err=err,message=mess)
        return
    endif
    ! Does not find root, cycle
    if(fcalls>itmax) then
        feas=.false.
        call cleanPointers(NewtonInfo,what=1,err=err,message=mess)
        return
    endif
    newton(2:3)=kappa(2)
    kappa(2)=kappa(2)+b(2)
endif

if(newton(1)==1 .or. newton(1)==4) then 
    call znewton(evalFunc=fNewton,dataIN=NewtonInfo,x1=max(0._mrk,delta),x2=y0,&
             tolX=0.00001_mrk,tolF=0.0000000001_mrk,xscale=5._mrk,&
             fscale=1._mrk,itmax=itmax,&
             xroot=newton(2),froot=fkappa,fcalls=fcalls,err=err,message=mess)
    if(err>0 .or. fcalls>itmax) then
        newton(2)=undefRN
    endif
    call znewton(evalFunc=fNewton,dataIN=NewtonInfo,x1=y0,x2=bigH,&
             tolX=0.00001_mrk,tolF=0.0000000001_mrk,xscale=5._mrk,&
             fscale=1._mrk,itmax=itmax,&
             xroot=newton(3),froot=fkappa,fcalls=fcalls,err=err,message=mess)
    if(err>0) then
        feas=.false.
        newton(3)=undefRN
        call cleanPointers(NewtonInfo,what=1,err=err,message=mess)
        return
    endif
    ! Does not find root, cycle
    if(fcalls>itmax) then
        feas=.false.
        newton(3)=undefRN
        call cleanPointers(NewtonInfo,what=1,err=err,message=mess)
        return
    endif
    kappa(2)=newton(3)+b(2)
endif

! Check feasability
if(any(a<=0._mrk) .or. any(c<=0._mrk)) then
    feas=.false.;return
endif
if(kappa(1)>kappa(2) .or. kappa(2)>kappa(3)) then 
    feas=.false.;return
endif
do t=1,nobs
    if(h(t)<=kappa(1)) then;Q(t)=0._mrk;cycle;endif
    if(h(t)<=kappa(2) .and. h(t)>kappa(1)) then
        if(h(t)-b(1)<=0._mrk) then;feas(t)=.false.;cycle;endif
        Q(t)=a(1)*(h(t)-b(1))**c(1);cycle
    endif
    if(h(t)<=kappa(3) .and. h(t)>kappa(2)) then
        if(h(t)-b(2)<=0._mrk) then;feas(t)=.false.;cycle;endif
        Q(t)=a(2)*(h(t)-b(2))**c(2);cycle
    endif
    if(h(t)>kappa(3)) then
        if(h(t)-b(2)<=0._mrk .or. h(t)-b(3)<=0) then;feas(t)=.false.;cycle;endif
        Q(t)=a(2)*(h(t)-b(2))**c(2) + a(3)*(h(t)-b(3))**c(3);cycle
    endif
enddo
end subroutine SPDpaper_Apply
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!==============!
! Private subs !
!==============!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine fNewton(dataIN,dataOUT,x,feas,fx,dfdxV,err,message)
! function to nullify by Newton method for SFD_Val_General
use kinds_dmsl_kit
use types_dmsl_kit,only:data_ricz_type
implicit none
type(data_ricz_type),intent(in),optional::dataIN
type(data_ricz_type),intent(inout),optional::dataOUT
real(mrk),intent(in)::x
logical(mlk),intent(out)::feas
real(mrk),intent(out),optional::fx,dfdxV(:)
integer(mik),intent(out)::err
character(*),intent(out)::message
    
feas=.true.;err=0;message='';
if(present(fx)) fx=undefRN
if(present(dfdxV)) dfdxV=undefRN 
if(present(fx))then
    if(dataIN%rp0(4)==0._mrk) then
        fx=(dataIN%rp0(1)-dataIN%rp0(3))*log(x)-log(dataIN%rp0(2))
    else
        fx=dataIN%rp0(1)*log(x)-log(dataIN%rp0(2))-dataIN%rp0(3)*log(x-dataIN%rp0(4))
    endif
endif
if(present(dfdxV))then
    if(dataIN%rp0(4)==0._mrk) then
        dfdxV(1)=(dataIN%rp0(1)-dataIN%rp0(3))/x
    else
        dfdxV(1)=dataIN%rp0(1)/x-dataIN%rp0(3)/(x-dataIN%rp0(4))
    endif
endif

end subroutine fNewton
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module SPDpaper_model