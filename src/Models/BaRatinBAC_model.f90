module BaRatinBAC_model

!~**********************************************************************
!~* Purpose: BaRatin rating curve parameterized in terms of b-a-c
!~*          parameters (in lieu of k-a-c in original BaRatin)
!~**********************************************************************
!~* Programmer: Ben Renard, Irstea Lyon
!~*             Modified from Valentin Mansanarez's code
!~*             With great assistance from Matteo Darienzo
!~**********************************************************************
!~* Last modified: 10/01/2019
!~**********************************************************************
!~* Comments:
!~**********************************************************************
!~* References: Mansanarez, V. (2016), Non-unique stage-discharge
!~*             relations: Bayesian analysis of complex rating curves
!~*             and their uncertainties, PhD thesis
!~*             Mansanarez et al. (2019) Shift happens! Adjusting
!~*             stage-discharge rating curves to riverbed morphological
!~*             changes at known times. Water Resources Research.
!~**********************************************************************
!~* 2Do List:
!~**********************************************************************
!~* Quick description of public procedures:
!~*     1. BaRatinBAC_GetParNumber, number of parameters of the RC
!~*     2. BaRatinBAC_Apply, compute Q=f(H|theta)
!~*     3. BaRatinBAC_XtraRead, read Bonnifait matrix
!~**********************************************************************

use kinds_dmsl_kit ! numeric kind definitions from DMSL

implicit none
Private
public :: BaRatinBAC_GetParNumber,BaRatinBAC_Apply,BaRatinBAC_XtraRead

Contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
pure subroutine BaRatinBAC_GetParNumber(ControlMatrix,npar,err,mess)

!^**********************************************************************
!^* Purpose: number of parameters of the RC
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified: 10/01/2019
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*     1. ControlMatrix, control Matrix
!^* OUT
!^*     1. npar, par. number
!^*     2.err, error code; <0:Warning, ==0:OK, >0: Error
!^*     3.mess, error message
!^**********************************************************************
use BaRatin_model, only:BaRatin_GetParNumber
integer(mik), intent(in)::ControlMatrix(:,:)
integer(mik), intent(out)::npar,err
character(*),intent(out)::mess

call BaRatin_GetParNumber(ControlMatrix,npar,err,mess)

end subroutine BaRatinBAC_GetParNumber
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine BaRatinBAC_Apply(H,theta,ControlMatrix,hmax,Q,k,ktype,feas,err,mess)

!^**********************************************************************
!^* Purpose: Apply rating curve equation Q=f(H|theta)
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^*             Modified from Valentin Mansanarez's code
!^**********************************************************************
!^* Last modified: 10/01/2019
!^**********************************************************************
!^* Comments: ktype is the type of solution of the continuity equation
!^*           * 0: control addition => simple solution k=b
!^*           * 1 and above: control substitution => need to solve the
!^*             continuity equation numerically; the number gives the
!^*             type of numerical solution, as described in the
!^*             analytical analysis of the continuity equation described
!^*             in the WRR paper below
!^**********************************************************************
!^* References: Mansanarez, V. (2016), Non-unique stage-discharge
!^*             relations: Bayesian analysis of complex rating curves
!^*             and their uncertainties, PhD thesis
!^*             Mansanarez et al. (2019) Shift happens! Adjusting
!^*             stage-discharge rating curves to riverbed morphological
!^*             changes at known times. Water Resources Research.
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*     1. H, water stage
!^*     2. theta, parameters (b-a-c parameterization)
!^*     3. ControlMatrix, control Matrix in Laurent's framework
!^*     4. hmax, higher bound of feasability range for the rating curve
!^* OUT
!^*     1. Q, RC-computed runoff
!^*     2. k, transition stages deduced by continuity
!^*     3. ktype, type of solution of the continuity equation (see comment above)
!^*     4. feas, feasible?
!^*     5. err, error code; <0:Warning, ==0:OK, >0: Error
!^*     6. mess, error message
!^**********************************************************************
use BaRatin_model, only:BaRatin_computeQ
real(mrk), intent(in)::H(:),theta(:),hmax
integer(mik), intent(in)::ControlMatrix(:,:)
real(mrk), intent(out)::Q(:),k(:),ktype(:)
logical, intent(out)::feas(:)
integer(mik), intent(out)::err
character(*),intent(out)::mess
!locals
character(250),parameter::procname='BaRatinBAC_Apply'
integer(mik)::ncontrol,i,p
real(mrk), allocatable::b(:),a(:),c(:)
logical::ok

err=0;mess='';Q=undefRN;k=undefRN;feas=.true.;ktype=undefIN

! Get number of controls (=number of ranges)
ncontrol=size(ControlMatrix,1)

! check size of theta
if(size(theta)/=3*ncontrol) then
    err=1;mess=trim(procname)//':Fatal:Incorrect size [theta]';feas=.false.;return
endif

! check size of k
if(size(k)/=ncontrol) then
    err=1;mess=trim(procname)//':Fatal:Incorrect size [k]';feas=.false.;return
endif

! Allocate and assign values to b,a,c
allocate(b(ncontrol),a(ncontrol),c(ncontrol))
do i=1,ncontrol
    p=3*(i-1)
    b(i)=theta(p+1);a(i)=theta(p+2);c(i)=theta(p+3)
enddo

! Apply Continuity Constraint & check feasability
call ApplyContinuity(b,a,c,ControlMatrix,hmax,k,ktype,ok,err,mess)
if(.not.ok) then;feas=.false.;return;endif

! Everything checked OK, can finally apply RC equation!
call BaRatin_computeQ(H,b,a,c,k,ControlMatrix,Q,feas)

end subroutine BaRatinBAC_Apply
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine BaRatinBAC_XtraRead(file,xtra,err,mess)

!^**********************************************************************
!^* Purpose: Read Xtra information for RC: Bonnifait matrix
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified: 23/07/2024
!^**********************************************************************
!^* Comments: Additional restriction compared with original BaRatin:
!^*           a new control cannot take over A SUM of controls.
!^*           Reason: analytical analysis of the continuity equation has
!^*           not been thoroughly performed yet, and the continuity
!^*           function may well have more than 2 roots max.
!^*           July 2024: Added a few additional checks, possibly partly
!^*           redundant with the previous restriction, but better safe than sorry
!^**********************************************************************
!^* References: Mansanarez et al. (2019) Shift happens! Adjusting
!^*             stage-discharge rating curves to riverbed morphological
!^*             changes at known times. Water Resources Research.
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*     1. file, Xtra file
!^* OUT
!^*     1. xtra, xtra information
!^*     2. err, error code; <0:Warning, ==0:OK, >0: Error
!^*     3. mess, error message
!^**********************************************************************
use BaRatin_model,only:CheckControlMatrix
use utilities_dmsl_kit,only:getSpareUnit,getNumItemsInFile
use types_dmsl_kit, only:data_ricz_type
character(*), intent(in)::file
type(data_ricz_type),intent(out):: xtra
integer(mik), intent(out)::err
character(*),intent(out)::mess
! locals
character(250),parameter::procname='BaRatinBAC_XtraRead'
integer(mik)::unt,nrange,i,n
character(250)::foo
logical::feas

err=0;mess=''
! open file and get number of items
call getSpareUnit(unt,err,mess)
if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
open(unit=unt,file=trim(file), status='old', iostat=err)
if(err>0) then;mess=trim(procname)//':problem opening file '//trim(file);return;endif
call getNumItemsInFile(unt=unt,preRwnd=.false.,nskip=0,nitems=n,postPos=0,&
                       jchar=foo,err=err,message=mess)
if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
! read control matrix
if(associated(xtra%rpm1)) nullify(xtra%rpm1);allocate(xtra%rpm1(n-1,n-1))
do i=1,(n-1)
    read(unt,*,iostat=err) xtra%rpm1(i,:)
    if(err/=0) then;mess=trim(procname)//':problem reading file '//trim(file);return;endif
enddo
! read hmax
read(unt,*,iostat=err) xtra%rs1
if(err/=0) then;mess=trim(procname)//':problem reading file '//trim(file);return;endif
close(unt)
! check control matrix
call CheckControlMatrix(nint(xtra%rpm1),feas,err,mess) ! same as original BaRatin
if(err/=0 .or. (.not.feas) ) then
    mess=trim(procname)//':invalid control matrix';
    err=1;return
endif
! Check additional constraint (see comment)
nrange=size(xtra%rpm1,dim=1)
if(nrange>1) then
    do i=2,nrange
        if(sum(xtra%rpm1(i,:))<sum(xtra%rpm1(i-1,:))) then ! decreasing number of active controls
            err=1
            mess=trim(procname)//':FATAL: the number of active controls is not allowed to decrease'&
                                 &'when stage increases. Use original BaRatin for such configuration.'
            return
        endif
        if(all( xtra%rpm1(i,1:(i-1)) - xtra%rpm1(i-1,1:(i-1)) == 0 )&
          .and. xtra%rpm1(i,i)==1 .and. xtra%rpm1(i-1,i)==0) then
            ! control addition case : OK, do nothing
        elseif( xtra%rpm1(i-1,i-1)==1 .and. xtra%rpm1(i-1,i)==0&
          .and. xtra%rpm1(i,i-1)==0 .and. xtra%rpm1(i,i)==1 ) then
            ! control replacement case : check that only the latest control is replaced
            if(i>2) then
                if( .not.all( xtra%rpm1(i,1:(i-2)) - xtra%rpm1(i-1,1:(i-2)) == 0 ) ) then
                    err=2
                    mess=trim(procname)//':FATAL: unsupported control matrix: '&
                                         &'unsupported control replacement.'
                    return
                endif
            endif
        else
            ! Should never arrive here?
            err=3
            mess=trim(procname)//':FATAL: unsupported control matrix.'
            return
        endif
    enddo
endif
end subroutine BaRatinBAC_XtraRead
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!==============!
! Private subs !
!==============!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine ApplyContinuity(b,a,c,ControlMatrix,hmax,k,ktype,feas,err,mess)

!^**********************************************************************
!^* Purpose: Apply Continuity constraints & evaluate feasability
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^*             Modified from Valentin Mansanarez's code
!^*             With great assistance from Matteo Darienzo
!^**********************************************************************
!^* Last modified:10/01/2019
!^**********************************************************************
!^* Comments: ktype is the type of solution of the continuity equation
!^*           * 0: control addition => simple solution k=b
!^*           * 1 and above: control substitution => need to solve the
!^*             continuity equation numerically; the number gives the
!^*             type of numerical solution, as described in the
!^*             analytical analysis of the continuity equation described
!^*             in the WRR paper below
!^**********************************************************************
!^* References: Mansanarez et al. (2019) Shift happens! Adjusting
!^*             stage-discharge rating curves to riverbed morphological
!^*             changes at known times. Water Resources Research.
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*     1.b, offsets
!^*     2.a, coefficients
!^*     3.c, exponents
!^*     4.ControlMatrix
!^*     5.hmax, higher bound of feasability range for the rating curve
!^* OUT
!^*     1.k, activation stages deduced by continuity
!^*     2.ktype, type of solution of the continuity equation (see comment above)
!^*     3.feas, feasability
!^*     4.err, error code; <0:Warning, ==0:OK, >0: Error
!^*     5.mess, error message
!^**********************************************************************

real(mrk), intent(in)::b(:),a(:),c(:),hmax
integer(mik),intent(in):: ControlMatrix(:,:)
real(mrk), intent(out)::k(:),ktype(:)
logical(mlk), intent(out)::feas
integer(mik), intent(out)::err
character(*),intent(out)::mess
! locals
integer(mik)::ncontrol,i,j,nroot,troot
real(mrk)::lbound,hbound,root

err=0;mess='';feas=.true.;k=undefRN;ktype=undefIN
ncontrol=size(a)
! check feasability
if(any(a<=0._mrk) .or. any(c==0._mrk)) then
    feas=.false.;return
endif
! first control
k(1)=b(1);ktype(1)=0
if(ncontrol<=1) return ! nothing more to do
! other controls
do j=ncontrol,2,-1 ! the loop starts with the higher control
    ! determine lower/higher bounds for the current activation stage
    lbound=max(b(j),b(j-1))
    if(j==ncontrol) then;hbound=hmax;else;hbound=k(j+1);endif
    ! compute current activation stage depending on the case: addition or replacement
    if( all( ControlMatrix(j,1:(j-1)) - ControlMatrix(j-1,1:(j-1)) == 0 ) ) then ! addition
        k(j)=b(j);ktype(j)=0
    else ! replacement
        call solveContinuity(a1=a(j-1),b1=b(j-1),c1=c(j-1),&
                             a2=a(j),b2=b(j),c2=c(j),&
                             lbound=lbound,hbound=hbound,&
                             root=root,nroot=nroot,troot=troot)
        if(nroot==0) then
            feas=.false.;return
        else
            k(j)=root;ktype(j)=troot
        endif
    endif
    ! feasability checks
    if(k(j)>=hbound) then
        feas=.false.;return
    endif
    if( any( k(j)<b(1:(j-1)) ) ) then
        feas=.false.;return
    endif
enddo

end subroutine ApplyContinuity
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine SolveContinuity(a1,b1,c1,a2,b2,c2,lbound,hbound,root,nroot,troot)

!^**********************************************************************
!^* Purpose: solve the continuity equation
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^*             Modified from Valentin Mansanarez's code
!^*             With great assistance from Matteo Darienzo
!^**********************************************************************
!^* Last modified:10/01/2019
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References: Mansanarez et al. (2019) Shift happens! Adjusting
!^*             stage-discharge rating curves to riverbed morphological
!^*             changes at known times. Water Resources Research.
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*     1. a1,b1,c1: coefficients of the first control
!^*     2. a2,b2,c2: coefficients of the second control
!^*     3. lbound,hbound: lower/higher bounds between which the solution is sought
!^* OUT
!^*     1.root, solution
!^*     2.nroot, number of solutions
!^*     3.troot, type of solution (see comments in the code and WRR paper)
!^**********************************************************************
real(mrk), intent(in)::a1,b1,c1,a2,b2,c2,lbound,hbound
real(mrk), intent(out)::root
integer(mik), intent(out)::nroot,troot
!locals
real(mrk),parameter::eps=epsilon(1.)
real(mrk)::x0,fx0
real(mrk)::foo1,foo2,f1,f2
logical::exist

root=undefRN;nroot=0;troot=undefIN

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CASE I
! c1=c2, easy: either a unique explicit root, or no root at all
if(abs(c1-c2)<=eps) then
    if(abs(a1-a2)<=eps) then
        ! no solution when a1=a2 AND c1=c2 (unless b1=b2 too, but then
        ! it's twice the same control which is not acceptable either).
        return
    else ! unique and explicit solution
        foo1=a1**(1/c1)
        foo2=a2**(1/c1)
        root=(foo2*b2-foo1*b1)/(foo2-foo1)
        if(root>lbound .and. root<hbound) then ! root in bounds, all good!
            nroot=1;troot=1
        endif
        return
    endif
endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CASE II
! b1=b2, easy as well: explicit root
if(abs(b1-b2)<=eps) then
    root=b1+(a1/a2)**(1/(c2-c1))
    if(root>lbound .and. root<hbound) then ! root in bounds, all good!
        nroot=1;troot=2
    endif
    return
endif

! In all other cases, need to compute x0 (extremum abscissa)
x0=(c2*b1-c1*b2)/(c2-c1)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CASE III
! x0 <= lbound
! => f is monotonic within the bounds
! => at most 1 solution
if(x0<=lbound) then
    f1=fContinuity(lbound+eps,a1,b1,c1,a2,b2,c2)
    f2=fContinuity(hbound,a1,b1,c1,a2,b2,c2)
    if(f1*f2>0._mrk) then ! function does not change sign => no solution
        return
    else ! search for unique root with Newton
        call findRoot(a1,b1,c1,a2,b2,c2,lbound+eps,hbound,f1,f2,root,exist)
        if(exist) then
            nroot=1;troot=3
        endif
        return
    endif
endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CASES IV to VI
! x0 > lbound
! => f is not monotonic between lbound and +infinity
! => 0, 1 or 2 solutions
if(x0>lbound) then
    ! evaluate f(x0)
    fx0=fContinuity(x0,a1,b1,c1,a2,b2,c2)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! CASE V
    ! f(x0)=0 => x0 is the unique root
    if(abs(fx0)<=eps) then
        root=x0;nroot=1;troot=5
        return
    endif
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! CASE VI
    ! cases with no roots
    if( (c2>c1 .and. fx0>0._mrk) .or. (c2<c1 .and. fx0<0._mrk) ) then
        nroot=0;troot=6
        return
    endif
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! CASE IV - the double-root case.
    ! Search for root numerically right of x0
    ! (decision rule 'keep the largest', see WRR paper).
    troot=4
    if(x0>=hbound) then ! the solution is out of the rating curve range
        return
    else ! search for root numerically
        f2=fContinuity(hbound,a1,b1,c1,a2,b2,c2)
        if(f2*fx0>0._mrk) then ! function does not change sign => no solution in RC range
            return
        else ! search for root numerically
            call findRoot(a1,b1,c1,a2,b2,c2,x0,hbound,fx0,f2,root,exist)
            if(exist) then
                nroot=1
            endif
            return
        endif
    endif
endif
end subroutine SolveContinuity
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
pure function fContinuity(x,a1,b1,c1,a2,b2,c2)

!^**********************************************************************
!^* Purpose: Continuity function
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^*             Modified from Valentin Mansanarez's code
!^*             With great assistance from Matteo Darienzo
!^**********************************************************************
!^* Last modified:10/01/2019
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References: Mansanarez et al. (2019) Shift happens! Adjusting
!^*             stage-discharge rating curves to riverbed morphological
!^*             changes at known times. Water Resources Research.
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*     1. x: activation stage to be found
!^*     2. a1,b1,c1: coefficients of the first control
!^*     3. a2,b2,c2: coefficients of the second control
!^* OUT
!^*     1.fContinuity:: value of the continuity function
!^**********************************************************************
real(mrk), intent(in)::x,a1,b1,c1,a2,b2,c2
real(mrk)::fContinuity

fContinuity=c2*log(x-b2)-c1*log(x-b1)-log(a1/a2)

end function fContinuity
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
pure function dfContinuity(x,a1,b1,c1,a2,b2,c2)

!^**********************************************************************
!^* Purpose: Derivative of the continuity function
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^*             Modified from Valentin Mansanarez's code
!^*             With great assistance from Matteo Darienzo
!^**********************************************************************
!^* Last modified:10/01/2019
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References: Mansanarez et al. (2019) Shift happens! Adjusting
!^*             stage-discharge rating curves to riverbed morphological
!^*             changes at known times. Water Resources Research.
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*     1. x: activation stage to be found
!^*     2. a1,b1,c1: coefficients of the first control
!^*     3. a2,b2,c2: coefficients of the second control
!^* OUT
!^*     1.dfContinuity:: value of the derivative of the continuity function
!^**********************************************************************
real(mrk), intent(in)::x,a1,b1,c1,a2,b2,c2
real(mrk)::dfContinuity

dfContinuity= (c2/(x-b2)) - (c1/(x-b1))

end function dfContinuity
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine fNewton(dataIN,dataOUT,x,feas,fx,dfdxV,err,message)

!^**********************************************************************
!^* Purpose: Wrapper for continuity function and its derivative, complying
!^*          with the interface required by DMSL's Newton algo (znewton)
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^*             Modified from Valentin Mansanarez's code
!^**********************************************************************
!^* Last modified:10/01/2019
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References: see DMSL
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*     1. dataIN: contains a1,b1,c1,a2,b2,c2 in dataIN%rp0(1:6)
!^*     2. dataOUT: not used
!^*     3. x: where functions f and df are evaluated
!^* OUT
!^*     1.feas: feasible?
!^*     2.fx: f(x)
!^*     3.dfdxV: df(x)
!^*     4.err: err code
!^*     5.message: message
!^**********************************************************************
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
    fx=fContinuity(x=x,a1=dataIN%rp0(1),b1=dataIN%rp0(2),&
                   c1=dataIN%rp0(3),a2=dataIN%rp0(4),&
                   b2=dataIN%rp0(5),c2=dataIN%rp0(6))
endif
if(present(dfdxV))then
    dfdxV(1)=dfContinuity(x=x,a1=dataIN%rp0(1),b1=dataIN%rp0(2),&
                          c1=dataIN%rp0(3),a2=dataIN%rp0(4),&
                          b2=dataIN%rp0(5),c2=dataIN%rp0(6))
endif
end subroutine fNewton
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine findRoot(a1,b1,c1,a2,b2,c2,x1,x2,f1,f2,root,exist)

!^**********************************************************************
!^* Purpose: Find root of the continuity function using DMSL's znewton
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^*             Modified from Valentin Mansanarez's code
!^**********************************************************************
!^* Last modified:10/01/2019
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References: see DMSL
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*     1. a1,b1,c1: coefficients of the first control
!^*     2. a2,b2,c2: coefficients of the second control
!^*     3. x1,x2: lower/higher bounds between which the solution is sought
!^*     4. f1,f2: f values at x1/x2
!^* OUT
!^*     1.root, solution
!^*     2.exist, does solution exist?
!^**********************************************************************
use numerix_dmsl_kit, only:znewton
use utilities_dmsl_kit, only:cleanPointers
use types_dmsl_kit,only:data_ricz_type
real(mrk), intent(in)::a1,b1,c1,a2,b2,c2,x1,x2,f1,f2
real(mrk), intent(out)::root
logical, intent(out)::exist
!locals
real(mrk),parameter::tolX=0.0001_mrk,tolF=0.000001_mrk
integer(mik),parameter::itmax=10000
real(mrk)::froot
integer(mik)::fcalls,err
character(250)::message
type(data_ricz_type)::NewtonInfo

root=undefRN;exist=.false.

! Fill in object NewtonInfo (see type data_ricz_type in types_dmsl_kit) with a1,b1,c1,a2,b2,c2
if(associated(NewtonInfo%rp0)) nullify(NewtonInfo%rp0);allocate(NewtonInfo%rp0(6))
NewtonInfo%rp0(1)=a1
NewtonInfo%rp0(2)=b1
NewtonInfo%rp0(3)=c1
NewtonInfo%rp0(4)=a2
NewtonInfo%rp0(5)=b2
NewtonInfo%rp0(6)=c2

! Call Newton solver
call znewton(evalFunc=fNewton,dataIN=NewtonInfo,x1=x1,x2=x2,f1=f1,f2=f2,&
             tolX=tolX,tolF=tolF,xscale=1._mrk,fscale=1._mrk,itmax=itmax,&
             xroot=root,froot=froot,fcalls=fcalls,err=err,message=message)
! conclude
if(err==0) exist=.true.
! clean memory
call cleanPointers(NewtonInfo,what=1,err=err,message=message)

end subroutine findRoot
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module BaRatinBAC_model
