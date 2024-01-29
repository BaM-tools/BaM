module HydraulicControl_section_model

!~**********************************************************************
!~* Purpose: General equation Q=f(H) for a section control, based on bathymetry
!~**********************************************************************
!~* Programmer: Ben Renard, INRAE
!~**********************************************************************
!~* Last modified: 29/01/2024
!~**********************************************************************
!~* Comments:
!~**********************************************************************
!~* References:
!~**********************************************************************
!~* 2Do List:
!~**********************************************************************
!~* Quick description of public procedures:
!~*     1. HydraulicControl_section_GetParNumber, number of parameters of the RC
!~*     2. HydraulicControl_section_Apply, compute Q=f(H|theta)
!~*     3. HydraulicControl_section_XtraRead, read h - A(h) - w(h) matrix
!~**********************************************************************

use kinds_dmsl_kit ! numeric kind definitions from DMSL

implicit none
Private
public :: HydraulicControl_section_GetParNumber,HydraulicControl_section_Apply,&
          HydraulicControl_section_XtraRead,HydraulicControl_section_Setup

Contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
pure subroutine HydraulicControl_section_GetParNumber(npar,err,mess)

!^**********************************************************************
!^* Purpose: number of parameters of the RC
!^**********************************************************************
!^* Programmer: Ben Renard, INRAE
!^**********************************************************************
!^* Last modified: 29/01/2024
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

err=0;mess='';npar=undefIN

npar=2

end subroutine HydraulicControl_section_GetParNumber
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine HydraulicControl_section_Apply(H,theta,hAwMatrix,Q,feas,err,mess)

!^**********************************************************************
!^* Purpose:
!^**********************************************************************
!^* Programmer: Ben Renard, INRAE
!^**********************************************************************
!^* Last modified:
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*     1. H, water stage
!^*     2. theta, parameters
!^*     3. hAwMatrix, Matrix h - A(h) - w(h)
!^* OUT
!^*     1. Q, RC-computed runoff
!^*     2. feas, feasible?
!^*     3.err, error code; <0:Warning, ==0:OK, >0: Error
!^*     4.mess, error message
!^**********************************************************************

real(mrk), intent(in)::H(:), theta(:)
real(mrk), intent(in)::hAwMatrix(:,:)
real(mrk), intent(out)::Q(:)
logical, intent(out)::feas(:)
integer(mik), intent(out)::err
character(*),intent(out)::mess
!locals
character(250),parameter::procname='HydraulicControl_section_Apply'
real(mrk)::coeff,g,x1,x2,y1,y2,hc,Ahc
real(mrk)::aow(size(hAwMatrix,dim=1)),lhs(size(hAwMatrix,dim=1)),fx(size(hAwMatrix,dim=1))
integer(mik)::i,n,p,k1,k2

err=0;mess='';Q=undefRN;feas=.true.

! Make sense of theta
coeff=theta(1)    ! coefficient for bathy multiplicative errors (1 if no bathy error)
g=theta(2)        ! gravity acceleration

! Check feasability
if(coeff<=0._mrk .or. g<=0._mrk) then
    feas=.false.;return
endif

n=size(H)
p=size(hAwMatrix,dim=1)
! A over w
aow(1)=0._mrk
aow(2:p)=hAwMatrix(2:p,2)/hAwMatrix(2:p,3)
! Left-hand site of critical stage equation
lhs=hAwMatrix(:,1) + 0.5_mrk*aow

do i=1,n
    if(H(i)<=hAwMatrix(1,1)) then
        Q(i)=0._mrk
    else if (H(i)>=hAwMatrix(p,1)) then
        feas(i)=.false.
    else
        ! solve critical stage equation
        fx=lhs-H(i)
        k1=maxloc(fx,dim=1,mask=fx<=0._mrk)
        k2=k1+1
        if( k1<=0 .or. k2>p) then
            feas(i)=.false.
        else
            x1=hAwMatrix(k1,1);x2=hAwMatrix(k2,1)
            y1=fx(k1);y2=fx(k2)
            hc=(x1*y2-x2*y1)/(y2-y1)
            ! Get corresponding A(hc)
            y1=hAwMatrix(k1,2);y2=hAwMatrix(k2,2)
            Ahc=(y1*(x2-hc)+y2*(hc-x1))/(x2-x1)
            ! get Q
            Q(i)=coeff*sqrt(2._mrk*g*(H(i)-hc))*Ahc
        endif
    endif
enddo


end subroutine HydraulicControl_section_Apply
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine HydraulicControl_section_XtraRead(file,xtra,err,mess)

!^**********************************************************************
!^* Purpose: Read Xtra information for RC: h - A(h) - w(h) matrix
!^**********************************************************************
!^* Programmer: Ben Renard, INRAE
!^**********************************************************************
!^* Last modified:29/01/2024
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
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
use utilities_dmsl_kit,only:getSpareUnit,getNumItemsInFile
use types_dmsl_kit, only:data_ricz_type
character(*), intent(in)::file
type(data_ricz_type),intent(out):: xtra
integer(mik), intent(out)::err
character(*),intent(out)::mess
! locals
character(250),parameter::procname='HydraulicControl_section_XtraRead'
character(250)::foo
integer(mik)::unt,i,n

err=0;mess=''
call getSpareUnit(unt,err,mess)
if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
open(unit=unt,file=trim(file), status='old', iostat=err)
if(err>0) then;mess=trim(procname)//':problem opening file '//trim(file);return;endif
call getNumItemsInFile(unt=unt,preRwnd=.false.,nskip=0,nitems=n,postPos=0,&
                       jchar=foo,err=err,message=mess)
if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
if(associated(xtra%rpm1)) nullify(xtra%rpm1);allocate(xtra%rpm1(n,3))
do i=1,n
    read(unt,*,iostat=err) xtra%rpm1(i,:)
    if(err/=0) then;mess=trim(procname)//':problem reading file '//trim(file);return;endif
enddo
close(unt)

end subroutine HydraulicControl_section_XtraRead

subroutine HydraulicControl_section_Setup(hAwMatrix,err,mess)

!^**********************************************************************
!^* Purpose: Check h - A(h) - w(h) matrix is valid
!^**********************************************************************
!^* Programmer: Ben Renard, INRAE
!^**********************************************************************
!^* Last modified:29/01/2024
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*     1. hAwMatrix, Xtra file
!^* OUT
!^*     1. err, error code; <0:Warning, ==0:OK, >0: Error
!^*     2. mess, error message
!^**********************************************************************
real(mrk), intent(in)::hAwMatrix(:,:)
integer(mik), intent(out)::err
character(*),intent(out)::mess
! locals
character(250),parameter::procname='HydraulicControl_section_Setup'
integer(mik)::n
real(mrk)::aow(size(hAwMatrix,dim=1))

err=0;mess=''
n=size(hAwMatrix,dim=1)
if(n<3) then
    err=-1
    mess=trim(procname)//':invalid h-A-w matrix: at least 3 points are required'
    return
endif

if(size(hAwMatrix,dim=2)<3) then
    err=-1
    mess=trim(procname)//':invalid h-A-w matrix: at least 3 columns are required'
    return
endif

if(any(hAwMatrix(2:n,1)-hAwMatrix(1:(n-1),1)<0)) then
    err=-1
    mess=trim(procname)//':invalid h-A-w matrix: h should be increasing'
    return
endif

if(any(hAwMatrix(2:n,2)<=0)) then
    err=-1
    mess=trim(procname)//':invalid h-A-w matrix: A should be positive'
    return
endif

if(any(hAwMatrix(2:n,3)<=0)) then
    err=-1
    mess=trim(procname)//':invalid h-A-w matrix: w should be positive'
    return
endif

if(any(hAwMatrix(2:n,2)-hAwMatrix(1:(n-1),2)<0)) then
    err=-1
    mess=trim(procname)//':invalid h-A-w matrix: A should be increasing'
    return
endif

if(any(hAwMatrix(2:n,3)-hAwMatrix(1:(n-1),3)<0)) then
    err=-1
    mess=trim(procname)//':invalid h-A-w matrix: w should be increasing'
    return
endif

aow(1)=0._mrk
aow(2:n)=hAwMatrix(2:n,2)/hAwMatrix(2:n,3)
if(any(aow(2:n)-aow(1:(n-1))<0)) then
    err=-1
    mess=trim(procname)//':invalid h-A-w matrix: A/w should be increasing'
    return
endif

end subroutine HydraulicControl_section_Setup

end module HydraulicControl_section_model
