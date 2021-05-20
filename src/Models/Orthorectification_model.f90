module Orthorectification_model

!~**********************************************************************
!~* Purpose: Orthorectification
!~**********************************************************************
!~* Programmer: Ben Renard, Irstea Lyon
!~**********************************************************************
!~* Last modified: 27/07/2015
!~**********************************************************************
!~* Comments:
!~**********************************************************************
!~* References:
!~**********************************************************************
!~* 2Do List: 
!~**********************************************************************
!~* Quick description of public procedures:
!~*		1. Ortho_GetParNumber, number of parameters 
!~*		2. Ortho_Apply, compute (c,l)=f(x,y,z|theta)
!~*		3. Ortho_XtraRead, read distortion function
!~**********************************************************************

use kinds_dmsl_kit ! numeric kind definitions from DMSL

implicit none
Private
public :: Ortho_GetParNumber,Ortho_Apply,Ortho_XtraRead

Character(100), parameter, PUBLIC:: & ! Catalogue of available distortion function
                    Disto_None="None",& ! no distortion
                    Disto_Radial="Radial" ! radial distortion

Contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
pure subroutine Ortho_GetParNumber(Disto_npar,npar,err,mess)

!^**********************************************************************
!^* Purpose: number of parameters
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified: 27/07/2015
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*		1. Disto_npar, number of parameters the distortion function
!^* OUT
!^*		1. npar, par. number
!^*		2.err, error code; <0:Warning, ==0:OK, >0: Error
!^*		3.mess, error message
!^**********************************************************************

integer(mik), intent(in)::Disto_npar
integer(mik), intent(out)::npar,err
character(*),intent(out)::mess
! locals

err=0;mess='';npar=11+Disto_npar

end subroutine Ortho_GetParNumber
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Ortho_Apply(Disto_ID,Disto_npar,IN,theta,OUT,Dpar,feas,err,mess)

!^**********************************************************************
!^* Purpose: Orthorectify
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified:27/07/2015
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*		1. Disto_ID, ID of the distortion function
!^*		2. Disto_npar, number of parameters of the distortion function
!^*		3. IN, real-world coordinate (x,y,z)
!^*		4. theta, parameters
!^* OUT
!^*		1. OUT, image coordinate (c,l)
!^*		2. Dpar, derived parameters (rotation matrix + ortho parameters)
!^*		3. feas, feasible?
!^*		4.err, error code; <0:Warning, ==0:OK, >0: Error
!^*		5.mess, error message
!^**********************************************************************

character(*), intent(in)::Disto_ID
integer(mik), intent(in)::Disto_npar
real(mrk), intent(in)::IN(:,:), theta(:)
real(mrk), intent(out)::OUT(:,:),Dpar(:)
logical, intent(out)::feas(:)
integer(mik), intent(out)::err
character(*),intent(out)::mess
!locals
character(250),parameter::procname='Ortho_Apply'
real(mrk)::x0,y0,z0,azim,roll,tilt ! intrinsic parameters
real(mrk)::f,c0,l0,lambda,k ! extrinsic parameters
real(mrk)::r11,r12,r13,r21,r22,r23,r31,r32,r33 ! rotation matrix
real(mrk)::a1,a2,a3,a4,b1,b2,b3,b4,c1,c2,c3 ! ortho parameters
real(mrk)::pc,pl ! pixel sizes
real(mrk)::x,y,z,c,l ! real-world / image coordinate
real(mrk)::denom,ac,al,alpha
integer(mik)::i,nobs
!real(mrk)::M(3,3),check(3,3)

err=0;mess='';feas=.true.;OUT=undefRN;Dpar=undefRN

! xtract parameters
x0=theta(1)
y0=theta(2)
z0=theta(3)
azim=theta(4)
roll=theta(5)
tilt=theta(6)
f=theta(7)
c0=theta(8)
l0=theta(9)
lambda=theta(10)
k=theta(11)
pc=lambda
pl=k*pc
if(pc<=0._mrk .or. pl<=0._mrk .or. k<=0._mrk) then;feas=.false.;return;endif
ac=f/pc
al=f/pl

! rotation matrix
r11=cos(azim)*cos(roll)+sin(azim)*cos(tilt)*sin(roll)
r12=-1._mrk*sin(azim)*cos(roll)+cos(azim)*cos(tilt)*sin(roll)
r13=sin(tilt)*sin(roll)
r21=-1._mrk*cos(azim)*sin(roll)+sin(azim)*cos(tilt)*cos(roll)
r22=sin(azim)*sin(roll)+cos(azim)*cos(tilt)*cos(roll)
r23=sin(tilt)*cos(roll)
r31=sin(azim)*sin(tilt)
r32=cos(azim)*sin(tilt)
r33=-1._mrk*cos(tilt)
! ortho-pars
alpha=r31*x0+r32*y0+r33*z0
if(alpha==0) then
    feas=.false.;return
else
    a1=(ac*r21-r31*c0)/alpha
    a2=(ac*r22-r32*c0)/alpha
    a3=(ac*r23-r33*c0)/alpha
    a4=-1._mrk*(a1*x0+a2*y0+a3*z0)
    b1=(al*r11-r31*l0)/alpha
    b2=(al*r12-r32*l0)/alpha
    b3=(al*r13-r33*l0)/alpha
    b4=-1._mrk*(b1*x0+b2*y0+b3*z0)
    c1=-1._mrk*r31/alpha
    c2=-1._mrk*r32/alpha
    c3=-1._mrk*r33/alpha
endif
! Assemble derived parameters
Dpar=(/r11,r12,r13,r21,r22,r23,r31,r32,r33,a1,a2,a3,a4,b1,b2,b3,b4,c1,c2,c3/)

! Compute image coordinate
nobs=size(IN,dim=1)
do i=1,nobs
    x=IN(i,1);y=IN(i,2);z=IN(i,3)
    denom=r31*(x-x0)+r32*(y-y0)+r33*(z-z0)
    if(denom==0._mrk) then;feas(i)=.false.;cycle;endif
    c=c0-ac*(r21*(x-x0)+r22*(y-y0)+r23*(z-z0))/denom
    l=l0-al*(r11*(x-x0)+r12*(y-y0)+r13*(z-z0))/denom
    ! Apply distortion
    call Disto_Apply(Disto_ID,Disto_npar,(/c,l/),theta,OUT(i,:),feas(i),err,mess)
    if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
enddo
end subroutine Ortho_Apply
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Ortho_XtraRead(file,xtra,err,mess)

!^**********************************************************************
!^* Purpose: Read Xtra information for Othorectif model: distortion
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
!^*		1. file, Xtra file
!^* OUT
!^*		1. xtra, xtra information
!^*		2. err, error code; <0:Warning, ==0:OK, >0: Error
!^*		3. mess, error message
!^**********************************************************************use types_dmsl_kit, only:data_ricz_type
use utilities_dmsl_kit,only:getSpareUnit
use types_dmsl_kit, only:data_ricz_type

character(*), intent(in)::file
type(data_ricz_type),intent(out):: xtra
integer(mik), intent(out)::err
character(*),intent(out)::mess
! locals
character(250),parameter::procname='Ortho_XtraRead'
integer(mik)::unt

err=0;mess=''
call getSpareUnit(unt,err,mess)
if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
open(unit=unt,file=trim(file), status='old', iostat=err)
if(err>0) then;mess=trim(procname)//':problem opening file '//trim(file);return;endif
read(unt,*,iostat=err) xtra%cs1 ! distortion ID
if(err/=0) then;mess=trim(procname)//':problem reading file '//trim(file);return;endif
read(unt,*,iostat=err) xtra%is1 ! distortion npar
if(err/=0) then;mess=trim(procname)//':problem reading file '//trim(file);return;endif
end subroutine Ortho_XtraRead
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!==============!
! Private subs !
!==============!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Disto_Apply(Disto_ID,Disto_npar,IN,theta,OUT,feas,err,mess)

!^**********************************************************************
!^* Purpose: Apply distortion function
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified:27/07/2015
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*		1. Disto_ID, ID of the distortion function
!^*		2. Disto_npar, number of parameters of the distortion function
!^*		2. IN, (c',l')
!^*		3. theta, parameters
!^* OUT
!^*		1. OUT, image coordinate (c,l)
!^*		2. feas, feasible?
!^*		3.err, error code; <0:Warning, ==0:OK, >0: Error
!^*		4.mess, error message
!^**********************************************************************

character(*), intent(in)::Disto_ID
integer(mik), intent(in)::Disto_npar
real(mrk), intent(in)::IN(2), theta(:)
real(mrk), intent(out)::OUT(2)
logical, intent(out)::feas
integer(mik), intent(out)::err
character(*),intent(out)::mess
!locals
character(250),parameter::procname='Disto_Apply'
real(mrk)::c0,l0,r,d
integer(mik)::i

err=0;mess='';feas=.true.;OUT=undefRN

select case (Disto_ID)
case(Disto_NONE)
    OUT=IN
case(Disto_Radial)
    c0=theta(8)
    l0=theta(9)
    r=sqrt( (IN(1)-c0)**2 + (IN(2)-l0)**2)
    d=1._mrk
    do i=1,Disto_npar
        d=d+theta(11+i)*r**(2*i)
    enddo
    OUT(1)=c0+d*(IN(1)-c0)
    OUT(2)=l0+d*(IN(2)-l0)
case default
    err=1;mess=trim(procname)//':Fatal:Unavailable [Disto_ID]';return
end select

end subroutine Disto_Apply
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module Orthorectification_model