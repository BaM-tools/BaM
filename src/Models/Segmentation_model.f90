module Segmentation_model
!~**********************************************************************
!~* Purpose: Segmentation of a time series (known the number of segments)
!~**********************************************************************
!^* Programmer: Ben Renard & Matteo Darienzo, Irstea Lyon
!^**********************************************************************
!^* Last modified: 23/10/2018
!~**********************************************************************


!~**********************************************************************
!~* References:
!~**********************************************************************
!~* 2Do List:
!~**********************************************************************
!~* Quick description of public procedures:

!~**********************************************************************

use kinds_dmsl_kit ! numeric kind definitions from DMSL
implicit none
Private
public :: Segmentation_GetParNumber, Segmentation_Apply, Segmentation_XtraRead
Contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
pure subroutine Segmentation_GetParNumber(nS,npar,err,mess)
!^**********************************************************************
!^* Purpose: number of parameters of the Segmentation model
!^**********************************************************************
!^* Programmer: Ben Renard & Matteo Darienzo, Irstea Lyon
!^**********************************************************************
!^* Last modified: 23/10/2018
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1. nS, number of segments
!^*    2. tmin, tmin between segments
!^* OUT
!^*    1. npar, par. number
!^*    2.err, error code; <0:Warning, ==0:OK, >0: Error
!^*    3.mess, error message
!^**********************************************************************
integer(mik), intent(in)::nS
integer(mik), intent(out)::npar,err
character(*),intent(out)::mess
! locals
err=0;mess=''
npar= nS + nS-1
end subroutine Segmentation_GetParNumber





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Segmentation_Apply(time,nS,tmin,theta,Y,feas,err,mess)
!^**********************************************************************
!^* Purpose: apply the
!^**********************************************************************
!^* Programmer: Ben Renard & Matteo Darienzo, Irstea Lyon
!^**********************************************************************
!^* Last modified: 23/10/2018
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*    1. nS, number of sources
!
!
!^* OUT
!^*    1. Y, output vector
!^*    2.
!^*    3.feas, feasibility (all theta between 0 and 1)
!^*    4.err, error code; <0:Warning, ==0:OK, >0: Error
!^*    5.mess, error message
!^**********************************************************************
integer(mik), intent(in)::nS
real(mrk), intent(in):: tmin, time(:), theta(:)
real(mrk), intent(out)::Y(:)
logical, intent(out)::feas
integer(mik), intent(out)::err
character(*),intent(out)::mess
!locals
integer(mik)::i,nt,indx
real(mrk)::tau(nS-1),mu(nS), taumin, taumax

err=0;
mess='';
Y=undefRN
feas = .TRUE.
!if(a(1)<=0._mrk .or. c(1)==0._mrk) then
!    feas=.false.;return
!endif

if (size(theta) /= (2*nS -1)) then
    feas =.FALSE.
    return
end if

mu =theta(1:nS)

if (nS > 1) then
   tau =theta((nS+1):)
endif

nt=size(time)

if (nS >= 3) then
  do i= 1, (nS-2)
    if (tau(i+1) <= tau(i)) then
        feas = .FALSE.
        return
    end if
    if ((tau(i+1) - tau(i))< tmin) then
        feas = .FALSE.
        return
    end if
    enddo
endif



if (nS > 1) then
    taumin = minval(time)
    taumax = maxval(time)
   if (any(tau <= (taumin+tmin)) .or. any(tau >= (taumax-tmin))) then
    feas = .FALSE.
    return
   end if
endif




do i=1,nt
    if (nS==1) then
        Y(i) = mu(1)
    else
        indx=GettRange(time(i),tau)
        Y(i) = mu(indx)
    endif
enddo




end subroutine Segmentation_Apply





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Segmentation_XtraRead(file,xtra,err,mess)
!^**********************************************************************
!^* Purpose: Read Xtra information (nS,tmin)
!^**********************************************************************
!^* Programmer: Ben Renard & Matteo Darienzo, Irstea Lyon
!^**********************************************************************
!^* Last modified: 23/10/2018
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
character(250),parameter::procname='Segmentation_XtraRead'
integer(mik)::unt
err=0;mess=''
call getSpareUnit(unt,err,mess)
if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
open(unit=unt,file=trim(file), status='old', iostat=err)
if(err>0) then;mess=trim(procname)//':problem opening file '//trim(file);return;endif
read(unt,*,iostat=err) xtra%is1 ! nS  number of segments
if(err/=0) then;mess=trim(procname)//':problem reading file '//trim(file);return;endif
read(unt,*,iostat=err) xtra%rs1 ! tmin between two segments
if(err/=0) then;mess=trim(procname)//':problem reading file '//trim(file);return;endif
close(unt)
end subroutine Segmentation_XtraRead
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



function GettRange(t,tau)
!#**********************************************************************
!#* Purpose: find the range where H belongs
!#**********************************************************************
!^* Programmer: Ben Renard & Matteo Darienzo, Irstea Lyon
!^**********************************************************************
!^* Last modified: 23/10/2018
!#**********************************************************************
!#* Comments:
!#**********************************************************************
!#* References:
!#**********************************************************************
!#* 2Do List:
!#**********************************************************************
!#* IN
!#*		1. t
!#*		2. k
!#* OUT
!#*		1.General_GetRange
!#**********************************************************************
real(mrk), intent(in)::t,tau(:)
integer(mik)::GettRange
!locals
integer(mik)::ntau,i
ntau=size(tau)
do i=1,ntau
    If(t<tau(i)) then
        GettRange=i;return
    else
        cycle
    endif
enddo
GettRange=ntau+1

end function GettRange
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


end module Segmentation_model
