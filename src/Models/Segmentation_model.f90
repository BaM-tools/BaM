module Segmentation_model
!~**********************************************************************
!~* Purpose: Segmentation of a time series (known the number of segments)
!~**********************************************************************
!^* Programmer: Ben Renard & Matteo Darienzo, Irstea Lyon
!^**********************************************************************
!^* Last modified: 23/10/2019 - allow using inter-shift durations
!^*                rather than shift times + both tmin and nmin constraints
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
subroutine Segmentation_Apply(time,nS,tmin,nmin,theta,option,Y,feas,err,mess)
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
!^*    1. time
!^*    2. nS, number of segments
!^*    3. tmin, minimum segment duration
!^*    4. nmin, minimum number of points per segment
!^*    5. theta, parameter vector
!^*    6. option. 1: shift times, 2: inter-shift durations
!^* OUT
!^*    1. Y, output vector
!^*    2.feas, feasibility (all theta between 0 and 1)
!^*    3.err, error code; <0:Warning, ==0:OK, >0: Error
!^*    4.mess, error message
!^**********************************************************************
integer(mik), intent(in)::nS,nmin,option
real(mrk), intent(in):: tmin, time(:), theta(:)
real(mrk), intent(out)::Y(:)
logical, intent(out)::feas
integer(mik), intent(out)::err
character(*),intent(out)::mess
!locals
character(250),parameter::procname='Segmentation_Apply'
integer(mik)::i,nt,indx
real(mrk)::duration(nS-1),tau(nS-1),mu(nS)
integer(mik)::segment(size(time))

err=0;
mess='';
Y=undefRN
feas = .TRUE.
nt=size(time)

!size check
if (size(theta) /= (2*nS -1)) then
    feas =.FALSE.
    return
end if

! Single-segment case is trivial, treat it and get out
if(nS==1) then;Y=theta(1);return;endif

! From now on, multi-segment case
! Get means
mu =theta(1:nS)
! Get inter-shift durations and shift times
select case (option)
case(1) ! parameters are shift times
    tau =theta((nS+1):)
    duration(1)=tau(1)
    if(nS>2) then
        do i=2,nS-1
            duration(i)=tau(i)-tau(i+1)
        enddo
    endif
case(2) ! parameters are inter-shift durations
    duration =theta((nS+1):)
    do i=1,nS-1
        tau(i)=sum(duration(1:i))
    enddo
case default
    err=1;mess=trim(procname)//': Fatal: Unavailable [option]';return
end select

! Feasability check
if(any(duration<=0.)) then;feas=.false.;return;endif
if(any(duration<=tmin)) then;feas=.false.;return;endif
if(any(tau>=maxval(time))) then;feas=.false.;return;endif

! Determine segment of each point
do i=1,nt
    segment(i)=GettRange(time(i),tau)
enddo

! Check that there are at least nmin values for each segment
do i=1,nS
    if(count(segment==i)<nmin) then;feas=.false.;return;endif
enddo

! Get mean of each segment
do i=1,nt
    Y(i) = mu(segment(i))
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
!^*     1. file, Xtra file
!^* OUT
!^*     1. xtra, xtra information
!^*     2. err, error code; <0:Warning, ==0:OK, >0: Error
!^*     3. mess, error message
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
read(unt,*,iostat=err) xtra%is2 ! nmin minimum number of values per segment
if(err/=0) then;mess=trim(procname)//':problem reading file '//trim(file);return;endif
read(unt,*,iostat=err) xtra%is3 ! option
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
!#*     1. t
!#*     2. k
!#* OUT
!#*     1.General_GetRange
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
