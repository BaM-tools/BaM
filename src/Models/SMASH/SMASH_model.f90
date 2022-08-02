module SMASH_model

!~**********************************************************************
!~* Purpose: Interface to SMASH model
!~**********************************************************************
!~* Programmer: Ben Renard, INRAE Aix
!~**********************************************************************
!~* Last modified: 02/08/2022
!~**********************************************************************
!~* Comments:
!~**********************************************************************
!~* References:
!~**********************************************************************
!~* 2Do List:
!~**********************************************************************
!~* Quick description of public procedures:
!~*     1. SMASH_GetParNumber, number of inferred parameters
!~*     2. SMASH_Apply, apply model
!~*     3. SMASH_XtraRead, read Xtra model information
!~**********************************************************************

use kinds_dmsl_kit ! numeric kind definitions from DMSL

implicit none
Private
public :: SMASH_GetParNumber, SMASH_Apply, SMASH_XtraRead

Contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine SMASH_GetParNumber(npar,err,mess)
!^**********************************************************************
!^* Purpose: number of parameters
!^**********************************************************************
!^* Programmer: Ben Renard, INRAE Aix
!^**********************************************************************
!^* Last modified: 02/08/2022
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*     Nothing
!^* OUT
!^*     1. npar, par. number
!^*     2.err, error code; <0:Warning, ==0:OK, >0: Error
!^*     3.mess, error message
!^**********************************************************************
integer(mik), intent(out)::npar,err
character(*),intent(out)::mess
! locals
character(250),parameter::procname='SMASH_GetParNumber'

err=0;mess='';npar=undefIN

end subroutine SMASH_GetParNumber

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine SMASH_Apply(X,& ! inputs:
                       theta,& ! parameters
                       Y,& ! outputs
                       feas,err,mess)
!^**********************************************************************
!^* Purpose: Apply SMASH
!^**********************************************************************
!^* Programmer: Ben Renard, INRAE Aix
!^**********************************************************************
!^* Last modified:
!^* - 02/08/2022
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*     1. XXX
!^* OUT
!^*     1. XXX
!^*     3. feas, feasible?
!^*     4.err, error code; <0:Warning, ==0:OK, >0: Error
!^*     5.mess, error message
!^**********************************************************************
real(mrk), intent(in)::X(:,:),theta(:)
real(mrk), intent(out)::Y(:,:)
logical, intent(out)::feas(:)
integer(mik), intent(out)::err
character(*),intent(out)::mess
!locals
character(250),parameter::procname='SMASH_Apply'

! Init
err=0;mess='';feas=.true.;Y=undefRN

end subroutine SMASH_Apply

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine SMASH_XtraRead(file,xtra,err,mess)
!^**********************************************************************
!^* Purpose: Read Xtra information
!^**********************************************************************
!^* Programmer: Ben Renard, INRAE Aix
!^**********************************************************************
!^* Last modified: 02/08/2022
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
character(250),parameter::procname='SMASH_XtraRead'
integer(mik)::unt

err=0;mess=''
call getSpareUnit(unt,err,mess)
if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
open(unit=unt,file=trim(file), status='old', iostat=err)
if(err>0) then;mess=trim(procname)//':problem opening file '//trim(file);return;endif
read(unt,*,iostat=err) ! xtra%rp1(1) ! Newton - xscale
if(err/=0) then;mess=trim(procname)//':problem reading file '//trim(file);return;endif
close(unt)

end subroutine SMASH_XtraRead
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!==============!
! Private subs !
!==============!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module SMASH_model
