module SedimentTransport_model

!~**********************************************************************
!~* Purpose: Sediment transport formulas
!~**********************************************************************
!~* Programmer: Ben Renard, Irstea Lyon
!~**********************************************************************
!~* Last modified: 30/07/2015
!~**********************************************************************
!~* Comments:
!~**********************************************************************
!~* References:
!~**********************************************************************
!~* 2Do List: Improve handling of feas & error in t-loop of Sediment_Apply.
!~*           But even more generally, this sediment model is a mess
!~*           with all these options and should probably be simplified
!~*           (make several sub-models?)
!~**********************************************************************
!~* Quick description of public procedures:
!~*		1. Sediment_GetParNumber, number of parameters 
!~*		2. Sediment_Apply, apply formula
!~*		3. Sediment_XtraRead, read model formulas & options
!~**********************************************************************

use kinds_dmsl_kit ! numeric kind definitions from DMSL

implicit none
Private
public :: Sediment_GetParNumber, Sediment_Apply, Sediment_XtraRead

Character(100), parameter:: & ! Catalogue of Sediment Transport formulas
                    Sediment_MPM="MPM",& ! Meyer-Peter & Mueller
                    Sediment_Camenen="Camenen",& ! Best formula in the world, ever
                    Sediment_Nielsen="Nielsen",& ! Nielsen
                    Sediment_Smart="Smart" ! Smart

Character(100), parameter:: & ! Catalogue of Critical Shear Stress formulas
                    CSS_Constant="Constant",& ! CSS = cte
                    CSS_Soulsby="Soulsby",& ! CSS = 0.3/(1+1.2D) + 0.55*(1-exp(-0.02D))
                    CSS_SW="SoulsbyWhitehouse" ! CSS = 0.24/D + 0.55*(1-exp(-0.02D))

Character(100), parameter:: & ! Catalogue of options for pseudo-parameters d,s,v and g
                    PPO_FIXED="FIX",& ! fixed value
                    PPO_INPUT="IN",& ! input data
                    PPO_PAR="PAR" ! parameter

Character(100), parameter:: & ! Catalogue of coeff options
                    CO_FIX="FIX",& ! coeff of empirical formulas are hard-coded
                    CO_PAR="PAR" ! coeff of empirical formulas are parameters

integer(mik),parameter::nINPUT_base=2
integer(mik),parameter::npseudopar=4
integer(mik),parameter::d_ix=1,s_ix=2,v_ix=3,g_ix=4

Contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Sediment_GetParNumber(ID,CSS_ID,PPO,CO,npar,err,mess)

!^**********************************************************************
!^* Purpose: number of parameters
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified: 30/07/2015
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*		1. ID, ID of base formula
!^*		2. CSS_ID, ID of CSS formula
!^*		3. PPO, pseudo-par options
!^*		4. CO, coeff option
!^* OUT
!^*		1. npar, par. number
!^*		2.err, error code; <0:Warning, ==0:OK, >0: Error
!^*		3.mess, error message
!^**********************************************************************
character(*), intent(in)::ID,CSS_ID,PPO(npseudopar),CO
integer(mik), intent(out)::npar,err
character(*),intent(out)::mess
! locals
character(250),parameter::procname='Sediment_GetParNumber'
integer(mik)::nb,nc,np

err=0;mess='';npar=undefIN
call GetNpar_base(ID,nb,err,mess)
if(err>0) then;mess=trim(procname)//':'//trim(mess);return;endif
call GetNpar_pseudo(PPO,np,err,mess)
if(err>0) then;mess=trim(procname)//':'//trim(mess);return;endif
call GetNpar_coeff(CSS_ID,CO,nc,err,mess)
if(err>0) then;mess=trim(procname)//':'//trim(mess);return;endif
npar=nb+np+nc
if(trim(CSS_ID)==CSS_Constant)npar=npar+1
end subroutine Sediment_GetParNumber
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Sediment_Apply(ID,CSS_ID,PPO,CO,IN,theta,pseudoval,OUT,feas,err,mess)

!^**********************************************************************
!^* Purpose: apply sediment transport formula
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified:30/07/2015
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*		1. ID, ID of base formula
!^*		2. CSS_ID, ID of CSS formula
!^*		3. PPO, pseudo-par options
!^*		4. CO, coeff option
!^*		5. IN, input vector
!^*		6. theta, parameter vector
!^*		7. pseudoval, values of pseudo-parameters (used only for those treated as FIX)
!^* OUT
!^*		1. OUT, sediment discharge (DIMENSIONAL)
!^*		2. feas, feasible?
!^*		3. err, error code; <0:Warning, ==0:OK, >0: Error
!^*		4. mess, error message
!^**********************************************************************

character(*), intent(in)::ID,CSS_ID,PPO(npseudopar),CO
real(mrk), intent(in)::IN(:,:), theta(:),pseudoval(npseudopar)
real(mrk), intent(out)::OUT(:)
logical, intent(out)::feas(:)
integer(mik), intent(out)::err
character(*),intent(out)::mess
!locals
character(250),parameter::procname='Sediment_Apply'
integer(mik)::nbase,ncoeff,t,nobs
real(mrk)::cte,css,pseudo(npseudopar),d,s,v,g,tau,dim
real(mrk),allocatable::base(:),coeff(:)

err=0;mess='';feas=.true.;OUT=undefRN

nobs=size(IN,1)

! get size of base & coeff 
call GetNpar_base(ID,nbase,err,mess)
if(err>0) then;mess=trim(procname)//':'//trim(mess);return;endif
call GetNCoeff(CSS_ID,ncoeff,err,mess)
if(err>0) then;mess=trim(procname)//':'//trim(mess);return;endif
if(allocated(base)) deallocate(base);allocate(base(nbase))
if(allocated(coeff)) deallocate(coeff);allocate(coeff(ncoeff))

do t=1,nobs
    if(any(IN(t,:)<=0._mrk)) then
        feas(t)=.false.;cycle
    endif

    ! Make sense of what's in theta/IN depending on various options
    call GetAllPar(ID,CSS_ID,PPO,CO,nbase,ncoeff,theta,IN(t,:),pseudoval,&
                   cte,base,pseudo,coeff,err,mess)
    if(err>0) then;feas(t)=.false.;cycle;endif
    d=pseudo(d_ix);s=pseudo(s_ix);v=pseudo(v_ix);g=pseudo(g_ix)

    ! compute CSS
    call CSS_Apply(CSS_ID,cte,pseudo,coeff,css,feas(t),err,mess)
    if(err>0) then;feas(t)=.false.;cycle;endif
    if(.not.feas(t)) cycle

    ! compute tau_star & dimensionality coeff
    if(s<=1._mrk) then;feas(t)=.false.;cycle;endif
    tau=(IN(t,1)*IN(t,2))/((s-1._mrk)*d)
    if(tau<=0._mrk) then;feas(t)=.false.;cycle;endif
    if(tau<=css) then;OUT(t)=0._mrk;cycle;endif
    dim=sqrt(g*(s-1._mrk)*(d**3))

    !apply base formula
    select case (trim(ID))
    case(Sediment_MPM)
        OUT(t)=dim*base(1)*( (tau-css)**base(2) )
    case(Sediment_Camenen)
        OUT(t)=dim*base(1)*( tau**base(2) )*exp(-1._mrk*base(3)*(css/tau))
    case(Sediment_Nielsen)
        OUT(t)=dim*base(1)*( tau**base(2) )*(tau-css)
    case(Sediment_Smart)
        OUT(t)=dim * base(1) * (IN(t,2)**base(2)) * ((s-1._mrk)**(1._mrk/6._mrk)) *&
            &(d**(1._mrk/6._mrk)) * (g**(-0.5_mrk)) * ( tau**base(3) ) * (tau-css)
    case default
        err=1;mess=trim(procname)//':Fatal:Unavailable [ID]';return
    end select
enddo
end subroutine Sediment_Apply
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Sediment_XtraRead(file,xtra,err,mess)

!^**********************************************************************
!^* Purpose: Read Xtra information for Sediment model: formulas & options
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
!^**********************************************************************
use types_dmsl_kit, only:data_ricz_type
use utilities_dmsl_kit,only:getSpareUnit

character(*), intent(in)::file
type(data_ricz_type),intent(out):: xtra
integer(mik), intent(out)::err
character(*),intent(out)::mess
! locals
character(250),parameter::procname='Sediment_XtraRead'
integer(mik)::unt

err=0;mess=''
call getSpareUnit(unt,err,mess)
if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
open(unit=unt,file=trim(file), status='old', iostat=err)
if(err>0) then;mess=trim(procname)//':problem opening file '//trim(file);return;endif
read(unt,*,iostat=err) xtra%cs1 ! base formula ID
if(err/=0) then;mess=trim(procname)//':problem reading file '//trim(file);return;endif
read(unt,*,iostat=err) xtra%cs2 ! Critical Shear Stress formula ID
if(err/=0) then;mess=trim(procname)//':problem reading file '//trim(file);return;endif
read(unt,*,iostat=err) xtra%cs3 ! Coefficient option
if(err/=0) then;mess=trim(procname)//':problem reading file '//trim(file);return;endif
if(associated(xtra%cp1)) nullify(xtra%cp1);allocate(xtra%cp1(npseudopar))
if(associated(xtra%rp1)) nullify(xtra%rp1);allocate(xtra%rp1(npseudopar))
read(unt,*,iostat=err) xtra%cp1 ! pseudo-par options
if(err/=0) then;mess=trim(procname)//':problem reading file '//trim(file);return;endif
read(unt,*,iostat=err) xtra%rp1 ! pseudo-par values
if(err/=0) then;mess=trim(procname)//':problem reading file '//trim(file);return;endif
close(unt)
end subroutine Sediment_XtraRead
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!==============!
! Private subs !
!==============!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine CSS_Apply(ID,cte,pseudo,coeff,css,feas,err,mess)

!^**********************************************************************
!^* Purpose: Apply Critical Shear Stress formula
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified:30/07/2015
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*		1. ID, ID of the CSS formula function
!^*		2. cte, constante
!^*		3. pseudo, pseudo-parameters d,s,g,v
!^*		4. coeff, coefficients of CSS formulas
!^* OUT
!^*		1. css, critical shear stress
!^*		2. feas, feasible?
!^*		3.err, error code; <0:Warning, ==0:OK, >0: Error
!^*		4.mess, error message
!^**********************************************************************

character(*), intent(in)::ID
real(mrk), intent(in)::cte,pseudo(npseudopar),coeff(:)
real(mrk), intent(out)::css
logical, intent(out)::feas
integer(mik), intent(out)::err
character(*),intent(out)::mess
!locals
character(250),parameter::procname='CSS_Apply'
real(mrk)::SD,d,s,v,g

err=0;mess='';feas=.true.;css=undefRN

if(any(pseudo<=0._mrk)) then;feas=.false.;return;endif

select case (trim(ID))
case(CSS_Constant)
    css=cte;return
case(CSS_Soulsby,CSS_SW)
    ! Compute Sediment-Diameter SD
    d=pseudo(d_ix);s=pseudo(s_ix);v=pseudo(v_ix);g=pseudo(g_ix)
    SD=( (g*(s-1._mrk)/(v**2))**(1._mrk/3._mrk) ) *d
    ! apply formula
    select case(trim(ID))
    case(CSS_Soulsby)
        css=(coeff(1)/(1._mrk+coeff(2)*SD)) * coeff(3)*(1._mrk-exp(-1._mrk*coeff(4)*SD))
    case(CSS_SW)
        css=(coeff(1)/SD) * coeff(2)*(1._mrk-exp(-1._mrk*coeff(3)*SD))
    end select
case default
    err=1;mess=trim(procname)//':Fatal:Unavailable [ID]';return
end select

end subroutine CSS_Apply
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine CSS_GetCoeff(ID,coeff,err,mess)

!^**********************************************************************
!^* Purpose: get default coefficient for Critical Shear Stress formula
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified:30/07/2015
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*		1. ID, ID of the CSS formula function
!^* OUT
!^*		1. coeff, coefficients of CSS formula
!^*		2.err, error code; <0:Warning, ==0:OK, >0: Error
!^*		3.mess, error message
!^**********************************************************************

character(*), intent(in)::ID
real(mrk), intent(out)::coeff(:)
integer(mik), intent(out)::err
character(*),intent(out)::mess
!locals
character(250),parameter::procname='CSS_GetCoeff'

err=0;mess='';coeff=undefRN

select case (trim(ID))
case(CSS_Constant)
    return
case(CSS_Soulsby)
    coeff=(/0.3_mrk,1.2_mrk,0.055_mrk,0.02_mrk/)
case(CSS_SW)
    coeff=(/0.24_mrk,0.055_mrk,0.02_mrk/)
case default
    err=1;mess=trim(procname)//':Fatal:Unavailable [ID]';return
end select

end subroutine CSS_GetCoeff
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine GetNpar_base(ID,npar,err,mess)

!^**********************************************************************
!^* Purpose: number of base parameters
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified:30/07/2015
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*		1. ID, ID of the base formula
!^* OUT
!^*		1. npar
!^*		2.err, error code; <0:Warning, ==0:OK, >0: Error
!^*		3.mess, error message
!^**********************************************************************

character(*), intent(in)::ID
integer(mik), intent(out)::npar,err
character(*),intent(out)::mess
!locals
character(250),parameter::procname='GetNpar_base'

err=0;mess='';npar=undefIN

select case (trim(ID))
case(Sediment_MPM)
    npar=2
case(Sediment_Camenen)
    npar=3
case(Sediment_Nielsen)
    npar=2
case(Sediment_Smart)
    npar=3
case default
    err=1;mess=trim(procname)//':Fatal:Unavailable [ID]';return
end select

end subroutine GetNpar_base
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine GetNpar_pseudo(PPO,npar,err,mess)

!^**********************************************************************
!^* Purpose: number of pseudo-parameters treated as parameters
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified:30/07/2015
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*		1. PPO, pseudo-parameters options
!^* OUT
!^*		1. npar
!^*		2.err, error code; <0:Warning, ==0:OK, >0: Error
!^*		3.mess, error message
!^**********************************************************************

character(*), intent(in)::PPO(npseudopar)
integer(mik), intent(out)::npar,err
character(*),intent(out)::mess
!locals
character(250),parameter::procname='GetNpar_pseudo'

err=0;mess=''
npar=count(PPO==PPO_PAR)
end subroutine GetNpar_pseudo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine GetNpar_coeff(CSS_ID,CO,npar,err,mess)

!^**********************************************************************
!^* Purpose: number of coefficients of CSS formula treated as parameters
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified:30/07/2015
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*		1. CSS_ID, ID of CSS formula
!^*		2. CO, coefficient option
!^* OUT
!^*		1. npar
!^*		2.err, error code; <0:Warning, ==0:OK, >0: Error
!^*		3.mess, error message
!^**********************************************************************

character(*), intent(in)::CSS_ID,CO
integer(mik), intent(out)::npar,err
character(*),intent(out)::mess
!locals
character(250),parameter::procname='GetNpar_coeff'

err=0;mess='';npar=undefIN
if(CO==CO_FIX) then
    npar=0
else
    select case (trim(CSS_ID))
    case(CSS_Constant)
        npar=1
    case(CSS_Soulsby)
        npar=4
    case(CSS_SW)
        npar=3
    case default
        err=1;mess=trim(procname)//':Fatal:Unavailable [CSS_ID]';return
    end select
endif
end subroutine GetNpar_coeff
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine GetNCoeff(CSS_ID,n,err,mess)

!^**********************************************************************
!^* Purpose: number of coefficients of CSS formula 
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified:30/07/2015
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*		1. CSS_ID, ID of CSS formula
!^* OUT
!^*		1.n
!^*		2.err, error code; <0:Warning, ==0:OK, >0: Error
!^*		3.mess, error message
!^**********************************************************************

character(*), intent(in)::CSS_ID
integer(mik), intent(out)::n,err
character(*),intent(out)::mess
!locals
character(250),parameter::procname='GetNCoeff'

err=0;mess='';n=undefIN
    select case (trim(CSS_ID))
    case(CSS_Constant)
        n=0
    case(CSS_Soulsby)
        n=4
    case(CSS_SW)
        n=3
    case default
        err=1;mess=trim(procname)//':Fatal:Unavailable [CSS_ID]';return
    end select
end subroutine GetNCoeff
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine GetAllPar(ID,CSS_ID,PPO,CO,nbase,ncoeff,theta,IN,pseudoval,&
                     cte,base,pseudo,coeff,err,mess)

!^**********************************************************************
!^* Purpose: Extract all parameters
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified:30/07/2015
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*		1. ID, ID of base formula
!^*		2. CSS_ID, ID of CSS formula
!^*		3. PPO, pseudo-par options
!^*		4. CO, coeff option
!^*		5. nbase, number of base parameters
!^*		6. ncoeff, number of coefficients
!^*		7. theta, packed parameter vector
!^*		8. IN, input vector
!^*		9. pseudoval, values of pseudo-parameters (used only for those treated as FIX)
!^* OUT
!^*		1. cte, cte used if CSS_ID="Constant"
!^*		2. base, base parameters
!^*		3. pseudo, pseudo-parameters
!^*		4. coeff, CSS coefficients
!^*		5. err, error code; <0:Warning, ==0:OK, >0: Error
!^*		6. mess, error message
!^**********************************************************************

character(*), intent(in)::ID,CSS_ID,PPO(npseudopar),CO
integer(mik),intent(in)::nbase,ncoeff
real(mrk),intent(in)::theta(:),IN(:),pseudoval(npseudopar)
real(mrk),intent(out)::cte,base(nbase),pseudo(npseudopar),coeff(ncoeff)
integer(mik), intent(out)::err
character(*),intent(out)::mess
!locals
character(250),parameter::procname='GetAllpar'
integer(mik)::ktheta,kin,i

err=0;mess='';cte=undefRN;base=undefRN;pseudo=undefRN;coeff=undefRN

ktheta=0;kin=nINPUT_base
if(nbase>0) then
    base=theta((ktheta+1):(ktheta+nbase))
    ktheta=ktheta+nbase
endif
if(trim(CSS_ID)==CSS_Constant) then
    cte=theta(ktheta+1)
    ktheta=ktheta+1
endif
do i=1,npseudopar
    select case (trim(PPO(i)))
    case(PPO_FIXED)
        pseudo(i)=pseudoval(i)
    case(PPO_PAR)
        pseudo(i)=theta(ktheta+1)
        ktheta=ktheta+1
    case(PPO_INPUT)
        pseudo(i)=IN(kin+1)
        kin=kin+1
    case default
        err=1;mess=trim(procname)//':Fatal:Unavailable [PPO]';return
    end select
enddo
if(trim(CSS_ID)/=CSS_Constant) then
    select case(trim(CO))
    case(CO_PAR)
        coeff=theta((ktheta+1):)
        ktheta=ktheta+ncoeff
    case(CO_FIX)
        call CSS_GetCoeff(CSS_ID,coeff,err,mess)
        if(err>0) then;mess=trim(procname)//':'//trim(mess);return;endif
    case default
        err=1;mess=trim(procname)//':Fatal:Unavailable [CO]';return
    end select
endif
end subroutine GetAllpar
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module SedimentTransport_model
