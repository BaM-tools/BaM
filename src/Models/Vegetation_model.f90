module Vegetation_model

!~**********************************************************************
!~* Purpose: Rating curve affected by vegetation
!~**********************************************************************
!~* Programmer: Ben Renard, Irstea Lyon
!~**********************************************************************
!~* Last modified: 08/10/2019 by Emeline Perret
!~**********************************************************************
!~* Comments:
!~      - correction of the problem related to the bending
!~**********************************************************************
!~* References:
!~**********************************************************************
!~* 2Do List:
!~**********************************************************************
!~* Quick description of public procedures:
!~*     1. Vegetation_GetParNumber, number of parameters
!~*     2. Vegetation_Apply, apply formula
!~*     3. Vegetation_XtraRead, read growth submodel
!~**********************************************************************

use kinds_dmsl_kit ! numeric kind definitions from DMSL

implicit none
Private
public :: Vegetation_GetParNumber, Vegetation_Apply, Vegetation_XtraRead,Growth_Yin1_temporal,Vegetation_fNewton

Character(100), parameter:: & ! Catalogue of growth formulas
                    Growth_Yin1="Growth_Yin1",& ! parameterization based on tm (time of maximum growth rate)
                    Growth_Yin2="Growth_Yin2",& ! parameterization based on alpha (exponent of rising limb)
                    Growth_Yin3="Growth_Yin3",& ! parameterization based on tf (end of vegetation cycle)
                    Growth_Yin1_temporal="Growth_Yin1_temporal" ! same as Yin1, but based on times only (no temperature used)

type,public:: Vegetation_NewtonOptionType ! options for newton algorithm
    real(mrk)::xscale=undefRN,fscale=undefRN,tolX=undefRN,tolF=undefRN
    integer(mik)::itmax=undefIN
end type Vegetation_NewtonOptionType

integer(mik),parameter:: npar_Q0=5,& ! B, S0, n1, b1, c1
                         npar_bend=2 ! chi, U_chi

Contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Vegetation_GetParNumber(growthFormula,nYear,npar,err,mess)
!^**********************************************************************
!^* Purpose: number of parameters
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified: 09/02/2018
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*     1. growthFormula, ID of the growth formula submodel
!^*     2. nYear, number of distinct years in record
!^* OUT
!^*     1. npar, par. number
!^*     2.err, error code; <0:Warning, ==0:OK, >0: Error
!^*     3.mess, error message
!^**********************************************************************
character(*), intent(in)::growthFormula
integer(mik), intent(in)::nYear
integer(mik), intent(out)::npar,err
character(*),intent(out)::mess
! locals
character(250),parameter::procname='Vegetation_GetParNumber'
integer(mik)::nGrowth,nGmax

err=0;mess='';npar=undefIN
call GetNpar_growth(growthFormula,nGrowth,err,mess)
if(err>0) then;mess=trim(procname)//':'//trim(mess);return;endif
if(trim(growthFormula)==trim(Growth_Yin1_temporal)) then
    nGmax=1
else
    nGmax=nYear
endif
npar=npar_Q0+npar_bend+nGrowth+nGmax
end subroutine Vegetation_GetParNumber
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Vegetation_Apply(growthFormula,nYear,IN,theta,NewtonOption,OUT,Dpar,feas,err,mess)
!^**********************************************************************
!^* Purpose: apply Vegetation formula
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified: 08/10/2019 by Emeline Perret
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*     1. growthFormula
!^*     2. nYear, number of years (have to be contiguous and non-empty)
!^*     3. IN, matrix of input variables (time, h, Traw, Tsmooth)
!^*     4. theta, parameter vector
!^*     5. NewtonOption
!^* OUT
!^*     1. OUT, matrix of output variables (Q,cosg0, sing0, g0)
!^*     2. Dpar, derived parameters (B.sqrt(S0)/n + optimal temperature for each year)
!^*     3. feas, feasible?
!^*     4. err, error code; <0:Warning, ==0:OK, >0: Error
!^*     5. mess, error message
!^**********************************************************************
use types_dmsl_kit,only:data_ricz_type
use utilities_dmsl_kit, only:cleanPointers
use numerix_dmsl_kit, only:znewton

character(*), intent(in)::growthFormula
integer(mik), intent(in)::nYear
real(mrk),intent(in)::IN(:,:),theta(:)
type(Vegetation_NewtonOptionType), intent(in)::NewtonOption
real(mrk),intent(out)::OUT(:,:),Dpar(:)
logical, intent(out)::feas(:)
integer(mik), intent(out)::err
character(*),intent(out)::mess
! locals
character(250),parameter::procname='Vegetation_Apply'
real(mrk),parameter::gravity=9.81_mrk
real(mrk)::time(size(IN,dim=1)),h(size(IN,dim=1)),Traw(size(IN,dim=1)),&
           Tsmooth(size(IN,dim=1)),g(size(IN,dim=1)),g0(size(IN,dim=1)),&
           cosg(size(IN,dim=1)),sing(size(IN,dim=1)),Q(size(IN,dim=1))
real(mrk)::B,S0,n1,b1,c1,chi,U_chi,Tm(nYear)
integer(mik)::npar,nGrowthPar
real(mrk),allocatable::growthPar(:)
real(mrk),allocatable::gmax(:)
integer(mik)::t,nt,fcalls
real(mrk)::y,Q0,gamma,gamma2,froot
type(data_ricz_type)::NewtonInfo

err=0;mess='';feas=.true.;OUT=undefRN
!----------------------------
! Make sense of input matrix
nt=size(IN,dim=1)
time=IN(:,1)
h=IN(:,2)
if(growthFormula/=Growth_Yin1_temporal) then
    Traw=IN(:,3)
    Tsmooth=IN(:,4)
endif
!----------------------------
! Make sense of parameter vector
call Vegetation_GetParNumber(growthFormula,nYear,npar,err,mess)
if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
if(size(theta)/=npar) then;err=1;mess=trim(procname)//':wrong size[theta]';return;endif
! base parameters
B=theta(1)
S0=theta(2)
n1=theta(3)
b1=theta(4)
c1=theta(5)
if(S0>=0._mrk) then
    Dpar(1)=B*sqrt(S0)/n1
else
    Dpar(1)=undefRN
endif
! bending
chi=theta(npar_Q0+1)
U_chi=theta(npar_Q0+2)    ! NOTE: modification here!
! growth function parameters
call GetNpar_growth(growthFormula,nGrowthPar,err,mess)
if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
if(allocated(growthPar)) deallocate(growthPar);allocate(growthPar(nGrowthPar))
growthPar=theta( (npar_Q0+npar_bend+1):(npar_Q0+npar_bend+nGrowthPar) )
! gmax
if(allocated(gmax)) deallocate(gmax)
if(growthFormula==Growth_Yin1_temporal) then
    allocate(gmax(1))
else
    allocate(gmax(nYear))
endif
gmax=theta( (npar_Q0+npar_bend+nGrowthPar+1):(npar_Q0+npar_bend+nGrowthPar+size(gmax)) )
! Check parameter feasability
if(B<=0._mrk .or. S0<=0._mrk .or. n1<=0._mrk .or. c1<=0._mrk .or.&
   chi<-2._mrk .or. chi>0._mrk .or. U_chi<=0._mrk .or. any(gmax<0._mrk)) then
    feas=.false.;return
endif
!----------------------------
! construct g(t)
call Apply_growth(growthFormula,time,Traw,Tsmooth,growthPar,gmax,&
                  g,g0,cosg,sing,Tm,feas,err,mess)
if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
if(all(.not.feas)) return
if(growthFormula/=Growth_Yin1_temporal) then
    Dpar(2:nYear+1)=Tm
endif
OUT(:,2)=g
OUT(:,3)=cosg
OUT(:,4)=sing
OUT(:,5)=g0
!----------------------------
! compute discharge
do t=1,nt
    if(.not.feas(t)) cycle
    y=h(t)-b1
    ! Q=0
    if(y<=0._mrk) then;Q(t)=0._mrk;cycle;endif
    ! growth=0
    Q0=((B*sqrt(S0))/n1)*(y**c1)
    if(g(t)<=0._mrk) then;Q(t)=Q0;cycle;endif
    ! general case
    gamma=(g(t)*y**(1._mrk/3._mrk))/(8._mrk*gravity*(n1**2))
    gamma2=(B**chi)*(y**chi)*(U_chi**chi)   ! NOTE: modification here!
    ! Pack parameters Q0, gamma and chi into RICZ structure
    if(associated(NewtonInfo%rp0)) nullify(NewtonInfo%rp0);allocate(NewtonInfo%rp0(4)) ! NOTE: modification here!
    NewtonInfo%rp0(1)=Q0
    NewtonInfo%rp0(2)=gamma
    NewtonInfo%rp0(3)=gamma2
    NewtonInfo%rp0(4)=chi
    ! Apply numerical resolution
    call znewton(evalFunc=Vegetation_fNewton,dataIN=NewtonInfo,x1=0._mrk,x2=Q0,&
                 tolX=NewtonOption%tolX,tolF=NewtonOption%tolF,xscale=NewtonOption%xscale,&
                 fscale=NewtonOption%fscale,itmax=NewtonOption%itmax,&
                 xroot=Q(t),froot=froot,fcalls=fcalls,err=err,message=mess)
    if(err>0) then
        feas(t)=.false.
        call cleanPointers(NewtonInfo,what=1,err=err,message=mess)
        if(err>0) then;mess=trim(procname)//':'//trim(mess);return;endif
        cycle
    endif
    ! Does not find root, cycle
    if(fcalls>NewtonOption%itmax) then
        feas(t)=.false.
        call cleanPointers(NewtonInfo,what=1,err=err,message=mess)
        if(err>0) then;mess=trim(procname)//':'//trim(mess);return;endif
        cycle
    endif
    ! clean RICZ structure to avoid memory leak
    call cleanPointers(NewtonInfo,what=1,err=err,message=mess)
    if(err>0) then;mess=trim(procname)//':'//trim(mess);return;endif
enddo
! return results
OUT(:,1)=Q
end subroutine Vegetation_Apply
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Vegetation_XtraRead(file,xtra,err,mess)
!^**********************************************************************
!^* Purpose: Read Xtra information (growth formula, nYear, Newton pars)
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified:09/02/2018
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
character(250),parameter::procname='Vegetation_XtraRead'
integer(mik), parameter::nNewtonPar=4
integer(mik)::unt

err=0;mess=''
call getSpareUnit(unt,err,mess)
if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
open(unit=unt,file=trim(file), status='old', iostat=err)
if(err>0) then;mess=trim(procname)//':problem opening file '//trim(file);return;endif
read(unt,*,iostat=err) xtra%cs1 ! growth formula ID
if(err/=0) then;mess=trim(procname)//':problem reading file '//trim(file);return;endif
read(unt,*,iostat=err) xtra%is1 ! nYear
if(err/=0) then;mess=trim(procname)//':problem reading file '//trim(file);return;endif
if(associated(xtra%rp1)) nullify(xtra%rp1);allocate(xtra%rp1(nNewtonPar))
read(unt,*,iostat=err) xtra%rp1(1) ! Newton - xscale
if(err/=0) then;mess=trim(procname)//':problem reading file '//trim(file);return;endif
read(unt,*,iostat=err) xtra%rp1(2) ! Newton - fscale
if(err/=0) then;mess=trim(procname)//':problem reading file '//trim(file);return;endif
read(unt,*,iostat=err) xtra%rp1(3) ! Newton - xtol
if(err/=0) then;mess=trim(procname)//':problem reading file '//trim(file);return;endif
read(unt,*,iostat=err) xtra%rp1(4) ! Newton - ftol
if(err/=0) then;mess=trim(procname)//':problem reading file '//trim(file);return;endif
read(unt,*,iostat=err) xtra%is2    ! Newton - maxiter
if(err/=0) then;mess=trim(procname)//':problem reading file '//trim(file);return;endif
close(unt)

end subroutine Vegetation_XtraRead
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!==============!
! Private subs !
!==============!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine GetNpar_growth(growthFormula,npar,err,mess)
!^**********************************************************************
!^* Purpose: number of parameters for the growth formula
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified:09/02/2018
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*     1. growthFormula, ID of the growth formula
!^* OUT
!^*     1. npar
!^*     2.err, error code; <0:Warning, ==0:OK, >0: Error
!^*     3.mess, error message
!^**********************************************************************
character(*), intent(in)::growthFormula
integer(mik), intent(out)::npar,err
character(*),intent(out)::mess
!locals
character(250),parameter::procname='GetNpar_growth'

err=0;mess='';npar=undefIN

select case (trim(growthFormula))
case(Growth_Yin1_temporal)
    npar=3
case(Growth_Yin1,Growth_Yin2,Growth_Yin3)
    npar=4
case default
    err=1;mess=trim(procname)//':Fatal:Unavailable [growthFormula]';return
end select

end subroutine GetNpar_growth
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Apply_growth(growthFormula,time,Traw,Tsmooth,par,&
                        gmax,g,g0,cosg,sing,Tm,feas,err,mess)
!^**********************************************************************
!^* Purpose: Apply growth formula
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified:09/02/2018
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*     1. growthFormula, ID of the growth formula
!^*     2. time
!^*     3. Traw, temperature
!^*     4. Tsmooth, smoothed temperature
!^*     5. par, parameters of the adimensional growth formula
!^*     6. gmax, gmax parameters for the total growth
!^* OUT
!^*     1. g, total growth series
!^*     2. g0, adimensional growth series
!^*     3. cosg, cosine of angular g0
!^*     4. sing, sine of angular g0
!^*     5. Tm, optimal temperature of each year
!^*     6. feas, feasible?
!^*     7. err, error code; <0:Warning, ==0:OK, >0: Error
!^*     8.  mess, error message
!^**********************************************************************
use utilities_dmsl_kit, only:ifirstTrueLoc,pi
character(*), intent(in)::growthFormula
real(mrk),intent(in)::time(:),Traw(:),Tsmooth(:),par(:),gmax(:)
real(mrk),intent(out)::g(:),g0(:),cosg(:),sing(:),Tm(:)
logical, intent(out)::feas(:)
integer(mik), intent(out)::err
character(*),intent(out)::mess
! locals
character(250),parameter::procname='Apply_growth'
real(mrk)::Tmin,Tb,Te,time_b(size(gmax)),time_e(size(gmax)),&
           time_m(size(gmax)),time_f(size(gmax)),alpha(size(gmax))
integer(mik)::i,j,k,nYear,n
real(mrk)::cumul
logical::mask(size(time)),IsGnull
real(mrk),allocatable::tim(:),temp(:),cum(:),stemp(:),foog0(:)
real(mrk)::t0,ti,tmax,tf,alf,tt

err=0;mess='';feas=.true.;g=undefRN;g0=undefRN;cosg=undefRN;sing=undefRN;Tm=undefRN
if(growthFormula==Growth_Yin1_temporal) then ! easy
    n=size(time)
    t0=par(1);ti=par(2);tmax=par(3)
    if(t0<=0._mrk .or. ti<=0._mrk .or. tmax<=0._mrk .or. &
       t0>=1._mrk .or. ti>=1._mrk .or. tmax>=1._mrk .or. &
       t0>=ti .or. ti>=tmax) then
        feas=.false.;return
    endif
    alf=(tmax-t0)/(tmax-ti)
    tf=2._mrk*tmax-ti
    do i=1,n
        tt=time(i)-real(floor(time(i)),mrk)
        if(tt>=t0 .and. tt<=tf) then
            g0(i)=( ((tt-t0)/(tmax-t0))**alf ) * (tf-tt)/(tf-tmax)
        else
            g0(i)=0._mrk
        endif
        cosg(i)=cos(g0(i)*pi)
        if(tt<=tmax) then
            sing(i)=sin(g0(i)*pi)
        else
            sing(i)=-1._mrk*sin(g0(i)*pi)
        endif
    enddo
    g=gmax(1)*g0
else ! much more convoluted...
    ! make sense of parameters common to all formulas
    Tmin=par(1);Tb=par(2);Te=par(3)
    if(Te<=Tb) then;feas=.false.;return;endif
    ! Loop over all years
    nYear=size(gmax)
    do i=1,nYear
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! determine number of points in this year and allocate stuff
        mask=(floor(time)==i-1)
        n=count(mask)
        if(n==0) then
            err=1
            mess=trim(procname)//':Fatal:empty year, currently unsupported'
            return
        endif
        if(allocated(tim)) deallocate(tim)
        if(allocated(temp)) deallocate(temp)
        if(allocated(cum)) deallocate(cum)
        if(allocated(stemp)) deallocate(stemp)
        if(allocated(foog0)) deallocate(foog0)
        allocate(tim(n),temp(n),cum(n),stemp(n),foog0(n))
        tim=pack(time,mask)
        temp=pack(Traw,mask)
        stemp=pack(Tsmooth,mask)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! compute cumulated temperature
        cumul=0._mrk
        if(temp(1)>=Tmin) cumul=cumul+365._mrk*(tim(1)-real(i-1,mrk))*(temp(1)-Tmin)
        cum(1)=cumul
        do j=2,n
            if(temp(j)>=Tmin) cumul=cumul+365._mrk*(tim(j)-tim(j-1))*(temp(j)-Tmin)
            cum(j)=cumul
        enddo
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! determine times of begining and peak of vegetation cycle
        IsGnull=.false.
        k=ifirstTrueLoc(cum>=Tb)
        if(k<=n) then
            time_b(i)=tim(k)
        else
            IsGnull=.true.
        endif
        k=ifirstTrueLoc(cum>=Te)
        if(k<=n) then
            time_e(i)=tim(k)
        else
            feas=.false.;return ! Do not allow vegetation continues from one year to the next
        endif
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! compute year-specific parameters and g0
        if(IsGnull) then
            Tm(i)=undefRN;alpha(i)=undefRN;time_f(i)=undefRN;
            foog0=0._mrk
        else
            ! formula-specific handling of other parameters
            select case(growthFormula)
            case(Growth_Yin1)
                Tm(i)=par(4)
                k=ifirstTrueLoc(stemp>=Tm(i))
                if(k<=n) then
                    time_m(i)=tim(k)
                else
                    feas=.false.;return ! inflexion point has to be reached within the year
                endif
                if(time_m(i)<=time_b(i)) time_m(i)=time_b(i)
                if(time_m(i)>=time_e(i)) then;feas=.false.;return;endif
                alpha(i)=(time_e(i)-time_b(i))/(time_e(i)-time_m(i))
                time_f(i)=2._mrk*time_e(i)-time_m(i)
            case(Growth_Yin2)
                if(par(4)<=0._mrk) then;feas=.false.;return;endif
                alpha(i)=par(4)
                time_f(i)=time_e(i)+(time_e(i)-time_b(i))/alpha(i)
                time_m(i)=((alpha(i)-1._mrk)*time_f(i)+2._mrk*time_b(i))/(alpha(i)+1._mrk)
            case(Growth_Yin3)
                if(par(4)<=0._mrk) then;feas=.false.;return;endif
                time_f(i)=time_e(i)+par(4)
                alpha(i)=(time_e(i)-time_b(i))/(time_f(i)-time_e(i))
                time_m(i)=((alpha(i)-1._mrk)*time_f(i)+2._mrk*time_b(i))/(alpha(i)+1._mrk)
            case default
                err=1;mess=trim(procname)//':Fatal:Unavailable [growthFormula]';return
            end select
            ! check feasability and compute optimal temperature
            if(time_m(i)<time_b(i)) time_m(i)=time_b(i)
            if(alpha(i)<=0._mrk .or. time_m(i)>=time_e(i) .or. time_f(i)>tim(n)) then
                feas=.false.;return
            endif
            Tm(i)=stemp(ifirstTrueLoc(tim>=time_m(i)))
            ! compute g0
            do j=1,n
                if(tim(j)<=time_b(i)) then;foog0(j)=0._mrk;cycle;endif
                if(tim(j)>=time_f(i)) then;foog0(j)=0._mrk;cycle;endif
                foog0(j)=((tim(j)-time_b(i))/(time_e(i)-time_b(i)))**(alpha(i)) &
                         * (time_f(i)-tim(j))/(time_f(i)-time_e(i))
            enddo
        endif
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Fill in g, g0 and cos/sin of g0
        k=0
        do j=1,size(g0)
            if(mask(j)) then
                k=k+1
                g0(j)=foog0(k)
                g(j)=g0(j)*gmax(i)
                cosg(j)=cos(foog0(k)*pi)
                if(tim(k)<=time_e(i)) then
                    sing(j)=sin(foog0(k)*pi)
                else
                    sing(j)=-1._mrk*sin(foog0(k)*pi)
                endif
            endif
        enddo
    enddo
endif
end subroutine Apply_growth
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Vegetation_fNewton(dataIN,dataOUT,x,feas,fx,dfdxV,err,message)  ! NOTE: modification here!
! function to nullify by Newton method
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
! locals
real(mrk)::gamma,chi,Q0,gamma2

feas=.true.;err=0;message='';
if(present(fx)) fx=undefRN
if(present(dfdxV)) dfdxV=undefRN
Q0=dataIN%rp0(1)
gamma=dataIN%rp0(2)
gamma2=dataIN%rp0(3)
chi=dataIN%rp0(4)
if(present(fx))then
    fx=x**2+(gamma*x**(2._mrk)*min(1.,x**(chi)/gamma2))-Q0**2
endif
if(present(dfdxV))then
    if(min(1.,x**(chi)/gamma2)==1.) then
        dfdxV(1)=2._mrk*x*(1+gamma)
        else
        dfdxV(1)=2._mrk*x+(gamma/gamma2*(2._mrk+chi)*x**(1._mrk+chi))
    end if
endif

end subroutine Vegetation_fNewton
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module Vegetation_model
