module Vegetation_SH_model

!~**********************************************************************
!~* Purpose: Rating curve affected by vegetation
!~**********************************************************************
!~* Programmer: Martin Morlot, INRAE Lyon
!~**********************************************************************
!~* Comments:
!~      - Comes back to stability before the second control
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
   public :: Vegetation_SH_GetParNumber,Vegetation_SH_XtraRead,Vegetation_SH_Apply

   integer(mik),parameter:: npar_Q0=5 ! B, S0, n, b, c
   ! integer(mik),parameter:: npar_SH=3 ! lambda0, lambda1, alpha
   integer(mik),parameter:: npar_SH=1 ! lambda0, (lambda1 and alpha are static)

   Character(100), parameter:: & ! Catalogue of growth formulas
      Growth_Yin1="Growth_Yin1",& ! parameterization based on tm (time of maximum growth rate)
      Growth_Yin2="Growth_Yin2",& ! parameterization based on alpha (exponent of rising limb)
      Growth_Yin3="Growth_Yin3",& ! parameterization based on tf (end of vegetation cycle)
      Growth_Yin1_temporal="Growth_Yin1_temporal" ! same as Yin1, but based on times only (no temperature used)


Contains

   subroutine Vegetation_SH_GetParNumber(growthFormula,nYear,ControlMatrix,npar,err,mess)
!^**********************************************************************
!^* Purpose: number of parameters
!^**********************************************************************
!^* Programmer: Martin Morlot, INRAE Lyon
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
!^*     3. ControlMatrix
!^* OUT
!^*     1. npar, par. number
!^*     2.err, error code; <0:Warning, ==0:OK, >0: Error
!^*     3.mess, error message
!^**********************************************************************
      use Vegetation_model, only:GetNpar_growth
      character(*), intent(in)::growthFormula
      integer(mik), intent(in)::nYear
      integer(mik), intent(in)::ControlMatrix(:,:)

      integer(mik), intent(out)::npar,err
      character(*),intent(out)::mess
! locals
      character(250),parameter::procname='Vegetation_SH_GetParNumber'
      integer(mik)::nGrowth,nGmax,ncontrol

      err=0;mess='';npar=undefIN
      call GetNpar_growth(growthFormula,nGrowth,err,mess)
      if(err>0) then;mess=trim(procname)//':'//trim(mess);return;endif
      if(trim(growthFormula)==trim(Growth_Yin1_temporal)) then
         nGmax=1
      else
         nGmax=nYear
      endif

      ncontrol=size(ControlMatrix,1) ! it is assumed that the size of ControlMatrix has been checked outside of this sub, one control varies and is not a fixed parameter, so it is excluded.
      if(ncontrol<1) then
         err=1;mess=trim(procname)//':ControlMatrix too small';return
      endif


      npar=npar_Q0+3*(ncontrol-1)+nGrowth+nGmax+npar_SH
   end subroutine Vegetation_SH_GetParNumber



   subroutine Vegetation_SH_XtraRead(file,xtra,err,mess)
!^**********************************************************************
!^* Purpose: Read Xtra information (growth formula, nYear, Matrix)
!^**********************************************************************
!^* Programmer: Martin Morlot, INRAE Lyon
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
      use utilities_dmsl_kit,only:getSpareUnit,getNumItemsInFile
      use BaRatin_model, only:CheckControlMatrix
      character(*), intent(in)::file
      type(data_ricz_type),intent(out):: xtra
      integer(mik), intent(out)::err
      character(*),intent(out)::mess
! locals
      character(250),parameter::procname='Vegetation_SH_XtraRead'
      character(250)::foo
      integer(mik)::unt,i,n
      logical::feas

      err=0;mess=''
      call getSpareUnit(unt,err,mess)
      if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
      open(unit=unt,file=trim(file),status='old',iostat=err)
      if(err/=0) then;mess=trim(procname)//':problem opening file '//trim(file);return;endif

! line 1: growth formula ID
      read(unt,*,iostat=err) xtra%cs1
      if(err/=0) then;mess=trim(procname)//':problem reading growth formula (line 1) of '//trim(file);close(unt);return;endif
! line 2: nYear
      read(unt,*,iostat=err) xtra%is1
      if(err/=0) then;mess=trim(procname)//':problem reading nYear (line 2) of '//trim(file);close(unt);return;endif
! line 3: index of the matrix column that varies
      read(unt,*,iostat=err) xtra%cs2
      if(err/=0) then;mess=trim(procname)//':problem reading varying parameter (line 3) of '//trim(file);close(unt);return;endif

! matrix size n, counted from the first matrix row
      call getNumItemsInFile(unt=unt,preRwnd=.true.,nskip=3,nitems=n,postPos=3,&
         jchar=foo,err=err,message=mess)
      if(err/=0) then;mess=trim(procname)//':'//trim(mess);close(unt);return;endif

! matrix, n x n
      if(associated(xtra%rpm1)) nullify(xtra%rpm1)
      allocate(xtra%rpm1(n,n))
      do i=1,n
         read(unt,*,iostat=err) xtra%rpm1(i,:)
         if(err/=0) then;mess=trim(procname)//':problem reading matrix row '//trim(file);close(unt);return;endif
      enddo
      close(unt)

      call CheckControlMatrix(nint(xtra%rpm1),feas,err,mess)
      if(err/=0 .or. (.not.feas)) then
         mess=trim(procname)//':invalid control matrix'
         err=1;return
      endif
   end subroutine Vegetation_SH_XtraRead


   subroutine Vegetation_SH_Apply(growthFormula,parAffect,nYear,IN,theta,ControlMatrix,OUT,Dpar,feas,err,mess)

!^**********************************************************************
!^* Purpose:
!^**********************************************************************
!^* Programmer: Martin Morlot, INRAE Lyon
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
!^*     3. ControlMatrix, control Matrix in Laurent's framework
!^* OUT
!^*     1. OUT, col 1 = Q, 2 = g, 3 = cos(g), 4 = sin(g), 5 = g0,
!^*             (col 6 = S(h) if available)
!^*     2. Dpar, derived parameters: a1, then b offset
!^*     3. feas, feasible?
!^*     4. err, error code; <0:Warning, ==0:OK, >0: Error
!^*     5. mess, error message
!^**********************************************************************

      use types_dmsl_kit,only:data_ricz_type
      use utilities_dmsl_kit,only:cleanPointers
      use Vegetation_model,only:GetNpar_growth,Apply_growth
      use BaRatin_model , only:ApplyContinuity

      character(*), intent(in)::growthFormula,parAffect
      integer(mik), intent(in)::nYear
      real(mrk),intent(in)::IN(:,:),theta(:)
      integer(mik), intent(in)::ControlMatrix(:,:)
      real(mrk),intent(out)::OUT(:,:),Dpar(:)
      logical, intent(out)::feas(:)
      integer(mik), intent(out)::err
      character(*),intent(out)::mess
!locals
      character(250),parameter::procname='Vegetation_SH_Apply'
      integer(mik)::nH,ncontrol,i,p
      real(mrk), allocatable::a(:),b(:),c(:),k(:)
      logical::ok
      real(mrk)::time(size(IN,dim=1)),h(size(IN,dim=1)),Traw(size(IN,dim=1)),&
         Tsmooth(size(IN,dim=1)),g(size(IN,dim=1)),g0(size(IN,dim=1)),&
         cosg(size(IN,dim=1)),sing(size(IN,dim=1)),Q(size(IN,dim=1))
      real(mrk)::Width,S0,n1,b1,c1,alpha,lambda0,lambda1,Tm(nYear)
      integer(mik)::npar,nGrowthPar
      real(mrk),allocatable::growthPar(:)
      real(mrk),allocatable::gmax(:)
      real(mrk),allocatable::S(:)
      integer(mik)::nt,nVeg


      err=0;mess='';Q=undefRN;feas=.true.
! Make sense of input matrix
      nt=size(IN,dim=1)
      time=IN(:,1)
      h=IN(:,2)
      if(growthFormula/=Growth_Yin1_temporal) then
         Traw=IN(:,3)
         Tsmooth=IN(:,4)
      endif

! Get number of H values
      nH=size(H)
! Get number of controls (=number of ranges)
      ncontrol=size(ControlMatrix,1) ! it is assumed that the size of ControlMatrix has been checked outside of this sub

! Make sense of parameter vector
      call Vegetation_SH_GetParNumber(growthFormula,nYear,ControlMatrix,npar,err,mess)
      if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
      if(size(theta)/=npar) then;err=1;mess=trim(procname)//':wrong size[theta]';return;endif
! base parameters
      Width=theta(1)
      S0=theta(2)
      n1=theta(3)
      b1=theta(4)
      c1=theta(5)
      lambda0=theta(6)
      !lambda1=theta(7)
      !alpha=theta(8)
      if(S0>=0._mrk) then
         Dpar(1)=Width*sqrt(S0)/n1
      else
         Dpar(1)=undefRN
      endif

      call GetNpar_growth(growthFormula,nGrowthPar,err,mess)
      if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
      if(allocated(growthPar)) deallocate(growthPar);allocate(growthPar(nGrowthPar))
      growthPar=theta( (npar_Q0+npar_SH+1):(npar_Q0+npar_SH+nGrowthPar) )
! gmax
      if(allocated(gmax)) deallocate(gmax)
      if(growthFormula==Growth_Yin1_temporal) then
         allocate(gmax(1))
         nVeg=1
      else
         allocate(gmax(nYear))
         nVeg=nYear+1
      endif
      gmax=theta( (npar_Q0+nGrowthPar+npar_SH+1):(npar_Q0+nGrowthPar+npar_SH+size(gmax)) )
! Check parameter feasability
      if(Width<=0._mrk .or. S0<=0._mrk .or. n1<=0._mrk .or. c1<=0._mrk .or.&
         any(gmax<0._mrk)) then
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

      ! Allocate and assign values to k,a,c
      allocate(a(ncontrol),c(ncontrol),k(ncontrol),b(ncontrol))
      a(1)=Dpar(1)
      k(1)=b1
      c(1)=c1
      do i=1,ncontrol-1
         p=npar_Q0+nGrowthPar+size(gmax)+npar_SH+3*(i-1)
         k(i+1)=theta(p+1);a(i+1)=theta(p+2);c(i+1)=theta(p+3)
      enddo

! Apply Continuity Constraint & check feasability
      call ApplyContinuity(a,c,k,ControlMatrix,b,ok,err,mess)
      if(.not.ok) then;feas=.false.;return;endif
      Dpar(nVeg+1:nVeg+ncontrol)=b

      lambda1=k(2)
      alpha=1._mrk
      call VerifyFeasibility(lambda0,lambda1,alpha,k,ok)
      if(.not.ok) then;feas=.false.;return;endif

! S(h)
      allocate(S(nt))
      call SH_compute(h,lambda0,lambda1,alpha,S,ok)
      if(.not.ok) then;feas=.false.;return;endif
      OUT(:,6)=S

! Everything checked OK, can finally apply RC equation!
      call Vegeation_SH_computeQ(H,g,S,n1,b,a,c,k,ControlMatrix,parAffect,Q,ok)
      if(.not.ok) then;feas=.false.;return;endif
      OUT(:,1)=Q

   end subroutine Vegetation_SH_Apply


   subroutine Vegeation_SH_computeQ(H,g,S,n1,b,a,c,k,ControlMatrix,parAffect,Q,feas)

!^**********************************************************************
!^* Purpose: Apply RC equation. All checks on parameters are assumed
!^*          already done, and continuity equation has already been applied
!^**********************************************************************
!^* Programmer: Martin Morlot, INRAE Lyon
!^**********************************************************************
!^* Last modified: 10/01/2018
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*     1. H, water stages
!^*     2. b, offset parameters
!^*     3. a, coefficient parameters
!^*     4. c, exponent parameters
!^*     5. k, activation stages
!^*     6. ControlMatrix, control Matrix in Laurent's framework
!^* OUT
!^*     1. Q, RC-computed runoff
!^*     2. feas, feasible?
!^**********************************************************************
      use BaRatin_model, only:GetHRange


      real(mrk), intent(in)::H(:),g(:),S(:),b(:),a(:),c(:),k(:),n1
      integer(mik), intent(in)::ControlMatrix(:,:)
      character(*), intent(in)::parAffect
      real(mrk), intent(out)::Q(:)
      logical, intent(out)::feas
!locals
      integer(mik)::t,i,range
      real(mrk)::aEff,bEFF,fveg,top,bottom
      real(mrk),parameter::gravity=9.81_mrk


      Q=undefRN;feas=.true.

      do t=1,size(H)
         ! Retrieve range corresponding to stage H
         range=GetHRange(H(t),k)
         ! Apply RC equation
         Q(t)=0._mrk
         range_loop: do i=1,range ! if range == 0 (h smaller than first activation), loop is inactive and Q=0
            if(ControlMatrix(range,i)==1) then
               if((H(t)-b(i))<0._mrk) then
                  Q(t)=0._mrk;exit range_loop
               endif
               aEff=a(i)
               bEff=b(i)
               if(i==1) then
                  ! vegetation-affected control: roughness scaled by (1+g*S)
                  if(parAffect=="a") then
                     top=g(t)*S(t)*((H(t)-b(i))**(((1._mrk/3._mrk))))
                     bottom=8*gravity*(n1**2)
                     fveg=1._mrk+(top)/(bottom)
                     aEff=a(i)*fveg**(-1._mrk/2._mrk)
                  else if(parAffect=="b") then
                     fveg=g(t)*S(t)
                     bEff=b(i)+fveg
                  else
                     feas=.false.;exit range_loop
                  endif
               endif
               Q(t)=Q(t)+aEff*(H(t)-bEff)**c(i)
            else
               if((H(t)-b(i))<0._mrk) then
                  write(*,*) 'WARNING: You should NEVER EVER arrive there!!!'
                  write(*,*) 'If you see this message, please contact developpers !!!'
                  feas=.false.;exit range_loop
               endif
            endif
         enddo range_loop
      enddo

   end subroutine


   subroutine SH_compute(H,lambda0,lambda1,alpha,S,feas)

!^**********************************************************************
!^* Purpose: Compute the stage-dependent weighting S(h):
!^*             S = 1               if h < lambda1
!^*             S = 0               if h > lambda2
!^*             S = 1 - x^alpha / (x^alpha + (1-x)^alpha) otherwise,
!^*                 x = (h-lambda1)/(lambda2-lambda1)
!^*          alpha = 1 gives a linear transition.
!^**********************************************************************
!^* Programmer: Martin Morlot, INRAE Lyon
!^**********************************************************************
!^* Last modified: 20/10/2025
!^**********************************************************************
!^* IN
!^*     1. H, water stages
!^*     2. lambdas, (lambda1, lambda2) with lambda1<lambda2
!^*     3. alpha, shape of the transition (>=1)
!^* OUT
!^*     1. S, S(h) values
!^*     2. feas, feasible?
!^**********************************************************************

      real(mrk), intent(in)::H(:),lambda0,lambda1,alpha
      real(mrk), intent(out)::S(:)
      logical, intent(out)::feas
!locals
      integer(mik)::t
      real(mrk)::x,xa,ya

      S=undefRN;feas=.true.

      do t=1,size(H)

         ! Apply RC equation
         S(t)=0._mrk

         if(H(t)<lambda0) then
            S(t) = 1._mrk
         else if(H(t)>lambda1) then
            S(t) = 0._mrk
         else
            x=(H(t)-lambda0)/(lambda1-lambda0)
            xa=x**alpha
            ya=(1._mrk-x)**alpha
            S(t)=1._mrk-xa/(xa+ya)
         endif
      enddo
   end subroutine

   subroutine VerifyFeasibility(lambda0,lambda1,alpha,k,feas)

!^**********************************************************************
!^* Purpose: Check feasibility of the parameters for S
!^**********************************************************************
!^* Programmer: Martin Morlot, INRAE Lyon
!^**********************************************************************
!^* Last modified: 20/10/2025
!^**********************************************************************
!^* IN
!^*     1. lambdas, (lambda0, lambda1) with lambda1<lambda2
!^*     2. k, shape of the transition (>=1)
!^* OUT
!^*     1. feas, feasible?
!^**********************************************************************

      real(mrk), intent(in)::k(:),lambda0,lambda1,alpha
      logical, intent(out)::feas

      feas=.true.

      if(lambda0>=lambda1 .or. alpha<1._mrk) then
         feas=.false.;return
      endif

      if(lambda0<k(1)) then
         feas=.false.;return
      endif

      if(lambda1>k(2)) then
         feas=.false.;return
      endif

   end subroutine

end module Vegetation_SH_model
