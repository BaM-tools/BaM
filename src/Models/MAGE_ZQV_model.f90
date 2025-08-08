module MAGE_ZQV_model

!~**********************************************************************
!~* Purpose: Interface to MAGE model
!~**********************************************************************
!~* Programmer: Ben Renard and Felipe Mendez and Theophile Terraz, INRAE Aix & INRAE Lyon
!~**********************************************************************
!~* Last modified: 29/07/2025
!~**********************************************************************
!~* Comments:
!~**********************************************************************
!~* References:
!~**********************************************************************
!~* 2Do List:
!~**********************************************************************
!~* Quick description of public procedures:
!~*     1. MAGE_GetParNumber, number of inferred parameters
!~*     2. MAGE_Run, run MAGE
!~*     3. MAGE_XtraRead, read Xtra model information
!~*     4. MAGE_Setup, complete model setup after reading Xtra information
!~**********************************************************************

use kinds_dmsl_kit ! numeric kind definitions from DMSL
use, intrinsic :: iso_fortran_env, only: real32, real64

implicit none
Private
public :: MAGE_ZQV_GetParNumber, MAGE_ZQV_Run, MAGE_ZQV_XtraRead, MAGE_ZQV_SetUp, mage_res

! Module variables
character(len_uLongStr)::RUGfile='' ! .RUG file specifying rugosities
real(mrk), allocatable::RUGx_start(:),RUGx_end(:) ! x's defining constant-rugosity patches
integer(mik), allocatable::RUGb(:) ! "Bief" (reach?) index
real(mrk), allocatable::ZxKmin(:,:) ! covariates Z's of the regression Kmin(x)=a1Z1(x)+...+apZp(x). Dimension size(RUGb)*p.
real(mrk), allocatable::ZxKmoy(:,:) ! covariates Z's of the regression Kmoy(x)=a1Z1(x)+...+apZp(x). Dimension size(RUGb)*p.
logical::applyExpKmin,applyExpKmoy ! apply exponential transformation to computed K's to ensure positivity?

type mage_res
  integer :: ibmax, ismax ! nombre de biefs et nombre de sections
  integer, allocatable :: is1(:), is2(:) ! pointeurs sur première et dernière sections du bief
  real(kind=real64), dimension(:), allocatable :: tQ, tZ, tV ! pas de temps
  real(kind=real32), dimension(:), allocatable :: pk ! pk
  real(kind=real32), dimension(:,:), allocatable :: valuesQ, valuesZ, valuesV ! valeurs
  contains
    procedure :: read_bin
    procedure :: get
end type mage_res

Contains

subroutine read_bin(self, filename, err, mess)

  implicit none
  ! prototype
  class(mage_res), intent(inout) :: self
  character(len=*), intent(in) :: filename
  ! variables locales
  integer :: lu ! unité fichier
  character :: line*80
  integer :: mage_version, ntq, ntz, ntv, np, i
  integer :: nmaxq, nmaxz, nmaxv, qptr, zptr, vptr
  integer, parameter :: dnmax = 1000
  real(kind=real64) :: t
  real(kind=real64), dimension(:), allocatable :: tmp1 ! tableau temporaire
  real(kind=real32), dimension(:), allocatable :: val ! tableau temporaire
  real(kind=real32), dimension(:,:), allocatable :: tmp2 ! tableau temporaire
  character :: a*1
  integer(mik), intent(out)::err ! i/o status
  character(*),intent(out)::mess
  character(8),parameter::procname='read_bin'

  ! ouverture du fichier bin
  open(newunit=lu,file=trim(filename),status='old',form='unformatted',iostat=err)
  if (err /= 0) then
    mess=trim(procname)//': Problem opening file '//trim(filename)
    return
  endif

  ! lecture première ligne
  read(lu) self%ibmax, self%ismax, mage_version
  if (mage_version <= 80) then
    ! on a un fichier BIN avec l'ancien entête
    err=1
    mess=trim(procname)//': FATAL: MAGEv7 is not supported'
    return
  endif

  ! allocation des tableaux
  nmaxq = 1000
  nmaxz = 1000
  nmaxv = 1000
  qptr = 0
  zptr = 0
  vptr = 0
  if(allocated(self%tQ)) deallocate(self%tQ)
  if(allocated(self%tZ)) deallocate(self%tZ)
  if(allocated(self%tV)) deallocate(self%tV)
  if(allocated(self%valuesQ)) deallocate(self%valuesQ)
  if(allocated(self%valuesZ)) deallocate(self%valuesZ)
  if(allocated(self%valuesV)) deallocate(self%valuesV)
  if(allocated(self%is1)) deallocate(self%is1)
  if(allocated(self%is2)) deallocate(self%is2)
  if(allocated(self%pk)) deallocate(self%pk)
  if(allocated(val)) deallocate(val)
  allocate (self%tQ(nmaxq), self%tZ(nmaxz), self%tV(nmaxz))
  allocate (self%valuesQ(self%ismax, nmaxq), self%valuesZ(self%ismax, nmaxz), self%valuesV(self%ismax, nmaxv))
  allocate (val(self%ismax))
  allocate (self%is1(self%ibmax), self%is2(self%ibmax))
  ! suite de l'entête
  read(lu) (self%is1(i),self%is2(i),i=1,self%ibmax)
  allocate(self%pk(self%ismax))
  read(lu) (self%pk(i),i=1,self%ismax)
  read(lu)

  do
    read(lu,iostat=err) np,t,a,(val(i),i=1,np)
    if (err < 0) exit ! fin du fichier
    if (a == "Q") then ! débit
      qptr = qptr + 1
      if (qptr > nmaxq) then  ! re-allocation mémoire
        allocate(tmp2(self%ismax, nmaxq+dnmax)) ! on ajoute dnmax à la taille allouée
        tmp2(1:self%ismax, 1:nmaxq) = self%valuesQ(1:self%ismax, 1:nmaxq) !on copie les anciennes valeurs de valuesQ
        call move_alloc(tmp2, self%valuesQ) ! valuesQ pointe sur l'espace mémoire de tmp2, tmp2 est désalloué
        allocate(tmp1(nmaxq+dnmax))
        tmp1(1:nmaxq) = self%tQ(1:nmaxq)
        call move_alloc(tmp1, self%tQ)
        nmaxq = nmaxq+dnmax
      endif
      self%tQ(qptr) = t
      self%valuesQ(:, qptr) = val(:)
    elseif (a == "Z") then ! cote
      zptr = zptr + 1
      if (zptr > nmaxz) then  ! re-allocation mémoire
        allocate(tmp2(self%ismax, nmaxz+dnmax)) ! on ajoute dnmax à la taille allouée
        tmp2(1:self%ismax, 1:nmaxz) = self%valuesZ(1:self%ismax, 1:nmaxz) !on copie les anciennes valeurs de valuesZ
        call move_alloc(tmp2, self%valuesZ) ! valuesZ pointe sur l'espace mémoire de tmp2, tmp2 est désalloué
        allocate(tmp1(nmaxz+dnmax))
        tmp1(1:nmaxz) = self%tZ(1:nmaxz)
        call move_alloc(tmp1, self%tZ)
        nmaxz = nmaxz+dnmax
      endif
      self%tZ(zptr) = t
      self%valuesZ(:, zptr) = val(:)
    elseif (a == "V") then ! vitesse
      vptr = vptr + 1
      if (vptr > nmaxv) then  ! re-allocation mémoire
        allocate(tmp2(self%ismax, nmaxv+dnmax)) ! on ajoute dnmax à la taille allouée
        tmp2(1:self%ismax, 1:nmaxv) = self%valuesV(1:self%ismax, 1:nmaxv) !on copie les anciennes valeurs de valuesV
        call move_alloc(tmp2, self%valuesV) ! valuesV pointe sur l'espace mémoire de tmp2, tmp2 est désalloué
        allocate(tmp1(nmaxv+dnmax))
        tmp1(1:nmaxv) = self%tV(1:nmaxv)
        call move_alloc(tmp1, self%tV)
        nmaxv = nmaxv+dnmax
      endif
      self%tV(vptr) = t
      self%valuesV(:, vptr) = val(:)
    else
      cycle
    endif
  enddo
  err = 0
  close(lu)

  ! réallocation à la bonne taille
  ! Q
  allocate(tmp2(self%ismax, qptr))
  tmp2(1:self%ismax, 1:qptr) = self%valuesQ(1:self%ismax, 1:qptr)
  call move_alloc(tmp2, self%valuesQ)
  allocate(tmp1(qptr))
  tmp1(1:qptr) = self%tQ(1:qptr)
  call move_alloc(tmp1, self%tQ)
  ! Z
  allocate(tmp2(self%ismax, zptr))
  tmp2(1:self%ismax, 1:zptr) = self%valuesZ(1:self%ismax, 1:zptr)
  call move_alloc(tmp2, self%valuesZ)
  allocate(tmp1(zptr))
  tmp1(1:zptr) = self%tZ(1:zptr)
  call move_alloc(tmp1, self%tZ)
  ! V
  allocate(tmp2(self%ismax, vptr))
  tmp2(1:self%ismax, 1:vptr) = self%valuesV(1:self%ismax, 1:vptr)
  call move_alloc(tmp2, self%valuesV)
  allocate(tmp1(vptr))
  tmp1(1:vptr) = self%tV(1:vptr)
  call move_alloc(tmp1, self%tV)

end subroutine read_bin

subroutine get(self,reach,pk,timesteps,variable,res,err,mess)

  implicit none
  ! prototype
  class(mage_res), target, intent(in) :: self
  integer(mik), dimension(:), intent(in) :: reach
  real(kind=real64), dimension(:), intent(in) :: pk, timesteps
  character, intent(in) :: variable*1
  integer(mik), intent(out)::err
  character(*),intent(out)::mess
  real(kind=real64), dimension(:), intent(out) :: res
  ! variables locales
  character(250),parameter::procname='MAGE_ZQV_get'
  real(kind=real64), dimension(:), pointer :: t ! pas de temps
  real(kind=real32), dimension(:,:), allocatable :: tmp2 ! tableau temporaire
  integer :: i, j, ptr, pkptr, tvalptr, pkvalptr
  real(kind=real32), dimension(:,:), pointer :: values
  real(kind=real32) :: v1, v2, ratio_pk
  real(kind=real64) :: ratio_t
  real(kind=real32), parameter :: eps = 0.0001 ! max percision
  logical :: trouvePk, trouvet

  if (size(reach) .ne. size(timesteps)) stop
  if (size(pk) .ne. size(timesteps)) stop
  if (size(res) .ne. size(timesteps)) stop

  if (variable == "Q" .or. variable == "q") then ! débit
    values => self%valuesQ
    t => self%tQ
  elseif (variable == "Z" .or. variable == "z") then ! cote
    values => self%valuesZ
    t => self%tZ
  elseif (variable == "V" .or. variable == "v") then ! vitesse
    values => self%valuesV
    t => self%tV
  endif

  do ptr = 1, size(res)
    trouvePk = .false.
    trouvet = .false.
    do i=1, size(t)-1
       if (t(i) < timesteps(ptr) .and. t(i+1) >= timesteps(ptr)) then
          ! trouvé !
          trouvet = .true.
          do j=self%is1(reach(ptr)), self%is2(reach(ptr))-1
             if (self%pk(j) < pk(ptr) .and. self%pk(j+1) >= pk(ptr)) then
                ! trouvé !
                trouvePk = .true.
                ratio_pk = (pk(ptr) - self%pk(j))/(self%pk(j+1) - self%pk(j))
                ratio_t = (timesteps(ptr) - t(i))/(t(i+1) - t(i))
                if (ratio_pk < eps) then
                    v1 = values(j, i)
                    v2 = values(j, i+1)
                else if (ratio_pk > 1.0 - eps) then
                    v1 = values(j+1, i)
                    v2 = values(j+1, i+1)
                else
                    v1 = values(j, i)*(1.0 - ratio_pk) + values(j+1, i)*(ratio_pk)
                    v2 = values(j, i+1)*(1.0 - ratio_pk) + values(j+1, i+1)*(ratio_pk)
                endif
                if (ratio_t < eps) then
                    res(ptr) = v1
                else if (ratio_t > 1.0 - eps) then
                    res(ptr) = v2
                else
                    res(ptr) = v1*(1.0 - ratio_t) + v2*(ratio_t)
                endif
                exit
             endif
          enddo
          if(.not.trouvePk)then
             mess=trim(procname)//" value not found in Mage result file (bad Pk)"
             err=1;return
          endif
          exit
       endif
    enddo
    if(.not.trouvet)then
       mess=trim(procname)//" value not found in Mage result file (bad timestamp)"
       err=1;return
    endif
 enddo

end subroutine get

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine MAGE_ZQV_GetParNumber(npar,err,mess)
!^**********************************************************************
!^* Purpose: number of inferred parameters
!^*          = number of "troncons" * 2 (main chanel and floodway)
!^**********************************************************************
!^* Programmer: Ben Renard and Felipe Mendez, INRAE Aix & INRAE Lyon
!^**********************************************************************
!^* Last modified:29/08/2022
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*     nothing
!^* OUT
!^*     1. npar, par. number
!^*     2.err, error code; <0:Warning, ==0:OK, >0: Error
!^*     3.mess, error message
!^**********************************************************************
integer(mik), intent(out)::npar,err
character(*),intent(out)::mess
! locals
character(250),parameter::procname='MAGE_ZQV_GetParNumber'

err=0;mess=''
npar=size(ZxKmin,dim=2)+size(ZxKmoy,dim=2)

end subroutine MAGE_ZQV_GetParNumber

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine MAGE_ZQV_Run(exeFile,version,projectDir,REPfile,&
                        X,theta,Y,feas,err,mess)
!^**********************************************************************
!^* Purpose: Run MAGE
!^**********************************************************************
!^* Programmer: Ben Renard and Felipe Mendez, INRAE Aix & INRAE Lyon
!^**********************************************************************
!^* Last modified: 29/08/2022
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*     1. exeFile, MAGE executable
!^*     2. version, MAGE version (e.g. 7 or 8)
!^*     3. projectDir, MAGE project directories (1-D array)
!^*     4. REPfile, MAGE .REP file
!^*     5. X, inputs: event number, reach number, x, t
!^*     6. theta, parameter vector
!^* OUT
!^*     1. Y, outputs: (Z,Q,V) computed at each location/time requested through X
!^*     3. feas, feasible?
!^*     4. err, error code; <0:Warning, ==0:OK, >0: Error
!^*     5. mess, error message
!^**********************************************************************
use utilities_dmsl_kit,only:getSpareUnit,replacedString,number_string,getNumItemsInFile,countSubstringInString
character(*), intent(in)::exeFile,version,projectDir(:),REPfile
real(mrk), intent(in)::X(:,:),theta(:)
real(mrk), intent(out)::Y(:,:)
logical, intent(out)::feas
integer(mik), intent(out)::err
character(*),intent(out)::mess
!locals
character(250),parameter::procname='MAGE_ZQV_Run'
integer(mik)::i,j,n,unt,k,ok,nitems,nevents,p,m
character(250)::forma
character(len_uLongStr)::cmdString,line
real(mrk)::Kmin(size(RUGb)),Kmoy(size(RUGb))
character(:), allocatable :: project
real(mrk), allocatable :: tim(:),loc(:)
integer(mik), allocatable :: reach(:)
type(mage_res) :: res
logical::mask(size(X,dim=1))

! Init
err=0;mess='';feas=.true.;Y=undefRN
n=size(RUGb)
nevents=size(projectDir)
forma='(A1,I3,6x,2f10.3,2f10.2)'
project=replacedString(trim(REPfile),'.REP','',.true.)

! Size check
k=size(theta)
if(k /= (size(ZxKmin,dim=2)+size(ZxKmoy,dim=2))) then
    mess=trim(procname)//': length of theta should be equal to the total number of columns in Zfiles'
    err=1;return
endif

! Compute Kmin and Kmoy from regression
Kmin=matmul(ZxKmin, theta( 1:size(ZxKmin,dim=2) ) )
Kmoy=matmul(ZxKmoy, theta( (size(ZxKmin,dim=2)+1):k ) )
if(applyExpKmin) Kmin=exp(Kmin)
if(applyExpKmoy) Kmoy=exp(Kmoy)

m=0
do j=1,nevents
    ! Determine X rows related to this event
    mask=nint(X(:,1),mik)==j
    p=count(mask)
    if(p==0) then
        cycle
    else
        if(allocated(tim)) deallocate(tim)
        if(allocated(loc)) deallocate(loc)
        allocate(tim(p),loc(p))
        reach=pack(X(:,2),mask)
        loc=pack(X(:,3),mask)
        tim=pack(X(:,4),mask)
    endif

    ! Write Kmin and Kmoy in .RUG file
    call getSpareUnit(unt,err,mess)
    if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
    open(unit=unt,file=trim(projectDir(j))//trim(RUGfile),status='replace',iostat=err)
    if(err/=0) then;mess=trim(procname)//':problem creating RUG file';return;endif
    select case(version)
    case('7','v7')
        write(unt,'(A)') '**************** STRICKLERS **********************'
        write(unt,'(A)') '*Bief        x deb     x fin     K min     K moy  '
        write(unt,'(A)') '*---======----------==========----------=========='
    case('8','v8')
        write(unt,'(A)') "* This file is generated by PAMHYR, please don't modify"
    case default
        if(err/=0) then;mess=trim(procname)//':FATAL: MAGE version is not supported';return;endif
    end select

    do i=1,n
        write(unt,forma) 'K',RUGb(i),RUGx_start(i),RUGx_end(i),Kmin(i),Kmoy(i)
    enddo
    close(unt)

    ! Call mage
    cmdString='cd '//trim(projectDir(j))//'&&'//trim(exeFile)//' '//trim(REPfile)
    call execute_command_line (trim(cmdString),wait=.true.,exitStat=err,cmdMsg=mess)
    if(err/=0) then
        err = 0
        feas=.false.
        mess=trim(procname)//':'//trim(mess);return
    endif

    ! Verify MAGE ran correctly
    call getSpareUnit(unt,err,mess)
    if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
    open(unit=unt,file=trim(projectDir(j))//trim(project)//trim('.TRA'),status='old',iostat=err)
    if(err/=0) then;mess=trim(procname)//':problem opening TRA file';return;endif
    call getNumItemsInFile(unt=unt,preRwnd=.true.,nskip=0,nitems=nitems,postPos=0,jchar=line,err=err,message=mess)
    if(err/=0) then;mess=trim(procname)//':problem reading TRA file';return;endif
    close(unt)
    ok = countSubstringInString(line,'FIN NORMALE DE MAGE')
    if(ok == 0) then
        feas=.false.;return
    endif

    ! read .BIN
    call res%read_bin(trim(projectDir(j))//project//'.BIN', err, mess) ! on remplit l'objet

    call res%get(reach,loc,tim,"Z",Y((m+1):(m+p),1),err,mess)
    call res%get(reach,loc,tim,"Q",Y((m+1):(m+p),2),err,mess)
    call res%get(reach,loc,tim,"V",Y((m+1):(m+p),3),err,mess)
    m=m+p
enddo

end subroutine MAGE_ZQV_Run

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine MAGE_ZQV_XtraRead(file,xtra,err,mess)
!^**********************************************************************
!^* Purpose: Read Xtra information
!^**********************************************************************
!^* Programmer: Ben Renard and Felipe Mendez, INRAE Aix & INRAE Lyon
!^**********************************************************************
!^* Last modified: 29/08/2022
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
use utilities_dmsl_kit,only:getSpareUnit,countSubstringInString
character(*), intent(in)::file
type(data_ricz_type),intent(out):: xtra
integer(mik), intent(out)::err
character(*),intent(out)::mess
! locals
character(250),parameter::procname='MAGE_XtraRead',sep=','
integer(mik)::unt,s,n
character(len_uLongStr)::foo

err=0;mess=''

call getSpareUnit(unt,err,mess)
if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
open(unit=unt,file=trim(file), status='old', iostat=err)
if(err>0) then;mess=trim(procname)//':problem opening file '//trim(file);return;endif
read(unt,*,iostat=err) xtra%cs1 ! 1. exe file
if(err/=0) then;mess=trim(procname)//':problem reading file '//trim(file);return;endif
read(unt,*,iostat=err) xtra%cs2 ! 2. MAGE version
if(err/=0) then;mess=trim(procname)//':problem reading file '//trim(file);return;endif
! Project directories (one per event): First read to get number of elements
read(unt,'(A)',iostat=err) foo
if(err/=0) then;mess=trim(procname)//':problem reading file '//trim(file);return;endif
n=countSubstringInString(foo,trim(sep))+1
if(associated(xtra%cp1)) nullify(xtra%cp1)
allocate(xtra%cp1(n))
read(foo,*,iostat=err) xtra%cp1 ! 3. Project Directories
if(err/=0) then;mess=trim(procname)//':problem reading file '//trim(file);return;endif
read(unt,*,iostat=err) xtra%cs4 ! 4. REP file
if(err/=0) then;mess=trim(procname)//':problem reading file '//trim(file);return;endif
! K-Z regressions
read(unt,*,iostat=err) xtra%cs5    ! 5. Z file (covariates for Kmin-Z regression)
if(err/=0) then;mess=trim(procname)//':problem reading file '//trim(file);return;endif
read(unt,*,iostat=err) xtra%is1    ! 6. Number of columns in Z file
if(err/=0) then;mess=trim(procname)//':problem reading file '//trim(file);return;endif
read(unt,*,iostat=err) xtra%ls1    ! 7. apply exponential transformation to computed  Kmin's?
if(err/=0) then;mess=trim(procname)//':problem reading file '//trim(file);return;endif
read(unt,*,iostat=err) xtra%cs6    ! 8. Z file (covariates for Kmoy-Z regression)
if(err/=0) then;mess=trim(procname)//':problem reading file '//trim(file);return;endif
read(unt,*,iostat=err) xtra%is2    ! 9. Number of columns in Z file
if(err/=0) then;mess=trim(procname)//':problem reading file '//trim(file);return;endif
read(unt,*,iostat=err) xtra%ls2    ! 10. apply exponential transformation to computed  Kmoy's?
if(err/=0) then;mess=trim(procname)//':problem reading file '//trim(file);return;endif
close(unt)

end subroutine MAGE_ZQV_XtraRead

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine MAGE_ZQV_setUp(version,projectDir,REPfile,Zfile,ZnCol,doExp,err,mess)
!^**********************************************************************
!^* Purpose: Finalize MAGE setup after having read Xtra information. In
!^*          particular this sub defines module variables (RUG stuff and outFiles)
!^**********************************************************************
!^* Programmer: Ben Renard and Felipe Mendez, INRAE Aix & INRAE Lyon
!^**********************************************************************
!^* Last modified: 24/02/2023
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*     1. version, MAGE version
!^*     2. projectDir, MAGE project directory
!^*     3. REPfile, MAGE .REP file
!^*     4. Zfile, containing covariates for K-Z regression
!^*     5. doExp, apply exponential transformation to computed  K's?
!^*     6. mage_extraire_args, arguments to mage_extraire, defines output variables
!^* OUT
!^*     1.err, error code; <0:Warning, ==0:OK, >0: Error
!^*     2.mess, error message
!^**********************************************************************
use utilities_dmsl_kit,only:getNumItemsInFile,getSpareUnit,replacedString
character(*), intent(in)::version,projectDir(:),REPfile,Zfile(2)
integer(mik), intent(in)::ZnCol(2)
logical, intent(in)::doExp(2)
integer(mik), intent(out)::err
character(*),intent(out)::mess
!locals
character(250),parameter::procname='MAGE_ZQV_setUp'
integer(mik)::unt,i,n,k,j,nskip
character(250)::foo,prefix,var,ib,pk,mode,dt0,forma,project
character(len_uLongStr)::line,fname
character(len_uLongStr),allocatable::fnames(:)
real(mrk)::c1,c2
real(mrk),allocatable::Zx(:,:)

! Init
err=0;mess=''

! Open REP file
call getSpareUnit(unt,err,mess)
if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
open(unit=unt,file=trim(projectDir(1))//trim(REPfile), status='old', iostat=err)
if(err>0) then;mess=trim(procname)//':problem opening file '//trim(projectDir(1))//trim(REPfile);return;endif
call getNumItemsInFile(unt=unt,preRwnd=.true.,nskip=0,nitems=n,postPos=0,jchar=foo,err=err,message=mess)
if(err>0) then;mess=trim(procname)//trim(mess);return;endif

! Interpret REP file
k=0
if(allocated(fnames)) deallocate(fnames);allocate(fnames(n))
do i=1,n
    read(unt,'(A)',iostat=err) line
    if(err>0) then;mess=trim(procname)//':problem reading file '//trim(REPfile);return;endif
    if(line(1:3)=='RUG') RUGfile=trim(line(5:))
enddo
close(unt)

! Read RUG file
select case(version)
case('7','v7')
    nskip=3
case('8','v8')
    nskip=1
case default
    if(err/=0) then;mess=trim(procname)//':FATAL: MAGE version is not supported';return;endif
end select
if(RUGfile=='') then;err=1;mess=trim(procname)//'FATAL: no .RUG file found';return;endif
call getSpareUnit(unt,err,mess)
if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
open(unit=unt,file=trim(projectDir(1))//trim(RUGfile), status='old', iostat=err)
if(err>0) then;mess=trim(procname)//':problem opening file '//trim(projectDir(1))//trim(RUGfile);return;endif
call getNumItemsInFile(unt=unt,preRwnd=.true.,nskip=nskip,nitems=n,postPos=0,jchar=foo,err=err,message=mess)
if(err>0) then;mess=trim(procname)//trim(mess);return;endif
if(allocated(RUGx_start)) deallocate(RUGx_start)
if(allocated(RUGx_end)) deallocate(RUGx_end)
if(allocated(RUGb)) deallocate(RUGb)
allocate(RUGx_start(n),RUGx_end(n),RUGb(n))
rewind(unt)
do i=1,nskip;read(unt,*);enddo
forma='(1x,I3,6x,4f10.0)'
do i=1,n
    read(unt,forma,iostat=err) RUGb(i),RUGx_start(i),RUGx_end(i),c1,c2
    if(err/=0) then;mess=trim(procname)//':problem reading file '//trim(projectDir(1))//trim(RUGfile);return;endif
enddo
close(unt)

! Read Z file
n=size(RUGb)
do j=1,2
    if(Zfile(j)=='') then
        if(allocated(Zx)) deallocate(Zx);allocate(Zx(n,n))
        Zx(:,:)=0._mrk
        forall(i=1:n)Zx(i,i)=1._mrk
    else
        call getSpareUnit(unt,err,mess)
        if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
        open(unit=unt,file=trim(Zfile(j)), status='old', iostat=err)
        if(err>0) then;mess=trim(procname)//':problem opening file '//trim(Zfile(j));return;endif
        call getNumItemsInFile(unt=unt,preRwnd=.true.,nskip=0,nitems=k,postPos=0,jchar=foo,err=err,message=mess)
        if(err>0) then;mess=trim(procname)//trim(mess);return;endif
        if(k /= (n+1)) then
            mess=trim(procname)//':Zfile should have n+1 rows (n rows of the RUG table + one header line). Zfile='//trim(Zfile(j))
            err=1;return
        endif
        if(allocated(Zx)) deallocate(Zx);allocate(Zx(n,ZnCol(j)))
        rewind(unt)
        read(unt,*) ! header
        do i=1,n
            read(unt,*,iostat=err) Zx(i,:)
            if(err/=0) then;mess=trim(procname)//':problem reading file '//trim(Zfile(j));return;endif
        enddo
        close(unt)
    endif
    if(j==1) then ! Applies to Kmin
        if(allocated(ZxKmin)) deallocate(ZxKmin);allocate(ZxKmin(size(Zx,dim=1),size(Zx,dim=2)))
        ZxKmin=Zx
    else ! Applies to Kmoy
        if(allocated(ZxKmoy)) deallocate(ZxKmoy);allocate(ZxKmoy(size(Zx,dim=1),size(Zx,dim=2)))
        ZxKmoy=Zx
    endif
enddo
if(allocated(Zx)) deallocate(Zx)

! exp. transformation
applyExpKmin=doExp(1)
applyExpKmoy=doExp(2)

end subroutine MAGE_ZQV_setUp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!==============!
! Private subs !
!==============!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module MAGE_ZQV_model
