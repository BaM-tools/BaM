program testSMASH

use kinds_dmsl_kit
use SMASH_model
use utilities_dmsl_kit,only:getSpareUnit
use DataRW_tools,only:ReadRasterASC, WriteRasterASC
use Geodesy_tools,only:rasterGridType,rasterToPoints,interpolateOnGrid

implicit none

character(len_vLongStr), parameter:: setupFile='../../tests/BaM_SMASH/Config_setup.txt'
logical, parameter:: reload=.true.
integer(mik):: err,unt,i
character(250):: mess,startDate,endDate
character(len_vLongStr)::loadScript,runScript,projectDir,precipDir,petDir,thetaDir
character(len_uLongStr)::cmdString
type(rasterGridType)::grid
real(mrk),allocatable::M(:,:)
real(mrk),allocatable::points(:,:)
real(mrk)::anchor(1,3),anchors(4,3)

! Get grid
call ReadRasterASC(file='/home/brenard/Desktop/2022SMASH/V5004030/PARAM/FLOW.asc',gridOnly=.true.,grid=grid,M=M,err=err,mess=mess)
if(err/=0) then;write(*,*) trim(mess);read(*,*);endif

! Paramter cp
anchors(1,:)=(/grid%xllcorner,grid%yllcorner,500._mrk/)
anchors(2,:)=(/grid%xllcorner+(grid%ncols+1)*grid%cellsize,grid%yllcorner,800._mrk/)
anchors(3,:)=(/grid%xllcorner,grid%yllcorner+(grid%nrows+1)*grid%cellsize,1000._mrk/)
anchors(4,:)=(/grid%xllcorner+(grid%ncols+1)*grid%cellsize,grid%yllcorner+(grid%nrows+1)*grid%cellsize,100._mrk/)
call interpolateOnGrid(pts=anchors,grid=grid,method='IDW',par=(/1._mrk/),M=M,err=err,mess=mess)
if(err/=0) then;write(*,*) trim(mess);read(*,*);endif
call WriteRasterASC(grid=grid,M=M,file='/home/brenard/Desktop/2022SMASH/V5004030/THETA/cp.asc',err=err,mess=mess)
if(err/=0) then;write(*,*) trim(mess);read(*,*);endif

! Paramter ctr
anchor(1,:)=(/grid%xllcorner,grid%yllcorner,100._mrk/)
call interpolateOnGrid(pts=anchor,grid=grid,method='IDW',par=(/1._mrk/),M=M,err=err,mess=mess)
if(err/=0) then;write(*,*) trim(mess);read(*,*);endif
call WriteRasterASC(grid=grid,M=M,file='/home/brenard/Desktop/2022SMASH/V5004030/THETA/ctr.asc',err=err,mess=mess)
if(err/=0) then;write(*,*) trim(mess);read(*,*);endif

! Paramter cr
anchor(1,:)=(/grid%xllcorner,grid%yllcorner,600._mrk/)
call interpolateOnGrid(pts=anchor,grid=grid,method='IDW',par=(/1._mrk/),M=M,err=err,mess=mess)
if(err/=0) then;write(*,*) trim(mess);read(*,*);endif
call WriteRasterASC(grid=grid,M=M,file='/home/brenard/Desktop/2022SMASH/V5004030/THETA/cr.asc',err=err,mess=mess)
if(err/=0) then;write(*,*) trim(mess);read(*,*);endif

! Paramter cr
anchor(1,:)=(/grid%xllcorner,grid%yllcorner,0._mrk/)
call interpolateOnGrid(pts=anchor,grid=grid,method='IDW',par=(/1._mrk/),M=M,err=err,mess=mess)
if(err/=0) then;write(*,*) trim(mess);read(*,*);endif
call WriteRasterASC(grid=grid,M=M,file='/home/brenard/Desktop/2022SMASH/V5004030/THETA/ml.asc',err=err,mess=mess)
if(err/=0) then;write(*,*) trim(mess);read(*,*);endif

! test rasterToPoints sub
!call rasterToPoints(grid=grid,M=M,pts=points,err=err,mess=mess)
!if(err/=0) then;write(*,*) trim(mess);read(*,*);endif
!do i=1,44 !size(points,dim=1)
!    write(*,*) points(i,:)
!enddo
!read(*,*)

call getSpareUnit(unt,err,mess)
if(err/=0) then;write(*,*) trim(mess);read(*,*);endif
open(unit=unt,file=trim(setupFile), status='old', iostat=err)
if(err>0) then;write(*,*) trim(mess);read(*,*);endif
read(unt,*,iostat=err) loadScript
if(err>0) then;write(*,*) trim(mess);read(*,*);endif
read(unt,*,iostat=err) runScript
if(err>0) then;write(*,*) trim(mess);read(*,*);endif
read(unt,*,iostat=err) projectDir
if(err>0) then;write(*,*) trim(mess);read(*,*);endif
read(unt,*,iostat=err) precipDir
if(err>0) then;write(*,*) trim(mess);read(*,*);endif
read(unt,*,iostat=err) petDir
if(err>0) then;write(*,*) trim(mess);read(*,*);endif
read(unt,*,iostat=err) thetaDir
if(err>0) then;write(*,*) trim(mess);read(*,*);endif
read(unt,*,iostat=err) startDate
if(err>0) then;write(*,*) trim(mess);read(*,*);endif
read(unt,*,iostat=err) endDate
if(err>0) then;write(*,*) trim(mess);read(*,*);endif

if(reload) then
    cmdString='python3'//' '//trim(loadScript)//' '//trim(projectDir)&
              &//' '//trim(precipDir)//' '//trim(petDir)&
              &//' '//trim(startDate)//' '//trim(endDate)
    call execute_command_line (trim(cmdString),exitStat=err,cmdMsg=mess)
endif

cmdString='python3'//' '//runScript//' '//projectDir//' '//thetaDir
call execute_command_line (trim(cmdString),exitStat=err,cmdMsg=mess)

write(*,*) 'Exit status:', err
write(*,*) 'Exit message:', mess

write(*,*) 'FINITO'
read(*,*)

end program testSMASH
