program testSMASH

use kinds_dmsl_kit
use SMASH_model
use utilities_dmsl_kit,only:getSpareUnit

implicit none

character(len_vLongStr), parameter:: setupFile='../../tests/BaM_SMASH/Config_setup.txt'
logical, parameter:: reload=.false.

integer(mik):: err,unt
character(250):: mess,startDate,endDate
character(len_vLongStr)::loadScript,runScript,projectDir,precipDir,petDir
character(len_uLongStr)::cmdString

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

cmdString='python3'//' '//runScript//' '//projectDir
call execute_command_line (trim(cmdString),exitStat=err,cmdMsg=mess)

write(*,*) 'Exit status:', err
write(*,*) 'Exit message:', mess

write(*,*) 'FINITO'
read(*,*)

end program testSMASH
