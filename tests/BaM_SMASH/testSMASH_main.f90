program testSMASH

use kinds_dmsl_kit
use SMASH_model
use types_dmsl_kit,only:data_ricz_type

implicit none

character(len_vLongStr), parameter:: setupFile='../../tests/BaM_SMASH/Config_setup.txt'
integer(mik):: err
character(250):: mess
type(data_ricz_type):: xtra
real(mrk)::Y(1460,1)
logical::feas(1)

call SMASH_XtraRead(file=trim(setupFile),xtra=xtra,err=err,mess=mess)
if(err/=0) then;write(*,*) trim(mess);read(*,*);endif

call SMASH_Load(loadScript=xtra%cp1(1),projectDir=xtra%cp1(3),&
                precipDir=xtra%cp1(4),petDir=xtra%cp1(5),err=err,mess=mess)
if(err/=0) then;write(*,*) trim(mess);read(*,*);endif

call SMASH_Run(runScript=xtra%cp1(2), projectDir=xtra%cp1(3),thetaDir=xtra%cp1(6),&
               theta=1._mrk*(/1,500,800,200,100,1200,1,100,1,600,1,0/),&
               Y=Y,feas=feas,err=err,mess=mess)
if(err/=0) then;write(*,*) trim(mess);read(*,*);endif

write(*,*) 'Exit status:', err
write(*,*) 'Exit message:', mess

write(*,*) 'FINITO'
read(*,*)

end program testSMASH
