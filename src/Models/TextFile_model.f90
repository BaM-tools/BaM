!~**********************************************************************
!~* Purpose: Model from a formula written in a text file
!~**********************************************************************
!~* Programmer: Ben Renard, Irstea Lyon
!~**********************************************************************
!~* Created: 18/06/2019
!~**********************************************************************
!~* Comments: Built on the fortran expression parser by Roland Schmehl.
!~*           See module fparser below for copyright and license info.
!~*
!~*           This f90 files starts with Roland Schmehl's modules
!~*           parameters and fparser. I just made a few modifications
!~*           (ctrl+F BEN RENARD to see them)
!~*
!~*           The module TextFile_model itself is at the end of this
!~*           file (ctrl+F module TextFile_model).
!~**********************************************************************
!~* References:
!~**********************************************************************
!~* 2Do List:
!~**********************************************************************
!~* Quick description of public procedures:
!~*     1. XXX
!~*     2. XXX
!~*     3. XXX
!~**********************************************************************

!oooooooooooooooooooooooooooooooooooooooooooooo
!o!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!o
!o ROLAND SCHMEHL'S FORTRAN EXPRESSION PARSER o
!o!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!o
!oooooooooooooooooooooooooooooooooooooooooooooo

! Roland Schmehl parser dowloaded here: https://sourceforge.net/projects/fparser/
! object-oriented versions of this can be found at:
! https://github.com/jacobwilliams/fortran_function_parser
! https://github.com/jacopo-chevallard/FortranParser (memory-leak? see reported issues)

! Copyright (c) 2000-2008, Roland Schmehl. All rights reserved.
!
! This software is distributable under the BSD license. See the terms of the
! BSD license in the documentation provided with this software.
!
!Copyright (c) 2000-2008, Roland Schmehl.
!
!All rights reserved.
!
!* Redistribution and use in source and binary forms, with or without
!modification, are permitted provided that the following conditions are
!met:
!
!* Redistributions of source code must retain the above copyright notice,
!this list of conditions and the following disclaimer.
!
!* Redistributions in binary form must reproduce the above copyright
!notice, this list of conditions and the following disclaimer in the
!documentation and/or other materials provided with the distribution.
!
!* Neither the name of the copyright holder nor the names of its
!contributors may be used to endorse or promote products derived from
!this software without specific prior written permission.
!
!
!THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
!"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
!LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
!A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT
!OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
!SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
!LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
!DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
!THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
!(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
!OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!

MODULE parameters
  !--------- -------- --------- --------- --------- --------- --------- --------- -----
  ! Specify data types
  !--------- -------- --------- --------- --------- --------- --------- --------- -----
  use kinds_dmsl_kit ! ADDED BY BEN RENARD
  IMPLICIT NONE
  INTEGER, PARAMETER :: rn = mrk ! MODIFIED BY BEN RENARD ! KIND(0.0d0)          ! Precision of real numbers
  INTEGER, PARAMETER :: is = mik ! MODIFIED BY BEN RENARD ! SELECTED_INT_KIND(1) ! Data type of bytecode
END MODULE parameters

!
! Copyright (c) 2000-2008, Roland Schmehl. All rights reserved.
!
! This software is distributable under the BSD license. See the terms of the
! BSD license in the documentation provided with this software.
!
MODULE fparser
  !------- -------- --------- --------- --------- --------- --------- --------- -------
  ! Fortran 90 function parser v1.1
  !------- -------- --------- --------- --------- --------- --------- --------- -------
  !
  ! This function parser module is intended for applications where a set of mathematical
  ! fortran-style expressions is specified at runtime and is then evaluated for a large
  ! number of variable values. This is done by compiling the set of function strings
  ! into byte code, which is interpreted efficiently for the various variable values.
  !
  ! The source code is available from http://fparser.sourceforge.net
  !
  ! Please send comments, corrections or questions to the author:
  ! Roland Schmehl <roland.schmehl@alumni.uni-karlsruhe.de>
  !
  !------- -------- --------- --------- --------- --------- --------- --------- -------
  ! The function parser concept is based on a C++ class library written by  Juha
  ! Nieminen <warp@iki.fi> available from http://warp.povusers.org/FunctionParser/
  !------- -------- --------- --------- --------- --------- --------- --------- -------
  USE parameters, ONLY: rn,is               ! Import KIND parameters
  IMPLICIT NONE
  !------- -------- --------- --------- --------- --------- --------- --------- -------
  PUBLIC                     :: initf,    & ! Initialize function parser for n functions
                                parsef,   & ! Parse single function string
                                evalf,    & ! Evaluate single function
                                EvalErrMsg  ! Error message (Use only when EvalErrType>0)
  INTEGER, PUBLIC            :: EvalErrType ! =0: no error occured, >0: evaluation error
  !------- -------- --------- --------- --------- --------- --------- --------- -------
  PRIVATE
  SAVE
  ! MODIFIED BY BEN RENARD: Added mathematical functions (from is0 onwards)
  INTEGER(is),                              PARAMETER :: cImmed   = 1,          &
                                                         cNeg     = 2,          &
                                                         cAdd     = 3,          &
                                                         cSub     = 4,          &
                                                         cMul     = 5,          &
                                                         cDiv     = 6,          &
                                                         cPow     = 7,          &
                                                         cAbs     = 8,          &
                                                         cExp     = 9,          &
                                                         cLog10   = 10,         &
                                                         cLog     = 11,         &
                                                         cSqrt    = 12,         &
                                                         cSinh    = 13,         &
                                                         cCosh    = 14,         &
                                                         cTanh    = 15,         &
                                                         cSin     = 16,         &
                                                         cCos     = 17,         &
                                                         cTan     = 18,         &
                                                         cAsin    = 19,         &
                                                         cAcos    = 20,         &
                                                         cAtan    = 21,         &
                                                         cIs0     = 22,         &
                                                         cIsPos   = 23,         &
                                                         cIsNeg   = 24,         &
                                                         cIsSPos  = 25,         &
                                                         cIsSNeg  = 26,         &
                                                         VarBegin = 27
  CHARACTER (LEN=1), DIMENSION(cAdd:cPow),  PARAMETER :: Ops      = (/ '+',     &
                                                                       '-',     &
                                                                       '*',     &
                                                                       '/',     &
                                                                       '^' /)

  ! MODIFIED BY BEN RENARD: defined begin/end of math functions
  integer(is),parameter::mathFirst=cAbs,mathLast=cIsSNeg
  CHARACTER (LEN=6), DIMENSION(mathFirst:mathLast), PARAMETER :: Funcs   = (/&
                                                                       'abs   ', &
                                                                       'exp   ', &
                                                                       'log10 ', &
                                                                       'log   ', &
                                                                       'sqrt  ', &
                                                                       'sinh  ', &
                                                                       'cosh  ', &
                                                                       'tanh  ', &
                                                                       'sin   ', &
                                                                       'cos   ', &
                                                                       'tan   ', &
                                                                       'asin  ', &
                                                                       'acos  ', &
                                                                       'atan  ', &
                                                                       'is0   ', &
                                                                       'ispos ', &
                                                                       'isneg ', &
                                                                       'isspos', &
                                                                       'issneg' &
                                                                       /)
  TYPE tComp
     INTEGER(is), DIMENSION(:), POINTER :: ByteCode
     INTEGER                            :: ByteCodeSize
     REAL(rn),    DIMENSION(:), POINTER :: Immed
     INTEGER                            :: ImmedSize
     REAL(rn),    DIMENSION(:), POINTER :: Stack
     INTEGER                            :: StackSize, &
                                           StackPtr
  END TYPE tComp
  TYPE (tComp),  DIMENSION(:),  POINTER :: Comp              ! Bytecode
  INTEGER,   DIMENSION(:),  ALLOCATABLE :: ipos              ! Associates function strings
  !
CONTAINS
  !
  SUBROUTINE initf (n)
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    ! Initialize function parser for n functions
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    IMPLICIT NONE
    INTEGER, INTENT(in) :: n                                 ! Number of functions
    INTEGER             :: i
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    ALLOCATE (Comp(n))
    DO i=1,n
       NULLIFY (Comp(i)%ByteCode,Comp(i)%Immed,Comp(i)%Stack)
    END DO
  END SUBROUTINE initf
  !
  SUBROUTINE parsef (i, FuncStr, Var)
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    ! Parse ith function string FuncStr and compile it into bytecode
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    IMPLICIT NONE
    INTEGER,                         INTENT(in) :: i         ! Function identifier
    CHARACTER (LEN=*),               INTENT(in) :: FuncStr   ! Function string
    CHARACTER (LEN=*), DIMENSION(:), INTENT(in) :: Var       ! Array with variable names
    CHARACTER (LEN=LEN(FuncStr))                :: Func      ! Function string, local use
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    IF (i < 1 .OR. i > SIZE(Comp)) THEN
       WRITE(*,*) '*** Parser error: Function number ',i,' out of range'
       STOP
    END IF
    ALLOCATE (ipos(LEN_TRIM(FuncStr)))                       ! Char. positions in orig. string
    Func = FuncStr                                           ! Local copy of function string
    CALL Replace ('**','^ ',Func)                            ! Exponent into 1-Char. format
    CALL RemoveSpaces (Func)                                 ! Condense function string
    CALL CheckSyntax (Func,FuncStr,Var)
    DEALLOCATE (ipos)
    CALL Compile (i,Func,Var)                                ! Compile into bytecode
  END SUBROUTINE parsef
  !
  FUNCTION evalf (i, Val) RESULT (res)
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    ! Evaluate bytecode of ith function for the values passed in array Val(:)
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    IMPLICIT NONE
    INTEGER,                INTENT(in) :: i                  ! Function identifier
    REAL(rn), DIMENSION(:), INTENT(in) :: Val                ! Variable values
    REAL(rn)                           :: res                ! Result
    INTEGER                            :: IP,              & ! Instruction pointer
                                          DP,              & ! Data pointer
                                          SP                 ! Stack pointer
    REAL(rn),                PARAMETER :: zero = 0._rn
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    DP = 1
    SP = 0
    DO IP=1,Comp(i)%ByteCodeSize
       SELECT CASE (Comp(i)%ByteCode(IP))

       CASE (cImmed); SP=SP+1; Comp(i)%Stack(SP)=Comp(i)%Immed(DP); DP=DP+1
       CASE   (cNeg); Comp(i)%Stack(SP)=-Comp(i)%Stack(SP)
       CASE   (cAdd); Comp(i)%Stack(SP-1)=Comp(i)%Stack(SP-1)+Comp(i)%Stack(SP); SP=SP-1
       CASE   (cSub); Comp(i)%Stack(SP-1)=Comp(i)%Stack(SP-1)-Comp(i)%Stack(SP); SP=SP-1
       CASE   (cMul); Comp(i)%Stack(SP-1)=Comp(i)%Stack(SP-1)*Comp(i)%Stack(SP); SP=SP-1
       CASE   (cDiv); IF (Comp(i)%Stack(SP)==0._rn) THEN; EvalErrType=1; res=zero; RETURN; ENDIF
                      Comp(i)%Stack(SP-1)=Comp(i)%Stack(SP-1)/Comp(i)%Stack(SP); SP=SP-1
       CASE   (cPow); Comp(i)%Stack(SP-1)=Comp(i)%Stack(SP-1)**Comp(i)%Stack(SP); SP=SP-1
       CASE   (cAbs); Comp(i)%Stack(SP)=ABS(Comp(i)%Stack(SP))
       CASE   (cExp); Comp(i)%Stack(SP)=EXP(Comp(i)%Stack(SP))
       CASE (cLog10); IF (Comp(i)%Stack(SP)<=0._rn) THEN; EvalErrType=3; res=zero; RETURN; ENDIF
                      Comp(i)%Stack(SP)=LOG10(Comp(i)%Stack(SP))
       CASE   (cLog); IF (Comp(i)%Stack(SP)<=0._rn) THEN; EvalErrType=3; res=zero; RETURN; ENDIF
                      Comp(i)%Stack(SP)=LOG(Comp(i)%Stack(SP))
       CASE  (cSqrt); IF (Comp(i)%Stack(SP)<0._rn) THEN; EvalErrType=3; res=zero; RETURN; ENDIF
                      Comp(i)%Stack(SP)=SQRT(Comp(i)%Stack(SP))
       CASE  (cSinh); Comp(i)%Stack(SP)=SINH(Comp(i)%Stack(SP))
       CASE  (cCosh); Comp(i)%Stack(SP)=COSH(Comp(i)%Stack(SP))
       CASE  (cTanh); Comp(i)%Stack(SP)=TANH(Comp(i)%Stack(SP))
       CASE   (cSin); Comp(i)%Stack(SP)=SIN(Comp(i)%Stack(SP))
       CASE   (cCos); Comp(i)%Stack(SP)=COS(Comp(i)%Stack(SP))
       CASE   (cTan); Comp(i)%Stack(SP)=TAN(Comp(i)%Stack(SP))
       CASE  (cAsin); IF ((Comp(i)%Stack(SP)<-1._rn).OR.(Comp(i)%Stack(SP)>1._rn)) THEN
                      EvalErrType=4; res=zero; RETURN; ENDIF
                      Comp(i)%Stack(SP)=ASIN(Comp(i)%Stack(SP))
       CASE  (cAcos); IF ((Comp(i)%Stack(SP)<-1._rn).OR.(Comp(i)%Stack(SP)>1._rn)) THEN
                      EvalErrType=4; res=zero; RETURN; ENDIF
                      Comp(i)%Stack(SP)=ACOS(Comp(i)%Stack(SP))
       CASE  (cAtan); Comp(i)%Stack(SP)=ATAN(Comp(i)%Stack(SP))
       ! MODIFIED BY BEN RENARD: Added functions
       CASE   (cIs0); IF (Comp(i)%Stack(SP)==zero) THEN
                          Comp(i)%Stack(SP)=1._rn
                      ELSE
                          Comp(i)%Stack(SP)=0._rn
                      ENDIF
       CASE (cIsPos); IF (Comp(i)%Stack(SP)>=zero) THEN
                          Comp(i)%Stack(SP)=1._rn
                      ELSE
                          Comp(i)%Stack(SP)=0._rn
                      ENDIF
       CASE (cIsNeg); IF (Comp(i)%Stack(SP)<=zero) THEN
                          Comp(i)%Stack(SP)=1._rn
                      ELSE
                          Comp(i)%Stack(SP)=0._rn
                      ENDIF
       CASE (cIsSPos); IF (Comp(i)%Stack(SP)>zero) THEN
                          Comp(i)%Stack(SP)=1._rn
                      ELSE
                          Comp(i)%Stack(SP)=0._rn
                      ENDIF
       CASE (cIsSNeg); IF (Comp(i)%Stack(SP)<zero) THEN
                          Comp(i)%Stack(SP)=1._rn
                      ELSE
                          Comp(i)%Stack(SP)=0._rn
                      ENDIF
       CASE  DEFAULT; SP=SP+1; Comp(i)%Stack(SP)=Val(Comp(i)%ByteCode(IP)-VarBegin+1)
       END SELECT
    END DO
    EvalErrType = 0
    res = Comp(i)%Stack(1)
  END FUNCTION evalf
  !
  SUBROUTINE CheckSyntax (Func,FuncStr,Var)
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    ! Check syntax of function string,  returns 0 if syntax is ok
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    IMPLICIT NONE
    CHARACTER (LEN=*),               INTENT(in) :: Func      ! Function string without spaces
    CHARACTER (LEN=*),               INTENT(in) :: FuncStr   ! Original function string
    CHARACTER (LEN=*), DIMENSION(:), INTENT(in) :: Var       ! Array with variable names
    INTEGER(is)                                 :: n
    CHARACTER (LEN=1)                           :: c
    REAL(rn)                                    :: r
    LOGICAL                                     :: err
    INTEGER                                     :: ParCnt, & ! Parenthesis counter
                                                   j,ib,in,lFunc
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    j = 1
    ParCnt = 0
    lFunc = LEN_TRIM(Func)
    step: DO
       IF (j > lFunc) CALL ParseErrMsg (j, FuncStr)
       c = Func(j:j)
       !-- -------- --------- --------- --------- --------- --------- --------- -------
       ! Check for valid operand (must appear)
       !-- -------- --------- --------- --------- --------- --------- --------- -------
       IF (c == '-' .OR. c == '+') THEN                      ! Check for leading - or +
          j = j+1
          IF (j > lFunc) CALL ParseErrMsg (j, FuncStr, 'Missing operand')
          c = Func(j:j)
          IF (ANY(c == Ops)) CALL ParseErrMsg (j, FuncStr, 'Multiple operators')
       END IF
       n = MathFunctionIndex (Func(j:))
       IF (n > 0) THEN                                       ! Check for math function
          j = j+LEN_TRIM(Funcs(n))
          IF (j > lFunc) CALL ParseErrMsg (j, FuncStr, 'Missing function argument')
          c = Func(j:j)
          IF (c /= '(') CALL ParseErrMsg (j, FuncStr, 'Missing opening parenthesis')
       END IF
       IF (c == '(') THEN                                    ! Check for opening parenthesis
          ParCnt = ParCnt+1
          j = j+1
          CYCLE step
       END IF
       IF (SCAN(c,'0123456789.') > 0) THEN                   ! Check for number
          r = RealNum (Func(j:),ib,in,err)
          IF (err) CALL ParseErrMsg (j, FuncStr, 'Invalid number format:  '//Func(j+ib-1:j+in-2))
          j = j+in-1
          IF (j > lFunc) EXIT
          c = Func(j:j)
       ELSE                                                  ! Check for variable
          n = VariableIndex (Func(j:),Var,ib,in)
          IF (n == 0) CALL ParseErrMsg (j, FuncStr, 'Invalid element: '//Func(j+ib-1:j+in-2))
          j = j+in-1
          IF (j > lFunc) EXIT
          c = Func(j:j)
       END IF
       DO WHILE (c == ')')                                   ! Check for closing parenthesis
          ParCnt = ParCnt-1
          IF (ParCnt < 0) CALL ParseErrMsg (j, FuncStr, 'Mismatched parenthesis')
          IF (Func(j-1:j-1) == '(') CALL ParseErrMsg (j-1, FuncStr, 'Empty parentheses')
          j = j+1
          IF (j > lFunc) EXIT
          c = Func(j:j)
       END DO
       !-- -------- --------- --------- --------- --------- --------- --------- -------
       ! Now, we have a legal operand: A legal operator or end of string must follow
       !-- -------- --------- --------- --------- --------- --------- --------- -------
       IF (j > lFunc) EXIT
       IF (ANY(c == Ops)) THEN                               ! Check for multiple operators
          IF (j+1 > lFunc) CALL ParseErrMsg (j, FuncStr)
          IF (ANY(Func(j+1:j+1) == Ops)) CALL ParseErrMsg (j+1, FuncStr, 'Multiple operators')
       ELSE                                                  ! Check for next operand
          CALL ParseErrMsg (j, FuncStr, 'Missing operator')
       END IF
       !-- -------- --------- --------- --------- --------- --------- --------- -------
       ! Now, we have an operand and an operator: the next loop will check for another
       ! operand (must appear)
       !-- -------- --------- --------- --------- --------- --------- --------- -------
       j = j+1
    END DO step
    IF (ParCnt > 0) CALL ParseErrMsg (j, FuncStr, 'Missing )')
  END SUBROUTINE CheckSyntax
  !
  FUNCTION EvalErrMsg () RESULT (msg)
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    ! Return error message
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    IMPLICIT NONE
    ! MODIFIED BY BEN RENARD: add error message for syntax error
    !CHARACTER (LEN=*), DIMENSION(4), PARAMETER :: m = (/ 'Division by zero                ', &
    !                                                     'Argument of SQRT negative       ', &
    !                                                     'Argument of LOG negative        ', &
    !                                                     'Argument of ASIN or ACOS illegal' /)
    CHARACTER (LEN=*), DIMENSION(5), PARAMETER :: m = (/ 'Division by zero                ', &
                                                         'Argument of SQRT negative       ', &
                                                         'Argument of LOG negative        ', &
                                                         'Argument of ASIN or ACOS illegal', &
                                                         'Syntax error                    '/)
    CHARACTER (LEN=LEN(m))                     :: msg
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    IF (EvalErrType < 1 .OR. EvalErrType > SIZE(m)) THEN
       msg = ''
    ELSE
       msg = m(EvalErrType)
    ENDIF
  END FUNCTION EvalErrMsg
  !
  SUBROUTINE ParseErrMsg (j, FuncStr, Msg)
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    ! Print error message and terminate program
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    IMPLICIT NONE
    INTEGER,                     INTENT(in) :: j
    CHARACTER (LEN=*),           INTENT(in) :: FuncStr       ! Original function string
    CHARACTER (LEN=*), OPTIONAL, INTENT(in) :: Msg
    INTEGER                                 :: k
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    IF (PRESENT(Msg)) THEN
       WRITE(*,*) '*** Error in syntax of function string: '//Msg
    ELSE
       WRITE(*,*) '*** Error in syntax of function string:'
    ENDIF
    WRITE(*,*)
    WRITE(*,'(A)') ' '//FuncStr
    DO k=1,ipos(j)
       WRITE(*,'(A)',ADVANCE='NO') ' '                       ! Advance to the jth position
    END DO
    WRITE(*,'(A)') '?'
    ! MODIFIED BY BEN RENARD: Do not stop but rather return an error message
    ! STOP
    EvalErrType=5
  END SUBROUTINE ParseErrMsg
  !
  FUNCTION OperatorIndex (c) RESULT (n)
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    ! Return operator index
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    IMPLICIT NONE
    CHARACTER (LEN=1), INTENT(in) :: c
    INTEGER(is)                   :: n,j
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    n = 0
    DO j=cAdd,cPow
       IF (c == Ops(j)) THEN
          n = j
          EXIT
       END IF
    END DO
  END FUNCTION OperatorIndex
  !
  FUNCTION MathFunctionIndex (str) RESULT (n)
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    ! Return index of math function beginnig at 1st position of string str
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    IMPLICIT NONE
    CHARACTER (LEN=*), INTENT(in) :: str
    INTEGER(is)                   :: n,j
    INTEGER                       :: k
    CHARACTER (LEN=LEN(Funcs))    :: fun
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    n = 0
    ! MODIFIED BY BEN RENARD: replaced cAbs/cAtan by mathFirst/mathLast to make includion of new functions easier
    DO j=mathFirst,mathLast                                  ! Check all math functions
       k = MIN(LEN_TRIM(Funcs(j)), LEN(str))
       CALL LowCase (str(1:k), fun)
       IF (fun == Funcs(j)) THEN                             ! Compare lower case letters
          n = j                                              ! Found a matching function
          EXIT
       END IF
    END DO
  END FUNCTION MathFunctionIndex
  !
  FUNCTION VariableIndex (str, Var, ibegin, inext) RESULT (n)
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    ! Return index of variable at begin of string str (returns 0 if no variable found)
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    IMPLICIT NONE
    CHARACTER (LEN=*),               INTENT(in) :: str       ! String
    CHARACTER (LEN=*), DIMENSION(:), INTENT(in) :: Var       ! Array with variable names
    INTEGER(is)                                 :: n         ! Index of variable
    INTEGER, OPTIONAL,              INTENT(out) :: ibegin, & ! Start position of variable name
                                                   inext     ! Position of character after name
    INTEGER                                     :: j,ib,in,lstr
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    n = 0
    lstr = LEN_TRIM(str)
    IF (lstr > 0) THEN
       DO ib=1,lstr                                          ! Search for first character in str
          IF (str(ib:ib) /= ' ') EXIT                        ! When lstr>0 at least 1 char in str
       END DO
       DO in=ib,lstr                                         ! Search for name terminators
          IF (SCAN(str(in:in),'+-*/^) ') > 0) EXIT
       END DO
       DO j=1,SIZE(Var)
          IF (str(ib:in-1) == Var(j)) THEN
             n = j                                           ! Variable name found
             EXIT
          END IF
       END DO
    END IF
    IF (PRESENT(ibegin)) ibegin = ib
    IF (PRESENT(inext))  inext  = in
  END FUNCTION VariableIndex
  !
  SUBROUTINE RemoveSpaces (str)
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    ! Remove Spaces from string, remember positions of characters in old string
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    IMPLICIT NONE
    CHARACTER (LEN=*), INTENT(inout) :: str
    INTEGER                          :: k,lstr
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    lstr = LEN_TRIM(str)
    ipos = (/ (k,k=1,lstr) /)
    k = 1
    DO WHILE (str(k:lstr) /= ' ')
       IF (str(k:k) == ' ') THEN
          str(k:lstr)  = str(k+1:lstr)//' '                  ! Move 1 character to left
          ipos(k:lstr) = (/ ipos(k+1:lstr), 0 /)             ! Move 1 element to left
          k = k-1
       END IF
       k = k+1
    END DO
  END SUBROUTINE RemoveSpaces
  !
  SUBROUTINE Replace (ca,cb,str)
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    ! Replace ALL appearances of character set ca in string str by character set cb
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    IMPLICIT NONE
    CHARACTER (LEN=*),       INTENT(in) :: ca
    CHARACTER (LEN=LEN(ca)), INTENT(in) :: cb                ! LEN(ca) must be LEN(cb)
    CHARACTER (LEN=*),    INTENT(inout) :: str
    INTEGER                             :: j,lca
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    lca = LEN(ca)
    DO j=1,LEN_TRIM(str)-lca+1
       IF (str(j:j+lca-1) == ca) str(j:j+lca-1) = cb
    END DO
  END SUBROUTINE Replace
  !
  SUBROUTINE Compile (i, F, Var)
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    ! Compile i-th function string F into bytecode
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    IMPLICIT NONE
    INTEGER,                         INTENT(in) :: i         ! Function identifier
    CHARACTER (LEN=*),               INTENT(in) :: F         ! Function string
    CHARACTER (LEN=*), DIMENSION(:), INTENT(in) :: Var       ! Array with variable names
    INTEGER                                     :: istat
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    IF (ASSOCIATED(Comp(i)%ByteCode)) DEALLOCATE ( Comp(i)%ByteCode, &
                                                   Comp(i)%Immed,    &
                                                   Comp(i)%Stack     )
    Comp(i)%ByteCodeSize = 0
    Comp(i)%ImmedSize    = 0
    Comp(i)%StackSize    = 0
    Comp(i)%StackPtr     = 0
    CALL CompileSubstr (i,F,1,LEN_TRIM(F),Var)               ! Compile string to determine size
    ALLOCATE ( Comp(i)%ByteCode(Comp(i)%ByteCodeSize), &
               Comp(i)%Immed(Comp(i)%ImmedSize),       &
               Comp(i)%Stack(Comp(i)%StackSize),       &
               STAT = istat                            )
    IF (istat /= 0) THEN
       WRITE(*,*) '*** Parser error: Memmory allocation for byte code failed'
       STOP
    ELSE
       Comp(i)%ByteCodeSize = 0
       Comp(i)%ImmedSize    = 0
       Comp(i)%StackSize    = 0
       Comp(i)%StackPtr     = 0
       CALL CompileSubstr (i,F,1,LEN_TRIM(F),Var)            ! Compile string into bytecode
    END IF
    !
  END SUBROUTINE Compile
  !
  SUBROUTINE AddCompiledByte (i, b)
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    ! Add compiled byte to bytecode
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    IMPLICIT NONE
    INTEGER,     INTENT(in) :: i                             ! Function identifier
    INTEGER(is), INTENT(in) :: b                             ! Value of byte to be added
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    Comp(i)%ByteCodeSize = Comp(i)%ByteCodeSize + 1
    IF (ASSOCIATED(Comp(i)%ByteCode)) Comp(i)%ByteCode(Comp(i)%ByteCodeSize) = b
  END SUBROUTINE AddCompiledByte
  !
  FUNCTION MathItemIndex (i, F, Var) RESULT (n)
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    ! Return math item index, if item is real number, enter it into Comp-structure
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    IMPLICIT NONE
    INTEGER,                         INTENT(in) :: i         ! Function identifier
    CHARACTER (LEN=*),               INTENT(in) :: F         ! Function substring
    CHARACTER (LEN=*), DIMENSION(:), INTENT(in) :: Var       ! Array with variable names
    INTEGER(is)                                 :: n         ! Byte value of math item
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    n = 0
    IF (SCAN(F(1:1),'0123456789.') > 0) THEN                 ! Check for begin of a number
       Comp(i)%ImmedSize = Comp(i)%ImmedSize + 1
       IF (ASSOCIATED(Comp(i)%Immed)) Comp(i)%Immed(Comp(i)%ImmedSize) = RealNum (F)
       n = cImmed
    ELSE                                                     ! Check for a variable
       n = VariableIndex (F, Var)
       IF (n > 0) n = VarBegin+n-1
    END IF
  END FUNCTION MathItemIndex
  !
  FUNCTION CompletelyEnclosed (F, b, e) RESULT (res)
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    ! Check if function substring F(b:e) is completely enclosed by a pair of parenthesis
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    IMPLICIT NONE
    CHARACTER (LEN=*), INTENT(in) :: F                       ! Function substring
    INTEGER,           INTENT(in) :: b,e                     ! First and last pos. of substring
    LOGICAL                       :: res
    INTEGER                       :: j,k
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    res=.false.
    IF (F(b:b) == '(' .AND. F(e:e) == ')') THEN
       k = 0
       DO j=b+1,e-1
          IF     (F(j:j) == '(') THEN
             k = k+1
          ELSEIF (F(j:j) == ')') THEN
             k = k-1
          END IF
          IF (k < 0) EXIT
       END DO
       IF (k == 0) res=.true.                                ! All opened parenthesis closed
    END IF
  END FUNCTION CompletelyEnclosed
  !
  RECURSIVE SUBROUTINE CompileSubstr (i, F, b, e, Var)
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    ! Compile i-th function string F into bytecode
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    IMPLICIT NONE
    INTEGER,                         INTENT(in) :: i         ! Function identifier
    CHARACTER (LEN=*),               INTENT(in) :: F         ! Function substring
    INTEGER,                         INTENT(in) :: b,e       ! Begin and end position substring
    CHARACTER (LEN=*), DIMENSION(:), INTENT(in) :: Var       ! Array with variable names
    INTEGER(is)                                 :: n
    INTEGER                                     :: b2,j,k,io
    CHARACTER (LEN=*),                PARAMETER :: calpha = 'abcdefghijklmnopqrstuvwxyz'// &
                                                            'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    ! Check for special cases of substring
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    IF     (F(b:b) == '+') THEN                              ! Case 1: F(b:e) = '+...'
!      WRITE(*,*)'1. F(b:e) = "+..."'
       CALL CompileSubstr (i, F, b+1, e, Var)
       RETURN
    ELSEIF (CompletelyEnclosed (F, b, e)) THEN               ! Case 2: F(b:e) = '(...)'
!      WRITE(*,*)'2. F(b:e) = "(...)"'
       CALL CompileSubstr (i, F, b+1, e-1, Var)
       RETURN
    ELSEIF (SCAN(F(b:b),calpha) > 0) THEN
       n = MathFunctionIndex (F(b:e))
       IF (n > 0) THEN
          b2 = b+INDEX(F(b:e),'(')-1
          IF (CompletelyEnclosed(F, b2, e)) THEN             ! Case 3: F(b:e) = 'fcn(...)'
!            WRITE(*,*)'3. F(b:e) = "fcn(...)"'
             CALL CompileSubstr(i, F, b2+1, e-1, Var)
             CALL AddCompiledByte (i, n)
             RETURN
          END IF
       END IF
    ELSEIF (F(b:b) == '-') THEN
       IF (CompletelyEnclosed (F, b+1, e)) THEN              ! Case 4: F(b:e) = '-(...)'
!         WRITE(*,*)'4. F(b:e) = "-(...)"'
          CALL CompileSubstr (i, F, b+2, e-1, Var)
          CALL AddCompiledByte (i, cNeg)
          RETURN
       ELSEIF (SCAN(F(b+1:b+1),calpha) > 0) THEN
          n = MathFunctionIndex (F(b+1:e))
          IF (n > 0) THEN
             b2 = b+INDEX(F(b+1:e),'(')
             IF (CompletelyEnclosed(F, b2, e)) THEN          ! Case 5: F(b:e) = '-fcn(...)'
!               WRITE(*,*)'5. F(b:e) = "-fcn(...)"'
                CALL CompileSubstr(i, F, b2+1, e-1, Var)
                CALL AddCompiledByte (i, n)
                CALL AddCompiledByte (i, cNeg)
                RETURN
             END IF
          END IF
       ENDIF
    END IF
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    ! Check for operator in substring: check only base level (k=0), exclude expr. in ()
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    DO io=cAdd,cPow                                          ! Increasing priority +-*/^
       k = 0
       DO j=e,b,-1
          IF     (F(j:j) == ')') THEN
             k = k+1
          ELSEIF (F(j:j) == '(') THEN
             k = k-1
          END IF
          IF (k == 0 .AND. F(j:j) == Ops(io) .AND. IsBinaryOp (j, F)) THEN
             IF (ANY(F(j:j) == Ops(cMul:cPow)) .AND. F(b:b) == '-') THEN ! Case 6: F(b:e) = '-...Op...' with Op > -
!               WRITE(*,*)'6. F(b:e) = "-...Op..." with Op > -'
                CALL CompileSubstr (i, F, b+1, e, Var)
                CALL AddCompiledByte (i, cNeg)
                RETURN
             ELSE                                                        ! Case 7: F(b:e) = '...BinOp...'
!               WRITE(*,*)'7. Binary operator',F(j:j)
                CALL CompileSubstr (i, F, b, j-1, Var)
                CALL CompileSubstr (i, F, j+1, e, Var)
                CALL AddCompiledByte (i, OperatorIndex(Ops(io)))
                Comp(i)%StackPtr = Comp(i)%StackPtr - 1
                RETURN
             END IF
          END IF
       END DO
    END DO
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    ! Check for remaining items, i.e. variables or explicit numbers
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    b2 = b
    IF (F(b:b) == '-') b2 = b2+1
    n = MathItemIndex(i, F(b2:e), Var)
!   WRITE(*,*)'8. AddCompiledByte ',n
    CALL AddCompiledByte (i, n)
    Comp(i)%StackPtr = Comp(i)%StackPtr + 1
    IF (Comp(i)%StackPtr > Comp(i)%StackSize) Comp(i)%StackSize = Comp(i)%StackSize + 1
    IF (b2 > b) CALL AddCompiledByte (i, cNeg)
  END SUBROUTINE CompileSubstr
  !
  FUNCTION IsBinaryOp (j, F) RESULT (res)
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    ! Check if operator F(j:j) in string F is binary operator
    ! Special cases already covered elsewhere:              (that is corrected in v1.1)
    ! - operator character F(j:j) is first character of string (j=1)
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    IMPLICIT NONE
    INTEGER,           INTENT(in) :: j                       ! Position of Operator
    CHARACTER (LEN=*), INTENT(in) :: F                       ! String
    LOGICAL                       :: res                     ! Result
    INTEGER                       :: k
    LOGICAL                       :: Dflag,Pflag
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    res=.true.
    IF (F(j:j) == '+' .OR. F(j:j) == '-') THEN               ! Plus or minus sign:
       IF (j == 1) THEN                                      ! - leading unary operator ?
          res = .false.
       ELSEIF (SCAN(F(j-1:j-1),'+-*/^(') > 0) THEN           ! - other unary operator ?
          res = .false.
       ELSEIF (SCAN(F(j+1:j+1),'0123456789') > 0 .AND. &     ! - in exponent of real number ?
               SCAN(F(j-1:j-1),'eEdD')       > 0) THEN
          Dflag=.false.; Pflag=.false.
          k = j-1
          DO WHILE (k > 1)                                   !   step to the left in mantissa
             k = k-1
             IF     (SCAN(F(k:k),'0123456789') > 0) THEN
                Dflag=.true.
             ELSEIF (F(k:k) == '.') THEN
                IF (Pflag) THEN
                   EXIT                                      !   * EXIT: 2nd appearance of '.'
                ELSE
                   Pflag=.true.                              !   * mark 1st appearance of '.'
                ENDIF
             ELSE
                EXIT                                         !   * all other characters
             END IF
          END DO
          IF (Dflag .AND. (k == 1 .OR. SCAN(F(k:k),'+-*/^(') > 0)) res = .false.
       END IF
    END IF
  END FUNCTION IsBinaryOp
  !
  FUNCTION RealNum (str, ibegin, inext, error) RESULT (res)
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    ! Get real number from string - Format: [blanks][+|-][nnn][.nnn][e|E|d|D[+|-]nnn]
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    IMPLICIT NONE
    CHARACTER (LEN=*),  INTENT(in) :: str                    ! String
    REAL(rn)                       :: res                    ! Real number
    INTEGER, OPTIONAL, INTENT(out) :: ibegin,              & ! Start position of real number
                                      inext                  ! 1st character after real number
    LOGICAL, OPTIONAL, INTENT(out) :: error                  ! Error flag
    INTEGER                        :: ib,in,istat
    LOGICAL                        :: Bflag,               & ! .T. at begin of number in str
                                      InMan,               & ! .T. in mantissa of number
                                      Pflag,               & ! .T. after 1st '.' encountered
                                      Eflag,               & ! .T. at exponent identifier 'eEdD'
                                      InExp,               & ! .T. in exponent of number
                                      DInMan,              & ! .T. if at least 1 digit in mant.
                                      DInExp,              & ! .T. if at least 1 digit in exp.
                                      err                    ! Local error flag
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    Bflag=.true.; InMan=.false.; Pflag=.false.; Eflag=.false.; InExp=.false.
    DInMan=.false.; DInExp=.false.
    ib   = 1
    in   = 1
    DO WHILE (in <= LEN_TRIM(str))
       SELECT CASE (str(in:in))
       CASE (' ')                                            ! Only leading blanks permitted
          ib = ib+1
          IF (InMan .OR. Eflag .OR. InExp) EXIT
       CASE ('+','-')                                        ! Permitted only
          IF     (Bflag) THEN
             InMan=.true.; Bflag=.false.                     ! - at beginning of mantissa
          ELSEIF (Eflag) THEN
             InExp=.true.; Eflag=.false.                     ! - at beginning of exponent
          ELSE
             EXIT                                            ! - otherwise STOP
          ENDIF
       CASE ('0':'9')                                        ! Mark
          IF     (Bflag) THEN
             InMan=.true.; Bflag=.false.                     ! - beginning of mantissa
          ELSEIF (Eflag) THEN
             InExp=.true.; Eflag=.false.                     ! - beginning of exponent
          ENDIF
          IF (InMan) DInMan=.true.                           ! Mantissa contains digit
          IF (InExp) DInExp=.true.                           ! Exponent contains digit
       CASE ('.')
          IF     (Bflag) THEN
             Pflag=.true.                                    ! - mark 1st appearance of '.'
             InMan=.true.; Bflag=.false.                     !   mark beginning of mantissa
          ELSEIF (InMan .AND..NOT.Pflag) THEN
             Pflag=.true.                                    ! - mark 1st appearance of '.'
          ELSE
             EXIT                                            ! - otherwise STOP
          END IF
       CASE ('e','E','d','D')                                ! Permitted only
          IF (InMan) THEN
             Eflag=.true.; InMan=.false.                     ! - following mantissa
          ELSE
             EXIT                                            ! - otherwise STOP
          ENDIF
       CASE DEFAULT
          EXIT                                               ! STOP at all other characters
       END SELECT
       in = in+1
    END DO
    err = (ib > in-1) .OR. (.NOT.DInMan) .OR. ((Eflag.OR.InExp).AND..NOT.DInExp)
    IF (err) THEN
       res = 0.0_rn
    ELSE
       READ(str(ib:in-1),*,IOSTAT=istat) res
       err = istat /= 0
    END IF
    IF (PRESENT(ibegin)) ibegin = ib
    IF (PRESENT(inext))  inext  = in
    IF (PRESENT(error))  error  = err
  END FUNCTION RealNum
  !
  SUBROUTINE LowCase (str1, str2)
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    ! Transform upper case letters in str1 into lower case letters, result is str2
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    IMPLICIT NONE
    CHARACTER (LEN=*),  INTENT(in) :: str1
    CHARACTER (LEN=*), INTENT(out) :: str2
    INTEGER                        :: j,k
    CHARACTER (LEN=*),   PARAMETER :: lc = 'abcdefghijklmnopqrstuvwxyz'
    CHARACTER (LEN=*),   PARAMETER :: uc = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    !----- -------- --------- --------- --------- --------- --------- --------- -------
    str2 = str1
    DO j=1,LEN_TRIM(str1)
       k = INDEX(uc,str1(j:j))
       IF (k > 0) str2(j:j) = lc(k:k)
    END DO
  END SUBROUTINE LowCase
  !
END MODULE fparser

!oooooooooooooooooooooooooooooo
!o!!!!!!!!!!!!!!!!!!!!!!!!!!!!o
!o MAIN MODULE TextFile_model o
!o!!!!!!!!!!!!!!!!!!!!!!!!!!!!o
!oooooooooooooooooooooooooooooo

module TextFile_model

use kinds_dmsl_kit ! numeric kind definitions from DMSL
use fparser,only: initf, parsef, evalf, EvalErrType, EvalErrMsg ! expression parser by Roland Schmehl

implicit none
Private
public :: TxtMdl_load,TxtMdl_Apply,TxtMdl_define

Contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine TxtMdl_load(file,err,mess)
!^**********************************************************************
!^* Purpose: Read and load model from a text file
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Created: 18/06/2019
!^**********************************************************************
!^* Comments: the text file should be structured as follows:
!^*           1. nIN ! number of input variables
!^*           2. X1,X2,etc. ! list of input variables, comma-separated
!^*           3. nPAR ! number of parameters
!^*           4. P1,P2,etc. ! list of parameters, comma-separated
!^*           5. nOUT ! number of output variables
!^*           6. Formula 1 ! Formula for 1st output variable
!^*           7. Formula 2 ! Formula for 2nd output variable
!^*           8. etc. until nOUT is reached
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*     1. file, text file
!^* OUT
!^*     1. err, error code; <0:Warning, ==0:OK, >0: Error
!^*     2. mess, error message
!^**********************************************************************
use utilities_dmsl_kit,only:getSpareUnit,countSubstringInString

character(*), intent(in)::file
integer(mik), intent(out)::err
character(*),intent(out)::mess
! locals
character(250),parameter::procname='TxtMdl_load'
integer(mik)::unt,nIN,nPAR,nOUT,i,j,n
character(250),allocatable:: nameIN(:),namePAR(:)
character(len_uLongStr)::foo
character(len_uLongStr),allocatable:: formula(:)

err=0;mess=''

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 1. Read content of text file
call getSpareUnit(unt,err,mess)
if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif
open(unit=unt,file=trim(file), status='old', iostat=err)
if(err>0) then;mess=trim(procname)//':problem opening file '//trim(file);return;endif
! INPUT VARIABLES
read(unt,*,iostat=err) nIN ! number of input variables
if(err/=0) then;mess=trim(procname)//':problem reading file '//trim(file);return;endif
if(nIN<0) then;err=1;mess=trim(procname)//':negative nIN not allowed';return;endif
if(allocated(nameIN)) deallocate(nameIN); allocate(nameIN(nIN))
if(nIN>0) then
    read(unt,*,iostat=err) nameIN ! names of input variables
    if(err/=0) then;mess=trim(procname)//':problem reading file '//trim(file);return;endif
else
    read(unt,*)
endif
! PARAMETERS
read(unt,*,iostat=err) nPAR ! number of parameters
if(err/=0) then;mess=trim(procname)//':problem reading file '//trim(file);return;endif
if(nPAR<0) then;err=1;mess=trim(procname)//':negative nPAR not allowed';return;endif
if(allocated(namePAR)) deallocate(namePAR); allocate(namePAR(nPAR))
if(nPAR>0) then
    read(unt,*,iostat=err) namePAR ! names of input variables
    if(err/=0) then;mess=trim(procname)//':problem reading file '//trim(file);return;endif
else
    read(unt,*)
endif
! OUTPUT VARIABLES
read(unt,*,iostat=err) nOUT ! number of output variables
if(err/=0) then;mess=trim(procname)//':problem reading file '//trim(file);return;endif
if(nOUT<=0) then;err=1;mess=trim(procname)//':nOUT should be at least one';return;endif
if(allocated(formula)) deallocate(formula); allocate(formula(nOUT))
do i=1,nOUT
    read(unt,'(A)',iostat=err) foo ! formula for ith output, comment included
    n=len_trim(foo)
    do j=1,n
        if(foo(j:j)=='!') exit
    enddo
    formula(i)=foo(1:(j-1)) ! formula for ith output with comments removed
    if(err/=0) then;mess=trim(procname)//':problem reading file '//trim(file);return;endif
enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 2. Define model
call TxtMdl_define(nameIN,namePAR,formula,err,mess)
if(err/=0) then;mess=trim(procname)//':'//trim(mess);return;endif

end subroutine TxtMdl_load
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine TxtMdl_define(nameIN,namePAR,formulas,err,mess)
!^**********************************************************************
!^* Purpose: define (Init and parse) the text model
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Last modified: 04/07/2019
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*     1. nameIN, character vector, name of input variables
!^*     2. namePAR, character vector, name of parameters
!^*     3. formulas, character vector, text formulas
!^* OUT
!^*     1. err, error code; <0:Warning, ==0:OK, >0: Error
!^*     2. mess, error message
!^**********************************************************************
character(*), intent(in)::nameIN(:),namePAR(:),formulas(:)
integer(mik), intent(out)::err
character(*),intent(out)::mess
! locals
character(250),parameter::procname='TxtMdl_define'
integer(mik)::nIN,nPAR,nOUT,i
character(250):: varname(size(nameIN)+size(namePAR))

err=0;mess=''
nIN=size(nameIN)
nPAR=size(namePAR)
nOUT=size(formulas)
if(nIN>0) varname(1:nIN)=nameIN
varname( (nIN+1):(nIN+nPAR) )=namePAR
call initf(nOUT)
do i=1,nOUT
    call parsef(i,formulas(i),varname)
    if(EvalErrType>0) then
        err=EvalErrType;mess=trim(procname)//':'//EvalErrMsg();return
    endif
enddo

end subroutine TxtMdl_define
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine TxtMdl_Apply(IN,theta,whichOUT,OUT,feas,err,mess)
!^**********************************************************************
!^* Purpose: apply text-file model
!^**********************************************************************
!^* Programmer: Ben Renard, Irstea Lyon
!^**********************************************************************
!^* Created:18/06/2019
!^**********************************************************************
!^* Comments:
!^**********************************************************************
!^* References:
!^**********************************************************************
!^* 2Do List:
!^**********************************************************************
!^* IN
!^*     1. IN, input vector
!^*     2. theta, parameter vector
!^*     3. [whichOUT], integer, optional, which output should be computed? If absent, all outputs are computed.
!^* OUT
!^*     1. OUT, output vector
!^*     2. feas, feasible?
!^*     3. err, error code; <0:Warning, ==0:OK, >0: Error
!^*     4. mess, error message
!^**********************************************************************
real(mrk), intent(in)::IN(:,:),theta(:)
integer(mik),intent(in),optional::whichOUT
real(mrk), intent(out)::OUT(:,:)
logical, intent(out)::feas(:)
integer(mik), intent(out)::err
character(*),intent(out)::mess
!locals
character(250),parameter::procname='TxtMdl_Apply'
integer(mik)::nOBS,nIN,nOUT,nPAR,t,k
real(mrk)::varval(size(IN,dim=2)+size(theta))

err=0;mess='';feas=.true.;OUT=undefRN
nOBS=size(IN,dim=1)
nIN=size(IN,dim=2)
nPAR=size(theta)
varval( (nIN+1):(nIN+nPAR) )=theta
if(present(whichOUT)) then
    do t=1,nOBS
        if(nIN>0) varval(1:nIN)=IN(t,:)
        OUT(t,1)=evalf(whichOUT,varval)
        if(EvalErrType>0 .or. OUT(t,1)/=OUT(t,1)) then ! the second condition aims at detecting NAs
            feas(t)=.false.
            OUT(t,1)=undefRN
        endif
    enddo
else
    nOUT=size(OUT,dim=2)
    do t=1,nOBS
        if(nIN>0) varval(1:nIN)=IN(t,:)
        do k=1,nOUT
            OUT(t,k)=evalf(k,varval)
            if(EvalErrType>0 .or. OUT(t,k)/=OUT(t,k)) then
                feas(t)=.false.
                OUT(t,k)=undefRN
            endif
        enddo
    enddo
endif

end subroutine TxtMdl_Apply
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module TextFile_model

