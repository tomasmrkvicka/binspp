@echo off
rem Checks valid SHA sums of binspp package
rem (c) 2022 Radim Remeš, Ladislav Beránek

SET filename=binspp_0.1.15.tar.gz
SET sha1=8dffc4b12edd5c609cd9ebab5e88f6707c3789f5

if not exist "%filename%" (
  echo File %filename% does not exists.
  echo Download it first.
  goto :EOF
rem ) ELSE (
rem   echo file %filename%: FOUND
)


SETLOCAL ENABLEDELAYEDEXPANSION
SET count=1
FOR /F "tokens=* USEBACKQ" %%F IN (`CertUtil -hashfile %filename% SHA1`) DO (
  SET line!count!=%%F
  SET /a count=!count!+1
)

if "%sha1%"=="%line2%" (
  echo file %filename% chekcsum is CORRECT
) ELSE (
  echo file %filename% chekcsum FAILED
)
ENDLOCAL
