
call "%SLREALTIME_QNX_SP_ROOT%\%SLREALTIME_QNX_VERSION%\qnxsdp-env.bat"

cd .

chcp 1252

if "%1"=="" (make  -f Continues_motor_control_v9_tensiontest.mk all) else (make  -f Continues_motor_control_v9_tensiontest.mk %1)
@if errorlevel 1 goto error_exit

exit /B 0

:error_exit
echo The make command returned an error of %errorlevel%
exit /B 1