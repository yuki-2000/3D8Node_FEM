cd /d %~dp0
call G:\anaconda3\Scripts\activate.bat
call activate base

rem do not know whether this runs
rem G:\anaconda3\python.exe %1

rem run only one file
rem python %1


for %%i in (%*) do (
python %%i
)


pause