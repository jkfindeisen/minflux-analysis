:: Variables for code location (change accordingly)
@echo off
set remote-path=\\s1602-fs01.bpc.mpg.de\DATA16020
set disk-name=T:
:: Change this path to the location of the notebooks
set path-to-startup-folder=T:\Antonio_Politi\Code\minflux-analysis\python
:: Path to conda installation (change accordingly)
set conda-activate=C:\Users\apoliti\Miniconda3\Scripts\activate.bat
:: If code is locally installed remove this comment
if not exist %disk-name% (
	net use %disk-name% %remote-path%
)
cd /D %path-to-startup-folder%
dir 
pause
@echo on
call %conda-activate%
start jupyter notebook

