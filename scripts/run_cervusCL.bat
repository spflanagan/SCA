for %%f in (B:\ubuntushare\SCA\results\parentsim\*.crv) do (
	start "" /B "C:\Program Files (x86)\Field Genetics\Cervus\Cervus CL\CervusCL.exe" %%f /ALF /O"
) 

for %%f in (B:\ubuntushare\SCA\results\parentsim\*.crv) do (
	start "" /B "C:\Program Files (x86)\Field Genetics\Cervus\Cervus CL\CervusCL.exe" %%f /SIM /O"
) 

for %%f in (B:\ubuntushare\SCA\results\parentsim\*.crv) do (
	start "" /B "C:\Program Files (x86)\Field Genetics\Cervus\Cervus CL\CervusCL.exe" %%f /PAR /O"
) 
