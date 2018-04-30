for %%f in (B:\ubuntushare\SCA\results\parentage_biallelic\dradPruned*.crv) do (
	start "" /B "C:\Program Files (x86)\Field Genetics\Cervus\Cervus CL\CervusCL.exe" %%f /ALF /O"
	start "" /B "C:\Program Files (x86)\Field Genetics\Cervus\Cervus CL\CervusCL.exe" %%f /SIM /O"
	start "" /B "C:\Program Files (x86)\Field Genetics\Cervus\Cervus CL\CervusCL.exe" %%f /PAR /O"
) 
pause