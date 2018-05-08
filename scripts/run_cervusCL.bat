for %%f in (C:\Users\sflan\Documents\GitHub\SCA\results\parentsim\*mat.crv) do (
	start "" /B "C:\Program Files (x86)\Field Genetics\Cervus\Cervus CL\CervusCL.exe" %%f /ALF /O"
) 

for %%f in (C:\Users\sflan\Documents\GitHub\SCA\results\parentsim\*mat.crv) do (
	start "" /B "C:\Program Files (x86)\Field Genetics\Cervus\Cervus CL\CervusCL.exe" %%f /SIM /O"
) 

for %%f in (C:\Users\sflan\Documents\GitHub\SCA\results\parentsim\*mat.crv) do (
	start "" /B "C:\Program Files (x86)\Field Genetics\Cervus\Cervus CL\CervusCL.exe" %%f /PAR /O"
) 

for %%f in (C:\Users\sflan\Documents\GitHub\SCA\results\parentsim\*pat.crv) do (
	start "" /B "C:\Program Files (x86)\Field Genetics\Cervus\Cervus CL\CervusCL.exe" %%f /SIM /O"
) 

for %%f in (C:\Users\sflan\Documents\GitHub\SCA\results\parentsim\*pat.crv) do (
	start "" /B "C:\Program Files (x86)\Field Genetics\Cervus\Cervus CL\CervusCL.exe" %%f /PAR /O"
) 
