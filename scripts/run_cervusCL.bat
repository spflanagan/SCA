for %%f in (C:\Users\sflan\Documents\GitHub\SCA\results\parentage_haplotypes\*mat.crv) do (
	start "" /B "C:\Program Files (x86)\Field Genetics\Cervus\Cervus CL\CervusCL.exe" %%f /ALF /O"
) 

for %%f in (C:\Users\sflan\Documents\GitHub\SCA\results\parentage_haplotypes\*mat.crv) do (
	start "" /B "C:\Program Files (x86)\Field Genetics\Cervus\Cervus CL\CervusCL.exe" %%f /SIM /O"
) 

for %%f in (C:\Users\sflan\Documents\GitHub\SCA\results\parentage_haplotypes\*mat.crv) do (
	start "" /B "C:\Program Files (x86)\Field Genetics\Cervus\Cervus CL\CervusCL.exe" %%f /PAR /O"
) 

for %%f in (C:\Users\sflan\Documents\GitHub\SCA\results\parentage_haplotypes\*pat.crv) do (
	start "" /B "C:\Program Files (x86)\Field Genetics\Cervus\Cervus CL\CervusCL.exe" %%f /SIM /O"
) 

for %%f in (C:\Users\sflan\Documents\GitHub\SCA\results\parentage_haplotypes\*pat.crv) do (
	start "" /B "C:\Program Files (x86)\Field Genetics\Cervus\Cervus CL\CervusCL.exe" %%f /PAR /O"
) 
