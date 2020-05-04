**************************************************************************
Project02
for CFD Course, Prof.Tang Huazhong, Spring Semester, 2020, Peking Univ.

Author : Chen Junlin, Dept.Mechanics, PKU
Contact: junlin.chen@pku.edu.cn 
Date   : Apr 30, 2020

Summary: Applications of Lax-Friedrichs scheme and MacCormack' scheme to
         solve 1D Euler equations.
**************************************************************************

--------------------------------------------------------------------------
Files
--------------------------------------------------------------------------
README.txt - what you read now
Flux_LF.m - Lax-Friedrich scheme for flux calculation
Flux_LW.m - Two-step Lax-Wendroff (MacCormack) scheme for calculation
Tstep.m - Compute the local time step with the CFL condition
RiemannExact.m - Exact solution of Riemann's 1-D problem by Virginia Notaro
CASE#_.m - cases from Prof.Tang's lecture notes, which apply the numerical schemes mentioned above
Project02_author name.pdf - reports of the results

--------------------------------------------------------------------------
How to run Project02
--------------------------------------------------------------------------
The program outputs numerical results case by case, followed by visualizations of the data. 
Example:
1. Open CASE1_stdRP.m in MATLAB;
2. Manipulate necessary parameters in the first block "user paremeters",
   frequently, you need to change "scheme" for alternative nemerical     
   schemes, and change "Tmax" to determine a the specific time for output;
3. Save the created figures to any directories, or using the auto-saved      
   ones;
5. Go back to 2. if you'd like to obtain results of other schemes or at  
   other time.
6*.Open Flux_LF.m or Flux_LW.m for editing if you want to develop the   
   numerical schemes. 

--------------------------------------------------------------------------
Many thanks for any advices that help improve the program!
--------------------------------------------------------------------------


     






 


