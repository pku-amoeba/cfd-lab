**************************************************************************
Project01 
for CFD Course, Prof.Tang Huazhong, Spring Semester, 2020, Peking Univ.

Author : Chen Junlin, Dept.Mechanics, PKU
Contact: junlin.chen@pku.edu.cn 
Date   : Apr 11, 2020

Summary: Applications of CIR/Upwind, Lax-Friedrichs, Lax-Wendroffs schemes
**************************************************************************

--------------------------------------------------------------------------
Files
--------------------------------------------------------------------------
README.txt - what you read now
UPWIND.m - upwind scheme for 1D convection equation
LAX_F.m - Lax-Friedrich scheme for 1D convection equation
LAX_W.m - Lax-Wendroff scheme for 1D convection equation
CASE#_.m - cases from Prof.Tang's lecture notes, which apply the numerical schemes referred above
PLOT.m - visualize the computational results
Project01_author name.pdf - reports of the results

--------------------------------------------------------------------------
How to run Project01
--------------------------------------------------------------------------
The program outputs numerical results in the form of ".dat", case by case, 
followed by the visualization of the data with the PLOT script provided or any other tools.
Example:
1. Open CASE1_wave_packet.m in MATLAB.
2. Manipulate necessary parameters in the first block "user paremeters",
   mostly, you need to change "ss" for alternative nemerical schemes,
   and change "T" to determine a specific time step for data print.
3. Open PLOT.m. Manipulate the settings according to your preference, if 
   your are familiar with the MATLAB visualization tools. Frequently, you
   have to change to file name(e.g. "wave_packet_upwind_2s.dat") of the 
   readily-inputed data, which can be found in the working directory.
4. Save the created figures to any directories, or using the auto-saved one 
5. Go back to 2. if you'd like to obtain results of other schemes or at  
   other time steps.
6*.Open UPWIND.m, LAX_F.m or LAX_W.m for editing if you want to develop the   
   numerical schemes. 

--------------------------------------------------------------------------
Many thanks for any advices for improvements of the program!
--------------------------------------------------------------------------


     






 


