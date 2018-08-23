To compile the program, command;
>> make

Then, to run the program, command orderly;
>> ./init N n M 
>> ./main

Here, 'N' is the number of stars in the cluster , 'n' and 'm' are the index parameters for the generalyzed polytropes to be used in distribution functions. First command will create "input.txt" which consisting of star parameters such as radius, radial velocity and tangential velocity. "main.cpp" will use that file to create cluster of stars.

After program returns, it'll create "radiis.txt" , "rminmax.txt" and "r_u_plane.txt" . Python files are using these files to plot.
*plot_radiis.py using --> "radiis.txt"
*plot_rminmax.py using --> "rminmax.txt"
*plot_ruplane.py using --> "r_u_plane.txt"

To use these, command;
>> phyton plot_radiis.py 
>> phyton plot_rminmax.py
>> phyton plot_ruplane.py

* "check.cpp" file includes control functions such as energy of the system and radii shells of the system.
* "init.cpp" file includes initialization routines for the generalyzed polytropes by using Simpson method for integration and Runge-Kutta for differentiation. General algorithm comes from Henon.
* "main.cpp" file includes step updating function and main loop that evolving the cluster itself.
