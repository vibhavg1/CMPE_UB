# CMPE_UB

Upper Bounding Methods for the constrained most probable explanation task.

Contents: Code and data for the paper:
    “Novel Upper Bounds for the Constrained Most Probable Explanation Task” 
    by Tahrima Rahman, Sara Rouhani and Vibhav Gogate, NeurIPS 2021 (To appear).


Requirements for Compilation:
* Access to CBC Mixed Integer Linear Programming (MILP) solver from Google. See: https://projects.coin-or.org/Cbc
* Google OR tools with C++ interface. See: https://developers.google.com/optimization
* To compile the code, you can use the provided CMakeLists.txt file as a reference
* The code will yield the following executable: CMPE_UB: Upper Bounding Algorithms for solving the CMPE task
	
Requirements for Running the Code:
* Once compiled, you can download the Markov network files from the UAI 2014 competition. See: http://www.hlt.utdallas.edu/~vgogate/uai14-competition/
* After compiling you can get help using: "<executable> -h". For example,

		./CMPE_UB -h
    


* The program will run the four bounding algorithms and output the results. On the standard output, for each algorithm, the program will print
	iteration#,time-in-seconds,current-bound

* Code is released under: MIT License.
