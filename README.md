If you use this code then please cite the following reference: 
Qi Xu, M.L. Sanyang, A.Kaban. Large Scale Continuous EDA Using Mutual Information. IEEE Congress on Evolutionary Computation 2016 (CEC-2016)


This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not, see <http://www.gnu.org/licenses/>.

"Experiments" folder contains files for all the experiments.
 
'batch_' files makes the experiments easy to run. it runs EDA-MCC-MI automatically by simply setting 4 parameters: (pD, pNp, pFun, pRuns), where

  pD    : Dimensionality
  pNp   : Population size
  pFun  : Function number
  pRuns : Number of runs.

The names of functions are displayed in benchmark.m.

EDA-MCC-MI stores the files for EDA-MCC-MI. Parameters: (D,Pop,FuncNum,runs). So this time you can control running times.

EDA-MCC-MI_THETA is used for analytical experiments. Parameters: (pD, pNp, pFun, pRuns, pTheta).


Result filename format:

EDA-MCC data: BFE_parameters
EDA-MCC-MI:   BFG_parameters
Theta/WEAK:   THETA/WEAK_parameters

=========================================================================

Plot Results: Trajectory.m shows an example on how to visualise the results. However, you may have to change the filenames or other parameters in order to show all kinds of results.

Note: you may also change the number of loops in Trajectory.m if the number of experiments is not 25.

=========================================================================

How to compile MI.cpp?

Compile_MI.m shows an example in windows system using Visual Studio. It makes a difference if you compile this file using some optimisation options, as well as multi-thread libs in C++11 standard library. If you use g++, you may have to google it, but it should not be hard to find.

