cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c    <Inference of signaling networks in biological systems using Belief Propagation Algorithm>
c    Copyright (C) 2013  MSKCC, Authors: Chris Sander, Evan Molinelli,  Anil Korkut
c	Andrea Pagnani, Alfredo Braunstein
c
c
c    This program is free software: you can redistribute it and/or modify
c    it under the terms of the GNU General Public License as published by
c    the Free Software Foundation, either version 3 of the License, or
c    (at your option) any later version.
c
c    This program is distributed in the hope that it will be useful,
c    but WITHOUT ANY WARRANTY; without even the implied warranty of
c    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
c    GNU General Public License for more details.
c
c    You should have received a copy of the GNU General Public License
c    along with this program.  If not, see <http://www.gnu.org/licenses/>.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       module mod_CONSTANTS
       implicit none

       double precision, dimension(:,:), allocatable :: XDATA,UDATA
       double precision, dimension(:), allocatable :: Omega !array of weight values
       double precision, dimension(:), allocatable :: Gamma !array of alpha/beta ratios
       double precision MaxWvals,thresh
       double precision lambda,beta


       character(len=25), dimension(:), allocatable:: ProtName	!character array of protein (or node) names
       integer, dimension(:), allocatable :: ProtObs	!Binary identifiers of whether protein node is measured or not
       integer, dimension(:), allocatable :: ObsArray	!Indices of the observed nodes/proteins

       integer Nexpts	!Number of experiments (M)
       integer Nnodes	!Number of model nodes (N)
       integer Nwvals	!Number of possible parameter assignments
       integer Nprior	!Number of Prior Knowledge Interactions provided by user
       integer Nobs		!Number of the model nodes that are experimentally observed/measured

       integer, dimension(:,:), allocatable :: PK !PK(:,3) must be ~=0


       !Constants for making directories and files
       character(len=25) SessionID	!User-provided label for results folder
       character(len=4) year
       character(len=2) month
       character(len=2) day
       integer*4 today(3)					!integer array containing the date
       character(len=33) SessionDirectory	!Directory for results
       character(len=50) MakeDirectory_1	!Command to make a directory (first level)
       character(len=100) aSystemCall		!text for any system call
       character(len=200) ModelDir			!Directory for Models
       character(len=300) ModelFileLoc		!location for a single model file

       double precision, dimension(:,:,:), allocatable :: FullNetMarg
       double precision, dimension(:,:), allocatable :: defaultMarg
       double precision, dimension(:,:), allocatable :: W_avg	!Full Matrix of average weights
       double precision, dimension(:), allocatable :: Wi_avg	!Single Row of average weights
       double precision, dimension(:,:), allocatable :: aModelW				!Full Matrix of weights (from Decimation, for gradient descent)

       logical decidone 	!TRUE when all models are completed
       integer deci_model
       integer Ndeci		!Total Number of Decimated Models
       integer opt_maxSkip	!Max number of fixes without recalculating BP


       integer NumRandNumbers				!Integer of Number of random numbers to call (Max Size)
       parameter (NumRandNumbers = 5000000)	!Arbitrarily large number of random numbers

       !Optional Inputs
       character(len=32) :: arg		!necessary for getarg() function
       integer arg_idx				!index over all optional arguments
       logical prntFactors			!FALSE -> normal execution, TRUE -> print all factors at the end of the iteration.
       logical randomIC				!FALSE -> start with uniform messages, TRUE -> start with random messages (rho)

       integer, dimension(:), allocatable :: testInds
       integer, dimension(:,:), allocatable :: testSubs
       logical doTestInd2Sub

       end module mod_CONSTANTS
