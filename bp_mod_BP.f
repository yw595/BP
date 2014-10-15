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
       module mod_BP
       ! BELIEF PROPAGATION VARIABLES
       ! 
       ! 'i' refers to the target node (the perceptron index)
       ! 'j' refers to non-cavity nodes
       ! 'k' refers to cavity nodes
       
          
             
       !ALLOCATABLE VARIABLES===========================================
       !BP method variables; the Total Field Approach
       double precision, dimension(:,:), allocatable :: Marginals
       double precision, dimension(:,:,:), allocatable :: rho
       double precision, dimension(:,:,:), allocatable :: rho_tm1
       double precision, dimension(:,:), allocatable :: Penalty
       double precision, dimension(:,:), allocatable :: TotalField
       integer, dimension(:), allocatable :: exp_permu
       integer, dimension(:), allocatable :: reg_permu
       double precision, dimension(:,:), allocatable :: PWij_mu
       double precision, dimension(:), allocatable :: mean_Wij
       double precision, dimension(:), allocatable :: mean_Wij2
       double precision, dimension(:), allocatable :: rho_update
       double precision, dimension(:,:,:), allocatable :: DeltaRho
       
       !Convenience variables for calculation
       double precision, dimension(:), allocatable :: sol_term1
       double precision, dimension(:,:), allocatable :: TFcorrection
       double precision, dimension(:,:), allocatable :: defaultPenalty
       double precision, dimension(:,:), allocatable :: positivePenalty
       double precision, dimension(:,:), allocatable :: negativePenalty
       
       !<DEBUG> Variables ********
       double precision, dimension(:,:), allocatable :: allRegrho 
       !================================================================
       
       !NON-ALLOCATABLE VARIABLES======================================================================================
       !BP option variables
       integer opt_uni_ic		!message initialization [=0:random intial messages][=1:uniform intitial messages]
       
       
       !BP method variables
       integer target_node		! 'i' the target node (perceptron index)
       integer target_idx		! MAYBE DELETE: index of the target node in the array of observed nodes
       integer mu				! the cavity factor
       integer factor_idx		! factor index
       integer varnode_idx		! variable node index
       integer reg				! the cavity variable
       integer reg_idx			! MAYBE DELETE: the index of the cavity variable in the random permutated array
       integer iteration_num	! current iteration number
       integer MaxIterations 	! the maximum allowable number of iterations to protect against non-convergence
       logical Converged
       double precision xi_mu	! the experimental expression of the target node in experiment mu
       double precision ui_mu	! the perturbation of the target node in experiment mu
       
       double precision S_ik	! Mean of the aggregate influence of the all-but-k nodes 
       double precision Svar_ik ! Varaince of the aggregate influence of the all-but-k nodes
       double precision Sdev_ik ! Standard deviation of the aggregate influence of the all-but-k nodes
       double precision MaxDeltaRho 	! Max of all the changes in messages for convergence
       
       !Convenience variables for indexing and calculations
       integer fxd_idx
       integer omega_idx
       integer pk_idx
       double precision S					! Mean Aggregate influence from ALL nodes
       double precision Svar				! Variance of the aggregate influence from ALL nodes
       double precision :: variance_correction
       double precision :: sol_term2
       double precision :: minusTerm		!the contribution to S that is replaced
       double precision :: plusTerm			!the replacement term
       !==============================================================================================================
       
       
       
       end module mod_BP