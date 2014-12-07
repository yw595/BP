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
	   !VVVVVVVVVVVVVVVVVVVVVVVV  ALLOCATE CONSTANTS VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV
       subroutine allocate_CONSTANTS()
       use mod_CONSTANTS
       integer targ
       integer reg
       integer widx
       integer i
       integer j

       allocate(XDATA(Nnodes,Nexpts))	!All Data points
       allocate(UDATA(Nnodes,Nexpts))	!All perturbations
       allocate(ProtName(Nnodes)) 		!Node Names must have fewer than 25 characters
       allocate(ProtObs(Nnodes))		!Binary identifiers of whether protein node is measured or not
       allocate(ObsArray(Nobs))			!The indices of the observed nodes
       allocate(Gamma(Nnodes))			!The ratio of beta/alpha (necessary to avoid NaN errors)
       allocate(Omega(Nwvals))			!The allowable weight values
       allocate(FullNetMarg(Nnodes,Nwvals,Nnodes))  !The third dimension indexes the Target
       allocate(defaultMarg(1,Nwvals))	!Default Marginal for assignment
       allocate(W_avg(Nnodes,Nnodes))	!Average Network from calculated Marginals (FullNetMarg)
       allocate(Wi_avg(Nnodes)) 		!Single Row of average weight values (Marginal)
       allocate(aModelW(Nnodes,Nnodes)) !A single weight matrix

       allocate(testInds(35))
       allocate(testSubs(35,2))
       allocate(outputIdxs(numOutputIdxs,3))

       !assign values to Omega
       do widx =1,Nwvals
          Omega(widx) = 2*MaxWvals*(widx-0.5*(Nwvals+1))
     :  		/(Nwvals-1)
       enddo


       !Default Marginal (Prob(wij=0) = 1)
       !The default is a definite empty model.
       do widx=1,Nwvals
          if (Omega(widx).ne.0) then
             defaultMarg(1,widx) = 0d0
          else
             defaultMarg(1,widx) = 1d0
          endif
       enddo


       !Assign DefaultMarg to FullNetMarg
       do targ=1,Nnodes
          do reg=1,Nnodes
             FullNetMarg(reg,:,targ) = defaultMarg(1,:)
          enddo
       enddo

       do i=1,35
          testInds(i)=i
       enddo

       end subroutine allocate_constants
       !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

       !VVVVVVVVVVVVVVVVVVVVVVVV  DE-ALLOCATE CONSTANTS VVVVVVVVVVVVVVVVVVVVVVVVVVVVVV
       subroutine deallocate_CONSTANTS()
       use mod_CONSTANTS

       deallocate(XDATA)
       deallocate(UDATA)
       deallocate(ProtName)
       deallocate(ProtObs)
       deallocate(ObsArray)
       deallocate(Gamma)
       deallocate(Omega)
       deallocate(FullNetMarg)
       deallocate(defaultMarg)
       deallocate(W_avg)
       deallocate(Wi_avg)
       deallocate(aModelW)
       if(doTestInd2Sub .EQ. 1) then
       else
          deallocate(PK)
       endif

       deallocate(testInds)
       deallocate(testSubs)
       deallocate(outputIdxs)

       end subroutine deallocate_CONSTANTS
       !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


       !VVVVVVVVVVVVVVV  ALLOCATE BELIEF PROPAGATION VARIABLES VVVVVVVVVVVVVVVVVVVVVVV
       subroutine allocate_BP()
       use mod_CONSTANTS		!Constants from user
       use mod_BP				!Variables for each call to BP method
       integer wvalidx 	!local index of parameter values


       !These variables are only allocated once in the whole procedure instead of before each call to BP
       !Nnodes = Number of Nodes in the Final Model (observed+unobserved)
       !Nwvals = Number of parameter values in Omega
       !Nexpts = Number of experiments in the data set
       !
       ! 'j' refers to the non-cavity variable nodes
       ! 'k' refers to the cavity variable node

       !BP method specific variables
       allocate(Marginals(Nnodes,Nwvals))			!Single Marginal Matrix for one target
       allocate(rho(Nnodes,Nwvals,Nexpts))			!The factor-to-variable messages (current iteration)
       allocate(rho_tm1(Nnodes,Nwvals,Nexpts))		!The factor-to-variable messages (previous iteration)

       allocate(Penalty(Nnodes,Nwvals))				!Matrix of penalties to all possible assignments
       													!*used for sparsity and prior knowledge constraints
       allocate(TotalField(Nnodes,Nwvals))			!Sum of all variable-to-factor messages (log space)
       													!*minimum value of epsilon to avoid "NaN errors"
       allocate(exp_permu(Nexpts))					!Array of randomly permutated order for the experiments
       allocate(reg_permu(Nnodes))					!Array of randomly permutated order of variable nodes
       													!*Already fixed nodes will be skipped
       allocate(PWij_mu(Nnodes,Nwvals))				!Single set of messages from ALL nodes (cavity AND non-cavity) to ONE Factor
       													!*used for the all-but-mu calculations
       allocate(mean_Wij(Nnodes))					!Mean Wij values for all j (non cavity)
       													!*used for calculating mean of P(s)
       allocate(mean_Wij2(Nnodes))					!Mean (Wij)^2 values for all j (non cavity)
       													!*used for calculating the standard deviation fo P(s)
       allocate(rho_update(Nwvals))					!The cavity factor-to-variable update message
       allocate(DeltaRho(Nnodes,Nwvals,Nexpts))		!The difference in the factor-to-variable messages for evaluating convergence


       !Convenience variable for calculations
       allocate(sol_term1(Nwvals))					!One part of the integral calculation used for the message update
       allocate(TFcorrection(Nnodes,Nwvals))		!Change to the TotalField to incorporate message updates
       allocate(defaultPenalty(1,Nwvals))			!Default Penalty (penalizing any non-zero parameter)
       allocate(positivePenalty(1,Nwvals))			!Positive Penalty (penalizing any non-positive parameter)
       allocate(negativePenalty(1,Nwvals))			!Negative Penalty (penalizing any non-negative parameter)

       !<DEBUG> Variables*************
       allocate(allRegrho(Nexpts,Nwvals)) 			!looking at all messages for a single regulator
       !END <DEBUG> Variables *******


       !Assign Default Penalty, Positive Penalty and Negative Penalty
	   !this will make accommodating prior knowledge easier
       defaultPenalty(1,:) = 1d0
       positivePenalty(1,:) = 1d0
       negativePenalty(1,:) = 1d0
       do wvalidx = 1,Nwvals
          if (Omega(wvalidx).EQ.0) then
             defaultPenalty(1,wvalidx) = 0d0
          elseif (Omega(wvalidx).GT.0) then
             positivePenalty(1,wvalidx) = 0d0
          elseif (Omega(wvalidx).LT.0) then
             negativePenalty(1,wvalidx) = 0d0
          endif
       enddo




       end subroutine allocate_BP
       !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

       !VVVVVVVVVVVVVVV  DEALLOCATE BELIEF PROPAGATION VARIABLES VVVVVVVVVVVVVVVVVVVVVVV
       subroutine deallocate_BP()
       use mod_CONSTANTS		!Constants from user
       use mod_BP				!Variables for each call to BP method

       deallocate(Marginals)
       deallocate(rho)
       deallocate(rho_tm1)

       deallocate(Penalty)
       deallocate(TotalField)
       deallocate(exp_permu)
       deallocate(reg_permu)
       deallocate(PWij_mu)
       deallocate(mean_Wij)
       deallocate(mean_Wij2)
       deallocate(rho_update)
       deallocate(DeltaRho)


       !Convenience variable for calculations
       deallocate(sol_term1)
       deallocate(TFcorrection)
       deallocate(defaultPenalty)
       deallocate(positivePenalty)
       deallocate(negativePenalty)

       !<DEBUG> Variables*************
       deallocate(allRegrho) 			!looking at all messages for a single regulator
       !END <DEBUG> Variables *******

       end subroutine deallocate_BP
       !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^




       !VVVVVVVVVVVVVVV  ALLOCATE DECIMATION VARIABLES VVVVVVVVVVVVVVVVVVVVVVV
       subroutine allocate_decimation()
       use mod_CONSTANTS
       use mod_BP
       use mod_DECIMATION


       NumToFix = Nnodes**2

       allocate(DeciMarginals(Nnodes,Nwvals,Nnodes))
       allocate(ModelW(Nnodes,Nnodes))
       allocate(FixedW(Nnodes,Nnodes))
       allocate(FixOrder(NumToFix,2))
       allocate(MCrandoms(NumToFix))
       allocate(Tracker(Nnodes,Nnodes))
       allocate(InteractionOrder(NumToFix))
       allocate(consecutive_fixes(Nnodes))

       !Assign DeciMarginals to be equal to those fit during BP
       DeciMarginals = FullNetMarg

       !Initialize FixedW to zero
       FixedW = 0
       Tracker = 0
       FixOrder = 0
       consecutive_fixes = 0

       end subroutine allocate_decimation
       !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


       !VVVVVVVVVVVVVVV  DE-ALLOCATE DECIMATION VARIABLES VVVVVVVVVVVVVVVVVVVVVVV
       subroutine deallocate_decimation()
       use mod_CONSTANTS
       use mod_BP
       use mod_DECIMATION


       deallocate(DeciMarginals)
       deallocate(ModelW)
       deallocate(FixedW)
       deallocate(FixOrder)
       deallocate(MCrandoms)
       deallocate(Tracker)
       deallocate(InteractionOrder)
       deallocate(consecutive_fixes)

       end subroutine deallocate_decimation
       !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

