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
       subroutine BeliefPropagation()
       use mod_CONSTANTS
       use mod_BP
       
       !Fixed Information Variables: Fix all diagonal elements
       integer SelfFB(1)
       double precision SelfFB_val(1)
       
       
       !!!Calculate Marginals for all observed Targets
       do target_node=1,Nnodes
          if (ProtObs(target_node).EQ.1) then !only continue if the node is observed
          	SelfFB(1) = target_node
          	SelfFB_val(1) = 0d0
          	call doBP(SelfFB,SelfFB_val,1) !Calculate Marginals for this Target with BP
          	FullNetMarg(:,:,target_node) = Marginals
          	write(*,*) 'Finished Target: ',ProtName(target_node)
          else
             call prepare_penalty(target_node)
             Marginals = lambda-Penalty
             call Normalize_Rows(Marginals,Nnodes,Nwvals)
             FullNetMarg(:,:,target_node) = Marginals
          endif
       enddo   

       
       end subroutine BeliefPropagation
       !--------------------------------------------------------------------------------------------------------
       
       
       !VVVVVVVVVVVVVVVVVVVVVV  DO THE BP ALGORITHM VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV
       subroutine doBP(FxdRegs,FxdVals,NumFxd)
       use mod_CONSTANTS
       use mod_BP
       
       !Variables derived from Decimation -------------------------------------
       integer NumFxd					!Number of fixed interactions
       integer FxdRegs(NumFxd)			!The regulators with fixed interactions
       double precision FxdVals(NumFxd)	!The strength of the fixed interactions
       !-----------------------------------------------------------------------
       
       
       
       ! (1) Prepare the Penalty Matrix
       !Accommodate Sparsity/Prior Knowledge
       call prepare_penalty(target_node)
                
       !(2) Initialize rho and rho_tm1 
       call initialize_messages(randomIC)
       if (NumFxd>0) then
          call Adjust_Messages(FxdRegs,FxdVals,NumFxd)
       endif
       
       
       !(3) Calculate Normalized Marginals from Messages in Rho & Penalty
       call calcMarg(rho,Nnodes,Nwvals,Nexpts,Marginals,Penalty)
       
       
       !(4) Start BP iteration
       Converged = .FALSE.		!convergence criterion
       iteration_num=1
       
             	
       !(5) Start BP
       do while (.not. Converged) 	!Start BP
          !(5a) Loop Through Factors (Experiments)
          call randperm(Nexpts,exp_permu)	!permute the experiments (factor nodes)
          do factor_idx=1,Nexpts
             mu = exp_permu(factor_idx)
             
             !Define Constants for Calculations involving this Factor Node
             xi_mu = XDATA(target_node,mu)	
             ui_mu = UDATA(target_node,mu)
             	      
             !Calculate Constants for Calculations involving this Factor Node
             TotalField = log(Marginals)
             TFcorrection = TotalField - log(rho(:,:,mu))
             PWij_mu = exp(TFcorrection)
             call Normalize_Rows(PWij_mu,Nnodes,Nwvals)
             
             mean_Wij = MATMUL(Omega,transpose(PWij_mu))
             mean_Wij2 = MATMUL((Omega**2),transpose(PWij_mu))
             S = sum(mean_Wij*XDATA(:,mu))+UDATA(target_node,mu)
             Svar = sum((mean_Wij2-(mean_Wij**2))*
     :		       XDATA(:,mu)**2)        	      
            
                 		
     		!Apply correction to 'S' and 'Svar' 
     		!-- remove incorrect influence from the FIXED interactions
     		!-- replace with correct influence from FIXED interactions
     		minusTerm = sum(mean_Wij(FxdRegs)*XDATA(FxdRegs,mu))
     		S = S- minusTerm
     		plusTerm = sum(FxdVals*XDATA(FxdRegs,mu))
     		S = S+plusTerm
     		    		
     		Svar = Svar - sum((mean_Wij2(FxdRegs)-
     :			(mean_Wij(FxdRegs)**2))*XDATA(FxdRegs,mu)**2)
            
     		
     		       
            !(5b) Loop Through Problem Variables (Regulator Edges)
            call randperm(Nnodes,reg_permu) !permute regulators (variable nodes)
            do varnode_idx =1,Nnodes
               reg = reg_permu(varnode_idx)	!reg = 'k' index
               
               if (minval(abs(FxdRegs-reg)).ne.0) then !skips all fixed edges
                  !Subtract reg's contribution to S
             	  S_ik = S-(mean_Wij(reg)*XDATA(reg,mu))
             	  Svar_ik = Svar - (mean_Wij2(reg)-
     :        	               mean_Wij(reg)**2)*
     :						   XDATA(reg,mu)**2
     			  Sdev_ik  = sqrt(Svar_ik)
     			      			  
     			  
     			  ! calculate_update(option,result)
     			  ! option = 0: sigmoid = tanh()
     			  ! option = 1: sigmoid = simple sigmoid
     			  call calculate_update(0,rho_update)
     			  
     			  rho(reg,:,mu) = rho_update
     			endif	
             enddo !end loop over variable nodes
             	
             !(5c) Update the Total Field after all-variable nodes have been updated             
             TotalField = TFcorrection + log(rho(:,:,mu))
             Marginals = exp(TotalField)
             call Normalize_Rows(Marginals,Nnodes,Nwvals)
             
             
          enddo   !end loop over factor nodes
             
          
          
          !(5d) Check convergence Criteria.
          DeltaRho = abs(rho-rho_tm1)
          MaxDeltaRho = maxval(DeltaRho)
          if (MaxDeltaRho.lt.thresh) then
             
             Converged = .TRUE.
             
             if(prntFactors) then
               call writeFactorsToFile(SessionDirectory,rho,
     :          Nnodes,Nwvals,Nexpts,ProtName,target_node)
             endif
             
          else
             iteration_num = iteration_num+1
             rho_tm1 = rho
          endif   
             
       enddo	!Stop Iterations
       
       
       end subroutine doBP
       
       
       
       