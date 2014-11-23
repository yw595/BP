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
c<<<<<<< local
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c=======
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c>>>>>>> other
       !VVVVVVVVVVVVVVVVVV   RANDOM VALUES VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV
       !!!! WARNING: Nvar*Nexpts*Nwvals must be less than 5Million
       subroutine random(s,SizeS)
       integer SizeS
       double precision s(SizeS)	! Declare the type of the rand() function
       integer icount                 	! Counts random numbers
       integer*4 timeArray(3)    	! Holds the hour, minute, and second

	   !--
	   ! In order to get a different sequence each time, we initialize the
       ! seed of the random number function with the sum of the current
	   !hour, minute, and second.
       !--
       call itime(timeArray)     ! Get the current time
       icount = rand( timeArray(1)+timeArray(2)+timeArray(3) )

	   !--
	   ! Calling rand() with an argument of zero generates the next number
	   ! in sequence.
	   !--
       do icount = 1,SizeS
          s(icount)=rand(0)
       enddo

       return
       end
       !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


       !VVVVVVVVVVVVVV   RANDOM VALUES WITH SEED ARGUMENT VVVVVVVVVVVVVVVVV
       subroutine random_withseed(s,SizeS,seed)
       integer SizeS
       integer seedlength
       integer, dimension(:), allocatable :: finalseed
       double precision s(SizeS)	!the returned array of random numbers
       integer seed	!has to be an integer
       integer icount

       call random_seed(size=seedlength)
       allocate(finalseed(seedlength))

       finalseed(:) = seed
       call random_seed(put=finalseed)	!this resets the seed for 'random_number()'

       call random_number(harvest=s)

       end subroutine random_withseed
       !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

       !VVVVVVVVVVVVVVVVVVVV  RANDOM PERMUTATION VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV
       subroutine randperm(num,myperm)
       !random ENOUGH, but could be better if we introduced a variable seed too

       integer  :: num
       integer myperm(num),timearray(8)
       integer :: mynumber, i, j, k
       real rand2(num),rand3(num)
       call date_and_time(values=timeArray)     ! Get the current time
       i = rand ( timeArray(8))

       do i=1,num
          rand2(i)=rand(0)
       enddo

       do i=1,num
          mynumber=1
          do j=1,num
             if (rand2(i) > rand2(j)) then
                mynumber = mynumber+1
             end if
          end do

          do k=1,i-1
             if (rand2(i) <= rand2(k) .and. rand2(i) >= rand2(k)) then
                mynumber = mynumber+1
             end if
          end do
          myperm(i) = mynumber
       end do

       return
       end
       !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

       !VVVVVVVVVVVVVVVVVVVV  NORMALIZE MATRIX BY ROW VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV
       !Normalizes the rows of the matrix
       !If the elements of the Rows are too small, there may be NaN Errors
       subroutine Normalize_Rows(Matrix,Nrows,Ncols)
       integer Nrows, row_idx
       integer Ncols, col_idx
       double precision Matrix(Nrows,Ncols)
       double precision MatCopy(Nrows,Ncols)
       MatCopy = Matrix
       do row_idx=1,Nrows
          Matrix(row_idx,:) = Matrix(row_idx,:)/
     :             sum(Matrix(row_idx,:))
       enddo

       !!-----Check For NaN Errors--------------------------------------
       do row_idx=1,Nrows
       do col_idx=1,Ncols
       if (Matrix(row_idx,col_idx).NE.Matrix(row_idx,col_idx)) then
          write(*,*) 'warning: NaN Error'
          write(*,*) 'Nrows = ',Nrows
          write(*,*) 'Ncols = ',Ncols
          write(*,*) MatCopy(row_idx,:)

          stop
       endif
       if (Matrix(row_idx,col_idx).LT.0) then
          write(*,*) 'Matrix contains negative numbers:'
          write(*,*) 'Fix in bo_calculations.f'
          stop
       endif
       enddo
       enddo
       !------------------------------------------------------------------

       end subroutine Normalize_Rows
       !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

       !VVVVVVVVVVVVVVVVVVVV  Calcualte Average Parameters VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV
       subroutine calculate_average_parameters(ProbDist,ParValues,
     :				  Matrix,dim1,dim2)
       integer dim1 	!dimension of Square Matrix
       integer dim2		!dimension of ParValues
       double precision ProbDist(dim1,dim2,dim1)
       double precision ParValues(dim2)
       double precision, intent(out) :: Matrix(dim1,dim1)
       integer RowIdx

       do RowIdx=1,dim1
          Matrix(RowIdx,:) = MATMUL(ParValues,transpose(
     :     		ProbDist(:,:,RowIdx)))
       enddo

       end subroutine calculate_average_parameters


       !VVVVVVVVVVVVVVVVVVVV  CONVERT INDICES TO SUBSCRIPTS VVVVVVVVVVVVVVVVVVVVVVVVVVV
       subroutine ind2sub(Rows,Cols,Indices,Subs,NumInd)
       integer Rows
       integer Cols
       integer NumInd
       integer Indices(NumInd)
       integer Subs(NumInd,2)	!By definition 2 subscripts per index (for a 2d matrix)
       integer j
       integer index


       do j=1,NumInd
          index = Indices(j)
          Subs(j,1) = mod(index,Rows)					!the row index
          Subs(j,2) = 1+ (index-(mod(index,Rows)))/Rows 	!column index
          if (Subs(j,1).eq.0) then
             Subs(j,1) = Rows	!mod=0 corresponds to last element
             Subs(j,2) = (index-(mod(index,Rows)))/Rows	!special case
          endif

       enddo

       end subroutine ind2sub
       !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

       !VVVVVVVVVVVVVVVVVVVV  CONVERT INDICES TO SUBSCRIPTS VVVVVVVVVVVVVVVVVVVVVVVVVVV
       subroutine ind2subprev(Rows,Cols,Indices,Subs,NumInd)
       integer Rows
       integer Cols
       integer NumInd
       integer Indices(NumInd)
       integer Subs(NumInd,2)	!By definition 2 subscripts per index (for a 2d matrix)
       integer j
       integer index


       do j=1,NumInd
          index = Indices(j)
          Subs(j,1) = mod(index,Cols)					!the row index
          Subs(j,2) = 1+ (index-(mod(index,Cols)))/Rows 	!column index
          if (Subs(j,1).eq.0) then
             Subs(j,1) = Rows	!mod=0 corresponds to last element
             Subs(j,2) = (index-(mod(index,Cols)))/Rows	!special case
          endif

       enddo

       end subroutine ind2subprev
       !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^





!--------------------------------------------------------------------------------------------------------------------------
!---------------------- 																	-------------------------------
!----------------------				SUBROUTINES FOR BELIEF PROPAGATION						-------------------------------
!----------------------																		-------------------------------
!--------------------------------------------------------------------------------------------------------------------------

       !(1)VVVVVVVVVVVVVV PREPARE PENALTY FOR TARGET NODE VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV
       subroutine prepare_penalty(TargIdx)
       use mod_CONSTANTS
       use mod_BP
       integer RegIdx		!local variable for the regulator index
       integer TargIdx		!local variable for the target index
       integer PriorIdx		!local variable for indexing through Prior Knowledge interactions
       integer LocsTargPK(Nprior)	!local variable for finding target nodes in PK

       !ASSIGN DEFAULT to all rows
       do RegIDX=1,Nnodes
          Penalty(RegIdx,:) = defaultPenalty(1,:)
       enddo

       !ASSIGN PRIOR KNOWLEDGE
       do PriorIdx=1,Nprior
          if (PK(PriorIdx,2).EQ.TargIdx) then
             RegIdx = PK(PriorIdx,1) 	!this is the regulator index

             if (PK(PriorIdx,3).GT.0) then
                Penalty(RegIdx,:) = positivePenalty(1,:)
             else
                Penalty(RegIdx,:) = negativePenalty(1,:)
        	 endif

    	  endif
       enddo

       Penalty = lambda*Penalty

       end subroutine prepare_penalty
       !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^




       !(2)VVVVVVVVVVVVVV  INITIALIZE MESSAGES FOR BELIEF PROPAGATION VVVVVVVVVVVVVVVVVVVVVVVV
       subroutine initialize_messages(useRandom)
       !Initialize messages into rho and rho_tm1
       use mod_CONSTANTS
       use mod_BP
       logical useRandom			!option for random (TRUE) or uniform (FALSE) initial messages for rho factors
       double precision InitialMessages(NumRandNumbers)	!array for storing random numbers
       integer NumMessagesNeeded
       integer RegIdx
       integer ExpIdx
       double precision tempMessage(1,Nwvals)		!temporary storage of messages for normalizing messages
       double precision NormConstant				!normalization constant for messages
       integer :: wvalIdx,mIdx		!indices for 3D matrix allocation

       !Total number of elements to be initialized (rho and rho_tm1)
       NumMessagesNeeded = (Nnodes*Nwvals*Nexpts)

       !ASSIGN
       if (useRandom) then
          !Retrieve seeded Random Numbers: seeded with hours,minutes,seconds of the current time
          call random(InitialMessages,NumRandNumbers)
          mIdx=1
          do RegIdx = 1,Nnodes
             do wvalIdx = 1,Nwvals
                do ExpIdx = 1,Nexpts
                   rho(RegIdx,wvalIdx,ExpIdx) = InitialMessages(mIdx)
                   mIdx=mIdx+1
                enddo
             enddo
          enddo

          !---- reshape() bug in mac os-10.7 and later ---
       	  !Put the random numbers in rho using |RESHAPE|
          !rho = reshape(InitialMessages(1:NumMessagesNeeded),
     	  !	:          (/ Nnodes,Nwvals,Nexpts /))

       elseif (.NOT.useRandom) then
          !Use Uniform Number for all messages
          rho = 1d0

          !---- reshape() bug in mac os-10.7 and later --
          !InitialMessages(1:NumMessagesNeeded) = 1d0
          !rho = reshape(InitialMessages(1:NumMessagesNeeded),
          !:          (/ Nnodes,Nwvals,Nexpts /))
       endif


       !NORMALIZE each row of each layer (sum(P) = 1)
       do ExpIdx = 1,Nexpts
          do RegIdx = 1,Nnodes
             tempMessage(1,:) = rho(RegIdx,:,ExpIdx)
             NormConstant = sum(tempMessage)
             rho(RegIdx,:,ExpIdx) = tempMessage(1,:)/NormConstant
          enddo
       enddo



       !Don't forget to initialize the reference messages for convergence criteria
       rho = rho+epsilon(1d0)	!Just to ensure there are no zero-value messages
       rho_tm1 = rho

       end subroutine initialize_messages
       !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

       !(3)VVVVVVVVVVVVVV  CALCULATE UPDATE MESSAGE VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV
       subroutine calculate_update(option,new_message)
       use mod_CONSTANTS
       use mod_BP
       integer option
       double precision arctanh_xi_mu 		!arctanh(xi_mu) for message calculation
       double precision xkWik(Nwvals)		!x_k_mu * Omega
       double precision new_message(Nwvals) !rho_update
       integer minIdx
       double precision minValue
       integer r

       !Execute the Proper Update Calculation
       if (option.EQ.0) then
       	  !Sigmoid = tanh()
       	  !Uses analytical integral: available in XXXX

          arctanh_xi_mu = atanh(xi_mu*Gamma(target_node))
          xkWik = XDATA(reg,mu)*Omega	!for all omega in Omega
          sol_term1 = xkWik+S_ik-arctanh_xi_mu+ui_mu
          sol_term1 = sol_term1**2

          sol_term2 = 2d0*beta*Svar_ik+1d0
          new_message = -beta*(sol_term1/sol_term2)
          minIdx=1
          do r=1,Nwvals
             if(new_message(r).LE.new_message(minIdx)) then
                minIdx=r
             endif
          enddo
          minValue=new_message(minIdx)
          do r=1,Nwvals
             new_message(r)=new_message(r)-minValue
          enddo
          do r=1,Nwvals-1
             !write(*,'(F6.2$)') new_message(r)
          enddo
          !write(*,'(F6.2)') new_message(r)
          new_message = exp(new_message)

          if (sum(new_message).eq.0) then
             write(*,*) new_message
             write(*,*) S_ik
             write(*,*) 'reg = ',reg
             write(*,*) 'FxdRegs=',FxdRegs
             write(*,*) 'target =',target_node
          endif
          call Normalize_Rows(new_message,1,Nwvals)

       endif

       end subroutine calculate_update
       !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

       !VVVVVVVVVVVVVVVVV  CREATE AUTO-MARGINAL FOR UNOBSERVED VVVVVVVVVVVV
       subroutine auto_marginal(TargetIdx)
       use mod_CONSTANTS
       use mod_BP
       integer TargetIdx
       integer RegIdx

       !The default marginals are P(wij=0) = 1
       !This function only needs to replace any constrained by prior
       ! knowledge.
       end subroutine auto_marginal
       !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^



       !VVVVVVVVVVVVVVVV   ADJUST_MESSAGES VVVVVVVVVVVVVVVVVVVVVVVVVVVVV
       !Takes fixed interactions and puts their individual messages to
       ! agree with those values.
       subroutine Adjust_Messages(FxdRegs,FxdVals,NumFxd)
       use mod_CONSTANTS
       use mod_BP
       integer NumFxd						!local
       integer FxdRegs(NumFxd)				!local
       double precision FxdVals(NumFxd)		!local

       integer j   							!local indexing variable
       integer k							!local indexing variable

       do j=1,NumFxd
          reg = FxdRegs(j)
          do k=1,Nexpts
             rho(reg,1:Nwvals,k) = epsilon(1d0)
             rho(reg,minloc(abs(FxdVals(j)-Omega)),k)
     :             = 1d0
          enddo
       enddo

       end subroutine Adjust_Messages
       !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


!--------------------------------------------------------------------------------------------------------------------------
!---------------------- 																	-------------------------------
!----------------------				SUBROUTINES FOR DECIMATION								-------------------------------
!----------------------																		-------------------------------
!--------------------------------------------------------------------------------------------------------------------------
       !VVVVVVVVVVVVVVVV  SET THE ORDER OF PARAMETERS TO FIX IN DECIMATION VVVVVVV
       subroutine set_param_fix_order(FinalOrder)
       use mod_CONSTANTS
       use mod_DECIMATION

       integer FinalOrder(NumToFix,2)

       call randperm(NumToFix,InteractionOrder)
       call ind2sub(Nnodes,Nnodes,InteractionOrder,
     : 		  FinalOrder,NumToFix)


       end subroutine set_param_fix_order
       !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^



       !VVVVVVVVVVVVVVVV  ASSIGN PARAMETER VALUE ACCORDING TO MARGINAL VVVVVVVVVVVVVV
       subroutine assign_parvalue(ProbDist,DistLength,
     :		MCprob,idx)
       use mod_DECIMATION		!I DON"T THINK WE NEED THIS MODULE HERE
       integer DistLength
       double precision ProbDist(DistLength)
       integer idx
       double precision MCprob
       double precision CumSum
       logical Found


       Found = .FALSE.
       CumSum = 0d0
       idx=1
       do while (Found.eqv..FALSE.)
          CumSum = CumSum + ProbDist(idx)
          if (MCprob.le.CumSum) then
             Found = .TRUE.
          else
       	     idx=idx+1
       	  endif
       enddo

       end subroutine assign_parvalue
       !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


       !VVVVVVVVVVVVVVVVVVVVVVVV Calculate Marginals from Messages      VVVVVVVVVVVVV
       !
       ! Built to handle large M (M>100) to avoid NaN errors
       !
       subroutine calcMarg(rho,N,Nwvals,M,Marginals,Penalty)
       integer N,Nwvals,M
       double precision, dimension(N,Nwvals,M) :: rho
       double precision, dimension(N,Nwvals) :: Marginals,Penalty


       double precision, dimension(N,Nwvals) :: TempMarg
       double precision, dimension(N,Nwvals) :: TotalField
       integer :: interval,numIntervals,remainder
       integer :: j,startIdx,finishIdx

       interval = 100
       numIntervals = ceiling(float(M)/float(interval))
       remainder = modulo(M,interval)

       !average numIntervals
       Marginals = 1d0	!initialize to 1
       do j=1,numIntervals
          startIdx = (j-1)*interval + 1
          if (j.EQ.numIntervals) then
             finishIdx = startIdx + remainder -1
          else
          	 finishIdx = j*interval
          endif


          !calculate contribution to marginal
          TotalField = sum(log(rho(:,:,startIdx:finishIdx)),3)	!sum along the dimension of M
          TempMarg = exp(TotalField)
          call Normalize_Rows(TempMarg,N,Nwvals)

          !factor into marginal
          Marginals = Marginals * TempMarg
          call Normalize_Rows(Marginals,N,Nwvals)


       enddo

       !factor in the penalty (penalty is positive, but exponent is negative)
       Marginals = Marginals * exp(-1d0*Penalty)
       call Normalize_Rows(Marginals,N,Nwvals)


       end subroutine calcMarg
       !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

!--------------------------------------------------------------------------------------------------------------------------
!---------------------- 																	-------------------------------
!----------------------				SUBROUTINES FOR SIMULATED ANNEALING						-------------------------------
!----------------------																		-------------------------------
!--------------------------------------------------------------------------------------------------------------------------
      subroutine doSA(Params,FixedSet,perc_idx,FxdOrder,Nunf) !VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV
      use mod_CONSTANTS
      double precision Params(Nnodes)
      integer FixedSet(Nnodes)
      integer perc_idx	!perceptron index (equivalent to target node)
      integer FxdOrder(Nnodes)
      integer Nunf

      !*** WARNING
      !Hardcoded to be Simulated Annealing over 5 Unfixed Values
      !************


      !SA variable declarations
      integer param_idx
      integer ToFix(Nunf)		!Nunf is the number of unfixed parameters
      double precision ToFixVals(Nunf)
      integer temp_idx		!temporary indexing variable
      integer mu			!index over the number of experimental conditions
      integer MCdone		!integer variable holding the state of the MC procedure
      integer rpermToFix(Nunf)	!integer array of randomly ordered parameters to fix
      double precision RandProbs(Nunf)		!random probabilities
      integer parVal_idx	!index for parameter value assignment
      double precision S(Nexpts) 	!Contribution of the Already Fixed parameters
      double precision S2(Nexpts)	!Contribution of ToFix Parameter Assignments
      double precision Error		!Total SoSq error for the ToFix Parameter Assignments
      double precision LastError	!Last Observed Error of the ToFix Parameter Assignments
      double precision LowestError	!Lowest Observed Error of the ToFix Parameter Assignments
      double precision BestSet(Nunf)		!Best accepted set of parameter assignments
      double precision ProbAccept	!Probability of accepting a non-lowest assignment
      double precision RandAccept(500000)	!Random probability for acceptance criteria
      integer NumRandCalls			!Number of times to evaluate acceptance criteria
      double precision Temperature	!Temperature in Simulated Annealing Sim
      integer iterations			!number of assignments tried
      integer temp_idx2				!indexing the FxdOrder
      !Prepare random values
      call random(RandAccept,500000)

      !identify the different variables
      temp_idx=1
      temp_idx2 = Nnodes-(Nunf-1)
      do param_idx=1,Nnodes
         if (FixedSet(param_idx).eq.0) then

            ToFix(temp_idx) = param_idx
            FxdOrder(temp_idx2) = param_idx
            temp_idx=temp_idx+1
            temp_idx2 = temp_idx2+1
         endif
      enddo


      !Calculate the constant contribution to the error of the fixed Parameters
      do mu=1,Nexpts
         S(mu) = sum(Params*XDATA(:,mu))+UDATA(perc_idx,mu)
      enddo

      !Start MC-Simulated Annealing
      LowestError  = 100000		!arbitrarily large initial value
      Temperature  = 20			!arbitrary starting point
      NumRandCalls = 1			!initialize the number of calls to 1
      iterations = 1			!start with the first iteration
      MCdone = 0
      do while (MCdone.eq.0)
         call randperm(Nunf,rpermToFix)
         call random(RandProbs,Nunf)

         !Fix configuration of ToFix parameters
         do temp_idx=1,Nunf
            param_idx = ToFix(rpermToFix(temp_idx))

            !write(*,*) 'param idx = ',param_idx
            !write(*,*) 'ProbDist  = ',
!     :       	FullNetMarg(param_idx,:,perc_idx)
            call assign_parvalue(FullNetMarg(param_idx,:,
     :		perc_idx),Nwvals,RandProbs(temp_idx),parVal_idx)

             ToFixVals(rpermToFix(temp_idx)) = Omega(parVal_idx)

          enddo
          !write(*,*) 'Newly Fixed Values: ',ToFixVals



         !Evaluate the Error of this new coniguration
         Error=0d0
         do mu=1,Nexpts
            S2(mu) = sum(ToFixVals(1:Nunf)*XDATA(ToFix,mu))
            Error = (Gamma(perc_idx)*XDATA(perc_idx,mu)
     :			- tanh(S(mu)+S2(mu)))**2 + Error
         enddo

        !accept if lowest error seen
         if (Error.lt.LowestError) then
            BestSet = ToFixVals
            LowestError = Error
         else !accept with probability as a function of error
            ProbAccept = exp(-Error/Temperature)
            if (RandAccept(NumRandCalls).le.ProbAccept) then
               BestSet = ToFixVals
               LowestError = Error
            endif
            NumRandCalls=NumRandCalls+1
            if (NumRandCalls.gt.500000) then
            	!get new random numbers
               call random(RandAccept,500000)
            endif
         endif

         !calculate termination criterion and alter temperature
         if (LowestError.lt.0.5) then
            MCdone=1
         endif
    	 if (Temperature.lt.1) then
    	    MCdone = 1
    	 endif
    	 if (mod(iterations,100).eq.0) then
    	    Temperature = .9*Temperature
    	 endif

        iterations = iterations+1
      enddo !end of while loop

      Params(ToFix) = BestSet
      FixedSet(:) = 1


      end subroutine doSA	!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^



