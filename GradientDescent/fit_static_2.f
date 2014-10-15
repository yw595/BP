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


C  [M_inferred,stats] = fit_static(z,p,omega,network_structure,options)
C
C  Zij     = experimental value  (node i, experiment j)
C  Pij     = perturbation matrix (node i, experiment j)
C  OMEGAij = observation matrix  (node i, experiment j)
C  STRUCTURErs = 1 iff node (s) may affect node (r)
C  OPTIONS is an object that contains
C          Sigmoid parameters (, , type)
C          Learning rate and momentum constant
C          Termination threshold and max iteration count
C          Initial guesses (not mandatory) for W, alpha, beta
C          Verbose and plot flags
C          Parameters used by structure learner
C
C  References:
C  Fernando J. Pineda: 
C  Generalization of Back-Propagation to Recurrent Neural Networks
C  Physical Review Letters 1987, VOL 59 NUMBER 19, PAGES 2229-22232

      subroutine fit_static_2(omega,topology,error)


	  !!******* EM: CHANGES ********************************************
	  ! (1) Also assumes that the observed nodes are listed before
	  !		the unobserved nodes.  I changed this reference the observed
	  !		matrix in data_o.txt inside the variable omega.  This 
	  !		change is made in the assignment of alpha and beta.
	  !		####WARNING####: This may affect things down the way.
	  !		The first nobs elements of ALPHA and BETA correspond to the 
	  !		observed alphas/betas.  Unobserved alphas/betas do not get
	  !		optimized.
	  !
	  !	(2) In call to Pineda ode, there is no check for testing steady
	  !		state.  Potentially simulating past steady state, but maybe
	  !		that case isn't a problem.  Potentially also NOT reaching
	  !		a steady state, which COULD be a problem.  Preliminary test
	  !		shows that we are reaching steady state.
	  !
	  !	(Q1)Make sure that z and x are consistently used. In this code
	  !		it looks like x stands for model nodes and z0 
	  !		stands for the convenience cariable and z for the measured
	  !		data.  Make changes to improve consistency
	  !
	  !	(3) There are some debugging variables also in sim.f
	  !
	  ! (4) Encode something to check for steady state check on both 
	  !		model simulations and pineda simulations.
	  !
	  !	(5) Make sure that we don't assume the fir n_obs nodes are 
	  !		observed. Instead check omega.  We allow alpha(i),beta(i) to
	  !		change even if (i) is not observed. Only the alpha(i) and
	  !		beta(i) with observed i get momentum change.  
	  !			*Perhaps this is something we should change.
	  !			*changes to parameter retrieval changed results, but for
	  !			the better.  If we allow all alpha and beta to change
	  !			the error got lower. 
	  !
	  !	(6) Removed ETAW: couldn't find a single use of it.
	  !
	  !****************************************************************


      use modmodel
      use modexp
      use modode
      use modoptions
      use modnonzero

      implicit none
      

      external ode_pineda_final, jac  
      
      
C     Sizes
      integer maxN,maxNexp,idx
      parameter (maxN=100)       ! max size of the network
      parameter (maxNexp=5000)       ! max number of pertubations  
      double precision topology(maxN,maxN),omega(maxN,maxNexp)
      
      integer nodes(maxN),n_nodes,n_interactions,PRIOR_PARAMETER,ETA
      double precision y(4*maxN+maxN*maxN+numnze),x0(maxN),z0(maxN)
      double precision w0(maxN*maxN),PRIOR(maxN*maxN)
      double precision yold(4*maxN+maxN*maxN+numnze)
      double precision omega2(maxN)
      double precision deltaError,lastError
      
      
      !**** EM*****
      double precision pinedaDiff 
      double precision tempDiff
      double precision simDiff
      !*******



C     Sim
      integer iexp
      double precision zpred(maxN)
      common/sim/iexp


C     Local
      double precision error
      integer t_end
      logical done
      integer icount, failcount
      integer i,ii,j,jc,k,ipert,iteration


      double precision CURRENT_ALPHA(maxN)
      double precision CURRENT_BETA(maxN)
      double precision CURRENT_W(maxN*maxN)

      double precision ALPHA_t_minus_1(maxN)
      double precision BETA_t_minus_1(maxN)
      double precision W_t_minus_1(maxN*maxN)
      
      double precision ALPHA_t_minus_2(maxN)
      double precision BETA_t_minus_2(maxN)
      double precision W_t_minus_2(maxN*maxN)


      double precision ALPHA_BP(maxN)
      double precision BETA_BP(maxN)
      double precision W_BP(maxN*maxN)


      double precision PREV_DELTA_ALPHA(maxN)
      double precision PREV_DELTA_BETA(maxN)
      double precision PREV_DELTA_W(maxN*maxN)



	  !ODE algorithm dlsode parameters
	  integer  iopt, iout, istate, itask, itol, liw, lrw
	  parameter (lrw = 22+ (4*maxN+maxN**2)*(5+1)+ 3*(4*maxN+maxN**2)
     $      + (4*maxN+maxN**2)**2 )
      parameter (liw = 20+(4*maxN+maxN**2))
      integer  mf,neq,iwork(liw),npts
      integer*4 time_Array(8), time_Array2(8)
      integer start_time,end_time
      double precision  atol, rtol, rwork(lrw),t,tout




	  !Parameters for termination
      call date_and_time(values=time_Array)
      start_time =60*24*30*time_array(2)+60*24*time_array(3)+ 
     :60*time_array(5) +time_array(6)
     


      t_end        = 4000                ! Integration right boundary.
      

C     1.5 Set initial values for the  model
      CURRENT_ALPHA(:) = 1.00
      CURRENT_BETA(:) = 1.00
      CURRENT_W(:) = 0d0
      k=0
      do j=1,modelN
         !if (omega(j,1).eq.1) then
            k=k+1
            CURRENT_ALPHA(k) = opt_ALPHAo(j)
            CURRENT_BETA(k) = opt_BETAo(j)
         !endif
      enddo
      
      do j=1,modelN             ! go column-wise
         do i=1,modelN
            CURRENT_W((j-1)*modelN+i)=opt_Wo(i,j)
         enddo
      enddo

      ALPHA_t_minus_1 = CURRENT_ALPHA
      BETA_t_minus_1  = CURRENT_BETA
      W_t_minus_1     = CURRENT_W
      
      ALPHA_t_minus_2 = CURRENT_ALPHA
      BETA_t_minus_2  = CURRENT_BETA
      W_t_minus_2     = CURRENT_W
      
      
C     1.6 Set initial values for the x and z variables
      do i=1,modelN
         x0(i)=0.d0 	!model variable    
         z0(i)=0.d0     !dummy variable in pineda (not be confused with z, the measured data) 
      enddo
      ode_ETA=opt_ETA
      
      ode_STRUCTURE=topology 	!think this is unused

      opt_iterations = opt_iterations
      
      
      
      
C     PART 2: MAIN LOOP 
      iteration=1
      done=.false.
      lastError = 100000	!Some Arbitrarily Large Error
      failcount=0
      do while (.not.done)      
         do ipert=1,Nexp
            
            iexp=ipert
            

C     2.1: RECURRENT BACKPROP SOLUTION FOR CURRENT EXAMPLE      
            do i=1,modelN
               ode_READOUTS(i)=omega(i,ipert)
            enddo
   
            !initialize ODE VECTOR
            do 40 i=1,modelN
               y(i)=CURRENT_ALPHA(i)			!first N elements are Alpha
               y(modelN+i)=CURRENT_BETA(i)		!second N elements are Beta
               y(modelN*2+i)=x0(i)				!third N elements are convenience variable 
               y(modelN*3+i)=z0(i)				!fourth N elements are model state variables
 40         continue
            do icount=1,numnze  
               i = nzrows(icount)
               j = nzcols(icount)
               idx = (j-1)*modelN+i
               y(4*modelN+icount) = CURRENT_W(idx)
             enddo           

      		!SIMUALTION STARTS for this experiment
            neq = 4*modelN+numnze
            t = 0.d0
            itol = 1
            rtol = 1.D-3
            atol = 1.D-8
            itask = 1
            istate = 1
            iopt = 1
            iwork(6) = 1000
            mf = 22
    
            do 70 iout = 1,t_end
               tout = iout
               pinedaDiff = 0d0
               tempDiff = 0d0
               do 65 ii=1,neq 
                  yold(ii)=y(ii)

                  !check steady-state convergence
                  tempDiff = abs(yold(ii)-y(ii))
                  if (tempDiff.gt.pinedaDiff) then
                     pinedaDiff = tempDiff
                  endif
                  
 65            continue

               call dlsode (ode_pineda_final, neq, y, t, tout, itol,
     $              rtol, atol,itask,istate, iopt, rwork, lrw, iwork,
     $              liw, jac, mf)

            
               if (istate.lt. 0) then !go to 900
                  write(*,*) 'old Eta = ',ode_ETA
                  write(*,*) 'changing eta and restarting iteration'
                  ode_ETA = .9*ode_ETA
                  write(*,*) 'new Eta = ',ode_ETA
                  go to 85
               endif

               if (done) then 
                  goto 90
               endif
 70         continue 
 			

 90         continue ! exit simulation, steady state reached
        
			!2.2.1 Retrieve Parameters
            do i=1,modelN
               if(i.le.nobs)then		!**** THIS MAY BE AN UNNECESSARY RESTRICTION ****
                  ALPHA_BP(i)=y(i)
                  BETA_BP(i)=y(modelN+i)
               else
               ALPHA_BP(i)=1
               BETA_BP(i)=1
               endif
            enddo
            do i=1,modelN*modelN
               W_BP(i)= 0.d0 !y(4*modelN+i)
            enddo
            do icount=1,numnze
               i = nzrows(icount)
               j = nzcols(icount)
               idx = (j-1)*modelN+i
               W_BP(idx) = y(4*modelN+icount)
            enddo   





C     2.2.2 previous update
            PREV_DELTA_ALPHA = ALPHA_t_minus_1 - ALPHA_t_minus_2
            PREV_DELTA_BETA  = BETA_t_minus_1  - BETA_t_minus_2
            PREV_DELTA_W     = W_t_minus_1     - W_t_minus_2
            
                   
             do i=1,modelN
                PREV_DELTA_ALPHA(i)=PREV_DELTA_ALPHA(i)*ode_READOUTS(i)
                PREV_DELTA_BETA(i)=PREV_DELTA_BETA(i)*ode_READOUTS(i)
               

			   !2.2.3 update parameters
      		   CURRENT_ALPHA(i)= ALPHA_BP(i)+ 
     :	 		   PREV_DELTA_ALPHA(i)*opt_momentum;
     		   CURRENT_BETA(i)= BETA_BP(i) + 
     :			   PREV_DELTA_BETA(i)*opt_momentum;
            enddo
            CURRENT_W   = W_BP  + PREV_DELTA_W*opt_momentum;
  
            
		    !2.2.4 book-keeping
            ALPHA_t_minus_2  = ALPHA_t_minus_1
            ALPHA_t_minus_1  = CURRENT_ALPHA
            BETA_t_minus_2   = BETA_t_minus_1
            BETA_t_minus_1   = CURRENT_BETA
            W_t_minus_2      = W_t_minus_1
            W_t_minus_1      = CURRENT_W

            write(*,'(A$)') '#'
            
            call date_and_time(values=time_Array2)
            end_time =60*24*30*time_array2(2)+60*24*time_array2(3)+ 
     :60*time_array2(5) + time_array2(6)
c           
            if((end_time-start_time).gt.30)then
               print*,"maximum time limit reached-terminated"
               done=.true.
               go to 999
            endif

         enddo                  ! end for all pert
         
	  
    
C     2.3 ERROR CALCULATION
    
    
C     2.3.1 Expansion. See comment 1.4. Here we'll expand the system with 4M+M*M
C     parameters to one with 4N+N*N parameters. If we didn't do this, then
C     we wouldn't be able to compare the error obtained for different
C     values of M.
         error=0.d0
         do j=1,Nexp
            ode_alpha=CURRENT_ALPHA
            ode_beta=CURRENT_BETA
            do jc=1,modelN       ! go column-wise
               do i=1,modelN
                  ode_WEIGHTS(i,jc)=CURRENT_W((jc-1)*modelN+i)
               enddo
            enddo
            ode_N_NODES=modelN
            do k=1,modelN
               ode_I(k)=p(k,j)  ! Perturbation bias 
            end do
            call forward_eval(zpred,istate)
            !write(*,*) (zpred(i),i=1,modelN) 
            if (istate.lt.0.d0) then
               write(*,*) 'istate',istate
            endif
            do i=1,modelN
               if (omega(i,j).gt.0.d0) then
                  error=error+omega(i,j)*(zpred(i)-z(i,j))**2
               endif
               !write(*,*) (z(i,j),i=1,modelN)
            enddo
         enddo
         
      	 !***** EVAN'S PRINT TO SCREEN           *********
      	 write(*,'(A11,F8.2)') 'SoSqError =',error
      	 write(*,'(A4)') 'done'
      	 !************************************************

		 !Update Current Values
         modelalpha=CURRENT_ALPHA		!TEMPORARY
         modelbeta=CURRENT_BETA			!TEMPORARY
         modelw=ode_WEIGHTS				!TEMPORARY

                
         
C     2.4 CHECK TERMINATION CRITERIA
 85      if (istate.lt.0) then
            if (failcount.gt.50) then
               done=.true.
            endif
            write(*,*) 'about to exit iteration loop...'
            istate=1
            write(*,*) 'failed iteration #',iteration
            iteration = iteration-1
            write(*,*) 'iteration-1 = ', iteration
            failcount=failcount+1



            
            write(*,*) 'fail count = ', failcount
            go to 980
         endif 
         deltaError = abs(error-lastError)/ode_ETA

         do ii=1,numnze
            i = nzrows(ii)
            j = nzcols(ii)
         enddo
         if (iteration.gt.opt_iterations) then
            done=.true.
         endif
         lastError = error
 980  iteration=iteration+1
      enddo
      
 
 999  continue
      
  
  
C Make a test prediction using the reverse engineered model  


      modelalpha=CURRENT_ALPHA
      modelbeta=CURRENT_BETA
      modelw=ode_WEIGHTS

         
      return
      
 910  format('dlsode, error halt. ISTATE =',I3)
      end

