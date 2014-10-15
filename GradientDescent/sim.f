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
	  !****** EM CHANGES ***********************************************
	  !	(1) Removed all 'debug' options as it was cluttering up the code in
	  !		forward_eval().
	  !
	  !	(2) Removed theta, theta2m, m, sigmoid_type, lambda, VM options.
	  !
	  !****************************************************************
	  
	  !SIMULATE MODEL -------------------------------------------------
      subroutine forward_eval(z,istate) 

      use modmodel
      external ode_pineda, jac  
      double precision z(modelN)

	  !ODE algorithm dlsode parameters
      integer  iopt, iout, istate, itask, itol, liw, lrw, maxN
      parameter (maxN=100)
      parameter (lrw = 22+ (4*maxN+maxN**2)*(5+1)+ 3*(4*maxN+maxN**2)
     $      + (4*maxN+maxN**2)**2 )
      parameter (liw = 20+(4*maxN+maxN**2))

      integer  mf,neq,iwork(liw),npts
      integer :: i,ii
      double precision  atol, rtol, rwork(lrw),t,tout, y(modelN),
     $     yold(modelN)
      
            
      logical done
      
	  !Check if arrays are sized appropriately
      if ((22+9*modelN+modelN**2).gt.lrw) then   ! check lrw
         write(*,*) 'forward-eval, lrw=',lrw,' required=',(22+9*modelN
     $        +modelN**2)
         stop
      endif
      if ((20+modelN).gt.liw) then   ! check liw
         write(*,*) 'forward-eval, liw=',liw,' required=',(20+modelN)
         stop
      endif

      
C     Initial values 
      do 30 i=1,modelN  
         y(i)=0.00
 30   continue
      
      neq = modelN
      t = 0.d0
      itol = 1
      rtol = 1.D-3
      atol = 1.D-8
      itask = 1
      istate = 1
      iopt = 0
      mf = 22
      

      npts=1000

      do 40 iout = 1,npts
         tout = iout
         tempDiff = 0d0
         simDiff = 0d0
         do 32 ii=1,neq 
            yold(ii)=y(ii)
 32      continue
 
         call dlsode (ode_pineda, neq, y, t, tout, itol, rtol, atol,
     $        itask,istate, iopt, rwork, lrw, iwork, liw, jac, mf)

         if (istate.lt.0)  then
            write(*,*) 'failure in forward_eval: stopping'
            stop
         endif
         
         !check if steady state
         done=.true.
         do 37 ii=1,neq 
            if (abs(yold(ii)-y(ii)).gt.atol/2) done=.false.
 37      continue
 
 		 
 
         if (done) then 
            goto 50
         endif
 40   continue !go to next time step
 50   continue ! exit simulation, steady state reached
      

      do 55 i=1,neq             ! save variable data
         z(i)=y(i)
 55   continue
      
      return
      end 

C ============================
C   Derivatives for steady state fitting
C ============================

      subroutine ode_pineda_final(neq,t,y,ydot)
      !***** EM CHANGES ********************************************
      !
      !	(1) removed all instances and dependencies on lambda
      !
      !*************************************************************
      
      use modmodel
      use modode
      use modexp
      use modoptions
      use modnonzero

      implicit none


C     ODE related
      integer neq
      double precision t,y(neq),ydot(neq)

C     Locals

      integer icount,iexp,i,j,ic,idx,idx2
      double precision wx_plus_p,x(modelN),zz(modelN)
      double precision ALPHA(modelN),BETA(modelN)
      double precision ww(modelN,modelN) 
      double precision deriv(modelN),JJ(modelN),sig_wx_eval(modelN)
      double precision sum(modelN),sum2(modelN)
      double precision V
      double precision CW(numnze)
      double precision nbr_pen
      


      common/sim/iexp

     
C =====================================================
C
C   ALPHA = y(1:N_NODES);
C   BETA  = y(N_NODES+1:2*N_NODES);
C   x     = y(2*N_NODES+1:3*N_NODES);
C   z     = y(3*N_NODES+1:4*N_NODES);
C   w = reshape(y(4*N_NODES+1:end),N_NODES,N_NODES);
C =====================================================

      do 10 i=1,modelN
         ALPHA(i)=y(i)
         BETA(i)=y(modelN+i)
         x(i)=y(2*modelN+i)
         zz(i)=y(3*modelN+i)
         sum(i) = 0.d0
         sum2(i) = 0.d0
 10   continue

      do 40 j=1,modelN
         do 30 i=1,modelN
            ww(i,j) = 0.d0
 30      continue
 40   continue
 
      do icount=1,numnze
         i = nzrows(icount)
         j = nzcols(icount)
         idx = (j-1)*modelN+i
         ww(i,j) = y(4*modelN+icount)
      enddo

      
      
      
C===============================================================
C 	  wx_eval     = w*x+I;         
C     sig_wx_eval = sigmoid(wx_eval);
C	  dxdt=-ALPHA.*x+BETA.*sig_wx_eval;
C===============================================================
      do 60 i=1,modelN
              
         wx_plus_p=p(i,iexp)                
         do 50 j=1,modelN
            if (ww(i,j).ne.0.d0) wx_plus_p=wx_plus_p+x(j)*ww(i,j)
 50      continue  
         
                  
         !for sigmoid = tanh()
         sig_wx_eval(i)   =  tanh(wx_plus_p)
         ydot(2*modelN+i) = -ALPHA(i)*x(i)+BETA(i)*sig_wx_eval(i)
         deriv(i)         =  (1-tanh(wx_plus_p)**2)
         
           
        
 60   continue
 
 

C==============================================================
C     equation 10 in Pineda
C J=zeros(N_NODES,1);
C J(READOUTS)=tau(READOUTS)-x(READOUTS);
C dzdt=-ALPHA.*z+J+((deriv.*z)'*w)';
C==============================================================

      do i=1,modelN
         if (ode_READOUTS(i).gt.0) then
            JJ(i)=z(i,iexp)-x(i)
         else
            JJ(i)=0.d0
         endif
      enddo
      

	  !calculate change to x ( model variable)
      do icount=1,numnze
         j = nzrows(icount)
         i = nzcols(icount)
         sum(i) = sum(i) +deriv(j)*ww(j,i)*zz(j)*BETA(j)
      enddo
      do i=1,modelN
         ydot(3*modelN+i) = -ALPHA(i)*zz(i)+sum(i)+JJ(i)
      enddo
          

C================================================================
C CALCULATE DERIVATIVE ON THE PARAMETERS Wij
C
C 
      do icount=1,numnze
         i = nzrows(icount)
         j = nzcols(icount)
         
         
         
         idx2 = 4*modelN
         
         nbr_pen = 1+20*abs((opt_Wo(i,j)-y(idx2+icount)))

         
         V = opt_L1*tanh(10*ww(j,i))
         ydot(idx2+icount) = ode_ETA*deriv(i)*zz(i)*x(j)*BETA(i)
         
         !Comment if you don't want to constrain changes to Wij
         !to be close to Wo (initial guess) values.
         ydot(idx2+icount) = ydot(idx2+icount) -ode_ETA*V
         ydot(idx2+icount) = ydot(idx2+icount)/nbr_pen
      enddo
C================================================================
      
      
     
     

C====================================================
C     KILLING WEIGHTS - SOME NEW CODE ADDED BY SVEN
								!f=find(STRUCTURE==0);
                                !dwdt(f)=-w(f)*0.1;
C====================================================

       do i=1,modelN 
          
          nbr_pen = 1+20*abs((opt_ALPHAo(i)-y(i)))
          ydot(i)=-ode_ETA*zz(i)*x(i)
          ydot(i) = ydot(i)/nbr_pen
          
          nbr_pen = 1+20*abs((opt_BETAo(i)-y(modelN+i)))
          ydot(modelN+i)=ode_ETA*zz(i)*sig_wx_eval(i)
          ydot(modelN+i) = ydot(modelN+i)/nbr_pen
       enddo


       end 
      




C     ODE model, called from dlsode to simulate the model
     
      subroutine ode_pineda(neq,t,y,ydot)
      
      !**** EM CHANGES ***********************************
      !
      !	(1) removed sigmoid type.  Assume it is tanh()
      !
      ! (2) removed ode_lambda.  It is not necessary or in any
      !		way a parameter of the model.
      !
      !****************************************************
      use modode
      use modmodel

C     ODE related
      integer neq
      double precision t,y(neq),ydot(neq)

C     Locals
      double precision decay_term,wy_plus_p
      double precision theta,m,wy2m
      integer iexp
      
    

      do 20 i=1,neq
         decay_term=-ode_alpha(i)*y(i)
         wy_plus_p=ode_I(i)        ! I = perturbation value
         do 10 j=1,neq
            if (ode_WEIGHTS(i,j).ne.0.d0) then
               wy_plus_p=wy_plus_p+y(j)*ode_WEIGHTS(i,j)
            endif
 10      continue  
                                
         !model equation for derivative   
        ydot(i)=decay_term+ode_beta(i)*tanh(wy_plus_p)
         
         
 20   continue

      end 



      
C     dummy subroutine for dlsode, not used
      subroutine jac (neq, t, y, ml, mu, pd, nrowpd)
      integer neq, ml, mu, nrowpd
      double precision t, y, pd
      dimension y(neq), pd(nrowpd,neq)
      return
      end
      