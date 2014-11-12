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
      subroutine doDecimation()
      use mod_CONSTANTS
      use mod_BP
      use mod_DECIMATION
      integer j

      !(1) -----------  SETUP ------------------------------------------------------------------
      !Allocate and Prepare variables and random numbers
      call allocate_decimation()			!Allocate Memory Space for this Decimation
      call set_param_fix_order(FixOrder)	!Predetermine random order of interactions to fix
      call random(MCrandoms,NumToFix)		!declare monte-carlo random variables

      do l = 1,Nnodes
         do m = 1,Nnodes
            ModelW(l,m)=0
         enddo
      enddo
      !Set diagonal elements to zero in ModelW,FixedW,Tracker
      !Tracker is an array that holds the information for the order
      !          of fixed regulators (column) for each target (row)
      do fxd_target=1,Nnodes
         ModelW(fxd_target,fxd_target) = 0d0
         FixedW(fxd_target,fxd_target) = 1
         Tracker(fxd_target,1) = fxd_target
         if (ProtObs(fxd_target).eq.0) then
            ModelW(fxd_target,:) = 0d0
            FixedW(fxd_target,:) = 1

         endif
      enddo
      !-----------------------------------------------------------------------------------------

      integer deciOption
      deciOption = 1

      !(2) ============  LOOP THROUGH UNFIXED ELEMENTS =============
      do interaction=1,NumToFix

         fxd_target = FixOrder(interaction,1)
         fxd_reg    = FixOrder(interaction,2)

         if (FixedW(fxd_target,fxd_reg).EQ.0) then
            if(deciOption.EQ.1) then
               call assign_parvalue(DeciMarginals(fxd_reg,:,
     :		  fxd_target),Nwvals,MCrandoms(interaction),widx)
            endif
            if((deciOption.EQ.2).OR.(deciOption.EQ.3)) then
               widx=1
               do l=1,Nwvals
                  if(DeciMarginals(fxd_reg,l,fxd_target).GE.
     :            DeciMarginals(fxd_reg,widx,fxd_target)) then
                     widx=l
                  endif
               enddo
            endif
            ModelW(fxd_target,fxd_reg) = Omega(widx)
            FixedW(fxd_target,fxd_reg) = 1


            !prepare fixed-interaction data for BP
            NumFixed = sum(FixedW(fxd_target,:))
            Tracker(fxd_target,NumFixed) = fxd_reg
            target_node = fxd_target

            !recalculate marginals if observed
            consecutive_fixes(target_node) =
     :	   consecutive_fixes(target_node)+1
            if(consecutive_fixes(target_node).ge.opt_maxSkip) then
               if((deciOption.EQ.1).OR.(deciOption.EQ.2)) then
                  opt_uni_ic = 0
                  call doBP(Tracker(target_node,1:NumFixed),
     :            ModelW(target_node,Tracker(target_node,1:NumFixed)),
     :            NumFixed)
                  DeciMarginals(:,:,target_node) = Marginals
                  consecutive_fixes(target_node)=0
               endif
            endif

            !do SIMULATED ANNEALING on last remaining unfixed parameters
            !  This is because BP may break down for only a few variables.
            if (NumFixed.GE.(Nnodes-5)) then
               NumUnfixed = Nnodes-NumFixed
               call doSA(ModelW(fxd_target,:),FixedW(fxd_target,:),
     :          	fxd_target,Tracker(fxd_target,:),NumUnfixed)
            endif

         endif

      enddo


      !(3)---------   RECORD RESULTS TO FILE ------------
      write(ModelIdx,'(i4)') ModelsDone
      ModelName = 'Model_'//trim(adjustl(ModelIdx))
      call write_matrix(ModelW,Nnodes,Nnodes,
     :		SessionDirectory,ModelName)
      call write_mat2sif(ModelW,Nnodes,Nnodes,
     :		SessionDirectory,ModelName)

 	  call deallocate_decimation()
 	  !--------------------------------------------------





      end subroutine doDecimation
