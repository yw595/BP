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
      program doBP_ALL

      !This program takes input for a series of perturbation experiments and infers marginal
      !probability distributions for each of N^2 possible parameters for a neural-network
      !like model.  The model form is described in Nelander et al, 2008.
      !HANDLE INPUTS--------------------------------------------------------------------------------------------
      use mod_CONSTANTS			!module containing variable declarations for problem constants
      use mod_BP   				!module containing variable declarations for BP calculation variables
      use mod_DECIMATION		!module containing variable declarations for Decimation calculation variables


      call read_input(21,'TestData/input.txt')				!read the information from input.txt with unit=5 "from keyboard"
      call allocate_constants()         !Allocate Memory Space for MultiDimensional Constants
      									!Assigns values to Omega
      									!creates the default Marginal P(wij=0) = 1

      call read_exprdata('TestData/data.txt') 	!read in experimental data
      call read_pertdata('TestData/pert.txt')	!read in perturbation data
      call read_ProtName('TestData/name.txt')    !read and display included protein names
      call write_gamma('TestData/gamma.txt')     !calculate and write alpha/beta ratios
      call read_PK('TestData/prior.txt') 		!create the Prior Knowledge Matrix

      !optional arguments
      do arg_idx=1,iargc()
         call getarg(arg_idx,arg)

         select case(arg)
            case('-f','--Factors')	!print Factors (individual contributions to marginal distributions)
               prntFactors = .TRUE.
               write(*,*) '-> Will print out a file with all factors'

            case('-r','--Random') 	!use random initial messages for rho and rho_tm1
            	randomIC = .TRUE.
            	write(*,*) '-> Will use random initial '//
     :       	'messages for all factors'
         end select
      enddo
      !-----------------------------------------------------------------------------------------------------------

      !(1) CREATE NECESSARY DIRECTORIES-------------------------------
      !get date information for labeling
      call idate(today)
      write(day,'(i2)') today(1)
      write(month,'(i2)') today(2)
      write(year,'(i4)') today(3)
      if (today(1).lt.10) then
         day = '0'//trim(adjustl(day))
      endif
      if (today(2).lt.10) then
         month = '0'//trim(adjustl(month))
      endif
      SessionDirectory = year//month//day//'_'//trim(SessionID)
      MakeDirectory_1 = 'mkdir -p '//SessionDirectory
      call system(MakeDirectory_1)
      ModelDir = SessionDirectory
      call system('rm -f -r '//trim(adjustl(SessionDirectory))//'/*')	!clear all files in SessionDirectory
      !------------------------------------------------------------------


      !(2) CALL BELIEF PROPAGATION ALGORITHM --------------------------------------------
      call allocate_BP()	!Allocate Memory Space for Belief Propagation
      MaxIterations = 200 	!Maximum Number of Iterations Before Quitting
       						!*Should be user-provided

      opt_uni_ic = 1		!use uniform initial messages

      call BeliefPropagation()	!Calculate the Marginals all observed nodes
      call write_all_marginals(SessionDirectory) !write contents of FullNetMarg to files
      call calculate_average_parameters(FullNetMarg,Omega,W_avg,
     : 		 Nnodes,Nwvals)
      call write_matrix(W_avg,Nnodes,Nnodes,
     : 		SessionDirectory,'W_average')

      call write_mat2sif(W_avg,Nnodes,Nnodes,
     :		SessionDirectory,'W_average')
      !------------------------------------------------------------------------------------


      !(3) PERFORM DECIMATION ALGORITHM inside loop for all requested models --------------------
      do ModelsDone = 1,Ndeci		!This loop can be parallelized
         call doDecimation(ModelsDone) !perform and record decimation
         write(*,*) 'Finished Model #',ModelsDone
      enddo
      call deallocate_BP()
      !--------------------------------------------------------------------------------------------




      !(4) DO OPTIMIZATION: PINEDA OR OTHER -------------------------------------------------------
      do ModelsDone = 1,Ndeci
         write(ModelIDX,'(I4)') ModelsDone

         ModelFileLoc = trim(adjustl(ModelDir))//'/Model_'//
     :	   trim(adjustl(ModelIDX))//'.txt'

         !(--a) read in model matrix
         call read_dpmatrix_fromfile(Nnodes,Nnodes,aModelW,
     :	   trim(adjustl(ModelFileLoc)))

         !(--b) reduce gradient descent training patterns
         call prep_data(XDATA,UDATA,20,Nnodes,Nexpts,Nobs,ProtObs)

         !(--c) prepare gradient descent
         call prep_cg(Nnodes,20,1,aModelW,gamma,Nobs,ModelDir)

         !(--d) send to gradient descent with pineda()
         call pineda(ModelDir,ModelIDX)
      enddo
      !---------------------------------------------------------------------------------------------





      end ! end of main
