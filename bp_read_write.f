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
c<<<<<<< local
c
c=======
c>>>>>>> other
      !VVVVVVVVVVVVVVVVVV   READ THE INPUT FILE   VVVVVVVVVVVVVVVVVVVVVVVVVVVVVV
      subroutine read_input(unit,filename)
      !to read in inputs from command: call read_input(5,'')
      !to read in inputs from a file:  call read_input(21,filename)

      use mod_CONSTANTS

      integer unit
      character(*) filename	!only for unit testing

      !for unit testing
      if (unit.ne.5) then
         open(unit,file=trim(adjustl(filename)),status='OLD')
      endif

      read(unit,"(I5)") doTestInd2Sub
      read(unit,"(A50)") inputDir
      read(unit,*) SessionID		!String label for Results Folder
      read(unit,"(i5)") Nexpts 		!number of experiments
      read(unit,"(i5)") Nnodes 		!number of nodes in the data set
      read(unit,"(i5)") Nwvals		!number of possible weight values
      read(unit,"(F5.2)") MaxWvals 	!the maximum magnitude of an allowed interaction
      read(unit,"(E8.2)") thresh 	!BP convergence threshold
      read(unit,"(F5.2)") lambda 	!Penalty on Non-Zero (or Prior Knowledge) edges
      read(unit,"(F5.2)") beta  	!Inverse temperature
      read(unit,"(i5)") Nprior 		!Number of Prior Knowledge interactions
      read(unit,"(i5)") Nobs  	 	!Number of Observed Proteins
      read(unit,"(i5)") Ndeci		!Number of Decimated Models to create
      read(unit,"(i5)") opt_maxSkip	!Decimation option to skip recalculating BP after every fixed value

      !close file if unit testing
      if (unit.ne.5) then
         close(unit)
      endif

      end subroutine read_input
      !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

       !VVVVVVVVVVVVVVVVVV   READ EXPERIMENTAL DATA VVVVVVVVVVVVVVVVVVV
       !reads in experimental data into XDATA variable
       subroutine read_exprdata(filename)
       use mod_CONSTANTS, only : XDATA,Nnodes,Nexpts
       character(*) filename
       integer node_idx,expr_idx

       open(21,file=trim(adjustl(filename)),status='OLD')
       do node_idx=1,Nnodes
          read(21,*) (XDATA(node_idx,expr_idx),expr_idx=1,Nexpts)
       enddo
       close(21)

       XDATA = XDATA + epsilon(1d0)

       end subroutine read_exprdata
       !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


       !VVVVVVVVVVVVVVVVVV   READ PERTURBATION DATA VVVVVVVVVVVVVVVVVVV
       !reads in perturbation data into UDATA variable
       subroutine read_pertdata(filename)
       use mod_constants
       character(*) filename
       integer node_idx,expr_idx

       open(21,file=trim(adjustl(filename)),status='OLD')
       do node_idx=1,Nnodes
          read(21,*) (UDATA(node_idx,expr_idx),expr_idx=1,Nexpts)
       enddo
       close(21)

       end subroutine read_pertdata
       !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^



       !VVVVVVVVVVVVVVVVVV   READ PROTEIN NAMES VVVVVVVVVVVVVVVVVVV
       !reads and displays protein names
       !reads and displays whether protein is observed
       !creates array of indices of the observed nodes
       !calculates the number of total observed nodes
       !Fixes marginals for unobserved nodes to default of P(W=0) = 1
       subroutine read_ProtName(filename)
       use mod_constants
       character(*) filename
       integer node_idx,obs_idx
       integer OmegaEq0(1)
       character*25 tempNodeName

       OmegaEq0 = minloc(abs(Omega))

	   obs_idx=1
       open(21,file=trim(adjustl(filename)),status='OLD')
       do node_idx=1,Nnodes
          read(21,"(i2,A25)") ProtObs(node_idx),ProtName(node_idx)
          if (ProtObs(node_idx).eq.1) then
             ObsArray(obs_idx) = node_idx
             obs_idx = obs_idx+1
          endif
       enddo
       close(21)

       write(*,*) 'The Observed Proteins Are-----------------------'
       do obs_idx=1,Nobs
          write(*,*) ProtName(ObsArray(obs_idx))
       enddo
       write(*,*) '------------------------------------------------'
       write(*,*) 'The Unobserved Proteins Are --------------------'
       do node_idx=1,Nnodes
          if (ProtObs(node_idx).ne.1) then
             write(*,*) ProtName(node_idx)
             !-- Fix these marginals to 0or1 for completeness
             FullNetMarg(:,:,node_idx) = 0d0
             FullNetMarg(:,OmegaEq0,node_idx) = 1d0
          endif
       enddo
       write(*,*) '-------------------------------------------------'




       end subroutine read_ProtName
       !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

       !VVVVVVVVVVVVVVVVVV   WRITE GAMMA RATIOS VVVVVVVVVVVVVVVVVVV
       !calculates and writes gamma ratios (alpha/beta)
       subroutine write_gamma(filename)
       use mod_constants
       character(*) filename
       integer node_idx
       double precision Xdiabs(Nnodes,Nexpts)

       Xdiabs = abs(XDATA)
       open(21,file=trim(adjustl(filename)))
       do node_idx =1,Nnodes
          Gamma(node_idx) = 0.91/maxval(Xdiabs(node_idx,1:Nexpts))
          write(21,'(F5.2)') Gamma(node_idx)
       enddo
       close(21)

       end subroutine write_gamma
       !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

       !VVVVVVVVVVVVVVVVVV   READ IN AND CREATE PRIOR KNOWLEDGE VVVVVVVVVVV
       subroutine read_PK(filename)
       use mod_constants
       integer prior_idx
       integer tempPrior(3)
       integer targ,reg,PKsign,reg_idx
       integer OmegaEq1(1)
       integer OmegaEq0(1)
       character(*) filename

       allocate(PK(Nprior,3))
       write(*,*) 'The PRIOR information IS --------------------'
       write(*,*) 'number of Prior Knowledge Edges =',Nprior
       write(*,*) '*********************************'

       open(21,file=trim(adjustl(filename)),status='old')
       do prior_idx=1,Nprior
          read(21,*) tempPrior
          reg  = tempPrior(1)
          targ = tempPrior(2)
          PKsign = tempPrior(3)

          PK(prior_idx,1) = reg 	!the regulator
          PK(prior_idx,2) = targ  	!the target
          PK(prior_idx,3) = PKsign  !the sign

          write(*,*) 'the target = ',ProtName(targ)
          write(*,*) 'the reg    = ',ProtName(reg)
          write(*,*) 'the sign   = ',PKsign
          write(*,*) '*********************************'

          ! (below) If we have Prior Knowledge on unobserved nodes, put them into the
          ! FullNetMarg variable, holding the final marginals.  This has
          ! no effect on BP results, but it may be useful for completeness.
          ! Used to assign Fixed Edges in Decimation allocation
          if (ProtObs(targ).ne.1) then
             OmegaEq1 = minloc(abs(Omega-PKsign))
             OmegaEq0 = minloc(abs(Omega))
             do reg_idx=1,Nnodes
                if (reg_idx.eq.reg) then
                   FullNetMarg(reg_idx,OmegaEq0,targ) = 0d0
                   FullNetMarg(reg_idx,OmegaEq1,targ) = 1d0
                endif
             enddo
          endif
          !End Handle unobserved targets with PK (obove)


       enddo
       write(*,*) '------------------------------------------------'

       end subroutine read_PK
       !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


       !VVVVVVVVVVVVVVVV  WRITE ALL MARGINALS VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV
       ! write to Marginals Subdirectory
       ! FileName = [index]_[Obs]_[NodeName]
       subroutine write_all_marginals(Directory)
       use mod_CONSTANTS
       character(len=100) MargFileName	!sufficiently large length
       character(len=3) TargIdx_string
       character(len=1) Obs_string
       character(*) Directory
       integer RegIdx
       integer ExpIdx
       integer TargIdx
       integer Widx
       integer FileId

       !Make 'Marginals' Directory if one does not exist
       aSystemCall = 'mkdir '//trim(adjustl(Directory))//
     :			'\Marginals'
       call system(aSystemCall)

       !writes the data contained in FullNetMarg
       do TargIdx=1,Nnodes
          write(TargIdx_string,'(i3)') TargIdx		!convert integer to string
          write(Obs_string,'(i1)') ProtObs(TargIdx)	!convert integer to string

          !FileName = '[Directory]/Marginals/[index]_[Obs]_[NodeName].txt
          MargFileName = trim(adjustl(Directory))//'\'//'Marginals\'
     :		    //trim(adjustl(TargIdx_string))//
     : 			'_'//Obs_string//'_'//
     :			trim(adjustl(ProtName(TargIdx)))//'.txt'

          open(21,file=MargFileName,action='write',status='replace')
          do RegIdx = 1,Nnodes
             do Widx =1,(Nwvals-1)
                write(21,'(F6.2$)') FullNetMarg(RegIdx,Widx,TargIdx)
             enddo
             write(21,'(F6.2)') FullNetMarg(RegIdx,Widx,TargIdx)
          enddo
          close(21)
       enddo
       end subroutine write_all_marginals
       !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


       !VVVVVVVVVVVVVVVV  WRITE A MATRIX TO FILE VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV
       subroutine write_matrix(Matrix,dim1,dim2,Directory,fname)
       integer dim1
       integer dim2
       character(len=*) Directory
       character(len=*) fname
       character(len=100) Full_fname
       double precision Matrix(dim1,dim2)
       integer dim1_idx
       integer dim2_idx
       parameter(cutoff=0.25)

       Full_fname = trim(adjustl(Directory))//'\'
     : 				//trim(adjustl(fname))//'.txt'

       open(21,file=trim(adjustl(Full_fname))
     :  	,action='write',status='replace')
       do dim1_idx=1,dim1
          do dim2_idx=1,(dim2-1)
             write(21,'(F6.2$)') Matrix(dim1_idx,dim2_idx)
          enddo
          write(21,'(F6.2)') Matrix(dim1_idx,dim2_idx)
       enddo
       close(21)

       end subroutine write_matrix
       !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


      !VVVVVVVVVVVVVVVVV  WRITE A MATRIX TO A SIF FILE VVVVVVVVVVVVVVVVVVVVVVVVV
      subroutine write_mat2sif(Matrix,dim1,dim2,Directory,fname)
      use mod_CONSTANTS
      integer dim1
      integer dim2
      double precision cutoff
      parameter(cutoff = 0.25)
      double precision Matrix(dim1,dim2)
      integer dim1_idx
      integer dim2_idx
      character(len=99) LineContent
      character(len=*) Directory
      character(len=*) fname
      character(len=60) Full_fname
      Full_fname = trim(adjustl(Directory))//'\'//
     :			trim(adjustl(fname))//'.sif'
      write(*,*) Full_fname


       open(21,file=trim(adjustl(Full_fname)),action
     :       ='write',status='replace')
       do dim1_idx =1,dim1
          do dim2_idx=1,dim2

             if (Matrix(dim1_idx,dim2_idx).gt.cutoff) then 	!activates
                LineContent = adjustl(trim(ProtName(dim2_idx)))//
     :           ' activates '//trim(adjustl(ProtName(dim1_idx)))
     			write(21,'(A)') adjustl(LineContent)
     		 elseif (Matrix(dim1_idx,dim2_idx).lt.(-cutoff)) then !inhibits
     			LineContent = adjustl(trim(ProtName(dim2_idx)))//
     :           ' inhibits '//trim(adjustl(ProtName(dim1_idx)))
     			write(21,'(A)') adjustl(LineContent)
             endif
          enddo
       enddo
       close(21)


       end subroutine write_mat2sif
       !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


       !VVVVVVVVVVVVVVVVVVV  WRITE FACTORS TO SINGLE FILE VVVVVVVVVVVVVVVVVVVVVVVV
       subroutine writeFactorsToFile(dir,Factors,Nregs,Nvals,Npats,
     :	ProtNames,targ_idx)
       character(*) dir				!directory
       integer :: targ_idx
       character(len=75) filepath	!full path of the factors.txt file
       character(len=25) :: ProtNames(Nregs)	!Protein Names
       integer :: Nregs,Nvals,Npats	!number of regulators,w-values,training-patterns
       double precision :: Factors(Nregs,Nvals,Npats)
       logical :: file_exists
       integer :: row_idx,col_idx,pat_idx	!for writing data to file
       character(len=100) :: factorLabel
       character(len=5) :: pat_label

       !check if file exists: file = [SessionDirectory]/Verbose/factors.txt'
       filepath = trim(adjustl(dir))//'\Verbose\factors.txt'
       inquire(file=trim(adjustl(filepath)),exist=file_exists)
       if(file_exists) then
          open(21,file=trim(adjustl(filepath)),STATUS='OLD',
     :      POSITION='APPEND')
       else
          call system('mkdir '//trim(adjustl(dir))//'\Verbose')
          open(21,file=trim(adjustl(filepath)),STATUS='NEW')
       endif

       !write factors to file
       do row_idx=1,Nregs
          do pat_idx =1,Npats
          	write(pat_label,'(i5)') pat_idx
          	factorLabel = trim(adjustl(ProtNames(targ_idx))) //
     :			'/' // trim(adjustl(ProtNames(row_idx))) // '/' //
     :			trim(adjustl(pat_label))

          	 write(21,'(A$)') trim(adjustl(factorLabel)) // ' '
             do col_idx=1,Nvals
                if (col_idx.lt.Nvals) then
                   write(21,'(F6.3$)') Factors(row_idx,col_idx,pat_idx)
                else
                   write(21,'(F6.3)') Factors(row_idx,col_idx,pat_idx)
                endif
             enddo
          enddo


       enddo


       close(21)


       end subroutine
       !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


       !VVVVVVVVVVVVVVVVV  READ A DOUBLE PRECISION MATRIX FROM FILE VVVVVVVVVVVVVV
       subroutine read_dpmatrix_fromfile(Nrows,Ncols,theMatrix,filename)
       !assumes format = F6.2
       integer :: Nros,Ncols
       integer :: ridx,cidx
       double precision theMatrix(Nrows,Ncols)
       character(*) filename


       open(21,file=trim(adjustl(filename)),STATUS='OLD')
       do ridx=1,Nrows
          do cidx=1,(Ncols-1)
             read(21,'(F6.2$)') theMatrix(ridx,cidx)
          enddo
          read(21,'(F6.2)') theMatrix(ridx,Ncols)
       enddo

       end subroutine read_dpmatrix_fromfile
       !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


       !VVVVVVVVVVVVVVVVV  READ A REAL PRECISION MATRIX FROM FILE VVVVVVVVVVVVVV
       subroutine read_realmatrix_fromfile(Nrows,Ncols,
     :  	theMatrix,filename)
       !assumes format = F6.2
       integer :: Nros,Ncols
       integer :: ridx,cidx
       real theMatrix(Nrows,Ncols)
       character(*) filename


       open(21,file=trim(adjustl(filename)),STATUS='OLD')
       do ridx=1,Nrows
          do cidx=1,(Ncols-1)
             read(21,'(F6.2$)') theMatrix(ridx,cidx)
          enddo
          read(21,'(F6.2)') theMatrix(ridx,Ncols)
       enddo

       end subroutine read_realmatrix_fromfile
