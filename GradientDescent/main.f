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

      subroutine pineda(dir,tag)
      !implicit none
      use modexp
      use modmodel
      use modoptions
      integer :: i,j
      character*150 fname
      double precision error
      double precision omega(expmaxN,expmaxNexp)	!consider putting in modexp and assigning in prep_data or prep_cg
      character(len=*) tag
      character(len=*) dir
      ! This reads arguments from file "command.txt"
      
      !******* EM: CHANGES ********************************************
	  ! (1) Removed all read and write arguments.  Kept only
	  !			(1a) read_perturb 
	  !			(1b) read_expdata
	  !		These were included so we don't have to fit_external_problem.
	  !
	  ! (2) Removed output_tag because right now the tag is not used.
	  !
	  ! (3) SKIPPED call to fit_external_program and went straight to
	  !		fit_static_2. This is because when pineda is called within
	  !		the bp_main all of the necessary arguments/parameters/
	  !		structures are stored internally and we don't need the 
	  !		fit_external_problem function to prepare and store the data.
	  !
	  !****************************************************************
      
      !read trimmed data: pert matrix and data matrix
      fname = "data_u.txt"
	  call read_perturb(fname)
	  fname = "data_x.txt"
	  call read_expdata(fname)
	  fname = "data_o.txt"
	  call read_omega(omega,fname) 
	  
	  
	  call fit_static_2(omega,opt_topology,error)
	  call print_model(error,dir,tag)
	  
	  
      
      
      
      end
