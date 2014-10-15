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
      subroutine prep_cg(n_node,n_Exp,round,wij,alb,n_obs,dir)
      use modoptions	!to assign the options  
      use modnonzero	!to assign numnze
      
      character index*100,char1,options*100,topology*100
      character(*) dir	!directory
      character(len=200) CommandFileName
      double precision alb(n_node)
      double precision :: wij(n_node,n_node),gamma(n_node)
      integer topo(n_node,n_node),n_node,n_exp,round,n_obs
      
      integer :: i,j
      
      character(len=5) char_numnze
      character(len=5) char_nnode
      character(len=5) char_nexp
      character(len=5) char_round
      character(len=5) char_nobs
      
c      write a temporary command file

	  !n_node = number of nodes
	  !n_Exp  = number of experiments
	  !round  = integer index of which decimated model
	  !wij    = Weight matrix before (input) GD
	  !alb    = alpha over beta (gamma?)
	  !n_obs  = number of observed nodes
	  
	  !******* EM: CHANGES ********************************************
	  ! (1) Removed swapping options betwen files and replaced with
	  !    	hardcoded parameters and direct assignments. This requires
	  !	   	the removal of the read_options call in fit_external_program
	  !
	  ! (2) Changed real wij to type 'double precision'.  This is ok
	  !    because it has to be turned into 'double precision' anyway
	  !    in the assignment to opt_Wo.
	  !
	  ! (3) Added opt_topology (instead of writing topology to file). 
	  !		This requires (a) the removal of the read_topology call
	  !		in fit_external_program, (b) adding opt_topology to 
	  !		modoptions.f
	  !
	  !	(4) Removed writing 'options_#.txt' file and 'topology_#.txt'
	  !		files. This requires removing them from the main.f program.
	  !
	  ! (5) Added 'dir' as an input argument, to hold the path to the 
	  ! 	'command.txt' file because now pineda is called from a 
	  !		different directory than than the source directory.
	  !
	  !	(6) Removed declaration of numnze from explicit declaration in 
	  !		the preamble.  Used modnonzero
	  !
	  !	(7) Removed writing command.txt as the data is already stored.
	  !
	  ! (8) Removed opt_lambda from here, fit_static_2,sim.f and modode
	  !
	  !	(9) Removed opt_sigmoid_type from here, fit_static_2, sim.f
	  !		and modode
	  !
	  !	(1) Removed opt_VM from here, fit_static_2, sim.f & modode
	  !****************************************************************
      write(index,*)round
       
      !HARDCODED PARAMETERS FOR GRADIENT DESCENT
      opt_ETA = 1.0000e-06	!learning rate for gradient descent
      opt_momentum = 0.50	!momentum rate for gradient descent
      opt_iterations = 20	!maximum number of iterations to take
      opt_termination_threshold = 0.010
      opt_L1 = 0.0010
        
      !ASSIGN OTHER OPTIONS(Wo,Alphao,Betao)
      do i=1,n_node
         opt_Wo(i,1:n_node) = wij(i,:)
         opt_ALPHAo(i) = 1.0000
         opt_BETAo(i) = 1/alb(i)
      enddo   
      
      
      !ASSIGN TOPOLOGY: SIGN() function is insufficient
      numnze=0
      do i=1,n_node
         do j=1,n_node
            if(wij(i,j).eq.0) then
               opt_topology(i,j) = 0
            else
               opt_topology(i,j) = 1
               numnze = numnze+1
               nzrows(numnze) = i
               nzcols(numnze) = j
            endif
         enddo
      enddo
      
      end !end of subroutine
      
        
       
       
       
       


       
       
