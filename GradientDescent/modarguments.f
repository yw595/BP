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
	   module modarguments
       character*150 arg_N_nodes
       character*150 arg_N_expts
       character*150 arg_N_observed       
       character*150 perturbfile    
       character*150 datafile
       character*150 omegafile
       character*150 topologyfile
       character*150 optionsfile
       character*150 output_tag
       !character*150 predfile
       character*150 arg_numnze
c       character*150 nze_file
       integer argument_N_nodes
       integer argument_N_expts
       integer argument_tag
       integer argument_numnze
       integer argument_n_observed      
       end module modarguments
