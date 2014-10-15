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
      module mod_DECIMATION
      implicit none
      
      integer ModelsDone	!The number of decimated models completed
      double precision, dimension(:,:,:), allocatable :: DeciMarginals  	!A copy of the full Marginals
      double precision, dimension(:,:), allocatable :: ModelW	!Matrix of real-valued weights
      integer, dimension(:,:), allocatable 			:: FixedW	!Matrix of integer indicators (0) unfixed, (1) fixed
      integer fxd_target  	!the target of the fixed interaction
      integer fxd_reg		!the regulator of the fixed interaction
      integer NumToFix		!the number of interactions to fix at any point throughout the decimation
      
      integer, dimension(:,:), allocatable 			:: FixOrder	!order in matrix subscripts
      integer, dimension(:), allocatable :: InteractionOrder   !order in matrix indices
      integer interaction
      
      integer, dimension(:), allocatable :: consecutive_fixes	!keeps track of the number of consecutive fixes inside decimation, without recalculating BP
      
      double precision, dimension(:), allocatable :: MCrandoms
      integer widx
      
      integer, dimension(:,:), allocatable :: Tracker
      integer NumFixed	
      integer NumUnfixed
      character(len=10) ModelName
      character(len=4) ModelIDX
            
            
      end module mod_DECIMATION