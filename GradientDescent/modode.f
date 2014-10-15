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
      module modode
      implicit none
      save
      
      !module containing variables for simulatin a single model given
      !an experimental condition.
      
      integer ode_maxN
      parameter (ode_maxN=100)      ! size of the network
      
      integer ode_N_NODES           ! Number of nodes after reduction 
      integer ode_READOUTS(ode_maxN)   ! Observed nodes (for the current experiment)
      integer ode_STRUCTURE(ode_maxN,ode_maxN)  ! w_ij's allowed a non-zero value.
      
      
      !for pineda_ode (model equation)
      double precision ode_I(ode_maxN)      ! Perturbation bias 
      double precision ode_tau(ode_maxN) ! Training example 
      double precision ode_WEIGHTS(ode_maxN,ode_maxN) ! Interaction weights
      double precision ode_alpha(ode_maxN) ! Half-life parameter
      double precision ode_beta(ode_maxN) ! Response amplitude parameter
      
      !for pineda_ode_final
      double precision ode_ETA      !  Learning rate
      double precision ode_Y0(ode_maxN)      !  Initial values
      
      
      
      end module modode
      
