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
      module modoptions
      implicit none
      save
      integer opt_maxN
      parameter (opt_maxN=100)      ! size of the network
      
      double precision  opt_ETA
      double precision  opt_momentum
      double precision  opt_termination_threshold
      integer opt_iterations
      double precision  opt_Wo(opt_maxN,opt_maxN)
      double precision  opt_topology(opt_maxN,opt_maxN)
      double precision  opt_ALPHAo(opt_maxN)
      double precision  opt_BETAo(opt_maxN)
      double precision  opt_prior(opt_maxN)!?
      double precision  opt_L1
      integer opt_verbose
      integer opt_plot
      double precision opt_temperature
      double precision opt_gamma
      character opt_resultpath
      
    
      
      end module modoptions
