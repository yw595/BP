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
C     Corresponds to M structure in Matlab 
C     (but with vectors of alpha and beta instead for more generality)
      module modmodel
      implicit none
      save
      integer modelmaxN
      parameter (modelmaxN=100)      
      integer modelN                 ! size of the network
      double precision  modelW(modelmaxN,modelmaxN)
      double precision  modelalpha(modelmaxN)
      double precision  modelbeta(modelmaxN)
      double precision  modelVM
      double precision  modellambda
      double precision  modelpred(modelmaxN)
      integer modelsigmoidtype
      integer nobs
      integer modeltype
      end module modmodel


