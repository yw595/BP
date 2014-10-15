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

C     Print perturbation matrix
      subroutine print_perturb()
      
      use modexp
      use modmodel

      do 10 i= 1,modelN
c         write(*,20) (p(i,j),j=1,Nexp)
 10   continue

 20   format(15F4.0)
      end

C     Print prediction matrix
      subroutine print_pred()
      
      use modexp
      use modmodel

      do 10 i= 1,modelN
c         write(*,20) pred(i,1)
 10   continue

 20   format(15F4.0)
      end

C     Print data matrix
      subroutine print_data()
      
      use modexp
      use modmodel
      
      do 10 i= 1,modelN
c         write(*,20) (z(i,j),j=1,Nexp)
 10   continue
 20   format(15F6.2)
      end
   
C     Print omega matrix
      subroutine print_omega(omega)
      
      use modexp
      use modmodel
      double precision omega(expmaxN,expmaxNexp)
            
      do 10 i= 1,modelN
c         write(*,20) (omega(i,j),j=1,Nexp)
 10   continue
 20   format(15F4.0)
      end
   

C     Read perturbation matrix form file
	  !**** KEPT IN EVANS CHANGES ***
      subroutine read_perturb(fname)
      use modmodel		!contains modelN
      use modexp		!contains Nexp

      character*150 fname
      integer u
      parameter (u=21)
      
      open (u,FILE=fname,STATUS='OLD') ! Open the file
      do 10 i= 1,modelN             
         read(u,*) (p(i,j),j=1,Nexp)
 10   continue

     
      end
      
C     Read prediction matrix form file
      subroutine read_pred(fname)
      use modexp
      use modmodel

      character*150 fname
      integer u
      parameter (u=21)
      
      open (u,FILE=fname,STATUS='OLD') ! Open the file
      !write(*,*)'read_perturb file: ',fname,' Nexp=',Nexp
      do 10 i= 1,modelN             
      !   write(*,*) i
          read(u,*) pred(i,1)
      !   write(*,*)(p(i,j),j=1,Nexp) 
 10   continue

     
      end
            
      
C     Read experimental data from file
      subroutine read_expdata(fname)
      use modexp
      use modmodel

      character*150 fname
      integer u
      parameter (u=21)
      
           
      open (u,FILE=fname,STATUS='OLD') ! Open the file
      !write(*,*)'read_expdata filename: ',fname,' Nexp=',Nexp
      do 10 i= 1,modelN             
      !   write(*,*) i
         read(u,*) (z(i,j),j=1,Nexp)
       !  write(*,*)(z(i,j),j=1,Nexp) 
      
 10   continue
      
      end
      
      
C           Read omega matrix form file
      subroutine read_omega(omega,fname)
      use modexp
      use modmodel

      double precision omega(expmaxN,expmaxNexp)
      character*150 fname
      integer u
      parameter (u=21)
      
           
      open (u,FILE=fname,STATUS='OLD') ! Open the file
      do 10 i= 1,modelN     
         read(u,*) (omega(i,j),j=1,Nexp)
      
 10   continue


      end
      
      
         
      
C     Read options from file
      subroutine read_options(fname)

      use modoptions
      use modmodel
      character*150 fname    
      integer u,x
      parameter (u=21)
      
      opt_sigmoid_type=1
      write(*,*) 'the model N is ',modelN
      open (u,FILE=fname,STATUS='OLD') ! Open the file
c      read(u,*) opt_resultspath
      read(u,*) opt_eta
      read(u,*) opt_momentum
      read(u,*) opt_iterations
      read(u,*) opt_termination_threshold
      read(u,*) opt_L1
      do 10 i=1,modelN
         read(u,*) (opt_Wo(i,j),j=1,modelN)
 10   continue

      do 40 i=1,modelN
         read(u,*) opt_ALPHAo(i)
 40   continue
      do 50 i=1,modelN
         read(u,*) opt_BETAo(i)
 50   continue
 
         
      
      
      write(*,30) 'Learning rate    : ',opt_eta
      write(*,20) 'Learning momentum: ',opt_momentum
      write(*,*) 'Max iterations   : ',opt_iterations
      write(*,20) 'Termination thres: ',opt_termination_threshold
      write(*,20) 'L1 coefficient   : ',opt_L1
      write(*,20) 'opt_Wo           : '
      do i=1,modelN
c         write(*,*)  (opt_Wo(i,j),j=1,modelN)
      enddo
    
 20   format(A,30F10.3)
 30   format(A,30F10.7)
      end
      
      
