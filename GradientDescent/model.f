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
      subroutine print_model(error,dir,tag)
      use modmodel
      use modarguments
      integer u1,u2,u3,u4,u5 
      character(*) tag
      character(*) dir
      integer i,j
      character*150 nameW,nameWx
      character*150 nameA,nameAx
      character*150 nameB,nameBx
      character*150 nameE,nameEx
      character*150 nameP
      double precision error   
       
       
      nameW=trim(adjustl(dir))//
     : trim(adjustl('/resultW_'))//trim(adjustl(tag))//'.txt'
      
      nameA=trim(adjustl(dir))//
     : trim(adjustl('/resulta_'))//trim(adjustl(tag))//'.txt'
     
      nameB=trim(adjustl(dir))//
     : trim(adjustl('/resultb_'))//trim(adjustl(tag))//'.txt'
     
      nameE=trim(adjustl(dir))//
     : trim(adjustl('/resultE_'))//trim(adjustl(tag))//'.txt'
      !write(nameP,50) 'result_P_',tag
      
 50   format (a9,i6)   
       
      u1=10
      u2=11
      u3=12
      u4=13
      !u5=14
      open (u1,FILE=nameW,STATUS='REPLACE') ! Open the file 
      open (u2,FILE=nameA,STATUS='REPLACE') ! Open the file 
      open (u3,FILE=nameB,STATUS='REPLACE') ! Open the file 
      open (u4,FILE=nameE,STATUS='REPLACE') ! Open the file 
      !open (u5,FILE=nameP,STATUS='REPLACE') ! Open the file 
      !write(*,*) 'N=',modelN
      !write(*,*) 'Model.type=',modeltype
      !write(*,*) 'Lambda=',modellambda
      !write(*,*) 'VM=',modelVM
      !write(*,*) 'Sigm.type=',modelsigmoidtype
      !write(*,30) 'alpha=',(modelalpha(i),i=1,modelN)
      !write(*,30) 'beta=',(modelbeta(i),i=1,modelN)
      do 10 i=1,modelN
         !write(*,20) 'W(',i,',:)=',(modelW(i,j),j=1,modelN)
         do j=1,(modelN-1)
            write(u1,'(F6.2$)') modelW(i,j)
         enddo   
         write(u1,'(F6.2)') modelW(i,modelN)
         write(u2,'(F6.2)') modelalpha(i)
         write(u3,'(F6.2)') modelbeta(i)
         !write(u5,*) ' ',(modelpred(i))
 10   continue
 20   format(A,I2,A,60F10.2)
 30   format(A,60F10.2)
 	  write(u4,'(F8.2)') error
      end
      
C     Print topology matrix
      subroutine print_topology(topology,N)
      
      use modexp
      use modmodel
      integer N
      double precision topology(modelmaxN,modelmaxN)

      do 10 i= 1,N
c         write(*,20) (topology(i,j),j=1,N)
 10   continue

 20   format(15F4.0)
      end
      

C     Read weight matrix from file
      subroutine read_weight(fname)
      use modmodel
      character*6 fname
      integer u
      parameter (u=21)
      
      open (u,FILE=fname,STATUS='OLD') ! Open the file
      
      do 10 i= 1,modelN             
         read(u,*) (modelW(i,j),j=1,modelN)
 10   continue
      
      
      end
      

      
C     Read topology matrix from file
      subroutine read_topology(topology,fname,N)
      use modexp
      use modmodel
	  double precision topology(modelmaxN,modelmaxN)

      character*150 fname
      integer u
      integer N
      parameter (u=21)
           
      open (u,FILE=fname,STATUS='OLD') ! Open the file
      write(*,*)'read_expdata filename: ',fname,' Nexp=',Nexp
      do 10 i= 1,modelN             
         read(u,*) (topology(i,j),j=1,N)
 10   continue

      write(*,*)'DONE read_topology file: ',fname
      
      end


C     This script generates a random model. 
C     Example:
C     M=random_model('hopfield2',5,2)
C     generates a 5-gene network with indegree of max 2
      
      subroutine random_model(indegree);
      use modmodel

      integer indegree
C     Model representation
C     type=0 'linear1'
C     type=1 'linear2'
C     type=2 'hopfield1'     
C     type=3 'hopfield2'      
C     type=4 'ssystem'
      double precision w(modelN,modelN)


      call randomnet(w,modelN,indegree) 
      
      if (modeltype.eq.0) then
         write(*,*) 'random_model, modeltype ',modeltype,' not impl.'
      elseif (modeltype.eq.1) then
         write(*,*) 'random_model, modeltype ',modeltype,' not impl.'
      elseif (modeltype.eq.2) then
         write(*,*) 'random_model, modeltype ',modeltype,' not impl.'
      elseif (modeltype.eq.3) then
         write(*,*) 'random_model, modeltype ',modeltype,' impl.'
         DO i=1,modelN
            DO j=1,modelN
               modelW(i,j)= w(i,j)
            END DO
            modelalpha(i)=1
            modelbeta(i)=1
         END DO
         modelVM=1.d0
         modellambda=1.d0
         modelsigmoidtype=1
      elseif (modeltype.eq.4) then
         write(*,*) 'random_model, modeltype ',modeltype,' not impl.'
      endif

      end
C     W=randomnet(n,k)
C     calculates a n node network with indegrees in the range 1-k
C     where nonzero W entries are either +1 or -1 with equal probability
C
C     See also: perturbationset, forward_eval
      subroutine randomnet(W,n,k) 
     
      integer k,n
      double precision W(n,n)
      integer r(n),v(n),temp(2)
      integer i,j,n_inputs
      debug=0
      
      do 20, i=1,n
         do 10, j=1,n
            W(i,j)=0;
 10      continue
 20   continue
      
      do 40, i=1,n
         call randperm(k,r) 
         n_inputs=r(1);
         call randperm(n,v) 
         do 30, j=1,n_inputs 
            W(i,v(j))=1.d0
            call randperm(2,temp)
            if (temp(1).gt.1) W(i,v(j))=-1.d0
 30      continue
 40   continue
      
      
      end
