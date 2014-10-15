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
         subroutine prep_data(xdata,udata,nexpcg,n_node,n_Exp,n_obs
     :		,ProtObs)
c		compute correlation data and select based on this
         use modexp			!to declare Nexp
         use modmodel		!to declare modelN
         integer n_exp,n_node,nexpcg,n_obs
         double precision xdata(n_node,n_exp),udata(n_node,n_exp)
         double precision xdatat(n_node,n_exp),udatat(n_node,n_exp)
         double precision xdata_X(n_node,n_exp),udata_x(n_node,n_exp)                
         integer    oDATACG(n_node,n_exp)
         integer ProtObs(n_node)
         
                  
         !************  EM: CHANGES  ***********************************
         !
         !	(1) I changed the format from auto (*) to F6.2, as it is 
         !		consistent with the way all other formatted data files 
         !		are read and written.
         !
         ! 	(2) There was a problem here, in that the subroutine assumes
         !		that the nodes are sorted with observed nodes preceeding
         ! 		the unobserved nodes; this is not necessarily true and we
         !		can't assume that the user will have sorted their model
         !		nodes in this way.  Included ProtObs argument. 
         !
         !	(3) Changed the dimension of xdata and udata to be of correct
         !		size. This fixed the problem with incorrect selection of
         !		the least correlated experiment indices.
         !		
         !	(4) Hardcoded to only accept EVEN NUMBER for the number of
         !		experiments to keep (nexpcg).  This has to be change for
         !		flexibility.  Right now, the program terminates. I was
         !		not able to add one to nexpcg without running into a 
         !		"bus error".
         !
         !	(5) Instead of passing the number of reduced experiments as
         !		an argument in "command.txt" I just assign to the 
         !		the variable defined in the relevant module (modexp).
         !
         !  (6) Use modmodel to declare the number of model nodes
         !
         !**************************************************************
         
         !Assign the new number of experiments (trimmed data) into 
         !modexp so that pineda knows how many experiments it is fitting to.
         Nexp = nexpcg 		!Nexp is stored in modexp
         modelN = n_node	!modelN is stored in modmodel
         nobs = n_obs		!nobs is stored in modmodel
         
         if (modulo(nexpcg,2).ne.0) then
            write(*,*) "-----------WARNING-------------------"
            write(*,*) "You must enter an EVEN number of" 
     :       //" experiments to keep."
     		write(*,*) ""
     		write(*,*) "PROGRAM TERMINATED."
     		write(*,*) "--------------------------------------"
     		stop
         endif   
         
         open(1,file="data_x.txt")
         open(2,file="data_u.txt")
         open(3,file="data_o.txt")
c                  do j=1,n_exp
c                  write(44,*)(xdata(i,j),i=1,n_node)
c                  enddo
         do i=1,n_node
         do j=1,n_exp
         xdatat(i,j)=xdata(i,j)
         udatat(i,j)=udata(i,j)
         enddo
         enddo
         
         
         call 
     :corr_sel(xdatat,udatat,nexpcg,n_node,n_Exp,xdata_x,udata_x,n_obs)

         do i=1,n_node
         do j=1,nexpcg

         if(ProtObs(i).eq.1)then	
         	oDATACG(i,j)=1
         else
            oDATACG(i,j)=0
         endif
         enddo
         enddo
         
         
         do i=1,n_node
            do j=1,(nexpcg-1)
               write(1,'(F6.2$)')xdata_x(i,j)
          	   write(2,'(F6.2$)')udata_x(i,j)
          	   write(3,"(I2$)")oDATACG(i,j)
          	enddo
          	write(1,'(F6.2)')xdata_x(i,nexpcg)
          	write(2,'(F6.2)')udata_x(i,nexpcg)
          	write(3,'(I2)')oDATACG(i,nexpcg)
         enddo 	
         close(1)
         close(2)
         close(3)
         end
         		
         subroutine corr_sel(xdata,
     :udata,nexpcg,n_node,n_Exp,xdata_x,udata_x,n_obs)
         integer n_exp,n_node,nexpcg,n_obs
         double precision xdata(n_node,n_exp),udata(n_node,n_exp)
         double precision xdatat(n_node,n_exp),udatat(n_node,n_exp)  
         double precision xdata_x(n_node,n_exp),udata_x(n_node,n_exp)            
         double precision s_xy(n_exp,n_exp),sx_sy(n_exp,n_exp)
         double precision sx_2(n_exp),s_x2(n_exp),sx(n_exp),std(n_exp)
         double precision corr(n_exp,n_exp),xy(n_node,n_exp,n_exp)
         double precision cormin
         integer l1(n_exp),l2(n_exp),t,l
         character(len=50) Fname_SubsetIdx
         character(len=2) string_nexpcg
         
		 !start calculating co-variance: x(i,j)*x(i,k), j.ne.k
         do k=1,n_exp
         do j=1,n_exp
         do i=1,n_node
         xy(i,k,j)=xdata(i,j)*xdata(i,k)
c         print*,xy(i,k,j),i,k,j
         enddo
         enddo
         enddo        
         
         !continue calculating co-variance: sum over all nodes
         do k=1,n_exp
         do j=1,n_exp
         s_xy(k,j)=sum(xy(1:n_obs,j,k))
         sx_sy(k,j)=sum(xdata(1:n_obs,k))*sum(xdata(1:n_obs,j))
         enddo
         
         s_x2(k)=sum(xdata(1:n_obs,k)**2)
         sx_2(k)=(sum(xdata(1:n_obs,k)))**2
         sx(k)=sum(xdata(1:n_obs,k))
         std(k)=sqrt(n_obs*s_x2(k)-sx_2(k))
c                  print*,s_x2(k),sx_2(k),sx(k),std(k)
         enddo
         
         do k=1,n_exp
         do j=1,n_exp
         corr(k,j)=(n_obs*s_xy(k,j)-sx(k)*sx(j))/(std(k)*std(j))
         if (corr(k,j).ne.corr(k,j)) then
            print*,k,j,std(k),std(j)
         endif
         enddo
         enddo
         
         
         
         do l=1,int(nexpcg*0.5)
         cormin=minval(abs(corr(1:n_exp,1:n_exp)))
         !print*,'******',cormin,l
         do k=1,n_exp
          do j=1,n_exp

          if(cormin.eq.abs(corr(k,j)))then
                   !print*,k,j
          l1(2*l-1)=k
          l1(2*l)=j

c	Exclude the already selected k & j from correlation analysis         
          do t=1,n_exp
          corr(t,k)=1
          corr(k,t)=1
          corr(t,j)=1
          corr(j,t)=1
          enddo
          go to 23
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc          
          endif
         enddo         
         enddo
23       continue         
         enddo 
         
         !shuffle the order depending on correlation
         write(string_nexpcg,"(I2)") nexpcg 
         
         Fname_SubsetIdx = "LeastCorrelatedIndices_"//
     :		trim(adjustl(string_nexpcg))//".txt"
         open(5,file=Fname_SubsetIdx,status="REPLACE")
         do j=1,nexpcg
            if (l1(j).lt.10) then
               write(5,'(I1A1$)') l1(j),' '
            else   
               write(5,'(I2A1$)') l1(j),' '
            endif   
         enddo
         close(5)
         
         
         do j=1,nexpcg
          do i=1,n_node
         xdata_x(i,j)=xdata(i,l1(j))
         udata_x(i,j)=udata(i,l1(j))
         enddo
         enddo
         end
               
             
