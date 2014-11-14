       subroutine testInd2Sub

       use mod_CONSTANTS

       integer i

       call ind2sub(5,7,testInds,testSubs,35)

       write(*,"(A20)") "current ind2sub:"

       do i=1,35
          write(*,"(I4,I4,I4)") testInds(i),
     :    testSubs(i,1),testSubs(i,2)
       enddo

       write(*,"(A20)") "previous ind2sub:"

       call ind2subprev(5,7,testInds,testSubs,35)

       do i=1,35
          write(*,"(I4,I4,I4)") testInds(i),
     :    testSubs(i,1),testSubs(i,2)
       enddo

       end
