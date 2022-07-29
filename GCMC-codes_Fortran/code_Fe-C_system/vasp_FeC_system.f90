      !Primary verison of 20210112 for VASP (Only relaxation; No AIMD yet)
      !Keep fixed atoms at the begining of structure file (bottom). Keep Fe goes before C in all structures file (generally (q/n)atom1 refer to Fe and atom2 refer to C)
      !Prepare different POTCARs in head directory (for subroutine runVasp)
      !First line of POSCAR: "Fe C". So that Ccount is the indicator of whether there is C in the system (if Ccount>1)
      
      !assume unit cell to be orthogonal
      !Only Fe-C system (with both Fe/C atoms as MC_atoms) in current version
      !Only the use of relaxation (but not AIMD) is considered in current version
      
      !orthogonal cell is assumed in subroutine(getRandomCoord)
      
      program GCMC_Vasp
      implicit none
      integer, parameter :: KREAL = kind(0.d0)
      integer, parameter :: nmcmax=10000  !!! this number should be consistent with the parameters in the two subrouthine 
      integer, parameter :: ntypemax=2
      integer, save   :: gSeed
      real(KREAL)  :: mass(ntypemax), wavelength(ntypemax), mu(ntypemax)
      real(KREAL) :: cmc(nmcmax,3), cmc_System(nmcmax,3), xaxismc_System(3), yaxismc_System(3), zaxismc_System(3),&
          & randCoord(3), &
          & cmc_trial(nmcmax,3), dist(3), &
          & xaxismc(3), yaxismc(3), zaxismc(3), abcrossproduct(3) 
      real(KREAL) :: temp, Rmax(2), Rmin(2), Pressure, &
          & emc, beta, &
          & distmin, E_System, rand, vacuum, &
          & volume, volume_System, prob, prob2, MCatomRemoved, deltaV, &
          & MCatomMoved, dista, zmin, zmax, disp_prob
      character*2 :: qatom(ntypemax), qadd(ntypemax)
      character*1 :: tfmc(nmcmax,3), tfmc_trial(nmcmax,3), tfmc_System(nmcmax,3)   
      character*5 :: filename
      character*11 :: filename2, rtime, rdate
      integer :: nMCatoms_System(ntypemax), nMCatoms(ntypemax)
      integer :: nmciter, iRX, imc, nfail, types, iter, &
         & ndumb, natoms, natoms1, natoms2, natoms1_System, natoms2_System, natoms_System, niterMC, nPermanentAtoms_System, &
         & iMCflag, icheck, nMCatomRemoved,iMCpick, iAtom, ipass, &
         & nMCatomMoved, i1, i2, i3, i4, iMCcheck, iTail, itype, nAtomTypes, imagetest_x, imagetest_y,&
         & imagetest_z, idebug, iter_cmc, type_check, natoms_trial, natoms1_trial, natoms2_trial, iz_place, Ccount
      logical, save :: debug
     
     debug = .TRUE.

!****************************************
!*     initialize MC parameters
!****************************************
!    Read Control_MC parameters
      if (debug) print*, 'start reading and writing controls'
      open (1,file='herefile',access='append')
           write (1,*) "start reading and writing controls"
      close (1)
      open (22,file='control_MC',status='old')
      read (22,'(i7)')   idebug        !debug switch
      read (22,'(i7)')   nmciter       !Number of MC trial iterations before stopping
      read (22,'(f8.2)') temp          !system temperature
      read (22,'(f8.2)') vacuum        !Vacuum space for slab/cluster geometries
      read (22,'(i7)')   nAtomTypes    !nuber of GC atom types
      do itype = 1,nAtomTypes
         read (22,'(a2)') qadd(itype)  !Atom types
      end do
      do itype = 1,nAtomTypes
         read (22,'(f8.2)')   mass(itype)  !atomic mass of each atom type
         wavelength(itype) = sqrt(1/(mass(itype)*temp))*17.45801643 !thermal de Broglie wavelength in angstroms
      end do
      do itype = 1,nAtomTypes
         read (22,'(f8.2)') mu(itype)  !chemical potential of each atom type
      end do 
      do itype = 1,nAtomTypes
         read (22,'(i7)')   nMCatoms_System(itype)  !initial number of each atom type (excluding fixed ones)
      end do   
      do itype = 1,nAtomTypes
         read (22,'(f8.2)')   Rmin(itype)  !Min radius for atom placement
      end do   
      do itype = 1,nAtomTypes
         read (22,'(f8.2)')   Rmax(itype)  !Max radius for atom placement   
      end do
      read (22,'(i7)')     iz_place        !Place atoms in a range of z? 0=no 1=yes
      read (22,'(f8.2)')   zmin            !min. z from where placement of atoms starts
      read (22,'(f8.2)')   zmax            !max. z where placement of atoms ends
      close (22)
      
      !Turn on debug switch for additional output if option selected
      if (idebug.eq.1) debug = .TRUE. 
      
!initalize random number generator seed
      call init_random_seed()

!initialize global variables 
     niterMC = 0 !iteration counters
     nfail = 0
     imc = 0
     beta = 1 / (0.00008617333262*temp) !Boltzmann Constant: 8.617333262×10−5 eV/K
     
!arrange files in directory and run initial VASP on System 
      
      if (debug) print*, 'start system calls'
      
      !call system
      !call system ("sed -i 's/NSW.*/NSW   = 0/g' INCAR")  !!run single point calculation
      !call system ('srun vasp_std')

      call ReadCoordinates(emc,cmc,tfmc,natoms1,natoms2,natoms,xaxismc,yaxismc,zaxismc)   !!start from a calculation result directly
      

!initialize system variables (store them from the output of ReadCoordinates)
      E_System = emc                      !initial Reax energy
      cmc_System = cmc                    !initial atom coordinates
      tfmc_System = tfmc
      natoms1_System = natoms1
      natoms2_System = natoms2
      natoms_System = natoms              !inital number of atoms
      xaxismc_System = xaxismc            !initial cell axis
      yaxismc_System = yaxismc
      zaxismc_System = zaxismc 
      
      nPermanentAtoms_System = natoms_System - sum(nMCatoms_System) !number of Permanent atoms that can't be removed
      
      abcrossproduct(1) = xaxismc(2)*yaxismc(3)-xaxismc(3)*yaxismc(2)
      abcrossproduct(2) = xaxismc(1)*yaxismc(3)-xaxismc(3)*yaxismc(1)
      abcrossproduct(3) = xaxismc(1)*yaxismc(2)-xaxismc(2)*yaxismc(1) 
      volume_System = sqrt(dot_product(abcrossproduct, abcrossproduct))*zaxismc(3)   !Non-orthogonal but z-axis has to be normal
      
!create log file (MClog)
      open (17, file='MClog')
          write (17,'(5a7)')'MCit','I','nAt1','nAt2','E'
          write (17,'(4i7,f14.8)')niterMC,imc,nMCatoms_System(:), E_System 
      close (17)

! (geolog)       
      open (16, file='geolog')
      write (16,'(4i7,f14.8)')niterMC,imc,nMCatoms_System(:), E_System 
      close (16)
      call system ('cat CONTCAR >> geolog')

! (fail-log)       
      open (15, file='faillog')
      write (15,'(4i7,f14.8)')nfail,imc,nMCatoms_System(:), E_System 
      close (15)
      call system ('cat CONTCAR >> faillog')
      
      call system ('rm CHG* CONTCAR EIGENVAL IBZKPT logfile OSZICAR OUTCAR PCDAT POSCAR REPORT vasprun.xml WAVECAR XDATCAR')  !! clean up files after the OUTCAR/CONTCAR have been strored
      
       open (1,file='herefile',access='append')
           call date_and_time(rdate,rtime)
           write (1,'(a50,i6,a15)') 'Initial vasp ran, starting MC loops',natoms_System,rtime
       close (1)
     
!****************************************
!*     MC LOOP
!****************************************
Do imc=1,nmciter              !START MC LOOP
     if (debug) print *,'Iteration:', imc
     
     ! select MC move
     !***************

     call Random_Number(rand)

     if (3*rand.lt.1.0) iMCflag = 0                         !! displace
     if ((3*rand.gt.1.0).and.(3*rand.lt.2.0)) iMCflag = 1     !! Remove
     if (3*rand.gt.2.0) iMCflag = 2                         !! Add
     
!*****************************************
!*     Displace Atom Block 
!*****************************************
     call Random_Number(disp_prob)     !! Determine the (artificial) ratio of accepted displacement.  
     if ((iMCflag.eq.0).and.((natoms_System-nPermanentAtoms_System).gt.0).and.(disp_prob.gt.0)) then     !!half of displacements are accepted                       
         if (debug) print*, 'displace atom'
         open (1,file='herefile',access='append')
           call date_and_time(rdate,rtime)
           write (1,'(i10,a40,i6,a15)') imc,"Displace", natoms_System, rtime
         close (1)

         !Derive trial geometry by displacement of one atom position
         !***********************************************************        
         cmc_trial = cmc_System

         call Random_Number(rand)
         nMCatomMoved = floor(rand*(natoms_System-nPermanentAtoms_System)+ nPermanentAtoms_System)+1   !! select random MC atom to move (rand0 -> select first unfixed atom)
         
         
         if (nMCatomMoved.gt.natoms1_System) then                   !determine which atom type was selectedn
               iAtom = 2                    
         else
               iAtom = 1
         end if         
         !! figure out what atom is displaced
         
call GetRandomCoords(natoms_System,cmc_System,Rmin,Rmax,xaxismc,yaxismc,zaxismc,nmcmax,randCoord,iz_place,zmin,zmax,iAtom, ipass)
         
     if (ipass.eq.1) then 
         cmc_trial(nMCatomMoved,:) = randCoord(:) ! replace with random coord
         tfmc_trial = tfmc_System
         natoms1_trial = natoms1_System
         natoms2_trial = natoms2_System
         natoms_trial = natoms_System
         
         !Run VASP
         !**********
         call Writecoordinates(xaxismc,yaxismc,zaxismc,cmc_trial,tfmc_trial,natoms1_trial,natoms2_trial,natoms_trial)
         call Writemagmom(natoms1_trial, natoms2_trial)
         call RunVasp(natoms_trial) ! run Vasp with new geo
         call ReadCoordinates(emc,cmc,tfmc,natoms1,natoms2,natoms,xaxismc,yaxismc,zaxismc) !read in new Reax data
         
         !Determine Acceptance Probablity
         !*******************************
         prob = Min(1.0,Exp(-beta*(emc-E_System)))         
         
     else if(ipass.eq.0) then
         prob = 0
     end if
     
         call Random_Number(prob2)         
         
         !Apply Acceptance Updates if Accepted
         !*************************************
         if (prob.lt.prob2) then
             nfail = nfail + 1
             
             open (15, file='faillog', access='append')
             write (15,'(4i7,f14.8)')nfail,imc,nMCatoms_System(:), E_System 
             close (15)             
             call system ('cat CONTCAR >> faillog') !log failed geometry step              
             
             call system ('rm CHG* CONTCAR EIGENVAL IBZKPT logfile OSZICAR OUTCAR PCDAT POSCAR REPORT vasprun.xml WAVECAR XDATCAR')  !! clean up files 
             
             open (1,file='herefile',access='append')
               call date_and_time(rdate,rtime)
               write (1,'(i10,a30,f10.4,i6,a15)') imc, "Displace fail",prob, natoms_System, rtime
             close (1)
         else
             niterMC = niterMC + 1
             
             open (16, file='geolog', access='append')
             write (16,'(4i7,f14.8)')niterMC,imc,nMCatoms_System(:), E_System 
             close (16)             
             call system ('cat CONTCAR >> geolog') !log accepted geometry step
             
             call system ('rm CHG* CONTCAR EIGENVAL IBZKPT logfile OSZICAR OUTCAR PCDAT POSCAR REPORT vasprun.xml WAVECAR XDATCAR')  !! clean up files
             
             E_System = emc                       !! Store new E/geometric informations (accept the change)              
             cmc_System = cmc                 
             tfmc_System = tfmc
             natoms1_System = natoms1
             natoms2_System = natoms2
             natoms_System = natoms   

             open (17, file='MClog', access='append')  !! Write log
                 write (17,'(4i7,f14.8)')niterMC,imc,nMCatoms_System(:), E_System 
             close (17)

             open (1,file='herefile',access='append')
               call date_and_time(rdate,rtime)
               write (1,'(i10,a30,f10.4,i6,a15)') imc, "Displace pass", prob, natoms_System, rtime
             close (1)
         end if
         
!*  End Displace Atom Block   
!********************************************

!********************************************
!*     Remove Atom Block
!********************************************
      elseif ( (iMCflag.eq.1).and.((natoms_System-nPermanentAtoms_System).gt.0) ) then     !skip if no exchangeable atoms left to remove
         if (debug) print*, 'remove atom'
         
         open (1,file='herefile',access='append')
           call date_and_time(rdate,rtime)
           write (1,'(i10,a40,i6,a15)') imc,"Remove", natoms_System, rtime
         close (1)
            
         !Remove atom 
         !************        
         call Random_Number(rand)

         nMCatomRemoved = floor(rand*(natoms_System-nPermanentAtoms_System))+nPermanentAtoms_System+1    !randomly select a MC atom to remove
                                             
         if (nMCatomRemoved.gt.natoms1_System) then                   !determine which atom type was selected for deletion
               natoms2_trial = natoms2_System -1
               natoms1_trial = natoms1_System
               iAtom = 2                    
         else
               natoms1_trial = natoms1_System -1
               natoms2_trial = natoms2_System
               iAtom = 1
         end if        
         
         cmc_trial = cmc_System                            
         tfmc_trial = tfmc_System
         natoms_trial = natoms_System - 1
         
         do  i3 = nMCatomRemoved,natoms_trial              !!change coordinates of following atoms
               cmc_trial(i3,:) = cmc_System(i3+1,:)
               tfmc_trial(i3,:) = tfmc_System(i3+1,:)
         end do 
         
         
         !Run VASP
         !**********
         call Writecoordinates(xaxismc,yaxismc,zaxismc,cmc_trial,tfmc_trial,natoms1_trial,natoms2_trial,natoms_trial)
         call Writemagmom(natoms1_trial, natoms2_trial)
         call RunVasp(natoms_trial) ! run Vasp with new geo
         call ReadCoordinates(emc,cmc,tfmc,natoms1,natoms2,natoms,xaxismc,yaxismc,zaxismc) !read in new Reax data
                
         !Determine Acceptance Probablity
         !*******************************
prob = Min(1.0,((nMCatoms_System(iAtom)*(wavelength(iAtom)**3))/(volume_System-vacuum))*Exp(-beta*((emc-E_System)+mu(iAtom))))         
         
         call Random_Number(prob2)
         
         !Apply Acceptance Updates if Accepted
         !*************************************
         if (prob.lt.prob2) then
             nfail = nfail + 1
             
             open (15, file='faillog', access='append')
             write (15,'(4i7,f14.8)')nfail,imc,nMCatoms_System(:), E_System 
             close (15)                          
             call system ('cat CONTCAR >> faillog') !log failed geometry step               
             
             call system ('rm CHG* CONTCAR EIGENVAL IBZKPT logfile OSZICAR OUTCAR PCDAT POSCAR REPORT vasprun.xml WAVECAR XDATCAR')  !! clean up files
             
             open (1,file='herefile',access='append')
               call date_and_time(rdate,rtime)
               write (1,'(i10,a30,f10.4,i6,a15)') imc, "Remove fail",prob, natoms_System, rtime
             close (1)
         else
             niterMC = niterMC + 1
             
             open (16, file='geolog', access='append')
             write (16,'(4i7,f14.8)')niterMC,imc,nMCatoms_System(:), E_System 
             close (16)                          
             call system ('cat CONTCAR >> geolog') !log accepted geometry step
             
             call system ('rm CHG* CONTCAR EIGENVAL IBZKPT logfile OSZICAR OUTCAR PCDAT POSCAR REPORT vasprun.xml WAVECAR XDATCAR')  !! clean up files
             
             E_System = emc                       !! Store new E/geometric informations (accept the change)             
             cmc_System = cmc                 
             tfmc_System = tfmc
             natoms1_System = natoms1
             natoms2_System = natoms2
             natoms_System = natoms   
             nMCatoms_System(iAtom) = nMCatoms_System(iAtom) - 1

             !write updated log files
             open (17, file='MClog', access='append')
                 write (17,'(4i7,f14.8)')niterMC,imc,nMCatoms_System(:), E_System 
             close (17)

             open (1,file='herefile',access='append')
               call date_and_time(rdate,rtime)
               write (1,'(i10,a30,f10.4,i6,a15)') imc, "Remove pass",prob, natoms_System, rtime
             close (1)
         end if
!*  End Remove Atom Block   
!********************************************

!********************************************
!*     Add Atom Block
!********************************************         
         elseif (iMCflag.eq.2) then
         
         if (debug) print*, '3 add atom'
         open (1,file='herefile',access='append')
           call date_and_time(rdate,rtime)
           write (1,'(i10,a40,i6,a15)') imc,"Add", natoms_System, rtime
         close (1)
         
         !Add atom to system
         !********************
         call Random_Number(rand) !randomly choose atom type
         iAtom = floor(rand*nAtomTypes) + 1
         iAtom = 2
         
         call GetRandomCoords(natoms_System,cmc_System,Rmin,Rmax,xaxismc,yaxismc,zaxismc,nmcmax,randCoord,iz_place,zmin,zmax,iAtom, ipass)

     if (ipass.eq.1) then
     
         cmc_trial = cmc_System
         tfmc_trial = tfmc_System
         natoms_trial = natoms_System + 1
         
         if (iAtom.eq.1) then     !!Add Fe
               natoms1_trial = natoms1_System +1
               natoms2_trial = natoms2_System
               if (natoms2_System.ne.0) then
                    do i4 = natoms1_trial+1,natoms_trial
                         cmc_trial(i4,:) = cmc_System(i4-1,:)
                         tfmc_trial(i4,:) = tfmc_System(i4-1,:)
                    end do
                    
                    cmc_trial(natoms1_trial,:) = randCoord(:) ! add random MC atom coords after Fe atoms
                    tfmc_trial(natoms1_trial,1) = 'T'
                    tfmc_trial(natoms1_trial,2) = 'T'
                    tfmc_trial(natoms1_trial,3) = 'T'
                    
               elseif (natoms2_System.eq.0) then
                    cmc_trial(natoms1_trial,:) = randCoord(:) ! add random MC atom coords to the end
                    tfmc_trial(natoms1_trial,1) = 'T'
                    tfmc_trial(natoms1_trial,2) = 'T'
                    tfmc_trial(natoms1_trial,3) = 'T'
                    
               end if     
         elseif (iAtom.eq.2) then   !!Add C
               natoms1_trial = natoms1_System
               natoms2_trial = natoms2_System +1
               cmc_trial(natoms_trial,:) = randCoord(:) ! add random MC atom coords to the end
                    tfmc_trial(natoms_trial,1) = 'T'
                    tfmc_trial(natoms_trial,2) = 'T'
                    tfmc_trial(natoms_trial,3) = 'T'
                    
         end if
         

         !Run VASP
         !**********
         call Writecoordinates(xaxismc,yaxismc,zaxismc,cmc_trial,tfmc_trial,natoms1_trial,natoms2_trial,natoms_trial)
         call Writemagmom(natoms1_trial, natoms2_trial)
         call RunVasp(natoms_trial) ! run Vasp with new geo
         call ReadCoordinates(emc,cmc,tfmc,natoms1,natoms2,natoms,xaxismc,yaxismc,zaxismc) !read in new Reax data
         
         !Determine Acceptance Probablity
         !*******************************
prob = Min(1.0,(((volume_System-vacuum)/((nMCatoms_System(iAtom)+1)*(wavelength(iAtom)**3)))*Exp(-beta*((emc-E_System)-mu(iAtom)))))

     else if(ipass.eq.0) then
         prob = 0
     end if
         
         call Random_Number(prob2)
         
         !Apply Acceptance Updates if Accepted
         !*************************************
         if (prob.lt.prob2) then
             nfail = nfail + 1
             
             open (15, file='faillog', access='append')
             write (15,'(4i7,f14.8)')nfail,imc,nMCatoms_System(:), E_System 
             close (15)                          
             call system ('cat CONTCAR >> faillog') !log failed geometry step              
             
             call system ('rm CHG* CONTCAR EIGENVAL IBZKPT logfile OSZICAR OUTCAR PCDAT POSCAR REPORT vasprun.xml WAVECAR XDATCAR')  !! clean up files
             
             open (1,file='herefile',access='append')
               call date_and_time(rdate,rtime)
               write (1,'(i10,a30,f10.4,i6,a15)') imc, "Add fail",prob, natoms_System, rtime
             close (1)
         else
             niterMC = niterMC + 1            
             !Update Phase Information

             E_System = emc                       !! Store new E/geometric informations (accept the change)             
             cmc_System = cmc                 
             tfmc_System = tfmc
             natoms1_System = natoms1
             natoms2_System = natoms2
             natoms_System = natoms   
             nMCatoms_System(iAtom) = nMCatoms_System(iAtom) + 1             

             !write updated log files
             open (17, file='MClog', access='append')
                 write (17,'(4i7,f14.8)')niterMC,imc,nMCatoms_System(:), E_System 
             close (17)

             open (16, file='geolog', access='append')
             write (16,'(4i7,f14.8)')niterMC,imc,nMCatoms_System(:), E_System 
             close (16)                          
             call system ('cat CONTCAR >> geolog') !log accepted geometry step
             
             call system ('rm CHG* CONTCAR EIGENVAL IBZKPT logfile OSZICAR OUTCAR PCDAT POSCAR REPORT vasprun.xml WAVECAR XDATCAR')  !! clean up files
             
             open (1,file='herefile',access='append')
               call date_and_time(rdate,rtime)
               write (1,'(i10,a30,f10.4,i6,a15)') imc, "Add pass",prob, natoms_System, rtime
             close (1)
         end if
!*  End Add Atom Block   
!*****************************************
     end if

    
!**********************************************
end do !END of MC Loop
stop
end !END of GEMC Program
!**********************************************



!SUBROUTINES: RunVASP, GetRandomCoords, ReadCoordinates, WriteCoordinates, Writemagmom, init_random_seed
!********************************************************************** 

   subroutine RunVASP(natoms) 
! Runs VASP
! Use Ccount to determine atom types
!**********************************************************************************
     implicit none
     integer, parameter :: KREAL = kind(0.d0)
     integer, intent(in) :: natoms
     integer :: Ccount
     character*5 :: filename
     character*11 :: rtime, rdate
     
     !!INCAR the MAGMOM should be updated with Writmagmom
     !call system ('sed -i "s/NSW.*/NSW   = 20/g" INCAR')           !!  NSW value is preset by tests
     
     !!KPOINTS is unchanged
     
     !!POSCAR is supplied from subroutine WriteCoordinates
     
     !!get updated POTCAR
      !call system ('grep C POSCAR  |wc -l>count')
      !open (2, file='count', status='unknown')
      !     read (2,'(i9)')Ccount       !Ccount indicating if there are C in the system (the head line of POSCAR is set as "C Fe" so Ccount is at least 1)
      !close (2)
      !call system ('rm count')
      !
      !
      !if ( Ccount> 1) then
      !    call system ('cp ~/potcars/Fe-C/POTCAR ./POTCAR')
      !else
      !    call system ('cp ~/potcars/Fe/POTCAR ./POTCAR')
      !end if
      !!!! Use a POSCAR of Fe-C (there are always Fe)
      
     call system ('srun -n 12 vasp_std')        ! run vasp_std with new geometry... mpi-parallel

     
     open (1,file='herefile',access='append')
          call date_and_time(rdate,rtime)
          write (1,'(a50,i6,a15)') "RunVasp done", natoms, rtime
     close (1)
     return
     
   end

!********************************************************************** 

   subroutine GetRandomCoords(natoms,cmc,Rmin,Rmax,xaxismc,yaxismc,zaxismc,nmcmax,randCoord,iz_place,zmin,zmax, iAtom, ipass)
! selects random coordiantes and checks distances to all atoms and images in the system 
!**********************************************************************
      implicit none
      integer, parameter :: KREAL = kind(0.d0)
      real(KREAL),intent(in) :: xaxismc(3), yaxismc(3), zaxismc(3), cmc(nmcmax,3), Rmin(2), Rmax(2), zmin, zmax
      integer, intent(in) :: natoms, nmcmax, iAtom 
      real(KREAL),intent(out):: randCoord(3)
      real(KREAL)  :: rand, distmin, dista, dist(3), distminCN
      integer :: imagetest_x, imagetest_y, imagetest_z, icheck, iz_place, ipass
      character*11 :: rtime, rdate
      
      
      distmin = 1000.0 !set high distmin value to ensure while loop is enetered
         
         distmin = 1000.0                 
         
         call Random_Number(rand)                          !select random coordinates
         randCoord(1) = rand
         call Random_Number(rand)
         randCoord(2) = rand
         call Random_Number(rand)
         if (iz_place.eq.1) then
          randCoord(3) = zmin+(zmax-zmin)*rand
         else
          randCoord(3) = rand
         end if

         
         do imagetest_x=-1,1            !calculate distance to all periodic images
             do imagetest_y=-1,1
                 do imagetest_z=-1,1 
                     do icheck=1,natoms       !calculate distance to all atoms in the system 
                         !check distances to -x,x,+x images, save lowest
                         dist(1) = abs((cmc(icheck,1)+imagetest_x-randCoord(1))*xaxismc(1))       !assume orthogonal cell
                         dist(2) = abs((cmc(icheck,2)+imagetest_y-randCoord(2))*yaxismc(2))
                         dist(3) = abs((cmc(icheck,3)+imagetest_z-randCoord(3))*zaxismc(3))
                         !save shortest distance
                         dista = sqrt(dot_product(dist, dist))
                         if (dista.lt.distmin) then
                              distmin = dista
                         end if
                     end do    
                 end do
             end do        
         end do    
      
      
      if ((distmin.gt.Rmin(iAtom)) .and. (distmin.lt.Rmax(iAtom))) then
       ipass = 1
      else
       ipass = 0
      end if 
      
      open (1,file='herefile',access='append')
          call date_and_time(rdate,rtime)
          write (1,'(a23,4f6.2,a15)') "GetRandomCoords done", distmin, randCoord(1),randCoord(2),randCoord(3),rtime
      close (1)
      
      return
      
  end
!**********************************************************************
!********************************************************************** 

subroutine ReadCoordinates(emc,cmc,tfmc,natoms1,natoms2,natoms,xaxismc,yaxismc,zaxismc)
! reads in 'CONTCAR' data
! Use Ccount from CONTCAR to determine atom types
!********************************************************************** 
      implicit none
      integer, parameter :: KREAL = kind(0.d0)
      integer :: natoms, natoms1, natoms2, i1, Ccount
      real(KREAL) :: emc, ndumb
      integer,parameter :: nmcmax=10000
      real(KREAL) :: cmc(nmcmax,3), xaxismc(3), yaxismc(3), zaxismc(3)
      character*1 :: tfmc(nmcmax,3)
      character*2 :: qatom(2), qdumb1, qdumb2, qdumb3, qdumb4
      character*9 :: dumb1, dumb2, dumb
      character*11 :: rtime, rdate
      
      call system("grep 'energy  without' OUTCAR |tail -1 >energytemp")
      call system("sed -i 's/energy  without entropy=/Fe/g' energytemp")
      call system("sed -i 's/energy(sigma->0) =/C/g' energytemp")
      
      open (96, file='energytemp', status='old')
      read (96,*)qdumb1, ndumb, qdumb2, emc
      close (96)
      
      call system('rm energytemp')
      
      call system ('grep C CONTCAR  |wc -l>count')   !!determine if there is C (Ccount==2 means there is C)
      open (2, file='count', status='old')
           read (2,'(i9)')Ccount
      close (2)
      call system ('rm count')
      
      open (97, file='CONTCAR', status='old')
      read (97,*)qdumb3,qdumb4
      read (97,*)ndumb
      read (97,*)xaxismc(1), yaxismc(1), zaxismc(1)
      read (97,*)xaxismc(2), yaxismc(2), zaxismc(2)
      read (97,*)xaxismc(3), yaxismc(3), zaxismc(3)
      
      if ( Ccount > 1) then    
          read (97,*)qatom(1), qatom(2)
          read (97,*)natoms1, natoms2
      else
          read (97,*)qatom(1)
          read (97,*)natoms1
          natoms2 = 0
      end if
      natoms = natoms1 + natoms2
      
      read (97,*)dumb1, dumb2
      read (97,*)dumb
      
      do i1=1,natoms
         read (97,*)cmc(i1,1),cmc(i1,2),cmc(i1,3),tfmc(i1,1),tfmc(i1,2),tfmc(i1,3)
      end do 
      close(97)

     
      open (1,file='herefile',access='append')
          call date_and_time(rdate,rtime)
          write (1,'(a50,i6,a15)') "ReadCoordinates done", natoms, rtime
      close (1)

      return
   
    end subroutine ReadCoordinates
!**********************************************************************
!********************************************************************** 

  
    subroutine Writecoordinates(xaxismc,yaxismc,zaxismc,cmc,tfmc,natoms1,natoms2,natoms)
! writes geometry data in POSCAR format
! Use natoms2 to determine atom types
!********************************************************************** 
      implicit none
      integer, parameter :: KREAL = kind(0.d0)
      integer :: natoms, natoms1, natoms2, i1
      integer,parameter :: nmcmax=10000
      real(KREAL) :: cmc(nmcmax,3), xaxismc(3), yaxismc(3), zaxismc(3), ndumb
      character*1 :: tfmc(nmcmax,3)
      character*2 :: qatom(2), qdumb1, qdumb2, qdumb3, qdumb4
      character*9 :: dumb1, dumb2, dumb
      character*11 :: rtime, rdate
      
      open (98, file='POSCAR', status='new')
      write (98,'(a5)')"Fe  C"
      write (98,'(a19)')"   1.00000000000000"
      write (98,'(3f20.16)')xaxismc(1), yaxismc(1), zaxismc(1)
      write (98,'(3f20.16)')xaxismc(2), yaxismc(2), zaxismc(2)
      write (98,'(3f20.16)')xaxismc(3), yaxismc(3), zaxismc(3)   
      
      if ( natoms2 > 0) then                           !!use natoms2 to determine if there is C, since C number might change before write coordinates 
          write (98,'(2a6)')"Fe","C"
          write (98,'(2i6)')natoms1, natoms2
      else
          write (98,'(a6)')"Fe"
          write (98,'(i6)')natoms1
      end if
      
      write (98,'(a18)')"Selective dynamics"
      write (98,'(a6)')"Direct"
      
      do i1=1,natoms
         write (98,'(3f23.16, 3a4)')cmc(i1,1),cmc(i1,2),cmc(i1,3),tfmc(i1,1),tfmc(i1,2),tfmc(i1,3)
      end do
      close (98)
      
      open (1,file='herefile',access='append')
          call date_and_time(rdate,rtime)
          write (1,'(a50,i6,a15)') "WriteCoordinates done", natoms, rtime
      close (1)

      return
   
    end subroutine Writecoordinates
!*********************************************************************     
!**********************************************************************



    subroutine Writemagmom(natoms1, natoms2)
! writes MAGMOM tag in INCAR based on number of atoms
! Use before run VASP
! Give MAGMOM tag at the last line of INCAR
! Use natoms2 (C) to determine atom types
!********************************************************************** 

      implicit none
      integer, parameter :: KREAL = kind(0.d0)
      integer :: natoms1, natoms2    !1 refer to C, 2 refer to Fe
      real(KREAL) :: magC, magFe

      magC = 0
      magFe = 2.60
      
      if (natoms2 > 0) then
           open (76,file='tmpmagmom',access='append')
               write (76,'(a9,i4,a1,f4.2,i4,a1,f4.2)') "MAGMOM = ", natoms1, "*", magFe, natoms2, "*", magC
           close (76)
      else
           open (76,file='tmpmagmom',access='append')
               write (76,'(a9,i4,a1,f4.2)') "MAGMOM  = ", natoms1, "*", magFe
           close (76)
      end if

      call system ("sed '/MAGMOM/d' INCAR > tmpincar")
      call system ('rm INCAR')
      call system ('cat tmpincar tmpmagmom > INCAR')
      call system ('rm tmpincar tmpmagmom')
      
    end subroutine Writemagmom
!*********************************************************************     
!**********************************************************************



  subroutine init_random_seed() 
! generate a random seed for initializing the random number generator
!********************************************************************** 
   implicit none
   integer, parameter :: KREAL = kind(0.d0)
   integer, allocatable :: seed(:)
   integer :: i, n, un, istat, dt(8), pid, t(2), s, getpid
   integer(8) :: count, tms
   real(8) :: rtemp
   real(KREAL)  :: rand
 
   call random_seed(size = n)
   allocate(seed(n))
      call system_clock(count)
      if (count /= 0) then
         t = transfer(count, t)
      else
         call date_and_time(values=dt)
         tms = (dt(1) - 1970) * 365_8 * 24 * 60 * 60 * 1000 &
              + dt(2) * 31_8 * 24 * 60 * 60 * 1000 &
              + dt(3) * 24 * 60 * 60 * 60 * 1000 &
              + dt(5) * 60 * 60 * 1000 &
              + dt(6) * 60 * 1000 + dt(7) * 1000 &
              + dt(8)
         t = transfer(tms, t)
      end if
      s = ieor(t(1), t(2))
      pid = getpid() + 1099279 ! Add a prime
      s = ieor(s, pid)
      if (n >= 3) then
         seed(1) = t(1) + 36269
         seed(2) = t(2) + 72551
         seed(3) = pid
         if (n > 3) then
            seed(4:) = s + 37 * (/ (i, i = 0, n - 4) /)
         end if
      else
         seed = s + 37 * (/ (i, i = 0, n - 1 ) /)
      end if
   call random_seed(put=seed)
   rtemp = rand(seed(1))
   end subroutine init_random_seed
