program nuclmass
implicit none
integer,parameter  :: dprec=8

!Values of fundamental constants and conversion factors
real(dprec),parameter ::   electron_mass=5.48579909065E-04_dprec  !electron mass in atomic mass units, u
real(dprec),parameter :: d_electron_mass=0.00000000016E-04_dprec
real(dprec),parameter ::   hartree2u=2.92126232205E-08_dprec !conversion factor for hartree to atomic mass unit, u
real(dprec),parameter :: d_hartree2u=0.00000000088E-08_dprec
real(dprec),parameter ::   ev2u=1.07354410233E-09_dprec !conversion factor  for eV to atomic mass unit, u
real(dprec),parameter :: d_ev2u=0.00000000032E-09_dprec !conversion factor  for eV to atomic mass unit, u


integer              :: Z
real(dprec)         :: atomic_mass,d_atomic_mass
real(dprec)         :: electron_binding_energy,d_electron_binding_energy
real(dprec)         :: nuclear_mass,d_nuclear_mass


!To get the nuclear mass accurate to 4 decimal points, the total electron binding
!energy only needs to be accurate to 0.1 hartree or even less.
!Ideally, the total electron binding energy should include finite-nuclear-mass and
!relativistic effects. But, again, in reality one might ignore these small effects
!without consequences.

write(*,'(a)') 'The calculation of masses of light nuclei (in the electron masses) including'
write(*,'(a)') 'error estimates. The values of the fundamental constants were taken from' 
write(*,'(a)') 'CODATA 2018 recommended values [Tiesinga et al., Rev. Mod. Phys 93, 025010 (2021)]'
write(*,'(a)') '  https://dx.doi.org/10.1103/RevModPhys.93.025010'
write(*,'(a)') '  https://physics.nist.gov/cuu/constants'
write(*,'(a)') 'Atomic weights were taken from AME2020 atomic mass evaluation'
write(*,'(a)') '[Wang et al., Chin. Physics C 45, 030003 (2021)]'
write(*,'(a)') '  https://dx.doi.org/10.1088/1674-1137/abddaf'
write(*,'(a)') 'Most common isotopic compositions are avaialble here'
write(*,'(a)') '  https://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl'

Z=1
  atomic_mass=1007825.031898E-06_dprec
d_atomic_mass=0.000014E-06_dprec
  electron_binding_energy=0.49973449614E0_dprec*hartree2u
d_electron_binding_energy=0.000002E0_dprec*hartree2u
nuclear_mass=(atomic_mass+electron_binding_energy)/electron_mass-Z
d_nuclear_mass=(d_atomic_mass+d_electron_binding_energy)/electron_mass+ &
               (atomic_mass+electron_binding_energy)*d_electron_mass/electron_mass**2
write(*,*)
write(*,'(a,f14.8)') 'M(1H) =   ',nuclear_mass
write(*,'(a,f10.8)') '             ±',d_nuclear_mass
write(*,'(a)') 'M(1H) =    1836.15267343(11)  [CODATA 2018]'

Z=1
  atomic_mass=2014101.777844E-06_dprec
d_atomic_mass=0.000015E-06_dprec
  electron_binding_energy=0.49987047167E0_dprec*hartree2u
d_electron_binding_energy=0.000002E0_dprec*hartree2u
nuclear_mass=(atomic_mass+electron_binding_energy)/electron_mass-Z
d_nuclear_mass=(d_atomic_mass+d_electron_binding_energy)/electron_mass+ &
               (atomic_mass+electron_binding_energy)*d_electron_mass/electron_mass**2
write(*,*)
write(*,'(a,f14.8)') 'M(2H) =   ',nuclear_mass
write(*,'(a,f10.8)') '             ±',d_nuclear_mass
write(*,'(a)') 'M(2H) =    3670.48296788(13)  [CODATA 2018]'

Z=1
  atomic_mass=3016049.281320E-06_dprec
d_atomic_mass=0.000081E-06_dprec
  electron_binding_energy=0.49991571296E0_dprec*hartree2u
d_electron_binding_energy=0.000002E0_dprec*hartree2u
nuclear_mass=(atomic_mass+electron_binding_energy)/electron_mass-Z
d_nuclear_mass=(d_atomic_mass+d_electron_binding_energy)/electron_mass+ &
               (atomic_mass+electron_binding_energy)*d_electron_mass/electron_mass**2
write(*,*)
write(*,'(a,f14.8)') 'M(3H) =   ',nuclear_mass
write(*,'(a,f10.8)') '             ±',d_nuclear_mass
write(*,'(a)') 'M(3H) =    5496.92153573(27)  [CODATA 2018]'

Z=2
  atomic_mass=3016029.321967E-06_dprec
d_atomic_mass=0.000060E-06_dprec
  electron_binding_energy=2.903271162695E0_dprec*hartree2u !includes rel. corr., source: M. Stanke et al. JCP 126, 194312 (2007)
d_electron_binding_energy=0.00002E0_dprec*hartree2u
nuclear_mass=(atomic_mass+electron_binding_energy)/electron_mass-Z
d_nuclear_mass=(d_atomic_mass+d_electron_binding_energy)/electron_mass+ &
               (atomic_mass+electron_binding_energy)*d_electron_mass/electron_mass**2
write(*,*)
write(*,'(a,f14.8)') 'M(3He) =  ',nuclear_mass
write(*,'(a,f10.8)') '             ±',d_nuclear_mass
write(*,'(a)') 'M(3He) =   5495.88528007(24)  [CODATA 2018]' 

Z=2
  atomic_mass=4002603.254130E-06_dprec
d_atomic_mass=0.000158E-06_dprec
  electron_binding_energy=2.903408504575E0_dprec*hartree2u !includes rel. corr., source: M. Stanke et al. JCP 126, 194312 (2007)
d_electron_binding_energy=0.00002E0_dprec*hartree2u
nuclear_mass=(atomic_mass+electron_binding_energy)/electron_mass-Z
d_nuclear_mass=(d_atomic_mass+d_electron_binding_energy)/electron_mass+ &
               (atomic_mass+electron_binding_energy)*d_electron_mass/electron_mass**2
write(*,*)
write(*,'(a,f14.8)') 'M(4He) =  ',nuclear_mass
write(*,'(a,f10.8)') '             ±',d_nuclear_mass
write(*,'(a)') 'M(4He) =   7294.29954142(24)  [CODATA 2018]'

Z=3
  atomic_mass=6015122.8874E-06_dprec
d_atomic_mass=0.0015E-06_dprec
  electron_binding_energy=7.47799240823E0_dprec*hartree2u !includes rel. corr., source:  M. Stanke et. al., PRA 78, 052507 (2008)
d_electron_binding_energy=0.0001E0_dprec*hartree2u
nuclear_mass=(atomic_mass+electron_binding_energy)/electron_mass-Z
d_nuclear_mass=(d_atomic_mass+d_electron_binding_energy)/electron_mass+ &
               (atomic_mass+electron_binding_energy)*d_electron_mass/electron_mass**2
write(*,*)
write(*,'(a,f14.8)') 'M(6Li) =  ',nuclear_mass
write(*,'(a,f10.8)') '             ±',d_nuclear_mass

Z=3
  atomic_mass=7016003.434E-06_dprec
d_atomic_mass=0.004E-06_dprec
  electron_binding_energy=7.47809366029E0_dprec*hartree2u !includes rel. corr., source:  M. Stanke et. al., PRA 78, 052507 (2008)
d_electron_binding_energy=0.0001E0_dprec*hartree2u
nuclear_mass=(atomic_mass+electron_binding_energy)/electron_mass-Z
d_nuclear_mass=(d_atomic_mass+d_electron_binding_energy)/electron_mass+ &
               (atomic_mass+electron_binding_energy)*d_electron_mass/electron_mass**2
write(*,*)
write(*,'(a,f14.8)') 'M(7Li) =  ',nuclear_mass
write(*,'(a,f10.8)') '             ±',d_nuclear_mass

Z=4
  atomic_mass=9012183.06E-06_dprec
d_atomic_mass=0.08E-06_dprec
  electron_binding_energy=14.668795768E0_dprec*hartree2u !includes rel. corr., source:  M. Stanke et. al., PRA 80, 022514 (2009)
d_electron_binding_energy= 0.0005E0_dprec*hartree2u
nuclear_mass=(atomic_mass+electron_binding_energy)/electron_mass-Z
d_nuclear_mass=(d_atomic_mass+d_electron_binding_energy)/electron_mass+ &
               (atomic_mass+electron_binding_energy)*d_electron_mass/electron_mass**2
write(*,*)
write(*,'(a,f14.8)') 'M(9Be) =  ',nuclear_mass
write(*,'(a,f10.8)') '             ±',d_nuclear_mass

Z=5
  atomic_mass=10012936.862E-06_dprec
d_atomic_mass=0.016E-06_dprec
  electron_binding_energy=(24.65250024E0_dprec+0.005E0_dprec)*hartree2u !nonrel. energy:  S. Bubin and L. Adamowicz, PRA 83, 022505 (2011); rel. corr. P. Seth et. al., JCP 134, 084105 (2011)
d_electron_binding_energy= 0.001E0_dprec*hartree2u
nuclear_mass=(atomic_mass+electron_binding_energy)/electron_mass-Z
d_nuclear_mass=(d_atomic_mass+d_electron_binding_energy)/electron_mass+ &
               (atomic_mass+electron_binding_energy)*d_electron_mass/electron_mass**2
write(*,*)
write(*,'(a,f14.8)') 'M(10B) =  ',nuclear_mass
write(*,'(a,f10.8)') '             ±',d_nuclear_mass

Z=5
  atomic_mass=11009305.167E-06_dprec
d_atomic_mass=0.013E-06_dprec
  electron_binding_energy=(24.65262387E0_dprec+0.005E0_dprec)*hartree2u !nonrel. energy:  S. Bubin and L. Adamowicz, PRA 83, 022505 (2011); rel. corr. P. Seth et. al., JCP 134, 084105 (2011)
d_electron_binding_energy= 0.001E0_dprec*hartree2u
nuclear_mass=(atomic_mass+electron_binding_energy)/electron_mass-Z
d_nuclear_mass=(d_atomic_mass+d_electron_binding_energy)/electron_mass+ &
               (atomic_mass+electron_binding_energy)*d_electron_mass/electron_mass**2
write(*,*)
write(*,'(a,f14.8)') 'M(11B) =  ',nuclear_mass
write(*,'(a,f10.8)') '             ±',d_nuclear_mass

Z=6
  atomic_mass=12.0E0_dprec
d_atomic_mass=0.000000000E0_dprec
  electron_binding_energy=(37.84446E0_dprec-0.0017166E0_dprec+0.015E0_dprec)*hartree2u
d_electron_binding_energy= 0.005E0_dprec*hartree2u
nuclear_mass=(atomic_mass+electron_binding_energy)/electron_mass-Z
d_nuclear_mass=(d_atomic_mass+d_electron_binding_energy)/electron_mass+ &
               (atomic_mass+electron_binding_energy)*d_electron_mass/electron_mass**2
write(*,*)
write(*,'(a,f14.8)') 'M(12C) =  ',nuclear_mass
write(*,'(a,f10.8)') '             ±',d_nuclear_mass

Z=6
  atomic_mass=13003354.83534E-06_dprec
d_atomic_mass=0.00025E-06_dprec
  electron_binding_energy=(37.84446E0_dprec-0.00158412E0_dprec+0.015E0_dprec)*hartree2u
d_electron_binding_energy= 0.005E0_dprec*hartree2u
nuclear_mass=(atomic_mass+electron_binding_energy)/electron_mass-Z
d_nuclear_mass=(d_atomic_mass+d_electron_binding_energy)/electron_mass+ &
               (atomic_mass+electron_binding_energy)*d_electron_mass/electron_mass**2
write(*,*)
write(*,'(a,f14.8)') 'M(13C) =  ',nuclear_mass
write(*,'(a,f10.8)') '             ±',d_nuclear_mass

Z=6
  atomic_mass=14003241.989E-06_dprec
d_atomic_mass=0.004E-06_dprec
  electron_binding_energy=(37.84446E0_dprec-0.00147098E0_dprec+0.015E0_dprec)*hartree2u
d_electron_binding_energy=0.005E0_dprec*hartree2u
nuclear_mass=(atomic_mass+electron_binding_energy)/electron_mass-Z
d_nuclear_mass=(d_atomic_mass+d_electron_binding_energy)/electron_mass+ &
               (atomic_mass+electron_binding_energy)*d_electron_mass/electron_mass**2
write(*,*)
write(*,'(a,f14.8)') 'M(14C) =  ',nuclear_mass
write(*,'(a,f10.8)') '             ±',d_nuclear_mass

Z=7
  atomic_mass=14003074.00425E-06_dprec
d_atomic_mass=0.00024E-06_dprec
  electron_binding_energy=(54.58867E0_dprec+0.03E0_dprec)*hartree2u
d_electron_binding_energy= 0.01E0_dprec*hartree2u
nuclear_mass=(atomic_mass+electron_binding_energy)/electron_mass-Z
d_nuclear_mass=(d_atomic_mass+d_electron_binding_energy)/electron_mass+ &
               (atomic_mass+electron_binding_energy)*d_electron_mass/electron_mass**2
write(*,*)
write(*,'(a,f14.8)') 'M(14N) =  ',nuclear_mass
write(*,'(a,f10.8)') '             ±',d_nuclear_mass

Z=7
  atomic_mass=15000108.8983E-06_dprec
d_atomic_mass=0.0006E-06_dprec
  electron_binding_energy=(54.58867E0_dprec+0.03E0_dprec)*hartree2u
d_electron_binding_energy= 0.01E0_dprec*hartree2u
nuclear_mass=(atomic_mass+electron_binding_energy)/electron_mass-Z
d_nuclear_mass=(d_atomic_mass+d_electron_binding_energy)/electron_mass+ &
               (atomic_mass+electron_binding_energy)*d_electron_mass/electron_mass**2
write(*,*)
write(*,'(a,f14.8)') 'M(15N) =  ',nuclear_mass
write(*,'(a,f10.8)') '             ±',d_nuclear_mass

Z=8
  atomic_mass=15994914.6193E-06_dprec
d_atomic_mass=0.0003E-06_dprec
  electron_binding_energy=(75.0673E0_dprec+0.053E0_dprec)*hartree2u
d_electron_binding_energy= 0.02E0_dprec*hartree2u
nuclear_mass=(atomic_mass+electron_binding_energy)/electron_mass-Z
d_nuclear_mass=(d_atomic_mass+d_electron_binding_energy)/electron_mass+ &
               (atomic_mass+electron_binding_energy)*d_electron_mass/electron_mass**2
write(*,*)
write(*,'(a,f14.8)') 'M(16O) =  ',nuclear_mass
write(*,'(a,f10.8)') '             ±',d_nuclear_mass

Z=8
  atomic_mass=16999131.7560E-06_dprec
d_atomic_mass=0.0007E-06_dprec
  electron_binding_energy=(75.0673E0_dprec+0.053E0_dprec)*hartree2u
d_electron_binding_energy= 0.02E0_dprec*hartree2u
nuclear_mass=(atomic_mass+electron_binding_energy)/electron_mass-Z
d_nuclear_mass=(d_atomic_mass+d_electron_binding_energy)/electron_mass+ &
               (atomic_mass+electron_binding_energy)*d_electron_mass/electron_mass**2
write(*,*)
write(*,'(a,f14.8)') 'M(17O) =  ',nuclear_mass
write(*,'(a,f10.8)') '             ±',d_nuclear_mass

Z=8
  atomic_mass=17999159.6121E-06_dprec
d_atomic_mass=0.0007E-06_dprec
  electron_binding_energy=(75.0673E0_dprec+0.053E0_dprec)*hartree2u
d_electron_binding_energy= 0.02E0_dprec*hartree2u
nuclear_mass=(atomic_mass+electron_binding_energy)/electron_mass-Z
d_nuclear_mass=(d_atomic_mass+d_electron_binding_energy)/electron_mass+ &
               (atomic_mass+electron_binding_energy)*d_electron_mass/electron_mass**2
write(*,*)
write(*,'(a,f14.8)') 'M(18O) =  ',nuclear_mass
write(*,'(a,f10.8)') '             ±',d_nuclear_mass

end program nuclmass