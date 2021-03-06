C.
      SUBROUTINE uhinit
C.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C         Defines HBOOK histogram/scatterplot definitions              C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C.
      IMPLICIT none
C.
      include 'gcflag.inc'          !geant
      include 'gcunit.inc'          !geant
C.
      include 'geometry.inc'        !local
      include 'u_geom.inc'          !local
      include 'gbox_info.inc'       !local
C.
      include 'uggeom.inc'          !local
      include 'higamcoinc.inc'      !local
      include 'res.inc'             !local
      include 'beamcom.inc'         !local
      include 'rescom.inc'          !local
      include 'history.inc'
      include 'gammahit.inc'

C.
      REAL sig, lm1, lm2, angdist, weight, ener
      EXTERNAL sig, angdist !, alphapos

C.
      INTEGER i, n, istat, nstrip, nlines
      CHARACTER*32 wfile
      character*11 strip
      character*2 num(16)
      character*8 chtags(NVAR10)
      data chtags/ 'GML0', 'GML1','GML2','GML3','GML4','GML5','GML6',
     &   'GML7','GML9','GML10','GML11','GML12','GML13','GML14','GML15',
     &  'GML16','GML18','GML19','GML20','GML21','GML22','GML23','GML24',
     &  'GML25', 'GML27','GML28','GML29',
     &  'HML0','HML1','HML2','HML3','HML4','HML5','HML6',
     &  'HML7','HML8','HML9'/
      data num/'01','02','03','04','05','06','07','08','09',
     &         '10','11','12','13','14','15','16'/
C.
C -->   Open a HBOOK direct access file
C.
      CALL namfil('dragon',idrun,'.hbook',wfile)
C
      CALL HROPEN(lunits(4),'HBOOK',wfile,'N',1024,istat)
C.
      If(istat.ne.0)then
        WRITE(lout,*)' Error: Bad return from HROPEN! '
        STOP
      Endif
C.
C -->   Initialize user HBOOK histograms and scatterplots
C.
      n = 0
C.
      CALL hbook1(n+ 1,' Initial - x - ',100,-2.0,2.0,0.0)
      CALL hbook1(n+ 2,' Initial - y - ',100,-2.0,2.0,0.0)
C.
      CALL hbook1(n+ 3,' Initial - dx - ',100,-100.0,100.0,0.0)
      CALL hbook1(n+ 4,' Initial - dy - ',100,-100.0,100.0,0.0)
C.
      CALL hbook1(n+ 5,' IniFin - x - Stops ',100,-2.0,2.0,0.0)
      CALL hbook1(n+ 6,' IniFin - y - Stops ',100,-2.0,2.0,0.0)
C.
      CALL hbook1(n+ 7,' IniFin - dx - Stops',100,-100.0,100.0,0.0)
      CALL hbook1(n+ 8,' IniFin - dy - Stops',100,-100.0,100.0,0.0)
C.
      CALL hbook1(n+ 9,' Initial Momentum ',100,-5.0,5.0,0.0)
      CALL hbook1(n+10,' Momentum spread (%) ',100,-5.0,5.0,0.0)
C.
      CALL hbook1(n+11,' Final - x - ',16, -2.4, 2.4,0.0)
      CALL hbook1(n+12,' Final - y - ',16, -2.4, 2.4,0.0)

      CALL hbook1(n+13,' Final - dx - ',100,-100.0,100.0,0.0)
      CALL hbook1(n+14,' Final - dy - ',100,-50.0,50.0,0.0)
C.
      CALL hbook1(n+15,' Final Energy ',4000,0.,20.0,0.0)
C.
      CALL hbook1(n+16,' Stop Length (cm) ',2000,0.0,2000.0,0.0)
      CALL hbook2(n+17, ' X vs Stop Length (cm) ',50,-20.,20.,
     &            200,0.,2000.,0.) 
C.
      n = 20
C.
      CALL hbook1( n+1,' True photon energy     ', 200, 0.,  20., 0.)
      CALL hbook1( n+2,' True photon pol. angle ', 200, 0., 200., 0.)
      CALL hbook1( n+3,' Photon conv. module   ',   29, 1.,  30., 0.)
      CALL hbook2( n+4,' Photon creation time vs z_react', 
     &    600,0.,300., 300,-15.,15.,0.)
      CALL hbook2( n+5,' Photon detection time vs z_react', 
     &    600,0.,300., 300,-15.,15.,0.)
C.
      CALL hbook1(n+11,' No. of Modules hit  ',  10,  0., 10., 0.)
      CALL hbook1(n+12,' Total energy dep.   ', 200,  0., 20., 0.)
      CALL hbook1(n+13,' x-coordinates of hit',  60, -15., 15., 0.)
      CALL hbook1(n+14,' y-coordinates of hit',  80, -20., 20., 0.)
      CALL hbook1(n+15,' z-coordinates of hit', 100, -20., 20., 0.)
      CALL hbook1(n+16,' Energy dep. in module     ', 200, 0., 20., 0.)
      CALL hbook1(n+17,' Energy dep. in 1. module  ', 200, 0., 20., 0.)
      CALL hbook1(n+18,' Energy dep. in 2. module  ', 200, 0., 20., 0.)
      CALL hbook1(n+19,' Energy dep. in 3. module  ', 200, 0., 20., 0.)
      CALL hbook1(n+20,' Energy dep. in 4. module  ', 200, 0., 20., 0.)
C.
      CALL hbook1(n+21,' Energy dep. in crystal    ', 200, 0., 20., 0.)
C.
      CALL hbook1(n+26,' True conversion z ', 100, -20., 20., 0.)
      CALL hbook1(n+27,' Energy weighted z ', 100, -20., 20., 0.)
C.
      CALL hbook1(n+28,' Distance: conv. and max-energy dep.    (xy)  '
     &               , 100, 0., 1., 0.)
      CALL hbook1(n+29,' Distance: conv. and max-energy dep. (  xyz)  '
     &               , 100, 0., 1., 0.)
      CALL hbook1(n+30,' Distance: PMT and max-energy dep. (xy)  '
     &               , 100, 0., 20., 0.)
C.
      CALL hbook1(n+31,' Number of photons detected in PMT '
     &               , 200, 10., 10000., 0.)
      CALL hbook1(n+32,' Photons in 1. PMT '
     &               , 2000, 0., 10000., 0.)
      CALL hbook1(n+33,' Photons in 2. PMT '
     &               ,  100, 0.,  1000., 0.)
      CALL hbook1(n+34,' Photons in 3. PMT '
     &               ,  100, 0.,  1000., 0.)
      CALL hbook1(n+35,' Photons in 4. PMT '
     &               ,  100, 0.,  1000., 0.)
C.
      CALL hbook1(n+36,' Number of PMTs hit above threshold'
     &               ,   20, 0.,    20., 0.)
      CALL hbook1(n+37,' Total Number of photons det. in PMTs > thrsld '
     &               , 2000, 0., 10000., 0.)
C.
      CALL hbook1(n+40,' Reconstructed x-position '
     &               , 120, -30., 30., 0.)
      CALL hbook1(n+41,' Reconstructed y-position '
     &               , 120, -30., 30., 0.)
      CALL hbook1(n+42,' Reconstructed z-position '
     &               , 120, -30., 30., 0.)
C.
      CALL hbook1(n+50,' Number of photons max. '
     &               , 100, 0., 15000., 0.)
      CALL hbook1(n+51,' Number of photons generated '
     &               , 100, 0.,  5000., 0.)
      CALL hbook1(n+52,' Number of photons lost LABS '
     &               , 100, 0.,  5000., 0.)
      CALL hbook1(n+53,' Number of photons lost REFL '
     &               , 100, 0.,  5000., 0.)
      CALL hbook1(n+54,' Number of photons lost ds < e '
     &               , 100, 0.,  1000., 0.)
      CALL hbook1(n+55,' Number of photons lost N > 1000 '
     &               , 100, 0.,  1000., 0.)
      CALL hbook1(n+56,' Number of photons unable to reflect '
     &               , 100, 0.,  1000., 0.)
      CALL hbook1(n+57,' Number of photons with error from GLISUR '
     &               , 100, 0.,  1000., 0.)
C.
      CALL hbook1(n+61,' Number of steps taken to PMT '
     &               , 100, 0., 200., 0.)
      CALL hbook1(n+62,' Total track length to PMT '
     &               , 100, 0., 100., 0.)
C.
      CALL hbook1(n+70,' Number of photon clusters ',10, 0., 10., 0.)
      CALL hbook1(n+71,' Energy of Cluster #1 ', 200, 0., 20., 0.)
      CALL hbook1(n+72,' Energy of Cluster #2 ', 200, 0., 20., 0.)
      CALL hbook1(n+73,' Energy of Cluster #3 ', 200, 0., 20., 0.)
      CALL hbook1(n+74,' Energy difference #1 ',  80, -2.,  2., 0.)
      CALL hbook1(n+75,' Energy difference #2 ',  80, -2.,  2., 0.)
      CALL hbook1(n+76,' Energy difference #3 ',  80, -2.,  2., 0.)
      CALL hbook1(n+77,' Dir. diff. [deg] #1 ', 90, 0., 90., 0.)
      CALL hbook1(n+78,' Dir. diff. [deg] #2 ', 90, 0., 90., 0.)
      CALL hbook1(n+79,' Dir. diff. [deg] #3 ', 90, 0., 90., 0.)
C.
      n = 100
C.
      CALL hbook2( n+1,' Initial -  y - vs -  x - ',
     &                100,-2.5,2.5,100,-2.5,2.5,0.0)
      CALL hbook2( n+2,' Initial - dy - vs - dx - ',
     &                100,-100.0,100.0,100,-100.0,100.0,0.0)
      CALL hbook2( n+3,' IniFin -  y - vs -  x - ',
     &                100,-2.0,2.0,100,-2.0,2.0,0.0)
      CALL hbook2( n+4,' IniFin - dy - vs - dx - ',
     &                100,-100.0,100.0,100,-100.0,100.0,0.0)
      CALL hbook2( n+5,' Initial - dx - vs - x -',
     &                100,-2.0,2.0,100,-100.0,100.0,0.0)
      CALL hbook2( n+6,' Initial - dy - vs - y -',
     &                100,-2.0,2.0,100,-100.0,100.0,0.0)
      CALL hbook2( n+7,' IniFin - dx - vs - x -',
     &                100,-2.0,2.0,100,-100.0,100.0,0.0)
      CALL hbook2( n+8,' IniFin - dy - vs - y -',
     &                100,-2.0,2.0,100,-100.0,100.0,0.0)
C.
      CALL hbook2(n+11,' Final -  y - vs -  x - ',
     &                16,-2.4,2.4,16,-2.4,2.4,0.0)
      CALL hbook2(n+12,' Final - dy - vs - dx - ',
     &                100,-100.0,100.0,100,-50.0,50.0,0.0)
      CALL hbook2(n+13,' Final - dx - vs - x ',
     &                100,-2.5,2.5,100,-100.0,100.0,0.0)
      CALL hbook2(n+14,' Final - dy - vs - y ',
     &                100,-3.0,3.0,100,-50.0,50.0,0.0)
C.
      CALL hbook2(n+15,' Final theta vs radius ',
     &                100,0.0,2.0,100,0.0,100.0,0.0)
C.
      n = 120
C.
      CALL hbook2( n+1,' True conversion position '
     &               , 60, -15., 15., 80, -20., 20., 0.)
      CALL hbook2( n+2,' Energy weighted position '
     &               , 60, -15., 15., 80, -20., 20., 0.)
      CALL hbook2( n+3,' Reconstructed position   '
     &               , 60, -15., 15., 80, -20., 20., 0.)
      CALL hbook2( n+4,' True conversion xy fngr-coordinates '
     &               , 60, -15., 15., 80, -20., 20., 0.)
      CALL hbook2( n+5,' True conversion zx fngr-coordinates '
     &               , 60, -15., 15., 80, -20., 20., 0.)
      CALL hbook2( n+6,' True conversion zy fngr-coordinates '
     &               , 60, -15., 15., 80, -20., 20., 0.)
C.
      CALL hbook2(n+11,' Max loop =  1 ',29, 1., 30., 100, 0., 10., 0.)  
      CALL hbook2(n+12,' Max loop =  2 ',29, 1., 30., 100, 0., 10., 0.)  
      CALL hbook2(n+13,' Max loop =  3 ',29, 1., 30., 100, 0., 10., 0.)  
      CALL hbook2(n+14,' Max loop =  4 ',29, 1., 30., 100, 0., 10., 0.)  
      CALL hbook2(n+15,' Max loop =  5 ',29, 1., 30., 100, 0., 10., 0.)  
      CALL hbook2(n+16,' Max loop =  6 ',29, 1., 30., 100, 0., 10., 0.)  
      CALL hbook2(n+17,' Max loop =  7 ',29, 1., 30., 100, 0., 10., 0.)  
      CALL hbook2(n+18,' Max loop =  8 ',29, 1., 30., 100, 0., 10., 0.)  
      CALL hbook2(n+19,' Max loop =  9 ',29, 1., 30., 100, 0., 10., 0.)  
      CALL hbook2(n+20,' Max loop = 10 ',29, 1., 30., 100, 0., 10., 0.)
C.
      n = 200
C.
      CALL hbook1(n+ 1,'Z-Stops in all col',1000,-TLrms,TLrms,0.)
      CALL hbook1(n+ 2,'R-Stops in targ entrance col',50,0.,Rrms/2.,0.)
      CALL hbook1(n+ 3,'R-Stops in targ exit     col',50,0.,Rrms/2.,0.)
      CALL hbook1(n+ 4,'TOF to TEND',200,0.,2.,0.)
      CALL hbook1(n+ 5,'Reaction z-pos',1000,-3*targetl,3*targetl,0.)
      CALL hbook1(n+ 6,'Beam Stops in target',1000,-TLrms,TLrms,0.)
      CALL hbook2(n+11,'stop/exit dist',
     &            200,-TLrms,TLrms,20,0.,Rrms/2.,0.)
      CALL hbook2(n+12,'stop dist     ',
     &            200,-TLrms,TLrms,20,0.,Rrms/2.,0.)
      CALL hbook2(n+13, 'Reaction position R-pos vs z-pos ',
     &                   50,-3*targetl,3*targetl,50,0.,1.,0.)
      CALL hbook2(n+14,'Exit spot',
     &            50,-Rrms/2.,Rrms/2.,50,-Rrms/2.,Rrms/2.,0.)
      CALL hbook2(n+15,'cosx vs X',
     &            100,-Rrms/2.,Rrms/2.,100,-.02,.02,0.)
      CALL hbook2(n+16,'cosy vs Y',
     &            100,-Rrms/2.,Rrms/2.,100,-.02,.02,0.)
      CALL hbook2(n+17,'cosx vs Xtarg',
     &            100,-1.,1.,100,-.02,.02,0.)
      CALL hbook2(n+18,'cosy vs Ytarg',
     &            100,-1.,1.,100,-.02,.02,0.)
C.
C.
      CALL hbook1(n+ 21,' Ini - x - Recoils ',100,-2.0,2.0,0.0)
      CALL hbook1(n+ 22,' Ini - y - Recoils ',100,-2.0,2.0,0.0)
C.
      CALL hbook1(n+ 23,' Ini - dx - Recoils',100,-100.0,100.0,0.0)
      CALL hbook1(n+ 24,' Ini - dy - Recoils',100,-100.0,100.0,0.0)
      CALL hbook1(n+ 25,' Momentum spread (%)-Recoils',100,-5.0,5.0,0.0)



c  MT adds histograms
      n = 300
      CALL hbook1(n+0,' Recoil Stopping Length (cm) ',
     &            2000, 0.0, 2000.0, 0.0)
      CALL hbook1(n+1,' Beam Particle Stopping Length (cm) ',
     &            2000, 0.0, 2000.0, 0.0)
C      CALL hbook2(n+10, ' DSSSD Energy (MeV) vs. x-strip position ',
C     &     16, -2.4, 2.4, 4000, 0.0, 20.0, 0.0)
C      CALL hbook2(n+11, ' DSSSD x-strip position vs. Energy (MeV) ',
C     &     4000, 0.0, 20.0, 16, -2.4, 2.4, 0.0) 

C.
C.    Strip detector separate spectra
C.      n = 300
C.      nstrip = 16
C.      do i = 1, nstrip
C.       strip = 'x-strip('//num(i)//')'
C.       CALL hbook1(n+i,strip,16,0.,16.,0.0)
C.       strip = 'y-strip('//num(i)//')'
C.       CALL hbook1(n+nstrip+i,strip,16,0.,16.,0.0)
C.      enddo
C.
C.    Strip detector hit patterns
C.        n = 400
C.        nstrip = 16
C.        CALL hbook1(n+1,'x-strip hit pattern',16,0.,16.,0.0)
C.        CALL hbook1(n+2,'y-strip hit pattern',16,0.,16.,0.0)

C.    user defined angular distribution

          CALL hbfun1(250,'ang. dist.',1000,-1.,1.,angdist)

C.    cross-section function and probability density
        lm1 = (1.-0.005)*beamo*m2/(m1+m2)
        lm2 = (1.+3.*(emax/beamenerg))*beamenerg*m2/(m1+m2)
        print*, m1, m2, beamenerg, beamo, el
        CALL hbfun1(500,'capture cross-section',1000.,lm1,
     +                   lm2,sig)
        CALL hcopy(500,501,'')

C.    CM energy distribution
       CALL hbook1(502,'CM energy distribution',1000,beamo*(1.-0.01),
     +              beamenerg*(1.+0.01),0.0)

C.    Beam caught after D1, scaler
       CALL hbook1(503,'Caught beam',1,1.,2.,0.0)

C.    Recoils which don't make it to ENDV
       CALL hbook2(504,'Stopped recoil pos.',100,-2000.,1000.,100,
     +             -1000.,100,0.0)
C.
C.    A quick way to keep track of how many recoils hit the end detector
      CALL hbook1(510, '# Recoils to end detector', 1, 0.0, 1.0, 0.0)
C.
C.    Number of BGOs trigger in coincedence with recoil detection
      CALL hbook1(511, '# BGOs triggered, in coin.', 10, 0.0, 10.0, 0.0)
C.
C.    Sum of all energy deposited in BGOs in coincedence with recoil detection
      CALL hbook1(512, 'Energy in BGOs, in coin.', 250, 0.0, 25.0, 0.0)
C.
C.    Sum of all energy deposited in BGOs vrs. Number of BGOs triggered,
C.    in coincedence with recoil detection
      CALL hbook2(513, 'Energy in BGOs vs # BGOs, in coin.', 10, 0.0,
     +            10.0, 250, 0.0, 25.0, 0.0)
C.
C.    Energy deposited in BGO with second greatest energy vrs energy deposited
C.    in BGO with greatest energy, in coincedence with recoil detection
      CALL hbook2(514, 'Energy of 2nd Energetic BGO vs 1st Energergetic
     +            BGO, in coin.', 250, 0.0, 25.0, 250, 0.0, 25.0, 0.0)
C.
C.    # of BGO with second greatest energy vrs # of BGO with greatest
C.    energy, in coincedence with recoil detection
      CALL hbook2(515, '# 2nd Energetic BGO vrs # 1st Energetic BGO, in
     +            coin.', 30, 1.0, 30.0, 30, 1.0, 30.0, 0.0)
C    
      CALL hbook1(516, 'Detected Energy (with P.H.D included)',4000,
     +     0.0,20.0,0.0)
C
      CALL hbook1(517,'Local t.o.f (MCP-DSSSD)',1000,0.,200.,0.)
      CALL hbook1(518,'Local t.o.f (MCP1-MCP0)',1000,0.,200.,0.)
C
      CALL hbook1(520,'Reaction x-position',1000,-5.,5.,0.)
      CALL hbook1(521,'reaction y-position',1000,-5.,5.,0.)

      CALL hbook2(522,'Q1 start x-y hit pattern',16,-2.45,2.45,16,
     +     -2.45,2.45,0.0)
C
C     Alpha source histogram(s)
      CALL hbook1(523,'Alpha source energy',512,0.,4096.,0)
C      CALL hbfun2(524,'Alpha position',100,-0.6,0.6,100,-0.6,0.6,
C     +     alphapos)
c      CALL hcopy (524,525,' ')

      ! uses a 0.6% detector efficiency convolution with 2 exponentials and a fermi function to fit the observed spectrum
      open(unit=10,file='hist6p2e.txt',status='OLD')


      nlines = 0
      do i = 1, 4096
       read(10,*,end=99) ener, weight
       nlines = nlines + 1
c     print*, ener, weight
       CALL hfill(523,ener,0,weight)
      enddo
 99   continue
      close(10)
      

C.
C.--> New ntuples 07.07.03
C.
C  'History' ntuple
      CALL HBNT(1000,'HISTORY',' ')
      CALL HBNAME(1000,'HISTORY',E_int,'E_int:R,E_rec:R,E_g(15):R,' //
     +                   'E_gp(15):R,cost_g(15):R,phi_g(15):R,' //
     +                   'cost_gp(15):R,cost_r:R,cosp_r:R,' //
     +                   'Nodec:I,' //
     +                   'react:I,recdet:I,x_r:R,y_r:R,'//
     +                   'z_r:R,thet_r:R,xstop:R,ystop:R,zstop:R,'//
     +                   'xint:R,yint:R,zint:R,'//
     +    'x:R,y:R,xp:R,yp:R,z:R,xtest(10):R,ytest(10):R,etest(10):R,'//
     +     'dsssdpos:I,binitenerg:R,half_thickness:R,'//
     +     'stop_power:R,beamtof:R' )
C.
C.--> New nutple 20040407
C.
C 'Gammahit' ntuple
      CALL HBNT(1001,'GAMMAHITS',' ')
      CALL HBNAME(1001,'GAMMAHIT',recoil_hit_ENDV,'recoil_hit_ENDV:I,'//
     + 'num_bgos_hit:I,e_bgos_total:R,num_bgo_first:I,'//
     + 'e_bgo_first:R,num_bgo_second:I,e_bgo_second:R,'//
     + 'pair_productions:I,num_bgos_hit_ab:I,'//
     + 'num_bgo_first_ab:I,e_bgo_first_ab:R,num_bgo_second_ab:I,'//
     + 'e_bgo_second_ab:R,gammatof:R,gammatofp:R,e0_conv:R,'//
     + 'isBGO:I,isLaBr3:I')


C -->   Define other ntuples
C.
      If(iswit(9).eq.2)then
C Gamma -HI coincidence ntuple
        CALL HBOOKN (100,'Gamma-HI',nvar10,'//HBOOK',1024,CHTAGS)
C.
C.-->   Setup parameters - filled once at the end of UGINIT
C.
         CALL HBNT(998,'GBOX Geant Setup Ntuple',' ')
         CALL HBNAME(998,'U_GEOM',x1_fngr,CH_U_GEOM)
C.
C.-->   Event variables - filled every event at the end of GUDIGI
C.
         CALL HBNT(999,'GBOX Geant Event Ntuple',' ')
         CALL HBNAME(999,'GBOX_INF',melem_gbox,CH_GBOX_INFO)
C.
      Endif
C.
      RETURN
      END
C.
      SUBROUTINE namfil(wstart,inum,wend,wfile)
C.
C-----------------------------------------------------------------------
C       Subroutine to create a file name containing a number inbedded
C
C       Input arg :     wstart  - Character string to be placed at the
C                                 start of the file name.
C
C                       inum    - I*4 to be converted to ASCII character
C                                 and appended to 'wstart' in file name.
C                                 No blanks or zeroes will be placed
C                                 before the number.
C
C                       wend    - Character string to terminate the file
C                                 name.
C
C       Output arg:     wfile   - Character string containing the file
C                                 name. The calling program must have
C                                 defined 'wfile' big enough to contain
C                                 all characters.
C
C-----------------------------------------------------------------------
C.
      IMPLICIT none
C.
      INTEGER inum, ifin, indexn
C.
      CHARACTER *(*) wstart, wend, wfile
      CHARACTER*10   wnum
C.
      Write(wnum,10)inum
   10 Format(I10)
C.
      ifin = indexn(wnum)
C.
      wfile = wstart//wnum(ifin:10)//wend
C.
      RETURN
      END
C.
