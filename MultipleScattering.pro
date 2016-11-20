PRO MultipleScattering
;--------------------------------------------------------------------------------------;
;                                   Parameters Initialize
;--------------------------------------------------------------------------------------;
    ; Lidar
    nPhotons = 1E4   ; the number of generated photons
    lambda = 532.0   ; the wavelength of the photons. Unit: nm
    rBeam = 5.0    ; the radius of the laser beam at z=0km. Unit: mm
    divBeam = 0.05   ; the divergence of the laser beam. Unit: mrad
    FOV = 1.0   ; the field of view(Full angle). Unit: mrad
    rTel = 150   ; the radius of the telescope. Unit: mm
    BinWidth = 20.0   ; the width of each bin. Unit: ns
    nBins = 200   ; the number of the recording bins
    distBeamTel = 300   ; the distance between laser beam and the centre of the 
                        ; telescope. Unit: mm
    SBeam = [1, 1, 0, 0]   ; the stokes vector of the incident laser beam.
    
    ; Medium
    clBase = 1.0   ; the cloud base. Unit: km
    nClLayers = 30.0   ; the number of the layers of the simulated cloud.
    clDh = 0.03   ; the delta h of each layer. Unit: km
    gamma = Fltarr(nClLayers)+6.0   ; gamma
    Reff = Fltarr(nClLayers)+0.677   ; effective radius. Unit: micros
    N0 = Fltarr(nClLayers)+1E10   ; droplets numbers. Unit: m^{-3}
    nAngs = 1800L   ; the number of scattering angles for Mie scattering
    relM = DComplex(1.33, 0)   ; the relative refractive index of the medium.     
    fileMie = 'MieScattering_x8.h5'   ; the h5 file containing the information about Mie scattering
    
    ; Monte-Carlo parameters
    thresholdForAlive = 0.01   ; the threshold for determining the photon condition
    chanceForAlive = 0.1   ; the chance in the roulette
    
    ALIVE = 1 & DEAD = 0
    photonStatus = DEAD   ; the status of each photon
;--------------------------------------------------------------------------------------;

;--------------------------------------------------------------------------------------;
;                                   Variables Initialize
;--------------------------------------------------------------------------------------;
    clTop = clBase + nClLayers*clDh   ; the cloud top. Unit: km
    ; Receiving Stokes vector. (200 Bins since penetrating the cloud)
    SVReturn = Fltarr(nBins, 4)
    len = Fltarr(nPhotons)   ; the length in each move. Unit: m
    nScatterings = Lonarr(nPhotons)   ; the number of the scatterings
;--------------------------------------------------------------------------------------;

;--------------------------------------------------------------------------------------;
;                                   Clouds Paramters
;--------------------------------------------------------------------------------------;
;    WaterMieScattering, nClLayers, gamma, Reff, N0, nAngs, relM, lambda, $
;                        FILE = fileMie
;--------------------------------------------------------------------------------------;

;--------------------------------------------------------------------------------------;
;                                   Read data
;--------------------------------------------------------------------------------------;
    S1Rel = ReadH5(fileMie, '/S1Rel')
    S1Img = ReadH5(fileMie, '/S1Img')
    S2Rel = ReadH5(fileMie, '/S2Rel')
    S2Img = ReadH5(fileMie, '/S2Img')
    scaAngs = ReadH5(fileMie, '/scaAngs')
    muExt = ReadH5(fileMie, '/Extinction')   ; extinction. Unit: m^{-1}
    muSca = ReadH5(fileMie, '/Scattering')   ; scattering. Unit: m^{-1}
;--------------------------------------------------------------------------------------;
    ; Usually bound corresponding to an optical depth of 2
    ; DM Winker and LR Poole, "Monte-Carlo calculations of 
    ; cloud returns for ground-based and space-based lidars," 
    ; Applied Physics B 60 (4), 341-344 (1995).   
    ; Note: we set 4. 
    clBoundX = 2.0/Mean(muExt)/1000.0   ; the boundary of the cloud 
                                        ; in X direction. Unit: km
    clBoundY = 2.0/Mean(muExt)/1000.0   ; the boundary of the cloud 
                                        ; in Y direction. Unit: km
                                        
;--------------------------------------------------------------------------------------;
;                                   Simulation
;--------------------------------------------------------------------------------------;
    dAng = scaAngs[1] - scaAngs[0]   ; the delta angle. Unit: rad
    
    ; Scattering parameters at each angle
    S1 = DComplex(S1Rel, S1Img) & S2 = DComplex(S2Rel, S2Img)
    s11 = Real_Part(0.5*(Abs(S1)^2 + Abs(S2)^2))
    s12 = Real_Part(0.5*(Abs(S2)^2 - Abs(S1)^2))
    s33 = Real_Part(0.5*(Conj(S2)*S1 + S2*Conj(S1)))
    s34 = Real_Part(DComplex(0, -0.5)*(S1*Conj(S2) - S2*Conj(S1)))
    albedo = muSca/muExt
                 
    seed = Ptr_New(100L)   ; Initial the seed pointer for random number generator
    
    FOR iPhoton = 0, nPhotons-1 DO BEGIN
    
        ; Launch
        Launch, 'Planar', (rBeam+divBeam/2.0*clBase*1000.0)/1000.0, seed, X = x, Y = y
        z = DOUBLE(clBase*1000.0) & x = DOUBLE(x) & y = DOUBLE(y)  ; Unit: m
        len[iPhoton] = Distance_Measure([[x,y,z], [0,0,0]], /DOUBLE)
        phDirCos = [x, y, z]/len[iPhoton]
        len[iPhoton] = 0.0
        SVIn = LaunchStokes(SBeam, [x, y, z])
        photonStatus = ALIVE
        iLayer = 0   ; the number of the layer that the photon is.
        
        WHILE (photonStatus) DO BEGIN
            
            ; Move to the next point
            rnd = RandomU(*seed)
            rnd = (rnd EQ 0.0)? 1-rnd : rnd
            ds = -Alog(rnd)/muExt[iLayer]
            x = x + phDirCos[0]*ds
            y = y + phDirCos[1]*ds
            z = z + phDirCos[2]*ds
            len[iPhoton] = len[iPhoton] + ds
            iLayer = FIX((z-clBase*1000.0) / (clDh*1000.0))

            ; Probability scattering into the FOV
            ; still in the Medium?
            IF (IsInMedium([x,y,z], [clBase, clTop, clBoundX, clBoundY]*1000.0)) $
                THEN BEGIN
                ; absorb?
                IF IsAbsorb(muSca[iLayer], muExt[iLayer], seed) THEN BEGIN
                    photonStatus = DEAD
                    BREAK
                ENDIF
                ; still in the FOV?
                IF (IsInFOV([x,y,z], FOV, rTel, distBeamTel)) THEN BEGIN   
                    ; Calculate the probability scattered into the Lidar

                    ; Incident and scattered direction Cosine
                    phIncDir = phDirCos
                    phScaDir = [distBeamTel/1000.0, 0, 0]-[x, y, z]
                    phScaDir = phScaDir/Sqrt(Total(phScaDir^2))
                    
                    ; scattering angle and rotation angles. Unit: rad
                    ScaRotAng, phIncDir, phScaDir, $
                               SCAANG = scaAng, ROTANG1 = rotAng1, ROTANG2 = rotAng2
                    
                    ; Solid angle of the telescope. Unit:Sr
                    solAng = 4.0*!PI*(rTel/1000.0)^2/ $ 
                             Distance_Measure([[x,y,z],[distBeamTel/1000.0,0,0]], $
                                              /DOUBLE)

                    ; probability without scattering again until entering the Lidar
                    iAng = Round(scaAng/dAng)   ; the index of the scattering angle
                    probEnter = solAng*Exp(-((Total(clDh*muExt[0:iLayer])- $ 
                               (-z/1000.0+clBase+clDh*(iLayer+1))*muExt[iLayer])/ $ 
                               Abs(phScaDir[2]))*1000.0)

                    ; Stokes vector after scattering
                    SVTemp1 = RotSphi(SVIn, rotAng1)
                    SVTemp2 = [[s11[iLayer, iAng], s12[iLayer, iAng], 0, 0], $
                               [s12[iLayer, iAng], s11[iLayer, iAng], 0, 0], $
                               [0, 0, s33[iLayer, iAng], -s34[iLayer, iAng]], $
                               [0, 0, -s34[iLayer, iAng], s33[iLayer, iAng]]] $
                               ## SVTemp1
                    SVSca = RotSphi(SVTemp2, -rotAng2)
                    SVSca = SVSca/s11[0, 0]

                    tReturn = (len[iPhoton] + (z - clBase*1000.0)/Abs(phScaDir[2]))/ $
                              !constant.c0 * 1E9   ; the time when Receiving. Unit: ns
                    iTReturn = FIX(tReturn / BinWidth)
                    SVReturn[iTReturn, *] = SVReturn[iTReturn, *] + $
                                            RotSphi(SVSca, Atanxoy(phScaDir)) * $
                                            probEnter
                    IF ~FINITE(TOTAL(SVReturn[iTReturn, *])) THEN STOP
                    ;IF scaAng LT 30.0/180.0*!PI THEN STOP                   
                ENDIF

                ; Randomly scattering

                ; Rejection method(Jessica D, Optics Express, 2005)
                REPEAT BEGIN 
                    theta = ACOS(2*RandomU(*seed) - 1.0)
                    phi = RandomU(*seed)*2.0*!PI
                    PhaseFunc0 = s11[iLayer, 0]*SVIn[0] + $
                                 s12[iLayer, 0]* $
                                 (SVIn[1]*COS(2.0*phi)+SVIn[2]*SIN(2.0*phi))
                    iAng = Round(theta / dAng)
                    PhaseFuncSca = s11[iLayer, iAng]*SVIn[0] + $
                                   s12[iLayer, iAng]* $
                                   (SVIn[1]*COS(2.0*phi)+SVIn[2]*SIN(2.0*phi))
                ENDREP UNTIL (RandomU(*seed)*PhaseFunc0 GE PhaseFuncSca)
                
                ; the Stokes Vector of the scattered photon
                phIncDirCos = phDirCos   ; saving the direction cosine of the incident
                                         ; photon
                phDirCos = UpdateDir(phDirCos, phi, theta)
                SVTemp1 = RotSphi(SVIn, phi)
                SVTemp2 = [[s11[iLayer, iAng], s12[iLayer, iAng], 0, 0], $
                           [s12[iLayer, iAng], s11[iLayer, iAng], 0, 0], $
                           [0, 0, s33[iLayer, iAng], -s34[iLayer, iAng]], $
                           [0, 0, -s34[iLayer, iAng], s33[iLayer, iAng]]] ## SVTemp1
                temp = Sqrt(((1.0-COS(theta)^2) * (1.0-phDirCos[2]^2)))
                IF (temp EQ 0.0) THEN BEGIN
                    gammaAng = !PI/2.0
                ENDIF ELSE BEGIN
                    gammaAngCos = ((phi GT !PI) AND (phi LT 2.0*!PI) ? 1.0: -1.0)* $
                                  (phDirCos[2]*COS(theta) - phIncDirCos[2]) / temp 
                    IF gammaAngCos LE -1.0 THEN gammaAngCos = -1.0 ELSE $
                    IF gammaAngCos GE 1.0 THEN gammaAngCos = 1.0
                    gammaAng = ACOS(gammaAngCos)
                ENDELSE
                SVIn = RotSphi(SVTemp2, -gammaAng)
                SVIn = SVIn/s11[0, 0]    ; normalized
                nScatterings[iPhoton] = nScatterings[iPhoton] + 1

            ENDIF ELSE BEGIN
                ; not in the medium. DEAD
                photonStatus = DEAD
                        Print, 'Photon '+ String(iPhoton, FORMAT='(I5)')+' is Dead!'
                        Print, String(len[iPhoton], FORMAT='(F6.2)')
                        Print, String(nScatterings[iPhoton], FORMAT='(I5)')
                BREAK
            ENDELSE

        ENDWHILE

    ENDFOR
;--------------------------------------------------------------------------------------;

;--------------------------------------------------------------------------------------;
;                                    Data visualization
;--------------------------------------------------------------------------------------;

    ; Cloud droplet scattering phase function
    W1 = Window(DIMENSION=[400,600])
    p1 = Plot(scaAngs*180.0/!PI, s11[0, *]/s11[0, 0], /CURRENT, $
              XTITLE='Scattering Angles($\deg$)', YTITLE='Phase Function(Normalized)', $
              XRANGE=[0,180], YRANGE=[1E-4, 1], /YLOG)
    t = Text(140, 0.9, $
             '$\sigma$='+String(muExt[0], FORMAT='(F4.1)')+ '$m^{-1}$'+ $
             '\n$R_eff$='+String(Reff[0], FORMAT='(F3.1)')+'$\mum$'+ $
             '\n$\gamma$='+String(gamma[0], FORMAT='(I02)'), /DATA, FONT_SIZE=10) 
                
    ; Stokes vector
    W2 = Window(DIMENSION=[1256,907])
    p1 = Plot(SVReturn[*, 0], Findgen(nBins)*BinWidth*!Constant.C0/1E9/2.0/1000.0, $
              /CURRENT, TITLE='Stokes Vector', $
              COLOR='r', $
              NAME='I', LAYOUT=[3, 1, 1])
    p2 = Plot(SVReturn[*, 1], Findgen(nBins)*BinWidth*!Constant.C0/1E9/2.0/1000.0, $
              /OVERPLOT, $
              COLOR='g', NAME='Q')
    p3 = Plot(SVReturn[*, 2], Findgen(nBins)*BinWidth*!Constant.C0/1E9/2.0/1000.0, $
              /OVERPLOT, $
              COLOR='b', NAME='U')
    p4 = Plot(SVReturn[*, 3], Findgen(nBins)*BinWidth*!Constant.C0/1E9/2.0/1000.0, $
              /OVERPLOT, $
              COLOR='cyan', NAME='V')
    l1 = Legend(TARGET=[p1, p2, p3, p4], POSITION=[0.2, 0.8])
    
    ; Polarization
    deRatio = (SVReturn[*, 0] - SVReturn[*,1])/(SVReturn[*, 0] + SVReturn[*,1])
    p5 = Plot(deRatio, Findgen(nBins)*BinWidth*!Constant.C0/1E9/2.0/1000.0, $
              /CURRENT, LAYOUT=[3,1,2], $
              XRANGE=[0, 1.0], YRANGE=[0, 0.8], $
              XTITLE='Depolarization Ratio', YTITLE='Height(km)')
    
    ; Number of scatterings
    nNScatterings = HistoGram(nScatterings, BINSIZE=1, MAX=30, MIN=1, $
                              REVERSE_INDICES=rIndex, LOCATIONS=temp)
    p6 = BarPlot(Findgen(30)+1, nNScatterings, /CURRENT, LAYOUT=[3,1,3], $
                 XRANGE=[-1,32], $
                 XTITLE='N Scatterings', YTITLE='Frequency', $
                 FILL_COLOR='r')                          
;--------------------------------------------------------------------------------------;
    
END
