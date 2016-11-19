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
    BinWidth = 200.0   ; the width of each bin. Unit: ns
    distBeamTel = 300   ; the distance between laser beam and the centre of the 
                        ; telescope. Unit: mm
    SBeam = [1, 1, 0, 0]   ; the stokes vector of the incident laser beam.
    
    ; Medium
    clBase = 1.0   ; the cloud base. Unit: km
    nClLayers = 30.0   ; the number of the layers of the simulated cloud.
    clDh = 0.03   ; the delta h of each layer. Unit: km
    gamma = Fltarr(nClLayers)+6.0   ; gamma
    Reff = Fltarr(nClLayers)+8.0   ; effective radius. Unit: micros
    N0 = Fltarr(nClLayers)+1E8   ; droplets numbers. Unit: m^{-3}
    nAngs = 1800L   ; the number of scattering angles for Mie scattering
    relM = DComplex(1.33, 0)   ; the relative refractive index of the medium.     
    fileMie = 'MieScattering.h5'   ; the h5 file containing the information about Mie scattering
    
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
    SVReturn = Fltarr(200, 4)
    len = 0.0   ; the length in each move. Unit: m
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
        len = Distance_Measure([[x,y,z], [0,0,0]], /DOUBLE)
        phDirCos = [x, y, z]/len
        len = 0.0
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
            len = len + ds
            iLayer = FIX((z-clBase*1000.0) / (clDh*1000.0))

            ; absorb?
            IF IsAbsorb(muSca[iLayer], muExt[iLayer], seed) THEN BEGIN
                photonStatus = DEAD
                BREAK
            ENDIF

            ; Probability scattering into the FOV
            ; still in the Medium?
            IF (IsInMedium([x,y,z], [clBase, clTop, clBoundX, clBoundY]*1000.0)) $
                THEN BEGIN
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
                               COS(phScaDir[2]))*1000.0)

                    ; Stokes vector after scattering
                    SVTemp1 = RotSphi(SVIn, rotAng1)
                    SVTemp2 = [[s11[iLayer, iAng], s12[iLayer, iAng], 0, 0], $
                               [s12[iLayer, iAng], s11[iLayer, iAng], 0, 0], $
                               [0, 0, s33[iLayer, iAng], -s34[iLayer, iAng]], $
                               [0, 0, -s34[iLayer, iAng], s33[iLayer, iAng]]] $
                               ## SVTemp1
                    SVSca = RotSphi(SVTemp2, -rotAng2)

                    tReturn = (len + (z - clBase/1000.0)/COS(phScaDir[2]))/ $
                              !constant.c0 * 1E9   ; the time when Receiving. Unit: ns
                    iTReturn = FIX(tReturn / BinWidth)
                    SVReturn[iTReturn, *] = SVReturn[iTReturn, *] + $
                                            RotSphi(SVSca, -Atanxoy(phScaDir)) * $
                                            probEnter                    
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
                    gammaAng = ACOS(((phi GT !PI) AND (phi LT 2.0*!PI) ? 1.0: -1.0)* $
                                    (phDirCos[2]*COS(theta) - phDirCos[2]) / temp )
                ENDELSE
                SVIn = RotSphi(SVTemp2, -gammaAng)

            ENDIF ELSE BEGIN
                ; not in the medium. DEAD
                photonStatus = DEAD
                BREAK
            ENDELSE

        ENDWHILE

    ENDFOR
;--------------------------------------------------------------------------------------;

    ; Data visualization
    
END
