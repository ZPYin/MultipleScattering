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
    nClLayers = 100.0   ; the number of the layers of the simulated cloud.
    clDh = 0.003   ; the delta h of each layer. Unit: km
    thresholdForAlive = 0.01   ; the threshold for determining the photon condition
    chanceForAlive = 0.1   ; the chance in the roulette
    fileMie = ''   ; the h5 file containing the information about Mie scattering
    
    ALIVE = 1 & DEAD = 0
    photonStatus = DEAD   ; the status of each photon
;--------------------------------------------------------------------------------------;

;--------------------------------------------------------------------------------------;
;                                   Variables Initialize
;--------------------------------------------------------------------------------------;
    clTop = clBase + nClLayers*clDh   ; the cloud top. Unit: km
    ; Receiving Stokes vector. (200 Bins since penetrating the cloud)
    IR = Fltarr(200)
    QR = Fltarr(200)
    UR = Fltarr(200)
    VR = Fltarr(200)
    len = 0.0   ; the length in each move. Unit: m
;--------------------------------------------------------------------------------------;

;--------------------------------------------------------------------------------------;
;                                   Read data
;--------------------------------------------------------------------------------------;
    S1Rel = ReadH5(fileMie, '/S1Rel')
    S1Img = ReadH5(fileMie, '/S1Img')
    S2Rel = ReadH5(fileMie, '/S2Rel')
    S2Img = ReadH5(fileMie, '/S2Img')
    scaAngs = ReadH5(fileMie, '/scaAngs')
    muExt = ReadH5(fileMie, '/Extinction')   ; extinction. Unit: km^{-1}
    muSca = ReadH5(fileMie, '/Scattering')   ; scattering. Unit: km^{-1}
;--------------------------------------------------------------------------------------;
    ; Usually bound corresponding to an optical depth of 2
    ; DM Winker and LR Poole, "Monte-Carlo calculations of 
    ; cloud returns for ground-based and space-based lidars," 
    ; Applied Physics B 60 (4), 341-344 (1995).    
    clBoundX = 1.0/Mean(muExt)   ; the boundary of the cloud in X direction. Unit: km
    clBoundY = 1.0/Mean(muSca)   ; the boundary of the cloud in Y direction. Unit: km
;--------------------------------------------------------------------------------------;
;                                   Simulation
;--------------------------------------------------------------------------------------;
    nAngs = Size(scaAngs, /DIM)[1]   ; the number of the angles of the phase function
    
    ; Scattering parameters at each angle
    S1 = DComplex(S1Rel, S1Img) & S2 = DComplex(S2Rel, S2Img)
    s11 = 0.5*(Abs(S1)^2 + Abs(S2)^2)
    s12 = 0.5*(Abs(S2)^2 - Abs(S1)^2)
    s33 = 0.5*(Conj(S2)*S1 + S2*Conj(S1))
    s34 = DComplex(0, -0.5)*(S1*Conj(S2) - S2*Conj(S1))
    albedo = muSca/(muSca+muExt)
    
    seed = Ptr_New(100L)   ; Initial the seed pointer for random number generator
    
    FOR iPhoton = 0, nPhotons-1 DO BEGIN
    
        ; Launch
        Launch, 'Planar', (rBeam+divBeam/2.0*clBase*1000.0)/1000.0, seed, X = x, Y = y
        z = clBase*1000.0   ; Unit: m
        len = Distance_Measure([[x,y,z], [0,0,0]], /DOUBLE)
        phDirCos = [x, y, z]/len
        SVIn = LaunchStokes(SBeam, [x, y, z])
        photonStatus = ALIVE
        iLayer = 0   ; the number of the layer that the photon is.
        
        WHILE (photonStatus) DO BEGIN
            
            ; Move to the next point
            rnd = RandomU(*seed)
            rnd = (rnd EQ 0.0)? 1-rnd : rnd
            len = len-Alog(rnd)/(muSca[iLayer]+muExt[iLayer])
            x = x + phDirCos[0]*len
            y = y + phDirCos[1]*len
            z = z + phDirCos[2]*len
            iLayer = FIX((z-clBase*1000.0) / (clDh*1000.0))

            ; absorb?
            IF IsAbsorb(muSca[iLayer], muExt[iLayer], seed) THEN BEGIN
                photonStatus = DEAD
                BREAK
            ENDIF

            ; Probability scattering into the FOV
            ; still in the FOV?
            IF (IsInFOV([x,y,z], FOV, rTel, distBeamTel)) THEN BEGIN
                ; still in the Medium?
                IF (IsInMedium([x,y,z], [clBase, clTop, clBoundX, clBoundY]*1000.0)) $
                    THEN BEGIN
                    it test
                ENDIF
            ENDIF


        ENDWHILE

    ENDFOR
;--------------------------------------------------------------------------------------;

END
