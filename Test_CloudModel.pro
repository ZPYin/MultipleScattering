PRO Test_CloudModel
    Rarr = Sequence(0.1, 8.0, 0.01)*1E-6    
    P1 = Plot(Rarr, CloudModel_LGD('lognm', $
        BETASCA=1e-6, GAMMA=1, REFF=2e-6, Radius=rarr), XRANGE=[Rarr[0], Rarr[-1]])
END