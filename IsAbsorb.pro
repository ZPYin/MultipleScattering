FUNCTION IsAbsorb, muSca, muExt, seed
    
    IF RandomU(*seed) LT muSca/muExt THEN Return, 0 ELSE Return, 1
    
END