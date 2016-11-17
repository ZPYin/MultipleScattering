FUNCTION IsAbsorb, muSca, muExt, seed
    
    IF muExt/(muSca+muExt) LT RandomU(*seed) THEN Return, 0 ELSE Return, 1
    
END