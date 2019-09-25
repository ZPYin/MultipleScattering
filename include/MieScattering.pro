;+
;
; :Author:
;  Wang Wei
;
; :Description:
;    mie散射程序，输出s1，s2，Qext ，Qsca. 
;    还包括不同散射角退偏比计算，但没有输出
;
; :Params:
;    m: DCOMPLEX(1.33,1)  ;相对折射率
;    x: 相对尺度，相对波长。 x = 2*!PI*r/lambda
;    theta: scattering angle. Unit: deg
; :Output:
;   [s1,s2,Qext,Qsca]
;   
; :Examples:
;   data=MieScattering(DCOMPLEX(1.33,1),8, 60)  
;   
; :History:
;  2016-11-16
;  
;-
;
;
;
FUNCTION MieScattering, m, x, theta 
      mx=m*x
      nmax = ROUND(x+4*X^1/3+2) ;最大n值
      
      ;递推求phi，xi ，xi=phi+i*chi
      phiX = DCOMPLEX(DBLARR(nmax+1))
      phimx = DCOMPLEX(DBLARR(nmax+1))
      phidx = DCOMPLEX(DBLARR(nmax+1))
      phidmx = DCOMPLEX(DBLARR(nmax+1)) 
      xix = DCOMPLEX(DBLARR(nmax+1))
      xidx = DCOMPLEX(DBLARR(nmax+1)) 
      chix = DCOMPLEX(DBLARR(nmax+1))
      chidx = DCOMPLEX(DBLARR(nmax+1)) 
      
      phix[0] = SIN(X) 
      phix[1] = 1./x*phix[0]-COS(x)
      phimx[0] = SIN(mX) 
      phimx[1] = 1./mx*phimx[0]-COS(mx)
      chix[0] = cos(X) 
      chix[1] = 1./x*chix[0]+sin(x)
            
      FOR i=2, nmax DO BEGIN
        phix[i]=(2.*i-1)/x*phix[i-1]-phix[i-2]
        phimx[i]=(2.*i-1)/mx*phimx[i-1]-phimx[i-2]
        chix[i]=(2.*i-1)/x*chix[i-1]-chix[i-2]
      ENDFOR
      
      FOR i=1, nmax DO BEGIN
        phidx[i]=-1.*i/x*phix[i]+phix[i-1]
        phidmx[i]=-1.*i/mx*phimx[i]+phimx[i-1]
        chidx[i]=-1.*i/x*chix[i]+chix[i-1]
      ENDFOR
      
      xix=DCOMPLEX(phix,-1.*chix)  ;符号 是加是减？
      xidx=DCOMPLEX(phidx,-1.*chidx)
      
      anup=phix*phidmx-m*phidx*phimx
      andown = xix*phidmx-m*xidx*phimx
      an=anup/andown     
      bnup=m*phix*phidmx-phidx*phimx
      bndown = m*xix*phidmx-xidx*phimx
      bn=bnup/bndown 
      

;      输出an，bn以测试
;      print,'an=',an  ;10.26 和matlab算出来的虚数部分符号相反，an[1]数值不对
;                10.27 更改后正常，xix符号变了 
;      print,'bn=',bn
      
      ;--------------计算pi，tao  算出来结果可靠----------------------
      tao = DBLARR(nmax+1) ;tao[0]并没有作用
      pi = DBLARR(nmax+1)  ;pi[0]并没有作用
      pi[0]=0.
      pi[1]=1.
      tao[0]=0.
      tao[1]=cos(theta*!DTOR)
      FOR i=2 ,nmax DO BEGIN
        pi[i] = (2.*i-1.)/(i-1.)*cos(theta*!DTOR)*pi[i-1]  $
                 -i/(i-1.)*pi[i-2]
        tao[i] = i*cos(theta*!DTOR)*pi[i]-(i+1.)*pi[i-1]       
      ENDFOR
      
      
      ;-----------计算各种参数-------------------------

      ;**********计算振幅函数s1，s2 ***************
      S1n = DCOMPLEX(DBLARR(nmax+1))
      S2n = DCOMPLEX(DBLARR(nmax+1))
      FOR N=1, nmax DO BEGIN
        s1n[n] = (2*n+1.)/(n*(n+1.))*(an[n]*pi[n]+bn[n]*tao[n])
        s2n[n] = (2*n+1.)/(n*(n+1.))*(bn[n]*pi[n]+an[n]*tao[n])
      ENDFOR
      s1 = TOTAL(s1n[1:nmax])
      s2 = TOTAL(s2n[1:nmax])
      
      
      ;**********消光系数Qext ,散射系数Qsca*********
      Qext_N = DBLARR(nmax+1)
      Qsca_N = DBLARR(nmax+1)
      FOR i=1,nmax DO BEGIN
        Qext_N(i) = (2.*i+1.)*REAL_PART(an(i)+bn(i))
        Qsca_N(i) = (2.*i+1.)*(ABS(an(i))^2+ABS(bn(i))^2)
      ENDFOR
      Qext=2./x^2*TOTAL(Qext_N(1:nmax))
      Qsca=2./x^2*TOTAL(Qsca_N(1:nmax))
      
      
;      ;********去极化率delta****************
;      F1=S2*CONJ(S2)
;      F2=S1*CONJ(S1)
;      F3=(S1*CONJ(S2)+S2*CONJ(S1))/2.
;      delta = (F2*COS(theta*!DTOR)^2+F1-2*F3*COS(theta*!DTOR))/ $
;              (3*F2*COS(theta*!DTOR)^2+3*F1+2*F3*COS(theta*!DTOR))
;      
;      
;      ;********i1,i2:散射强度函数*******
;      i1=ABS(S1)^2
;      i2=ABS(S2)^2
      
;      data = COMPLEX(DBLARR(4))
      data = [s1,s2,Qext,Qsca]
      RETURN , data
       
 END    
