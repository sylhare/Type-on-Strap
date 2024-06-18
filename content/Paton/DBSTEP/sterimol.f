C ----
C ---- Minute changes to sterimol.for to make program compiled under GNU g77
C ---- by Jan Labanowski, Feb. 28, 2009
C ----
C  *****   STERIMOL  PROGRAM  (DATE: SEPT. 1976)      *****
C  *****   INPUT  CARDS   *****                                     
C  THE FIRST CARD CONTAINS THE  33  USED SYMBOLS; STARTING WITH 1       
C  IN COLUMN  1.  SYMBOL NR  27  IS A SPACE.                        
C      1234567890ABCDEFHINOPRSTXZ (&*),=                      
C  NEXT FOLLOW (IF WANTED) THE FORMULAE   FOR THE 'RADICALS',             
C  OF THE FOLLOWING GENERAL TYPE:                             
C       ZKK=FORMULA*                                                
C  WHERE  KK  IS A TWO-DIGIT INTEGER NUMBER  (01 - 99).                   
C  A RADICAL MAY BE ANY SYMBOL STRING WHICH IS PART OF                    
C  A VALID FORMULA.  IF NECESSARY, IT MAY BE CONTINUED                    
C  LIKE A FORMULA CARD.  CLOSE OR REPLACE THE                       
C  RADICAL-CARD SET BY A BLANK CARD.                                
C    FOR EACH MOLECULE THE FOLLOWING CARDS ARE GIVEN:               
C  (1)  THE INTEGER  IPR  (I4): THE OUTPUT INDEX (1 OR 2);                
C  (2)  THE FORMULA.  IF ONE CARD IS NOT SUFFICIENT, END                  
C  THIS CARD WITH A  &  SIGN AND CONTINUE ON THE NEXT CARD.               
C  (3)  NUMBER OF THE TORSION ANGLES (I4).  IF NONE ARE TO BE             
C  GIVEN, REPLACE THIS CARD BY A BLANK ONE.                         
C  (4)  THE TORSION ANGLES (1X,F7.2; DECIMAL POINTS AT                    
C  COLUMNS  6, 14, 22, ....) .                                      
C  AFTER THE SETS OF MOLECULE-INPUT CARDS:  ONE BLANK CARD.               
C  *****   OUTPUT   *****                                           
C       FOR  IPR = 1 :                                        
C  THE RADICAL FORMULAE (IF PRESENT);                               
C  THE MOLECULE FORMULA AS IT WAS INPUT;                      
C  SERIAL NUMBER, NUMBER OF ATOMS, OUTPUT INDEX;                    
C  FULLY WRITTEN FORMULA, WITH 'TRANSLATION' OF THE RADICALS;             
C  LIST OF THE TORSION ANGLES;                                      
C  TABLE OF THE VANDERWAALS RADII AND COORDINATES OF THE ATOMS;           
C  LIST OF THE STERIC PARAMETERS: L; B(1) - B(4)  AND               
C  B(5) (MAXIMUM WIDTH)                                             
C  A WARNING IF NEIGHBOURING ATOMS HAVE   OVERLAPPING VANDERWAALS SPHERES 
C       FOR  IPR = 2  IN ADDITION:                            
C  AN  X, Y  AND A  Y, Z  PROJECTION OF   THE MOLECULE;                   
C  A LIST OF THE ATOMS, ORDERED BY INCREASING  X  COORDINATES,            
C  SERVING TO RAPIDLY IDENTIFY THE PLOTTED ATOMS.                   
C                                                             
      INTEGER DD
      COMMON/ALG/NR,N,IPR,DD(99,100),NA(99),KH(500)                 
      COMMON/CHAR/ICH(33)
      COMMON/ZDD/NNVZ,ZDAT(30)
      NR=0
C  READ THE RADICAL FORMULAE, IF PRESENT.
    1 CALL RECHA(KH,IA)                                             
      IF (IA.EQ.0) GOTO 5                                           
      J1=MOD(KH(2),10)                                              
      J2=MOD(KH(3),10)                                              
      J=10*J1+J2                                                    
      IF(J.GT.0) GOTO 3
C  THE RADICAL Z00 REPRESENTS A SPECIAL ATOM;
C  THIS IS CALLED Z, WHICH SYMBOL IS FOLLOWED BY THE NUMBER OF BONDS.
      NNVZ=KH(6)
      K=4*NNVZ+1
C  READ THE V.D.WAALS RADIUS (*100) AND THE BOND VALUES,
C  EACH FOLLOWED BY THE CORRESPONDING BOND VECTOR.
      READ(5,*) (ZDAT(I), I=1,K)
      GOTO 1
   3  CONTINUE
      N=IA-4                                                        
      NA(J)=N                                                       
      DO 2 I=1,N                                                    
    2 DD(J,I)=KH(I+4)                                               
      GOTO 1                                                        
    5 CALL STRML                                                    
      CALL PARAM                                                    
      IF(IPR.EQ.2) CALL TRANSP                                      
      END
       BLOCK DATA
C   THE ARRAYS  IDAT  AND  VDAT CONTAIN   THE ATOMIC DATA.                
C   IDAT CONSISTS OF SEQUENCES OF THE FOLLOWING KIND:               
C        100 + K , W , N , N * (B, -J) ,  IN WHICH:                 
C   K = ATOM-TYPE NUMBER:                                           
C   W = 100 * VANDERWAALS RADIUS;                             
C   N = NUMBER OF BONDS;                                            
C   B = BOND VALUE;                                                 
C   J = FIRST SUBSCRIPT NUMBER OF A TRIPLE OF ARRAY ELEMENTS              
C     OF ARRAY  VDAT, GIVING THE COORDINATES OF THE PARTIAL BOND.       
       COMMON/CHAR/ICH(33)
       COMMON/DAT1/IDAT(184),VDAT(177)
       DATA IDAT/101,150,4,1,-1,1,-4,1,-7,1,-10,103,160,3,2,-13,1,-16,  
     1 1,-19,104,160,2,3,-22,1,-25,105,150,3,4,-28,2,-31,1,-34,106,     
     2 170,3,3,-37,3,-40,1,-43,107,170,3,2,-46,2,-49,1,-52,108,170,3,   
     3 3,-55,2,-58,1,-61,109,150,4,2,-64,2,-67,1,-70,1,-73,110,100,1,1, 
     4 -76,111,150,4,1,-79,1,-82,1,-85,1,-88,113,170,3,2,-166,2,-169,   
     5 2,-172,115,145,3,5,-91,4,-94,1,-97,120,135,2,1,-100,1,-103,122,  
     6 135,1,2,-106,123,140,4,1,-109,1,-112,1,-115,1,-118,124,170,2,1,  
     7 -121,1,-124,125,100,6,1,-127,1,-130,1,-133,1,-136,1,-139,1,-142, 
     8 126,135,1,1,-145,127,180,1,1,-148,128,140,4,1,-151,1,-154,1,     
     9 -157,1,-160,129,195,1,1,-163,130,215,1,1,-175/               
       DATA VDAT/-.77,0.,0.,.256667,.725963,0.,.256667,-.362981,        
     1 .628702,.256667,-.362981,-.628702,-.67,0.,0.,.421,.596,0.,       
     5 .421,-.596,0.,                                               
     2 -.6,0.,0.,.7,.0001,0.,-.72,0.,0.,.384,.549,0.,.317,-.68,0.,      
     3 -.69,0.,0.,.213222,.656229,0.,.429083,-.590582,0.,-.69,0.,0.,    
     4 .345,.597558,0.,.365,-.632119,0.,-.69,0.,0.,.4617,.51277,0.,     
     5 .213222,-.656229,0.,-.76,0.,0.,-.38,-.193,-.629,.348,.686,0.,    
     6 .348,-.568,.384,-.33,0.,0.,-.70,0.,0.,.23333,.659966,0.,           
     7 .233333,-.329983,.571548,.233333,-.329983,-.571548,-.7,0.,0.,    
     8 .327,.503,0.,.285,-.639,0.,-.66,0.,0.,.226,.62,0.,-.57,0.,0.,    
     9 -.96,0.,0.,.32,.905097,0.,.32,-.452548,.783837,.32,-.452548,     
     1  -.783837,-1.04,0.,0.,.252,1.009,0.,-1.0,0.,                 
     2 0.,0.,1.0,0.,0.,0.,1.0,0.,-1.0,0.,0.,0.,-1.0,1.0,0.,0.,-.57,     
     3 0.,0.,-.99,0.,0.,-.77,0.,0.,.256667,.725963,0.,.256667,            
     4 -.362981,.628702,.256667,-.362981,-.628702,-1.14,0.,0.,            
     5 -.69,0.,0.,.345,.597558,0.,.345,-.597558,0.,-1.35,0.,0./
       CHARACTER*4 ICH
       DATA ICH/'1','2','3','4','5','6','7','8','9','10',
     1 'A','B','C','D','E','F','H','I','N','O','P','R','S',
     2 'T','X','Z',' ','(','&','*',')',',','='/
       END
       SUBROUTINE DATPR(JS)                                         
C   THIS SUBROUTINE EXTRACTS, FOR A PARTICULAR ATOM TYPE  JS,             
C   THE ATOMIC DATA FROM THE BLOCK-DATA VALUES AND PUTS                 
C   THEM INTO THE VARIABLES OF THE COMMON/DAT2/.                    
       INTEGER ZZV,H
       COMMON/DAT1/IDAT(184),VDAT(177)                              
       COMMON/DAT2/NDAT,IAN,NNV,ZZV(6),VV(6,3)                      
       JSS=JS
       IS=0                                                   
       IF(JSS.EQ.14.OR.JSS.EQ.16.OR.JSS.EQ.17) IS=1                 
       IF(IS.EQ.1) JSS=JSS-10                                       
       DO 30 K=1,NDAT                                               
       IF(IDAT(K).NE.100+JSS) GOTO 30                               
       K1=K+1                                                       
       GOTO 40                                                      
   30  CONTINUE                                                     
   40  IAN=IDAT(K1)                                                 
     
       NNV=IDAT(K1+1)                                               
      K3=K1+2                                                       
      DO 60 J=1,NNV                                                 
      ZZV(J)=IDAT(K3)                                               
      H=-1-IDAT(K3+1)                                               
      DO 65 I=1,3                                                   
   65 VV(J,I)=VDAT(H+I)                                             
   60 K3=K3+2                                                       
      IF(JSS.NE.10) GOTO 70                                         
      NNV=2                                                   
      ZZV(2)=5                                                      
      VV(2,1)=-.33                                                  
      VV(2,2)=0.                                                    
      VV(2,3)=0.                                                    
   70 CONTINUE                                                      
      RETURN                                                        
      END                                                     
      SUBROUTINE ORIENT(IU,NNV,VV,BV,CV,PHI)                        
C  THIS SUBROUTINE IS CALLED BY SUBROUTINE DRAAI.                   
      INTEGER H
      DIMENSION VV(6,3),BV(3),CV(3)                           
      DIMENSION X(5),Y(5),Z(5),G(3,3)                               
      P= .017453293
      COF= COS(P*PHI)                                               
      SIF= SIN(P*PHI)                                               
      SS= 0.0                                                       
      DO 230 KM= 1,3                                                
  230 SS= SS+BV(KM)**2                                              
      SS= SQRT(SS)                                                  
      SSS= 0.0                                                      
      DO 240 KM= 1,3                                                
      X(KM)= BV(KM)/SS                                              
  240 SSS= SSS+X(KM)*CV(KM)                                         
      SS= 0.0                                                       
      DO 250 KM= 1,3                                                
      Y(KM)= CV(KM)-SSS*X(KM)                                       
  250 SS= SS+Y(KM)**2                                               
      SS= SQRT(SS)                                                  
      DO 260 KM= 1,3                                                
  260 Y(KM)= Y(KM)/SS                                               
      X(4)= X(1)                                                    
      X(5)= X(2)                                                    
      Y(4)= Y(1)                                                    
      Y(5)= Y(2)                                                    
      DO 270 KM= 1,3                                                
      Z(KM)= X(KM+1)*Y(KM+2)-Y(KM+1)*X(KM+2)                        
      IF(IU.EQ.-1) GOTO 265                                         
      G(KM,1)= X(KM)                                                
      G(KM,2)= Y(KM)                                                
      G(KM,3)= Z(KM)                                                
      GOTO 270                                                      
  265 G(1,KM)=X(KM)                                                 
      G(2,KM)=Y(KM)                                                 
      G(3,KM)=Z(KM)                                                 
  270 CONTINUE                                                      
      DO 280 H= 1,NNV                                               
      SS= COF*VV(H,2)-SIF*VV(H,3)                             
      SSS= COF*VV(H,3)+SIF*VV(H,2)                            
      VV(H,2)= SS                                                   
      VV(H,3)= SSS                                                  
      DO 290 KM= 1,3                                                
      X(KM)=0.0                                                     
       DO 300 J= 1,3                                                
  300  X(KM)=X(KM)+G(KM,J)*VV(H,J)                            
  290  CONTINUE                                                     
       DO 280 KM= 1,3                                               
  280  VV(H,KM)=X(KM)                                               
       RETURN                                                       
        END                                                   
       SUBROUTINE DRAAI(N,A,B,C,BA,CA,PHI)                          
C   THIS SUBROUTINE GIVES TO A N-VALUED   VECTOR SYSTEM   A  SUCH         
C   AN ORIENTATION, THAT                                            
C   (1)  A VECTOR BA, CONNECTED TO THIS   SYSTEM, WILL POINT IN AN        
C   OPPOSITE DIRECTION AS A GIVEN VECTOR  B  ;                      
C   (2)  ANOTHER VECTOR  CA, CONNECTED TO THE SYSTEM, WILL MAKE A       
C   DIHEDRAL ANGLE  PHI  WITH  B  AND A   GIVEN VECTOR  C ;               
C   CA  IS TURNED CLOCK-WISE, LOOKING IN THE DIRECTION OF  B .            
       DIMENSION A(6,3),B(3),C(3),BA(3),CA(3)                       
       DIMENSION BAA(3)                                             
       DO 1 K=1,3                                                   
    1  BAA(K)=-BA(K)                                                
       CALL ORIENT(-1,N,A,BAA,CA,0.0)                               
       CALL ORIENT(1,N,A,B,C,PHI)                             
       RETURN                                                       
       END                                                    
       SUBROUTINE STRML                                             
C   THIS IS THE CENTRAL SUBROUTINE, WHICH TRANSLATES THE GIVEN            
C   INPUT FORMULA INTO THE MOLECULAR STRUCTURE, GIVEN BY                  
C   THE ATOMIC COORDINATES OF ARRAY  C.                             
CCCCC
CCCCC The following notes are discoveries made while attempting to figure out
CCCCC how this mess works.  They are here so that you, poor soul that you are,
CCCCC don't have to reinvent the wheel (well, at least not the whole thing).
CCCCC
CCCCC Variable    Description
CCCCC --------    -----------
CCCCC B(x,y)            Number of bonds between atom number x and atom number
CCCCC             y.
CCCCC F()         Stack for keeping track of the first atom of the
CCCCC             chain (this does not include the first atom of the
CCCCC             main chain--only side chains).
CCCCC H           Stack pointer to F().
CCCCC IC          Atom number of last atom in this chain.
CCCCC IZ          Number of bonds from current atom to previous atom.
CCCCC K           Number of the current atom we are looking at.
CCCCC J           Position in the input buffer (SY(J) is the character
CCCCC             we are currently looking at).
CCCCC R(x,y)            Torsion angle between atom number x and atom number y.
CCCCC SY()        Input buffer.
CCCCC
CCCCC This is it.  Good luck--MLR
CCCCC
       IMPLICIT INTEGER(A,B,D,F,H,S,Z)                              
       COMMON/CHAR/ICH(33)
       COMMON/ALG/NR,N,IPR,DD(99,100),NA(99),SY(500)                
       COMMON AN(65),S(65),C(65,3),B(65,65),R(65,65)                  
       COMMON/DAT2/NDAT,IAN,NNV,ZZV(6),VV(6,3)                      
       COMMON/ZDD/NNVZ,ZDAT(30)
       DIMENSION REST(50),BBB(9),
     1 BB(9),IR(65,3),BH(65),F(20),NV(65),ZV(6,65),V(6,3,65),BK(3),     
     2 CK(3),BI(3),CI(3)                                            
       REAL BI,BK                                                   
       NDAT=184                                                     
       READ(5,*,END=9999) IPR
       IF(IPR.EQ.0) GOTO 9999
       WRITE(6,600)                                                 
       NR= NR+1
C   THE MOLECULAR FORMULA IS READ.                            
       CALL RECHA(SY,IA)                                            
       IF(IA.EQ.0) GOTO 9999                                        
       A=0                                                    
    4  IF(A.EQ.IA) GOTO 8                                           
       A=A+1                                                        
       K=SY(A)                                                      
       IF(K.NE.26) GOTO 4                                           
C   A RADICAL SYMBOL (Z) IS ENCOUNTERED   AND MUST BE TRANSLATED.         
       J1=MOD(SY(A+1),10)                                           
       J2=MOD(SY(A+2),10)                                           
       J=10*J1+J2                                                   
       IF(J.GT.0) GOTO 710
C   THE RADICAL REPRESENTS A SPECIAL ATOM, WHICH IS
C   INDICATED BY THE SYMBOL Z.
       IA=IA-2
       J1=A
       J2=A+2
       NJ=IA-A
       DO 700 I=1,NJ
  700  SY(J1+I)=SY(J2+I)
       GOTO 4
  710  CONTINUE
       HA=NA(J)                                                     
       IF(HA.EQ.3) GOTO 6                                           
       J1=IA+HA-2                                                   
       J2=IA+1                                                      
       NJ=IA-A-2                                                    
       DO 5 I=1,NJ                                                  
    5  SY(J1-I)=SY(J2-I)                                            
    6 A=A-1                                                   
      DO 7 I=1,HA                                                   
      A=A+1                                                   
    7 SY(A)=DD(J,I)                                                 
      IA=J1-1                                                       
      GOTO 4                                                        
    8 CONTINUE                                                      
C  THE MOLECULAR SYMBOL SEQUENCE HAS NOW BEEN STORED IN SY.               
      DO 15 K=1,65                                                  
      DO 15 J=1,65                                                  
      B(K,J)=0                                                      
      R(K,J)=0                                                      
   15 CONTINUE                                                      
   10 FORMAT(//,13H MOLECULE NR ,I3,//,18H NUMBER OF ATOMS =,I3,        
     118H    OUTPUT INDEX =,I2,//)                            
      READ(5,*) KM
      IF(KM.EQ.0) GOTO 18
C.....                                                              
C....  THE INPUT OF THE ANGLES 'REST' HAS BEEN CHANGED TO FREE             
C....  FORMAT.  ***** NOTE **** THE NUMBER OF ANGLES TO BE                
C....  INPUT IS A FUNCTION OF THE VALUES. THE INPUT IS FREE FORM        
C....  SO ALL THAT NEED BE ENTERED IS THE VALUES SEPARATED BY             
C....  SPACES OR COMMAS. EXAMPLE:----S                              
C...   180,0,65.3,180 180 0 0 60 60 45.67 33.33,76.5,0,0                   
C....                                                         
C....  THERE ARE 14 ANGLES ENTERED IN THE ABOVE                    
C....  EXAMPLE. OBVIOUSLY, THE NUMBER OF ANGLES ENTERED IS A              
C....  FUNCTION OF THE LENGTH. *****  NOTE **** ONLY THE FIRST 72        
C...   COLS OF THE CARD MAY BE USED. IF THIS IS NOT SUFFICIENT         
C....  THEN THE PROGRAM SHOULD BE CHANGED TO GO BACK TO THE               
C....  ORIGINAL METHOD OF INPUT WHICH IS GIVEN IN                   
C....  LINES 23200 AND 23300......  PAUL D GOODWIN.....                   
C.....                                                              
      READ(5,*)(REST(I),I=1,KM)                                     
   18 HH= 0
      K= 0                                                    
      H= 0                                                    
      IC= 1                                                   
      IZ=1                                                    
      IA= 0                                                   
      DO 20 J=1,9                                                   
   20 BBB(J)= 0                                                     
      N=65                                                    
      DO 22 J=1,65                                                  
      DO 21 I=1,3                                                   
   21 IR(J,I)=0                                                     
   22 BH(J)=0                                                       
      J= 1                                                    
   25 IF(J.GT.A) GOTO 35                                            
C  INTERPRETATION OF THE SYMBOL SEQUENCE  SY.                       
C  THE ARRAY  B(I,J)  BECOMES THE CONNECTIVITY MATRIX (0 = NO BOND);    
C  ARRAY  IR(K,J)  WILL CONTAIN FOR EACH ATOM NR  K:                
C    FOR  J = 1: THE NUMBER OF THE PRECEDING MAIN-CHAIN ATOM;             
C    FOR  J = 2: THE NUMBER OF THE NEXT   MAIN-CHAIN ATOM;                
C    FOR  J = 3: THE NUMBER OF THE FIRST SIDE-CHAIN ATOM.                 
      SJ= SY(J)                                                     
      IF(SJ.LT.12.OR.SJ.EQ.14.OR.SJ.EQ.15.OR.SJ.EQ.22.OR.SJ.GT.23)      
     1 GOTO 40                                                      
      K= K+1                                                        
      IF(SJ.NE.12) GOTO 41                                          
      S(K)= 28                                                      
      GOTO 50                                                       
   41 IF(SJ.NE.13) GOTO 42                                          
      S(K)= 1                                                       
      GOTO 50                                                       
   42 IF(SJ.NE.16) GOTO 43                                          
      S(K)= 26                                                      
      GOTO 50                                                       
   43 IF(SJ.NE.17) GOTO 44                                          
      S(K)= 10                                                      
      GOTO 50                                                       
   44 IF(SJ.NE.19) GOTO 45                                          
      S(K)= 11                                                      
      GOTO 50                                                       
   45 IF(SJ.NE.20) GOTO 46                                          
      S(K)= 20                                                      
      GOTO 50                                                       
   46 IF(SJ.NE.21) GOTO 47                                          
      S(K)= 23                                                      
      GOTO 50                                                       
   47 IF(SJ.NE.23) GOTO 48                                          
      S(K)=24                                                       
      GOTO 50                                                       
   48 IF(SJ.NE.18) S(K)=0                                           
      S(K)=30                                                       
   50 IF(IA.EQ.0) GOTO 60                                           
      HH= HH+1                                                      
      R(IC,K)= REST(HH)                                             
      IA= 0                                                   
   60 IF(K.EQ.1) GOTO 30                                            
      B(IC,K)= IZ                                                   
      IR(K,1)=IC                                                    
      IF(BH(IC).EQ.0) GOTO 400                                      
      BH(IC)=0                                                      
      IR(IC,3)=K                                                    
  400 IF(IR(IC,2).NE.0) GOTO 410                              
      IR(IC,2)=K                                                    
  410 CONTINUE                                                      
      IZ= 1                                                   
      IC= K                                                   
      GOTO 30                                                       
   40 IF(SJ.GT.10) GOTO 55                                          
      S(K)=S(K)+SJ                                                  
      IF(S(K).EQ.2) S(K)= 27                                        
      GOTO 30                                                       
   55 IF(SJ.NE.28) GOTO 70                                          
      H= H+1                                                        
      F(H)= K                                                       
      IC= K                                                   
      BH(IC)=1                                                      
      DO 420 JJ=1,IC                                                
      I=IR(JJ,2)                                                    
      IF(I.LT.0) IR(JJ,2)=I-1                                       
  420 CONTINUE                                                      
      IR(IC,2)=-1                                                   
      GOTO 30                                                       
   70 IF(SJ.NE.31) GOTO 80                                          
      IC= F(H)                                                      
      H= H-1                                                        
      DO 430 JJ=1,IC                                                
      I=IR(JJ,2)                                                    
      IF(I.LT.0) IR(JJ,2)=I+1                                       
  430 CONTINUE                                                      
      GOTO 30                                                       
   80 IF(SJ.NE.32) GOTO 85                                          
      IC= F(H)                                                      
      GOTO 30                                                       
   85 IF(SJ.NE.25) GOTO 90                                          
C  RING CLOSURE  X---X                                              
      J= J+1                                                        
      SJ=SY(J)                                                      
      IF(BBB(SJ).EQ.1) GOTO 100                                     
      BB(SJ)= IC                                                    
      BBB(SJ)= 1                                                    
      IF(BH(IC).EQ.0) GOTO 440                                      
       BH(IC)=0                                                     
       IR(IC,3)=-SJ                                                 
  440  CONTINUE                                                     
       GOTO 30                                                      
  100  I=BB(SJ)                                                     
       B(I,IC)=IZ                                                   
       IF(BH(IC).EQ.0) GOTO 450                                     
       BH(IC)=0                                                     
       IR(IC,3)=I                                                   
  450  CONTINUE                                                     
       IF(IR(I,3).EQ.-SJ) IR(I,3)=IC                                
       IZ= 1                                                        
       IF(IA.EQ.0) GOTO 110                                         
       HH= HH+1                                                     
       R(BB(SJ),IC)= REST(HH)                                       
       IA= 0                                                        
  110  BBB(SJ)= 0                                                   
       GOTO 30                                                      
   90  IF(SJ.EQ.22) IA=1                                            
C   SPECIAL BONDS  D, T, A OR E                                     
       IF(SJ.EQ.14) IZ=2                                            
       IF(SJ.EQ.24) IZ=3                                            
       IF(SJ.EQ.11) IZ=4                                            
       IF(SJ.EQ.15) IZ=5                                            
       IF(SJ.NE.26) GOTO 30
C   THE SYMBOL Z INDICATES A SPECIAL ATOM, WHICH
C   THE DATA ARE NOW PUT INTO THE APPROPRIATE ARRAYS.
C   ITS ATOM-TYPE NUMBER S(K) IS TAKEN =0  .
       K=K+1
       S(K)=0
       NV(K)=NNVZ
       AN(K)=ZDAT(1)
       DO 720 I=1,NNVZ
       IK=4*I-2
       ZV(I,K)=ZDAT(IK)
       DO 720 II=1,3
  720  V(I,II,K)=ZDAT(IK+II)
       GOTO 50
   30  J= J+1                                                       
       GOTO 25                                                      
   35  N=K                                                    
       IF(KM.EQ.0) GOTO 37                                          
       WRITE(6,10) NR,N,IPR                                         
       WRITE(6,12) (ICH(SY(I)), I=1,A)
   12  FORMAT(' ',120A1,/)
       WRITE(6,36) (REST(I), I= 1,KM)                                
   36  FORMAT(//,18H   TORSION ANGLES:,10(1X,F7.2),//)
C   - THE ATOMS ARE NOW CONNECTED -                           
C   IN THE NOW FOLLOWING FIRST SCAN OVER ALL ATOM PAIRS, EVERY            
C   POSITIVE  B(K,I) VALUE IS SUBSTITUTED BY THE REVERSE OF THE           
C   SEQUENCE NUMBER OF THE BOND VECTOR OF ATOM  K , WHICH                 
C   CONSTITUTES THE BOND  K--I.                                     
   37  DO 120 K= 1,N                                                
       DO 120 J= 1,3                                                
  120  C(K,J)= 0.0                                                  
       DO 130 K=1,N                                                 
C   FOR THE SPECIAL ATOM Z (WITH S(K)=0), DATPR IS NOT CALLED.
       IF (S(K).EQ.0) GOTO 130
       CALL DATPR(S(K))                                             
       NV(K)=NNV                                                    
       AN(K)=IAN                                                    
       DO 140 J=1,NNV                                               
       ZV(J,K)=ZZV(J)                                               
       DO 140 I=1,3                                                 
  140  V(J,I,K)=VV(J,I)                                             
  130  CONTINUE                                                     
       DO 150 K=1,N                                                 
       NNV=NV(K)                                                    
       IF(K.EQ.N) GOTO 220                                          
       KM=K+1                                                       
       DO 160 I=KM,N                                                
       BKI=B(K,I)                                                   
       IF(BKI.LE.0) GOTO 160                                        
       NVI=NV(I)                                                    
       DO 170 J=1,NNV                                               
       IF(ZV(J,K).NE.BKI) GOTO 170                            
       JB=1                                                   
      DO 180 H=1,N                                                  
      IF(B(K,H).EQ.-J) JB=0                                         
  180 CONTINUE                                                      
      IF(JB.EQ.0) GOTO 170                                          
      B(K,I)=-J                                                     
      GOTO 190                                                      
  170 CONTINUE                                                      
  190 DO 200 J=1,NVI                                                
      IF(ZV(J,I).NE.BKI) GOTO 200                             
      JB=1                                                    
      DO 210 H=1,N                                                  
      IF(B(I,H).EQ.-J) JB=0                                         
  210 CONTINUE                                                      
      IF(JB.EQ.0) GOTO 200                                          
      B(I,K)=-J                                                     
      GOTO 160                                                      
  200 CONTINUE                                                      
  160 CONTINUE                                                      
  150 CONTINUE                                                      
C   NOW FOLLOWS THE SECOND SCAN OVER ALL ATOM PAIRS, IN WHICH             
C   THE ARRAYS  B  AND  IR  PROVIDE THE   DATA FROM WHICH                 
C   THE BINDING VECTOR AND THE ORIENTATION VECTOR OF BOTH                 
C   ATOMS OF A BOUND PAIR  K, I  ARE BEING CALCULATED               
C   (BK, CK  AND BI, CI, RESPECTIVELY).    TOGETHER WITH THE              
C   KNOWN BOND-VECTOR SET  VV(1:NVI)  OF ATOM  I  THESE ARE               
C   FED INTO SUBROUTINE  DRAAI, WHICH GIVES THE PROPER                    
C   ORIENTATION TO THIS ATOM.                                       
  220 DO 230 K=1,N                                                  
      IF(K.EQ.N) GOTO 335                                           
      KM=K+1                                                        
      DO 240 I=KM,N                                                 
      IF(B(K,I).EQ.0) GOTO 240                                      
      NVI=NV(I)                                                     
      DO 250 J=1,3                                                  
      IRK=IR(K,J)                                                   
      IF(IRK.NE.0.AND.IRK.NE.I) GOTO 260                      
  250 CONTINUE                                                      
  260 DO 270 J=2,3                                                  
      IRI=IR(I,J)                                                   
      IF(IRI.NE.0.AND.IRI.NE.K) GOTO 280                      
  270 CONTINUE                                                      
  280 DO 290 J=1,3                                                  
      DO 300 H=1,NVI                                                
  300 VV(H,J)=V(H,J,I)                                              
      BK(J)=V(-B(K,I),J,K)                                          
      IF(IRK.EQ.0.OR.IRK.EQ.I) GOTO 302                             
      CK(J)=V(-B(K,IRK),J,K)                                        
      GOTO 304                                                      
  302 CK(J)=1.                                                      
  304 CONTINUE                                                      
      BI(J)=VV(-B(I,K),J)                                           
      IF(IRI.EQ.0.OR.IRI.EQ.K) GOTO 310                             
      CI(J)=VV(-B(I,IRI),J)                                         
      GOTO 290                                                      
  310 CI(J)=1.                                                      
  290 CONTINUE                                                      
      CALL DRAAI(NVI,VV,BK,CK,BI,CI,R(K,I))                         
      DO 312 H=1,NVI                                                
      DO 312 J=1,3                                                  
  312 V(H,J,I)=VV(H,J)                                              
      H= 1                                                    
       DO 320 J= 1,3                                                
       IF(ABS(C(I,J)).GT.0.0001) H= 2                               
  320  CONTINUE                                                     
C   CALCULATION OF THE POSITION OF THE NEW ATOM.                    
       DO 330 J= 1,3                                                
  330  C(I,J)=(C(I,J)+C(K,J)+BK(J)-VV(-B(I,K),J))/H                 
  240  CONTINUE                                                     
  230  CONTINUE                                                     
  335  CONTINUE                                                     
       WRITE(6,340)                                                 
  340  FORMAT(43H   VANDERWAALS RADII (*100) AND COORDINATES,//)
       DO 345 K= 1,N                                                
  345  WRITE(6,350) K,AN(K),(C(K,I),I= 1,3)                         
  350  FORMAT(2(3X,I3),3(3X,F6.2))
       IF(IPR.LE.1) GOTO 355                                        
C   TWO PROJECTIONS OF THE MOLECULE ARE   PLOTTED BY                
C   SUBROUTINE  PLOTAT .                                            
       CALL PLOTAT(1,2)                                             
       CALL PLOTAT(2,3)                                             
  355  WRITE(6,360)                                                 
  360  FORMAT(//)
  600  FORMAT(1H1)                                                  
       RETURN
 9999  CONTINUE
       END                                                    
       SUBROUTINE PLOTAT(A,B)                                       
       INTEGER A,B,DD,S,SO,SY,Y
       COMMON/CHAR/ICH(33)
       COMMON/ALG/NR,N,IPR,DD(99,100),NA(99),SY(500)                
       COMMON/DAT1/IDAT(184),VDAT(177)                              
       COMMON JAN(65),S(65),C(65,3),JB(65,65)                   
       DIMENSION Y(65,3)                                            
       REAL MI,MA
C   THIS SUBROUTINE GIVES A PROJECTION OF THE                       
C   MOLECULE IN THE PLAIN (A,B) (1, 2, 3 = X, Y, Z, RESPECTIVELY).      
C   THE SCALE IS INDICATED IN  MM  PER ANGSTROM UNIT.               
C   HALOGEN ATOMS ARE INDICATED BY THE LETTER  X.                   
C   FIRST THE SCALE AND THE SCALED (INTEGER) ATOMIC                 
C   COORDINATES  Y  ARE CALCULATED.                           
       J= A                                                   
   10  MI= 1000.                                                    
       MA= -MI                                                      
       DO 20 K= 1,N                                                 
       XX= C(K,J)                                                   
       IF(XX.GT.MA) MA= XX                                          
       IF(XX.LT.MI) MI= XX                                          
   20  CONTINUE                                                     
       IF(J.NE.A) GOTO 30                                           
       XO= MI                                                       
       BR= MA-MI                                                    
       J= B                                                   
       GOTO 10                                                      
   30  YO= MA                                                       
       HO= MA-MI                                                    
   50 IF (ABS(BR*HO).LT..1) GOTO 160                                
       MA= 300/BR                                                   
       MI= 246/HO                                                   
       IF(MA.GT.MI) MA= MI                                          
       WRITE(6,55) MA                                               
   55  FORMAT(8H1   1A = ,F6.1,3H MM,/)                             
       DO 60 K= 1,N                                                 
       Y(K,A)= 1+INT(.5+(C(K,A)-XO)*MA/2.54)                        
   60 Y(K,B)= 1+INT(.5+(YO-C(K,B))*MA/4.23)                         
C  FIND, FOR EACH LINE. IF AND WHERE ATOMS SHOULD BE PRINTED              
CCCCC The following line was formerly:  MA=ICH(27)
      MA=32
      DO 140 I= 1,59                                                
C  GENERATE AN EMPTY STRING  P.                                     
      DO 150 II=1,128                                               
  150 SY(II)=MA                                                     
      DO 85 K= 1,N                                                  
      IF(Y(K,B).NE.I) GOTO 85                                       
C  BACK TRANSLATION OF ATOMIC INDICES INTO THE EBCDIC CODE                
      SO= S(K)                                                      
      IF(SO.GT.0) GOTO 69
      IV=26
      GOTO 70
   69 CONTINUE
      IF(SO.GE.10.AND.SO.NE.13) GOTO 71                             
      IV=13                                                   
      GOTO 70                                                       
   71 IF(SO.NE.10) GOTO 72                                          
      IV=17                                                   
      GOTO 70                                                       
   72 IF(SO.GE.18) GOTO 73                                          
      IV=19                                                   
      GOTO 70                                                       
   73 IF(SO.GE.23) GOTO 74                                          
      IV=20                                                   
      GOTO 70                                                       
   74 IF(SO.NE.23) GOTO 75                                          
      IV=21                                                   
      GOTO 70                                                       
   75 IF(SO.GE.26) GOTO 76                                          
      IV=23                                                   
      GOTO 70                                                       
   76 IF(SO.GE.31) GOTO 77                                          
      IV=25                                                   
      GOTO 70                                                       
   77 IV= 0                                                   
      IF(IV.EQ.0) GOTO 85                                           
C  THE ATOMIC SYMBOLS  IV  ARE BEING PUT AT THE CORRECT                   
C  POSITIONS  Y(K,A)  OF THE STRING  P.                             
   70 SY(Y(K,A))=ICH(IV)                                            
   85 CONTINUE                                                      
C  THE FINISHED LINE IS PRINTED.                              
      WRITE(6,80)  (SY(II), II=1,128)                               
   80 FORMAT(1X,128A1)                                              
  140 CONTINUE                                                      
  160 RETURN                                                        
      END                                                     
      SUBROUTINE PARAM                                              
C  THIS SUBROUTINE CALCULATES THE STERIC PARAMETERS OF THE                
C  MOLECULE, USING THE ATOMIC COORDINATES  C  AND THE               
C  VANDERWAALS RADII AN.
      INTEGER D,DD,S
      COMMON/ALG/NR,N,IPR,DD(99,100),NA(99),SY(500)                 
      COMMON JAN(65),S(65),C(65,3),JB(65,65),D(65,65)                 
      DIMENSION XAN(65),B(4),BB(4)
C*    THE FIRST STEP IS TO MEASURE THE LENGTH -- THIS IS DONE
C*    BY GOING THROUGH THE COORDINATES IN THE X DIRECTION 
C*    MULTIPLY BY THE VAN DER WAALS RADIUS (JAN(K)*0.01) THEN
C*    ASSIGNING EACH AS ONE GOES TO AN ADDRESS AMA OR AMI AND 
C*    COMPARING UNTIL THE LARGEST ONE IS FOUND. ACTUAL DISTANCE
C*    IS ASSIGNED TO AL
      DO 2 K=1,N                                                    
    2 XAN(K)=JAN(K)*0.01                                                  
      WRITE(6,5)                                                    
    5 FORMAT(1X)                                                    
      AMI=1000.                                                     
      AMA=-1000.                                                    
      DO 10 K=2,N                                                   
      CZ=C(K,1)                                               
      RW=XAN(K)                                                     
       A=CZ+RW                                                      
       IF(A.GT.AMA) AMA=A                                           
       A=CZ-RW                                                      
       IF(A.LT.AMI) AMI=A                                           
   10  CONTINUE                                                     
C   AL  IS THE LENGTH ALONG THE  X  AXIS, MEASURED FROM THE               
C   CENTRE OF THE FIRST ATOM, WHICH SHOULD BE A                     
C   HYDROGEN.  THE LENGTH IS CORRECTED BY ADDING  0.40 A TO               
C   BECOME THE ORIGINAL STANDARD VALUE,   WHICH WAS TAKEN FROM A C6       
C   ATOM AS A STARTING POINT.                                       
       AL=-AMI+0.4                                                  
       AMI=1000.                                                    
       AMA=-1000.                                                   
C   ROTATION OF THE SUBSTITUENT AROUND THE  X AXIS,                 
C   WITH STEPS OF 1.72 DEGREES, OVER SOMEWHAT MORE                  
C   THAN  90 DEGREES.                                               
       DO 50 I=1,53                                                 
       SI=SIN(0.03*I)                                               
       CO=COS(0.03*I)                                               
       B(1)=1000.                                                   
       B(2)=-1000.                                                  
       B(3)=1000.                                                   
       B(4)=-1000.                                                  
C   AT EACH POSITION, THE FOUR WIDTHS IN PERPENDICULAR                    
C   DIRECTIONS ARE MEASURED; AT THE POSITION WHERE                  
C   THE SMALLEST WIDTH IS ENCOUNTERED, THE FOUR VALUES ARE                
C   TAKEN TO BECOME THE B(1 - 4) VALUES.                      
       DO 20 K=2,N                                                  
       RW=XAN(K)                                                    
       CZ=C(K,2)*CO + C(K,3)*SI                                   
       A=CZ+RW                                                      
       IF(A.GT.B(2)) B(2)=A                                         
       A=CZ-RW                                                      
       IF(A.LT.B(1)) B(1)=A                                         
       CZ=C(K,3)*CO - C(K,2)*SI                                    
       A=CZ+RW                                                      
       IF(A.GT.B(4)) B(4)=A                                         
       A=CZ-RW                                                      
       IF(A.LT.B(3)) B(3)=A                                         
   20  CONTINUE                                                     
C   THE LARGEST WIDTH IS ALSO CALCULATED (AMA).  IT IS PRINTED            
C   AS THE PARAMETER B(5).                                          
       DO 40 K=1,4                                                  
       ABSBK=ABS(B(K))                                              
       IF(ABSBK.GT.AMA) AMA=ABSBK                             
       IF(ABSBK.GE.AMI) GOTO 40                                     
       AMI=ABSBK                                                    
       DO 30 KK=1,4                                                 
   30  BB(KK)=ABS(B(KK))                                            
   40  CONTINUE                                                     
   50  CONTINUE                                                     
       WRITE(6,5)                                                   
       WRITE(6,5)                                                   
       WRITE(6,60) AL,(BB(J), J=1,4),AMA                      
   60 FORMAT ('  STERIC PARAMETERS:   L=',F6.2//
     1 '  B(1) - B(4): ' ,4(3X,F6.2),'   B(5):',3X,F6.2)
       WRITE(6,5)                                                   
C   FINALLY, IN ORDER TO AVOID IMPOSSIBLE CONFORMATIONS,
C   THE DISTANCES BETWEEN ATOMS WITH THREE                          
C   OR MORE BONDS BETWEEN THEM ARE CALCULATED.  IF OVERLAPS               
C   LARGER THAN 0.05 A ARE FOUND, A WARNING IS PRINTED.                   
      WRITE(6,5)                                                    
      WRITE(6,5)                                                    
      INF=10000                                                     
      DO 110 I=1,N                                                  
      DO 100 J=1,N                                                  
  100 D(I,J)=INF                                                    
  110 D(I,I)=0                                                      
      DO 120 K=1,N                                                  
      JB(K,K)=0                                                     
      DO 120 I=K,N                                                  
      IF(JB(K,I).EQ.0) GOTO 120                                     
      D(I,K)=1                                                      
      D(K,I)=1                                                      
  120 CONTINUE                                                      
      DO 131 I=1,N
      DO 131 J=1,N
      IF(D(I,J).GE.INF) GOTO 131
      IS=D(I,J)      
      DO 132 K=1,N
      ISS=IS+D(J,K)                                                 
      IF(ISS.GE.D(I,K)) GOTO 132
      D(K,I)=ISS                                                    
      D(I,K)=ISS                                                    
  132 CONTINUE
  131 CONTINUE
      IB=0                                                    
      DO 90 K=1,N                                                   
      DO 90 I=K,N                                                   
      IF(D(K,I).LT.3) GOTO 90                                       
      AMA=0.                                                        
      DO 70 J=1,3                                                   
   70 AMA=AMA+(C(K,J)-C(I,J))**2                              
      AMA=SQRT(AMA)                                                 
      AMI=XAN(K)+XAN(I)-AMA                                       
      IF(AMI.LT.0.05) GOTO 90                                       
      IF(IB.EQ.1) GOTO 80                                           
      WRITE(6,75)                                                   
   75 FORMAT(51H  TAKE CARE; MUTUAL PENETRATION OF FOLLOWING ATOMS:,/)
      IB=1                                                    
   80 WRITE(6,85) K,I,AMI                                           
   85 FORMAT(' NRS. ',I2,' AND ',I2,'   ; OVERLAP IS ',F6.2)
   90 CONTINUE                                                      
      RETURN                                                        
      END                                                     
      SUBROUTINE TRANSP                                             
C   THIS SUBROUTINE TRANSPOSES THE ATOMS, ORDERING THEM ACCORDING       
C   TO INCREASING  X COORDINATE VALUES.    THIS MAKES THE IDENTIFICATION  
C   OF THE SYMBOLS OF THE MOLECULE PLOTS EASIER.                    
      INTEGER AN,S,DD
      COMMON/ALG/NR,N,IPR,DD(99,100),NA(99),SY(500)                 
      COMMON AN(65),S(65),C(65,3),B(65,65),R(65,65)                   
      J=0
    2 IF(J.EQ.N-1) GOTO 4                                           
      J=J+1                                                   
      SS=1000.                                                      
      DO 1 K=J,N                                                    
      IF(J.EQ.1) S(K)=K                                             
      IF(C(K,1).GE.SS) GOTO 1                                       
      KK=K                                                    
      SS=C(K,1)                                                     
    1 CONTINUE                                                      
      IF(KK.EQ.J) GOTO 2                                            
      K=AN(J)                                                       
      AN(J)=AN(KK)                                                  
      AN(KK)=K                                                      
      K=S(J)                                                        
      S(J)=S(KK)                                                    
      S(KK)=K                                                       
      DO 3 I=1,3                                                    
      SS=C(J,I)                                                     
      C(J,I)=C(KK,I)                                                
    3 C(KK,I)=SS                                                    
      GOTO 2                                                        
    4 WRITE(6,600)                                                  
      DO 5 K=1,N                                                    
    5 WRITE(6,601) K,S(K),AN(K),(C(K,I), I=1,3)                     
      WRITE(6,602)                                                  
      RETURN                                                        
  600 FORMAT(44H1  ATOMS ORDERED BY INCREASING X COORDINATES,///)       
  601 FORMAT(3(3X,I3),3(3X,F6.2))                             
  602 FORMAT(1H1)                                                   
      END                                                     
      SUBROUTINE RECHA(KH,II)                                       
      COMMON/CHAR/ICH(33)
      DIMENSION KH(500),KA(72)                                      
C..........................XX.....                            
      II=0                                                    
C.....                                                              
C....   THIS INPUT ROUTINE HAS BEEN CHANGED SO THAT                 
C....   THE FORMULAE CAN ONLY BE ENTERED IN THE                     
C....   FIRST 72 COLUMNS OF THE DATA CARDS.   THE LAST 8                  
C....   COLUMNS ARE USED TO HOLD THE LINE NUMBERS                   
C....   WITH THE RESULT BEING AN EASIER   TO EDIT DATA DECK.              
C.....          THE CHANGE TO THE ORIGINAL PROGRAM ENTAIL                 
C....   CHANGING THE 80 TO A 72 IN THE RIGHT PLACES.                
C.....         THESE CHANGES HAVE BEEN INDICATED                          
C....   WITH A XX UNDER THE VALUES....                              
C....                                                         
    1 READ(5,500,END=6) (KA(J), J=1,72)                             
C...................................XX........                      
      DO 5 J=1,72                                                   
C..............XX..........                                         
      KAJ=KA(J)                                                     
      DO 2 I=1,33                                                   
      IF(KAJ.EQ.ICH(I)) GOTO 3                                      
    2 CONTINUE                                                      
    3 IF(I.EQ.30) GOTO 6                                            
      IF(I.EQ.29) GOTO 1                                            
      IF(I.EQ.27) GOTO 4                                            
      II=II+1                                                       
      KH(II)=I                                                      
      GOTO 5                                                        
    4 IF(J.EQ.72) GOTO 6                                            
C.............XX...............                                     
    5 CONTINUE                                                      
    6 IF(II.EQ.0) RETURN                                            
      IF(I.EQ.27) GOTO 7                                            
      WRITE(6,600) (ICH(KH(I)), I=1,II)                             
      RETURN                                                        
    7 WRITE(6,601)                                                  
  500 FORMAT(72A1)
C............XX........                                             
  600 FORMAT(1H ,120A1)                                             
  601 FORMAT(//,13H  WRONG INPUT)                             
      END
