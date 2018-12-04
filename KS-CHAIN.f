C         AN N-Body PROGRAM WITH CHAIN INTEGRATION. External potential can be included (look EXTERNALU and EXTERNAL ACCELERATIONS) 
          IMPLICIT REAL*8 (A-H,M,O-Z)
          COMMON/DIAGNOSTICS/IWR,IWK,H,GAMMA
           REAL*8 G0(3),G(3),X(100),V(100),M(25)
          COMMON/INTEGR/JMAX
         LOGICAL NEWREG
          TIME=0.0
          III=123456789
C         Initial values
C        iwr =writing index (normally =0, if iwr>0, then some diagno output)
c        N= number of bodies
c        DELTAT = output interval (approximate, may be longer if DELTAT is small)
c        TMAX= maximum time
c        EPS= error tolerance (1.e-13 recommended for accurate computation)
c        KSMX=MAximum number of steps between outputs (use = 1 for frequent output)

         READ(5,*)IWR, N,DELTAT,TMAX,EPS,KSMX
         JMAX=MIN(12,3+KSMX) !if KSMX small (1 or 2 or..) => short steps because of this [modify this as U want)
                             !if U want longer steps even with small KSMX, set JMAX=10  
         OPEN(66,file='C') ! output of time and x y z to this file
          MASS=0
          DO I=1,N ! read mass, x, y ,z, vx, vy, vz for body I
          K0=3*I-3
          READ(5,*)M(I),(X(K0+K),K=1,3),(V(K0+K),K=1,3)
         MASS=MASS+M(I)
          END DO
          call   CONSTANTS OF MOTION(X,V,M,N,ENER0,G0,AL)
          TCR=MASS**2.5/(2.*ABS(ENER0))**1.5
           ANGMO=ABS(ENER0)*TCR
           WRITE(6,*)' ENERGY',ENER0
          NEWREG=.TRUE.
100       CONTINUE
          call CHAIN(N,X,V,M,TIME,DELTAT,NEWREG,KSMX,EPS)
c        Diagnostics
         call   CONSTANTS OF MOTION(X,V,M,N,ENER1,G,AL)
         WRITE(6,123)TIME,(ENER0-ENER1)/AL,((G(J)-G0(J))/ANGMO,J=1,3)! energy check & angular momentum check
123      FORMAT(1x,'T:',F12.2,' dE/L:',1pg9.1,' dRxV:',1p3g9.1)
         write(66,166)TIME,(X(I),I=1,3*N)
166      format(1x,1p,600g12.4)
        IF(TIME.LT.TMAX)GOTO 100
        END
        SUBROUTINE CHAIN(NN,XX,VX,MX,TIME,DELTAT,NEWREG,KSMX,EPS)
C        CHAIN INTEGRATION. Perturbations & CM-motion included.
         INCLUDE 'COMMON1.CH'
         INCLUDE 'COMMON2.CH'
         COMMON/DIAGNOSTICS/IWR,IWK,H,GAMMA
        REAL*8 G0(3),Y(NMX8),XX(1),VX(1),MX(1)
         LOGICAL MUSTSWITCH,NEWREG
        EXTERNAL DERIVATIVES
        SAVE
c          Initial constants of motion
          NEQ=8*NN
          CHTIME=0.0
           Y(NEQ)=CHTIME
           IF(NEWREG)THEN
           N=NN
          DO I=1,N
          M(I)=MX(I)
           END DO
          DO I=1,3*N
           X(I)=XX(I)
          V(I)=VX(I)
          END DO
           call CONSTANTS OF MOTION(X,V,M,N,ENER0,G0,ALAG)
           call FIND CHAIN INDICES
           call EVALUATE Q and P
C          Remove XTRNLU & UG if external potential is absent.
           call TAKE Y FROM COMMON (Y)
           IF(STEP.EQ.0.0)STEP=MASS**2.5/ALAG**0.5*EPS**0.2
          STIME=0.0
           NEWREG=.FALSE.
          END IF
          KSTEPS=0
777        KSTEPS=KSTEPS+1
           oldstep=step
c*********************************
c                          IWK=0
c*********************************
           call  DIFSY1(NEQ,DERIVATIVES,EPS,STEP,STIME,Y)
           if(step.gt.10.*oldstep)step=10.*oldstep
           if(step.eq.0.0)then
          write(6,*)' stepsize=0!',char(7)
                                          STOP
           end if
            call CHECK SWITCHING CONDITIONS(MUST SWITCH)
            IF(MUST SWITCH)THEN
            call SWITCH(Y)
                         IF(IWR.GT.0)THEN
                         WRITE(6,1232)(INAME(KW),KW=1,N)
1232                        FORMAT(1X,' I-CHAIN',20I3)
                         END IF
            END IF
           ! IF(IWR.GT.0)WRITE(6,1299)KSTEPS,KSMX,STEP,Y(NEQ),GAMMA
1299           FORMAT(1X,' STEPS:',2I11,2F9.1,1PE10.1)
                 CHTIME=Y(NEQ) 
                 IF((CHTIME.GT.DELTAT).OR.(KSTEPS.GT.KSMX))THEN
                call PUT Y TO COMMON(Y)
                 call EVALUATE X AND V
                 DO I=1,3*N
                 XX(I)=X(I)
                VX(I)=V(I)
                 END DO
                 TIME=TIME+CHTIME
                 RETURN
                ELSE
                 GOTO 777
                 END IF
        END
        SUBROUTINE CHECK SWITCHING CONDITIONS(MUSTSWITCH)
        INCLUDE 'COMMON1.CH'
        LOGICAL MUSTSWITCH
        DATA Ncall,NSWITCH/0,200/
        SAVE
        MUST SWITCH=.FALSE.
        Ncall=Ncall+1
C       Switch anyway after every NSWITCHth step.
        IF(Ncall.GE.NSWITCH)THEN
        Ncall=0
        MUST SWITCH=.TRUE.
        RETURN
        END IF
C       Inspect the structure of the chain.
C       NOTE: Inverse values 1/r are used instead of r's itself.
        ADISTI=0.5*(N-1)/RSUM
        LRI=N-1
        DO I=1,N-2
       DO J=I+2,N
       LRI=LRI+1
C       Do not inspect if 1/r is small.
        IF(RINV(LRI).GT.ADISTI)THEN
         IF(J-I.GT.2)THEN
C        Check for a dangerous long loop.
C          RINVMX=MAX(RINV(I-1),RINV(I),RINV(J-1),RINV(J))
           IF(I.GT.1)THEN
           RINVMX=MAX(RINV(I-1),RINV(I))
           ELSE
           RINVMX=RINV(1)
           END IF
           RINVMX=MAX(RINVMX,RINV(J-1))
           IF(J.LT.N)RINVMX=MAX(RINVMX,RINV(J))
           IF(RINV(LRI).GT.RINVMX)THEN
           MUST SWITCH=.TRUE.
           Ncall=0
           RETURN
           END IF
         ELSE
C        Is this a triangle with smallest size not regularised?
           IF( RINV(LRI).GT.MAX(RINV(I),RINV(I+1)))THEN
           MUST SWITCH=.TRUE.
          Ncall=0
           RETURN
          END IF
         END IF
        END IF
       END DO
       END DO
        RETURN
        END
        SUBROUTINE SWITCH(Y)
        INCLUDE 'COMMON1.CH'
       INCLUDE 'COMMON2.CH'
       REAL*8 Y(1)
       SAVE
       call PUT Y TO COMMON (Y)
       call CHAIN TRANSFORMATION
        call TAKE Y FROM COMMON (Y)
       RETURN
       END
       SUBROUTINE TAKE Y FROM COMMON (Y)
       INCLUDE 'COMMON1.CH'
        INCLUDE 'COMMON2.CH'
       REAL*8 Y(1)
       SAVE
       NC=N-1
        L=0
        DO I=1,4*NC
        L=L+1
       Y(L)=Q(I)
       END DO
        DO I=1,3
       L=L+1
        Y(L)=CMX(I)
        END DO
        L=L+1
       Y(L)=ENERGY
       DO I=1,4*NC
        L=L+1
       Y(L)=P(I)
       END DO
       DO I=1,3
        L=L+1
       Y(L)=CMV(I)
       END DO
        L=L+1
       Y(L)=CHTIME
        RETURN
       END
         SUBROUTINE FIND CHAIN INDICES
         INCLUDE 'COMMON1.CH'
        INCLUDE 'COMMON2.CH'
        REAL*8 RIJ2(NMXM)
       INTEGER IC(NMX2),IJ(NMXM,2),IND(NMXM)
       LOGICAL USED(NMXM),SUC,LOOP
       SAVE
        L=0
        DO I=1,N-1
       DO J=I+1,N
       L=L+1
       RIJ2(L)=SQUARE(X(3*I-2),X(3*J-2))
        IJ(L,1)=I
       IJ(L,2)=J
       USED(L)=.FALSE.
        END DO
       END DO
        call HEAPSORT(L,RIJ2,IND)
        LMIN=1+NMX
       LMAX=2+NMX
       IC(LMIN)=IJ(IND(1),1)
       IC(LMAX)=IJ(IND(1),2)
       USED(IND(1))=.TRUE.
1       DO I=2,L
       LI=IND(I)
       IF( .NOT.USED(LI))THEN
       call CHECK CONNECTION(IC,LMIN,LMAX,IJ,LI,SUC,LOOP)
        IF(SUC)THEN
       USED(LI)=.TRUE.
       GOTO 2
       ELSE
        USED(LI)=LOOP
       END IF
       END IF
       END DO
2       IF(LMAX-LMIN+1.LT.N)GO TO 1
       L=0
       DO I=LMIN,LMAX
       L=L+1
       INAME(L)=IC(I)
       END DO
       RETURN
       END
       SUBROUTINE CHECK CONNECTION(IC,LMIN,LMAX,IJ,LI,SUC,LOOP)
        INCLUDE 'COMMON1.CH'
        INCLUDE 'COMMON2.CH'
       INTEGER IC(1),ICC(2),IJ(NMXM,2)
       LOGICAL SUC,LOOP
       SAVE
       SUC=.FALSE.
        LOOP=.FALSE.
       ICC(1)=IC(LMIN)
       ICC(2)=IC(LMAX)
       DO I=1,2
       DO J=1,2
       IF(ICC(I).EQ.IJ(LI,J))THEN
        JC=3-J
        LOOP=.TRUE.
       DO L=LMIN,LMAX
       IF(IC(L).EQ.IJ(LI,JC))RETURN
       END DO
       SUC=.TRUE.
        LOOP=.FALSE.
       IF(I.EQ.1)THEN
       LMIN=LMIN-1
       IC(LMIN)=IJ(LI,JC)
       RETURN
       ELSE
       LMAX=LMAX+1
       IC(LMAX)=IJ(LI,JC)
       RETURN
       END IF
       END IF
       END DO
       END DO
       RETURN
       END
        SUBROUTINE HEAPSORT(N,Array,Indx)
       implicit real*8 (a-h,o-z)
        dimension Array(*),Indx(*)
        SAVE
        do 11 j=1,N
        Indx(j)=J
11      continue
        if(N.lt.2)RETURN
        l=N/2+1
        ir=N
10      CONTINUE
        IF(l.gt.1)THEN
        l=l-1
        Indxt=Indx(l)
        q=Array(Indxt)
        ELSE
        Indxt=Indx(ir)
        q=Array(Indxt)
        Indx(ir)=Indx(1)
        ir=ir-1
        IF(ir.eq.1)THEN
        Indx(1)=Indxt
        RETURN
        END IF
        END IF
        i=l
        j=l+l
20      IF(j.le.ir)THEN
            IF(j.lt.ir)THEN
               IF(Array(Indx(j)).lt.Array(Indx(j+1)))j=j+1
            END IF
            IF(q.lt.Array(Indx(j)))THEN
               Indx(i)=Indx(j)
               i=j
               j=j+j
            ELSE
               j=ir+1
            END IF
         GOTO 20
         END IF
         Indx(i)=Indxt
         GO TO 10
         END
        SUBROUTINE EVALUATE Q AND P
        INCLUDE 'COMMON1.CH'
        INCLUDE 'COMMON2.CH'
        REAL*8 GA(3)
        SAVE
C       Center of mass
        DO K=1,3
       CMX(K)=0.0
       CMV(K)=0.0
       END DO
       MASS=0.0
        DO I=1,N
        L=3*(I-1)
       MC(I)=M(INAME(I))
       MASS=MASS+MC(I)
        DO K=1,3
        CMX(K)=CMX(K)+M(I)*X(L+K)
        CMV(K)=CMV(K)+M(I)*V(L+K)
        END DO
       END DO
       DO K=1,3
       CMX(K)=CMX(K)/MASS
       CMV(K)=CMV(K)/MASS
       END DO
C       AUXILIARY QUANTITIES
        DO I=1,N-1
        TKK(I)=.5D0*(1./MC(I)+1./MC(I+1))
       TK1(I)=-1./MC(I)
       MKK(I)=MC(I)*MC(I+1)
       DO J=I+1,N
        MIJ(I,J)=MC(I)*MC(J)
       MIJ(J,I)=MIJ(I,J)
       END DO
       END DO
        call CONSTANTS OF MOTION(X,V,M,N,ENER0,GA,ALAG)
        call EXTERNALU(X,M,N,UG)
          UG=0.0
           ENERGY=ENER0
C     &  -.5D0*MASS*(CMV(1)**2+CMV(2)**2+CMV(3)**2)+UG
C       Pysical momenta
       DO I=1,N
       L=3*(I-1)
        LF=3*INAME(I)-3
        DO K=1,3
         XI(L+K)=X(LF+K)
        PI(L+K)=M(INAME(I))*(V(LF+K)-CMV(K))
        END DO
       END DO
C       Chain momenta
       L=3*(N-2)
        DO K=1,3
        WC(K)=-PI(K)
       WC(L+K)=PI(L+K+3)
       END DO
       DO I=2,N-2
       L=3*(I-1)
       DO K=1,3
       WC(L+K)=WC(L+K-3)-PI(L+K)
       END DO
       END DO
C       Chain coordinates
        DO I=1,N-1
       L=3*(I-1)
       DO K=1,3
       XC(L+K)=XI(L+K+3)-XI(L+K)
       END DO
       END DO
C       KS-transformations
       DO I=1,N-1
       L1=3*(I-1)+1
        KS1=4*(I-1)+1
       call PH TO KS (XC(L1),WC(L1),Q(KS1),P(KS1))
       END DO
       RETURN
       END
       SUBROUTINE EVALUATE X AND V
       INCLUDE 'COMMON1.CH'
        INCLUDE 'COMMON2.CH'
       REAL*8 Q0(3)
       SAVE
C       Transformation of the KS-variables to the physical ones.
C       First transform to chain coordinates.
       DO I=1,N-1
       L1=3*(I-1)+1
       KS1=4*(I-1)+1
       call KS TO PH(Q(KS1),P(KS1),XC(L1),WC(L1))
       END DO
C       Obtain physical variables from chain quantities.
       L=3*(N-2)
       DO K=1,3
       PI(K)=-WC(K)
       PI(L+K+3)=WC(L+K)
       END DO
       DO I=2,N-1
       L=3*(I-1)
       DO K=1,3
       PI(L+K)=WC(L+K-3)-WC(L+K)
       END DO
       END DO
       DO K=1,3
       XI(K)=0.0
        Q0(K)=0.0
       END DO
       DO I=1,N-1
       L=3*(I-1)
       DO K=1,3
       XI(L+3+K)=XI(L+K)+XC(L+K)
       END DO
       END DO
        DO I=1,N
        L=3*(I-1)
       DO K=1,3
       Q0(K)=Q0(K)+XI(L+K)*MC(I)/MASS
       END DO
       END DO
C       Rearrange according to INAME(i) and add CM.
       DO I=1,N
       L=3*(I-1)
       LF=3*(INAME(I)-1)
       DO K=1,3
       X(LF+K)=XI(L+K)-Q0(K)+CMX(K)
       V(LF+K)=PI(L+K)/MC(I)+CMV(K)
       END DO
       END DO
       RETURN
       END
       SUBROUTINE CHAIN TRANSFORMATION
       INCLUDE 'COMMON1.CH'
        INCLUDE 'COMMON2.CH'
       REAL*8 XCNEW(NMX3)
        INTEGER IOLD(NMX)
        SAVE
C       Transformation of the chain.
C       First transform to chain coordinates.
       DO I=1,N-1
       L1=3*(I-1)+1
       KS1=4*(I-1)+1
       call KS TO PH(Q(KS1),P(KS1),XC(L1),WC(L1))
       END DO
       L2=3*(INAME(1)-1)
       DO K=1,3
       X(L2+K)=0.0
       END DO
C       X's are needed when determining new chain indices.
       DO I=1,N-1
       L=3*(I-1)
        L1=L2
        L2=3*(INAME(I+1)-1)
       DO K=1,3
       X(L2+K)=X(L1+K)+XC(L+K)
       END DO
       END DO
C       Store the old chain indices.
        DO I=1,N
       IOLD(I)=INAME(I)
       END DO
C       Find new ones.
       call FIND CHAIN INDICES
C       Transform chain momenta
       L1=3*(IOLD(1)-1)
        LN=3*(IOLD(N)-1)
        L=3*(N-2)
       DO K=1,3
       PI(L1+K)=-WC(K)
       PI(LN+K)=WC(L+K)
       END DO
       DO I=2,N-1
       L=3*(I-1)
        LI=3*(IOLD(I)-1)
       DO K=1,3
       PI(LI+K)=WC(L+K-3)-WC(L+K)
       END DO
       END DO
        L1=3*(INAME(1)-1)
        LN=3*(INAME(N)-1)
       L=3*(N-2)
        DO K=1,3
        WC(K)=-PI(L1+K)
       WC(L+K)=PI(LN+K)
       END DO
       DO I=2,N-2
       L=3*(I-1)
        LI=3*(INAME(I)-1)
       DO K=1,3
       WC(L+K)=WC(L+K-3)-PI(LI+K)
       END DO
       END DO
C       Construct new chain coordinates. Transformation matrix
C       (from old to new) has only coefficients -1, 0 or +1.
        DO I=1,3*(N-1)
       XCNEW(I)=0.0
       END DO
        DO ICNEW=1,N-1
C       Obtain K0 &  K1 such that iold(k0)=iname(icnew)
c                                 iold(k1)=iname(icnew+1)
       LNEW=3*(ICNEW-1)
        DO I=1,N
        IF(IOLD(I).EQ.INAME(ICNEW))K0=I
        IF(IOLD(I).EQ.INAME(ICNEW+1))K1=I
       END DO
        DO ICOLD=1,N-1
        LOLD=3*(ICOLD-1)
        IF( (K1.GT.ICOLD).AND.(K0.LE.ICOLD))THEN
C       ADD
        DO K=1,3
        XCNEW(LNEW+K)=XCNEW(LNEW+K)+XC(LOLD+K)
       END DO
        ELSEIF( (K1.LE.ICOLD).AND.(K0.GT.ICOLD) )THEN
C       SUBTRACT
        DO K=1,3
       XCNEW(LNEW+K)=XCNEW(LNEW+K)-XC(LOLD+K)
       END DO
        END IF
       END DO
        END DO
C       KS-transformations
       DO I=1,N-1
       L1=3*(I-1)+1
        KS1=4*(I-1)+1
       call PH TO KS (XCNEW(L1),WC(L1),Q(KS1),P(KS1))
       END DO
C       Auxiliary quantities.
       MASS=0.0
        DO I=1,N
        L=3*(I-1)
       MC(I)=M(INAME(I))
       MASS=MASS+MC(I)
        END DO
        DO I=1,N-1
        TKK(I)=.5D0*(1./MC(I)+1./MC(I+1))
       TK1(I)=-1./MC(I)
       MKK(I)=MC(I)*MC(I+1)
       DO J=I+1,N
        MIJ(I,J)=MC(I)*MC(J)
        MIJ(J,I)=MIJ(I,J)
       END DO
       END DO
       RETURN
       END
       SUBROUTINE PUT Y TO COMMON (Y)
       INCLUDE 'COMMON1.CH'
        INCLUDE 'COMMON2.CH'
       REAL*8 Y(*)
       SAVE
       NC=N-1
        L=0
        DO I=1,4*NC
       L=L+1
       Q(I)=Y(L)
       END DO
        DO I=1,3
       L=L+1
        CMX(I)=Y(L)
        END DO
       L=L+1
       ENERGY=Y(L)
       DO I=1,4*NC
       L=L+1
       P(I)=Y(L)
       END DO
       DO I=1,3
       L=L+1
       CMV(I)=Y(L)
       END DO
       L=L+1
       CHTIME=Y(L)
        RETURN
       END
        SUBROUTINE DERIVATIVES(Y,D)
        INCLUDE 'COMMON1.CH'
        SAVE
        REAL*8 Y(*),D(*)
        NC=N-1
        LQ=1
       LX=4*NC+1
        LE=4*NC+4
        LP=4*NC+5
        LV=8*NC+5
        LT=8*NC+8
        call DERIVATIVES OF CHAIN VARIABLES
     &            (Y(LQ),Y(LX),Y(LE),Y(LP),Y(LV),Y(LT)
     &            ,D(LQ),D(LX),D(LE),D(LP),D(LV),D(LT))
        RETURN
       END
       SUBROUTINE DERIVATIVES OF CHAIN VARIABLES
     &  (Q,CMX,ENERGY,P,CMV,CHTIME,DQ,DX,DE,DP,DV,DT)
        INCLUDE 'COMMON1.CH'
C       INCLUDE 'COMMON2.CH'
        COMMON/DIAGNOSTICS/IWR,IWK,H,GAMMA
        SAVE
        REAL*8 Q(*),CMX(*),P(*),CMV(*)
       REAL*8 DQ(*),DX(3),DP(*),DV(3)
       REAL*8 W(NMX4),AK(NMX4),DK(NMX),FNC(NMX3),FXTNL(NMX3)
        REAL*8 TP(NMX4),TQ(NMX4),UQ(NMX4),FAUX(4),AUX(-1:1),XAUX(3)
        REAL*8 FCM(3)
C        M&A 1990 EQS. (77)->(80)
         UC=0.0
        RSUM=0.0
C        (77)
         DO I=1,N-1
         L=4*(I-1)
         call QxP(Q(L+1),P(L+1),W(L+1))
         RIJL=Q(L+1)**2+Q(L+2)**2+Q(L+3)**2+Q(L+4)**2  ! +eps*eps
C        Evaluate RSUM for decisionmaking.
         RSUM=RSUM+RIJL
         RINV(I)=1./RIJL
         A=.5D0*RINV(I)
         UC=UC+MKK(I)*RINV(I)
         DO K=1,4
         W(L+K)=A*W(L+K)
        END DO
         END DO
         LRI=N-1
C        (78)
         TKIN=.5D0*MASS*(CMV(1)**2+CMV(2)**2+CMV(3)**2)
         DO I=1,N-1
        J1=-1
        J2=+1
         IF(I.EQ.1)J1=0
         IF(I.EQ.N-1)J2=0
         AUX(-1)=0.5d0*TK1(I)
        AUX( 0)=TKK(I)
         AUX(+1)=0.5d0*TK1(I+1)
         L=4*(I-1)
         DK(I)=0.0
           DO K=1,4
          AA=0.0
             DO J=J1,J2
             LJ=L+4*J
            AA=AA+AUX(J)*W(LJ+K)
            END DO
           AK(L+K)=AA
C           (79)
           DK(I)=DK(I)+AA*W(L+K)
          END DO
C           (80)
       TKIN=TKIN+DK(I)
       END DO
C       Physical coordinates
        DO K=1,3
       XI(K)=0.0
       END DO
       DO I=1,N-1
        L=3*(I-1)
       L1=3*(I-1)+1
       KS1=4*(I-1)+1
       call QxQ(Q(KS1),XC(L1))
       DO K=1,3
       XI(L+3+K)=XI(L+K)+XC(L+K)
       END DO
       END DO
C       External force (FXTNL = W`, not the 'usual' force!)
        call EXTERNAL CHAINFORCE(FXTNL,FCM,CMX,CHTIME,UG)
C        UG=0.0
C       Non-chained contribution
       UNC=0.0
        DO I=1,3*(N-1)
       FNC(I)=FXTNL(I)
       END DO
       DO I=1,N-2
       LI=3*(I-1)
       DO J=I+2,N
       LJ=3*(J-1)
       RIJ2=0.0 ! +eps*eps
       DO K=1,3
       XAUX(K)=XI(LI+K)-XI(LJ+K)
        RIJ2=RIJ2+XAUX(K)**2
       END DO
       RIJ2INV=1./RIJ2
C       Store the inverse distances.
        LRI=LRI+1
       RINV(LRI)=SQRT(RIJ2INV)
       FM=MIJ(I,J)*RINV(LRI)
        UNC=UNC+FM
       FM=FM*RIJ2INV
C       Fij atraction
       DO K=1,3
       FAUX(K)=FM*XAUX(K)
       END DO
C       Add the contribution to interactions depending on Rij
       DO IK=I,J-1
       L=3*(IK-1)
       DO K=1,3
       FNC(L+K)=FNC(L+K)+FAUX(K)
       END DO
       END DO
       END DO
       END DO
C       Evaluate UQ & TP
        DO I=1,N-1
       L1=3*(I-1)+1
       KS1=4*(I-1)+1
       KS=4*(I-1)
        call QFORCE(Q(KS1),FNC(L1),UQ(KS1))
        call QTxP(Q(KS1),AK(KS1),TP(KS1))
C       The * operation (85)
         AK(KS+4)=-AK(KS+4)
       call QTxP(P(KS1),AK(KS1),TQ(KS1))
       DO K=1,4
       UQ(KS+K)=UQ(KS+K)-2.0D0*MKK(I)*Q(KS+K)*RINV(I)**2  ! RINV**2-> Q^2 Rinv^3
        TQ(KS+K)=TQ(KS+K)-4.D0*DK(I)*Q(KS+K)
       END DO! K
       END DO! I
C       NOTE:Above the division by R (in TP & TQ) is delayed.
C       Proceed to final evaluation of derivatives (90)->(94)
       UPOT=UC+UNC+UG
        G=1./(TKIN+UPOT)
       H=TKIN-UPOT
**********************
c        IF(IWK.EQ.0)THEN
c        WRITE(6,987)' H,E,H-E,UG',H,ENERGY,H-ENERGY,UG
c987       FORMAT(A11,1P4G14.7)
c        IWK=1
c        END IF
**********************
       GAMMA=(H-ENERGY)*G
       GT= (1.-GAMMA)*G
       GU=-(1.+GAMMA)*G
       DO I=1,N-1
       KS=4*(I-1)
C       Apply here the division by R (to TP & TQ)
C       NOTE: TP & TQ never get 'rigth' values. In fact TP=R*TPtrue ...
        GToverR=GT*RINV(I)
       DO K=1,4
       DQ(KS+K)=GToverR*TP(KS+K)
       DP(KS+K)=-GToverR*TQ(KS+K)-GU*UQ(KS+K)
       END DO
       END DO
       DT=G
       DO K=1,3
       DX(K)=CMV(K)*G
        DV(K)=FCM(K)*G
       END DO
C       Evaluate E`
        DE=0.0
        IF(UG.NE.0.0)RETURN
        DO I=1,N-1
       L1=3*(I-1)+1
       KS1=4*(I-1)+1
       KS=4*(I-1)
        call QFORCE(Q(KS1),FXTNL(L1),FAUX)
       DO K=1,4
        DE=DE+DQ(KS+K)*FAUX(K)
       END DO! K
       END DO! I
        RETURN
        END
        SUBROUTINE EXTERNAL CHAINFORCE (FW,FV,XCM,CHTIME,UG)
        INCLUDE 'COMMON1.CH'
*************
*        IMPLICIT REAL*8 (A-H,M,O-Z)
*        PARAMETER (NMX=25,NMX2=2*NMX,NMX3=3*NMX,NMX4=4*NMX,
*     &  NMX8=8*NMX,NMXm=NMX*(NMX-1)/2)
*         COMMON/DataForChainRoutinesOne/X(NMX3),V(NMX3),M(NMX),
*     &   XC(NMX3),WC(NMX3),MC(NMX)
*     &  ,XI(NMX3),PI(NMX3),MASS,RINV(NMXm),RSUM,MKK(NMX)
*     &  ,MIJ(NMX,NMX),TKK(NMX),TK1(NMX),INAME(NMX),N
************
        REAL*8 FW(*),FV(3),XCM(3)
        REAL*8 ACC(NMX3),DDP(NMX3),Q0(3)
        SAVE
       DO K=1,3
        Q0(K)=0.0
       END DO
        DO I=1,N
        L=3*(I-1)
       DO K=1,3
       Q0(K)=Q0(K)+XI(L+K)*MC(I)/MASS
       END DO
       END DO
C       Rearrange according to INAME(i) and add CM.
       DO I=1,N
       L=3*(I-1)
       LF=3*(INAME(I)-1)
       DO K=1,3
       X(LF+K)=XI(L+K)-Q0(K)+XCM(K)
       END DO
       END DO
        call EXTERNALU(X,M,N,UG)
        call EXTERNAL ACCELERATIONS(ACC,XCM,CHTIME)
C       Center of mass force
        DO K=1,3
       FV(K)=0.0
       END DO
        DO I=1,N
        L=3*(I-1)
        DO K=1,3
        FV(K)=FV(K)+M(I)/MASS*ACC(L+K)
        END DO
       END DO
C       Physical chain forces
       DO I=1,N
       L=3*(I-1)
        LF=3*(INAME(I)-1)
        DO K=1,3
        DDP(L+K)=MC(I)*(ACC(LF+K)-FV(K))
        END DO
       END DO
       L=3*(N-2)
        DO K=1,3
        FW(K)=-DDP(K)
       FW(L+K)=DDP(L+K+3)
       END DO
       DO I=2,N-2
       L=3*(I-1)
       DO K=1,3
       FW(L+K)=FW(L+K-3)-DDP(L+K)
       END DO
       END DO
        RETURN
       END
       SUBROUTINE EXTERNALU(X,M,NB,UG)
        IMPLICIT REAL*8 (A-H,M,O-Z)
       REAL*8 X(*),M(*)
       SAVE
C      EXTERNAL potential (if needed). The force must be given in 
c      The subroutine EXTERNAL ACCELERATIONS. 
       UG=0.0
               IF(1.eq.1)RETURN! remove this if U need external potential 
C       THIS IS JUST AN EXAMPLE (replace by what U need)
        L=0
       DO I=1,NB
        DO K=1,3
       L=L+1
       UG=UG-5.0*X(L)**2*M(I)
       END DO
       END DO
       RETURN
       END

       SUBROUTINE EXTERNAL ACCELERATIONS(ACC,XCM,CHTIME)
        INCLUDE 'COMMON1.CH'
C       INCLUDE 'COMMON2.CH'
       REAL*8 ACC(*),XCM(3)
       SAVE
       if(1.eq.1)return  ! remove this if U need external acceration
C         A TEST (modify accordingly also EXTERNALU) 
        do i=1,n
        l=3*(i-1)
        do k=1,3  ! For now set to 0
        acc(l+k)=-10.0*x(l+k) ! real thing must be provided by the user, also in EXTERNALU
        end do !                ! if external forces are wanted
       end do
       RETURN
       END
C************* AUXILIARY ROUTINES FOLLOW *******************
       REAL*8 FUNCTION SQUARE(X,Y)
       REAL*8 X(3),Y(3)
       SAVE
       SQUARE=(X(1)-Y(1))**2+(X(2)-Y(2))**2+(X(3)-Y(3))**2
       RETURN
       END
      SUBROUTINE QTxP (Q,P,W)
C      W = L(Q)P ; KS-matrix of Q times P
      IMPLICIT REAL*8 (A-H,O-Z)
       REAL*8 Q(4),P(4),w(4)
       SAVE
      W(1)=(+Q(1)*P(1)+Q(2)*P(2)+Q(3)*P(3)+Q(4)*P(4))
      W(2)=(-Q(2)*P(1)+Q(1)*P(2)+Q(4)*P(3)-Q(3)*P(4))
      W(3)=(-Q(3)*P(1)-Q(4)*P(2)+Q(1)*P(3)+Q(2)*P(4))
      W(4)=(+Q(4)*P(1)-Q(3)*P(2)+Q(2)*P(3)-Q(1)*P(4))
      RETURN
      END
      SUBROUTINE QxP (Q,P,W)
c      W=L'(Q)P ; L'=transpose of L times P
      IMPLICIT REAL*8 (A-H,O-Z)
       REAL*8 Q(4),P(4),W(4)
       SAVE
      W(1)=(Q(1)*P(1)-Q(2)*P(2)-Q(3)*P(3)+Q(4)*P(4))
      W(2)=(Q(2)*P(1)+Q(1)*P(2)-Q(4)*P(3)-Q(3)*P(4))
      W(3)=(Q(3)*P(1)+Q(4)*P(2)+Q(1)*P(3)+Q(2)*P(4))
      W(4)=(Q(4)*P(1)-Q(3)*P(2)+Q(2)*P(3)-Q(1)*P(4))
      RETURN
      END
        SUBROUTINE CONSTANTS OF MOTION(X,XD,M,NB,ENERGY,G,Alag)
        IMPLICIT real*8 (A-H,m,O-Z)
        DIMENSION X(*),XD(*),G(3),M(*)
        SAVE
        T=0.0
        UG=0.0
        call EXTERNALU(X,M,NB,UG) 
        U=UG
        G(1)=0.
        G(2)=0.
        G(3)=0.
        RMIN=1.D30
        DO 10 I=1,NB
        K1=(I-1)*3+1
        K2=K1+1
        K3=K2+1
        T=T+.5d0*M(I)*(XD(K1)**2+XD(K2)**2+XD(K3)**2)
        G(1)=G(1)+M(I)*(X(K2)*XD(K3)-X(K3)*XD(K2))
        G(2)=G(2)-M(I)*(X(K1)*XD(K3)-X(K3)*XD(K1))
        G(3)=G(3)+M(I)*(X(K1)*XD(K2)-X(K2)*XD(K1))
        IF(I.EQ.Nb)GO TO 10
        J1=I+1
        DO 9 J=J1,Nb
        KI=(I-1)*3
        KJ=(J-1)*3
        R2=0. !+eps*eps
        DO 8 K=1,3
        KI=KI+1
        KJ=KJ+1
8       R2=R2+(X(KI)-X(KJ))**2
        U=U+M(I)*M(J)/SQRT(R2)
9       CONTINUE
10      CONTINUE
        ENERGY=T-U
        Alag=T+U
        RETURN
        END
       SUBROUTINE PH TO KS (XR,PR,Q,P)
        IMPLICIT REAL*8 (A-H,O-Z)
        REAL*8 XR(*),PR(*),Q(*),P(*)
        SAVE
      R2=XR(1)**2+XR(2)**2+XR(3)**2
      R=SQRT(R2)
      Q(4)=0.
      A=R+ABS(XR(1))
      Q(1)=SQRT(.5D0*A)
      B=1./(2.0D0*Q(1))
      Q(2)=XR(2)*B
      Q(3)=XR(3)*B
        IF(XR(1).lt.0.)THEN
        U1=Q(1)
        Q(1)=Q(2)
        Q(2)=U1
        U3=Q(3)
        Q(3)=Q(4)
        Q(4)=U3
       END IF
      P(1)=2.D0*(+Q(1)*PR(1)+Q(2)*PR(2)+Q(3)*PR(3))
      P(2)=2.D0*(-Q(2)*PR(1)+Q(1)*PR(2)+Q(4)*PR(3))
      P(3)=2.D0*(-Q(3)*PR(1)-Q(4)*PR(2)+Q(1)*PR(3))
      P(4)=2.D0*(+Q(4)*PR(1)-Q(3)*PR(2)+Q(2)*PR(3))
      RETURN
      END
      SUBROUTINE KS TO PH (Q,P,XR,PR)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 Q(*),P(*),XR(*),PR(*)
      SAVE
      XR(1)=Q(1)**2-Q(2)**2-Q(3)**2+Q(4)**2
      XR(2)=2.D0*(Q(1)*Q(2)-Q(3)*Q(4))
      XR(3)=2.D0*(Q(1)*Q(3)+Q(2)*Q(4))
      R=Q(1)**2+Q(2)**2+Q(3)**2+Q(4)**2
      A=0.5D0/R
      PR(1)=(Q(1)*P(1)-Q(2)*P(2)-Q(3)*P(3)+Q(4)*P(4))*A
      PR(2)=(Q(2)*P(1)+Q(1)*P(2)-Q(4)*P(3)-Q(3)*P(4))*A
      PR(3)=(Q(3)*P(1)+Q(4)*P(2)+Q(1)*P(3)+Q(2)*P(4))*A
      RETURN
      END
        SUBROUTINE QFORCE(Q,F,qf)
        IMPLICIT REAL*8 (a-h,o-z)
        REAL*8 Q(4),F(3),qf(4)
        SAVE
c       KS transpose times 2F
        qf(1)=2.D0*(+Q(1)*F(1)+Q(2)*F(2)+Q(3)*F(3))
        qf(2)=2.D0*(-Q(2)*F(1)+Q(1)*F(2)+Q(4)*F(3))
        qf(3)=2.D0*(-Q(3)*F(1)-Q(4)*F(2)+Q(1)*F(3))
        qf(4)=2.D0*(+Q(4)*F(1)-Q(3)*F(2)+Q(2)*F(3))
        RETURN
        END
      SUBROUTINE QxQ(Q,XR)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 Q(*),XR(*)
      SAVE
      XR(1)=Q(1)**2-Q(2)**2-Q(3)**2+Q(4)**2
      XR(2)=2.D0*(Q(1)*Q(2)-Q(3)*Q(4))
      XR(3)=2.D0*(Q(1)*Q(3)+Q(2)*Q(4))
       RETURN
       END
      SUBROUTINE DIFSY1(N,F,EPS,H,X,Y)
      IMPLICIT REAL*8 (A-H,O-Z)
C          BULIRSCH-STOER INTEGRATOR.
C          --------------------------
C        Works if Gamma=(H-E)/L, for other time transformations 'eps'
C        must be scaled appropriately such that the test is esentially
C        of the form   (H-E)/L<eps.
C        Convergence test used here is: Abs(Q'dP)<eps  & Abs(P'dQ)<eps
c        see: Mikkola (1987) In 'The Few Body Problem'  p.311 (bottom)
C        NMX=Maximum number of equations; modify if necessary
c    This works only if eqs. are canonical and we have
c    first the P:s then Q:s (or the other way).
c    One additional eq. is allowed (e.g. for t'=...??) but not checked.
      PARAMETER(NMX=250)
      REAL*8 Y(N),YA(NMX),YL(NMX),YM(NMX),DY(NMX),DZ(NMX),DT(NMX,7),
     1D(7),X,XN,H,G,B,B1,U,V,C,TA,W,DABS
      DIMENSION EP(4)
      LOGICAL KONV,BO,KL,GR,FYBAD
      COMMON/INTEGR/JMAX
      DATA EP/0.4D-1,0.16D-2,0.64D-4,0.256D-5/
      save
      NHALF2=(N/2)*2
      JTI=0
      FY=1.
       DO I=1,N
       YA(I)=Y(I)
       END DO
      call F(Y,DZ)
   10 XN=X+H
      BO=.FALSE.
      M=1
      JR=2
      JS=3
      DO  J=1,MAX(4,JMAX)
       IF(BO)THEN
       D(2)=16D0/9.d0
       D(4)=64.d0/9.d0
       D(6)=256.D0/9.d0
       ELSE
       D(2)=2.25D0
       D(4)=9.D0
       D(6)=3.6D1
       END IF
       IF(J.GT.7)THEN
       L=7
       D(7)=6.4D1
       ELSE
       L=J
       D(L)=M*M
       END IF
      KONV=L.GT.3
      M=M+M
      G=H/(M)
      B=G+G
      M=M-1
       DO I=1,N
       YL(I)=YA(I)
       YM(I)=YA(I)+G*DZ(I)
       END DO
      DO K=1,M
      call F(YM(1),DY(1))
       DO I=1,N
       U=YL(I)+B*DY(I)
       YL(I)=YM(I)
       YM(I)=U
       END DO
      END DO
      call F(YM(1),DY(1))
      KL=L.LT.2
      GR=L.GT.5
      FS=0.
      DO I=1,N
      V=DT(I,1)
      C=(YM(I)+YL(I)+G*DY(I))*0.5D0
      DT(I,1)=C
      TA=C
      IF(.NOT.KL)THEN
      DO K=2,L
      B1=D(K)*V
      B=B1-C
      W=C-V
      U=V
       IF(B.NE.0.D0)THEN
       B=W/B
       U=C*B
       C=B1*B
       END IF
      V=DT(I,K)
      DT(I,K)=U
      TA=U+TA
      END DO
      is=I+N/2
      IF(IS.GT.NHALF2)IS=I-(N/2)
      DYIS=DABS(DY(IS))
      IF( I.GT.NHALF2)DYIS=0.0
        IF(KONV)THEN
        TEST=DABS( (Y(I)-TA)*DYIS )
        IF(TEST.GT.EPS) KONV=.FALSE.
        END IF
      IF(.NOT.GR)THEN
      FV=DABS(W)*DYIS
      IF(FS.LT.FV) FS=FV
      END IF
      END IF
      Y(I)=TA
      END DO
       IF(FS.NE.0.D0)THEN
       FA=FY
       K=L-1
       FY=(EP(K)/FS)**(1./(L+K))
       FA7=0.7*FA
       IF(L.EQ.2)FA7=0.0
       FYBAD=.NOT.((FA7.GT.FY).OR.(FY.GT.0.7))
        IF(FY BAD)THEN
        H=H*FY
        JTI=JTI+1
         IF(JTI.GT.5)THEN
         H=0.0
          DO I=1,N
          Y(I)=YA(I)
          END DO
         RETURN
         END IF
        GOTO 10
        END IF
       END IF
      IF(KONV)THEN
      X=XN
      H=H*FY
      RETURN
      END IF
      D(3)=4.D0
      D(5)=1.6D1
      BO=.NOT.BO
      M=JR
      JR=JS
      JS=M+M
      END DO
      H=0.5*H
      GOTO 10
      END
