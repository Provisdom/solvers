//
// Translated by CS2J (http://www.cs2j.com): 2/18/2016 3:06:47 PM
//

package java.Provisdom.Optimization;

import java.Provisdom.LinearAlgebra.BLAS;

public class DoNlp2   
{
    public DoNlp2(int mTotalConstraints, int mEqualityConstraints, int nVariables) throws Exception {
        if (nVariables <= 0)
        {
            Object[] objArray1 = new Object[]{ nVariables };
            throw new IllegalArgumentException("nVariables <= 0");
        }
        else if (mTotalConstraints < 0)
        {
            Object[] objArray2 = new Object[]{ mTotalConstraints };
            throw new IllegalArgumentException("mTotalConstraints < 0");
        }
        else if (mEqualityConstraints < 0)
        {
            Object[] objArray3 = new Object[]{ mEqualityConstraints };
            throw new IllegalArgumentException("mEqualityConstraints < 0");
        }
           
        double[] numArray1 = new double[nVariables];
        this.XGUESS = numArray1;
        double[] numArray2 = new double[nVariables];
        for (int num1 = 0;num1 < nVariables;num1++)
        {
            numArray2[num1] = 1;
        }
        this._scale = numArray2;
        double[] numArray3 = new double[nVariables];
        double[] numArray4 = new double[nVariables];
        this._totalConstraints = mTotalConstraints;
        this._numberOfEqualityConstraints = mEqualityConstraints;
        this._numberOfVariables = nVariables;
        this._lowerBounds = numArray3;
        this._upperBounds = numArray4;
        BLAS.set(nVariables,-1.7976931348623157E+308,this._lowerBounds,0,1);
        BLAS.set(nVariables,1.7976931348623157E+308,this._upperBounds,0,1);
        this._maxIterations = 200;
        this.IPRINT = 0;
        this.PENALTY = 1;
        this.DELTA0 = 0.5;
        this.RELPRE = 2.2204460492503131E-16;
        this.FEPS = 2.2204460492503131E-16;
        this.IDTYPE = 1;
        this.TBND = 1;
        this.SMALL_W = Math.exp(2 * Math.log(7.4014868308343765E-17));
        this.DEL_MIN = Math.min(this.DELTA0 / 10, Math.max(1E-06 * this.DELTA0, this.SMALL_W));
        this.INITIAL_DEL_MIN = this.DEL_MIN;
        this.SCF_MAX = 10000;
        this.NSTEP = 40;
    }

    private void dN10PF() throws Exception {
        for (int num1 = 1;num1 <= this.N;num1++)
        {
            for (int num2 = 1;num2 <= this.N;num2++)
            {
                this.A[num1][num2] = 0;
            }
            this.A[num1][num1] = this.MATSC;
            this.DIAG0[num1] = this.MATSC;
        }
        this.ACCINF[0x1b][this.ITSTEP] = -1;
        this.ACCINF[14][this.ITSTEP] = 1;
        if (!this.SILENT)
        {
            this.dN44PF(20);
        }
         
    }

    private void dN11PF() throws Exception {
        this.DIRDER = this.dN18PF(1,this.N,this.GRADF,this.D) * this.SCF;
        for (int num1 = 1;num1 <= this.NRES;num1++)
        {
            double num2 = this.dN20PF(1,this.N,num1,this.GRES,this.NX,this.D) * this.W[num1];
            double num3 = this.RES[num1];
            if (num1 <= this.NH)
            {
                if ((num3 / this.GRESN[num1]) < (-1000 * this.EPSMAC))
                {
                    this.DIRDER -= num2;
                }
                else if ((num3 / this.GRESN[num1]) > (1000 * this.EPSMAC))
                {
                    this.DIRDER += num2;
                }
                else
                {
                    this.DIRDER += Math.abs(num2);
                }  
            }
            else if (this.BIND[num1] == 1)
            {
                if ((Math.abs(num3) / this.GRESN[num1]) <= (1000 * this.EPSMAC))
                {
                    this.DIRDER -= Math.min(0, num2);
                }
                else if ((num3 / this.GRESN[num1]) < (-1000 * this.EPSMAC))
                {
                    if (num2 > 0)
                    {
                        num2 = Math.min(num2, -this.RES[num1] * this.W[num1]);
                    }
                     
                    this.DIRDER -= num2;
                }
                  
            }
              
        }
    }

    private void dN12PF() throws Exception {
        int num1;
        this.XNORM = this.dN31PF(0,this.N - 1,this.X);
        double num2 = this.BETA * (this.XNORM + 1);
        this.DNORM = this.dN31PF(1,this.N,this.D);
        this.e0NORM = this.dN31PF(1,this.N,this.E0);
        this.DSCAL = 1;
        if ((this.DNORM * this.e0NORM) == 0)
        {
            this.COSPHI = 0;
        }
        else
        {
            this.COSPHI = this.dN18PF(1,this.N,this.D,this.E0) / (this.e0NORM * this.DNORM);
        } 
        if (this.DNORM > num2)
        {
            double num3 = num2 / this.DNORM;
            this.DNORM = num2;
            this.DSCAL = num3;
            for (num1 = 1;num1 <= this.N;num1++)
            {
                this.D[num1] *= num3;
                this.DD[num1] = (this.DD[num1] * num3) * num3;
            }
        }
         
        for (num1 = 1;num1 <= this.N;num1++)
        {
            if (this.LLOW[num1] && ((this.X[num1 - 1] + (this.SIGSM * this.D[num1])) <= this.UG[num1]))
            {
                this.D[num1] = 0;
                this.DD[num1] = Math.max(0, this.DD[num1]);
            }
             
            if (this.LUP[num1] && ((this.X[num1 - 1] + (this.SIGSM * this.D[num1])) >= this.OG[num1]))
            {
                this.D[num1] = 0;
                this.DD[num1] = Math.min(0, this.DD[num1]);
            }
             
        }
        this.DNORM = this.dN31PF(1,this.N,this.D);
    }

    private void dN13PF() throws Exception {
        boolean flag1 = true;
        int num1 = 1;
        while (num1 <= this.N)
        {
            flag1 = flag1 && (((this.D[num1] == 0) || (this.LUP[num1] && (this.D[num1] > 0))) || (this.LLOW[num1] && (this.D[num1] < 0)));
            num1++;
        }
        if (flag1)
        {
            this.STMAXL = this.SIGSM;
            for (num1 = 1;num1 <= this.N;num1++)
            {
                if (this.LLOW[num1] && (this.D[num1] < 0))
                {
                    if ((-this.D[num1] * this.SIGLA) >= (this.X[num1 - 1] - this.UG[num1]))
                    {
                        this.STMAXL = Math.max(this.STMAXL, (this.X[num1 - 1] - this.UG[num1]) / -this.D[num1]);
                    }
                    else
                    {
                        this.STMAXL = this.SIGLA;
                    } 
                }
                 
                if (this.LUP[num1] && (this.D[num1] > 0))
                {
                    if ((this.D[num1] * this.SIGLA) >= (this.OG[num1] - this.X[num1 - 1]))
                    {
                        this.STMAXL = Math.max(this.STMAXL, (this.OG[num1] - this.X[num1 - 1]) / this.D[num1]);
                    }
                    else
                    {
                        this.STMAXL = this.SIGLA;
                    } 
                }
                 
            }
        }
        else
        {
            this.STMAXL = this.SIGLA;
        } 
        this.STMAXL = Math.min(this.SIGLA, this.STMAXL);
    }

    private void dN14PF() throws Exception {
        this.PHI1 = this.PHIMIN;
        this.PSI1 = this.PSIMIN;
        this.UPSI1 = this.UPSIM;
        this.SIG = this.SIGMIN;
        this.FX1 = this.FMIN;
        int num1 = 0;
        while (num1 < this.N)
        {
            this.X1[num1] = this.XMIN[num1];
            num1++;
        }
        for (num1 = 1;num1 <= this.NRES;num1++)
        {
            this.RES1[num1] = this.RESMIN[num1];
        }
    }

    private void dN15PF() throws Exception {
        this.PHIMIN = this.PHI1;
        this.UPSIM = this.UPSI1;
        this.PSIMIN = this.PSI1;
        this.FMIN = this.FX1;
        this.SIGMIN = this.SIG;
        int num1 = 0;
        while (num1 < this.N)
        {
            this.XMIN[num1] = this.X1[num1];
            num1++;
        }
        for (num1 = 1;num1 <= this.NRES;num1++)
        {
            this.RESMIN[num1] = this.RES1[num1];
        }
    }

    private void dN16PF(double SIGACT, double[] SIGRES, boolean[] REJECT, boolean[] ERROR, java.Provisdom.Optimization.DoNlp2.IFunction FCN, double[] XLB, double[] XUB) throws Exception {
        this.SIG = SIGACT;
        int num1 = 1;
        while (num1 <= this.N)
        {
            this.X1[num1 - 1] = this.X[num1 - 1] + (this.SIG * (this.D[num1] + (this.SIG * this.DD[num1])));
            if (this.LLOW[num1])
            {
                this.X1[num1 - 1] = Math.max(this.X1[num1 - 1], this.UG[num1]);
            }
             
            if (this.LUP[num1])
            {
                this.X1[num1 - 1] = Math.min(this.X1[num1 - 1], this.OG[num1]);
            }
             
            num1++;
        }
        REJECT[0] = false;
        ERROR[0] = false;
        SIGRES[0] = this.SIG;
        this.UPSI1 = 0;
        this.PSI1 = 0;
        for (int num2 = 1;num2 <= this.NRES;num2++)
        {
            double num3;
            num1 = this.SORT[num2];
            if (num1 <= this.NH)
            {
                this.CFUERR[num1] = false;
                this.RES1[num1] = this.dN40PF(num1, this.X1, FCN, XLB, XUB);
                if (this.CFUERR[num1])
                {
                    ERROR[0] = true;
                    return ;
                }
                 
                num3 = Math.abs(this.RES1[num1]);
            }
            else
            {
                this.CFUERR[num1] = false;
                this.RES1[num1] = this.dN42PF(num1 - this.NH, this.X1, FCN, XLB, XUB);
                if (this.CFUERR[num1])
                {
                    ERROR[0] = true;
                    return ;
                }
                 
                num3 = -Math.min(0, this.RES1[num1]);
                if ((this.RES1[num1] < -this.DELMIN) && (this.BIND[num1] == 0))
                {
                    this.VIOLIS[0]++;
                    this.VIOLIS[this.VIOLIS[0]] = num1;
                }
                 
            } 
            this.UPSI1 += num3;
            if ((this.UPSI1 > this.TAU0) && (this.PHASE != -1))
            {
                REJECT[0] = true;
                return ;
            }
             
            this.PSI1 += num3 * this.W[num1];
            if ((((this.RES1[num1] * this.RES[num1]) < 0) && (this.SIG <= 1)) && ((this.BIND[num1] == 0) || ((this.BIND[num1] == 1) && (((Math.abs(this.RES[num1]) / this.GRESN[num1]) >= (1000 * this.EPSMAC)) || ((Math.abs(this.RES1[num1]) / this.GRESN[num1]) >= (1000 * this.EPSMAC))))))
            {
                SIGRES[0] = Math.min(SIGRES[0], (this.SIG * this.RES[num1]) / (this.RES[num1] - this.RES1[num1]));
            }
             
        }
        if (this.PHASE == -1)
        {
            this.FX1 = 0;
        }
        else
        {
            this.FFUERR = false;
            this.FX1 = this.dN38PF(this.X1, FCN, XLB, XUB);
            if (this.FFUERR)
            {
                ERROR[0] = true;
                return ;
            }
             
        } 
        this.PHI1 = (this.SCF * this.FX1) + this.PSI1;
    }

    private void dN17PF(double SIG1TH, java.Provisdom.Optimization.DoNlp2.IFunction FCN, double[] XLB, double[] XUB) throws Exception {
        double[] numArray1 = new double[this.NSTEP + 1];
        double[] numArray2 = new double[1];
        boolean[] flagArray1 = new boolean[1];
        boolean[] flagArray2 = new boolean[1];
        double num6 = 0;
        boolean flag5 = false;
        numArray1[1] = 0.5;
        numArray1[2] = 0.25;
        int num1 = 3;
        while (num1 <= 40) {
            numArray1[num1] = 0.1;
            num1++;
        }
        int num2 = 0;
        boolean flag6 = false;
        this.PHI = (this.SCF * this.FX) + this.PSI;
        this.SIG = SIG1TH;
        this.VIOLIS[0] = 0;
        if (!this.SILENT) {
            this.dN44PF(8);
        }

        Label_00AC:
        {
            num2++;
            if (num2 > this.NSTEP) {
                this.STPTRM = -1;
                this.SIG = 0;
            } else {
                numArray2[0] = num6;
                flagArray1[0] = flag5;
                flagArray2[0] = flag6;
                this.dN16PF(this.SIG, numArray2, flagArray1, flagArray2, FCN, XLB, XUB);
                num6 = numArray2[0];
                flag5 = flagArray1[0];
                if (flagArray2[0]) {
                    if (this.SIG > 1) {
                        this.dN14PF();
                        this.FX0 = this.FX;
                        this.FX = this.FX1;
                        this.UPSI0 = this.UPSI;
                        this.UPSI = this.UPSI1;
                        this.PSI0 = this.PSI;
                        this.PSI = this.PSI1;
                        this.STPTRM = 1;
                        this.SIG0 = this.SIG;
                        num1 = 1;
                        while (num1 <= this.N) {
                            this.X0[num1 - 1] = this.X[num1 - 1];
                            this.E0[num1] = this.D[num1];
                            this.X[num1 - 1] = this.X1[num1 - 1];
                            this.DifX[num1] = this.X[num1 - 1] - this.X0[num1 - 1];
                            num1++;
                        }
                        this.e0NORM = this.DNORM;
                        this.X0NORM = this.XNORM;
                        for (num1 = 1; num1 <= this.NRES; num1++) {
                            this.RES[num1] = this.RES1[num1];
                        }
                        return;
                    }

                    this.SIG = numArray1[num2] * this.SIG;
                    break Label_00AC;
                }

                if (!flag5) {
                    if (!this.SILENT) {
                        this.dN44PF(9);
                    }

                    if (this.SIG <= 1) {
                        double num7;
                        if (((this.LASTCH >= (this.ITSTEP - 3)) || (this.PHASE != 2)) || this.SINGUL) {
                            num7 = this.PHI - this.PHI1;
                        } else {
                            double num5 = this.PHI;
                            for (int num3 = 1; num3 <= 3; num3++) {
                                num5 = Math.max((this.SCF * this.ACCINF[2, this.ITSTEP - num3])
                                +this.ACCINF[4, this.ITSTEP - num3],num5);
                            }
                            num7 = num5 - this.PHI1;
                        }
                        boolean flag1 = num7 >= Math.min((-this.SIG * this.DELTA) * this.DIRDER, this.LEVEL);
                        boolean flag2 = (this.UPSI - this.UPSI1) >= ((((this.SIG * this.DELTA) * this.DELTA) * this.UPSI) / this.TAUQP);
                        boolean flag3 = (this.UPSI <= (this.TAU0 * 0.5)) && (this.UPSI1 <= this.TAU0);
                        boolean flag4 = this.UPSI > (this.TAU0 * 0.5);
                        if (flag1 && (flag3 || (flag4 && flag2))) {
                            if (((this.SIG == 1) && ((((this.COSPHI >= this.THETA) && (this.SIG0 >= 1)) && ((((this.PHASE + 1) * (this.PHASE - 2)) != 0) && !this.SINGUL)) || (num7 >= ((-this.SIG * this.DELTA1) * this.DIRDER)))) && ((this.STMAXL > 1) && (this.UPSI < (this.TAU0 * 0.5)))) {
                                this.dN15PF();
                                this.SIG = Math.min(this.STMAXL, this.SIG + this.SIG);
                                break Label_00AC;
                            }

                            if (((this.SIG > 1) || (this.UPSI <= (this.TAU0 * 0.5))) || (this.UPSI1 <= this.UPSI)) {
                                this.FX0 = this.FX;
                                this.FX = this.FX1;
                                this.UPSI0 = this.UPSI;
                                this.UPSI = this.UPSI1;
                                this.PSI0 = this.PSI;
                                this.PSI = this.PSI1;
                                this.STPTRM = 1;
                                this.SIG0 = this.SIG;
                                num1 = 1;
                                while (num1 <= this.N) {
                                    this.X0[num1 - 1] = this.X[num1 - 1];
                                    this.E0[num1] = this.D[num1];
                                    this.X[num1 - 1] = this.X1[num1 - 1];
                                    this.DifX[num1] = this.X[num1 - 1] - this.X0[num1 - 1];
                                    num1++;
                                }
                                this.e0NORM = this.DNORM;
                                this.X0NORM = this.XNORM;
                                for (num1 = 1; num1 <= this.NRES; num1++) {
                                    this.RES[num1] = this.RES1[num1];
                                }
                                return;
                            }

                        }

                        if (num6 < this.SIG) {
                            this.SIG = Math.min(0.5 * this.SIG, Math.max(numArray1[num2] * this.SIG, num6));
                        } else {
                            double num4 = (num7 - (this.DIRDER * this.SIG)) * 2;
                            if (num4 > (this.EPSMAC * ((this.SCF * Math.abs(this.FX)) + this.PSI))) {
                                this.SIG = Math.min(0.5 * this.SIG, Math.max((double) (numArray1[num2] * this.SIG), (double) (((-this.DIRDER * this.SIG) * this.SIG) / num4)));
                            } else {
                                this.SIG = numArray1[num2] * this.SIG;
                            }
                        }
                        if ((this.SIG * Math.max(1, this.DNORM)) < this.SIGSM) {
                            this.STPTRM = -1;
                            this.SIG = 0;
                            return;
                        }

                        break Label_00AC;
                    }

                    if (this.PHI1 >= this.PHIMIN) {
                        this.dN14PF();
                        this.FX0 = this.FX;
                        this.FX = this.FX1;
                        this.UPSI0 = this.UPSI;
                        this.UPSI = this.UPSI1;
                        this.PSI0 = this.PSI;
                        this.PSI = this.PSI1;
                        this.STPTRM = 1;
                        this.SIG0 = this.SIG;
                        num1 = 1;
                        while (num1 <= this.N) {
                            this.X0[num1 - 1] = this.X[num1 - 1];
                            this.E0[num1] = this.D[num1];
                            this.X[num1 - 1] = this.X1[num1 - 1];
                            this.DifX[num1] = this.X[num1 - 1] - this.X0[num1 - 1];
                            num1++;
                        }
                        this.e0NORM = this.DNORM;
                        this.X0NORM = this.XNORM;
                        for (num1 = 1; num1 <= this.NRES; num1++) {
                            this.RES[num1] = this.RES1[num1];
                        }
                    } else {
                        if (this.SIG < this.STMAXL) {
                            this.dN15PF();
                            this.SIG = Math.min(this.STMAXL, this.SIG + this.SIG);
                            break Label_00AC;
                        }

                        this.FX0 = this.FX;
                        this.FX = this.FX1;
                        this.UPSI0 = this.UPSI;
                        this.UPSI = this.UPSI1;
                        this.PSI0 = this.PSI;
                        this.PSI = this.PSI1;
                        this.STPTRM = 1;
                        this.SIG0 = this.SIG;
                        num1 = 1;
                        while (num1 <= this.N) {
                            this.X0[num1 - 1] = this.X[num1 - 1];
                            this.E0[num1] = this.D[num1];
                            this.X[num1 - 1] = this.X1[num1 - 1];
                            this.DifX[num1] = this.X[num1 - 1] - this.X0[num1 - 1];
                            num1++;
                        }
                        this.e0NORM = this.DNORM;
                        this.X0NORM = this.XNORM;
                        for (num1 = 1; num1 <= this.NRES; num1++) {
                            this.RES[num1] = this.RES1[num1];
                        }
                    }
                } else if (this.SIG > 1) {
                    this.dN14PF();
                    this.FX0 = this.FX;
                    this.FX = this.FX1;
                    this.UPSI0 = this.UPSI;
                    this.UPSI = this.UPSI1;
                    this.PSI0 = this.PSI;
                    this.PSI = this.PSI1;
                    this.STPTRM = 1;
                    this.SIG0 = this.SIG;
                    num1 = 1;
                    while (num1 <= this.N) {
                        this.X0[num1 - 1] = this.X[num1 - 1];
                        this.E0[num1] = this.D[num1];
                        this.X[num1 - 1] = this.X1[num1 - 1];
                        this.DifX[num1] = this.X[num1 - 1] - this.X0[num1 - 1];
                        num1++;
                    }
                    this.e0NORM = this.DNORM;
                    this.X0NORM = this.XNORM;
                    for (num1 = 1; num1 <= this.NRES; num1++) {
                        this.RES[num1] = this.RES1[num1];
                    }
                } else {
                    this.SIG = numArray1[num2] * this.SIG;
                    break Label_00AC;
                }
            }
        }
    }

    private double dN18PF(int i, int j, double[] A, double[] B) throws Exception {
        if (i > j)
        {
            return 0;
        }
         
        double num3 = 0;
        for (int num1 = i;num1 <= j;num1++)
        {
            num3 += A[num1] * B[num1];
        }
        return num3;
    }

    private double dN19PF(int N, int M, int j, double[][] A, int LDA, double[] B) throws Exception {
        double num2 = 0;
        for (int num3 = N;num3 <= M;num3++)
        {
            num2 += A[num3][j] * B[num3];
        }
        return num2;
    }

    private double dN20PF(int N, int M, int j, double[][] A, int LDA, double[] B) throws Exception {
        double num2 = 0;
        for (int num3 = N;num3 <= M;num3++)
        {
            num2 += A[j][num3] * B[num3];
        }
        return num2;
    }

    private void dN21PG(java.Provisdom.Optimization.DoNlp2.IGradient GRAD, int FTYPE, int I, double[] XTR, double[] VALUES) throws Exception {
        this._totalConstraints = (this.NG + this.NH) + this.NBOUNDS;
        for (int num1 = 0;num1 < this.N;num1++)
        {
            VALUES[num1] = 0;
        }
        if (FTYPE == 1)
        {
            this.ICGF++;
            GRAD.gradient(XTR, 0, VALUES);
        }
        else if (FTYPE == 2)
        {
            this.CGRES[I]++;
            GRAD.gradient(XTR, I, VALUES);
        }
        else
        {
            if (FTYPE != 3)
            {
                return ;
            }
             
            this.CGRES[I + this.NH]++;
            GRAD.gradient(XTR, I + this.NH, VALUES);
        }  
    }

    private void dN22PF(double[] GPHI) throws Exception {
        for (int num1 = 1;num1 <= this.N;num1++)
        {
            GPHI[num1] = this.GRADF[num1] * this.SCF;
            int num2 = 1;
            while (num2 <= this.NH)
            {
                GPHI[num1] -= this.U[num2] * this.GRES[num2][num1];
                num2++;
            }
            for (num2 = this.NH + 1;num2 <= this.ALIST[0];num2++)
            {
                int num3 = this.ALIST[num2];
                if (this.U[num3] > 0)
                {
                    GPHI[num1] -= this.GRES[num3][num1] * this.U[num3];
                }
                 
            }
        }
    }

    private void dN23PF(int NLOW, int NRL) throws Exception {
        double[] numArray1 = new double[this.NX + 1];
        double[] numArray2 = new double[this.NX + 1];
        double num10 = 0;
        if (NLOW <= NRL)
        {
            int num1;
            int num2;
            int num4;
            int num5;
            int num7;
            int num9;
            double num11;
            if (NLOW == 1)
            {
                this.RANK = 0;
            }
             
            double num15 = 1 / ((double)((this.N + this.N) + this.N));
            int num3 = NLOW;
            while (num3 <= NRL)
            {
                this.DIAG[num3] = 0;
                this.BETAQ[num3] = 0;
                this.COLNO[num3] = num3;
                num4 = 1;
                while (num4 <= this.N)
                {
                    numArray1[num4] = this.GRES[this.ALIST[num3]][num4];
                    num4++;
                }
                num10 = this.dN30PF(this.A, numArray1, numArray2, this.N);
                if (num10 == 0)
                {
                    this.CSCAL[num3] = 1;
                    this.COLLE[num3] = 0;
                    for (num4 = 1;num4 <= this.N;num4++)
                    {
                        this.QR[num3][num4] = 0;
                    }
                }
                else
                {
                    num4 = 1;
                    while (num4 <= this.N)
                    {
                        numArray1[num4] = numArray2[this.PERM[num4]];
                        num4++;
                    }
                    num11 = 1 / Math.sqrt(Math.max(num10, this.RHO * this.RHO));
                    this.CSCAL[num3] = num11;
                    if (NLOW > 1)
                    {
                        this.dN24PF(1, 0, 1, this.RANK, this.N, this.QR, this.BETAQ, numArray1, numArray2);
                        for (num4 = 1;num4 <= this.N;num4++)
                        {
                            numArray1[num4] = numArray2[num4];
                        }
                    }
                     
                    for (num4 = 1;num4 <= this.N;num4++)
                    {
                        this.QR[num3][num4] = numArray1[num4] * num11;
                    }
                    this.COLLE[num3] = (this.dN31PF(this.RANK + 1, this.N, numArray1) * num11) * (this.dN31PF(this.RANK + 1, this.N, numArray1) * num11);
                } 
                num3++;
            }
            if ((NLOW > 1) && (this.RANK < (NLOW - 1)))
            {
                num7 = (NLOW - 1) - this.RANK;
                int num8 = (NRL - NLOW) + 1;
                for (num3 = 1;num3 <= Math.min(num7, num8);num3++)
                {
                    num9 = this.RANK + num3;
                    num5 = (NRL - num3) + 1;
                    num11 = this.BETAQ[num5];
                    this.BETAQ[num5] = this.BETAQ[num9];
                    this.BETAQ[num9] = num11;
                    num4 = this.COLNO[num5];
                    this.COLNO[num5] = this.COLNO[num9];
                    this.COLNO[num9] = num4;
                    num11 = this.COLLE[num5];
                    this.COLLE[num5] = this.COLLE[num9];
                    this.COLLE[num9] = num11;
                    for (num4 = 1;num4 <= this.N;num4++)
                    {
                        num11 = this.QR[num5, num4];
                        this.QR[num5][num4] = this.QR[num9][num4];
                        this.QR[num9][num4] = num11;
                    }
                }
            }
             
            if (NLOW > 1)
            {
                num1 = this.RANK + 1;
                num2 = (num1 + NRL) - NLOW;
            }
            else
            {
                num1 = NLOW;
                num2 = NRL;
            } 
            for (num3 = num1;num3 <= num2;num3++)
            {
                num9 = num3;
                double num16 = this.COLLE[num3];
                num4 = num3 + 1;
                while (num4 <= num2)
                {
                    if (this.COLLE[num4] > num16)
                    {
                        num16 = this.COLLE[num4];
                    }
                     
                    num4++;
                }
                num4 = num2;
                while (num4 >= num3)
                {
                    if (this.COLLE[num4] >= (num16 / 3))
                    {
                        num9 = num4;
                    }
                     
                    num4--;
                }
                if (num9 != num3)
                {
                    num4 = this.COLNO[num3];
                    this.COLNO[num3] = this.COLNO[num9];
                    this.COLNO[num9] = num4;
                    num11 = this.COLLE[num3];
                    this.COLLE[num3] = this.COLLE[num9];
                    this.COLLE[num9] = num11;
                    for (num5 = 1;num5 <= this.N;num5++)
                    {
                        num11 = this.QR[num3][num5];
                        this.QR[num3][num5] = this.QR[num9][num5];
                        this.QR[num9][num5] = num11;
                    }
                }
                 
                num10 = 0;
                num4 = num3;
                while (num4 <= this.N)
                {
                    num11 = this.QR[num3][num4];
                    numArray1[num4] = num11;
                    num10 = (num11 * num11) + num10;
                    num4++;
                }
                if (num10 <= (this.RHO * this.RHO))
                {
                    for (num4 = num3;num4 <= num2;num4++)
                    {
                        this.COLLE[num4] = 0;
                        for (num5 = num3;num5 <= this.N;num5++)
                        {
                            this.QR[num4][num5] = 0;
                        }
                    }
                    this.RANK = num3 - 1;
                    return ;
                }
                 
                double num14 = numArray1[num3];
                double num12 = -Math.sqrt(num10);
                if (Math.abs(num14) <= (-num12 * num15))
                {
                    num11 = 0;
                    int num6 = num3 + 1;
                    for (num4 = num3 + 1;num4 <= this.N;num4++)
                    {
                        if (Math.abs(numArray1[num4]) > num11)
                        {
                            num11 = Math.abs(numArray1[num4]);
                            num6 = num4;
                        }
                         
                    }
                    num5 = this.PERM1[num3];
                    this.PERM1[num3] = this.PERM1[num6];
                    this.PERM1[num6] = num5;
                }
                 
                if (num14 < 0)
                {
                    num12 = -num12;
                }
                 
                double num13 = 1 / (num10 - (num14 * num12));
                this.DIAG[num3] = num12;
                this.BETAQ[num3] = num13;
                numArray1[num3] = num14 - num12;
                this.QR[num3][num3] = numArray1[num3];
                this.RANK = num3;
                num7 = num3 + 1;
                for (num4 = num7;num4 <= num2;num4++)
                {
                    num10 = num13 * this.dN20PF(num3, this.N, num4, this.QR, this.NX, numArray1);
                    for (num5 = num3;num5 <= this.N;num5++)
                    {
                        this.QR[num4][num5] -= num10 * numArray1[num5];
                    }
                    this.COLLE[num4] -= this.QR[num4][num3] * this.QR[num4][num3];
                }
            }
        }
         
    }

    private void dN24PF(int ID, int INCR, int IS1, int IS2, int M, double[][] A, double[] BETA, double[] B, double[] C) throws Exception {
        int num2 = 1;
        while (num2 <= M)
        {
            C[num2] = B[num2];
            num2++;
        }
        if ((IS1 <= M) && (IS2 >= IS1))
        {
            for (num2 = IS1;num2 <= IS2;num2++)
            {
                int num5 = num2;
                if (ID < 0)
                {
                    num5 = (IS2 - num5) + IS1;
                }
                 
                int num3 = num5 + INCR;
                double num1 = BETA[num5] * this.dN20PF(num3, M, num5, A, this.NX, C);
                for (int num4 = num3;num4 <= M;num4++)
                {
                    C[num4] -= num1 * A[num5][num4];
                }
            }
        }
         
    }

    private void dN25PF(int NLOW, int NUP, double[] B, double[] XX) throws Exception {
        double[] numArray1 = new double[this.NX + 1];
        int num2 = NUP;
        while (num2 >= NLOW)
        {
            double num1 = 0;
            for (int num3 = num2 + 1;num3 <= NUP;num3++)
            {
                num1 += this.QR[num3][num2] * numArray1[num3];
            }
            numArray1[num2] = (B[num2] - num1) / this.DIAG[num2];
            num2--;
        }
        for (num2 = NLOW;num2 <= NUP;num2++)
        {
            XX[num2] = numArray1[num2] * this.CSCAL[this.COLNO[num2]];
        }
    }

    private void dN26PF(int NLOW, int NUP, double[] B, double[] XX) throws Exception {
        int num1 = NLOW;
        while (num1 <= NUP)
        {
            XX[num1] = B[num1] * this.CSCAL[this.COLNO[num1]];
            num1++;
        }
        for (num1 = NLOW;num1 <= NUP;num1++)
        {
            double num3 = 0;
            for (int num2 = NLOW;num2 <= (num1 - 1);num2++)
            {
                num3 += this.QR[num1][num2] * XX[num2];
            }
            XX[num1] = (XX[num1] - num3) / this.DIAG[num1];
        }
    }

    private double dN27PF(double A, double B) throws Exception {
        double num1 = Math.abs(A);
        double num2 = Math.abs(B);
        if (num1 > num2)
        {
            return (num1 * Math.sqrt(1 + ((num2 / num1) * (num2 / num1))));
        }
         
        if (num2 > num1)
        {
            return (num2 * Math.sqrt(1 + ((num1 / num2) * (num1 / num2))));
        }
         
        return (num1 * Math.sqrt(2));
    }

    private boolean dN28PF(double[][] R, double[] Z, double[] Y, int N) throws Exception {
        int num2;
        double num8;
        double num9;
        double num10;
        double[] numArray1 = new double[this.NX + 1];
        double[] numArray2 = new double[this.NX + 1];
        double[] numArray3 = new double[this.NX + 1];
        boolean flag1 = false;
        double num6 = 0;
        int num1 = 1;
        while (num1 <= (N - 1))
        {
            numArray1[num1] = R[num1][num1 + 1];
            R[num1][num1 + 1] = 0;
            num1++;
        }
        double num5 = 0;
        num1 = 1;
        while (num1 <= N)
        {
            num5 += Z[num1] * Z[num1];
            num1++;
        }
        if (num5 != 0)
        {
            int num3;
            num6 = this.dN30PF(R, Z, numArray3, N);
            num6 = Math.sqrt(num6 + 1);
            num1 = N;
            while (num1 >= 2)
            {
                if (numArray3[num1] != 0)
                {
                    num3 = num1 - 1;
                    num8 = numArray3[num3];
                    num9 = numArray3[num1];
                    numArray3[num3] = this.dN27PF(num8,num9);
                    num8 /= numArray3[num3];
                    num9 = -num9 / numArray3[num3];
                    R[num3][num1] = num9 * R[num3][num3];
                    R[num3][num3] = num8 * R[num3][num3];
                    for (num2 = num1;num2 <= N;num2++)
                    {
                        num10 = (num8 * R[num2][num3]) - (num9 * R[num2][num1]);
                        R[num2][num1] = (num9 * R[num2][num3]) + (num8 * R[num2][num1]);
                        R[num2][num3] = num10;
                    }
                }
                 
                num1--;
            }
            num1 = 1;
            while (num1 <= N)
            {
                R[num1][1] *= num6;
                num1++;
            }
            for (num1 = 1;num1 <= (N - 1);num1++)
            {
                num3 = num1 + 1;
                num8 = R[num1][num1];
                num9 = -R[num1][num3];
                num10 = this.dN27PF(num8,num9);
                if (num10 != 0)
                {
                    num8 /= num10;
                    num9 /= num10;
                    R[num1][num1] = num10;
                    R[num1][num3] = 0;
                    for (num2 = num1 + 1;num2 <= N;num2++)
                    {
                        num10 = (num8 * R[num2][num1]) - (num9 * R[num2][num3]);
                        R[num2][num3] = (num9 * R[num2][num1]) + (num8 * R[num2][num3]);
                        R[num2][num1] = num10;
                    }
                }
                 
            }
        }
         
        double num4 = 0;
        num1 = 1;
        while (num1 <= N)
        {
            num4 += Y[num1] * Y[num1];
            num1++;
        }
        if (num4 != 0)
        {
            num6 = this.dN30PF(R, Y, numArray3, N);
            if (num6 >= 1)
            {
                flag1 = true;
            }
            else
            {
                num6 = Math.sqrt(Math.abs((double) (1 - num6)));
                double num7 = num6;
                for (num1 = N;num1 >= 1;num1--)
                {
                    num8 = num7;
                    num9 = numArray3[num1];
                    num7 = this.dN27PF(num8,num9);
                    if (num7 != 0)
                    {
                        num8 /= num7;
                        num9 /= num7;
                        numArray2[num1] = num9 * R[num1][num1];
                        R[num1][num1] = num8 * R[num1][num1];
                        for (num2 = num1 + 1;num2 <= N;num2++)
                        {
                            num10 = (num8 * R[num2][num1]) - (num9 * numArray2[num2]);
                            numArray2[num2] = (num9 * R[num2][num1]) + (num8 * numArray2[num2]);
                            R[num2][num1] = num10;
                        }
                    }
                     
                }
            } 
        }
         
        for (num1 = 1;num1 <= (N - 1);num1++)
        {
            R[num1][num1 + 1] = numArray1[num1];
        }
        return flag1;
    }

    private double dN29PF(double[][] A, double[] B, double[] Y, int N) throws Exception {
        double num4 = 0;
        for (int num1 = N;num1 >= 1;num1--)
        {
            double num3 = B[num1];
            for (int num2 = num1 + 1;num2 <= N;num2++)
            {
                num3 -= A[num2][num1] * Y[num2];
            }
            num3 /= A[num1][num1];
            Y[num1] = num3;
            num4 = (num3 * num3) + num4;
        }
        return num4;
    }

    private double dN30PF(double[][] A, double[] B, double[] Y, int N) throws Exception {
        double num4 = 0;
        for (int num1 = 1;num1 <= N;num1++)
        {
            double num3 = B[num1];
            for (int num2 = 1;num2 <= (num1 - 1);num2++)
            {
                num3 -= A[num1][num2] * Y[num2];
            }
            num3 /= A[num1][num1];
            Y[num1] = num3;
            num4 = (num3 * num3) + num4;
        }
        return num4;
    }

    private double dN31PF(int NL, int NM, double[] XX) throws Exception {
        if (NM < NL)
        {
            return 0;
        }
         
        double num3 = Math.abs(XX[NL]);
        int num1 = NL + 1;
        while (num1 <= NM)
        {
            num3 = Math.max(num3, Math.abs(XX[num1]));
            num1++;
        }
        if (num3 == 0)
        {
            return 0;
        }
         
        double num4 = 0;
        for (num1 = NL;num1 <= NM;num1++)
        {
            num4 += (XX[num1] / num3) * (XX[num1] / num3);
        }
        return (num3 * Math.sqrt(num4));
    }

    private void dN32PF() throws Exception {
        int num4;
        int num5;
        int num6;
        int num7;
        double num12;
        double num14;
        double num15;
        double num16;
        double num17;
        double num18;
        double num19;
        double num25;
        double num26;
        double num27;
        boolean flag1;
        double[] numArray1 = new double[(this.NX + this.NRESM) + 1];
        double[] numArray2 = new double[(2 * this.NRESM) + 1];
        double[] numArray3 = new double[(2 * this.NRESM) + 1];
        double[] numArray4 = new double[(this.NX + this.NRESM) + 1];
        double[] numArray5 = new double[(this.NX + this.NRESM) + 1];
        double[] numArray6 = new double[(2 * this.NRESM) + 1];
        double[] numArray7 = new double[(this.NX + this.NRESM) + 1];
        double[] numArray8 = new double[(2 * this.NRESM) + 1];
        double[] numArray9 = new double[(this.NX + this.NRESM) + 1];
        double[] numArray10 = new double[(2 * this.NRESM) + 1];
        double[] numArray11 = new double[(this.NX + this.NRESM) + 1];
        double[] numArray12 = new double[this.NRESM + 1];
        int[] numArray13 = new int[(2 * this.NRESM) + 1];
        int[] numArray14 = new int[(2 * this.NRESM) + 1];
        int[] numArray15 = new int[(2 * this.NRESM) + 1];
        int[] numArray16 = new int[(2 * this.NRESM) + 1];
        int[] numArray17 = new int[(2 * this.NRESM) + 1];
        double num28 = 1;
        int num1 = this.NX + this.NRESM;
        int num2 = this.NRESM * this.NRESM;
        this.NDUAL = this.N + this.NR;
        this._numberOfEqualityConstraints = this.NH;
        this.MI = (2 * this.NR) - this.NH;
        double num29 = this.EPSMAC / this.TOLMAC;
        if (this.ANALYT)
        {
            num12 = ((2 * this.NR) * this.EPSMAC) * 1000;
        }
        else
        {
            num12 = (2 * this.NR) * Math.max(this.EPSDif, this.EPSMAC * 1000);
        } 
        this.QPTERM = 0;
        int num3 = 1;
        while (num3 <= this.NR)
        {
            for (num4 = 1;num4 <= this.N;num4++)
            {
                numArray11[num4] = this.GRES[this.ALIST[num3]][num4];
            }
            numArray12[num3] = 1;
            if (this.dN31PF(1, this.N, numArray11) == 0)
            {
                numArray12[num3] = 0;
            }
             
            num3++;
        }
        Label_026E:num3 = 1;
        while (num3 <= this.NDUAL)
        {
            this.DDUAL[num3] = 0;
            for (num4 = 1;num4 <= this.NDUAL;num4++)
            {
                this.R[num3, num4] = 0;
                this.XJ[num3, num4] = 0;
            }
            num3++;
        }
        this.RNORM = 1;
        this.RLOW = 1;
        double num30 = 0;
        num3 = 1;
        while (num3 <= this.NRES)
        {
            this.U[num3] = 0;
            if (this.W[num3] > num30)
            {
                num30 = this.W[num3];
            }
             
            num3++;
        }
        this.ACCINF[0x13, this.ITSTEP] = this.CLOW;
        this.ACCINF[20, this.ITSTEP] = num30;
        this.ACCINF[0x1f, this.ITSTEP] = this.TAUQP;
        num3 = 1;
        while (num3 <= (this._numberOfEqualityConstraints + this.MI))
        {
            this.UD[num3] = 0;
            num3++;
        }
        double num22 = Math.Abs(this.A[1, 1]);
        num3 = 1;
        while (num3 <= this.N)
        {
            num22 = Math.Max(num22, Math.Abs(this.A[num3, num3]));
            num3++;
        }
        num22 *= 10;
        num3 = 1;
        while (num3 <= this.N)
        {
            if (Math.Abs(this.A[num3, num3]) < (Math.Sqrt(this.RHO1) * num22))
            {
                this.A[num3, num3] = Math.Sqrt(this.RHO1) * num22;
            }
             
            num3++;
        }
        this.dN37PF(this.NX,this.N,this.A,num1,this.NDUAL,this.XJ);
        num22 = Math.Abs(this.A[1, 1]);
        int num8 = this.NR;
        double num23 = Math.Abs(this.XJ[1 + num8, 1 + num8]);
        num3 = 1;
        while (num3 <= this.N)
        {
            num22 = Math.Max(num22, Math.Abs(this.A[num3, num3]));
            num23 = Math.Max(num23, Math.Abs(this.XJ[num3 + num8, num3 + num8]));
            num3++;
        }
        double num13 = 0;
        num3 = 1;
        while (num3 <= this.N)
        {
            for (num4 = 1;num4 <= this.N;num4++)
            {
                num13 += this.A[num4, num3] * this.A[num4, num3];
            }
            num3++;
        }
        num13 /= (double)this.N;
        double num24 = 1 / Math.Sqrt(num13);
        num3 = 1;
        while (num3 <= num8)
        {
            this.XJ[num3, num3] = num24;
            num3++;
        }
        num3 = 1;
        while (num3 <= this.NDUAL)
        {
            if (num3 > num8)
            {
                numArray1[num3] = this.GRADF[num3 - num8] * this.SCF;
            }
             
            if (num3 <= num8)
            {
                numArray1[num3] = this.TAUQP;
            }
             
            num3++;
        }
        num3 = 1;
        while (num3 <= this.N)
        {
            num16 = 0;
            for (num4 = 1;num4 <= (num3 - 1);num4++)
            {
                num16 += this.A[num3, num4] * numArray11[num4 + num8];
            }
            numArray11[num3 + num8] = (numArray1[num3 + num8] - num16) / this.A[num3, num3];
            num3++;
        }
        num3 = this.N;
        while (num3 >= 1)
        {
            num16 = 0;
            for (num4 = num3 + 1;num4 <= this.N;num4++)
            {
                num16 += this.A[num4, num3] * numArray11[num4 + num8];
            }
            numArray11[num3 + num8] = (numArray11[num3 + num8] - num16) / this.A[num3, num3];
            num3--;
        }
        num3 = 1;
        while (num3 <= num8)
        {
            numArray11[num3] = 0;
            num3++;
        }
        num3 = 1;
        while (num3 <= this.NDUAL)
        {
            numArray5[num3] = -numArray11[num3];
            num3++;
        }
        double num20 = 0.5 * this.DN18PF(1, this.NDUAL, numArray1, numArray5);
        this.IQ = this.NR;
        num3 = 1;
        while (num3 <= this.IQ)
        {
            numArray13[num3] = num3;
            this.R[num3, num3] = 1;
            this.UD[num3] = this.TAUQP;
            num3++;
        }
        this.RNORM = 1;
        this.RLOW = 1;
        num3 = 1;
        while (num3 <= this._numberOfEqualityConstraints)
        {
            num4 = 1;
            while (num4 <= this.IQ)
            {
                this.UD1[num4] = this.UD[num4];
                num4++;
            }
            Label_0753:this.UD1[this.IQ + 1] = 0;
            numArray13[this.IQ + 1] = -this.NR - num3;
            num4 = 1;
            while (num4 <= this.N)
            {
                numArray4[num4 + num8] = this.GRES[num3, num4];
                num4++;
            }
            num4 = 1;
            while (num4 <= num8)
            {
                numArray4[num4] = 0;
                num4++;
            }
            numArray4[num3] = 1;
            if (this.RES[num3] > 0)
            {
                numArray4[num3] = -1;
            }
             
            num4 = 1;
            while (num4 <= this.NDUAL)
            {
                this.NP[num4] = numArray4[num4];
                num4++;
            }
            this.DN33PF(numArray7);
            if (this.IQ != 0)
            {
                this.DN34PF(numArray8);
            }
             
            for (int num32 = 1;num32 <= 8;num32++)
            {
            }
            num14 = this.DN31PF(1, this.NDUAL, numArray7);
            num25 = this.DN18PF(1, this.NDUAL, numArray7, this.NP);
            if ((num14 >= (num12 * this.RNORM)) && (num25 > 0))
            {
                num19 = (-this.DN18PF(1, this.NDUAL, this.NP, numArray5) - this.RES[num3]) / num25;
            }
            else if ((-this.DN18PF(1, this.NDUAL, this.NP, numArray5) - this.RES[num3]) >= 0)
            {
                num19 = num29;
            }
            else
            {
                num19 = -num29;
            }  
            if (this.IQ != 0)
            {
                this.DN34PF(numArray8);
            }
             
            num7 = 0;
            if (num19 > 0)
            {
                num18 = num29;
                for (num5 = 1;num5 <= this.IQ;num5++)
                {
                    if (((numArray8[num5] > 0) && (numArray13[num5] > 0)) && ((this.UD1[num5] / numArray8[num5]) < num18))
                    {
                        num18 = this.UD1[num5] / numArray8[num5];
                    }
                     
                }
                num17 = Math.Min(num18, num19);
            }
            else
            {
                num18 = num29;
                for (num5 = 1;num5 <= this.IQ;num5++)
                {
                    if (((numArray8[num5] < 0) && (numArray13[num5] > 0)) && ((this.UD1[num5] / Math.Abs(numArray8[num5])) < num18))
                    {
                        num18 = this.UD1[num5] / Math.Abs(numArray8[num5]);
                    }
                     
                }
                num18 = -num18;
                num17 = Math.Max(num18, num19);
            } 
            if (Math.Abs(num17) >= num29)
            {
                this.QPTERM = -2;
                this.ACCINF[30, this.ITSTEP] = -2;
                this.ACCINF[13, this.ITSTEP] = num28;
                this.ACCINF[14, this.ITSTEP] = num22 * num23;
                num3 = 1;
                while (num3 <= this.N)
                {
                    this.D[num3] = numArray5[num3 + num8];
                    num3++;
                }
                this.DNORM = this.dN31PF(1,this.N,this.D);
                num26 = 0;
                num3 = 1;
                while (num3 <= num8)
                {
                    num26 += Math.Abs(numArray5[num3]);
                    num3++;
                }
                this.ACCINF[0x20, this.ITSTEP] = num26;
                this.INFEAS = num26;
                flag1 = false;
                num26 = 0;
                num27 = 0;
                num3 = 1;
                while (num3 <= this.IQ)
                {
                    if (numArray13[num3] < 0)
                    {
                        this.U[-(numArray13[num3] + this.NR)] = this.UD[num3];
                    }
                    else if (numArray13[num3] > this.NR)
                    {
                        this.U[this.ALIST[(numArray13[num3] - this.NR) + this.NH]] = this.UD[num3];
                    }
                      
                    num3++;
                }
                num30 = 0;
                num4 = 1;
                while (num4 <= this.N)
                {
                    this.NP[num4] = this.GRADF[num4] * this.SCF;
                    num4++;
                }
                num3 = 1;
                while (num3 <= this.NRES)
                {
                    for (num4 = 1;num4 <= this.N;num4++)
                    {
                        this.NP[num4] -= this.GRES[num3, num4] * this.U[num3];
                    }
                    this.W1[num3] = this.W[num3];
                    if (num3 <= this.NH)
                    {
                        if (this.SLACK[num3] > Math.Abs(this.RES[num3]))
                        {
                            this.W1[num3] = Math.Abs(this.U[num3]);
                        }
                         
                        if (this.SLACK[num3] <= Math.Abs(this.RES[num3]))
                        {
                            if ((this.W[num3] <= Math.Abs(this.U[num3])) && (Math.Abs(this.U[num3]) <= (this.W[num3] + this.TAU)))
                            {
                                this.W1[num3] = this.W[num3] + (2 * this.TAU);
                            }
                            else
                            {
                                this.W1[num3] = Math.Max(this.W[num3], (this.NY * Math.Abs(this.U[num3])) + this.TAU);
                            } 
                        }
                         
                        num26 += Math.Abs(this.RES[num3]) * this.W1[num3];
                        num27 += Math.Abs(this.RESST[num3]) * this.W1[num3];
                    }
                    else
                    {
                        if ((this.SLACK[num3] > -Math.Min(-num12, this.RES[num3])) && (this.BIND[num3] == 1))
                        {
                            this.W1[num3] = Math.Abs(this.U[num3]);
                        }
                        else if (((this.BIND[num3] == 1) && (this.SLACK[num3] <= -Math.Min(-num12, this.RES[num3]))) && ((this.U[num3] <= (this.W[num3] + this.TAU)) && (this.W[num3] >= this.U[num3])))
                        {
                            this.W1[num3] = this.W[num3] + (2 * this.TAU);
                        }
                        else if (this.BIND[num3] == 1)
                        {
                            this.W1[num3] = Math.Max(this.W[num3], (this.NY * Math.Abs(this.U[num3])) + this.TAU);
                        }
                           
                        num26 -= this.W1[num3] * Math.Min(0, this.RES[num3]);
                        num27 -= this.W1[num3] * Math.Min(0, this.RESST[num3]);
                    } 
                    if (this.W[num3] != this.W1[num3])
                    {
                        this.LASTCH = this.ITSTEP;
                    }
                     
                    this.W[num3] = this.W1[num3];
                    num30 = Math.Max(num30, this.W[num3]);
                    num3++;
                }
                this.PSIST = num27;
                this.PSI = num26;
                this.B2N = Math.Sqrt(this.dN18PF(1,this.N,this.NP,this.NP));
                if (this.SCF == 0)
                {
                    this.B2N = -1;
                }
                else
                {
                    this.B2N /= this.SCF;
                } 
                if (flag1)
                {
                    this.CLOW++;
                    this.LASTCH = this.ITSTEP;
                    this.LASTDW = this.ITSTEP;
                }
                 
                if (!this.SILENT)
                {
                    this.dN44PF(12);
                }
                 
                this.ACCINF[0x13, this.ITSTEP] = this.CLOW;
                this.ACCINF[20, this.ITSTEP] = num30;
                this.ACCINF[0x1f, this.ITSTEP] = this.TAUQP;
                this.dN12PF();
                this.dN11PF();
                if (this.DIRDER < 0)
                {
                    if (-this.DIRDER > ((this.EPSMAC * 100) * (((this.SCF * Math.Abs(this.FX)) + this.PSI) + 1)))
                    {
                        return ;
                    }
                     
                    if (this.INFEAS <= Math.Max(this.UPSI, this.NRES * this.DELMIN))
                    {
                        return ;
                    }
                     
                }
                 
                if (this.TAUQP > (this.TAUMAX / this.TAUFAC))
                {
                    this.QPTERM = -1;
                    this.ACCINF[30, this.ITSTEP] = this.QPTERM;
                    this.ACCINF[0x1f, this.ITSTEP] = this.TAUQP;
                    this.ACCINF[0x20, this.ITSTEP] = this.INFEAS;
                    for (num3 = 1;num3 <= this.N;num3++)
                    {
                        this.D[num3] = 0;
                    }
                    this.DNORM = 0;
                    this.DIRDER = 0;
                    return ;
                }
                 
                this.TAUQP *= this.TAUFAC;
                goto Label_026E
            }
             
            if (Math.Abs(num19) >= num29)
            {
                num5 = 1;
                while (num5 <= this.IQ)
                {
                    this.UD1[num5] += num17 * -numArray8[num5];
                    if ((this.UD1[num5] < 0) && (numArray13[num5] > 0))
                    {
                        this.UD1[num5] = 0;
                    }
                     
                    num5++;
                }
                this.UD1[this.IQ + 1] += num17;
                numArray17[0] = 0;
                for (num4 = 1;num4 <= this.IQ;num4++)
                {
                    if ((this.UD1[num4] <= num12) && (numArray13[num4] > 0))
                    {
                        numArray17[0]++;
                        numArray17[numArray17[0]] = numArray13[num4];
                    }
                     
                }
                for (num5 = 1;num5 <= numArray17[0];num5++)
                {
                    num7 = numArray17[num5];
                    numArray14[num7] = num7;
                    this.DN35PF(numArray13, num7);
                }
                goto Label_0753
            }
             
            num5 = 1;
            while (num5 <= this.NDUAL)
            {
                numArray5[num5] += num17 * numArray7[num5];
                num5++;
            }
            num5 = 1;
            while (num5 <= this.IQ)
            {
                this.UD1[num5] += num17 * -numArray8[num5];
                if ((this.UD1[num5] < 0) && (numArray13[num5] > 0))
                {
                    this.UD1[num5] = 0;
                }
                 
                num5++;
            }
            this.UD1[this.IQ + 1] += num17;
            num20 += (num17 * this.DN18PF(1, this.NDUAL, numArray7, this.NP)) * ((0.5 * num17) + this.UD1[this.IQ + 1]);
            if (Math.Abs((double)(num19 - num18)) <= num12)
            {
                numArray17[0] = 0;
                num4 = 1;
                while (num4 <= this.IQ)
                {
                    if ((this.UD1[num4] <= num12) && (numArray13[num4] > 0))
                    {
                        numArray17[0]++;
                        numArray17[numArray17[0]] = numArray13[num4];
                    }
                     
                    num4++;
                }
                for (num5 = 1;num5 <= numArray17[0];num5++)
                {
                    num7 = numArray17[num5];
                    numArray14[num7] = num7;
                    this.DN35PF(numArray13, num7);
                }
                numArray13[this.IQ + 1] = -num3 - this.NR;
                this.dN36PF();
                for (num4 = 1;num4 <= this.IQ;num4++)
                {
                    this.UD[num4] = this.UD1[num4];
                }
            }
            else
            {
                if (num17 != num19)
                {
                    numArray17[0] = 0;
                    for (num4 = 0;num4 <= this.IQ;num4++)
                    {
                        if ((this.UD1[num4] <= num12) && (numArray13[num4] > 0))
                        {
                            numArray17[0]++;
                            numArray17[numArray17[0]] = numArray13[num4];
                        }
                         
                    }
                    for (num5 = 1;num5 <= numArray17[0];num5++)
                    {
                        num7 = numArray17[num5];
                        numArray14[num7] = num7;
                        this.DN35PF(numArray13, num7);
                    }
                    goto Label_0753
                }
                 
                numArray13[this.IQ + 1] = -num3 - this.NR;
                this.dN36PF();
                for (num4 = 1;num4 <= this.IQ;num4++)
                {
                    this.UD[num4] = this.UD1[num4];
                }
            } 
            num3++;
        }
        num3 = 1;
        while (num3 <= this.MI)
        {
            numArray14[num3] = num3;
            num3++;
        }
        Label_13E3:num3 = 1;
        while (num3 <= this.IQ)
        {
            num6 = numArray13[num3];
            if (num6 > 0)
            {
                numArray14[num6] = 0;
            }
             
            num3++;
        }
        double num21 = 0;
        num3 = 1;
        while (num3 <= this.MI)
        {
            numArray15[num3] = 1;
            num16 = 0;
            if (num3 > this.NR)
            {
                num5 = this.ALIST[(num3 + this.NH) - num8];
                num4 = 1;
                while (num4 <= this.N)
                {
                    numArray3[num4 + num8] = this.GRES[num5, num4];
                    num4++;
                }
                for (num4 = 1;num4 <= num8;num4++)
                {
                    numArray3[num4] = 0;
                }
                numArray3[(this.NH + num3) - num8] = 1;
                numArray2[num3] = this.RES[num5];
            }
            else
            {
                for (num4 = 1;num4 <= this.NDUAL;num4++)
                {
                    numArray3[num4] = 0;
                }
                numArray2[num3] = 0;
                numArray3[num3] = 1;
            } 
            num16 = this.DN18PF(1, this.NDUAL, numArray3, numArray5) + numArray2[num3];
            numArray6[num3] = num16;
            num21 += Math.Min(0, num16);
            num3++;
        }
        num3 = 1;
        while (num3 <= this.IQ)
        {
            numArray10[num3] = this.UD[num3];
            numArray16[num3] = numArray13[num3];
            num3++;
        }
        num3 = 1;
        while (num3 <= this.NDUAL)
        {
            numArray9[num3] = numArray5[num3];
            num3++;
        }
        Label_1568:num15 = 0;
        num6 = 0;
        num3 = 1;
        while (num3 <= this.MI)
        {
            if (((numArray6[num3] < num15) && (numArray14[num3] != 0)) && (numArray15[num3] != 0))
            {
                num15 = numArray6[num3];
                num6 = num3;
            }
             
            num3++;
        }
        if (this.IQ > 1)
        {
            num28 = this.RNORM / this.RLOW;
        }
        else
        {
            num28 = 1;
        } 
        if ((Math.Abs(num21) <= (num12 * ((num22 * num23) + num28))) || (num6 == 0))
        {
            double num10;
            double num11;
            this.QPTERM = 1;
            this.ACCINF[30, this.ITSTEP] = 1;
            this.ACCINF[13, this.ITSTEP] = num28;
            this.ACCINF[14, this.ITSTEP] = num22 * num23;
            num3 = 1;
            while (num3 <= this.N)
            {
                this.D[num3] = numArray5[num3 + num8];
                num3++;
            }
            this.DNORM = this.dN31PF(1,this.N,this.D);
            this.INFEAS = 0;
            num3 = 1;
            while (num3 <= num8)
            {
                this.INFEAS += Math.Abs(numArray5[num3]);
                num3++;
            }
            this.ACCINF[0x1f, this.ITSTEP] = this.TAUQP;
            this.ACCINF[0x20, this.ITSTEP] = this.INFEAS;
            flag1 = false;
            num26 = 0;
            num27 = 0;
            num3 = 1;
            while (num3 <= this.IQ)
            {
                if (numArray13[num3] < 0)
                {
                    this.U[-(numArray13[num3] + this.NR)] = this.UD[num3];
                }
                else if (numArray13[num3] > this.NR)
                {
                    this.U[this.ALIST[(numArray13[num3] + this.NH) - this.NR]] = this.UD[num3];
                }
                  
                num3++;
            }
            num30 = 0;
            num4 = 1;
            while (num4 <= this.N)
            {
                this.NP[num4] = this.GRADF[num4] * this.SCF;
                num4++;
            }
            num3 = 1;
            while (num3 <= this.NRES)
            {
                for (num4 = 1;num4 <= this.N;num4++)
                {
                    this.NP[num4] -= this.GRES[num3, num4] * this.U[num3];
                }
                num3++;
            }
            this.B2N = this.dN31PF(1,this.N,this.NP);
            if (this.SCF != 0)
            {
                this.B2N /= this.SCF;
            }
             
            double num9 = 0;
            num3 = 1;
            while (num3 <= this.NR)
            {
                num9 += Math.Abs(numArray5[num3]) * numArray12[num3];
                num3++;
            }
            if (((this.UPSI <= (this.DELMIN * this.NRES)) && (this.B2N <= (((this.GFN + 1) * this.EPSX) * 100))) && ((this.PHASE >= 0) && (this.INFEAS <= (this.DELMIN * this.NRES))))
            {
                for (num3 = 1;num3 <= this.N;num3++)
                {
                    this.D[num3] = 0;
                }
                this.DNORM = 0;
                this.QPTERM = 1;
                this.DIRDER = 0;
                this.OPTITE = 3;
                return ;
            }
             
            if ((num9 > ((1 - (this.DELTA1 / this.TAUQP)) * this.UPSI)) && ((this.dN31PF(1,this.N,this.D) <= (Math.Min(num9, num9 * num9) * 10)) || (this.UPSI > (this.TAU0 * 0.5))))
            {
                num3 = 1;
                while (num3 <= this.NRES)
                {
                    this.U[num3] = 0;
                    this.SLACK[num3] = 0;
                    num3++;
                }
                if ((this.TAUQP * this.TAUFAC) > this.TAUMAX)
                {
                    this.QPTERM = -1;
                    this.ACCINF[30, this.ITSTEP] = this.QPTERM;
                    this.ACCINF[0x1f, this.ITSTEP] = this.TAUQP;
                    this.ACCINF[0x20, this.ITSTEP] = this.INFEAS;
                    for (num3 = 1;num3 <= this.N;num3++)
                    {
                        this.D[num3] = 0;
                    }
                    this.DNORM = 0;
                    return ;
                }
                 
                this.TAUQP *= this.TAUFAC;
                goto Label_026E
            }
             
            num3 = 1;
            while (num3 <= this.NRES)
            {
                this.SLACK[num3] = 0;
                num3++;
            }
            num3 = 1;
            while (num3 <= this.NR)
            {
                this.SLACK[this.ALIST[num3]] = numArray5[num3];
                num3++;
            }
            flag1 = false;
            num3 = 1;
            while (num3 <= this.NRES)
            {
                this.W1[num3] = this.W[num3];
                if (num3 <= this.NH)
                {
                    if (Math.Abs(this.SLACK[num3]) > (Math.Abs(this.RES[num3]) + num12))
                    {
                        this.W1[num3] = Math.Abs(this.U[num3]);
                    }
                    else
                    {
                        this.W1[num3] = (this.NY * Math.Abs(this.U[num3])) + this.TAU;
                    } 
                }
                else if (this.BIND[num3] == 0)
                {
                    this.W1[num3] = Math.Max(this.W[num3] * 0.8, this.TAU);
                }
                else
                {
                    if ((this.RES[num3] >= 0) && (this.SLACK[num3] <= num12))
                    {
                        this.W1[num3] = Math.Max((double)((this.NY * Math.Abs(this.U[num3])) + this.TAU), (double)((Math.Abs(this.U[num3]) + this.W1[num3]) * 0.5));
                    }
                     
                    if ((this.RES[num3] >= 0) && (this.SLACK[num3] > num12))
                    {
                        this.W1[num3] = Math.Abs(this.U[num3]);
                    }
                     
                    if ((this.RES[num3] < 0) && (this.SLACK[num3] <= (-this.RES[num3] + num12)))
                    {
                        this.W1[num3] = Math.Max((double)((this.NY * Math.Abs(this.U[num3])) + this.TAU), (double)((Math.Abs(this.U[num3]) + this.W1[num3]) * 0.5));
                    }
                     
                    if ((this.RES[num3] < 0) && (this.SLACK[num3] > (-this.RES[num3] + num12)))
                    {
                        this.W1[num3] = Math.Abs(this.U[num3]);
                    }
                     
                }  
                if (this.W1[num3] < this.W[num3])
                {
                    flag1 = true;
                }
                 
                num3++;
            }
            if (flag1)
            {
                num10 = 0;
                num11 = 0;
                num3 = 1;
                while (num3 <= this.NRES)
                {
                    if (num3 <= this.NH)
                    {
                        num10 += this.W1[num3] * Math.Abs(this.RESST[num3]);
                        num11 += this.W1[num3] * Math.Abs(this.RES[num3]);
                    }
                    else
                    {
                        num10 -= Math.Min(0, this.RESST[num3]) * this.W1[num3];
                        num11 -= Math.Min(0, this.RES[num3]) * this.W1[num3];
                    } 
                    num3++;
                }
                double num31 = ((this.FXST - this.FX) * this.SCF) + (num10 - num11);
                if ((num31 >= (this.ETA * this.CLOW)) && ((this.ITSTEP - this.LASTDW) >= Math.Max(5, Math.Min(this.N / 10, 20))))
                {
                    this.LASTDW = this.ITSTEP;
                    this.LASTCH = this.ITSTEP;
                    this.LEVEL = num31 / ((double)this.ITERMA);
                    this.PSIST = num10;
                    this.PSI = num11;
                    num3 = 1;
                    while (num3 <= this.NRES)
                    {
                        if (this.W1[num3] != this.W[num3])
                        {
                            this.LASTCH = this.ITSTEP;
                        }
                         
                        this.W[num3] = this.W1[num3];
                        num3++;
                    }
                    this.CLOW++;
                    if (this.CLOW > (this.ITSTEP / 10))
                    {
                        this.ETA *= 1.3;
                        if (!this.SILENT)
                        {
                            this.dN44PF(11);
                        }
                         
                    }
                     
                    if (!this.SILENT)
                    {
                        this.dN44PF(12);
                    }
                     
                    this.dN12PF();
                    this.dN11PF();
                    if (this.DIRDER < 0)
                    {
                        if (-this.DIRDER > ((this.EPSMAC * 100) * (((this.SCF * Math.Abs(this.FX)) + this.PSI) + 1)))
                        {
                            return ;
                        }
                         
                        if (this.INFEAS <= Math.Max(this.UPSI, this.NRES * this.DELMIN))
                        {
                            return ;
                        }
                         
                    }
                     
                    if (this.TAUQP > (this.TAUMAX / this.TAUFAC))
                    {
                        this.QPTERM = -1;
                        this.ACCINF[30, this.ITSTEP] = this.QPTERM;
                        this.ACCINF[0x1f, this.ITSTEP] = this.TAUQP;
                        this.ACCINF[0x20, this.ITSTEP] = this.INFEAS;
                        for (num3 = 1;num3 <= this.N;num3++)
                        {
                            this.D[num3] = 0;
                        }
                        this.DNORM = 0;
                        this.DIRDER = 0;
                        return ;
                    }
                     
                    this.TAUQP *= this.TAUFAC;
                    goto Label_026E
                }
                 
            }
             
            num3 = 1;
            while (num3 <= this.NRES)
            {
                this.W1[num3] = this.W[num3];
                if (num3 <= this.NH)
                {
                    if (this.SLACK[num3] > Math.Abs(this.RES[num3]))
                    {
                        this.W1[num3] = Math.Abs(this.U[num3]);
                    }
                     
                    if (this.SLACK[num3] <= Math.Abs(this.RES[num3]))
                    {
                        if ((this.W[num3] <= Math.Abs(this.U[num3])) && (Math.Abs(this.U[num3]) <= (this.W[num3] + this.TAU)))
                        {
                            this.W1[num3] = this.W[num3] + (2 * this.TAU);
                        }
                        else
                        {
                            this.W1[num3] = Math.Max(this.W[num3], (this.NY * Math.Abs(this.U[num3])) + this.TAU);
                        } 
                    }
                     
                }
                else if ((this.SLACK[num3] > -Math.Min(-num12, this.RES[num3])) && (this.BIND[num3] == 1))
                {
                    this.W1[num3] = Math.Abs(this.U[num3]);
                }
                else if (((this.BIND[num3] == 1) && (this.SLACK[num3] <= -Math.Min(-num12, this.RES[num3]))) && ((this.U[num3] <= (this.W[num3] + this.TAU)) && (this.W[num3] >= this.U[num3])))
                {
                    this.W1[num3] = this.W[num3] + (2 * this.TAU);
                }
                else if (this.BIND[num3] == 1)
                {
                    this.W1[num3] = Math.Max(this.W[num3], (this.NY * Math.Abs(this.U[num3])) + this.TAU);
                }
                    
                num3++;
            }
            num30 = 0;
            num3 = 1;
            while (num3 <= this.NRES)
            {
                if ((this.W1[num3] > this.W[num3]) || (this.W1[num3] < this.W[num3]))
                {
                    this.LASTCH = this.ITSTEP;
                }
                 
                if (this.W1[num3] > this.W[num3])
                {
                    this.LASTUP = this.ITSTEP;
                }
                 
                if (this.W1[num3] < this.W[num3])
                {
                    this.LASTDW = this.ITSTEP;
                }
                 
                this.W[num3] = this.W1[num3];
                num30 = Math.Max(num30, this.W[num3]);
                num3++;
            }
            num10 = 0;
            num11 = 0;
            num3 = 1;
            while (num3 <= this.NRES)
            {
                if (num3 <= this.NH)
                {
                    num10 += this.W[num3] * Math.Abs(this.RESST[num3]);
                    num11 += this.W[num3] * Math.Abs(this.RES[num3]);
                }
                else
                {
                    num10 -= this.W[num3] * Math.Min(0, this.RESST[num3]);
                    num11 -= this.W[num3] * Math.Min(0, this.RES[num3]);
                } 
                num3++;
            }
            this.PSIST = num10;
            this.PSI = num11;
            if (!this.SILENT)
            {
                this.dN44PF(12);
            }
             
            this.ACCINF[20, this.ITSTEP] = num30;
            this.ACCINF[0x13, this.ITSTEP] = this.CLOW;
            this.dN12PF();
            this.dN11PF();
            if (this.DIRDER < 0)
            {
                if (-this.DIRDER > ((this.EPSMAC * 100) * (((this.SCF * Math.Abs(this.FX)) + this.PSI) + 1)))
                {
                    return ;
                }
                 
                if (this.INFEAS <= Math.Max(this.UPSI, this.NRES * this.DELMIN))
                {
                    return ;
                }
                 
            }
             
            if (this.TAUQP > (this.TAUMAX / this.TAUFAC))
            {
                this.QPTERM = -1;
                this.ACCINF[30, this.ITSTEP] = this.QPTERM;
                this.ACCINF[0x1f, this.ITSTEP] = this.TAUQP;
                this.ACCINF[0x20, this.ITSTEP] = this.INFEAS;
                for (num3 = 1;num3 <= this.N;num3++)
                {
                    this.D[num3] = 0;
                }
                this.DNORM = 0;
                this.DIRDER = 0;
                return ;
            }
             
            this.TAUQP *= this.TAUFAC;
            goto Label_026E
        }
         
        if (num6 > this.NR)
        {
            num5 = this.ALIST[(num6 + this.NH) - this.NR];
            num4 = 0;
            while (num4 <= this.N)
            {
                numArray3[num4 + num8] = this.GRES[num5, num4];
                num4++;
            }
            for (num4 = 1;num4 <= num8;num4++)
            {
                numArray3[num4] = 0;
            }
            numArray3[(this.NH + num6) - this.NR] = 1;
            numArray2[num6] = this.RES[num5];
        }
        else
        {
            for (num4 = 1;num4 <= this.NDUAL;num4++)
            {
                numArray3[num4] = 0;
            }
            numArray2[num6] = 0;
            numArray3[num6] = 1;
        } 
        num3 = 1;
        while (num3 <= this.NDUAL)
        {
            this.NP[num3] = numArray3[num3];
            num3++;
        }
        num3 = 1;
        while (num3 <= this.IQ)
        {
            this.UD1[num3] = this.UD[num3];
            num3++;
        }
        this.UD1[this.IQ + 1] = 0;
        numArray13[this.IQ + 1] = num6;
        Label_2607:this.DN33PF(numArray7);
        if (this.IQ != 0)
        {
            this.DN34PF(numArray8);
        }
         
        num7 = 0;
        num18 = num29;
        num5 = 1;
        while (num5 <= this.IQ)
        {
            if (((numArray13[num5] > 0) && (numArray8[num5] > 0)) && ((this.UD1[num5] / numArray8[num5]) < num18))
            {
                num18 = this.UD1[num5] / numArray8[num5];
            }
             
            num5++;
        }
        num14 = this.DN31PF(1, this.NDUAL, numArray7);
        num25 = this.DN18PF(1, this.NDUAL, numArray7, this.NP);
        if ((num14 >= (num12 * this.RNORM)) && (num25 > 0))
        {
            num19 = -numArray6[num6] / num25;
        }
        else
        {
            num19 = num29;
        } 
        num17 = Math.Min(num18, num19);
        if (num17 >= num29)
        {
            this.QPTERM = -2;
            this.ACCINF[30, this.ITSTEP] = -2;
            this.ACCINF[13, this.ITSTEP] = num28;
            this.ACCINF[14, this.ITSTEP] = num22 * num23;
            num3 = 1;
            while (num3 <= this.N)
            {
                this.D[num3] = numArray5[num3 + num8];
                num3++;
            }
            this.DNORM = this.dN31PF(1,this.N,this.D);
            num26 = 0;
            num3 = 1;
            while (num3 <= num8)
            {
                num26 += Math.Abs(numArray5[num3]);
                num3++;
            }
            this.ACCINF[0x20, this.ITSTEP] = num26;
            this.INFEAS = num26;
            flag1 = false;
            num26 = 0;
            num27 = 0;
            num3 = 1;
            while (num3 <= this.IQ)
            {
                if (numArray13[num3] < 0)
                {
                    this.U[-(numArray13[num3] + this.NR)] = this.UD[num3];
                }
                else if (numArray13[num3] > this.NR)
                {
                    this.U[this.ALIST[(numArray13[num3] - this.NR) + this.NH]] = this.UD[num3];
                }
                  
                num3++;
            }
            num30 = 0;
            num4 = 1;
            while (num4 <= this.N)
            {
                this.NP[num4] = this.GRADF[num4] * this.SCF;
                num4++;
            }
            num3 = 1;
            while (num3 <= this.NRES)
            {
                for (num4 = 1;num4 <= this.N;num4++)
                {
                    this.NP[num4] -= this.GRES[num3, num4] * this.U[num3];
                }
                this.W1[num3] = this.W[num3];
                if (num3 <= this.NH)
                {
                    if (this.SLACK[num3] > Math.Abs(this.RES[num3]))
                    {
                        this.W1[num3] = Math.Abs(this.U[num3]);
                    }
                     
                    if (this.SLACK[num3] <= Math.Abs(this.RES[num3]))
                    {
                        if ((this.W[num3] <= Math.Abs(this.U[num3])) && (Math.Abs(this.U[num3]) <= (this.W[num3] + this.TAU)))
                        {
                            this.W1[num3] = this.W[num3] + (2 * this.TAU);
                        }
                        else
                        {
                            this.W1[num3] = Math.Max(this.W[num3], (this.NY * Math.Abs(this.U[num3])) + this.TAU);
                        } 
                    }
                     
                    num26 += Math.Abs(this.RES[num3]) * this.W1[num3];
                    num27 += Math.Abs(this.RESST[num3]) * this.W1[num3];
                }
                else
                {
                    if ((this.SLACK[num3] > -Math.Min(-num12, this.RES[num3])) && (this.BIND[num3] == 1))
                    {
                        this.W1[num3] = Math.Abs(this.U[num3]);
                    }
                    else if (((this.BIND[num3] == 1) && (this.SLACK[num3] <= -Math.Min(-num12, this.RES[num3]))) && ((this.U[num3] <= (this.W[num3] + this.TAU)) && (this.W[num3] >= this.U[num3])))
                    {
                        this.W1[num3] = this.W[num3] + (2 * this.TAU);
                    }
                    else if (this.BIND[num3] == 1)
                    {
                        this.W1[num3] = Math.Max(this.W[num3], (this.NY * Math.Abs(this.U[num3])) + this.TAU);
                    }
                       
                    num26 -= this.W1[num3] * Math.Min(0, this.RES[num3]);
                    num27 -= this.W1[num3] * Math.Min(0, this.RESST[num3]);
                } 
                if (this.W[num3] != this.W1[num3])
                {
                    this.LASTCH = this.ITSTEP;
                }
                 
                this.W[num3] = this.W1[num3];
                num30 = Math.Max(num30, this.W[num3]);
                num3++;
            }
            this.PSIST = num27;
            this.PSI = num26;
            this.B2N = Math.Sqrt(this.dN18PF(1,this.N,this.NP,this.NP));
            if (this.SCF == 0)
            {
                this.B2N = -1;
            }
            else
            {
                this.B2N /= this.SCF;
            } 
            if (flag1)
            {
                this.CLOW++;
                this.LASTCH = this.ITSTEP;
                this.LASTDW = this.ITSTEP;
            }
             
            if (!this.SILENT)
            {
                this.dN44PF(12);
            }
             
            this.ACCINF[0x13, this.ITSTEP] = this.CLOW;
            this.ACCINF[20, this.ITSTEP] = num30;
            this.ACCINF[0x1f, this.ITSTEP] = this.TAUQP;
            this.dN12PF();
            this.dN11PF();
            if (this.DIRDER < 0)
            {
                if (-this.DIRDER > ((this.EPSMAC * 100) * (((this.SCF * Math.Abs(this.FX)) + this.PSI) + 1)))
                {
                    return ;
                }
                 
                if (this.INFEAS <= Math.Max(this.UPSI, this.NRES * this.DELMIN))
                {
                    return ;
                }
                 
            }
             
            if (this.TAUQP > (this.TAUMAX / this.TAUFAC))
            {
                this.QPTERM = -1;
                this.ACCINF[30, this.ITSTEP] = this.QPTERM;
                this.ACCINF[0x1f, this.ITSTEP] = this.TAUQP;
                this.ACCINF[0x20, this.ITSTEP] = this.INFEAS;
                for (num3 = 1;num3 <= this.N;num3++)
                {
                    this.D[num3] = 0;
                }
                this.DNORM = 0;
                this.DIRDER = 0;
                return ;
            }
             
            this.TAUQP *= this.TAUFAC;
            goto Label_026E
        }
         
        if (num19 >= num29)
        {
            num5 = 1;
            while (num5 <= this.IQ)
            {
                this.UD1[num5] += num17 * -numArray8[num5];
                if ((this.UD1[num5] < 0) && (numArray13[num5] > 0))
                {
                    this.UD1[num5] = 0;
                }
                 
                num5++;
            }
            this.UD1[this.IQ + 1] += num17;
            numArray17[0] = 0;
            for (num3 = 1;num3 <= this.IQ;num3++)
            {
                if ((this.UD1[num3] <= num12) && (numArray13[num3] > 0))
                {
                    numArray17[0]++;
                    numArray17[numArray17[0]] = numArray13[num3];
                }
                 
            }
            for (num5 = 1;num5 <= numArray17[0];num5++)
            {
                num7 = numArray17[num5];
                numArray14[num7] = num7;
                this.DN35PF(numArray13, num7);
            }
            goto Label_2607
        }
         
        num5 = 1;
        while (num5 <= this.NDUAL)
        {
            numArray5[num5] += num17 * numArray7[num5];
            num5++;
        }
        num5 = 1;
        while (num5 <= this.IQ)
        {
            this.UD1[num5] += num17 * -numArray8[num5];
            if ((this.UD1[num5] < 0) && (numArray13[num5] > 0))
            {
                this.UD1[num5] = 0;
            }
             
            num5++;
        }
        this.UD1[this.IQ + 1] += num17;
        num20 += (num17 * this.DN18PF(1, this.NDUAL, numArray7, this.NP)) * ((0.5 * num17) + this.UD1[this.IQ + 1]);
        if (num19 <= (num18 - num12))
        {
            if (this.dN31PF(this.IQ + 1,this.NDUAL,this.DDUAL) >= (this.EPSMAC * this.RNORM))
            {
                this.dN36PF();
                numArray14[num6] = 0;
                for (num3 = 1;num3 <= this.IQ;num3++)
                {
                    this.UD[num3] = this.UD1[num3];
                }
                goto Label_13E3
            }
             
            this.IPTR = num6;
            this.IQTR = this.IQ;
            num3 = 1;
            while (num3 <= this.IQ)
            {
                this.AITR[num3] = numArray13[num3];
                num3++;
            }
            this.SSTR = num15;
            this.RIITR = this.dN31PF(this.IQ + 1,this.NDUAL,this.DDUAL);
            numArray15[num6] = 0;
            num3 = 1;
            while (num3 <= this.MI)
            {
                numArray14[num3] = num3;
                num3++;
            }
            num3 = 1;
            while (num3 <= this.IQ)
            {
                numArray13[num3] = numArray16[num3];
                if (numArray13[num3] > 0)
                {
                    numArray14[numArray13[num3]] = 0;
                }
                 
                this.UD1[num3] = numArray10[num3];
                num3++;
            }
            for (num3 = 1;num3 <= this.NDUAL;num3++)
            {
                numArray5[num3] = numArray9[num3];
            }
            goto Label_1568
        }
         
        num16 = 0;
        if (num6 > this.NR)
        {
            num5 = this.ALIST[(num6 + this.NH) - this.NR];
            num4 = 1;
            while (num4 <= this.N)
            {
                numArray3[num4 + num8] = this.GRES[num5, num4];
                num4++;
            }
            for (num4 = 1;num4 <= num8;num4++)
            {
                numArray3[num4] = 0;
            }
            numArray3[(this.NH + num6) - this.NR] = 1;
            numArray2[num6] = this.RES[num5];
            numArray6[num6] = this.DN18PF(1, this.NDUAL, numArray3, numArray5) + numArray2[num6];
        }
        else
        {
            numArray6[num6] = numArray5[num6];
        } 
        numArray17[0] = 0;
        num3 = 1;
        while (num3 <= this.IQ)
        {
            if ((this.UD1[num3] <= num12) && (numArray13[num3] > 0))
            {
                numArray17[0]++;
                numArray17[numArray17[0]] = numArray13[num3];
            }
             
            num3++;
        }
        for (num5 = 1;num5 <= numArray17[0];num5++)
        {
            num7 = numArray17[num5];
            numArray14[num7] = num7;
            this.DN35PF(numArray13, num7);
        }
        if (num19 > (num18 + num12))
        {
            goto Label_2607
        }
         
        if (this.dN31PF(this.IQ + 1,this.NDUAL,this.DDUAL) >= (this.EPSMAC * this.RNORM))
        {
            this.dN36PF();
            numArray14[num6] = 0;
            for (num3 = 1;num3 <= this.IQ;num3++)
            {
                this.UD[num3] = this.UD1[num3];
            }
            goto Label_13E3
        }
         
        this.IPTR = num6;
        this.IQTR = this.IQ;
        num3 = 1;
        while (num3 <= this.IQ)
        {
            this.AITR[num3] = numArray13[num3];
            num3++;
        }
        this.SSTR = num15;
        this.RIITR = this.dN31PF(this.IQ + 1,this.NDUAL,this.DDUAL);
        numArray15[num6] = 0;
        num3 = 1;
        while (num3 <= this.MI)
        {
            numArray14[num3] = num3;
            num3++;
        }
        num3 = 1;
        while (num3 <= this.IQ)
        {
            numArray13[num3] = numArray16[num3];
            if (numArray13[num3] > 0)
            {
                numArray14[numArray13[num3]] = 0;
            }
             
            this.UD1[num3] = numArray10[num3];
            num3++;
        }
        for (num3 = 1;num3 <= this.NDUAL;num3++)
        {
            numArray5[num3] = numArray9[num3];
        }
        goto Label_1568
    }

    private void dN33PF(double[] Z) throws Exception {
        int num2;
        int num1 = 1;
        while (num1 <= this.NDUAL)
        {
            double num3 = 0;
            for (num2 = 1;num2 <= this.NDUAL;num2++)
            {
                num3 += this.XJ[num1, num2] * this.NP[num2];
            }
            this.DDUAL[num1] = num3;
            num1++;
        }
        for (num1 = 1;num1 <= this.NDUAL;num1++)
        {
            Z[num1] = 0;
            for (num2 = this.IQ + 1;num2 <= this.NDUAL;num2++)
            {
                Z[num1] += this.XJ[num2, num1] * this.DDUAL[num2];
            }
        }
    }

    private void dN34PF(double[] RV) throws Exception {
        int num1 = this.NX + this.NRESM;
        int num2 = this.NRESM * 2;
        for (int num4 = this.IQ;num4 >= 1;num4--)
        {
            double num3 = 0;
            for (int num5 = num4 + 1;num5 <= this.IQ;num5++)
            {
                num3 += this.R[num5, num4] * RV[num5];
            }
            RV[num4] = (this.DDUAL[num4] - num3) / this.R[num4, num4];
        }
    }

    private void dN35PF(int[] AI, int L) throws Exception {
        int num5;
        int num6;
        double num7;
        double num8;
        double num9;
        double num10;
        double num11;
        double num12 = new double();
        double num13 = new double();
        double num14 = new double();
        int num1 = this.NX + this.NRESM;
        int num2 = this.NRESM * 2;
        int num3 = 1;
        int num4 = 1;
        while (num4 <= this.IQ)
        {
            if (AI[num4] == L)
            {
                num3 = num4;
                num4 = num3;
                while (num4 <= (this.IQ - 1))
                {
                    AI[num4] = AI[num4 + 1];
                    this.UD1[num4] = this.UD1[num4 + 1];
                    for (num5 = 1;num5 <= this.NDUAL;num5++)
                    {
                        this.R[num4, num5] = this.R[num4 + 1, num5];
                    }
                    num4++;
                }
                AI[this.IQ] = AI[this.IQ + 1];
                this.UD1[this.IQ] = this.UD1[this.IQ + 1];
                AI[this.IQ + 1] = 0;
                this.UD1[this.IQ + 1] = 0;
                num5 = 1;
                while (num5 <= this.IQ)
                {
                    this.R[this.IQ, num5] = 0;
                    num5++;
                }
                this.IQ--;
                if (this.IQ != 0)
                {
                    for (num5 = num3;num5 <= this.IQ;num5++)
                    {
                        num9 = this.R[num5, num5];
                        num10 = this.R[num5, num5 + 1];
                        num11 = this.dN27PF(num9,num10);
                        if (num11 != 0)
                        {
                            num12 = num9 / num11;
                            num13 = num10 / num11;
                            this.R[num5, num5 + 1] = 0;
                            if (num12 < 0)
                            {
                                this.R[num5, num5] = -num11;
                                num12 = -num12;
                                num13 = -num13;
                            }
                            else
                            {
                                this.R[num5, num5] = num11;
                            } 
                            num14 = num13 / (1 + num12);
                            num6 = num5 + 1;
                            while (num6 <= this.IQ)
                            {
                                num7 = this.R[num6, num5];
                                num8 = this.R[num6, num5 + 1];
                                this.R[num6, num5] = (num7 * num12) + (num8 * num13);
                                this.R[num6, num5 + 1] = (num14 * (num7 + this.R[num6, num5])) - num8;
                                num6++;
                            }
                            for (num6 = 1;num6 <= this.NDUAL;num6++)
                            {
                                num7 = this.XJ[num5, num6];
                                num8 = this.XJ[num5 + 1, num6];
                                this.XJ[num5, num6] = (num7 * num12) + (num8 * num13);
                                this.XJ[num5 + 1, num6] = (num14 * (this.XJ[num5, num6] + num7)) - num8;
                            }
                        }
                         
                    }
                }
                 
                this.RNORM = 1;
                this.RLOW = 1;
                if (this.IQ >= 1)
                {
                    this.RNORM = Math.Abs(this.R[1, 1]);
                    this.RLOW = Math.Abs(this.R[1, 1]);
                    num4 = 1;
                    while (num4 < this.IQ)
                    {
                        num4++;
                        this.RNORM = Math.Max(this.RNORM, Math.Abs(this.R[num4, num4]));
                        this.RLOW = Math.Min(this.RLOW, Math.Abs(this.R[num4, num4]));
                    }
                }
                 
                return ;
            }
             
            num4++;
        }
        num4 = num3;
        while (num4 <= (this.IQ - 1))
        {
            AI[num4] = AI[num4 + 1];
            this.UD1[num4] = this.UD1[num4 + 1];
            for (num5 = 1;num5 <= this.NDUAL;num5++)
            {
                this.R[num4, num5] = this.R[num4 + 1, num5];
            }
            num4++;
        }
        AI[this.IQ] = AI[this.IQ + 1];
        this.UD1[this.IQ] = this.UD1[this.IQ + 1];
        AI[this.IQ + 1] = 0;
        this.UD1[this.IQ + 1] = 0;
        num5 = 1;
        while (num5 <= this.IQ)
        {
            this.R[this.IQ, num5] = 0;
            num5++;
        }
        this.IQ--;
        if (this.IQ != 0)
        {
            for (num5 = num3;num5 <= this.IQ;num5++)
            {
                num9 = this.R[num5, num5];
                num10 = this.R[num5, num5 + 1];
                num11 = this.dN27PF(num9,num10);
                if (num11 != 0)
                {
                    num12 = num9 / num11;
                    num13 = num10 / num11;
                    this.R[num5, num5 + 1] = 0;
                    if (num12 < 0)
                    {
                        this.R[num5, num5] = -num11;
                        num12 = -num12;
                        num13 = -num13;
                    }
                    else
                    {
                        this.R[num5, num5] = num11;
                    } 
                    num14 = num13 / (1 + num12);
                    num6 = num5 + 1;
                    while (num6 <= this.IQ)
                    {
                        num7 = this.R[num6, num5];
                        num8 = this.R[num6, num5 + 1];
                        this.R[num6, num5] = (num7 * num12) + (num8 * num13);
                        this.R[num6, num5 + 1] = (num14 * (num7 + this.R[num6, num5])) - num8;
                        num6++;
                    }
                    for (num6 = 1;num6 <= this.NDUAL;num6++)
                    {
                        num7 = this.XJ[num5, num6];
                        num8 = this.XJ[num5 + 1, num6];
                        this.XJ[num5, num6] = (num7 * num12) + (num8 * num13);
                        this.XJ[num5 + 1, num6] = (num14 * (this.XJ[num5, num6] + num7)) - num8;
                    }
                }
                 
            }
        }
         
        this.RNORM = 1;
        this.RLOW = 1;
        if (this.IQ >= 1)
        {
            this.RNORM = Math.Abs(this.R[1, 1]);
            this.RLOW = Math.Abs(this.R[1, 1]);
            num4 = 1;
            while (num4 < this.IQ)
            {
                num4++;
                this.RNORM = Math.Max(this.RNORM, Math.Abs(this.R[num4, num4]));
                this.RLOW = Math.Min(this.RLOW, Math.Abs(this.R[num4, num4]));
            }
        }
         
    }

    private void dN36PF() throws Exception {
        int num1 = this.NX + this.NRESM;
        int num2 = this.NRESM * 2;
        for (int num4 = this.NDUAL;num4 >= (this.IQ + 2);num4--)
        {
            double num6 = this.DDUAL[num4 - 1];
            double num7 = this.DDUAL[num4];
            double num8 = this.dN27PF(num6,num7);
            if (num8 != 0)
            {
                this.DDUAL[num4] = 0;
                double num9 = num7 / num8;
                double num10 = num6 / num8;
                if (num10 < 0)
                {
                    num10 = -num10;
                    num9 = -num9;
                    this.DDUAL[num4 - 1] = -num8;
                }
                else
                {
                    this.DDUAL[num4 - 1] = num8;
                } 
                double num13 = num9 / (1 + num10);
                for (int num5 = 1;num5 <= this.NDUAL;num5++)
                {
                    double num11 = this.XJ[num4 - 1, num5];
                    double num12 = this.XJ[num4, num5];
                    this.XJ[num4 - 1, num5] = (num11 * num10) + (num12 * num9);
                    this.XJ[num4, num5] = (num13 * (num11 + this.XJ[num4 - 1, num5])) - num12;
                }
            }
             
        }
        this.IQ++;
        int num3 = 1;
        while (num3 <= this.IQ)
        {
            this.R[this.IQ, num3] = this.DDUAL[num3];
            num3++;
        }
        this.RNORM = 1;
        this.RLOW = 1;
        if (this.IQ >= 1)
        {
            this.RNORM = Math.Abs(this.R[1, 1]);
            this.RLOW = Math.Abs(this.R[1, 1]);
            num3 = 1;
            while (num3 < this.IQ)
            {
                num3++;
                this.RNORM = Math.Max(this.RNORM, Math.Abs(this.R[num3, num3]));
                this.RLOW = Math.Min(this.RLOW, Math.Abs(this.R[num3, num3]));
            }
        }
         
    }

    private void dN37PF(int NDIM, int N, double[][] A, int NDUALM, int NDUAL, double[][] XJ) throws Exception {
        int num4 = NDUAL - N;
        for (int num2 = N;num2 >= 1;num2--)
        {
            XJ[num2 + num4, num2 + num4] = 1 / A[num2, num2];
            for (int num3 = num2 - 1;num3 >= 1;num3--)
            {
                double num5 = 0;
                for (int num1 = num3 + 1;num1 <= num2;num1++)
                {
                    num5 += A[num1, num3] * XJ[num2 + num4, num1 + num4];
                }
                XJ[num2 + num4, num3 + num4] = -num5 / A[num3, num3];
            }
        }
    }

    private double dN38PF(double[] X, Provisdom.Optimization.DoNlp2.IFunction FCN, double[] XLB, double[] XUB) throws Exception {
        double num2 = 0;
        if (this.BLOC)
        {
            return num2;
        }
         
        for (int num1 = 0;num1 < this.N;num1++)
        {
            this.XTR[num1] = X[num1] * this.XSC[num1];
        }
        return this.Evaluate(FCN, 1, 0, this.XTR, XLB, XUB);
    }

    private void dN39PF(double[] X, double[] GRADF, Provisdom.Optimization.DoNlp2.IFunction FCN, double[] XLB, double[] XUB) throws Exception {
        int num2;
        double num1 = 45;
        double num18 = 0;
        double num9 = 0;
        double num10 = 0;
        double num11 = 0;
        double num12 = 0;
        double num13 = 0;
        double num14 = 0;
        double num15 = 0;
        double num17 = 0;
        double num16 = 0;
        if (this.GUNIT[0, 1] == 1)
        {
            for (num2 = 1;num2 <= this.N;num2++)
            {
                GRADF[num2] = 0;
            }
            GRADF[this.GUNIT[0, 2]] = this.GUNIT[0, 3] * this.XSC[this.GUNIT[0, 2] - 1];
            goto Label_0541
        }
         
        if (!this.BLOC)
        {
            num2 = 0;
            while (num2 < this.N)
            {
                this.XTR[num2] = this.XSC[num2] * X[num2];
                num2++;
            }
            if (!this.ANALYT)
            {
                if (this.DifFTYPE == 1)
                {
                    this.DELDif = Math.Min((double)(0.1 * Math.Sqrt(this.EPSFCN)), (double)0.01);
                    num18 = this.Evaluate(FCN, 1, 0, this.XTR, XLB, XUB);
                    for (num2 = 0;num2 < this.N;num2++)
                    {
                        num17 = this.XTR[num2];
                        num16 = Math.Min((this.DELDif * Math.Abs(num17)) + this.DELDif, this.TAUBND);
                        num16 = Math.Min(0.01, num16);
                        if (num17 >= 0)
                        {
                            this.XTR[num2] = num17 + num16;
                        }
                        else
                        {
                            this.XTR[num2] = num17 - num16;
                        } 
                        num9 = this.Evaluate(FCN, 1, 0, this.XTR, XLB, XUB);
                        GRADF[num2 + 1] = (num9 - num18) / (this.XTR[num2] - num17);
                        this.XTR[num2] = num17;
                    }
                    goto Label_051F
                }
                 
                if (this.DifFTYPE == 2)
                {
                    this.DELDif = Math.Min((double)(0.1 * Math.Pow(this.EPSFCN, 0.33333333333333331)), (double)0.01);
                    for (num2 = 0;num2 < this.N;num2++)
                    {
                        num17 = this.XTR[num2];
                        num16 = Math.Min((this.DELDif * Math.Abs(num17)) + this.DELDif, this.TAUBND);
                        num16 = Math.Min(0.01, num16);
                        this.XTR[num2] = num17 + num16;
                        num10 = this.Evaluate(FCN, 1, 0, this.XTR, XLB, XUB);
                        this.XTR[num2] = num17 - num16;
                        num11 = this.Evaluate(FCN, 1, 0, this.XTR, XLB, XUB);
                        GRADF[num2 + 1] = (num10 - num11) / (num16 + num16);
                        this.XTR[num2] = num17;
                    }
                    goto Label_051F
                }
                 
                this.DELDif = Math.Min((double)(0.1 * Math.Pow(this.EPSFCN, 0.14285714285714285)), (double)0.01);
                for (num2 = 0;num2 < this.N;num2++)
                {
                    num17 = this.XTR[num2];
                    num16 = Math.Min((double)((this.DELDif * Math.Abs(num17)) + this.DELDif), (double)(this.TAUBND / 4));
                    num16 = Math.Min(0.01, num16);
                    this.XTR[num2] = num17 - num16;
                    num10 = this.Evaluate(FCN, 1, 0, this.XTR, XLB, XUB);
                    this.XTR[num2] = num17 + num16;
                    num11 = this.Evaluate(FCN, 1, 0, this.XTR, XLB, XUB);
                    num16 += num16;
                    double num3 = num16;
                    this.XTR[num2] = num17 - num16;
                    num12 = this.Evaluate(FCN, 1, 0, this.XTR, XLB, XUB);
                    this.XTR[num2] = num17 + num16;
                    num13 = this.Evaluate(FCN, 1, 0, this.XTR, XLB, XUB);
                    num16 += num16;
                    double num4 = num16;
                    this.XTR[num2] = num17 - num16;
                    num14 = this.Evaluate(FCN, 1, 0, this.XTR, XLB, XUB);
                    this.XTR[num2] = num17 + num16;
                    num15 = this.Evaluate(FCN, 1, 0, this.XTR, XLB, XUB);
                    this.XTR[num2] = num17;
                    double num5 = num16 + num16;
                    double num6 = (num11 - num10) / num3;
                    double num7 = (num13 - num12) / num4;
                    double num8 = (num15 - num14) / num5;
                    num8 = num7 - num8;
                    num7 = num6 - num7;
                    num8 = num7 - num8;
                    GRADF[num2 + 1] = (num6 + (0.4 * num7)) + (num8 / num1);
                }
                goto Label_051F
            }
             
            this.DN21PG((Provisdom.Optimization.DoNlp2.IGradient)this.F, 1, 0, this.XTR, GRADF);
            for (num2 = this.N;num2 >= 1;num2--)
            {
                GRADF[num2] = this.XSC[num2 - 1] * GRADF[num2 - 1];
            }
        }
         
        return ;
        Label_051F:num2 = 1;
        while (num2 <= this.N)
        {
            GRADF[num2] = this.XSC[num2 - 1] * GRADF[num2];
            num2++;
        }
        Label_0541:    ;
    }

    private void solve(Provisdom.Optimization.DoNlp2.IFunction FCN, int IPRINT, int MAXITN, double[] XLB, double[] XUB) throws Exception {
        int num1;
        int num4;
        int num5;
        double num12 = new double();
        double num16 = new double();
        double num17 = new double();
        int num24 = new int();
        double num26 = new double();
        double[] numArray1 = new double[this.NX + 1];
        double[] numArray2 = new double[this.NX + 1];
        double[] numArray3 = new double[this.NX + 1];
        double[] numArray4 = new double[this.NX + 1];
        int[] numArray5 = new int[(this.NRESM + 1) + 1];
        int[] numArray6 = new int[this.NRESM + 1];
        int num3 = 1;
        while (num3 <= this.N)
        {
            this.D[num3] = 0;
            this.E0[num3] = 0;
            num3++;
        }
        double num14 = 0;
        double num30 = 0;
        double num28 = 0;
        double num27 = 0;
        boolean flag2 = false;
        this.ITSTEP = 0;
        this.ALIST[0] = this.NH;
        numArray5[0] = 0;
        this.VIOLIS[0] = 0;
        this.UPSI = 0;
        this.PSI = 0;
        this.PSI0 = 0;
        this.SIG0 = 0;
        this.e0NORM = 1;
        double num18 = 1;
        this.DNORM = 1;
        this.DEL = this.DEL0;
        int num8 = 0;
        int num9 = 0;
        int num23 = 0;
        int num7 = 0;
        int num6 = 0;
        int num10 = 0;
        this.MATSC = 1;
        this.TAUQP = 1;
        boolean flag1 = false;
        this.IDENT = false;
        boolean flag3 = false;
        if ((this.N > 100) || (this.NRES > 100))
        {
            this.TE3 = false;
        }
         
        num3 = 1;
        while (num3 <= this.N)
        {
            this.PERM[num3] = num3;
            this.PERM1[num3] = num3;
            num3++;
        }
        num3 = 1;
        while (num3 <= this.NH)
        {
            this.BINe0[num3] = 1;
            this.BIND[num3] = 1;
            this.ALIST[num3] = num3;
            num3++;
        }
        if (this.ANALYT)
        {
            num26 = Math.Min(this.EPSX, Math.Sqrt(this.EPSMAC));
        }
        else
        {
            num26 = this.EPSDif;
            if (this.EPSX < (this.EPSDif * this.EPSDif))
            {
                this.EPSX = this.EPSDif * this.EPSDif;
            }
             
        } 
        num26 = Math.Max(this.EPSMAC * 1000, Math.Min(0.001, num26));
        num3 = 1;
        while (num3 <= this.NRES)
        {
            if (num3 <= this.NH)
            {
                this.CFUERR[num3] = false;
                this.RES[num3] = this.DN40PF(num3, this.X, FCN, XLB, XUB);
                if (this.CFUERR[num3])
                {
                    this.OPTITE = -10;
                    this.throwError(this.OPTITE);
                    return ;
                }
                 
                num14 = Math.Abs(this.RES[num3]);
                if (!this.GCONST[num3])
                {
                    this.DN41PF(num3, this.X, numArray2, FCN, XLB, XUB);
                    this.VAL[num3] = true;
                    for (num4 = 1;num4 <= this.N;num4++)
                    {
                        this.GRES[num3, num4] = numArray2[num4];
                    }
                }
                 
            }
            else
            {
                this.BIND[num3] = 0;
                this.CFUERR[num3] = false;
                this.RES[num3] = this.DN42PF(num3 - this.NH, this.X, FCN, XLB, XUB);
                if (this.CFUERR[num3])
                {
                    this.OPTITE = -10;
                    this.throwError(this.OPTITE);
                    return ;
                }
                 
                num14 = -Math.Min(0, this.RES[num3]);
                if (this.RES[num3] <= this.DELMIN)
                {
                    this.BIND[num3] = 1;
                    this.ALIST[0]++;
                    this.ALIST[this.ALIST[0]] = num3;
                    if (!this.GCONST[num3])
                    {
                        this.VAL[num3] = true;
                        this.DN43PF(num3 - this.NH, this.X, numArray2, FCN, XLB, XUB);
                        for (num4 = 1;num4 <= this.N;num4++)
                        {
                            this.GRES[num3, num4] = numArray2[num4];
                        }
                    }
                     
                }
                 
            } 
            this.UPSI += num14;
            this.PSI += num14 * this.W[num3];
            if (this.VAL[num3] & !this.GCONST[num3])
            {
                this.GRESN[num3] = Math.Max(1, this.DN31PF(1, this.N, numArray2));
            }
             
            num3++;
        }
        if (this.UPSI >= this.TAU0)
        {
            this.SCF = 0;
            this.PHASE = -1;
        }
        else
        {
            this.FFUERR = false;
            this.FX = this.DN38PF(this.X, FCN, XLB, XUB);
            if (this.FFUERR)
            {
                this.OPTITE = -9;
                this.throwError(this.OPTITE);
                return ;
            }
             
            if (!this.VAL[0])
            {
                this.DN39PF(this.X, this.GRADF, FCN, XLB, XUB);
                this.VAL[0] = true;
            }
             
            this.SCF = 1;
            this.PHASE = 0;
            this.FXST = this.FX;
            this.PSIST = this.PSI;
            this.UPSIST = this.UPSI;
            for (num4 = 1;num4 <= this.NRES;num4++)
            {
                this.RESST[num4] = this.RES[num4];
            }
            this.ETA = 0;
        } 
        Label_05AA:if (!this.IDENT)
        {
            this.ITSTEP++;
            if (this.ITSTEP > MAXITN)
            {
                this.OPTITE = -3;
                this.ITSTEP = MAXITN;
                this.B2N = this.ACCINF[8, this.ITSTEP];
                this.B2N0 = this.ACCINF[7, this.ITSTEP];
                this.throwError(this.OPTITE);
                return ;
            }
             
            flag2 = false;
            this.QPTERM = 0;
            num27 = this.DEL;
            this.DEL = 0;
            this.B2N0 = -1;
            this.B2N = -1;
            this.SINGUL = false;
            flag1 = false;
            num3 = 1;
            while (num3 <= this.N)
            {
                flag1 = flag1 || (this.PERM[num3] != this.PERM1[num3]);
                this.PERM[num3] = this.PERM1[num3];
                num3++;
            }
            for (num3 = 1;num3 <= this.NRES;num3++)
            {
                this.DIAG[num3] = 0;
            }
        }
         
        this.NR = this.ALIST[0];
        int num25 = this.NR;
        num3 = 1;
        while (num3 <= this.NRES)
        {
            numArray6[num3] = this.BIND[num3];
            num3++;
        }
        num4 = 1;
        while (num4 <= 0x20)
        {
            this.ACCINF[num4, this.ITSTEP] = 0;
            num4++;
        }
        this.GFN = this.dN31PF(1,this.N,this.GRADF);
        if (((((this.NRES > 0) & (this.PHASE >= 0)) & !this.IDENT) & (this.ITSTEP > 1)) & (((this.ACCINF[10, this.ITSTEP - 1] == -1) & (this.SCF0 == 1)) || (this.ACCINF[10, this.ITSTEP - 1] == 1)))
        {
            num14 = 0;
            num3 = 1;
            while (num3 <= this.NRES)
            {
                if (this.GUNIT[num3, 1] != 1)
                {
                    num14 = Math.Max(num14, this.GRESN[num3]);
                }
                 
                num3++;
            }
            num17 = num14 / Math.Max(1 / this.SCFMAX, this.GFN);
            if (num17 < (1 / this.SCFMAX))
            {
                num17 = 1 / this.SCFMAX;
            }
             
            if (num17 > this.SCFMAX)
            {
                num17 = this.SCFMAX;
            }
             
            if (((((this.FXST - this.FX) * num17) + ((num17 / this.SCF) * (this.PSIST - this.PSI))) >= (((num17 / this.SCF) * this.ETA) * this.CLOW)) && ((this.LASTCH <= (this.ITSTEP - 4)) & ((num17 < (0.1 * this.SCF)) || (num17 > (10 * this.SCF)))))
            {
                this.CLOW++;
                num14 = num17 / this.SCF;
                this.PSI *= num14;
                this.PSIST *= num14;
                num3 = 1;
                while (num3 <= this.NRES)
                {
                    this.U[num3] *= num14;
                    num3++;
                }
                num18 *= num14;
                this.SCF = num17;
                this.LASTCH = this.ITSTEP;
                num14 = Math.Sqrt(num14);
                for (num3 = 1;num3 <= this.N;num3++)
                {
                    this.DIAG0[num3] = num14 * this.DIAG0[num3];
                    for (num4 = 1;num4 <= this.N;num4++)
                    {
                        this.A[num3, num4] *= num14;
                    }
                }
                this.MATSC *= num14;
            }
             
        }
         
        this.ACCINF[1, this.ITSTEP] = this.ITSTEP;
        this.ACCINF[2, this.ITSTEP] = this.FX;
        this.ACCINF[3, this.ITSTEP] = this.SCF;
        this.ACCINF[4, this.ITSTEP] = this.PSI;
        this.ACCINF[5, this.ITSTEP] = this.UPSI;
        if (!this.SILENT)
        {
            this.dN44PF(1);
        }
         
        if (this.NR >= 1)
        {
            this.dN23PF(1,this.NR);
        }
        else
        {
            this.RANK = 0;
        } 
        num14 = this.DN30PF(this.A, this.GRADF, numArray2, this.N);
        num3 = 1;
        while (num3 <= this.N)
        {
            this.QGF[num3] = numArray2[this.PERM[num3]];
            num3++;
        }
        this.DN24PF(1, 0, 1, this.RANK, this.N, this.QR, this.BETAQ, this.QGF, numArray4);
        num3 = 1;
        while (num3 <= this.N)
        {
            this.QGF[num3] = numArray4[num3];
            num3++;
        }
        this.B2N0 = this.dN31PF(this.RANK + 1,this.N,this.QGF);
        double num13 = 0;
        num3 = 1;
        while (num3 <= this.NRES)
        {
            if (num3 <= this.NH)
            {
                num13 += Math.Abs(this.RES[num3]) / this.GRESN[num3];
            }
            else
            {
                num13 -= Math.Min(0, this.RES[num3]) / this.GRESN[num3];
            } 
            num3++;
        }
        if (((this.ITSTEP > 1) && (this.ACCINF[8, this.ITSTEP - 1] >= 0)) && (!flag3 && (this.ACCINF[0x12, this.ITSTEP - 1] >= 0)))
        {
            flag3 = true;
            this.ETA = ((((this.ACCINF[8, this.ITSTEP - 1] / Math.Max(1, this.GFN)) + num13) + Math.Min(1, num30)) + Math.Min(1, Math.Abs(num28))) / ((double)Math.Min(30 * this.N, this.ITERMA));
            this.LEVEL = this.ETA;
        }
         
        if (this.ITSTEP > 1)
        {
            num12 = this.DELMIN;
        }
        else
        {
            num12 = this.DEL01;
        } 
        num14 = ((this.SCF * (this.FX0 - this.FX)) + this.PSI0) - this.PSI;
        if (((num14 > 0) & (this.SCF != 0)) & (this.ITSTEP > 1))
        {
            num12 = Math.Max(num12, Math.Exp(0.48999999999999994 * Math.Log(num14)));
        }
         
        if (this.SCF == 0)
        {
            num12 = Math.Min(this.DEL0 * 0.0001, Math.Max(num12, this.UPSI * 0.01));
        }
         
        double num11 = this.DELMIN;
        num3 = 1;
        while (num3 <= this.VIOLIS[0])
        {
            num4 = this.VIOLIS[num3];
            num11 = Math.Max(num11, ((this.RES[num4 - 1] / this.GRESN[num4 - 1]) / this.DELFAC[num4 - 1]) * 1.1);
            num3++;
        }
        this.DEL = Math.Min(this.DEL0, Math.Max(Math.Min(num11, 5 * num12), num12));
        if (this.VIOLIS[0] == 0)
        {
            this.DEL = Math.Min(this.DEL, this.DEL01);
        }
         
        if ((this.PHASE == 2) && (this.VIOLIS[0] == 0))
        {
            for (num3 = this.NH + 1;num3 <= this.NRES;num3++)
            {
                if (this.BINe0[num3] == 1)
                {
                    this.DEL = Math.Min(this.DEL01, Math.Max(this.DEL, Math.Abs(this.RES[num3]) / this.GRESN[num3]));
                }
                 
            }
        }
         
        num14 = this.DEL;
        num3 = 1;
        while (num3 <= numArray5[0])
        {
            num4 = numArray5[num3];
            num16 = ((this.RES[num4] / this.GRESN[num4]) * 0.99) / this.DELFAC[num4];
            if (num16 >= (this.DEL * 0.01))
            {
                num14 = Math.Min(num14, num16);
            }
             
            num3++;
        }
        this.DEL = num14;
        if (((this.ITSTEP > 1) && !this.IDENT) && (this.SCF != 0))
        {
            num14 = 0;
            num3 = this.NH + 1;
            while (num3 <= this.NRES)
            {
                num14 += (Math.Max((double)0, (double)((this.RES[num3] / this.GRESN[num3]) - this.DELMIN)) * Math.Max((double)0, (double)(this.U[num3] - this.SMALLW))) / this.GRESN[num3];
                num3++;
            }
            if (num14 > 0)
            {
                for (num3 = this.NH + 1;num3 <= this.NRES;num3++)
                {
                    if ((this.U[num3] > this.SMALLW) && ((this.RES[num3] / this.GRESN[num3]) > this.DELMIN))
                    {
                        this.DEL = Math.Max(this.DELMIN, Math.Min(this.DEL, ((this.RES[num3] / this.GRESN[num3]) * 0.99) / this.DELFAC[num3]));
                    }
                     
                }
            }
             
        }
         
        if ((this.ITSTEP > 1) && (this.ACCINF[30, this.ITSTEP - 1] < 0))
        {
            this.DEL = Math.Min(10 * num27, this.DEL0);
        }
         
        num3 = this.NH + 1;
        while (num3 <= this.NRES)
        {
            num14 = this.RES[num3] / this.GRESN[num3];
            if ((this.BIND[num3] == 0) && (num14 <= (this.DEL * this.DELFAC[num3])))
            {
                this.BIND[num3] = 1;
                this.ALIST[0]++;
                this.ALIST[this.ALIST[0]] = num3;
                if (!this.VAL[num3])
                {
                    this.VAL[num3] = true;
                    this.DN43PF(num3 - this.NH, this.X, numArray2, FCN, XLB, XUB);
                    for (num4 = 1;num4 <= this.N;num4++)
                    {
                        this.GRES[num3, num4] = numArray2[num4];
                    }
                    this.GRESN[num3] = Math.Max(1, this.DN31PF(1, this.N, numArray2));
                }
                 
            }
             
            num3++;
        }
        int num21 = this.RANK;
        int num22 = this.NR;
        this.NR = this.ALIST[0];
        this.dN23PF(num22 + 1,this.NR);
        this.DN24PF(1, 0, num21 + 1, this.RANK, this.N, this.QR, this.BETAQ, this.QGF, numArray4);
        num3 = 1;
        while (num3 <= this.N)
        {
            this.QGF[num3] = numArray4[num3];
            num3++;
        }
        num3 = 1;
        while (num3 <= this.N)
        {
            numArray2[num3] = -this.QGF[num3] * this.SCF;
            num3++;
        }
        num3 = 1;
        while (num3 <= this.NRES)
        {
            this.YU[num3] = 0;
            num3++;
        }
        this.DN25PF(1, this.RANK, numArray2, this.YU);
        double num15 = 0;
        num18 = 0;
        num3 = 1;
        while (num3 <= this.NRES)
        {
            this.U[num3] = 0;
            num3++;
        }
        int num20 = 0;
        num28 = 0;
        num3 = 1;
        while (num3 <= this.RANK)
        {
            num18 = Math.Max(num18, Math.Abs(this.YU[num3]));
            num5 = this.ALIST[this.COLNO[num3]];
            this.U[num5] = -this.YU[num3];
            if ((num5 > this.NH) && ((-this.YU[num3] / this.GRESN[num5]) < num28))
            {
                num20 = num5;
                num28 = -this.YU[num3] / this.GRESN[num5];
            }
             
            num3++;
        }
        if (this.SCF == 0)
        {
            this.B2N = -1;
        }
        else
        {
            for (num3 = 1;num3 <= this.N;num3++)
            {
                numArray3[num3] = this.SCF * this.GRADF[num3];
                for (num4 = 1;num4 <= this.NRES;num4++)
                {
                    numArray3[num3] -= this.GRES[num4, num3] * this.U[num4];
                }
            }
            this.B2N = this.DN31PF(1, this.N, numArray3) / this.SCF;
        } 
        if (!this.SILENT)
        {
            this.dN44PF(3);
        }
         
        double num19 = this.DEL;
        if (this.B2N >= 0)
        {
            if ((((Math.Abs(this.B2N) / (this.GFN + 1)) + Math.Max(0, -num28)) + num13) <= 0)
            {
                num19 = Math.Max(this.DEL, 0.1 * Math.Min(this.DEL0, Math.Exp(-9.7999999999999989)));
            }
            else
            {
                num19 = Math.Max(this.DEL, 0.1 * Math.Min(this.DEL0, Math.Exp(0.7 * Math.Log(((Math.Abs(this.B2N) / (this.GFN + 1)) + Math.Max(0, -num28)) + num13))));
            } 
        }
         
        num3 = 1;
        while (num3 <= numArray5[0])
        {
            num4 = numArray5[num3];
            num16 = ((this.RES[num4] / this.GRESN[num4]) * 0.99) / this.DELFAC[num4];
            if (num16 >= (num19 * 0.01))
            {
                num19 = Math.Max(this.DELMIN, Math.Min(num19, num16));
            }
             
            num3++;
        }
        num30 = 0;
        num3 = this.NH + 1;
        while (num3 <= this.NRES)
        {
            num30 += (Math.Max((double)0, (double)((this.RES[num3] / this.GRESN[num3]) - this.DELMIN)) * Math.Max((double)0, (double)(this.U[num3] - this.SMALLW))) / this.GRESN[num3];
            num3++;
        }
        if ((((this.UPSI <= this.DELMIN) && (this.B2N <= (this.EPSX * (this.GFN + 1)))) && ((this.B2N != -1) && (num28 >= -this.SMALLW))) && (num30 <= ((this.DELMIN * this.SMALLW) * this.NRES)))
        {
            this.OPTITE = 0;
            return ;
        }
         
        int num2 = this.ALIST[0];
        num3 = 1;
        while (num3 <= this.NRES)
        {
            num14 = (this.RES[num3] / this.GRESN[num3]) / this.DELFAC[num3];
            if (((num14 > this.DEL) && (num14 <= num19)) && (this.BIND[num3] == 0))
            {
                this.BIND[num3] = 1;
                this.ALIST[0]++;
                this.ALIST[this.ALIST[0]] = num3;
                if (!this.VAL[num3])
                {
                    this.VAL[num3] = true;
                    this.DN43PF(num3 - this.NH, this.X, numArray2, FCN, XLB, XUB);
                    for (num4 = 1;num4 <= this.N;num4++)
                    {
                        this.GRES[num3, num4] = numArray2[num4];
                    }
                    this.GRESN[num3] = Math.Max(1, this.DN31PF(1, this.N, numArray2));
                }
                 
            }
             
            num3++;
        }
        this.DEL = num19;
        this.ACCINF[6, this.ITSTEP] = this.DEL;
        this.ACCINF[7, this.ITSTEP] = this.B2N0;
        this.ACCINF[9, this.ITSTEP] = this.ALIST[0];
        this.ACCINF[10, this.ITSTEP] = -1;
        this.NR = this.ALIST[0];
        if (num2 != this.NR)
        {
            num21 = this.RANK;
            this.dN23PF(num2 + 1,this.NR);
            this.DN24PF(1, 0, num21 + 1, this.RANK, this.N, this.QR, this.BETAQ, this.QGF, numArray4);
            for (num3 = 1;num3 <= this.N;num3++)
            {
                this.QGF[num3] = numArray4[num3];
            }
        }
         
        if (!this.SILENT)
        {
            this.dN44PF(2);
        }
         
        if (this.RANK == this.NR)
        {
            num3 = 1;
            while (num3 <= this.N)
            {
                numArray2[num3] = -this.QGF[num3] * this.SCF;
                num3++;
            }
            num3 = 1;
            while (num3 <= this.NRES)
            {
                this.YU[num3] = 0;
                num3++;
            }
            this.DN25PF(1, this.RANK, numArray2, this.YU);
            num15 = 0;
            num18 = 0;
            num3 = 1;
            while (num3 <= this.NRES)
            {
                this.U[num3] = 0;
                num3++;
            }
            num20 = 0;
            num28 = 0;
            num3 = 1;
            while (num3 <= this.RANK)
            {
                num18 = Math.Max(num18, Math.Abs(this.YU[num3]));
                num5 = this.ALIST[this.COLNO[num3]];
                this.U[num5] = -this.YU[num3];
                if (num5 > this.NH)
                {
                    num15 = Math.Min(num15, -this.YU[num3]);
                    if ((-this.YU[num3] / this.GRESN[num5]) < num28)
                    {
                        num20 = num5;
                        num28 = -this.YU[num3] / this.GRESN[num5];
                    }
                     
                }
                 
                num3++;
            }
            if (this.SCF != 0)
            {
                for (num3 = 1;num3 <= this.N;num3++)
                {
                    numArray3[num3] = this.SCF * this.GRADF[num3];
                    for (num4 = 1;num4 <= this.NRES;num4++)
                    {
                        numArray3[num3] -= this.GRES[num4, num3] * this.U[num4];
                    }
                }
                this.B2N = this.DN31PF(1, this.N, numArray3) / this.SCF;
            }
             
            this.ACCINF[8, this.ITSTEP] = this.B2N;
            this.ACCINF[11, this.ITSTEP] = num15;
            this.dN46PF();
            if (!this.SILENT)
            {
                this.dN44PF(4);
            }
             
            numArray5[0] = 0;
            if (((this.PHASE >= 0) && (this.B2N != -1)) && (Math.Abs(num28) >= Math.Max(this.SMALLW, (Math.Abs(this.B2N) / (this.GFN + 1)) * this.C1D)))
            {
                for (num3 = this.NH + 1;num3 <= this.NR;num3++)
                {
                    num5 = this.ALIST[this.COLNO[num3]];
                    if ((-this.YU[num3] / this.GRESN[num5]) <= -this.SMALLW)
                    {
                        numArray5[0]++;
                        numArray5[numArray5[0]] = num5;
                    }
                     
                }
            }
             
            this.EQRES = true;
            num3 = 1;
            while (num3 <= this.NRES)
            {
                this.EQRES = this.EQRES && (this.BIND[num3] == this.BINe0[num3]);
                num3++;
            }
            if (this.NR > 1)
            {
                num14 = 0;
                num16 = 1;
                for (num3 = 1;num3 <= this.NR;num3++)
                {
                    num14 = Math.Max(num14, Math.Abs(this.DIAG[num3]));
                    num16 = Math.Min(num16, Math.Abs(this.DIAG[num3]));
                }
                this.ACCINF[13, this.ITSTEP] = num14 / num16;
            }
            else if (this.NR == 1)
            {
                this.ACCINF[13, this.ITSTEP] = 1;
            }
            else
            {
                this.ACCINF[13, this.ITSTEP] = -1;
            }  
            num14 = Math.Abs(this.A[1, 1]);
            num16 = Math.Abs(this.A[1, 1]);
            num3 = 2;
            while (num3 <= this.N)
            {
                num14 = Math.Max(num14, Math.Abs(this.A[num3, num3]));
                num16 = Math.Min(num16, Math.Abs(this.A[num3, num3]));
                num3++;
            }
            this.ACCINF[14, this.ITSTEP] = (num14 / num16) * (num14 / num16);
            if (!this.SILENT)
            {
                this.dN44PF(5);
            }
             
            num30 = 0;
            num3 = this.NH + 1;
            while (num3 <= this.NRES)
            {
                num30 += (Math.Max((double)0, (double)((this.RES[num3] / this.GRESN[num3]) - this.DELMIN)) * Math.Max((double)0, (double)(this.U[num3] - this.SMALLW))) / this.GRESN[num3];
                num3++;
            }
            if ((((num15 >= -this.SMALLW) && (num30 <= ((this.DELMIN * this.SMALLW) * this.NRES))) && ((this.UPSI <= (this.NRES * this.DELMIN)) && (this.UPSI0 <= (this.NRES * this.DELMIN)))) && (((Math.Abs((double)(this.FX - this.FX0)) <= (num26 * (Math.Abs(this.FX) + 1))) && (this.B2N != -1)) && (this.B2N <= ((100 * this.EPSX) * (this.GFN + 1)))))
            {
                num23++;
            }
            else
            {
                num23 = 0;
            } 
            if (((this.PHASE >= 0) && ((this.ACCINF[14, this.ITSTEP] > 1000) || !this.ANALYT)) && (num23 > this.N))
            {
                this.OPTITE = 4;
                this.throwError(this.OPTITE);
                return ;
            }
             
            this.SCF0 = 1;
            if ((this.PHASE >= 0) && (this.UPSI > (this.TAU0 * 0.5)))
            {
                this.SCF0 = Math.Max((double)(1 / this.SCFMAX), (double)(((((2 * (this.TAU0 - this.UPSI)) / this.TAU0) * this.UPSI) * 0.1) / Math.Max(1, this.GFN))) / this.SCF;
            }
             
            this.ACCINF[15, this.ITSTEP] = this.SCF0;
            num3 = this.NR + 1;
            while (num3 <= this.N)
            {
                numArray1[num3] = numArray2[num3] * this.SCF0;
                num3++;
            }
            double num29 = 1;
            if (((-num15 * this.C1D) > (this.B2N + this.UPSI)) && (this.B2N != -1))
            {
                num29 = this.C1D;
            }
             
            if (this.UPSI > (this.TAU0 * 0.5))
            {
                num29 = 0;
            }
             
            num3 = 1;
            while (num3 <= this.NR)
            {
                num5 = this.ALIST[this.COLNO[num3]];
                num14 = this.RES[num5];
                if (((num5 > this.NH) && (-this.YU[num3] < 0)) && (num14 > 0))
                {
                    num14 = -num14;
                }
                 
                if ((num5 > this.NH) && (-this.YU[num3] < 0))
                {
                    num14 -= this.YU[num3] * num29;
                }
                 
                numArray3[num3] = -num14;
                num3++;
            }
            this.DN26PF(1, this.NR, numArray3, numArray1);
            this.DN24PF(-1, 0, 1, this.NR, this.N, this.QR, this.BETAQ, numArray1, numArray3);
            num3 = 1;
            while (num3 <= this.N)
            {
                numArray1[this.PERM[num3]] = numArray3[num3];
                num3++;
            }
            num14 = this.DN29PF(this.A, numArray1, this.D, this.N);
            num24 = this.CLOW;
            if (this.PHASE >= 0)
            {
                this.dN8LPF();
            }
             
            if (this.CLOW > num24)
            {
                num14 = this.W[0];
                for (num3 = 1;num3 <= this.NRES;num3++)
                {
                    num14 = Math.Max(num14, this.W[num3]);
                }
                this.TAUQP = Math.Max(1, Math.Min(this.TAUQP, num14));
            }
             
            if (num28 < -this.SMALLW)
            {
                this.PHASE = Math.Min(1, this.PHASE);
            }
             
            if (!this.EQRES)
            {
                this.PHASE = Math.Min(0, this.PHASE);
            }
             
            if (this.EQRES && (this.UPSI < this.TAU0))
            {
                this.PHASE = Math.Max(1, this.PHASE);
            }
             
            this.dN12PF();
            this.dN11PF();
            if ((((this.DNORM <= (this.EPSX * (this.XNORM + this.EPSX))) && (this.UPSI <= this.DELMIN)) && ((this.B2N != -1) && (num28 >= -this.SMALLW))) && (this.B2N <= (this.EPSX * (this.GFN + 1))))
            {
                this.OPTITE = 1;
                return ;
            }
             
            goto Label_434F
        }
         
        this.SINGUL = true;
        this.PHASE = Math.Min(this.PHASE, 0);
        this.ACCINF[10, this.ITSTEP] = 1;
        this.SCF0 = 1;
        if ((this.PHASE >= 0) && (this.UPSI > (this.TAU0 * 0.5)))
        {
            num17 = Math.Max(1 / this.SCFMAX, Math.Min(this.SCFMAX, ((((2 * (this.TAU0 - this.UPSI)) / this.TAU0) * this.UPSI) * this.TAU) / Math.Max(1, this.GFN)));
            if ((((((this.FXST - this.FX) * num17) + ((num17 / this.SCF) * (this.PSIST - this.PSI))) >= (((num17 / this.SCF) * this.ETA) * this.CLOW)) && (this.LASTCH <= (this.ITSTEP - 4))) && ((num17 < (0.1 * this.SCF)) || (num17 > (10 * this.SCF))))
            {
                this.CLOW++;
                num14 = num17 / this.SCF;
                this.SCF0 = num14;
                this.PSI *= num14;
                this.PSIST *= num14;
                num3 = 1;
                while (num3 <= this.NRES)
                {
                    this.U[num3] *= num14;
                    num3++;
                }
                num18 *= num14;
                this.SCF = num17;
                this.LASTCH = this.ITSTEP;
                this.ACCINF[15, this.ITSTEP] = this.SCF;
                num14 = Math.Sqrt(num14);
                for (num3 = 1;num3 <= this.N;num3++)
                {
                    this.DIAG0[num3] = num14 * this.DIAG0[num3];
                    for (num4 = 1;num4 <= this.N;num4++)
                    {
                        this.A[num3, num4] *= num14;
                    }
                }
                this.MATSC *= num14;
            }
             
        }
         
        this.ACCINF[0x20, this.ITSTEP] = this.UPSI;
        if (!this.SILENT)
        {
            this.dN44PF(15);
        }
         
        this.ACCINF[13, this.ITSTEP] = -1;
        num14 = Math.Abs(this.A[1, 1]);
        num16 = num14;
        num3 = 2;
        while (num3 <= this.N)
        {
            num14 = Math.Max(num14, Math.Abs(this.A[num3, num3]));
            num16 = Math.Min(num16, Math.Abs(this.A[num3, num3]));
            num3++;
        }
        this.ACCINF[14, this.ITSTEP] = (num14 / num16) * (num14 / num16);
        if (!this.SILENT)
        {
            this.dN44PF(5);
        }
         
        num24 = this.CLOW;
        double num31 = this.TAUQP;
        num3 = 1;
        while (num3 <= this.NRES)
        {
            this.U[num3] = 0;
            num3++;
        }
        num3 = 1;
        while (num3 <= this.N)
        {
            this.DD[num3] = 0;
            num3++;
        }
        this.dN32PF();
        if (((this.DNORM == 0) && (this.QPTERM == 1)) && (this.OPTITE == 3))
        {
            return ;
        }
         
        if ((this.DNORM <= (this.EPSX * (Math.Min(this.XNORM, 1) + this.EPSX))) && (this.QPTERM < 0))
        {
            if ((this.UPSI >= (this.NRES * this.DELMIN)) && flag2)
            {
                this.OPTITE = -1;
                if (!this.SILENT)
                {
                    this.dN44PF(0x12);
                }
                 
                this.throwError(this.OPTITE);
                return ;
            }
             
            if (flag2)
            {
                this.OPTITE = this.QPTERM - 5;
                if (!this.SILENT)
                {
                    this.dN44PF(0x12);
                }
                 
                this.throwError(this.OPTITE);
                return ;
            }
             
            flag2 = true;
            num3 = 1;
            while (num3 <= this.NRES)
            {
                this.W[num3] = 1;
                num3++;
            }
            this.LASTCH = this.ITSTEP;
            this.dN10PF();
            this.IDENT = true;
            this.TAUQP = 1;
            this.ALIST[0] = num25;
            for (num3 = 1;num3 <= this.NRES;num3++)
            {
                this.BIND[num3] = numArray6[num3];
            }
            if (this.SCF == 0)
            {
                if (this.UPSI >= this.TAU0)
                {
                    this.SCF = 0;
                    this.PHASE = -1;
                }
                else
                {
                    this.FFUERR = false;
                    this.FX = this.DN38PF(this.X, FCN, XLB, XUB);
                    if (this.FFUERR)
                    {
                        this.OPTITE = -9;
                        this.throwError(this.OPTITE);
                        return ;
                    }
                     
                    if (!this.VAL[0])
                    {
                        this.DN39PF(this.X, this.GRADF, FCN, XLB, XUB);
                        this.VAL[0] = true;
                    }
                     
                    this.SCF = 1;
                    this.PHASE = 0;
                    this.FXST = this.FX;
                    this.PSIST = this.PSI;
                    this.UPSIST = this.UPSI;
                    for (num4 = 1;num4 <= this.NRES;num4++)
                    {
                        this.RESST[num4] = this.RES[num4];
                    }
                    this.ETA = 0;
                } 
            }
             
            goto Label_05AA
        }
         
        this.ACCINF[11, this.ITSTEP] = 0;
        numArray5[0] = 0;
        num15 = 0;
        this.dN46PF();
        if (((this.QPTERM >= 0) || (this.QPTERM == -3)) && (this.SCF != 0))
        {
            num18 = Math.Abs(this.U[1]);
            num3 = 2;
            while (num3 <= this.NRES)
            {
                num18 = Math.Max(num18, Math.Abs(this.U[num3]));
                num3++;
            }
            for (num3 = 1;num3 <= this.N;num3++)
            {
                numArray3[num3] = this.SCF * this.GRADF[num3];
                for (num4 = 1;num4 <= this.NRES;num4++)
                {
                    numArray3[num3] -= this.GRES[num4, num3] * this.U[num4];
                }
            }
            this.B2N = this.DN31PF(1, this.N, numArray3) / this.SCF;
        }
        else
        {
            this.B2N = -1;
        } 
        if (!this.SILENT)
        {
            this.dN44PF(2);
            this.dN44PF(0x10);
        }
         
        if (this.B2N == -1)
        {
            this.B2N = this.EPSMAC / this.TOLMAC;
        }
         
        if ((this.QPTERM >= 0) && (this.DNORM <= ((0.01 * this.EPSX) * (this.EPSX + Math.Min(1, this.XNORM)))))
        {
            if (this.UPSI <= (this.NRES * this.DELMIN))
            {
                this.OPTITE = 6;
            }
            else
            {
                this.OPTITE = -5;
            } 
            this.throwError(this.OPTITE);
            return ;
        }
         
        if (this.QPTERM < 0)
        {
        }
         
        if (this.CLOW > num24)
        {
            num14 = 1;
            for (num3 = 1;num3 <= this.NRES;num3++)
            {
                num14 = Math.Max(num14, this.W[num3]);
            }
            this.TAUQP = Math.Max(1, num14);
        }
         
        if (this.TAUQP > Math.Pow(this.TAUFAC, 3 * num31))
        {
            this.TAUQP = Math.Pow(this.TAUFAC, 3 * num31);
        }
         
        this.B2N0 = this.B2N;
        num15 = 0;
        num28 = 0;
        if (this.QPTERM >= 1)
        {
            num30 = 0;
            for (num3 = this.NH + 1;num3 <= this.NRES;num3++)
            {
                num30 += (Math.Max((double)0, (double)((this.RES[num3] / this.GRESN[num3]) - this.DELMIN)) * Math.Max((double)0, (double)(this.U[num3] - this.SMALLW))) / this.GRESN[num3];
            }
        }
        else
        {
            num30 = 1;
        } 
        Label_2088:this.ACCINF[0x10, this.ITSTEP] = this.XNORM;
        this.ACCINF[0x11, this.ITSTEP] = this.DNORM;
        this.ACCINF[0x12, this.ITSTEP] = this.PHASE;
        this.CFINCR = this.ICF;
        if (this.DIRDER >= 0)
        {
            this.STPTRM = -2;
            this.SIG = 0;
            goto Label_2525
        }
         
        if (-this.DIRDER <= ((this.EPSMAC * 100) * (((this.SCF * Math.Abs(this.FX)) + this.PSI) + 1)))
        {
            if (this.UPSI > (this.DELMIN * this.NRES))
            {
                this.OPTITE = -1;
                this.STPTRM = -1;
            }
            else
            {
                this.OPTITE = 2;
                this.STPTRM = 1;
            } 
            this.SIG = 0;
            goto Label_2525
        }
         
        if ((((this.PHASE >= 1) && (this.DNORM <= (this.SMALLD * (this.XNORM + this.SMALLD)))) && ((this.SCF0 == 1) && (num28 >= -this.SMALLW))) && !this.SINGUL)
        {
            this.PHASE = 2;
        }
         
        if ((this.PHASE == 2) && (this.DNORM > (this.XNORM + this.SMALLD)))
        {
            this.PHASE = 1;
        }
         
        this.dN13PF();
        num3 = 1;
        while (num3 <= this.N)
        {
            this.DD[num3] = 0;
            this.X1[num3 - 1] = this.X[num3 - 1] + this.D[num3];
            num3++;
        }
        boolean flag4 = false;
        num3 = 1;
        while (num3 <= this.N)
        {
            if ((this.LLOW[num3] && (this.X1[num3 - 1] < (this.UG[num3] - this.TAUBND))) || (this.LUP[num3] && (this.X1[num3 - 1] > (this.OG[num3] + this.TAUBND))))
            {
                flag4 = true;
            }
             
            num3++;
        }
        if (((this.PHASE == 2) && (this.DNORM > (this.XNORM * Math.Sqrt(this.EPSMAC)))) && (!this.SINGUL && !flag4))
        {
            num3 = 1;
            while (num3 <= this.ALIST[0])
            {
                numArray3[num3] = 0;
                if ((num3 <= this.NH) && !this.GCONST[this.ALIST[num3]])
                {
                    this.CFUERR[num3] = false;
                    numArray3[num3] = this.DN40PF(num3, this.X1, FCN, XLB, XUB);
                    if (!this.CFUERR[num3])
                    {
                        goto Label_23C8
                    }
                     
                    goto Label_24EB
                }
                 
                if (!this.GCONST[this.ALIST[num3]])
                {
                    this.CFUERR[this.ALIST[num3]] = false;
                    numArray3[num3] = this.DN42PF(this.ALIST[num3] - this.NH, this.X1, FCN, XLB, XUB);
                    if (this.CFUERR[this.ALIST[num3]])
                    {
                        goto Label_24EB
                    }
                     
                }
                 
                Label_23C8:numArray3[num3] = -numArray3[num3];
                num3++;
            }
            num3 = 1;
            while (num3 <= this.ALIST[0])
            {
                numArray2[num3] = numArray3[this.COLNO[num3]];
                num3++;
            }
            this.DN26PF(1, this.NR, numArray2, this.DD);
            num3 = this.NR + 1;
            while (num3 <= this.N)
            {
                this.DD[num3] = 0;
                num3++;
            }
            this.DN24PF(-1, 0, 1, this.NR, this.N, this.QR, this.BETAQ, this.DD, numArray3);
            num3 = 1;
            while (num3 <= this.N)
            {
                this.DD[this.PERM[num3]] = numArray3[num3];
                num3++;
            }
            num14 = this.dN29PF(this.A,this.DD,this.DD,this.N);
            if (Math.Sqrt(num14) > (0.5 * this.DNORM))
            {
                for (num3 = 1;num3 <= this.N;num3++)
                {
                    this.DD[num3] = 0;
                }
            }
             
        }
         
        Label_24EB:if (!this.SILENT)
        {
            this.dN44PF(7);
        }
         
        this.SIG = Math.Min(1, this.STMAXL);
        this.DN17PF(this.SIG, FCN, XLB, XUB);
        Label_2525:this.CFINCR = this.ICF - this.CFINCR;
        if (!this.SILENT)
        {
            this.dN44PF(10);
        }
         
        num14 = ((this.SCF * (this.FX0 - this.FX)) + this.PSI0) - this.PSI;
        if (Math.Abs(num14) <= (this.EPSPHI * ((this.SCF * Math.Abs(this.FX)) + this.PSI)))
        {
            num10++;
        }
        else
        {
            num10 = 0;
        } 
        if (num10 > this.NUMSM)
        {
            this.OPTITE = 7;
            this.throwError(this.OPTITE);
            return ;
        }
         
        if (this.SIG <= 0.05)
        {
            if (this.SIG0 <= 0.05)
            {
                num6++;
            }
             
        }
        else
        {
            num6 = 0;
        } 
        this.ACCINF[0x15, this.ITSTEP] = this.SIG;
        this.ACCINF[0x16, this.ITSTEP] = this.CFINCR;
        this.ACCINF[0x17, this.ITSTEP] = this.DIRDER;
        this.ACCINF[0x18, this.ITSTEP] = this.DSCAL;
        this.ACCINF[0x19, this.ITSTEP] = this.COSPHI;
        this.ACCINF[0x1a, this.ITSTEP] = this.VIOLIS[0];
        if (((this.SIG == 0) && (this.STPTRM == 1)) && (this.OPTITE == 2))
        {
            if (!this.SILENT)
            {
                this.dN44PF(0x11);
            }
             
            return ;
        }
         
        if ((((this.STPTRM == 1) && (this.SIG <= 0.0001)) && ((this.ACCINF[13, this.ITSTEP] > 10000) && !this.SINGUL)) && (this.NRES > 0))
        {
            if (this.ACCINF[14, this.ITSTEP] > 10000)
            {
                this.dN10PF();
            }
             
            this.IDENT = true;
            this.SINGUL = true;
        }
        else
        {
            if (this.STPTRM >= 0)
            {
                if ((((this.SINGUL && (this.ITSTEP > this.N)) && ((Math.Abs((double)(this.FX - this.FX0)) <= (num26 * (Math.Abs(this.FX) + 1))) && (this.PHASE >= 0))) && (((this.UPSI <= (this.NRES * this.DELMIN)) && (this.UPSI0 <= (this.NRES * this.DELMIN))) && ((num30 <= ((this.DELMIN * this.SMALLW) * this.NRES)) && (this.INFEAS <= this.UPSI)))) && !this.IDENT)
                {
                    this.OPTITE = 4;
                    this.throwError(this.OPTITE);
                    return ;
                }
                 
                if (((this.SINGUL && (this.UPSI <= (this.DELMIN * this.NRES))) && ((this.UPSI0 <= (this.DELMIN * this.NRES)) && (this.B2N != -1))) && (((this.B2N <= (((this.GFN + 1) * this.EPSX) * 100)) && (this.PHASE >= 0)) && ((num30 <= ((this.DELMIN * this.SMALLW) * this.NRES)) && (this.INFEAS <= this.UPSI))))
                {
                    this.OPTITE = 3;
                    return ;
                }
                 
                num5 = 0;
                num3 = 1;
                while (num3 <= this.N)
                {
                    if (Math.Abs(this.DifX[num3]) >= (this.EPSX * (Math.Abs(this.X[num3 - 1]) + 0.01)))
                    {
                        num5 = 1;
                    }
                     
                    num3++;
                }
                if (num5 == 0)
                {
                    num9++;
                }
                else
                {
                    num9 = 0;
                } 
                if ((num9 > this.NRESET) && this.SINGUL)
                {
                    this.OPTITE = 5;
                    this.throwError(this.OPTITE);
                    return ;
                }
                 
                this.XNORM = this.dN31PF(0,this.N - 1,this.X);
                this.IDENT = false;
                this.dN22PF(this.GPHI0);
                num3 = 0;
                while (num3 <= this.NRES)
                {
                    if (!this.GCONST[num3])
                    {
                        this.VAL[num3] = false;
                    }
                     
                    num3++;
                }
                if ((this.PHASE >= 0) && !this.GCONST[0])
                {
                    this.VAL[0] = true;
                    this.DN39PF(this.X, this.GRADF, FCN, XLB, XUB);
                }
                 
                num3 = 1;
                while (num3 <= this.ALIST[0])
                {
                    num1 = this.ALIST[num3];
                    if (!this.VAL[num1])
                    {
                        this.VAL[num1] = true;
                        if (num1 <= this.NH)
                        {
                            this.DN41PF(num1, this.X, numArray3, FCN, XLB, XUB);
                        }
                        else
                        {
                            this.DN43PF(num1 - this.NH, this.X, numArray3, FCN, XLB, XUB);
                        } 
                        for (num4 = 1;num4 <= this.N;num4++)
                        {
                            this.GRES[num1, num4] = numArray3[num4];
                        }
                        this.GRESN[num1] = Math.Max(1, this.DN31PF(1, this.N, numArray3));
                    }
                     
                    num3++;
                }
                this.dN22PF(this.GPHI1);
                num3 = 1;
                while (num3 <= this.N)
                {
                    numArray3[num3] = this.X[num3 - 1] - this.X0[num3 - 1];
                    numArray2[num3] = this.GPHI1[num3] - this.GPHI0[num3];
                    num3++;
                }
                num14 = Math.Sqrt(this.DN31PF(1, this.N, numArray2) / this.DN31PF(1, this.N, numArray3));
                if ((num14 != 0) && (this.PHASE >= 0))
                {
                    this.MATSC = Math.Max(1 / this.SCFMAX, Math.Min(this.SCFMAX, num14 / 2));
                }
                 
                this.ALIST[0] = 0;
                for (num3 = 1;num3 <= this.NRES;num3++)
                {
                    this.U0[num3] = this.U[num3];
                    this.BINe0[num3] = this.BIND[num3];
                    this.RES0[num3] = this.RES[num3];
                    if (num3 <= this.NH)
                    {
                        this.BIND[num3] = 1;
                    }
                    else
                    {
                        this.BIND[num3] = 0;
                        if ((this.RES[num3] / this.GRESN[num3]) <= this.DELMIN)
                        {
                            this.BIND[num3] = 1;
                        }
                         
                        if (!this.VAL[num3] && (this.BIND[num3] == 1))
                        {
                            this.VAL[num3] = true;
                            this.DN43PF(num3 - this.NH, this.X, numArray3, FCN, XLB, XUB);
                            for (num4 = 1;num4 <= this.N;num4++)
                            {
                                this.GRES[num3, num4] = numArray3[num4];
                            }
                            this.GRESN[num3] = Math.Max(1, this.DN31PF(1, this.N, numArray3));
                        }
                         
                    } 
                    if (this.BIND[num3] == 1)
                    {
                        this.ALIST[0]++;
                        this.ALIST[this.ALIST[0]] = num3;
                    }
                     
                }
                if (this.SCF != 0)
                {
                    if (((num7 > this.NRESET) || (num6 > this.NRESET)) || (num8 > this.NRESET))
                    {
                        num8 = 0;
                        num6 = 0;
                        num7 = 0;
                        this.dN10PF();
                    }
                    else
                    {
                        this.dN9LPF();
                        if (!this.SILENT)
                        {
                            this.dN44PF(14);
                        }
                         
                    } 
                }
                 
                if (this.ACCINF[0x1b, this.ITSTEP] == 1)
                {
                    if (((this.ITSTEP > 1) && (this.ACCINF[0x1d, this.ITSTEP - 1] != 0)) && (this.ACCINF[0x1d, this.ITSTEP] != 0))
                    {
                        num7++;
                    }
                    else
                    {
                        num7 = 0;
                    } 
                }
                 
                if (this.PHASE == -1)
                {
                    if (this.UPSI >= this.TAU0)
                    {
                        this.SCF = 0;
                        this.PHASE = -1;
                    }
                    else
                    {
                        this.FFUERR = false;
                        this.FX = this.DN38PF(this.X, FCN, XLB, XUB);
                        if (this.FFUERR)
                        {
                            this.OPTITE = -9;
                            this.throwError(this.OPTITE);
                            return ;
                        }
                         
                        if (!this.VAL[0])
                        {
                            this.DN39PF(this.X, this.GRADF, FCN, XLB, XUB);
                            this.VAL[0] = true;
                        }
                         
                        this.SCF = 1;
                        this.PHASE = 0;
                        this.FXST = this.FX;
                        this.PSIST = this.PSI;
                        this.UPSIST = this.UPSI;
                        for (num4 = 1;num4 <= this.NRES;num4++)
                        {
                            this.RESST[num4] = this.RES[num4];
                        }
                        this.ETA = 0;
                    } 
                }
                 
                goto Label_05AA
            }
             
            if (!this.IDENT)
            {
                this.IDENT = true;
                numArray5[0] = 0;
                this.VIOLIS[0] = 0;
                num8 = 0;
                num6 = 0;
                num7 = 0;
                this.dN10PF();
                this.ALIST[0] = num25;
                for (num3 = 1;num3 <= this.NRES;num3++)
                {
                    this.BIND[num3] = numArray6[num3];
                }
                if (this.UPSI >= this.TAU0)
                {
                    if (this.UPSI >= this.TAU0)
                    {
                        this.SCF = 0;
                        this.PHASE = -1;
                    }
                    else
                    {
                        this.FFUERR = false;
                        this.FX = this.DN38PF(this.X, FCN, XLB, XUB);
                        if (this.FFUERR)
                        {
                            this.OPTITE = -9;
                            this.throwError(this.OPTITE);
                            return ;
                        }
                         
                        if (!this.VAL[0])
                        {
                            this.DN39PF(this.X, this.GRADF, FCN, XLB, XUB);
                            this.VAL[0] = true;
                        }
                         
                        this.SCF = 1;
                        this.PHASE = 0;
                        this.FXST = this.FX;
                        this.PSIST = this.PSI;
                        this.UPSIST = this.UPSI;
                        for (num4 = 1;num4 <= this.NRES;num4++)
                        {
                            this.RESST[num4] = this.RES[num4];
                        }
                        this.ETA = 0;
                    } 
                }
                 
                goto Label_05AA
            }
             
            if ((!this.SINGUL && this.IDENT) && ((this.ACCINF[13, this.ITSTEP] > 10000) && (this.NRES > 0)))
            {
                this.SINGUL = true;
                this.IDENT = true;
            }
            else
            {
                if (this.STPTRM == -2)
                {
                    this.OPTITE = -4;
                    this.throwError(this.OPTITE);
                }
                else
                {
                    if ((this.SIG == 0) && (this.OPTITE == -1))
                    {
                        this.throwError(this.OPTITE);
                        return ;
                    }
                     
                    this.OPTITE = -2;
                    this.throwError(this.OPTITE);
                } 
                return ;
            } 
        } 
        this.SINGUL = true;
        this.PHASE = Math.Min(this.PHASE, 0);
        this.ACCINF[10, this.ITSTEP] = 1;
        this.SCF0 = 1;
        if ((this.PHASE >= 0) && (this.UPSI > (this.TAU0 * 0.5)))
        {
            num17 = Math.Max(1 / this.SCFMAX, Math.Min(this.SCFMAX, ((((2 * (this.TAU0 - this.UPSI)) / this.TAU0) * this.UPSI) * this.TAU) / Math.Max(1, this.GFN)));
            if ((((((this.FXST - this.FX) * num17) + ((num17 / this.SCF) * (this.PSIST - this.PSI))) >= (((num17 / this.SCF) * this.ETA) * this.CLOW)) && (this.LASTCH <= (this.ITSTEP - 4))) && ((num17 < (0.1 * this.SCF)) || (num17 > (10 * this.SCF))))
            {
                this.CLOW++;
                num14 = num17 / this.SCF;
                this.SCF0 = num14;
                this.PSI *= num14;
                this.PSIST *= num14;
                num3 = 1;
                while (num3 <= this.NRES)
                {
                    this.U[num3] *= num14;
                    num3++;
                }
                num18 *= num14;
                this.SCF = num17;
                this.LASTCH = this.ITSTEP;
                this.ACCINF[14, this.ITSTEP] = this.SCF;
                num14 = Math.Sqrt(num14);
                for (num3 = 1;num3 <= this.N;num3++)
                {
                    this.DIAG0[num3] = num14 * this.DIAG0[num3];
                    for (num4 = 1;num4 <= this.N;num4++)
                    {
                        this.A[num3, num4] *= num14;
                    }
                }
                this.MATSC *= num14;
            }
             
        }
         
        this.ACCINF[0x20, this.ITSTEP] = this.UPSI;
        if (!this.SILENT)
        {
            this.dN44PF(15);
        }
         
        this.ACCINF[13, this.ITSTEP] = -1;
        num14 = Math.Abs(this.A[1, 1]);
        num16 = num14;
        num3 = 2;
        while (num3 <= this.N)
        {
            num14 = Math.Max(num14, Math.Abs(this.A[num3, num3]));
            num16 = Math.Min(num16, Math.Abs(this.A[num3, num3]));
            num3++;
        }
        this.ACCINF[14, this.ITSTEP] = (num14 / num16) * (num14 / num16);
        if (!this.SILENT)
        {
            this.dN44PF(5);
        }
         
        num24 = this.CLOW;
        num31 = this.TAUQP;
        num3 = 1;
        while (num3 <= this.NRES)
        {
            this.U[num3] = 0;
            num3++;
        }
        num3 = 1;
        while (num3 <= this.N)
        {
            this.DD[num3] = 0;
            num3++;
        }
        this.dN32PF();
        if (((this.DNORM == 0) && (this.QPTERM == 1)) && (this.OPTITE == 3))
        {
            return ;
        }
         
        if ((this.DNORM <= (this.EPSX * (Math.Min(this.XNORM, 1) + this.EPSX))) && (this.QPTERM < 0))
        {
            if ((this.UPSI >= (this.NRES * this.DELMIN)) && flag2)
            {
                this.OPTITE = -1;
                if (!this.SILENT)
                {
                    this.dN44PF(0x12);
                }
                 
                this.throwError(this.OPTITE);
                return ;
            }
             
            if (flag2)
            {
                this.OPTITE = this.QPTERM - 5;
                if (!this.SILENT)
                {
                    this.dN44PF(0x12);
                }
                 
                this.throwError(this.OPTITE);
                return ;
            }
             
            flag2 = true;
            num3 = 1;
            while (num3 <= this.NRES)
            {
                this.W[num3] = 1;
                num3++;
            }
            this.LASTCH = this.ITSTEP;
            this.dN10PF();
            this.IDENT = true;
            this.TAUQP = 1;
            this.ALIST[0] = num25;
            for (num3 = 1;num3 <= this.NRES;num3++)
            {
                this.BIND[num3] = numArray6[num3];
            }
            if (this.SCF == 0)
            {
                if (this.UPSI >= this.TAU0)
                {
                    this.SCF = 0;
                    this.PHASE = -1;
                }
                else
                {
                    this.FFUERR = false;
                    this.FX = this.DN38PF(this.X, FCN, XLB, XUB);
                    if (this.FFUERR)
                    {
                        this.OPTITE = -9;
                        this.throwError(this.OPTITE);
                        return ;
                    }
                     
                    if (!this.VAL[0])
                    {
                        this.DN39PF(this.X, this.GRADF, FCN, XLB, XUB);
                        this.VAL[0] = true;
                    }
                     
                    this.SCF = 1;
                    this.PHASE = 0;
                    this.FXST = this.FX;
                    this.PSIST = this.PSI;
                    this.UPSIST = this.UPSI;
                    for (num4 = 1;num4 <= this.NRES;num4++)
                    {
                        this.RESST[num4] = this.RES[num4];
                    }
                    this.ETA = 0;
                } 
            }
             
            goto Label_05AA
        }
         
        this.ACCINF[11, this.ITSTEP] = 0;
        numArray5[0] = 0;
        num15 = 0;
        this.dN46PF();
        if (((this.QPTERM >= 0) || (this.QPTERM == -3)) && (this.SCF != 0))
        {
            num18 = Math.Abs(this.U[1]);
            num3 = 2;
            while (num3 <= this.NRES)
            {
                num18 = Math.Max(num18, Math.Abs(this.U[num3]));
                num3++;
            }
            for (num3 = 1;num3 <= this.N;num3++)
            {
                numArray3[num3] = this.SCF * this.GRADF[num3];
                for (num4 = 1;num4 <= this.NRES;num4++)
                {
                    numArray3[num3] -= this.GRES[num4, num3] * this.U[num4];
                }
            }
            this.B2N = this.DN31PF(1, this.N, numArray3) / this.SCF;
        }
        else
        {
            this.B2N = -1;
        } 
        if (!this.SILENT)
        {
            this.dN44PF(2);
            this.dN44PF(0x10);
        }
         
        if (this.B2N == -1)
        {
            this.B2N = this.EPSMAC / this.TOLMAC;
        }
         
        if ((this.QPTERM >= 0) && (this.DNORM <= ((0.01 * this.EPSX) * (this.EPSX + Math.Min(1, this.XNORM)))))
        {
            if (this.UPSI <= (this.NRES * this.DELMIN))
            {
                this.OPTITE = 6;
            }
            else
            {
                this.OPTITE = -5;
            } 
            this.throwError(this.OPTITE);
            return ;
        }
         
        if (this.QPTERM < 0)
        {
        }
         
        if (this.CLOW > num24)
        {
            num14 = 1;
            for (num3 = 1;num3 <= this.NRES;num3++)
            {
                num14 = Math.Max(num14, this.W[num3]);
            }
            this.TAUQP = Math.Max(1, num14);
        }
         
        if (this.TAUQP > Math.Pow(this.TAUFAC, 3 * num31))
        {
            this.TAUQP = Math.Pow(this.TAUFAC, 3 * num31);
        }
         
        this.B2N0 = this.B2N;
        num15 = 0;
        num28 = 0;
        if (this.QPTERM >= 1)
        {
            num30 = 0;
            for (num3 = this.NH + 1;num3 <= this.NRES;num3++)
            {
                num30 += (Math.Max((double)0, (double)((this.RES[num3] / this.GRESN[num3]) - this.DELMIN)) * Math.Max((double)0, (double)(this.U[num3] - this.SMALLW))) / this.GRESN[num3];
            }
            goto Label_2088
        }
         
        num30 = 1;
        goto Label_2088
        Label_434F:this.ACCINF[0x10, this.ITSTEP] = this.XNORM;
        this.ACCINF[0x11, this.ITSTEP] = this.DNORM;
        this.ACCINF[0x12, this.ITSTEP] = this.PHASE;
        this.CFINCR = this.ICF;
        if (this.DIRDER >= 0)
        {
            this.STPTRM = -2;
            this.SIG = 0;
            goto Label_47EC
        }
         
        if (-this.DIRDER <= ((this.EPSMAC * 100) * (((this.SCF * Math.Abs(this.FX)) + this.PSI) + 1)))
        {
            if (this.UPSI > (this.DELMIN * this.NRES))
            {
                this.OPTITE = -1;
                this.STPTRM = -1;
            }
            else
            {
                this.OPTITE = 2;
                this.STPTRM = 1;
            } 
            this.SIG = 0;
            goto Label_47EC
        }
         
        if ((((this.PHASE >= 1) && (this.DNORM <= (this.SMALLD * (this.XNORM + this.SMALLD)))) && ((this.SCF0 == 1) && (num28 >= -this.SMALLW))) && !this.SINGUL)
        {
            this.PHASE = 2;
        }
         
        if ((this.PHASE == 2) && (this.DNORM > (this.XNORM + this.SMALLD)))
        {
            this.PHASE = 1;
        }
         
        this.dN13PF();
        num3 = 1;
        while (num3 <= this.N)
        {
            this.DD[num3] = 0;
            this.X1[num3 - 1] = this.X[num3 - 1] + this.D[num3];
            num3++;
        }
        flag4 = false;
        num3 = 1;
        while (num3 <= this.N)
        {
            if ((this.LLOW[num3] && (this.X1[num3 - 1] < (this.UG[num3] - this.TAUBND))) || (this.LUP[num3] && (this.X1[num3 - 1] > (this.OG[num3] + this.TAUBND))))
            {
                flag4 = true;
            }
             
            num3++;
        }
        if (((this.PHASE == 2) && (this.DNORM > (this.XNORM * Math.Sqrt(this.EPSMAC)))) && (!this.SINGUL && !flag4))
        {
            num3 = 1;
            while (num3 <= this.ALIST[0])
            {
                numArray3[num3] = 0;
                if ((num3 <= this.NH) && !this.GCONST[this.ALIST[num3]])
                {
                    this.CFUERR[num3] = false;
                    numArray3[num3] = this.DN40PF(num3, this.X1, FCN, XLB, XUB);
                    if (!this.CFUERR[num3])
                    {
                        goto Label_468F
                    }
                     
                    goto Label_47B2
                }
                 
                if (!this.GCONST[this.ALIST[num3]])
                {
                    this.CFUERR[this.ALIST[num3]] = false;
                    numArray3[num3] = this.DN42PF(this.ALIST[num3] - this.NH, this.X1, FCN, XLB, XUB);
                    if (this.CFUERR[this.ALIST[num3]])
                    {
                        goto Label_47B2
                    }
                     
                }
                 
                Label_468F:numArray3[num3] = -numArray3[num3];
                num3++;
            }
            num3 = 1;
            while (num3 <= this.ALIST[0])
            {
                numArray2[num3] = numArray3[this.COLNO[num3]];
                num3++;
            }
            this.DN26PF(1, this.NR, numArray2, this.DD);
            num3 = this.NR + 1;
            while (num3 <= this.N)
            {
                this.DD[num3] = 0;
                num3++;
            }
            this.DN24PF(-1, 0, 1, this.NR, this.N, this.QR, this.BETAQ, this.DD, numArray3);
            num3 = 1;
            while (num3 <= this.N)
            {
                this.DD[this.PERM[num3]] = numArray3[num3];
                num3++;
            }
            num14 = this.dN29PF(this.A,this.DD,this.DD,this.N);
            if (Math.Sqrt(num14) > (0.5 * this.DNORM))
            {
                for (num3 = 1;num3 <= this.N;num3++)
                {
                    this.DD[num3] = 0;
                }
            }
             
        }
         
        Label_47B2:if (!this.SILENT)
        {
            this.dN44PF(7);
        }
         
        this.SIG = Math.Min(1, this.STMAXL);
        this.DN17PF(this.SIG, FCN, XLB, XUB);
        Label_47EC:this.CFINCR = this.ICF - this.CFINCR;
        if (!this.SILENT)
        {
            this.dN44PF(10);
        }
         
        num14 = ((this.SCF * (this.FX0 - this.FX)) + this.PSI0) - this.PSI;
        if (Math.Abs(num14) <= (this.EPSPHI * ((this.SCF * Math.Abs(this.FX)) + this.PSI)))
        {
            num10++;
        }
        else
        {
            num10 = 0;
        } 
        if (num10 > this.NUMSM)
        {
            this.OPTITE = 7;
            this.throwError(this.OPTITE);
        }
        else
        {
            if (this.SIG <= 0.05)
            {
                if (this.SIG0 <= 0.05)
                {
                    num6++;
                }
                 
            }
            else
            {
                num6 = 0;
            } 
            this.ACCINF[0x15, this.ITSTEP] = this.SIG;
            this.ACCINF[0x16, this.ITSTEP] = this.CFINCR;
            this.ACCINF[0x17, this.ITSTEP] = this.DIRDER;
            this.ACCINF[0x18, this.ITSTEP] = this.DSCAL;
            this.ACCINF[0x19, this.ITSTEP] = this.COSPHI;
            this.ACCINF[0x1a, this.ITSTEP] = this.VIOLIS[0];
            if (((this.SIG == 0) && (this.STPTRM == 1)) && (this.OPTITE == 2))
            {
                if (!this.SILENT)
                {
                    this.dN44PF(0x11);
                }
                 
            }
            else
            {
                if ((((this.STPTRM == 1) && (this.SIG <= 0.0001)) && ((this.ACCINF[13, this.ITSTEP] > 10000) && !this.SINGUL)) && (this.NRES > 0))
                {
                    if (this.ACCINF[14, this.ITSTEP] > 10000)
                    {
                        this.dN10PF();
                    }
                     
                    this.IDENT = true;
                    this.SINGUL = true;
                }
                else
                {
                    if (this.STPTRM >= 0)
                    {
                        if ((((this.SINGUL && (this.ITSTEP > this.N)) && ((Math.Abs((double)(this.FX - this.FX0)) <= (num26 * (Math.Abs(this.FX) + 1))) && (this.PHASE >= 0))) && (((this.UPSI <= (this.NRES * this.DELMIN)) && (this.UPSI0 <= (this.NRES * this.DELMIN))) && ((num30 <= ((this.DELMIN * this.SMALLW) * this.NRES)) && (this.INFEAS <= this.UPSI)))) && !this.IDENT)
                        {
                            this.OPTITE = 4;
                            this.throwError(this.OPTITE);
                        }
                        else if (((this.SINGUL && (this.UPSI <= (this.DELMIN * this.NRES))) && ((this.UPSI0 <= (this.DELMIN * this.NRES)) && (this.B2N != -1))) && (((this.B2N <= (((this.GFN + 1) * this.EPSX) * 100)) && (this.PHASE >= 0)) && ((num30 <= ((this.DELMIN * this.SMALLW) * this.NRES)) && (this.INFEAS <= this.UPSI))))
                        {
                            this.OPTITE = 3;
                        }
                        else
                        {
                            num5 = 0;
                            num3 = 1;
                            while (num3 <= this.N)
                            {
                                if (Math.Abs(this.DifX[num3]) >= (this.EPSX * (Math.Abs(this.X[num3 - 1]) + 0.01)))
                                {
                                    num5 = 1;
                                }
                                 
                                num3++;
                            }
                            if (num5 == 0)
                            {
                                num9++;
                            }
                            else
                            {
                                num9 = 0;
                            } 
                            if ((num9 > this.NRESET) && this.SINGUL)
                            {
                                this.OPTITE = 5;
                                this.throwError(this.OPTITE);
                            }
                            else
                            {
                                this.XNORM = this.dN31PF(0,this.N - 1,this.X);
                                this.IDENT = false;
                                this.dN22PF(this.GPHI0);
                                num3 = 0;
                                while (num3 <= this.NRES)
                                {
                                    if (!this.GCONST[num3])
                                    {
                                        this.VAL[num3] = false;
                                    }
                                     
                                    num3++;
                                }
                                if ((this.PHASE >= 0) && !this.GCONST[0])
                                {
                                    this.VAL[0] = true;
                                    this.DN39PF(this.X, this.GRADF, FCN, XLB, XUB);
                                }
                                 
                                num3 = 1;
                                while (num3 <= this.ALIST[0])
                                {
                                    num1 = this.ALIST[num3];
                                    if (!this.VAL[num1])
                                    {
                                        this.VAL[num1] = true;
                                        if (num1 <= this.NH)
                                        {
                                            this.DN41PF(num1, this.X, numArray3, FCN, XLB, XUB);
                                        }
                                        else
                                        {
                                            this.DN43PF(num1 - this.NH, this.X, numArray3, FCN, XLB, XUB);
                                        } 
                                        for (num4 = 1;num4 <= this.N;num4++)
                                        {
                                            this.GRES[num1, num4] = numArray3[num4];
                                        }
                                        this.GRESN[num1] = Math.Max(1, this.DN31PF(1, this.N, numArray3));
                                    }
                                     
                                    num3++;
                                }
                                this.dN22PF(this.GPHI1);
                                num3 = 1;
                                while (num3 <= this.N)
                                {
                                    numArray3[num3] = this.X[num3 - 1] - this.X0[num3 - 1];
                                    numArray2[num3] = this.GPHI1[num3] - this.GPHI0[num3];
                                    num3++;
                                }
                                num14 = Math.Sqrt(this.DN31PF(1, this.N, numArray2) / this.DN31PF(1, this.N, numArray3));
                                if ((num14 != 0) && (this.PHASE >= 0))
                                {
                                    this.MATSC = Math.Max(1 / this.SCFMAX, Math.Min(this.SCFMAX, num14 / 2));
                                }
                                 
                                this.ALIST[0] = 0;
                                for (num3 = 1;num3 <= this.NRES;num3++)
                                {
                                    this.U0[num3] = this.U[num3];
                                    this.BINe0[num3] = this.BIND[num3];
                                    this.RES0[num3] = this.RES[num3];
                                    if (num3 <= this.NH)
                                    {
                                        this.BIND[num3] = 1;
                                    }
                                    else
                                    {
                                        this.BIND[num3] = 0;
                                        if ((this.RES[num3] / this.GRESN[num3]) <= this.DELMIN)
                                        {
                                            this.BIND[num3] = 1;
                                        }
                                         
                                        if (!this.VAL[num3] && (this.BIND[num3] == 1))
                                        {
                                            this.VAL[num3] = true;
                                            this.DN43PF(num3 - this.NH, this.X, numArray3, FCN, XLB, XUB);
                                            for (num4 = 1;num4 <= this.N;num4++)
                                            {
                                                this.GRES[num3, num4] = numArray3[num4];
                                            }
                                            this.GRESN[num3] = Math.Max(1, this.DN31PF(1, this.N, numArray3));
                                        }
                                         
                                    } 
                                    if (this.BIND[num3] == 1)
                                    {
                                        this.ALIST[0]++;
                                        this.ALIST[this.ALIST[0]] = num3;
                                    }
                                     
                                }
                                if (this.SCF != 0)
                                {
                                    if (((num7 > this.NRESET) || (num6 > this.NRESET)) || (num8 > this.NRESET))
                                    {
                                        num8 = 0;
                                        num6 = 0;
                                        num7 = 0;
                                        this.dN10PF();
                                    }
                                    else
                                    {
                                        this.dN9LPF();
                                        if (!this.SILENT)
                                        {
                                            this.dN44PF(14);
                                        }
                                         
                                    } 
                                }
                                 
                                if (this.ACCINF[0x1b, this.ITSTEP] == 1)
                                {
                                    if (((this.ITSTEP > 1) && (this.ACCINF[0x1d, this.ITSTEP - 1] != 0)) && (this.ACCINF[0x1d, this.ITSTEP] != 0))
                                    {
                                        num7++;
                                    }
                                    else
                                    {
                                        num7 = 0;
                                    } 
                                }
                                 
                                if (this.PHASE != -1)
                                {
                                    goto Label_05AA
                                }
                                 
                                if (this.UPSI >= this.TAU0)
                                {
                                    this.SCF = 0;
                                    this.PHASE = -1;
                                    goto Label_05AA
                                }
                                 
                                this.FFUERR = false;
                                this.FX = this.DN38PF(this.X, FCN, XLB, XUB);
                                if (this.FFUERR)
                                {
                                    this.OPTITE = -9;
                                    this.throwError(this.OPTITE);
                                }
                                else
                                {
                                    if (!this.VAL[0])
                                    {
                                        this.DN39PF(this.X, this.GRADF, FCN, XLB, XUB);
                                        this.VAL[0] = true;
                                    }
                                     
                                    this.SCF = 1;
                                    this.PHASE = 0;
                                    this.FXST = this.FX;
                                    this.PSIST = this.PSI;
                                    this.UPSIST = this.UPSI;
                                    for (num4 = 1;num4 <= this.NRES;num4++)
                                    {
                                        this.RESST[num4] = this.RES[num4];
                                    }
                                    this.ETA = 0;
                                    goto Label_05AA
                                } 
                            } 
                        }  
                        return ;
                    }
                     
                    if (!this.IDENT)
                    {
                        this.IDENT = true;
                        numArray5[0] = 0;
                        this.VIOLIS[0] = 0;
                        num8 = 0;
                        num6 = 0;
                        num7 = 0;
                        this.dN10PF();
                        this.ALIST[0] = num25;
                        for (num3 = 1;num3 <= this.NRES;num3++)
                        {
                            this.BIND[num3] = numArray6[num3];
                        }
                        if (this.UPSI >= this.TAU0)
                        {
                            if (this.UPSI >= this.TAU0)
                            {
                                this.SCF = 0;
                                this.PHASE = -1;
                            }
                            else
                            {
                                this.FFUERR = false;
                                this.FX = this.DN38PF(this.X, FCN, XLB, XUB);
                                if (this.FFUERR)
                                {
                                    this.OPTITE = -9;
                                    this.throwError(this.OPTITE);
                                    return ;
                                }
                                 
                                if (!this.VAL[0])
                                {
                                    this.DN39PF(this.X, this.GRADF, FCN, XLB, XUB);
                                    this.VAL[0] = true;
                                }
                                 
                                this.SCF = 1;
                                this.PHASE = 0;
                                this.FXST = this.FX;
                                this.PSIST = this.PSI;
                                this.UPSIST = this.UPSI;
                                for (num4 = 1;num4 <= this.NRES;num4++)
                                {
                                    this.RESST[num4] = this.RES[num4];
                                }
                                this.ETA = 0;
                            } 
                        }
                         
                        goto Label_05AA
                    }
                     
                    if ((!this.SINGUL && this.IDENT) && ((this.ACCINF[13, this.ITSTEP] > 10000) && (this.NRES > 0)))
                    {
                        this.SINGUL = true;
                        this.IDENT = true;
                    }
                    else
                    {
                        if (this.STPTRM == -2)
                        {
                            this.OPTITE = -4;
                            this.throwError(this.OPTITE);
                        }
                        else
                        {
                            if ((this.SIG == 0) && (this.OPTITE == -1))
                            {
                                this.throwError(this.OPTITE);
                                return ;
                            }
                             
                            this.OPTITE = -2;
                            this.throwError(this.OPTITE);
                        } 
                        return ;
                    } 
                } 
                this.SINGUL = true;
                this.PHASE = Math.Min(this.PHASE, 0);
                this.ACCINF[10, this.ITSTEP] = 1;
                this.SCF0 = 1;
                if ((this.PHASE >= 0) && (this.UPSI > (this.TAU0 * 0.5)))
                {
                    num17 = Math.Max(1 / this.SCFMAX, Math.Min(this.SCFMAX, ((((2 * (this.TAU0 - this.UPSI)) / this.TAU0) * this.UPSI) * this.TAU) / Math.Max(1, this.GFN)));
                    if ((((((this.FXST - this.FX) * num17) + ((num17 / this.SCF) * (this.PSIST - this.PSI))) >= (((num17 / this.SCF) * this.ETA) * this.CLOW)) && (this.LASTCH <= (this.ITSTEP - 4))) && ((num17 < (0.1 * this.SCF)) || (num17 > (10 * this.SCF))))
                    {
                        this.CLOW++;
                        num14 = num17 / this.SCF;
                        this.SCF0 = num14;
                        this.PSI *= num14;
                        this.PSIST *= num14;
                        num3 = 1;
                        while (num3 <= this.NRES)
                        {
                            this.U[num3] *= num14;
                            num3++;
                        }
                        num18 *= num14;
                        this.SCF = num17;
                        this.LASTCH = this.ITSTEP;
                        this.ACCINF[this.ITSTEP, 14] = this.SCF;
                        num14 = Math.Sqrt(num14);
                        for (num3 = 1;num3 <= this.N;num3++)
                        {
                            this.DIAG0[num3] = num14 * this.DIAG0[num3];
                            for (num4 = 1;num4 <= this.N;num4++)
                            {
                                this.A[num3, num4] *= num14;
                            }
                        }
                        this.MATSC *= num14;
                    }
                     
                }
                 
                this.ACCINF[0x20, this.ITSTEP] = this.UPSI;
                if (!this.SILENT)
                {
                    this.dN44PF(15);
                }
                 
                this.ACCINF[13, this.ITSTEP] = -1;
                num14 = Math.Abs(this.A[1, 1]);
                num16 = num14;
                num3 = 2;
                while (num3 <= this.N)
                {
                    num14 = Math.Max(num14, Math.Abs(this.A[num3, num3]));
                    num16 = Math.Min(num16, Math.Abs(this.A[num3, num3]));
                    num3++;
                }
                this.ACCINF[14, this.ITSTEP] = (num14 / num16) * (num14 / num16);
                if (!this.SILENT)
                {
                    this.dN44PF(5);
                }
                 
                num24 = this.CLOW;
                num31 = this.TAUQP;
                num3 = 1;
                while (num3 <= this.NRES)
                {
                    this.U[num3] = 0;
                    num3++;
                }
                num3 = 1;
                while (num3 <= this.N)
                {
                    this.DD[num3] = 0;
                    num3++;
                }
                this.dN32PF();
                if (((this.DNORM != 0) || (this.QPTERM != 1)) || (this.OPTITE != 3))
                {
                    if ((this.DNORM <= (this.EPSX * (Math.Min(this.XNORM, 1) + this.EPSX))) && (this.QPTERM < 0))
                    {
                        if ((this.UPSI >= (this.NRES * this.DELMIN)) && flag2)
                        {
                            this.OPTITE = -1;
                            if (!this.SILENT)
                            {
                                this.dN44PF(0x12);
                            }
                             
                            this.throwError(this.OPTITE);
                            return ;
                        }
                         
                        if (flag2)
                        {
                            this.OPTITE = this.QPTERM - 5;
                            if (!this.SILENT)
                            {
                                this.dN44PF(0x12);
                            }
                             
                            this.throwError(this.OPTITE);
                            return ;
                        }
                         
                        flag2 = true;
                        num3 = 1;
                        while (num3 <= this.NRES)
                        {
                            this.W[num3] = 1;
                            num3++;
                        }
                        this.LASTCH = this.ITSTEP;
                        this.dN10PF();
                        this.IDENT = true;
                        this.TAUQP = 1;
                        this.ALIST[0] = num25;
                        for (num3 = 1;num3 <= this.NRES;num3++)
                        {
                            this.BIND[num3] = numArray6[num3];
                        }
                        if (this.SCF == 0)
                        {
                            if (this.UPSI >= this.TAU0)
                            {
                                this.SCF = 0;
                                this.PHASE = -1;
                            }
                            else
                            {
                                this.FFUERR = false;
                                this.FX = this.DN38PF(this.X, FCN, XLB, XUB);
                                if (this.FFUERR)
                                {
                                    this.OPTITE = -9;
                                    this.throwError(this.OPTITE);
                                    return ;
                                }
                                 
                                if (!this.VAL[0])
                                {
                                    this.DN39PF(this.X, this.GRADF, FCN, XLB, XUB);
                                    this.VAL[0] = true;
                                }
                                 
                                this.SCF = 1;
                                this.PHASE = 0;
                                this.FXST = this.FX;
                                this.PSIST = this.PSI;
                                this.UPSIST = this.UPSI;
                                for (num4 = 1;num4 <= this.NRES;num4++)
                                {
                                    this.RESST[num4] = this.RES[num4];
                                }
                                this.ETA = 0;
                            } 
                        }
                         
                        goto Label_05AA
                    }
                     
                    this.ACCINF[11, this.ITSTEP] = 0;
                    numArray5[0] = 0;
                    num15 = 0;
                    this.dN46PF();
                    if (((this.QPTERM >= 0) || (this.QPTERM == -3)) && (this.SCF != 0))
                    {
                        num18 = Math.Abs(this.U[1]);
                        num3 = 2;
                        while (num3 <= this.NRES)
                        {
                            num18 = Math.Max(num18, Math.Abs(this.U[num3]));
                            num3++;
                        }
                        for (num3 = 1;num3 <= this.N;num3++)
                        {
                            numArray3[num3] = this.SCF * this.GRADF[num3];
                            for (num4 = 1;num4 <= this.NRES;num4++)
                            {
                                numArray3[num3] -= this.GRES[num4, num3] * this.U[num4];
                            }
                        }
                        this.B2N = this.DN31PF(1, this.N, numArray3) / this.SCF;
                    }
                    else
                    {
                        this.B2N = -1;
                    } 
                    if (!this.SILENT)
                    {
                        this.dN44PF(2);
                        this.dN44PF(0x10);
                    }
                     
                    if (this.B2N == -1)
                    {
                        this.B2N = this.EPSMAC / this.TOLMAC;
                    }
                     
                    if ((this.QPTERM >= 0) && (this.DNORM <= ((0.01 * this.EPSX) * (this.EPSX + Math.Min(1, this.XNORM)))))
                    {
                        if (this.UPSI <= (this.NRES * this.DELMIN))
                        {
                            this.OPTITE = 6;
                        }
                        else
                        {
                            this.OPTITE = -5;
                        } 
                        this.throwError(this.OPTITE);
                    }
                    else
                    {
                        if (this.QPTERM < 0)
                        {
                        }
                         
                        if (this.CLOW > num24)
                        {
                            num14 = 1;
                            for (num3 = 1;num3 <= this.NRES;num3++)
                            {
                                num14 = Math.Max(num14, this.W[num3]);
                            }
                            this.TAUQP = Math.Max(1, num14);
                        }
                         
                        if (this.TAUQP > Math.Pow(this.TAUFAC, 3 * num31))
                        {
                            this.TAUQP = Math.Pow(this.TAUFAC, 3 * num31);
                        }
                         
                        this.B2N0 = this.B2N;
                        num15 = 0;
                        num28 = 0;
                        if (this.QPTERM >= 1)
                        {
                            num30 = 0;
                            for (num3 = this.NH + 1;num3 <= this.NRES;num3++)
                            {
                                num30 += (Math.Max((double)0, (double)((this.RES[num3] / this.GRESN[num3]) - this.DELMIN)) * Math.Max((double)0, (double)(this.U[num3] - this.SMALLW))) / this.GRESN[num3];
                            }
                            goto Label_434F
                        }
                         
                        num30 = 1;
                        goto Label_434F
                    } 
                }
                 
            } 
        } 
    }

    private double dN40PF(int i, double[] X, Provisdom.Optimization.DoNlp2.IFunction FCN, double[] XLB, double[] XUB) throws Exception {
        double num2 = 0;
        if (this.BLOC)
        {
            return num2;
        }
         
        for (int num1 = 0;num1 < this.N;num1++)
        {
            this.XTR[num1] = X[num1] * this.XSC[num1];
        }
        return this.Evaluate(FCN, 2, i, this.XTR, XLB, XUB);
    }

    private void dN41PF(int i, double[] X, double[] GRADHI, Provisdom.Optimization.DoNlp2.IFunction FCN, double[] XLB, double[] XUB) throws Exception {
        int num2;
        double num1 = 45;
        double num9 = 0;
        double num10 = 0;
        double num11 = 0;
        double num12 = 0;
        double num13 = 0;
        double num14 = 0;
        double num15 = 0;
        double num16 = 0;
        double num17 = 0;
        double num18 = 0;
        if (this.GUNIT[i, 1] == 1)
        {
            for (num2 = 1;num2 <= this.N;num2++)
            {
                GRADHI[num2] = 0;
            }
            GRADHI[this.GUNIT[i, 2]] = this.GUNIT[i, 3] * this.XSC[this.GUNIT[i, 2] - 1];
            goto Label_054D
        }
         
        if (!this.BLOC)
        {
            num2 = 0;
            while (num2 < this.N)
            {
                this.XTR[num2] = this.XSC[num2] * X[num2];
                num2++;
            }
            if (!this.ANALYT)
            {
                if (this.DifFTYPE == 1)
                {
                    this.DELDif = Math.Min((double)(0.1 * Math.Sqrt(this.EPSFCN)), (double)0.01);
                    num18 = this.Evaluate(FCN, 2, i, this.XTR, XLB, XUB);
                    for (num2 = 0;num2 < this.N;num2++)
                    {
                        num17 = this.XTR[num2];
                        num16 = Math.Min((this.DELDif * Math.Abs(num17)) + this.DELDif, this.TAUBND);
                        num16 = Math.Min(0.01, num16);
                        if (num17 >= 0)
                        {
                            this.XTR[num2] = num17 + num16;
                        }
                        else
                        {
                            this.XTR[num2] = num17 - num16;
                        } 
                        num9 = this.Evaluate(FCN, 2, i, this.XTR, XLB, XUB);
                        GRADHI[num2 + 1] = (num9 - num18) / (this.XTR[num2] - num17);
                        this.XTR[num2] = num17;
                    }
                    goto Label_052B
                }
                 
                if (this.DifFTYPE == 2)
                {
                    this.DELDif = Math.Min((double)(0.1 * Math.Pow(this.EPSFCN, 0.33333333333333331)), (double)0.01);
                    for (num2 = 0;num2 < this.N;num2++)
                    {
                        num17 = this.XTR[num2];
                        num16 = Math.Min((this.DELDif * Math.Abs(num17)) + this.DELDif, this.TAUBND);
                        num16 = Math.Min(0.01, num16);
                        this.XTR[num2] = num17 + num16;
                        num10 = this.Evaluate(FCN, 2, i, this.XTR, XLB, XUB);
                        this.XTR[num2] = num17 - num16;
                        num11 = this.Evaluate(FCN, 2, i, this.XTR, XLB, XUB);
                        GRADHI[num2 + 1] = (num10 - num11) / (num16 + num16);
                        this.XTR[num2] = num17;
                    }
                    goto Label_052B
                }
                 
                this.DELDif = Math.Min((double)(0.1 * Math.Pow(this.EPSFCN, 0.14285714285714285)), (double)0.01);
                for (num2 = 0;num2 < this.N;num2++)
                {
                    num17 = this.XTR[num2];
                    num16 = Math.Min((double)((this.DELDif * Math.Abs(num17)) + this.DELDif), (double)(this.TAUBND / 4));
                    num16 = Math.Min(0.01, num16);
                    this.XTR[num2] = num17 - num16;
                    num10 = this.Evaluate(FCN, 2, i, this.XTR, XLB, XUB);
                    this.XTR[num2] = num17 + num16;
                    num11 = this.Evaluate(FCN, 2, i, this.XTR, XLB, XUB);
                    num16 += num16;
                    double num3 = num16;
                    this.XTR[num2] = num17 - num16;
                    num12 = this.Evaluate(FCN, 2, i, this.XTR, XLB, XUB);
                    this.XTR[num2] = num17 + num16;
                    num13 = this.Evaluate(FCN, 2, i, this.XTR, XLB, XUB);
                    num16 += num16;
                    double num4 = num16;
                    this.XTR[num2] = num17 - num16;
                    num14 = this.Evaluate(FCN, 2, i, this.XTR, XLB, XUB);
                    this.XTR[num2] = num17 + num16;
                    num15 = this.Evaluate(FCN, 2, i, this.XTR, XLB, XUB);
                    this.XTR[num2] = num17;
                    double num5 = num16 + num16;
                    double num6 = (num11 - num10) / num3;
                    double num7 = (num13 - num12) / num4;
                    double num8 = (num15 - num14) / num5;
                    num8 = num7 - num8;
                    num7 = num6 - num7;
                    num8 = num7 - num8;
                    GRADHI[num2 + 1] = (num6 + (0.4 * num7)) + (num8 / num1);
                }
                goto Label_052B
            }
             
            this.DN21PG((Provisdom.Optimization.DoNlp2.IGradient)this.F, 2, i, this.XTR, GRADHI);
            for (num2 = this.N;num2 >= 1;num2--)
            {
                GRADHI[num2] = this.XSC[num2 - 1] * GRADHI[num2 - 1];
            }
        }
         
        return ;
        Label_052B:num2 = 1;
        while (num2 <= this.N)
        {
            GRADHI[num2] = this.XSC[num2 - 1] * GRADHI[num2];
            num2++;
        }
        Label_054D:    ;
    }

    private double dN42PF(int i, double[] X, Provisdom.Optimization.DoNlp2.IFunction FCN, double[] XLB, double[] XUB) throws Exception {
        double num2 = 0;
        if (this.BLOC)
        {
            return num2;
        }
         
        for (int num1 = 0;num1 < this.N;num1++)
        {
            this.XTR[num1] = X[num1] * this.XSC[num1];
        }
        return this.Evaluate(FCN, 3, i, this.XTR, XLB, XUB);
    }

    private void dN43PF(int i, double[] X, double[] GRADGI, Provisdom.Optimization.DoNlp2.IFunction FCN, double[] XLB, double[] XUB) throws Exception {
        int num2;
        double num1 = 45;
        double num9 = 0;
        double num10 = 0;
        double num11 = 0;
        double num12 = 0;
        double num13 = 0;
        double num14 = 0;
        double num15 = 0;
        double num16 = 0;
        double num17 = 0;
        double num18 = 0;
        if (this.GUNIT[i + this.NH, 1] == 1)
        {
            for (num2 = 1;num2 <= this.N;num2++)
            {
                GRADGI[num2] = 0;
            }
            GRADGI[this.GUNIT[i + this.NH, 2]] = this.GUNIT[i + this.NH, 3] * this.XSC[this.GUNIT[i + this.NH, 2] - 1];
            goto Label_0567
        }
         
        if (!this.BLOC)
        {
            num2 = 0;
            while (num2 < this.N)
            {
                this.XTR[num2] = X[num2] * this.XSC[num2];
                num2++;
            }
            if (!this.ANALYT)
            {
                if (this.DifFTYPE == 1)
                {
                    this.DELDif = Math.Min((double)(0.1 * Math.Sqrt(this.EPSFCN)), (double)0.01);
                    num18 = this.Evaluate(FCN, 3, i, this.XTR, XLB, XUB);
                    for (num2 = 0;num2 < this.N;num2++)
                    {
                        num17 = this.XTR[num2];
                        num16 = Math.Min((this.DELDif * Math.Abs(num17)) + this.DELDif, this.TAUBND);
                        num16 = Math.Min(0.01, num16);
                        if (num17 >= 0)
                        {
                            this.XTR[num2] = num17 + num16;
                        }
                        else
                        {
                            this.XTR[num2] = num17 - num16;
                        } 
                        num9 = this.Evaluate(FCN, 3, i, this.XTR, XLB, XUB);
                        GRADGI[num2 + 1] = (num9 - num18) / (this.XTR[num2] - num17);
                        this.XTR[num2] = num17;
                    }
                    goto Label_0545
                }
                 
                if (this.DifFTYPE == 2)
                {
                    this.DELDif = Math.Min((double)(0.1 * Math.Pow(this.EPSFCN, 0.33333333333333331)), (double)0.01);
                    for (num2 = 0;num2 < this.N;num2++)
                    {
                        num17 = this.XTR[num2];
                        num16 = Math.Min((this.DELDif * Math.Abs(num17)) + this.DELDif, this.TAUBND);
                        num16 = Math.Min(0.01, num16);
                        this.XTR[num2] = num17 + num16;
                        num10 = this.Evaluate(FCN, 3, i, this.XTR, XLB, XUB);
                        this.XTR[num2] = num17 - num16;
                        num11 = this.Evaluate(FCN, 3, i, this.XTR, XLB, XUB);
                        GRADGI[num2 + 1] = (num10 - num11) / (num16 + num16);
                        this.XTR[num2] = num17;
                    }
                    goto Label_0545
                }
                 
                this.DELDif = Math.Min((double)(0.1 * Math.Pow(this.EPSFCN, 0.14285714285714285)), (double)0.01);
                for (num2 = 0;num2 < this.N;num2++)
                {
                    num17 = this.XTR[num2];
                    num16 = Math.Min((double)((this.DELDif * Math.Abs(num17)) + this.DELDif), (double)(this.TAUBND / 4));
                    num16 = Math.Min(0.01, num16);
                    this.XTR[num2] = num17 - num16;
                    num10 = this.Evaluate(FCN, 3, i, this.XTR, XLB, XUB);
                    this.XTR[num2] = num17 + num16;
                    num11 = this.Evaluate(FCN, 3, i, this.XTR, XLB, XUB);
                    num16 += num16;
                    double num3 = num16;
                    this.XTR[num2] = num17 - num16;
                    num12 = this.Evaluate(FCN, 3, i, this.XTR, XLB, XUB);
                    this.XTR[num2] = num17 + num16;
                    num13 = this.Evaluate(FCN, 3, i, this.XTR, XLB, XUB);
                    num16 += num16;
                    double num4 = num16;
                    this.XTR[num2] = num17 - num16;
                    num14 = this.Evaluate(FCN, 3, i, this.XTR, XLB, XUB);
                    this.XTR[num2] = num17 + num16;
                    num15 = this.Evaluate(FCN, 3, i, this.XTR, XLB, XUB);
                    this.XTR[num2] = num17;
                    double num5 = num16 + num16;
                    double num6 = (num11 - num10) / num3;
                    double num7 = (num13 - num12) / num4;
                    double num8 = (num15 - num14) / num5;
                    num8 = num7 - num8;
                    num7 = num6 - num7;
                    num8 = num7 - num8;
                    GRADGI[num2 + 1] = (num6 + (0.4 * num7)) + (num8 / num1);
                }
                goto Label_0545
            }
             
            this.DN21PG((Provisdom.Optimization.DoNlp2.IGradient)this.F, 3, i, this.XTR, GRADGI);
            for (num2 = this.N;num2 >= 1;num2--)
            {
                GRADGI[num2] = this.XSC[num2 - 1] * GRADGI[num2 - 1];
            }
        }
         
        return ;
        Label_0545:num2 = 1;
        while (num2 <= this.N)
        {
            GRADGI[num2] = this.XSC[num2 - 1] * GRADGI[num2];
            num2++;
        }
        Label_0567:    ;
    }

    private void dN44PF(int ICASE) throws Exception {
        Object[] objArray1 = new Object[10];
        if (this.TE2)
        {
            int num3;
            if (ICASE == 2)
            {
                this.keyMessage = "DoNlp2.FinerOutput1";
                objArray1[0] = this.DEL.ToString();
                objArray1[1] = this.B2N0.ToString();
                objArray1[2] = this.B2N.ToString();
                objArray1[3] = this.GFN.ToString();
                this.keyMessage = "DoNlp2.BlankLine";
                if (this.ALIST[0] != 0)
                {
                    this.keyMessage = "DoNlp2.FinerOutput2";
                    this.keyMessage = "DoNlp2.FinerOutput2A";
                    for (num3 = 1;num3 <= this.ALIST[0];num3++)
                    {
                        objArray1[0] = String.Concat(new String[]{ "(   ", Convert.ToString(this.ALIST[num3]), "   ", this.RES[this.ALIST[num3]].ToString(), "   ", this.GRESN[this.ALIST[num3]].ToString(), ")" });
                    }
                    this.keyMessage = "DoNlp2.BlankLine";
                }
                 
                if ((this.ALIST[0] != 0) && !this.SINGUL)
                {
                    this.keyMessage = "DoNlp2.FinerOutput3";
                    this.keyMessage = "DoNlp2.FinerOutput3A";
                    for (num3 = 1;num3 <= this.ALIST[0];num3++)
                    {
                        objArray1[0] = this.DIAG[num3].ToString();
                    }
                    this.keyMessage = "DoNlp2.BlankLine";
                }
                 
                if ((this.ALIST[0] == 0) || !this.TE3)
                {
                    return ;
                }
                 
                for (int num1 = 1;num1 <= this.ALIST[0];num1++)
                {
                    int num2 = this.ALIST[num1];
                    this.keyMessage = "DoNlp2.FinerOutput4";
                    objArray1[0] = Convert.ToString(num2);
                    this.keyMessage = "DoNlp2.FinerOutput4A";
                    objArray1[0] = this.GRES[num2, 1].ToString();
                    for (num3 = 2;num3 <= this.N;num3++)
                    {
                        objArray1[0] = objArray1[0] + "   " + this.GRES[num2, num3].ToString();
                    }
                }
            }
            else if (ICASE == 3)
            {
                if ((this.NR == 0) || (this.PHASE == -1))
                {
                    return ;
                }
                 
                this.keyMessage = "DoNlp2.FinerOutput5";
                this.keyMessage = "DoNlp2.FinerOutput5A";
                this.keyMessage = "DoNlp2.FinerOutput5B";
                for (num3 = 1;num3 <= this.NR;num3++)
                {
                    objArray1[0] = "   " + Convert.ToString(this.ALIST[num3]) + "   " + this.U[this.ALIST[num3]].ToString();
                }
                this.keyMessage = "DoNlp2.BlankLine";
            }
            else if (ICASE == 4)
            {
                if ((this.NR == 0) || (this.PHASE == -1))
                {
                    return ;
                }
                 
                this.keyMessage = "DoNlp2.BlankLine";
                this.keyMessage = "DoNlp2.FinerOutput6";
                this.keyMessage = "DoNlp2.FinerOutput5A";
                this.keyMessage = "DoNlp2.FinerOutput5B";
                for (num3 = 1;num3 <= this.NR;num3++)
                {
                    objArray1[0] = "   " + Convert.ToString(this.ALIST[num3]) + "   " + this.U[this.ALIST[num3]].ToString();
                }
            }
            else if (ICASE == 5)
            {
                this.keyMessage = "DoNlp2.FinerOutput7";
                objArray1[0] = this.ACCINF[13, this.ITSTEP].ToString();
                this.keyMessage = "DoNlp2.BlankLine";
                if (this.PHASE == -1)
                {
                    return ;
                }
                 
                this.keyMessage = "DoNlp2.FinerOutput8";
                objArray1[0] = this.ACCINF[14, this.ITSTEP].ToString();
                this.keyMessage = "DoNlp2.BlankLine";
            }
            else
            {
                if (ICASE == 6)
                {
                    return ;
                }
                 
                if (ICASE == 7)
                {
                    this.keyMessage = "DoNlp2.FinerOutput9";
                    objArray1[0] = Convert.ToString(this.PHASE);
                    objArray1[1] = this.SCF0.ToString();
                    this.keyMessage = "DoNlp2.FinerOutput9A";
                    this.keyMessage = "DoNlp2.FinerOutput4A";
                    objArray1[0] = this.D[1].ToString();
                    num3 = 2;
                    while (num3 <= this.N)
                    {
                        objArray1[0] = objArray1[0] + "   " + this.D[num3].ToString();
                        num3++;
                    }
                    this.keyMessage = "DoNlp2.BlankLine";
                    if (this.PHASE != 2)
                    {
                        return ;
                    }
                     
                    this.keyMessage = "DoNlp2.FinerOutput10";
                    this.keyMessage = "DoNlp2.FinerOutput4A";
                    objArray1[0] = this.DD[1].ToString();
                    for (num3 = 2;num3 <= this.N;num3++)
                    {
                        objArray1[0] = objArray1[0] + "   " + this.DD[num3].ToString();
                    }
                }
                else if (ICASE == 8)
                {
                    double num4 = this.TAU0 * 0.5;
                    double num5 = (this.FX * this.SCF) + this.PSI;
                    this.keyMessage = "DoNlp2.FinerOutput11";
                    this.keyMessage = "DoNlp2.BlankLine";
                    this.keyMessage = "DoNlp2.FinerOutput11A";
                    objArray1[0] = num5.ToString();
                    objArray1[1] = this.DIRDER.ToString();
                    objArray1[2] = this.PSI.ToString();
                    objArray1[3] = num4.ToString();
                    this.keyMessage = "DoNlp2.FinerOutput11B";
                    objArray1[0] = this.FX.ToString();
                    objArray1[1] = this.DSCAL.ToString();
                    objArray1[2] = this.SCF.ToString();
                    objArray1[3] = this.UPSI.ToString();
                }
                else if (ICASE == 9)
                {
                    this.keyMessage = "DoNlp2.FinerOutput11C";
                    objArray1[0] = this.SIG.ToString();
                    objArray1[1] = this.FX1.ToString();
                    objArray1[2] = this.PSI1.ToString();
                    objArray1[3] = this.UPSI1.ToString();
                    this.keyMessage = "DoNlp2.BlankLine";
                }
                else
                {
                    if (ICASE == 10)
                    {
                        return ;
                    }
                     
                    if (ICASE == 11)
                    {
                        return ;
                    }
                     
                    if (ICASE == 12)
                    {
                        return ;
                    }
                     
                    if (ICASE == 13)
                    {
                        return ;
                    }
                     
                    if (ICASE == 14)
                    {
                        return ;
                    }
                     
                    if (ICASE == 15)
                    {
                        return ;
                    }
                     
                    if (ICASE == 0x10)
                    {
                        return ;
                    }
                     
                    if (ICASE == 0x11)
                    {
                        return ;
                    }
                     
                    if (ICASE == 0x12)
                    {
                        return ;
                    }
                     
                    if (ICASE == 0x13)
                    {
                        return ;
                    }
                     
                    if (ICASE == 20)
                    {
                        return ;
                    }
                     
                    if (ICASE == 0x15)
                    {
                    }
                     
                }   
            }    
        }
         
    }

    private void dN45PF(double[][] A, int ME, int NE, int MA, int NA, String HEAD, int LOGNUM, boolean FIX) throws Exception {
        Object[] objArray1 = new Object[NA + 1];
        int num5 = 4;
        int num3 = 0;
        while (num3 < NA)
        {
            int num4 = num3 + 1;
            num3 = Math.Min((num4 + num5) - 1, NA);
            this.keyMessage = "DoNlp2.FinestOutput2";
            objArray1[0] = "ROW/COLUMN      ";
            int num2 = num4;
            while (num2 <= num3)
            {
                objArray1[0] = objArray1[0] + Convert.ToString(num2) + "     ";
                num2++;
            }
            for (int num1 = 1;num1 <= MA;num1++)
            {
                if (FIX)
                {
                    this.keyMessage = "DoNlp2.FinestOutput3";
                    objArray1[0] = "   " + Convert.ToString(num1) + "      ";
                    for (num2 = num4;num2 <= num3;num2++)
                    {
                        objArray1[0] = objArray1[0] + A[num2, num1].ToString() + "  ";
                    }
                }
                else
                {
                    this.keyMessage = "DoNlp2.FinestOutput4";
                    objArray1[0] = "   " + Convert.ToString(num1) + "   ";
                    for (num2 = num4;num2 <= num3;num2++)
                    {
                        objArray1[0] = objArray1[0] + A[num2, num1].ToString() + "  ";
                    }
                } 
            }
        }
    }

    private void dN46PF() throws Exception {
        Object[] objArray1 = new Object[7];
        if (this.TE0)
        {
            double num1 = this.ACCINF[11, this.ITSTEP];
            this.keyMessage = "DoNlp2.BlankLine";
            this.keyMessage = "DoNlp2.ConfigOutput";
            objArray1[0] = Convert.ToString(this.ITSTEP);
            objArray1[1] = this.FX.ToString();
            objArray1[2] = this.UPSI.ToString();
            objArray1[3] = this.B2N.ToString();
            objArray1[4] = num1.ToString();
            objArray1[5] = Convert.ToString(this.NR);
            objArray1[6] = Convert.ToString((int)this.ACCINF[10, this.ITSTEP]);
            this.keyMessage = "DoNlp2.BlankLine";
        }
         
    }

    private double evaluate(Provisdom.Optimization.DoNlp2.IFunction FCN, int FTYPE, int i, double[] XTR, double[] XLB, double[] XUB) throws Exception {
        double num3 = 0;
        int num1 = (this.NG + this.NH) + this.NBOUNDS;
        if (FTYPE == 1)
        {
            this.ICF++;
            this.TMPBUL[0] = this.FFUERR;
            num3 = FCN.F(XTR, 0, this.TMPBUL);
            this.FFUERR = this.TMPBUL[0];
            return num3;
        }
         
        if (FTYPE == 2)
        {
            this.CRES[i]++;
            this.TMPBUL[0] = this.CFUERR[i];
            num3 = FCN.F(XTR, i, this.TMPBUL);
            this.CFUERR[i] = this.TMPBUL[0];
            return num3;
        }
         
        if (FTYPE == 3)
        {
            if (this.GUNIT[i + this.NH, 1] == -1)
            {
                this.CRES[i + this.NH]++;
            }
             
            if (i > (this.NG - this.NBOUNDS))
            {
                int num2 = this.GUNIT[i + this.NH, 2];
                if (this.GUNIT[i + this.NH, 3] > 0)
                {
                    return (XTR[num2 - 1] - XLB[num2 - 1]);
                }
                 
                return (XUB[num2 - 1] - XTR[num2 - 1]);
            }
             
            this.TMPBUL[0] = this.CFUERR[i + this.NH];
            num3 = FCN.F(XTR, i + this.NH, this.TMPBUL);
            this.CFUERR[i + this.NH] = this.TMPBUL[0];
        }
         
        return num3;
    }

    private void throwError(double OPTITE) throws Exception {
        int num1 = (int)OPTITE;
        if (num1 == -10)
        {
            throw new ConstraintEvaluationException();
        }
         
        if (num1 == -9)
        {
            throw new ObjectiveEvaluationException();
        }
         
        if (num1 == -8)
        {
            throw new WorkingSetSingularException();
        }
         
        if (num1 == -7)
        {
            throw new QPInfeasibleException();
        }
         
        if ((num1 == -6) || (num1 == -5))
        {
            throw new PenaltyFunctionPointInfeasibleException();
        }
         
        if (num1 == -4)
        {
            throw new LimitingAccuracyException();
        }
         
        if (num1 == -3)
        {
            throw new TooManyIterationsException(this.ITSTEP);
        }
         
        if (num1 == -2)
        {
            Object[] objArray1 = new Object[]{ this.SIGSM, this.SIGLA };
            throw new NoAcceptableStepsizeException(this.SIGSM, this.SIGLA);
        }
         
        if (num1 == -1)
        {
            throw new BadInitialGuessException();
        }
         
        if (num1 == 4)
        {
            throw new IllConditionedException();
        }
         
        if (num1 == 5)
        {
            throw new SingularException();
        }
         
        if (num1 == 6)
        {
            throw new LinearlyDependentGradientsException();
        }
         
        if (num1 == 7)
        {
            throw new TooManyIterationsException(this.NUMSM);
        }
         
    }

    private void dN6LPF() throws Exception {
        int num1;
        String[] textArray1 = new String[20];
        Object[] objArray1 = new Object[10];
        textArray1[1] = "CONSTRAINT EVALUATION }S ERROR WITH CURRENT POINT";
        textArray1[2] = "OBJECTIVE EVALUATION }S ERROR WITH CURRENT POINT";
        textArray1[3] = "QPSOLVER: WORKING SET SINGULAR IN DUAL EXTED QP ";
        textArray1[4] = "QPSOLVER: EXTED QP-PROBLEM SEEMINGLY INFEASIBLE ";
        textArray1[5] = "QPSOLVER: NO DESCENT DIRECTION FROM QP FOR TAU=TAU_MAX";
        textArray1[6] = "QPSOLVER: ON EXIT CORRECTION SMALL, INFEASIBLE POINT";
        textArray1[7] = "STEPSIZESELECTION: COMPUTED D NOT A DIRECTION OF DESCENT";
        textArray1[8] = "MORE THAN MAXIT ITERATION STEPS";
        textArray1[9] = "STEPSIZESELECTION: NO ACCEPTABLE STEPSIZE IN [SIGSM,SIGLA]";
        textArray1[10] = "STEPSIZESELECTION: DIRECTIONAL DERIV. VERY SMALL, INFEASIBLE";
        textArray1[11] = "KT-CONDITIONS SATISFIED, NO FURTHER CORRECTION COMPUTED";
        textArray1[12] = "KT-CONDITIONS SATISFIED, COMPUTED CORRECTION SMALL";
        textArray1[13] = "STEPSIZESELECTION: X (ALMOST) FEASIBLE, DIR. DERIV. VERY SMALL";
        textArray1[14] = "KT-CONDITIONS (RELAXED) SATISFIED, SINGULAR POINT";
        textArray1[15] = "VERY SLOW PRIMAL PROGRESS, SINGULAR OR ILLCONDITONED PROBLEM";
        textArray1[0x10] = "VERY SLOW PROGRESS IN X, SINGULAR PROBLEM";
        textArray1[0x11] = "CORRECTION VERY SMALL, ALMOST FEASIBLE BUT SINGULAR POINT";
        textArray1[0x12] = "NUMSM  SMALL DifFERENCES IN PENALTY FUNCTION,TERMINATE";
        if (this.SCF != 0)
        {
            for (num1 = 1;num1 <= this.NRES;num1++)
            {
                this.U[num1] /= this.SCF;
            }
        }
         
        if (!this.SILENT || this.INTAKT)
        {
            num1 = 0;
            int num2 = 0;
            double num12 = 0;
            int num3 = 1;
            while (num3 <= this.NRES)
            {
                num1 += this.CRES[num3];
                num2 += this.CGRES[num3];
                if (num3 > this.NH)
                {
                    num12 = Math.Min(num12, this.U[num3]);
                }
                 
                num3++;
            }
            int num15 = num1;
            int num16 = num2;
            int num14 = 0;
            int num19 = 0;
            int num17 = 0;
            int num18 = 0;
            num3 = 1;
            while (num3 <= this.ITSTEP)
            {
                if (this.ACCINF[10, num3] == 1)
                {
                    num14++;
                }
                 
                if (this.ACCINF[0x1b, num3] == 1)
                {
                    num18++;
                }
                 
                if ((this.ACCINF[0x1d, num3] == 0) && (this.ACCINF[0x1b, num3] == 1))
                {
                    num17++;
                }
                 
                if (this.ACCINF[0x1b, num3] == -1)
                {
                    num19++;
                }
                 
                num3++;
            }
            num3 = ((int)this.OPTITE) + 11;
            switch(num3)
            {
                case 1: 
                    {
                        this.keyMessage = "DoNlp2.FineOutput1";
                        break;
                    }
                case 2: 
                    {
                        this.keyMessage = "DoNlp2.FineOutput2";
                        break;
                    }
                case 3: 
                    {
                        this.keyMessage = "DoNlp2.FineOutput3";
                        break;
                    }
                case 4: 
                    {
                        this.keyMessage = "DoNlp2.FineOutput4";
                        break;
                    }
                case 5: 
                    {
                        this.keyMessage = "DoNlp2.FineOutput5";
                        break;
                    }
                case 6: 
                    {
                        this.keyMessage = "DoNlp2.FineOutput6";
                        break;
                    }
                case 7: 
                    {
                        this.keyMessage = "DoNlp2.FineOutput7";
                        break;
                    }
                case 8: 
                    {
                        this.keyMessage = "DoNlp2.FineOutput8";
                        break;
                    }
                case 9: 
                    {
                        this.keyMessage = "DoNlp2.FineOutput9";
                        break;
                    }
                case 10: 
                    {
                        this.keyMessage = "DoNlp2.FineOutput10";
                        break;
                    }
                case 11: 
                    {
                        this.keyMessage = "DoNlp2.FineOutput11";
                        break;
                    }
                case 12: 
                    {
                        this.keyMessage = "DoNlp2.FineOutput12";
                        break;
                    }
                case 13: 
                    {
                        this.keyMessage = "DoNlp2.FineOutput13";
                        break;
                    }
                case 14: 
                    {
                        this.keyMessage = "DoNlp2.FineOutput14";
                        break;
                    }
                case 15: 
                    {
                        this.keyMessage = "DoNlp2.FineOutput15";
                        break;
                    }
                case 0x10: 
                    {
                        this.keyMessage = "DoNlp2.FineOutput16";
                        break;
                    }
                case 0x11: 
                    {
                        this.keyMessage = "DoNlp2.FineOutput17";
                        break;
                    }
                case 0x12: 
                    {
                        this.keyMessage = "DoNlp2.FineOutput18";
                        break;
                    }
                default: 
                    {
                        this.keyMessage = "DoNlp2.FineOutput19";
                        break;
                    }
            
            }
            if (!this.SILENT)
            {
            }
             
            if ((this.ITSTEP > 1) && (this.OPTITE == 0))
            {
                this.ITSTEP--;
            }
             
            double num13 = this.ACCINF[14, this.ITSTEP];
            num1 = this.ITSTEP;
            while ((num1 > 1) && (num13 == 0))
            {
                num1--;
                num13 = this.ACCINF[14, num1];
            }
            if (this.SILENT)
            {
                this.TE1 = false;
            }
             
            if (this.TE1)
            {
                for (num1 = 1;num1 <= this.ITSTEP;num1++)
                {
                    int num4 = (int)this.ACCINF[1, num1];
                    int num5 = (int)this.ACCINF[9, num1];
                    int num6 = (int)this.ACCINF[10, num1];
                    int num7 = (int)this.ACCINF[0x12, num1];
                    int num8 = (int)this.ACCINF[0x13, num1];
                    int num9 = (int)this.ACCINF[0x16, num1];
                    int num10 = (int)this.ACCINF[0x1a, num1];
                    int num11 = (int)this.ACCINF[0x1b, num1];
                }
            }
             
        }
         
    }

    private void dN7LPF(Provisdom.Optimization.DoNlp2.IFunction FCN, double[] XLB, double[] XUB) throws Exception {
        int num2;
        double num5;
        double[] numArray1 = new double[this.NX];
        char[] chArray1 = new char[8];
        char[] chArray2 = new char[8];
        char[] chArray3 = new char[8];
        int num1 = 0;
        while (num1 < 8)
        {
            chArray2[num1] = 'X';
            num1++;
        }
        this.EPSMAC = 2.2204460492503131E-16;
        this.TOLMAC = 2.2250738585072014E-308;
        if (this.NRESET > this.N)
        {
            this.NRESET = this.N;
        }
         
        if (this.NRESET <= 4)
        {
            this.NRESET = 4;
        }
         
        this.LASTCH = 0;
        this.LASTDW = 0;
        this.LASTUP = 0;
        this.LEVEL = 1;
        this.TAU = 0.1;
        this.ITERMA = this.MAXIT;
        this.EPSX = 1E-05;
        this.SIGSM = Math.Sqrt(this.EPSMAC);
        this.SMALLD = 0.1;
        this.RHO = 1E-06;
        this.RHO1 = 1E-10;
        this.DEL01 = this.DEL0 / 10;
        this.C1D = 0.01;
        this.TAUFAC = 10;
        this.TAUMAX = this.SCFMAX * this.SCFMAX;
        double num4 = this.EPSMAC / this.TOLMAC;
        this.FX = 0;
        this.B2N = 0;
        this.B2N0 = 0;
        this.NRES = this.NG + this.NH;
        if (this.COLD)
        {
            for (num1 = 1;num1 <= this.NX;num1++)
            {
                for (num2 = 1;num2 <= this.NX;num2++)
                {
                    this.A[num1, num2] = 0;
                }
                this.A[num1, num1] = 1;
                this.DIAG0[num1] = 1;
            }
        }
         
        num1 = 1;
        while (num1 <= this.NRESM)
        {
            this.DIAG[num1] = 0;
            for (num2 = 1;num2 <= this.NX;num2++)
            {
                this.QR[num1, num2] = 0;
                this.GRES[num1, num2] = 0;
            }
            num1++;
        }
        num1 = 1;
        while (num1 <= this.NX)
        {
            numArray1[num1 - 1] = 0;
            this.UG[num1] = -num4;
            this.OG[num1] = num4;
            this.LLOW[num1] = false;
            this.LUP[num1] = false;
            num1++;
        }
        num1 = 1;
        while (num1 <= this.NH)
        {
            this.DELFAC[num1] = 1;
            num1++;
        }
        num1 = this.NH + 1;
        while (num1 <= this.NRES)
        {
            this.DELFAC[num1] = 1;
            if (this.GUNIT[num1, 1] == 1)
            {
                num5 = this.DN42PF(num1 - this.NH, numArray1, FCN, XLB, XUB);
                if (this.GUNIT[num1, 3] > 0)
                {
                    this.LLOW[this.GUNIT[num1, 2]] = true;
                    this.UG[this.GUNIT[num1, 2]] = -num5 / ((double)this.GUNIT[num1, 3]);
                }
                else
                {
                    this.OG[this.GUNIT[num1, 2]] = -num5 / ((double)this.GUNIT[num1, 3]);
                    this.LUP[this.GUNIT[num1, 2]] = true;
                } 
            }
             
            num1++;
        }
        num1 = this.NH + 1;
        while (num1 <= this.NRES)
        {
            if (this.GUNIT[num1, 1] == 1)
            {
                num2 = this.GUNIT[num1, 2];
                if ((this.OG[num2] < num4) && (this.UG[num2] > -num4))
                {
                    this.DEL0 = Math.Min(this.DEL0, ((this.OG[num2] - this.UG[num2]) * 0.1) * Math.Abs(this.GUNIT[num1, 3]));
                }
                 
            }
             
            num1++;
        }
        num1 = this.NH + 1;
        while (num1 <= this.NRES)
        {
            if (this.GUNIT[num1, 1] == 1)
            {
                num2 = this.GUNIT[num1, 2];
                if (this.GUNIT[num1, 3] > 0)
                {
                    this.DELFAC[num1] = Math.Max(this.DELFAC[num1], Math.Abs(this.UG[num2]) * 0.1);
                    if (this.OG[num2] < num4)
                    {
                        this.DELFAC[num1] = Math.Min(this.DELFAC[num1], (this.OG[num2] - this.UG[num2]) / (10 * this.DEL0));
                    }
                     
                }
                else
                {
                    this.DELFAC[num1] = Math.Max(this.DELFAC[num1], Math.Abs(this.OG[num2]) * 0.1);
                    if (this.UG[num2] > -num4)
                    {
                        this.DELFAC[num1] = Math.Min(this.DELFAC[num1], (this.OG[num2] - this.UG[num2]) / (10 * this.DEL0));
                    }
                     
                } 
            }
             
            num1++;
        }
        double num3 = num4;
        num1 = 1;
        while (num1 <= this.N)
        {
            if (this.UG[num1] > 0)
            {
                num3 = Math.Min(num3, this.OG[num1]);
            }
             
            if (this.OG[num1] < 0)
            {
                num3 = Math.Min(num3, -this.UG[num1]);
            }
             
            num1++;
        }
        if (this.GUNIT[0, 1] == 1)
        {
            this.GCONST[0] = true;
            this.VAL[0] = true;
            for (num1 = 1;num1 <= this.N;num1++)
            {
                this.GRADF[num1] = 0;
            }
            this.GRADF[this.GUNIT[0, 2]] = this.GUNIT[0, 3] * this.XSC[this.GUNIT[0, 2] - 1];
            this.GFN = Math.Abs(this.GUNIT[0, 3]);
        }
        else
        {
            this.VAL[0] = false;
            for (num1 = 1;num1 <= this.N;num1++)
            {
                this.GRADF[num1] = 0;
            }
        } 
        num1 = 1;
        while (num1 <= this.NH)
        {
            if (this.GUNIT[num1, 1] == 1)
            {
                this.VAL[num1] = true;
                this.GCONST[num1] = true;
                this.GRES[num1, this.GUNIT[num1, 2]] = this.GUNIT[num1, 3] * this.XSC[this.GUNIT[num1, 2] - 1];
                this.GRESN[num1] = Math.Abs(this.GUNIT[num1, 3]) * this.XSC[this.GUNIT[num1, 2] - 1];
                double num6 = this.DN40PF(num1, numArray1, FCN, XLB, XUB);
                double num7 = -num6 / ((double)this.GUNIT[num1, 3]);
                this.X[this.GUNIT[num1, 2] - 1] = num7;
            }
             
            num1++;
        }
        num1 = this.NH + 1;
        while (num1 <= this.NRES)
        {
            if (this.GUNIT[num1, 1] == 1)
            {
                num5 = this.DN42PF(num1 - this.NH, this.X, FCN, XLB, XUB);
                num5 = (2 * this.DELMIN) - num5;
                if (num5 > 0)
                {
                    this.X[this.GUNIT[num1, 2] - 1] += num5 / ((double)this.GUNIT[num1, 3]);
                }
                 
                this.GRES[num1, this.GUNIT[num1, 2]] = this.GUNIT[num1, 3] * this.XSC[this.GUNIT[num1, 2] - 1];
                this.GRESN[num1] = Math.Abs(this.GUNIT[num1, 3]) * this.XSC[this.GUNIT[num1, 2] - 1];
                this.VAL[num1] = true;
                this.GCONST[num1] = true;
            }
             
            num1++;
        }
        for (num1 = 1;num1 <= this.NRES;num1++)
        {
            this.BIND[num1] = 0;
            this.BINe0[num1] = 0;
            this.U[num1] = 0;
            this.U0[num1] = 0;
            this.CRES[num1] = 0;
            this.CGRES[num1] = 0;
            if (this.COLD)
            {
                this.W[num1] = 1;
            }
             
            this.SORT[num1] = num1;
        }
        this.CLOW = 1;
        this.NY = 2;
        this.SCF = 1;
        this.SCF0 = 1;
        this.SIGLA = 2048;
        this.BETA = 4;
        this.DELTA1 = 0.9;
        this.DELTA = 0.001;
        this.THETA = 0.9;
        this.ICF = 0;
        this.ICGF = 0;
    }

    private void dN8LPF() throws Exception {
        double num2;
        boolean flag1 = false;
        int num1 = 1;
        while (num1 <= this.NRES)
        {
            num2 = (this.NY * Math.Abs(this.U[num1])) + this.TAU;
            if (num2 > this.W[num1])
            {
                this.W1[num1] = num2 + this.TAU;
            }
            else
            {
                this.W1[num1] = this.W[num1];
                if ((num2 < (this.W[num1] * 0.5)) && (this.BIND[num1] == 1))
                {
                    this.W1[num1] = (num2 + this.W[num1]) * 0.5;
                }
                 
            } 
            if (this.W1[num1] < this.W[num1])
            {
                flag1 = true;
            }
             
            num1++;
        }
        double num3 = 0;
        double num4 = 0;
        num1 = 1;
        while (num1 <= this.NRES)
        {
            if (num1 <= this.NH)
            {
                num3 += this.W1[num1] * Math.Abs(this.RESST[num1]);
                num4 += this.W1[num1] * Math.Abs(this.RES[num1]);
            }
            else
            {
                num3 -= Math.Min(0, this.RESST[num1]) * this.W1[num1];
                num4 -= Math.Min(0, this.RES[num1]) * this.W1[num1];
            } 
            num1++;
        }
        double num5 = ((this.FXST - this.FX) * this.SCF) + (num3 - num4);
        if ((flag1 && (num5 >= (this.ETA * this.CLOW))) && ((this.ITSTEP - this.LASTDW) > Math.Max(5, Math.Min(20, this.N / 10))))
        {
            if (this.CLOW > (this.ITSTEP / 10))
            {
                this.ETA = 1.3 * this.ETA;
                if (!this.SILENT)
                {
                    this.dN44PF(11);
                }
                 
            }
             
            this.LASTCH = this.ITSTEP;
            this.LASTDW = this.ITSTEP;
            this.LEVEL = num5 / ((double)this.ITERMA);
            this.PSIST = num3;
            this.PSI = num4;
            for (num1 = 1;num1 <= this.NRES;num1++)
            {
                this.W[num1] = this.W1[num1];
            }
            this.CLOW++;
        }
        else
        {
            num3 = 0;
            num4 = 0;
            for (num1 = 1;num1 <= this.NRES;num1++)
            {
                if (this.W1[num1] > this.W[num1])
                {
                    this.LASTUP = this.ITSTEP;
                    this.LASTCH = this.ITSTEP;
                }
                 
                this.W[num1] = Math.Max(this.W[num1], this.W1[num1]);
                if (num1 <= this.NH)
                {
                    num3 += this.W[num1] * Math.Abs(this.RESST[num1]);
                    num4 += this.W[num1] * Math.Abs(this.RES[num1]);
                }
                else
                {
                    num3 -= this.W[num1] * Math.Min(0, this.RESST[num1]);
                    num4 -= this.W[num1] * Math.Min(0, this.RES[num1]);
                } 
            }
            this.PSIST = num3;
            this.PSI = num4;
        } 
        num2 = 0;
        if (this.NRES >= 1)
        {
            num2 = this.W[0];
        }
         
        for (num1 = 2;num1 <= this.NRES;num1++)
        {
            num2 = Math.Max(num2, this.W[num1]);
        }
        this.ACCINF[20, this.ITSTEP] = num2;
        this.ACCINF[0x13, this.ITSTEP] = this.CLOW;
        if (!this.SILENT)
        {
            this.dN44PF(12);
        }
         
    }

    private void dN9LPF() throws Exception {
        double[] numArray1 = new double[this.NX + 1];
        double[] numArray2 = new double[this.NX + 1];
        double[] numArray3 = new double[this.NX + 1];
        double[] numArray4 = new double[this.NRESM + 1];
        double[] numArray5 = new double[this.NX + 1];
        double[] numArray6 = new double[this.NX + 1];
        boolean flag1 = false;
        int num1 = 1;
        while (num1 <= this.N)
        {
            numArray3[num1] = this.dN19PF(num1,this.N,num1,this.A,this.NX,this.DifX);
            numArray1[num1] = this.GPHI1[num1] - this.GPHI0[num1];
            num1++;
        }
        if (this.DN31PF(1, this.N, numArray1) == 0)
        {
            this.ACCINF[0x1b, this.ITSTEP] = 0;
            this.ACCINF[0x1c, this.ITSTEP] = 0;
            this.ACCINF[0x1d, this.ITSTEP] = 0;
        }
        else
        {
            int num2;
            double num12 = new double();
            num1 = 1;
            while (num1 <= this.N)
            {
                numArray2[num1] = this.DN20PF(1, num1, num1, this.A, this.NX, numArray3);
                num1++;
            }
            num1 = 1;
            while (num1 <= this.ALIST[0])
            {
                numArray4[num1] = this.DN20PF(1, this.N, this.ALIST[num1], this.GRES, this.NX, this.DifX);
                numArray4[num1] /= this.GRESN[this.ALIST[num1]];
                num1++;
            }
            double num13 = this.dN31PF(1,this.N,this.DifX);
            double num7 = Math.Min((double)0.5, (double)(this.DNORM * this.DNORM));
            double num11 = 0;
            double num10 = Math.Abs(this.A[1, 1]);
            num11 = 0;
            num1 = 1;
            while (num1 <= this.N)
            {
                for (num2 = num1;num2 <= this.N;num2++)
                {
                    num11 += this.A[num2, num1] * this.A[num2, num1];
                }
                num10 = Math.Min(num10, Math.Abs(this.A[num1, num1]));
                num1++;
            }
            if (num10 == 0)
            {
                num12 = this.EPSMAC / this.TOLMAC;
            }
            else
            {
                num12 = (num11 / num10) * num10;
            } 
            double num3 = this.DN31PF(1, this.N, numArray3) * this.DN31PF(1, this.N, numArray3);
            double num4 = this.DN18PF(1, this.N, numArray1, this.DifX);
            if ((num3 <= (((this.RHO1 * num11) * num13) * num13)) || (num12 >= (1 / this.RHO1)))
            {
                this.dN10PF();
            }
            else
            {
                double num6;
                double num9;
                if (this.NRES == 0)
                {
                    num6 = 1;
                    if (num4 < (0.2 * num3))
                    {
                        num6 = (0.8 * num3) / (num3 - num4);
                        for (num1 = 1;num1 <= this.N;num1++)
                        {
                            numArray1[num1] = (num6 * numArray1[num1]) + ((1 - num6) * numArray2[num1]);
                        }
                        num4 = this.DN18PF(1, this.N, numArray1, this.DifX);
                    }
                     
                    num9 = 1 / Math.Sqrt(num4);
                    num1 = 1;
                    while (num1 <= this.N)
                    {
                        numArray1[num1] *= num9;
                        numArray6[num1] = numArray1[num1];
                        num1++;
                    }
                    num9 = 1 / Math.Sqrt(num3);
                    for (num1 = 1;num1 <= this.N;num1++)
                    {
                        numArray5[num1] = numArray2[num1] * num9;
                    }
                    this.ACCINF[0x1c, this.ITSTEP] = num4 / num3;
                    this.ACCINF[0x1d, this.ITSTEP] = num6;
                    this.ACCINF[0x1b, this.ITSTEP] = 2;
                    if (num6 != 1)
                    {
                        this.ACCINF[0x1b, this.ITSTEP] = 3;
                    }
                     
                }
                else
                {
                    double num8;
                    double num15 = new double();
                    double num14 = this.DN31PF(1, this.ALIST[0], numArray4);
                    num9 = 1 / Math.Sqrt(num3);
                    num1 = 1;
                    while (num1 <= this.N)
                    {
                        numArray5[num1] = numArray2[num1] * num9;
                        num1++;
                    }
                    if ((num4 >= (this.RHO1 * this.DN18PF(1, this.N, numArray1, numArray1))) && (this.DN31PF(1, this.N, numArray1) >= (Math.Sqrt(this.EPSMAC) * num13)))
                    {
                        num8 = 0;
                        for (num1 = 1;num1 <= this.N;num1++)
                        {
                            numArray6[num1] = numArray1[num1];
                        }
                        num15 = num4;
                    }
                    else
                    {
                        double num5 = ((num7 * num13) * num13) + (num14 * num14);
                        if (num4 >= (this.RHO1 * this.DN18PF(1, this.N, numArray1, numArray1)))
                        {
                            num8 = 1;
                        }
                        else
                        {
                            num8 = 1 + ((((num7 * num13) * num13) + Math.Abs(num4)) / num5);
                        } 
                        for (num1 = 1;num1 <= this.N;num1++)
                        {
                            num9 = 0;
                            for (num2 = 1;num2 <= this.ALIST[0];num2++)
                            {
                                num10 = this.GRES[this.ALIST[num2], num1] * numArray4[num2];
                                num10 /= this.GRESN[this.ALIST[num2]];
                                num9 += num10;
                            }
                            numArray6[num1] = numArray1[num1] + (num8 * ((num7 * this.DifX[num1]) + num9));
                        }
                        num15 = this.DN18PF(1, this.N, numArray6, this.DifX);
                    } 
                    num9 = 1 / Math.Sqrt(num15);
                    num1 = 1;
                    while (num1 <= this.N)
                    {
                        numArray6[num1] *= num9;
                        num1++;
                    }
                    num6 = 1;
                    if (num4 < (0.2 * num3))
                    {
                        num6 = (0.8 * num3) / (num3 - num4);
                        for (num1 = 1;num1 <= this.N;num1++)
                        {
                            numArray1[num1] = (num6 * numArray1[num1]) + ((1 - num6) * numArray2[num1]);
                        }
                        num4 = this.DN18PF(1, this.N, numArray1, this.DifX);
                    }
                     
                    num9 = 1 / Math.Sqrt(num4);
                    num1 = 1;
                    while (num1 <= this.N)
                    {
                        numArray1[num1] *= num9;
                        num1++;
                    }
                    if (this.DN31PF(1, this.N, numArray1) <= (0.001 * this.DN31PF(1, this.N, numArray6)))
                    {
                        for (num1 = 1;num1 <= this.N;num1++)
                        {
                            numArray6[num1] = numArray1[num1];
                        }
                        this.ACCINF[0x1c, this.ITSTEP] = num4 / num3;
                        this.ACCINF[0x1d, this.ITSTEP] = num6;
                        this.ACCINF[0x1b, this.ITSTEP] = 2;
                        if (num6 != 1)
                        {
                            this.ACCINF[0x1b, this.ITSTEP] = 3;
                        }
                         
                    }
                    else
                    {
                        this.ACCINF[0x1b, this.ITSTEP] = 1;
                        this.ACCINF[0x1c, this.ITSTEP] = num7;
                        this.ACCINF[0x1d, this.ITSTEP] = num8;
                    } 
                } 
                flag1 = this.DN28PF(this.A, numArray6, numArray5, this.N);
                num9 = Math.Abs(this.A[1, 1]);
                num10 = num9;
                num1 = 1;
                while (num1 < this.N)
                {
                    num1++;
                    num9 = Math.Max(num9, Math.Abs(this.A[num1, num1]));
                    num10 = Math.Min(num10, Math.Abs(this.A[num1, num1]));
                }
                if (!flag1 && ((num10 * num10) > ((this.RHO1 * num9) * num9)))
                {
                    return ;
                }
                 
                this.dN10PF();
            } 
        } 
    }

    public double[] getConstraintResiduals() throws Exception {
        if (this.RES == null)
        {
            return null;
        }
         
        double[] numArray1 = new double[this.NRESM - this.NBOUNDS];
        for (int num1 = 0;num1 < (this.NRESM - this.NBOUNDS);num1++)
        {
            numArray1[num1] = this.RES[num1 + 1];
        }
        return numArray1;
    }

    public double getProjectedGradientNorm() throws Exception {
        return B2N;
    }

    public double getConstraintNorm() throws Exception {
        return UPSI;
    }

    public double[] getLagrangeMultiplierEst() throws Exception {
        if (this.U == null)
        {
            return null;
        }
         
        double[] numArray1 = new double[this.NRESM - this.NBOUNDS];
        for (int num1 = 0;num1 < (this.NRESM - this.NBOUNDS);num1++)
        {
            numArray1[num1] = this.U[num1 + 1];
        }
        return numArray1;
    }

    public double[] getGuess() throws Exception {
        return (double[])XGUESS;
    }

    public void setGuess(double[] value) throws Exception {
        XGUESS = (double[])value.Clone();
    }

    public double[] getLowerBounds() throws Exception {
        return (double[])_lowerBounds.Clone();
    }

    public void setLowerBounds(double[] value) throws Exception {
        _lowerBounds = (double[])value.Clone();
        for (int i = 0;i < _lowerBounds.Length;i++)
        {
            _lowerBounds[i] = double.IsNegativeInfinity(_lowerBounds[i]) ? double.MinValue : _lowerBounds[i];
        }
    }

    public double[] getScale() throws Exception {
        return _scale;
    }

    public void setScale(double[] value) throws Exception {
        for (int num1 = 0;num1 < value.GetLength(0);num1++)
        {
            if (value[num1] <= 0)
            {
                Object[] objArray1 = new Object[]{ "" };
                ExceptionThrower.ThrowArgumentException("Numerical.Math", "XscaleNonPositive", objArray1);
            }
             
        }
        this._scale = (double[])value.Clone();
    }

    public double[] getUpperBounds() throws Exception {
        return _upperBounds;
    }

    public void setUpperBounds(double[] value) throws Exception {
        this._upperBounds = (double[])value.Clone();
        for (int i = 0;i < this._upperBounds.Length;i++)
        {
            _upperBounds[i] = double.IsPositiveInfinity(_upperBounds[i]) ? double.MaxValue : _upperBounds[i];
        }
    }

    public double[] solve(Provisdom.Optimization.DoNlp2.IFunction f) throws Exception {
        this.F = f;
        double[] numArray1 = new double[this._numberOfVariables + 1];
        double[] numArray2 = new double[this._numberOfVariables];
        boolean flag1 = false;
        Provisdom.Optimization.DoNlp2.IGradient gradient1 = (this.F instanceof Provisdom.Optimization.DoNlp2.IGradient) ? ((Provisdom.Optimization.DoNlp2.IGradient)this.F) : null;
        if (this.DEL_MIN != this.INITIAL_DEL_MIN)
        {
            flag1 = true;
        }
         
        if ((this.PENALTY != 1) && (this.DELTA0 == 0.5))
        {
            this.DELTA0 = 0.5 * this.PENALTY;
        }
         
        if (!flag1)
        {
            this.DEL_MIN = Math.Min(this.DELTA0 / 10, Math.Max(1E-06 * this.DELTA0, this.SMALL_W));
            if (gradient1 == null)
            {
                this.DEL_MIN = Math.Min(0.1 * this.DELTA0, Math.Max(this.RELPRE, this.DEL_MIN));
            }
             
        }
         
        this.N = this._numberOfVariables;
        this.NH = this._numberOfEqualityConstraints;
        this.NG = this._totalConstraints - this._numberOfEqualityConstraints;
        int num2 = 0;
        while (num2 < this._numberOfVariables)
        {
            if (this._lowerBounds[num2] >= this._upperBounds[num2])
            {
                Object[] objArray1 = new Object[]{ "IllegalBounds", num2, this._lowerBounds[num2], this._upperBounds[num2] };
                ExceptionThrower.ThrowArgumentException("Numerical.Math", "IllegalBounds", objArray1);
            }
             
            num2++;
        }
        this.NBOUNDS = 0;
        int num1 = 0;
        while (num1 < this._numberOfVariables)
        {
            if (this._lowerBounds[num1] != -1.7976931348623157E+308)
            {
                this.NBOUNDS++;
            }
             
            if (this._upperBounds[num1] != 1.7976931348623157E+308)
            {
                this.NBOUNDS++;
            }
             
            num1++;
        }
        this.NG += this.NBOUNDS;
        this.NX = this.N;
        this.MAXIT = this._maxIterations;
        this.NRESM = this.NH + this.NG;
        double[][] numArray3 = new double[0x21, this.MAXIT + 1];
        this.ACCINF = numArray3;
        double[] numArray4 = new double[this.NX];
        this.X = numArray4;
        double[] numArray5 = new double[this.NX];
        this.X0 = numArray5;
        double[] numArray6 = new double[this.NX];
        this.X1 = numArray6;
        double[] numArray7 = new double[this.NX];
        this.XMIN = numArray7;
        double[] numArray8 = new double[this.NRESM + 1];
        this.RESMIN = numArray8;
        double[] numArray9 = new double[this.NX + 1];
        this.D = numArray9;
        double[] numArray10 = new double[this.NX + 1];
        this.E0 = numArray10;
        double[] numArray11 = new double[this.NX + 1];
        this.DD = numArray11;
        double[] numArray12 = new double[this.NX + 1];
        this.DifX = numArray12;
        double[] numArray13 = new double[this.NX + 1];
        this.GRADF = numArray13;
        double[] numArray14 = new double[this.NX + 1];
        this.QGF = numArray14;
        double[][] numArray15 = new double[this.NRESM + 1, this.NX + 1];
        this.GRES = numArray15;
        double[] numArray16 = new double[this.NRESM + 1];
        this.GRESN = numArray16;
        double[] numArray17 = new double[this.NX + 1];
        this.GPHI0 = numArray17;
        double[] numArray18 = new double[this.NX + 1];
        this.GPHI1 = numArray18;
        double[][] numArray19 = new double[this.NRESM + 1, this.NX + 1];
        this.QR = numArray19;
        double[] numArray20 = new double[this.NRESM + 1];
        this.BETAQ = numArray20;
        double[] numArray21 = new double[this.NRESM + 1];
        this.DIAG = numArray21;
        double[] numArray22 = new double[this.NRESM + 1];
        this.CSCAL = numArray22;
        double[] numArray23 = new double[this.NRESM + 1];
        this.COLLE = numArray23;
        int[] numArray24 = new int[(2 * this.NRESM) + 1];
        this.COLNO = numArray24;
        int[] numArray25 = new int[this.NX + 1];
        this.PERM = numArray25;
        int[] numArray26 = new int[this.NX + 1];
        this.PERM1 = numArray26;
        boolean[] flagArray1 = new boolean[this.NRESM + 1];
        this.VAL = flagArray1;
        boolean[] flagArray2 = new boolean[1];
        this.TMPBUL = flagArray2;
        boolean[] flagArray3 = new boolean[this.NRESM + 1];
        this.GCONST = flagArray3;
        int[][] numArray27 = new int[this.NRESM + 1, 4];
        this.GUNIT = numArray27;
        boolean[] flagArray4 = new boolean[this.NX + 1];
        this.LLOW = flagArray4;
        boolean[] flagArray5 = new boolean[this.NX + 1];
        this.LUP = flagArray5;
        double[][] numArray28 = new double[this.NX + 1, this.NX + 1];
        this.A = numArray28;
        double[] numArray29 = new double[this.NX + 1];
        this.DIAG0 = numArray29;
        int[] numArray30 = new int[this.NRESM + 1];
        this.BIND = numArray30;
        int[] numArray31 = new int[this.NRESM + 1];
        this.BINe0 = numArray31;
        int[] numArray32 = new int[(this.NSTEP * this.NRESM) + 1];
        this.VIOLIS = numArray32;
        int[] numArray33 = new int[this.NRESM + 1];
        this.ALIST = numArray33;
        int[] numArray34 = new int[this.NRESM + 1];
        this.SORT = numArray34;
        double[] numArray35 = new double[this.NRESM + 1];
        this.RES = numArray35;
        double[] numArray36 = new double[this.NRESM + 1];
        this.RES0 = numArray36;
        double[] numArray37 = new double[this.NRESM + 1];
        this.RES1 = numArray37;
        double[] numArray38 = new double[this.NRESM + 1];
        this.RESST = numArray38;
        double[] numArray39 = new double[this.NRESM + 1];
        this.U = numArray39;
        double[] numArray40 = new double[this.NRESM + 1];
        this.U0 = numArray40;
        double[] numArray41 = new double[this.NRESM + 1];
        this.W = numArray41;
        double[] numArray42 = new double[this.NRESM + 1];
        this.W1 = numArray42;
        double[] numArray43 = new double[this.NRESM + 1];
        this.WORK = numArray43;
        double[] numArray44 = new double[this.NRESM + 1];
        this.YU = numArray44;
        double[] numArray45 = new double[this.NRESM + 1];
        this.SLACK = numArray45;
        double[] numArray46 = new double[this.NX + 1];
        this.UG = numArray46;
        double[] numArray47 = new double[this.NX + 1];
        this.OG = numArray47;
        double[] numArray48 = new double[this.NRESM + 1];
        this.DELFAC = numArray48;
        double[] numArray49 = new double[this.NX + 1];
        this.XST = numArray49;
        int[] numArray50 = new int[this.NRESM + 1];
        this.CRES = numArray50;
        int[] numArray51 = new int[this.NRESM + 1];
        this.CGRES = numArray51;
        boolean[] flagArray6 = new boolean[this.NRESM + 1];
        this.CFUERR = flagArray6;
        double[] numArray52 = new double[this.NX];
        this.XTR = numArray52;
        double[] numArray53 = new double[this.NX];
        this.XSC = numArray53;
        double[] numArray54 = new double[this.NRESM + 1];
        this.FU = numArray54;
        double[][] numArray55 = new double[this.NRESM + 1, this.NX + 1];
        this.FUGRAD = numArray55;
        double[][] numArray56 = new double[7, this.NRESM + 1];
        this.FUD = numArray56;
        double[][] numArray57 = new double[(this.NRESM + this.NX) + 1, (this.NRESM + this.NX) + 1];
        this.XJ = numArray57;
        double[] numArray58 = new double[(this.NX + this.NRESM) + 1];
        this.DDUAL = numArray58;
        double[][] numArray59 = new double[(this.NX + this.NRESM) + 1, (this.NX + this.NRESM) + 1];
        this.R = numArray59;
        double[] numArray60 = new double[(this.NX + this.NRESM) + 1];
        this.NP = numArray60;
        double[] numArray61 = new double[(2 * this.NRESM) + 1];
        this.UD = numArray61;
        double[] numArray62 = new double[(2 * this.NRESM) + 1];
        this.UD1 = numArray62;
        int[] numArray63 = new int[(2 * this.NRESM) + 1];
        this.AITR = numArray63;
        boolean flag2 = false;
        num1 = 0;
        while (num1 < this._numberOfVariables)
        {
            if (this.XGUESS[num1] != 0)
            {
                flag2 = true;
            }
             
            num1++;
        }
        if (flag2)
        {
            Array.Copy(this.XGUESS, 0, this.X, 0, this._numberOfVariables);
        }
        else
        {
            for (num1 = 0;num1 < this._numberOfVariables;num1++)
            {
                if ((this._lowerBounds[num1] <= 0) && (this._upperBounds[num1] >= 0))
                {
                    this.X[num1] = 0;
                }
                else if (this._lowerBounds[num1] >= 0)
                {
                    this.X[num1] = this._lowerBounds[num1];
                }
                else if (this._upperBounds[num1] <= 0)
                {
                    this.X[num1] = this._upperBounds[num1];
                }
                   
            }
        } 
        num2 = 0;
        while (num2 < this._numberOfVariables)
        {
            if (this._scale[num2] <= 0)
            {
            }
             
            num2++;
        }
        num2 = 0;
        while (num2 <= this.NRESM)
        {
            this.GUNIT[num2, 1] = -1;
            this.GUNIT[num2, 2] = 0;
            this.GUNIT[num2, 3] = 0;
            num2++;
        }
        if (this.NG > this._totalConstraints)
        {
            num2 = this._totalConstraints + 1;
            num1 = 1;
            while (num1 <= this._numberOfVariables)
            {
                if (this._lowerBounds[num1 - 1] != -1.7976931348623157E+308)
                {
                    this.GUNIT[num2, 1] = 1;
                    this.GUNIT[num2, 2] = num1;
                    this.GUNIT[num2, 3] = 1;
                    num2++;
                }
                 
                num1++;
            }
            for (num1 = 1;num1 <= this._numberOfVariables;num1++)
            {
                if (this._upperBounds[num1 - 1] != 1.7976931348623157E+308)
                {
                    this.GUNIT[num2, 1] = 1;
                    this.GUNIT[num2, 2] = num1;
                    this.GUNIT[num2, 3] = -1;
                    num2++;
                }
                 
            }
        }
         
        this.BLOC = false;
        this.ANALYT = false;
        if (this.F instanceof Provisdom.Optimization.DoNlp2.IGradient)
        {
            this.ANALYT = true;
        }
         
        this.EPSFCN = this.FEPS;
        this.DifFTYPE = this.IDTYPE;
        this.TAUBND = this.TBND;
        num1 = 0;
        while (num1 < this.NX)
        {
            this.XSC[num1] = 1;
            this.XTR[num1] = 0;
            num1++;
        }
        num1 = 0;
        while (num1 <= this.NRESM)
        {
            this.GCONST[num1] = false;
            this.VAL[num1] = false;
            if (num1 > 0)
            {
                this.GRESN[num1] = 1;
            }
             
            num1++;
        }
        num1 = 1;
        while (num1 <= this.NRESM)
        {
            this.CFUERR[num1] = false;
            num1++;
        }
        this.FFUERR = false;
        this.SILENT = false;
        this.INTAKT = false;
        this.TE0 = false;
        this.TE1 = false;
        this.TE2 = false;
        this.TE3 = false;
        this.COLD = true;
        Array.Copy(this._scale, 0, this.XSC, 0, this._numberOfVariables);
        this.TAU0 = this.PENALTY;
        this.DEL0 = this.DELTA0;
        this.DELMIN = this.DEL_MIN;
        this.SMALLW = this.SMALL_W;
        this.SCFMAX = this.SCF_MAX;
        this.EPSDif = this.RELPRE;
        this.dN7LPF(this.F,this._lowerBounds,this._upperBounds);
        this.NUMSM = Math.Max(this.N, 10);
        num1 = 1;
        while (num1 <= this.N)
        {
            this.XST[num1] = this.X[num1 - 1];
            this.X[num1 - 1] /= this.XSC[num1 - 1];
            num1++;
        }
        this.NRESET = this.N;
        this.EPSPHI = 1000 * this.EPSMAC;
        num1 = 1;
        while (num1 <= this.N)
        {
            if (this.LLOW[num1])
            {
                this.UG[num1] /= this.XSC[num1 - 1];
            }
             
            if (this.LUP[num1])
            {
                this.OG[num1] /= this.XSC[num1 - 1];
            }
             
            num1++;
        }
        num1 = 0;
        while (num1 <= this.NRES)
        {
            if ((this.GUNIT[num1, 1] != 1) && this.GCONST[num1])
            {
                if (num1 == 0)
                {
                    this.VAL[0] = true;
                    this.dN39PF(this.X,this.GRADF,this.F,this._lowerBounds,this._upperBounds);
                }
                else
                {
                    this.VAL[num1] = true;
                    if (num1 <= this.NH)
                    {
                        this.DN41PF(num1, this.X, numArray1, this.F, this._lowerBounds, this._upperBounds);
                    }
                    else
                    {
                        this.DN43PF(num1 - this.NH, this.X, numArray1, this.F, this._lowerBounds, this._upperBounds);
                    } 
                    for (num2 = 1;num2 <= this.N;num2++)
                    {
                        this.GRES[num1, num2] = numArray1[num2];
                    }
                    this.GRESN[num1] = Math.Max(1, this.DN31PF(1, this.N, numArray1));
                } 
            }
             
            num1++;
        }
        for (num1 = 0;num1 < this.N;num1++)
        {
            numArray2[num1] = double.NaN;
        }
        this.FVALUE = double.NaN;
        this.solve(this.F,this.IPRINT,this._maxIterations,this._lowerBounds,this._upperBounds);
        this.FVALUE = this.FX;
        Array.Copy(this.X, 0, numArray2, 0, this.N);
        if (this.IPRINT >= 1)
        {
            this.dN6LPF();
        }
         
        return numArray2;
    }

    public double getBindingThreshold() throws Exception {
        return this.DELTA0;
    }

    public void setBindingThreshold(double value) throws Exception {
        if (value <= 0)
        {
            Object[] objArray1 = new Object[]{ "Binding threshold", value };
            ExceptionThrower.ThrowArgumentException("Numerical.Math", "NotPositive", objArray1);
        }
        else
        {
            this.DELTA0 = value;
        } 
    }

    public double getBoundViolationBound() throws Exception {
        return this.TBND;
    }

    public void setBoundViolationBound(double value) throws Exception {
        if (value <= 0)
        {
            Object[] objArray1 = new Object[]{ "Bound violation bound", value };
            ExceptionThrower.ThrowArgumentException("Numerical.Math", "NotPositive", objArray1);
        }
        else
        {
            this.TBND = value;
        } 
    }

    public int getDifferentiationType() throws Exception {
        return this.IDTYPE;
    }

    public void setDifferentiationType(int value) throws Exception {
        if ((value <= 0) || (value >= 4))
        {
            Object[] objArray1 = new Object[]{ "Differentiation type", value, "1 <= idtype <= 3" };
            ExceptionThrower.ThrowArgumentException("Numerical.Math", "NotInInterval", objArray1);
        }
        else
        {
            this.IDTYPE = value;
        } 
    }

    public double getFunctionPrecision() throws Exception {
        return this.FEPS;
    }

    public void setFunctionPrecision(double value) throws Exception {
        if (value <= 0)
        {
            Object[] objArray1 = new Object[]{ "Function precision", value };
            ExceptionThrower.ThrowArgumentException("Numerical.Math", "NotPositive", objArray1);
        }
        else
        {
            this.FEPS = value;
        } 
    }

    public double getGradientPrecision() throws Exception {
        return this.RELPRE;
    }

    public void setGradientPrecision(double value) throws Exception {
        if (value <= 0)
        {
            Object[] objArray1 = new Object[]{ "Gradient precision", value };
            ExceptionThrower.ThrowArgumentException("Numerical.Math", "NotPositive", objArray1);
        }
        else
        {
            this.RELPRE = value;
        } 
    }

    public int getMaximumIterations() throws Exception {
        return this._maxIterations;
    }

    public void setMaximumIterations(int value) throws Exception {
        if (value <= 0)
        {
            Object[] objArray1 = new Object[]{ "MaximumIterations", value };
            ExceptionThrower.ThrowArgumentException("Numerical.Math", "NotPositive", objArray1);
        }
        else
        {
            this._maxIterations = value;
        } 
    }

    public double getMultiplierError() throws Exception {
        return this.SMALL_W;
    }

    public void setMultiplierError(double value) throws Exception {
        if (value <= 0)
        {
            Object[] objArray1 = new Object[]{ "Multiplier error", value };
            ExceptionThrower.ThrowArgumentException("Numerical.Math", "NotPositive", objArray1);
        }
        else
        {
            this.SMALL_W = value;
        } 
    }

    public double getPenaltyBound() throws Exception {
        return this.PENALTY;
    }

    public void setPenaltyBound(double value) throws Exception {
        if (value <= 0)
        {
            Object[] objArray1 = new Object[]{ "Penalty bound", value };
            ExceptionThrower.ThrowArgumentException("Numerical.Math", "NotPositive", objArray1);
        }
        else
        {
            this.PENALTY = value;
        } 
    }

    public double getScalingBound() throws Exception {
        return this.SCF_MAX;
    }

    public void setScalingBound(double value) throws Exception {
        if (value <= 0)
        {
            Object[] objArray1 = new Object[]{ "Scaling variable", value };
            ExceptionThrower.ThrowArgumentException("Numerical.Math", "NotPositive", objArray1);
        }
        else
        {
            this.SCF_MAX = value;
        } 
    }

    public double getViolationBound() throws Exception {
        return this.DEL_MIN;
    }

    public void setViolationBound(double value) throws Exception {
        if (value <= 0)
        {
            Object[] objArray1 = new Object[]{ "Violation bound", value };
            ExceptionThrower.ThrowArgumentException("Numerical.Math", "NotPositive", objArray1);
        }
        else
        {
            this.DEL_MIN = value;
        } 
    }

    private double[][] A = new double[][]();
    private double[][] ACCINF = new double[][]();
    private int[] AITR = new int[]();
    private int[] ALIST = new int[]();
    private boolean ANALYT = new boolean();
    private double B2N = new double();
    private double B2N0 = new double();
    private double BETA = new double();
    private double[] BETAQ = new double[]();
    private int[] BIND = new int[]();
    private int[] BINe0 = new int[]();
    private boolean BLOC = new boolean();
    private double C1D = new double();
    private int CFINCR = new int();
    private boolean[] CFUERR = new boolean[]();
    private int[] CGRES = new int[]();
    private int CLOW = new int();
    private boolean COLD = new boolean();
    private double[] COLLE = new double[]();
    private int[] COLNO = new int[]();
    private double COSPHI = new double();
    private int[] CRES = new int[]();
    private double[] CSCAL = new double[]();
    private double[] D = new double[]();
    private double[] DD = new double[]();
    private double[] DDUAL = new double[]();
    private double DEL = new double();
    private double DEL_MIN = new double();
    private double DEL0 = new double();
    private double DEL01 = new double();
    private double DELDif = new double();
    private double[] DELFAC = new double[]();
    private double DELMIN = new double();
    private double DELTA = new double();
    private double DELTA0 = new double();
    private double DELTA1 = new double();
    private double[] DIAG = new double[]();
    private double[] DIAG0 = new double[]();
    private int DifFTYPE = new int();
    private double[] DifX = new double[]();
    private double DIRDER = new double();
    private double DNORM = new double();
    private double DSCAL = new double();
    private double[] E0 = new double[]();
    private double e0NORM = new double();
    private double EPSDif = new double();
    private double EPSFCN = new double();
    private double EPSMAC = new double();
    private double EPSPHI = new double();
    private double EPSX = new double();
    private boolean EQRES = new boolean();
    private double ETA = new double();
    private Provisdom.Optimization.DoNlp2.IFunction F;
    private double FEPS = new double();
    private boolean FFUERR = new boolean();
    private static final double FIVE = 5;
    private double FMIN = new double();
    private static final double FOUR = 4;
    private double[] FU = new double[]();
    private double[][] FUD = new double[][]();
    private double[][] FUGRAD = new double[][]();
    private double FVALUE = new double();
    private double FX = new double();
    private double FX0 = new double();
    private double FX1 = new double();
    private double FXST = new double();
    private boolean[] GCONST = new boolean[]();
    private double GFN = new double();
    private double[] GPHI0 = new double[]();
    private double[] GPHI1 = new double[]();
    private double[] GRADF = new double[]();
    private double[][] GRES = new double[][]();
    private double[] GRESN = new double[]();
    private int[][] GUNIT = new int[][]();
    private int ICF = new int();
    private int ICGF = new int();
    private boolean IDENT = new boolean();
    private int IDTYPE = new int();
    private double INFEAS = new double();
    private double INITIAL_DEL_MIN = new double();
    private boolean INTAKT = new boolean();
    private int IPRINT = new int();
    private int IPTR = new int();
    private int IQ = new int();
    private int IQTR = new int();
    private int ITERMA = new int();
    private int ITSTEP = new int();
    public String keyMessage = new String();
    private int LASTCH = new int();
    private int LASTDW = new int();
    private int LASTUP = new int();
    private double LEVEL = new double();
    private boolean[] LLOW = new boolean[]();
    private boolean[] LUP = new boolean[]();
    private int _totalConstraints = new int();
    private double MATSC = new double();
    private int MAXIT = new int();
    private int _maxIterations = new int();
    private int _numberOfEqualityConstraints = new int();
    private int MI = new int();
    private int N = new int();
    private int NBOUNDS = new int();
    private int NDUAL = new int();
    private int NG = new int();
    private int NH = new int();
    private double[] NP = new double[]();
    private int NR = new int();
    private int NRES = new int();
    private int NRESET = new int();
    private int NRESM = new int();
    private int NSTEP = new int();
    private int NUMSM = new int();
    private int _numberOfVariables = new int();
    private int NX = new int();
    private double NY = new double();
    private double[] OG = new double[]();
    private static final double ONE = 1;
    private static final double ONEP3 = 1.3;
    private double OPTITE = new double();
    private static final double P2 = 0.2;
    private static final double P4 = 0.4;
    private static final double P5 = 0.5;
    private static final double P7 = 0.7;
    private static final double P8 = 0.8;
    private static final double P9 = 0.9;
    private double PENALTY = new double();
    private int[] PERM = new int[]();
    private int[] PERM1 = new int[]();
    private int PHASE = new int();
    private double PHI = new double();
    private double PHI1 = new double();
    private double PHIMIN = new double();
    private double PSI = new double();
    private double PSI0 = new double();
    private double PSI1 = new double();
    private double PSIMIN = new double();
    private double PSIST = new double();
    private double[] QGF = new double[]();
    private int QPTERM = new int();
    private double[][] QR = new double[][]();
    private double[][] R = new double[][]();
    private int RANK = new int();
    private double RELPRE = new double();
    private double[] RES = new double[]();
    private double[] RES0 = new double[]();
    private double[] RES1 = new double[]();
    private double[] RESMIN = new double[]();
    private double[] RESST = new double[]();
    private double RHO = new double();
    private double RHO1 = new double();
    private double RIITR = new double();
    private double RLOW = new double();
    private double RNORM = new double();
    private double SCF = new double();
    private double SCF_MAX = new double();
    private double SCF0 = new double();
    private double SCFMAX = new double();
    private static final double SEVEN = 7;
    private double SIG = new double();
    private double SIG0 = new double();
    private double SIGLA = new double();
    private double SIGMIN = new double();
    private double SIGSM = new double();
    private boolean SILENT = new boolean();
    private boolean SINGUL = new boolean();
    private static final double SIX = 6;
    private double[] SLACK = new double[]();
    private double SMALL_W = new double();
    private double SMALLD = new double();
    private double SMALLW = new double();
    private int[] SORT = new int[]();
    private double SSTR = new double();
    private double STMAXL = new double();
    private double STPTRM = new double();
    private double TAU = new double();
    private double TAU0 = new double();
    private double TAUBND = new double();
    private double TAUFAC = new double();
    private double TAUMAX = new double();
    private double TAUQP = new double();
    private double TBND = new double();
    private boolean TE0 = new boolean();
    private boolean TE1 = new boolean();
    private boolean TE2 = new boolean();
    private boolean TE3 = new boolean();
    private double THETA = new double();
    private static final double THREE = 3;
    private static final double TM1 = 0.1;
    private static final double TM10 = 1E-10;
    private static final double TM2 = 0.01;
    private static final double TM3 = 0.001;
    private static final double TM4 = 0.0001;
    private static final double TM5 = 1E-05;
    private static final double TM6 = 1E-06;
    private static final double TM8 = 1E-08;
    private boolean[] TMPBUL = new boolean[]();
    private static final double TOL = 2.2204460492503131E-16;
    private double TOLMAC = new double();
    private static final double TP1 = 10;
    private static final double TP2 = 100;
    private static final double TP3 = 1000;
    private static final double TP4 = 10000;
    private static final double TWO = 2;
    private static final double TWOM2 = 0.25;
    private static final double TWOP11 = 2048;
    private static final double TWOP4 = 16;
    private double[] U = new double[]();
    private double[] U0 = new double[]();
    private double[] UD = new double[]();
    private double[] UD1 = new double[]();
    private double[] UG = new double[]();
    private double UPSI = new double();
    private double UPSI0 = new double();
    private double UPSI1 = new double();
    private double UPSIM = new double();
    private double UPSIST = new double();
    private boolean[] VAL = new boolean[]();
    private int[] VIOLIS = new int[]();
    private double[] W = new double[]();
    private double[] W1 = new double[]();
    private double[] WORK = new double[]();
    private double[] X = new double[]();
    private double[] X0 = new double[]();
    private double X0NORM = new double();
    private double[] X1 = new double[]();
    private double[] XGUESS = new double[]();
    private double[][] XJ = new double[][]();
    private double[] _lowerBounds = new double[]();
    private double[] XMIN = new double[]();
    private double XNORM = new double();
    private double[] XSC = new double[]();
    private double[] _scale = new double[]();
    private double[] XST = new double[]();
    private double[] XTR = new double[]();
    private double[] _upperBounds = new double[]();
    private double[] YU = new double[]();
    private static final double ZERO = 0;
    public interface IFunction   
    {
        double f(double[] x, int iact, boolean[] ierr) throws Exception ;
    
    }

    public interface IGradient   extends Provisdom.Optimization.DoNlp2.IFunction
    {
        void gradient(double[] x, int iact, double[] result) throws Exception ;
    
    }

}


