//
// Translated by CS2J (http://www.cs2j.com): 2/18/2016 3:06:48 PM
//

package java.Provisdom.Optimization;

import provisdom.LinearAlgebra.BLAS;
import provisdom.Optimization.NlpLinearConstraintSolver;

public class NlpLinearConstraintSolver   
{
    public NlpLinearConstraintSolver(Provisdom.Optimization.NlpLinearConstraintSolver.IFunction fcn, int nvar, int ncon, int neq, double[] a, double[] b, double[] lowerBound, double[] upperBound) throws Exception {
        this._objective = fcn;
        if (nvar <= 0)
        {
            Object[] objArray1 = new Object[]{ "nvar", nvar };
            ExceptionThrower.ThrowArgumentException("Numerical.Math", "NotPositive", objArray1);
        }
         
        this._var = nvar;
        if (ncon < 0)
        {
            Object[] objArray2 = new Object[]{ "ncon", ncon };
            ExceptionThrower.ThrowArgumentException("Numerical.Math", "Negative", objArray2);
        }
         
        this._con = ncon;
        if (neq < 0)
        {
            Object[] objArray3 = new Object[]{ "neq", neq };
            ExceptionThrower.ThrowArgumentException("Numerical.Math", "Negative", objArray3);
        }
         
        this._eq = neq;
        if (neq > ncon)
        {
            Object[] objArray4 = new Object[]{ "neq", "ncon" };
            ExceptionThrower.ThrowArgumentException("Numerical.Math", "XLessOrEqualY", objArray4);
        }
         
        if (a.Length != (ncon * nvar))
        {
            Object[] objArray5 = new Object[]{ "a", a.Length, ncon * nvar };
            ExceptionThrower.ThrowArgumentException("Numerical.Math", "NotEqual", objArray5);
        }
         
        this._a = new double[a.Length];
        a.CopyTo(this._a, 0);
        if (b.Length != ncon)
        {
            Object[] objArray6 = new Object[]{ "b", b.Length, ncon };
            ExceptionThrower.ThrowArgumentException("Numerical.Math", "NotEqual", objArray6);
        }
         
        this._b = new double[b.Length];
        b.CopyTo(this._b, 0);
        if (lowerBound.Length != nvar)
        {
            Object[] objArray7 = new Object[]{ "lowerBound", lowerBound.Length, nvar };
            ExceptionThrower.ThrowArgumentException("Numerical.Math", "NotEqual", objArray7);
        }
         
        this._lowerBounds = new double[lowerBound.Length];
        lowerBound.CopyTo(this._lowerBounds, 0);
        if (upperBound.Length != nvar)
        {
            Object[] objArray8 = new Object[]{ "upperBound", upperBound.Length, nvar };
            ExceptionThrower.ThrowArgumentException("Numerical.Math", "NotEqual", objArray8);
        }
         
        this._upperBounds = new double[upperBound.Length];
        upperBound.CopyTo(this._upperBounds, 0);
        for (int num1 = 0;num1 < nvar;num1++)
        {
            if (lowerBound[num1] > upperBound[num1])
            {
                Object[] objArray9 = new Object[]{ "lowerBound[" + num1 + "]", "upperBound[" + num1 + "]" };
                ExceptionThrower.ThrowArgumentException("Numerical.Math", "XLessOrEqualY", objArray9);
            }
             
        }
        this._iactUser = new int[(2 * nvar) + ncon];
        this._alamdaUser = new double[nvar];
        this._guess = new double[nvar];
        this._gradient = (fcn instanceof Provisdom.Optimization.NlpLinearConstraintSolver.IGradient) ? ((Provisdom.Optimization.NlpLinearConstraintSolver.IGradient)fcn) : null;
        this._maxObjective = 0x7fffffff;
        this._tol = Math.Sqrt(2.2204460492503131E-16);
    }

    public int[] getFinalActiveConstraints() throws Exception {
        if (this._iactUser == null)
        {
            return null;
        }
         
        return (int[])this._iactUser.Clone();
    }

    public double[] getLagrangeMultiplierEstimate() throws Exception {
        if (this._act > 0)
        {
            double[] numArray1 = new double[this._act];
            BLAS.copy(this._act, this._alamdaUser, numArray1);
            return numArray1;
        }
         
        return null;
    }

    public double[] getSolution() throws Exception {
        if (this._value == null)
        {
            return null;
        }
         
        return (double[])this._value.Clone();
    }

    private void l_l10nf(int n, int m, double[] a, int ia, double[] b, double[] xl, double[] xu, double[] x, int[] iact, int[] nact, int[] info, double[] w, int z, int u, int xbig, double relacc, double tol, int meql) throws Exception {
        int num6 = 0;
        double num12 = 0;
        if (nact[0] != 0)
        {
            for (int num7 = 1;num7 <= nact[0];num7++)
            {
                int num1;
                int num2;
                double num9;
                double num10;
                double num11;
                double num14 = new double();
                int num8 = num7 - 1;
                int num5 = iact[num8];
                if (num5 <= m)
                {
                    num9 = b[num5 - 1];
                    num10 = Math.Abs(b[num5 - 1]);
                    num11 = num10;
                    for (num1 = 1;num1 <= n;num1++)
                    {
                        num2 = num1 - 1;
                        double num15 = a[(num2 * ia) + (num5 - 1)];
                        num14 = num15 * x[num2];
                        num9 -= num14;
                        num10 += Math.Abs(num14);
                        num11 += Math.Abs(num15) * w[xbig + num2];
                    }
                }
                else
                {
                    num6 = num5 - m;
                    if (num6 <= n)
                    {
                        num9 = x[num6 - 1] - xl[num6 - 1];
                        num10 = Math.Abs(x[num6 - 1]) + Math.Abs(xl[num6 - 1]);
                        num11 = w[(xbig + num6) - 1] + Math.Abs(xl[num6 - 1]);
                        num12 = xl[num6 - 1];
                    }
                    else
                    {
                        num6 -= n;
                        num9 = xu[num6 - 1] - x[num6 - 1];
                        num10 = Math.Abs(x[num6 - 1]) + Math.Abs(xu[num6 - 1]);
                        num11 = w[(xbig + num6) - 1] + Math.Abs(xu[num6 - 1]);
                        num12 = xu[num6 - 1];
                    } 
                } 
                if (num9 != 0)
                {
                    num14 = num9 / num10;
                    if (num7 <= meql)
                    {
                        num14 = -Math.Abs(num14);
                    }
                     
                    if ((tol == relacc) || ((num14 + relacc) < 0))
                    {
                        info[0] = 1;
                        double num13 = num9 * w[u + num8];
                        int num4 = num7;
                        for (num1 = 1;num1 <= n;num1++)
                        {
                            double[] numArray1 = new double[]();
                            IntPtr ptr1 = new IntPtr();
                            num2 = num1 - 1;
                            (numArray1 = x)[(int)(ptr1 = (IntPtr)num2)] = numArray1[(int)ptr1] + (num13 * w[(z + num4) - 1]);
                            num4 += n;
                            w[xbig + num2] = Math.Max(w[xbig + num2], Math.Abs(x[num2]));
                        }
                        if (num5 > m)
                        {
                            x[num6 - 1] = num12;
                        }
                         
                    }
                    else if ((num9 / num11) > tol)
                    {
                        iact[num8] = -iact[num8];
                    }
                      
                }
                 
            }
            int num3 = nact[0];
            do
            {
                if (iact[num3 - 1] < 0)
                {
                    iact[num3 - 1] = -iact[num3 - 1];
                    this.l_l15nf(n, m, a, ia, iact, nact, w, z, u, relacc, num3);
                }
                 
                num3--;
            }
            while (num3 > meql);
        }
         
    }

    private void l_l10ng(int n, int m, double[] a, int ia, double[] b, double[] xl, double[] xu, double[] x, int[] iact, int[] nact, int[] info, double[] w, int z, int u, int xbig, double relacc, double tol, int meql) throws Exception {
        int num6 = 0;
        double num12 = 0;
        if (nact[0] != 0)
        {
            for (int num7 = 1;num7 <= nact[0];num7++)
            {
                int num1;
                int num2;
                double num9;
                double num10;
                double num11;
                double num14 = new double();
                int num8 = num7 - 1;
                int num5 = iact[num8];
                if (num5 <= m)
                {
                    num9 = b[num5 - 1];
                    num10 = Math.Abs(b[num5 - 1]);
                    num11 = num10;
                    for (num1 = 1;num1 <= n;num1++)
                    {
                        num2 = num1 - 1;
                        double num15 = a[(num2 * ia) + (num5 - 1)];
                        num14 = num15 * x[num2];
                        num9 -= num14;
                        num10 += Math.Abs(num14);
                        num11 += Math.Abs(num15) * w[xbig + num2];
                    }
                }
                else
                {
                    num6 = num5 - m;
                    if (num6 <= n)
                    {
                        num9 = x[num6 - 1] - xl[num6 - 1];
                        num10 = Math.Abs(x[num6 - 1]) + Math.Abs(xl[num6 - 1]);
                        num11 = w[(xbig + num6) - 1] + Math.Abs(xl[num6 - 1]);
                        num12 = xl[num6 - 1];
                    }
                    else
                    {
                        num6 -= n;
                        num9 = xu[num6 - 1] - x[num6 - 1];
                        num10 = Math.Abs(x[num6 - 1]) + Math.Abs(xu[num6 - 1]);
                        num11 = w[(xbig + num6) - 1] + Math.Abs(xu[num6 - 1]);
                        num12 = xu[num6 - 1];
                    } 
                } 
                if (num9 != 0)
                {
                    num14 = num9 / num10;
                    if (num7 <= meql)
                    {
                        num14 = -Math.Abs(num14);
                    }
                     
                    if ((tol == relacc) || ((num14 + relacc) < 0))
                    {
                        info[0] = 1;
                        double num13 = num9 * w[u + num8];
                        int num4 = num7;
                        for (num1 = 1;num1 <= n;num1++)
                        {
                            double[] numArray1 = new double[]();
                            IntPtr ptr1 = new IntPtr();
                            num2 = num1 - 1;
                            (numArray1 = x)[(int)(ptr1 = (IntPtr)num2)] = numArray1[(int)ptr1] + (num13 * w[(z + num4) - 1]);
                            num4 += n;
                            w[xbig + num2] = Math.Max(w[xbig + num2], Math.Abs(x[num2]));
                        }
                        if (num5 > m)
                        {
                            x[num6 - 1] = num12;
                        }
                         
                    }
                    else if ((num9 / num11) > tol)
                    {
                        iact[num8] = -iact[num8];
                    }
                      
                }
                 
            }
            int num3 = nact[0];
            do
            {
                if (iact[num3 - 1] < 0)
                {
                    iact[num3 - 1] = -iact[num3 - 1];
                    this.l_l15ng(n, m, a, ia, iact, nact, w, z, u, relacc, num3);
                }
                 
                num3--;
            }
            while (num3 > meql);
        }
         
    }

    private void l_l11nf(int n, int m, double[] a, int ia, double[] b, double[] xl, double[] xu, double[] x, int[] iact, int[] nact, double[] par, double[] w, int g, int z, int u, int xbig, int bres, int d, int ztg, double relacc, double tol, double[] stepcb, double[] sumres, int meql, int[] msat, int mtot, int[] indxbd, int gm, int gmnew, int parnew, int cgrad) throws Exception {
        int num1;
        int num2;
        double[] numArray2 = new double[1];
        double[] numArray1 = numArray2;
        int num3 = mtot - msat[0];
        if (num3 > 0)
        {
            for (num1 = 1;num1 <= n;num1++)
            {
                num2 = num1 - 1;
                w[g + num2] = 0;
            }
            sumres[0] = 0;
        }
         
        int num10 = msat[0];
        int num9 = nact[0];
        msat[0] = nact[0];
        int num8 = meql + 1;
        for (int num6 = num8;num6 <= mtot;num6++)
        {
            int num5;
            double num12 = new double();
            double num13 = new double();
            int num7 = num6 - 1;
            int num4 = iact[num7];
            if (num4 <= m)
            {
                num12 = b[num4 - 1];
                num13 = Math.Abs(b[num4 - 1]);
                for (num1 = 1;num1 <= n;num1++)
                {
                    num2 = num1 - 1;
                    num12 += -x[num2] * a[(num2 * ia) + (num4 - 1)];
                    num13 += Math.Abs((double)(w[xbig + num2] * a[(num2 * ia) + (num4 - 1)]));
                }
            }
            else
            {
                num5 = num4 - m;
                if (num5 <= n)
                {
                    num12 = x[num5 - 1] - xl[num5 - 1];
                    num13 = Math.Abs(w[(xbig + num5) - 1]) + Math.Abs(xl[num5 - 1]);
                }
                else
                {
                    num5 -= n;
                    num12 = xu[num5 - 1] - x[num5 - 1];
                    num13 = Math.Abs(w[(xbig + num5) - 1]) + Math.Abs(xu[num5 - 1]);
                } 
            } 
            w[(bres + num4) - 1] = num12;
            double num15 = 0;
            if (num13 != 0)
            {
                num15 = num12 / num13;
            }
             
            if (((num6 > num10) && (num15 < 0)) && ((num15 + relacc) >= 0))
            {
                double num14 = new double();
                if (num4 <= m)
                {
                    num14 = Math.Abs(b[num4 - 1]);
                    for (num1 = 1;num1 <= n;num1++)
                    {
                        num2 = num1 - 1;
                        num14 += Math.Abs((double)(x[num2] * a[(num2 * ia) + (num4 - 1)]));
                    }
                }
                else
                {
                    num5 = num4 - m;
                    if (num5 <= n)
                    {
                        num14 = Math.Abs(x[num5 - 1]) + Math.Abs(xl[num5 - 1]);
                    }
                    else
                    {
                        num14 = Math.Abs(x[(num5 - n) - 1]) + Math.Abs(xu[(num5 - n) - 1]);
                    } 
                } 
                if (Math.Abs(num12) <= (num14 * relacc))
                {
                    num15 = 0;
                }
                 
            }
             
            int num11 = 0;
            if (num6 <= nact[0])
            {
                num11 = 1;
            }
             
            if (num11 == 0)
            {
                if ((num6 <= num10) || (num15 >= 0))
                {
                    int[] numArray3 = new int[]();
                    (numArray3 = msat)[0] = numArray3[0] + 1;
                    if (msat[0] < num6)
                    {
                        iact[num7] = iact[msat[0] - 1];
                    }
                     
                    if (num15 > tol)
                    {
                        iact[msat[0] - 1] = num4;
                    }
                    else
                    {
                        num9++;
                        iact[msat[0] - 1] = iact[num9 - 1];
                        iact[num9 - 1] = num4;
                    } 
                }
                else
                {
                    IntPtr ptr1 = new IntPtr();
                    if (num4 <= m)
                    {
                        for (num1 = 1;num1 <= n;num1++)
                        {
                            num2 = num1 - 1;
                            (numArray2 = w)[(int)(ptr1 = (IntPtr)(g + num2))] = numArray2[(int)ptr1] + a[(num2 * ia) + (num4 - 1)];
                        }
                    }
                    else
                    {
                        num4 -= m;
                        if (num4 <= n)
                        {
                            (numArray2 = w)[(int)(ptr1 = (IntPtr)((g + num4) - 1))] = numArray2[(int)ptr1] - 1;
                        }
                        else
                        {
                            (numArray2 = w)[(int)(ptr1 = (IntPtr)(((g + num4) - n) - 1))] = numArray2[(int)ptr1] + 1;
                        } 
                    } 
                    (numArray2 = sumres)[0] = numArray2[0] + Math.Abs(num12);
                } 
            }
             
        }
        stepcb[0] = 0;
        if ((num3 <= 0) || (msat[0] != mtot))
        {
            this.l_l16nf(n, m, a, ia, iact, nact, par, w, g, z, u, d, ztg, relacc, numArray1, meql, num9, gm, gmnew, parnew, cgrad);
            if (numArray1[0] < 0)
            {
                this.l_l17nf(n, m, a, ia, iact, w, bres, d, stepcb, numArray1, num9, msat, mtot, indxbd);
            }
             
            if (num3 == 0)
            {
                sumres[0] = numArray1[0];
            }
             
        }
         
    }

    private void l_l11ng(int n, int m, double[] a, int ia, double[] b, double[] xl, double[] xu, double[] x, int[] iact, int[] nact, double[] par, double[] w, int g, int z, int u, int xbig, int bres, int d, int ztg, double relacc, double tol, double[] stepcb, double[] sumres, int meql, int[] msat, int mtot, int[] indxbd, int gm, int gmnew, int parnew, int cgrad) throws Exception {
        int num1;
        int num2;
        double[] numArray2 = new double[1];
        double[] numArray1 = numArray2;
        int num3 = mtot - msat[0];
        if (num3 > 0)
        {
            for (num1 = 1;num1 <= n;num1++)
            {
                num2 = num1 - 1;
                w[g + num2] = 0;
            }
            sumres[0] = 0;
        }
         
        int num10 = msat[0];
        int num9 = nact[0];
        msat[0] = nact[0];
        int num8 = meql + 1;
        for (int num6 = num8;num6 <= mtot;num6++)
        {
            int num5;
            double num12 = new double();
            double num13 = new double();
            int num7 = num6 - 1;
            int num4 = iact[num7];
            if (num4 <= m)
            {
                num12 = b[num4 - 1];
                num13 = Math.Abs(b[num4 - 1]);
                for (num1 = 1;num1 <= n;num1++)
                {
                    num2 = num1 - 1;
                    num12 += -x[num2] * a[(num2 * ia) + (num4 - 1)];
                    num13 += Math.Abs((double)(w[xbig + num2] * a[(num2 * ia) + (num4 - 1)]));
                }
            }
            else
            {
                num5 = num4 - m;
                if (num5 <= n)
                {
                    num12 = x[num5 - 1] - xl[num5 - 1];
                    num13 = Math.Abs(w[(xbig + num5) - 1]) + Math.Abs(xl[num5 - 1]);
                }
                else
                {
                    num5 -= n;
                    num12 = xu[num5 - 1] - x[num5 - 1];
                    num13 = Math.Abs(w[(xbig + num5) - 1]) + Math.Abs(xu[num5 - 1]);
                } 
            } 
            w[(bres + num4) - 1] = num12;
            double num15 = 0;
            if (num13 != 0)
            {
                num15 = num12 / num13;
            }
             
            if (((num6 > num10) && (num15 < 0)) && ((num15 + relacc) >= 0))
            {
                double num14 = new double();
                if (num4 <= m)
                {
                    num14 = Math.Abs(b[num4 - 1]);
                    for (num1 = 1;num1 <= n;num1++)
                    {
                        num2 = num1 - 1;
                        num14 += Math.Abs((double)(x[num2] * a[(num2 * ia) + (num4 - 1)]));
                    }
                }
                else
                {
                    num5 = num4 - m;
                    if (num5 <= n)
                    {
                        num14 = Math.Abs(x[num5 - 1]) + Math.Abs(xl[num5 - 1]);
                    }
                    else
                    {
                        num14 = Math.Abs(x[(num5 - n) - 1]) + Math.Abs(xu[(num5 - n) - 1]);
                    } 
                } 
                if (Math.Abs(num12) <= (num14 * relacc))
                {
                    num15 = 0;
                }
                 
            }
             
            int num11 = 0;
            if (num6 <= nact[0])
            {
                num11 = 1;
            }
             
            if (num11 == 0)
            {
                if ((num6 <= num10) || (num15 >= 0))
                {
                    int[] numArray3 = new int[]();
                    (numArray3 = msat)[0] = numArray3[0] + 1;
                    if (msat[0] < num6)
                    {
                        iact[num7] = iact[msat[0] - 1];
                    }
                     
                    if (num15 > tol)
                    {
                        iact[msat[0] - 1] = num4;
                    }
                    else
                    {
                        num9++;
                        iact[msat[0] - 1] = iact[num9 - 1];
                        iact[num9 - 1] = num4;
                    } 
                }
                else
                {
                    IntPtr ptr1 = new IntPtr();
                    if (num4 <= m)
                    {
                        for (num1 = 1;num1 <= n;num1++)
                        {
                            num2 = num1 - 1;
                            (numArray2 = w)[(int)(ptr1 = (IntPtr)(g + num2))] = numArray2[(int)ptr1] + a[(num2 * ia) + (num4 - 1)];
                        }
                    }
                    else
                    {
                        num4 -= m;
                        if (num4 <= n)
                        {
                            (numArray2 = w)[(int)(ptr1 = (IntPtr)((g + num4) - 1))] = numArray2[(int)ptr1] - 1;
                        }
                        else
                        {
                            (numArray2 = w)[(int)(ptr1 = (IntPtr)(((g + num4) - n) - 1))] = numArray2[(int)ptr1] + 1;
                        } 
                    } 
                    (numArray2 = sumres)[0] = numArray2[0] + Math.Abs(num12);
                } 
            }
             
        }
        stepcb[0] = 0;
        if ((num3 <= 0) || (msat[0] != mtot))
        {
            this.l_l16ng(n, m, a, ia, iact, nact, par, w, g, z, u, d, ztg, relacc, numArray1, meql, num9, gm, gmnew, parnew, cgrad);
            if (numArray1[0] < 0)
            {
                this.l_l17ng(n, m, a, ia, iact, w, bres, d, stepcb, numArray1, num9, msat, mtot, indxbd);
            }
             
            if (num3 == 0)
            {
                sumres[0] = numArray1[0];
            }
             
        }
         
    }

    private void l_l12nf(int n, int m, double[] a, int ia, int[] iact, int nact, double[] par, double[] w, int g, int reskt, int z, int u, int bres, double[] relaxf, int meql, double[] ssqkt, int parw, int resktw) throws Exception {
        int num2;
        int num5;
        int num7;
        int num8 = new int();
        int num9;
        double[] numArray1 = new double[]();
        double num12 = 0;
        int num1 = 1;
        while (num1 <= n)
        {
            num2 = num1 - 1;
            w[reskt + num2] = w[g + num2];
            num1++;
        }
        if (nact <= 0)
        {
            ssqkt[0] = 0;
            for (num1 = 1;num1 <= n;num1++)
            {
                num2 = num1 - 1;
                (numArray1 = ssqkt)[0] = numArray1[0] + Math.Pow(w[g + num2], 2);
            }
            goto Label_02C0
        }
         
        int num3 = 0;
        Label_0038:num9 = 1;
        while (num9 <= nact)
        {
            IntPtr ptr1 = new IntPtr();
            int num10 = num9 - 1;
            num7 = (nact + 1) - num9;
            num5 = iact[num7 - 1];
            double num13 = 0;
            int num4 = num7;
            num1 = 1;
            while (num1 <= n)
            {
                num2 = num1 - 1;
                num13 += w[(z + num4) - 1] * w[reskt + num2];
                num4 += n;
                num1++;
            }
            num13 *= w[(u + num7) - 1];
            if (num3 == 0)
            {
                par[num7 - 1] = 0;
            }
             
            if ((num7 <= meql) || ((par[num7 - 1] + num13) < 0))
            {
                (numArray1 = par)[(int)(ptr1 = (IntPtr)(num7 - 1))] = numArray1[(int)ptr1] + num13;
            }
            else
            {
                num13 = -par[num7 - 1];
                par[num7 - 1] = 0;
            } 
            if (num13 != 0)
            {
                if (num5 <= m)
                {
                    for (num1 = 1;num1 <= n;num1++)
                    {
                        num2 = num1 - 1;
                        (numArray1 = w)[(int)(ptr1 = (IntPtr)(reskt + num2))] = numArray1[(int)ptr1] + (-num13 * a[(num2 * ia) + (num5 - 1)]);
                    }
                }
                else
                {
                    int num6 = num5 - m;
                    if (num6 <= n)
                    {
                        (numArray1 = w)[(int)(ptr1 = (IntPtr)((reskt + num6) - 1))] = numArray1[(int)ptr1] + num13;
                    }
                    else
                    {
                        (numArray1 = w)[(int)(ptr1 = (IntPtr)(((reskt + num6) - n) - 1))] = numArray1[(int)ptr1] - num13;
                    } 
                } 
            }
             
            num9++;
        }
        ssqkt[0] = 0;
        if (nact == n)
        {
            return ;
        }
         
        num1 = 1;
        while (num1 <= n)
        {
            num2 = num1 - 1;
            (numArray1 = ssqkt)[0] = numArray1[0] + Math.Pow(w[reskt + num2], 2);
            num1++;
        }
        if (num3 == 0)
        {
            num3 = 1;
            for (num7 = 1;num7 <= nact;num7++)
            {
                num8 = num7 - 1;
                w[parw + num8] = par[num8];
            }
            for (num1 = 1;num1 <= n;num1++)
            {
                num2 = num1 - 1;
                w[resktw + num2] = w[reskt + num2];
            }
            num12 = ssqkt[0];
            goto Label_0038
        }
         
        if (num12 < ssqkt[0])
        {
            for (num7 = 1;num7 <= nact;num7++)
            {
                num8 = num7 - 1;
                par[num8] = w[parw + num8];
            }
            for (num1 = 1;num1 <= n;num1++)
            {
                num2 = num1 - 1;
                w[reskt + num2] = w[resktw + num2];
            }
            ssqkt[0] = num12;
        }
         
        Label_02C0:relaxf[0] = 0;
        if (meql < nact)
        {
            int num11 = meql + 1;
            for (num7 = num11;num7 <= nact;num7++)
            {
                num8 = num7 - 1;
                num5 = iact[num8];
                if (w[(bres + num5) - 1] > 0)
                {
                    (numArray1 = relaxf)[0] = numArray1[0] + (-par[num8] * w[(bres + num5) - 1]);
                }
                 
            }
        }
         
    }

    private void l_l12ng(int n, int m, double[] a, int ia, int[] iact, int nact, double[] par, double[] w, int g, int reskt, int z, int u, int bres, double[] relaxf, int meql, double[] ssqkt, int parw, int resktw) throws Exception {
        int num2;
        int num5;
        int num7;
        int num8 = new int();
        int num9;
        double[] numArray1 = new double[]();
        double num12 = 0;
        int num1 = 1;
        while (num1 <= n)
        {
            num2 = num1 - 1;
            w[reskt + num2] = w[g + num2];
            num1++;
        }
        if (nact <= 0)
        {
            ssqkt[0] = 0;
            for (num1 = 1;num1 <= n;num1++)
            {
                num2 = num1 - 1;
                (numArray1 = ssqkt)[0] = numArray1[0] + Math.Pow(w[g + num2], 2);
            }
            goto Label_02C0
        }
         
        int num3 = 0;
        Label_0038:num9 = 1;
        while (num9 <= nact)
        {
            IntPtr ptr1 = new IntPtr();
            int num10 = num9 - 1;
            num7 = (nact + 1) - num9;
            num5 = iact[num7 - 1];
            double num13 = 0;
            int num4 = num7;
            num1 = 1;
            while (num1 <= n)
            {
                num2 = num1 - 1;
                num13 += w[(z + num4) - 1] * w[reskt + num2];
                num4 += n;
                num1++;
            }
            num13 *= w[(u + num7) - 1];
            if (num3 == 0)
            {
                par[num7 - 1] = 0;
            }
             
            if ((num7 <= meql) || ((par[num7 - 1] + num13) < 0))
            {
                (numArray1 = par)[(int)(ptr1 = (IntPtr)(num7 - 1))] = numArray1[(int)ptr1] + num13;
            }
            else
            {
                num13 = -par[num7 - 1];
                par[num7 - 1] = 0;
            } 
            if (num13 != 0)
            {
                if (num5 <= m)
                {
                    for (num1 = 1;num1 <= n;num1++)
                    {
                        num2 = num1 - 1;
                        (numArray1 = w)[(int)(ptr1 = (IntPtr)(reskt + num2))] = numArray1[(int)ptr1] + (-num13 * a[(num2 * ia) + (num5 - 1)]);
                    }
                }
                else
                {
                    int num6 = num5 - m;
                    if (num6 <= n)
                    {
                        (numArray1 = w)[(int)(ptr1 = (IntPtr)((reskt + num6) - 1))] = numArray1[(int)ptr1] + num13;
                    }
                    else
                    {
                        (numArray1 = w)[(int)(ptr1 = (IntPtr)(((reskt + num6) - n) - 1))] = numArray1[(int)ptr1] - num13;
                    } 
                } 
            }
             
            num9++;
        }
        ssqkt[0] = 0;
        if (nact == n)
        {
            return ;
        }
         
        num1 = 1;
        while (num1 <= n)
        {
            num2 = num1 - 1;
            (numArray1 = ssqkt)[0] = numArray1[0] + Math.Pow(w[reskt + num2], 2);
            num1++;
        }
        if (num3 == 0)
        {
            num3 = 1;
            for (num7 = 1;num7 <= nact;num7++)
            {
                num8 = num7 - 1;
                w[parw + num8] = par[num8];
            }
            for (num1 = 1;num1 <= n;num1++)
            {
                num2 = num1 - 1;
                w[resktw + num2] = w[reskt + num2];
            }
            num12 = ssqkt[0];
            goto Label_0038
        }
         
        if (num12 < ssqkt[0])
        {
            for (num7 = 1;num7 <= nact;num7++)
            {
                num8 = num7 - 1;
                par[num8] = w[parw + num8];
            }
            for (num1 = 1;num1 <= n;num1++)
            {
                num2 = num1 - 1;
                w[reskt + num2] = w[resktw + num2];
            }
            ssqkt[0] = num12;
        }
         
        Label_02C0:relaxf[0] = 0;
        if (meql < nact)
        {
            int num11 = meql + 1;
            for (num7 = num11;num7 <= nact;num7++)
            {
                num8 = num7 - 1;
                num5 = iact[num8];
                if (w[(bres + num5) - 1] > 0)
                {
                    (numArray1 = relaxf)[0] = numArray1[0] + (-par[num8] * w[(bres + num5) - 1]);
                }
                 
            }
        }
         
    }

    private void l_l13nf(Provisdom.Optimization.NlpLinearConstraintSolver.IFunction fcn, int n, double[] x, double[] w, int g, int d, int xs, int gs, double relacc, double stepcb, double ddotg, double[] f, double[] step, int[] nfvals, int nfmax, int gopt, double[] obj) throws Exception {
        int num2;
        double num21 = new double();
        int[] numArray2 = new int[1];
        int[] numArray1 = numArray2;
        double num11 = 0;
        double num5 = 0;
        double num15 = 0.9;
        numArray1[0] = 0;
        double num14 = -1;
        int num1 = 1;
        while (num1 <= n)
        {
            num2 = num1 - 1;
            w[xs + num2] = x[num2];
            w[gs + num2] = w[g + num2];
            w[gopt + num2] = w[g + num2];
            if (w[d + num2] != 0)
            {
                num21 = Math.Abs((double)(x[num2] / w[d + num2]));
                if ((num14 < 0) || (num21 < num14))
                {
                    num14 = num21;
                }
                 
            }
             
            num1++;
        }
        step[0] = Math.Min(1, stepcb);
        double num19 = Math.Max((double)(relacc * num14), (double)(1E-12 * step[0]));
        step[0] = Math.Max(num19, step[0]);
        double num16 = 0;
        double num10 = f[0];
        double num4 = ddotg;
        double num18 = 0;
        double num12 = f[0];
        double num7 = ddotg;
        double num17 = 0;
        double num20 = 0;
        double num13 = f[0];
        double num9 = Math.Abs(ddotg);
        Label_013D:num1 = 1;
        while (num1 <= n)
        {
            num2 = num1 - 1;
            x[num2] = w[xs + num2] + (step[0] * w[d + num2]);
            num1++;
        }
        f[0] = fcn.F(x);
        this.l_l21nf(fcn, n, x, f[0], w, g, numArray1);
        (numArray2 = numArray1)[0] = numArray2[0] + 1;
        double num8 = 0;
        num1 = 1;
        while (num1 <= n)
        {
            num2 = num1 - 1;
            num8 += w[d + num2] * w[g + num2];
            num1++;
        }
        if ((f[0] <= num13) && ((f[0] < num13) || (Math.Abs(num8) < num9)))
        {
            num20 = step[0];
            num13 = f[0];
            obj[0] = f[0];
            for (num1 = 1;num1 <= n;num1++)
            {
                num2 = num1 - 1;
                w[gopt + num2] = w[g + num2];
            }
            num9 = Math.Abs(num8);
        }
         
        if ((nfvals[0] + numArray1[0]) == nfmax)
        {
            if (step[0] != num20)
            {
                step[0] = num20;
                f[0] = num13;
                for (num1 = 1;num1 <= n;num1++)
                {
                    num2 = num1 - 1;
                    x[num2] = w[xs + num2] + (step[0] * w[d + num2]);
                    w[g + num2] = w[gopt + num2];
                }
            }
             
            (numArray2 = nfvals)[0] = numArray2[0] + numArray1[0];
        }
        else
        {
            int num3 = 0;
            if (f[0] >= (num10 + ((0.1 * (step[0] - num16)) * num4)))
            {
                if (((num17 > 0) || (f[0] > num10)) || (num8 > (0.5 * ddotg)))
                {
                    num17 = step[0];
                    num11 = f[0];
                    num5 = num8;
                    num3 = 1;
                }
                 
                if (num3 == 0)
                {
                    num16 = step[0];
                    num10 = f[0];
                    num4 = num8;
                }
                 
            }
             
            if (num3 == 0)
            {
                if (num8 >= (0.7 * num4))
                {
                    if (step[0] != num20)
                    {
                        step[0] = num20;
                        f[0] = num13;
                        for (num1 = 1;num1 <= n;num1++)
                        {
                            num2 = num1 - 1;
                            x[num2] = w[xs + num2] + (step[0] * w[d + num2]);
                            w[g + num2] = w[gopt + num2];
                        }
                    }
                     
                    (numArray2 = nfvals)[0] = numArray2[0] + numArray1[0];
                    return ;
                }
                 
                num18 = step[0];
                num12 = f[0];
                num7 = num8;
            }
             
            if ((num17 > 0) && (num18 >= (num15 * num17)))
            {
                if (step[0] != num20)
                {
                    step[0] = num20;
                    f[0] = num13;
                    for (num1 = 1;num1 <= n;num1++)
                    {
                        num2 = num1 - 1;
                        x[num2] = w[xs + num2] + (step[0] * w[d + num2]);
                        w[g + num2] = w[gopt + num2];
                    }
                }
                 
                (numArray2 = nfvals)[0] = numArray2[0] + numArray1[0];
            }
            else
            {
                double[] numArray3 = new double[]();
                if (num17 == 0)
                {
                    if (step[0] == stepcb)
                    {
                        if (step[0] != num20)
                        {
                            step[0] = num20;
                            f[0] = num13;
                            for (num1 = 1;num1 <= n;num1++)
                            {
                                num2 = num1 - 1;
                                x[num2] = w[xs + num2] + (step[0] * w[d + num2]);
                                w[g + num2] = w[gopt + num2];
                            }
                        }
                         
                        (numArray2 = nfvals)[0] = numArray2[0] + numArray1[0];
                        return ;
                    }
                     
                    num21 = 10;
                    if (num8 > (0.9 * ddotg))
                    {
                        num21 = ddotg / (ddotg - num8);
                    }
                     
                    step[0] = Math.Min(num21 * step[0], stepcb);
                    goto Label_013D
                }
                 
                if ((numArray1[0] == 1) || (num18 > 0))
                {
                    double num6 = ((2 * (num11 - num12)) / (num17 - num18)) - (0.5 * (num7 + num5));
                    if (num6 >= 0)
                    {
                        num14 = Math.Max((double)0.1, (double)((0.5 * num7) / (num7 - num6)));
                    }
                    else
                    {
                        num14 = ((0.5 * num5) - num6) / (num5 - num6);
                    } 
                    step[0] = num18 + (num14 * (num17 - num18));
                    goto Label_013D
                }
                 
                (numArray3 = step)[0] = numArray3[0] * 0.1;
                if (step[0] >= num19)
                {
                    goto Label_013D
                }
                 
                if (step[0] != num20)
                {
                    step[0] = num20;
                    f[0] = num13;
                    for (num1 = 1;num1 <= n;num1++)
                    {
                        num2 = num1 - 1;
                        x[num2] = w[xs + num2] + (step[0] * w[d + num2]);
                        w[g + num2] = w[gopt + num2];
                    }
                }
                 
                (numArray2 = nfvals)[0] = numArray2[0] + numArray1[0];
            } 
        } 
    }

    private void l_l13ng(Provisdom.Optimization.NlpLinearConstraintSolver.IFunction fcn, Provisdom.Optimization.NlpLinearConstraintSolver.IGradient grad, int n, double[] x, double[] w, int g, int d, int xs, int gs, double relacc, double stepcb, double ddotg, double[] f, double[] step, int[] nfvals, int nfmax, int gopt, double[] obj) throws Exception {
        int num2;
        double num22 = new double();
        int[] numArray3 = new int[]();
        double num12 = 0;
        double num6 = 0;
        double num16 = 0.9;
        int num3 = 0;
        double num15 = -1;
        int num1 = 1;
        while (num1 <= n)
        {
            num2 = num1 - 1;
            w[xs + num2] = x[num2];
            w[gs + num2] = w[g + num2];
            w[gopt + num2] = w[g + num2];
            if (w[d + num2] != 0)
            {
                num22 = Math.Abs((double)(x[num2] / w[d + num2]));
                if ((num15 < 0) || (num22 < num15))
                {
                    num15 = num22;
                }
                 
            }
             
            num1++;
        }
        step[0] = Math.Min(1, stepcb);
        double num20 = Math.Max((double)(relacc * num15), (double)(1E-12 * step[0]));
        step[0] = Math.Max(num20, step[0]);
        double num17 = 0;
        double num11 = f[0];
        double num5 = ddotg;
        double num19 = 0;
        double num13 = f[0];
        double num8 = ddotg;
        double num18 = 0;
        double num21 = 0;
        double num14 = f[0];
        double num10 = Math.Abs(ddotg);
        Label_0132:num1 = 1;
        while (num1 <= n)
        {
            num2 = num1 - 1;
            x[num2] = w[xs + num2] + (step[0] * w[d + num2]);
            num1++;
        }
        f[0] = fcn.F(x);
        double[] numArray1 = new double[n];
        BLAS.copy(n, w, g, 1, numArray1, 0, 1);
        grad.Gradient(x, numArray1);
        BLAS.copy(n, numArray1, 0, 1, w, g, 1);
        num3++;
        double num9 = 0;
        num1 = 1;
        while (num1 <= n)
        {
            num2 = num1 - 1;
            num9 += w[d + num2] * w[g + num2];
            num1++;
        }
        if ((f[0] <= num14) && ((f[0] < num14) || (Math.Abs(num9) < num10)))
        {
            num21 = step[0];
            num14 = f[0];
            obj[0] = f[0];
            for (num1 = 1;num1 <= n;num1++)
            {
                num2 = num1 - 1;
                w[gopt + num2] = w[g + num2];
            }
            num10 = Math.Abs(num9);
        }
         
        if ((nfvals[0] + num3) != nfmax)
        {
            int num4 = 0;
            if (f[0] >= (num11 + ((0.1 * (step[0] - num17)) * num5)))
            {
                if (((num18 > 0) || (f[0] > num11)) || (num9 > (0.5 * ddotg)))
                {
                    num18 = step[0];
                    num12 = f[0];
                    num6 = num9;
                    num4 = 1;
                }
                 
                if (num4 == 0)
                {
                    num17 = step[0];
                    num11 = f[0];
                    num5 = num9;
                }
                 
            }
             
            if (num4 == 0)
            {
                if (num9 >= (0.7 * num5))
                {
                    goto Label_03EA
                }
                 
                num19 = step[0];
                num13 = f[0];
                num8 = num9;
            }
             
            if ((num18 <= 0) || (num19 < (num16 * num18)))
            {
                double[] numArray2 = new double[]();
                if (num18 == 0)
                {
                    if (step[0] == stepcb)
                    {
                        goto Label_03EA
                    }
                     
                    num22 = 10;
                    if (num9 > (0.9 * ddotg))
                    {
                        num22 = ddotg / (ddotg - num9);
                    }
                     
                    step[0] = Math.Min(num22 * step[0], stepcb);
                    goto Label_0132
                }
                 
                if ((num3 == 1) || (num19 > 0))
                {
                    double num7 = ((2 * (num12 - num13)) / (num18 - num19)) - (0.5 * (num8 + num6));
                    if (num7 >= 0)
                    {
                        num15 = Math.Max((double)0.1, (double)((0.5 * num8) / (num8 - num7)));
                    }
                    else
                    {
                        num15 = ((0.5 * num6) - num7) / (num6 - num7);
                    } 
                    step[0] = num19 + (num15 * (num18 - num19));
                    goto Label_0132
                }
                 
                (numArray2 = step)[0] = numArray2[0] * 0.1;
                if (step[0] >= num20)
                {
                    goto Label_0132
                }
                 
            }
             
        }
         
        Label_03EA:if (step[0] != num21)
        {
            step[0] = num21;
            f[0] = num14;
            for (num1 = 1;num1 <= n;num1++)
            {
                num2 = num1 - 1;
                x[num2] = w[xs + num2] + (step[0] * w[d + num2]);
                w[g + num2] = w[gopt + num2];
            }
        }
         
        (numArray3 = nfvals)[0] = numArray3[0] + num3;
    }

    private void l_l14nf(int n, double[] x, int nact, double[] w, int g, int z, int ztg, int xs, int gs, double[] zznorm) throws Exception {
        int num2;
        int num3;
        int num7;
        double num9 = 0;
        double num10 = 0;
        double num12 = 0;
        int num1 = 1;
        while (num1 <= n)
        {
            num2 = num1 - 1;
            w[xs + num2] = x[num2] - w[xs + num2];
            num9 += Math.Pow(w[xs + num2], 2);
            num12 += w[gs + num2] * w[xs + num2];
            w[gs + num2] = w[g + num2] - w[gs + num2];
            num10 += w[gs + num2] * w[xs + num2];
            num1++;
        }
        if (num10 < (0.1 * Math.Abs(num12)))
        {
            return ;
        }
         
        int num4 = n;
        Label_00BC:num7 = num4;
        num4--;
        if (num4 > nact)
        {
            if (w[(ztg + num7) - 1] != 0)
            {
                num12 = Math.Abs(w[(ztg + num7) - 1]) * Math.Sqrt(1 + Math.Pow(w[(ztg + num4) - 1] / w[(ztg + num7) - 1], 2));
                double num13 = w[(ztg + num4) - 1] / num12;
                double num14 = w[(ztg + num7) - 1] / num12;
                w[(ztg + num4) - 1] = num12;
                num3 = num4;
                for (num1 = 1;num1 <= n;num1++)
                {
                    num2 = num1 - 1;
                    num12 = (num13 * w[z + num3]) - (num14 * w[(z + num3) - 1]);
                    w[(z + num3) - 1] = (num13 * w[(z + num3) - 1]) + (num14 * w[z + num3]);
                    w[z + num3] = num12;
                    num3 += n;
                }
            }
             
            goto Label_00BC
        }
         
        if (zznorm[0] < 0)
        {
            zznorm[0] = num9 / num10;
        }
        else
        {
            num12 = Math.Sqrt((zznorm[0] * num9) / num10);
            zznorm[0] = Math.Min(zznorm[0], num12);
            zznorm[0] = Math.Max(zznorm[0], 0.1 * num12);
        } 
        int num8 = nact + 1;
        num12 = Math.Sqrt(num10);
        num3 = num8;
        num1 = 1;
        while (num1 <= n)
        {
            num2 = num1 - 1;
            w[(z + num3) - 1] = w[xs + num2] / num12;
            num3 += n;
            num1++;
        }
        if (num8 < n)
        {
            int num6 = num8 + 1;
            for (num4 = num6;num4 <= n;num4++)
            {
                double[] numArray1 = new double[]();
                IntPtr ptr1 = new IntPtr();
                int num5 = num4 - 1;
                num12 = 0;
                num3 = num4;
                num1 = 1;
                while (num1 <= n)
                {
                    num2 = num1 - 1;
                    num12 += w[gs + num2] * w[(z + num3) - 1];
                    num3 += n;
                    num1++;
                }
                num12 /= num10;
                double num11 = 0;
                num3 = num4;
                num1 = 1;
                while (num1 <= n)
                {
                    num2 = num1 - 1;
                    (numArray1 = w)[(int)(ptr1 = (IntPtr)((z + num3) - 1))] = numArray1[(int)ptr1] + (-num12 * w[xs + num2]);
                    num11 += Math.Pow(w[(z + num3) - 1], 2);
                    num3 += n;
                    num1++;
                }
                if (num11 < zznorm[0])
                {
                    num12 = Math.Sqrt(zznorm[0] / num11);
                    num3 = num4;
                    for (num1 = 1;num1 <= n;num1++)
                    {
                        num2 = num1 - 1;
                        (numArray1 = w)[(int)(ptr1 = (IntPtr)((z + num3) - 1))] = numArray1[(int)ptr1] * num12;
                        num3 += n;
                    }
                }
                 
            }
        }
         
    }

    private void l_l14ng(int n, double[] x, int nact, double[] w, int g, int z, int ztg, int xs, int gs, double[] zznorm) throws Exception {
        int num2;
        int num3;
        int num7;
        double num9 = 0;
        double num10 = 0;
        double num12 = 0;
        int num1 = 1;
        while (num1 <= n)
        {
            num2 = num1 - 1;
            w[xs + num2] = x[num2] - w[xs + num2];
            num9 += Math.Pow(w[xs + num2], 2);
            num12 += w[gs + num2] * w[xs + num2];
            w[gs + num2] = w[g + num2] - w[gs + num2];
            num10 += w[gs + num2] * w[xs + num2];
            num1++;
        }
        if (num10 < (0.1 * Math.Abs(num12)))
        {
            return ;
        }
         
        int num4 = n;
        Label_00BC:num7 = num4;
        num4--;
        if (num4 > nact)
        {
            if (w[(ztg + num7) - 1] != 0)
            {
                num12 = Math.Abs(w[(ztg + num7) - 1]) * Math.Sqrt(1 + Math.Pow(w[(ztg + num4) - 1] / w[(ztg + num7) - 1], 2));
                double num13 = w[(ztg + num4) - 1] / num12;
                double num14 = w[(ztg + num7) - 1] / num12;
                w[(ztg + num4) - 1] = num12;
                num3 = num4;
                for (num1 = 1;num1 <= n;num1++)
                {
                    num2 = num1 - 1;
                    num12 = (num13 * w[z + num3]) - (num14 * w[(z + num3) - 1]);
                    w[(z + num3) - 1] = (num13 * w[(z + num3) - 1]) + (num14 * w[z + num3]);
                    w[z + num3] = num12;
                    num3 += n;
                }
            }
             
            goto Label_00BC
        }
         
        if (zznorm[0] < 0)
        {
            zznorm[0] = num9 / num10;
        }
        else
        {
            num12 = Math.Sqrt((zznorm[0] * num9) / num10);
            zznorm[0] = Math.Min(zznorm[0], num12);
            zznorm[0] = Math.Max(zznorm[0], 0.1 * num12);
        } 
        int num8 = nact + 1;
        num12 = Math.Sqrt(num10);
        num3 = num8;
        num1 = 1;
        while (num1 <= n)
        {
            num2 = num1 - 1;
            w[(z + num3) - 1] = w[xs + num2] / num12;
            num3 += n;
            num1++;
        }
        if (num8 < n)
        {
            int num6 = num8 + 1;
            for (num4 = num6;num4 <= n;num4++)
            {
                double[] numArray1 = new double[]();
                IntPtr ptr1 = new IntPtr();
                int num5 = num4 - 1;
                num12 = 0;
                num3 = num4;
                num1 = 1;
                while (num1 <= n)
                {
                    num2 = num1 - 1;
                    num12 += w[gs + num2] * w[(z + num3) - 1];
                    num3 += n;
                    num1++;
                }
                num12 /= num10;
                double num11 = 0;
                num3 = num4;
                num1 = 1;
                while (num1 <= n)
                {
                    num2 = num1 - 1;
                    (numArray1 = w)[(int)(ptr1 = (IntPtr)((z + num3) - 1))] = numArray1[(int)ptr1] + (-num12 * w[xs + num2]);
                    num11 += Math.Pow(w[(z + num3) - 1], 2);
                    num3 += n;
                    num1++;
                }
                if (num11 < zznorm[0])
                {
                    num12 = Math.Sqrt(zznorm[0] / num11);
                    num3 = num4;
                    for (num1 = 1;num1 <= n;num1++)
                    {
                        num2 = num1 - 1;
                        (numArray1 = w)[(int)(ptr1 = (IntPtr)((z + num3) - 1))] = numArray1[(int)ptr1] * num12;
                        num3 += n;
                    }
                }
                 
            }
        }
         
    }

    private void l_l15nf(int n, int m, double[] a, int ia, int[] iact, int[] nact, double[] w, int z, int u, double relacc, int idrop) throws Exception {
        int num8 = 0;
        int num5 = 0;
        int num12 = nact[0] - 1;
        if (idrop == nact[0])
        {
            nact[0] = num12;
        }
        else
        {
            int num6 = iact[idrop - 1];
            for (int num9 = idrop;num9 <= num12;num9++)
            {
                int num1;
                int num2;
                int num7;
                double num14 = new double();
                int num10 = num9 - 1;
                int num11 = num9 + 1;
                int num4 = iact[num11 - 1];
                iact[num10] = num4;
                if (num4 <= m)
                {
                    num14 = 0;
                    num7 = num9;
                    for (num1 = 1;num1 <= n;num1++)
                    {
                        num2 = num1 - 1;
                        num14 += w[(z + num7) - 1] * a[(num2 * ia) + (num4 - 1)];
                        num7 += n;
                    }
                }
                else
                {
                    int num3 = num4 - m;
                    if (num3 <= n)
                    {
                        num8 = (num3 * n) - n;
                        num14 = -w[(z + num8) + num10];
                    }
                    else
                    {
                        num3 -= n;
                        num8 = (num3 * n) - n;
                        num14 = w[(z + num8) + num10];
                    } 
                } 
                double num19 = w[(u + num11) - 1];
                double num16 = num14 * num19;
                double num13 = Math.Abs(num16);
                if ((num13 * relacc) < 1)
                {
                    num13 = Math.Sqrt(1 + (num13 * num13));
                }
                 
                double num20 = num16 / num13;
                double num22 = 1 / num13;
                num7 = num9;
                if (num4 > m)
                {
                    for (num1 = 1;num1 <= n;num1++)
                    {
                        num2 = num1 - 1;
                        num16 = (num20 * w[z + num7]) - (num22 * w[(z + num7) - 1]);
                        w[(z + num7) - 1] = (num20 * w[(z + num7) - 1]) + (num22 * w[z + num7]);
                        w[z + num7] = num16;
                        num7 += n;
                    }
                    w[((z + num8) + num11) - 1] = 0;
                }
                else
                {
                    double num21 = 0;
                    num1 = 1;
                    while (num1 <= n)
                    {
                        num2 = num1 - 1;
                        double num17 = num20 * w[z + num7];
                        double num18 = num22 * w[(z + num7) - 1];
                        num16 = Math.Abs(a[(num2 * ia) + (num4 - 1)]) * (Math.Abs(num17) + Math.Abs(num18));
                        if (num16 > num21)
                        {
                            num21 = num16;
                            num5 = num1;
                        }
                         
                        w[(z + num7) - 1] = (num20 * w[(z + num7) - 1]) + (num22 * w[z + num7]);
                        w[z + num7] = num17 - num18;
                        num7 += n;
                        num1++;
                    }
                    double num15 = 0;
                    num7 = num11;
                    for (num1 = 1;num1 <= n;num1++)
                    {
                        num2 = num1 - 1;
                        num15 += w[(z + num7) - 1] * a[(num2 * ia) + (num4 - 1)];
                        num7 += n;
                    }
                    if (num15 != 0)
                    {
                        double[] numArray1 = new double[]();
                        IntPtr ptr1 = new IntPtr();
                        num7 = ((num5 * n) - n) + num11;
                        (numArray1 = w)[(int)(ptr1 = (IntPtr)((z + num7) - 1))] = numArray1[(int)ptr1] + (-num15 / a[((num5 - 1) * ia) + (num4 - 1)]);
                    }
                     
                } 
                w[(u + num11) - 1] = -num13 * w[u + num10];
                w[u + num10] = num19 / num13;
            }
            iact[nact[0] - 1] = num6;
            nact[0] = num12;
        } 
    }

    private void l_l15ng(int n, int m, double[] a, int ia, int[] iact, int[] nact, double[] w, int z, int u, double relacc, int idrop) throws Exception {
        int num8 = 0;
        int num5 = 0;
        int num12 = nact[0] - 1;
        if (idrop == nact[0])
        {
            nact[0] = num12;
        }
        else
        {
            int num6 = iact[idrop - 1];
            for (int num9 = idrop;num9 <= num12;num9++)
            {
                int num1;
                int num2;
                int num7;
                double num14 = new double();
                int num10 = num9 - 1;
                int num11 = num9 + 1;
                int num4 = iact[num11 - 1];
                iact[num10] = num4;
                if (num4 <= m)
                {
                    num14 = 0;
                    num7 = num9;
                    for (num1 = 1;num1 <= n;num1++)
                    {
                        num2 = num1 - 1;
                        num14 += w[(z + num7) - 1] * a[(num2 * ia) + (num4 - 1)];
                        num7 += n;
                    }
                }
                else
                {
                    int num3 = num4 - m;
                    if (num3 <= n)
                    {
                        num8 = (num3 * n) - n;
                        num14 = -w[(z + num8) + num10];
                    }
                    else
                    {
                        num3 -= n;
                        num8 = (num3 * n) - n;
                        num14 = w[(z + num8) + num10];
                    } 
                } 
                double num19 = w[(u + num11) - 1];
                double num16 = num14 * num19;
                double num13 = Math.Abs(num16);
                if ((num13 * relacc) < 1)
                {
                    num13 = Math.Sqrt(1 + (num13 * num13));
                }
                 
                double num20 = num16 / num13;
                double num22 = 1 / num13;
                num7 = num9;
                if (num4 > m)
                {
                    for (num1 = 1;num1 <= n;num1++)
                    {
                        num2 = num1 - 1;
                        num16 = (num20 * w[z + num7]) - (num22 * w[(z + num7) - 1]);
                        w[(z + num7) - 1] = (num20 * w[(z + num7) - 1]) + (num22 * w[z + num7]);
                        w[z + num7] = num16;
                        num7 += n;
                    }
                    w[((z + num8) + num11) - 1] = 0;
                }
                else
                {
                    double num21 = 0;
                    num1 = 1;
                    while (num1 <= n)
                    {
                        num2 = num1 - 1;
                        double num17 = num20 * w[z + num7];
                        double num18 = num22 * w[(z + num7) - 1];
                        num16 = Math.Abs(a[(num2 * ia) + (num4 - 1)]) * (Math.Abs(num17) + Math.Abs(num18));
                        if (num16 > num21)
                        {
                            num21 = num16;
                            num5 = num1;
                        }
                         
                        w[(z + num7) - 1] = (num20 * w[(z + num7) - 1]) + (num22 * w[z + num7]);
                        w[z + num7] = num17 - num18;
                        num7 += n;
                        num1++;
                    }
                    double num15 = 0;
                    num7 = num11;
                    for (num1 = 1;num1 <= n;num1++)
                    {
                        num2 = num1 - 1;
                        num15 += w[(z + num7) - 1] * a[(num2 * ia) + (num4 - 1)];
                        num7 += n;
                    }
                    if (num15 != 0)
                    {
                        double[] numArray1 = new double[]();
                        IntPtr ptr1 = new IntPtr();
                        num7 = ((num5 * n) - n) + num11;
                        (numArray1 = w)[(int)(ptr1 = (IntPtr)((z + num7) - 1))] = numArray1[(int)ptr1] + (-num15 / a[((num5 - 1) * ia) + (num4 - 1)]);
                    }
                     
                } 
                w[(u + num11) - 1] = -num13 * w[u + num10];
                w[u + num10] = num19 / num13;
            }
            iact[nact[0] - 1] = num6;
            nact[0] = num12;
        } 
    }

    private void l_l16nf(int n, int m, double[] a, int ia, int[] iact, int[] nact, double[] par, double[] w, int g, int z, int u, int d, int ztg, double relacc, double[] ddotg, int meql, int mdeg, int gm, int gmnew, int parnew, int cgrad) throws Exception {
        int num1;
        int num2;
        double[] numArray2 = new double[1];
        double[] numArray1 = numArray2;
        int num6 = 0;
        int num7 = 1;
        Label_0017:if (num7 == 1)
        {
            for (num1 = 1;num1 <= n;num1++)
            {
                num2 = num1 - 1;
                w[gm + num2] = w[g + num2];
            }
            num6 = nact[0];
        }
         
        if (num6 > 0)
        {
            IntPtr ptr1 = new IntPtr();
            double num8 = 0;
            int num3 = num6;
            num1 = 1;
            while (num1 <= n)
            {
                num2 = num1 - 1;
                num8 += w[(z + num3) - 1] * w[gm + num2];
                num3 += n;
                num1++;
            }
            num8 *= w[(u + num6) - 1];
            if ((num6 > meql) && (num8 > 0))
            {
                this.l_l15nf(n, m, a, ia, iact, nact, w, z, u, relacc, num6);
                num7 = 1;
                goto Label_0017
            }
             
            int num4 = iact[num6 - 1];
            if (num4 <= m)
            {
                for (num1 = 1;num1 <= n;num1++)
                {
                    num2 = num1 - 1;
                    (numArray2 = w)[(int)(ptr1 = (IntPtr)(gm + num2))] = numArray2[(int)ptr1] + (-num8 * a[(num2 * ia) + (num4 - 1)]);
                }
            }
            else
            {
                int num5 = num4 - m;
                if (num5 <= n)
                {
                    (numArray2 = w)[(int)(ptr1 = (IntPtr)((gm + num5) - 1))] = numArray2[(int)ptr1] + num8;
                }
                else
                {
                    (numArray2 = w)[(int)(ptr1 = (IntPtr)(((gm + num5) - n) - 1))] = numArray2[(int)ptr1] - num8;
                } 
            } 
            par[num6 - 1] = num8;
            num6--;
            num7 = 0;
            goto Label_0017
        }
         
        ddotg[0] = 0;
        if (nact[0] < n)
        {
            this.l_l18nf(n, m, a, ia, iact, nact, par, w, z, u, d, ztg, gm, relacc, numArray1, meql, mdeg, gmnew, parnew, cgrad);
            if (numArray1[0] < 0)
            {
                for (num1 = 1;num1 <= n;num1++)
                {
                    num2 = num1 - 1;
                    (numArray2 = ddotg)[0] = numArray2[0] + (w[d + num2] * w[g + num2]);
                }
            }
             
        }
         
    }

    private void l_l16ng(int n, int m, double[] a, int ia, int[] iact, int[] nact, double[] par, double[] w, int g, int z, int u, int d, int ztg, double relacc, double[] ddotg, int meql, int mdeg, int gm, int gmnew, int parnew, int cgrad) throws Exception {
        int num1;
        int num2;
        double[] numArray2 = new double[1];
        double[] numArray1 = numArray2;
        int num6 = 0;
        int num7 = 1;
        Label_0017:if (num7 == 1)
        {
            for (num1 = 1;num1 <= n;num1++)
            {
                num2 = num1 - 1;
                w[gm + num2] = w[g + num2];
            }
            num6 = nact[0];
        }
         
        if (num6 > 0)
        {
            IntPtr ptr1 = new IntPtr();
            double num8 = 0;
            int num3 = num6;
            num1 = 1;
            while (num1 <= n)
            {
                num2 = num1 - 1;
                num8 += w[(z + num3) - 1] * w[gm + num2];
                num3 += n;
                num1++;
            }
            num8 *= w[(u + num6) - 1];
            if ((num6 > meql) && (num8 > 0))
            {
                this.l_l15ng(n, m, a, ia, iact, nact, w, z, u, relacc, num6);
                num7 = 1;
                goto Label_0017
            }
             
            int num4 = iact[num6 - 1];
            if (num4 <= m)
            {
                for (num1 = 1;num1 <= n;num1++)
                {
                    num2 = num1 - 1;
                    (numArray2 = w)[(int)(ptr1 = (IntPtr)(gm + num2))] = numArray2[(int)ptr1] + (-num8 * a[(num2 * ia) + (num4 - 1)]);
                }
            }
            else
            {
                int num5 = num4 - m;
                if (num5 <= n)
                {
                    (numArray2 = w)[(int)(ptr1 = (IntPtr)((gm + num5) - 1))] = numArray2[(int)ptr1] + num8;
                }
                else
                {
                    (numArray2 = w)[(int)(ptr1 = (IntPtr)(((gm + num5) - n) - 1))] = numArray2[(int)ptr1] - num8;
                } 
            } 
            par[num6 - 1] = num8;
            num6--;
            num7 = 0;
            goto Label_0017
        }
         
        ddotg[0] = 0;
        if (nact[0] < n)
        {
            this.l_l18ng(n, m, a, ia, iact, nact, par, w, z, u, d, ztg, gm, relacc, numArray1, meql, mdeg, gmnew, parnew, cgrad);
            if (numArray1[0] < 0)
            {
                for (num1 = 1;num1 <= n;num1++)
                {
                    num2 = num1 - 1;
                    (numArray2 = ddotg)[0] = numArray2[0] + (w[d + num2] * w[g + num2]);
                }
            }
             
        }
         
    }

    private void l_l17nf(int n, int m, double[] a, int ia, int[] iact, double[] w, int bres, int d, double[] stepcb, double[] ddotg, int mdeg, int[] msat, int mtot, int[] indxbd) throws Exception {
        int num10 = new int();
        int num11 = new int();
        double[] numArray1 = new double[]();
        int[] numArray2 = new int[]();
        int num4 = 0;
        double num13 = 0;
        int num3 = 0;
        stepcb[0] = 0;
        indxbd[0] = 0;
        int num6 = mdeg;
        int num9 = num10 = num11 = 1;
        Label_0033:if (num9 == 1)
        {
            num6++;
            if (num6 > mtot)
            {
                num9 = 0;
                num10 = 0;
                num11 = 1;
            }
             
        }
         
        if (num10 == 1)
        {
            num4 = iact[num6 - 1];
            if (num4 <= m)
            {
                num13 = 0;
                for (int num1 = 1;num1 <= n;num1++)
                {
                    int num2 = num1 - 1;
                    num13 += w[d + num2] * a[(num2 * ia) + (num4 - 1)];
                }
            }
            else
            {
                int num5 = num4 - m;
                if (num5 <= n)
                {
                    num13 = -w[(d + num5) - 1];
                }
                else
                {
                    num13 = w[((d + num5) - n) - 1];
                } 
            } 
            int num12 = 0;
            if (num3 == 1)
            {
                num11 = 0;
                num12 = 1;
            }
             
            if (num12 == 0)
            {
                if ((num13 * w[(bres + num4) - 1]) <= 0)
                {
                    w[(bres + num4) - 1] = 0;
                }
                else
                {
                    IntPtr ptr1 = new IntPtr();
                    (numArray1 = w)[(int)(ptr1 = (IntPtr)((bres + num4) - 1))] = numArray1[(int)ptr1] / num13;
                    if ((stepcb[0] == 0) || (w[(bres + num4) - 1] < stepcb[0]))
                    {
                        stepcb[0] = w[(bres + num4) - 1];
                        indxbd[0] = num6;
                    }
                     
                } 
                num9 = 1;
                num10 = 1;
                num11 = 1;
                goto Label_0033
            }
             
        }
         
        if (num11 == 1)
        {
            if (indxbd[0] <= msat[0])
            {
                return ;
            }
             
            num3 = 1;
            num6 = indxbd[0];
            num9 = 0;
            num10 = 1;
            num11 = 1;
            goto Label_0033
        }
         
        (numArray2 = msat)[0] = numArray2[0] + 1;
        iact[indxbd[0] - 1] = iact[msat[0] - 1];
        iact[msat[0] - 1] = num4;
        w[(bres + num4) - 1] = 0;
        indxbd[0] = msat[0];
        (numArray1 = ddotg)[0] = numArray1[0] - num13;
        if ((ddotg[0] < 0) && (msat[0] < mtot))
        {
            double num14 = 0;
            int num8 = mdeg + 1;
            for (num6 = num8;num6 <= mtot;num6++)
            {
                int num7 = num6 - 1;
                num4 = iact[num7];
                if ((w[(bres + num4) - 1] > 0) && ((num14 == 0) || (w[(bres + num4) - 1] < num14)))
                {
                    num14 = w[(bres + num4) - 1];
                    indxbd[0] = num6;
                }
                 
            }
            if (num14 > 0)
            {
                stepcb[0] = num14;
                num9 = 0;
                num10 = 0;
                num11 = 1;
                goto Label_0033
            }
             
        }
         
    }

    private void l_l17ng(int n, int m, double[] a, int ia, int[] iact, double[] w, int bres, int d, double[] stepcb, double[] ddotg, int mdeg, int[] msat, int mtot, int[] indxbd) throws Exception {
        int num10 = new int();
        int num11 = new int();
        double[] numArray1 = new double[]();
        int[] numArray2 = new int[]();
        int num4 = 0;
        double num13 = 0;
        int num3 = 0;
        stepcb[0] = 0;
        indxbd[0] = 0;
        int num6 = mdeg;
        int num9 = num10 = num11 = 1;
        Label_0033:if (num9 == 1)
        {
            num6++;
            if (num6 > mtot)
            {
                num9 = 0;
                num10 = 0;
                num11 = 1;
            }
             
        }
         
        if (num10 == 1)
        {
            num4 = iact[num6 - 1];
            if (num4 <= m)
            {
                num13 = 0;
                for (int num1 = 1;num1 <= n;num1++)
                {
                    int num2 = num1 - 1;
                    num13 += w[d + num2] * a[(num2 * ia) + (num4 - 1)];
                }
            }
            else
            {
                int num5 = num4 - m;
                if (num5 <= n)
                {
                    num13 = -w[(d + num5) - 1];
                }
                else
                {
                    num13 = w[((d + num5) - n) - 1];
                } 
            } 
            int num12 = 0;
            if (num3 == 1)
            {
                num11 = 0;
                num12 = 1;
            }
             
            if (num12 == 0)
            {
                if ((num13 * w[(bres + num4) - 1]) <= 0)
                {
                    w[(bres + num4) - 1] = 0;
                }
                else
                {
                    IntPtr ptr1 = new IntPtr();
                    (numArray1 = w)[(int)(ptr1 = (IntPtr)((bres + num4) - 1))] = numArray1[(int)ptr1] / num13;
                    if ((stepcb[0] == 0) || (w[(bres + num4) - 1] < stepcb[0]))
                    {
                        stepcb[0] = w[(bres + num4) - 1];
                        indxbd[0] = num6;
                    }
                     
                } 
                num9 = 1;
                num10 = 1;
                num11 = 1;
                goto Label_0033
            }
             
        }
         
        if (num11 == 1)
        {
            if (indxbd[0] <= msat[0])
            {
                return ;
            }
             
            num3 = 1;
            num6 = indxbd[0];
            num9 = 0;
            num10 = 1;
            num11 = 1;
            goto Label_0033
        }
         
        (numArray2 = msat)[0] = numArray2[0] + 1;
        iact[indxbd[0] - 1] = iact[msat[0] - 1];
        iact[msat[0] - 1] = num4;
        w[(bres + num4) - 1] = 0;
        indxbd[0] = msat[0];
        (numArray1 = ddotg)[0] = numArray1[0] - num13;
        if ((ddotg[0] < 0) && (msat[0] < mtot))
        {
            double num14 = 0;
            int num8 = mdeg + 1;
            for (num6 = num8;num6 <= mtot;num6++)
            {
                int num7 = num6 - 1;
                num4 = iact[num7];
                if ((w[(bres + num4) - 1] > 0) && ((num14 == 0) || (w[(bres + num4) - 1] < num14)))
                {
                    num14 = w[(bres + num4) - 1];
                    indxbd[0] = num6;
                }
                 
            }
            if (num14 > 0)
            {
                stepcb[0] = num14;
                num9 = 0;
                num10 = 0;
                num11 = 1;
                goto Label_0033
            }
             
        }
         
    }

    private void l_l18nf(int n, int m, double[] a, int ia, int[] iact, int[] nact, double[] par, double[] w, int z, int u, int d, int ztg, int gm, double relacc, double[] ddotgm, int meql, int mdeg, int gmnew, int parnew, int cgrad) throws Exception {
        int num1;
        int num2;
        int num6;
        int num10 = new int();
        double[] numArray1 = new double[]();
        IntPtr ptr1 = new IntPtr();
        int num4 = 0;
        int num9 = 0;
        int num3 = 0;
        int num12 = meql + 1;
        double num16 = 0;
        int num14 = 1;
        int num15 = 1;
        Label_0023:if (num14 == 1)
        {
            this.l_l19nf(n, nact[0], w, z, d, ztg, gm, relacc, ddotgm);
            if (ddotgm[0] == 0)
            {
                return ;
            }
             
            if (nact[0] == mdeg)
            {
                return ;
            }
             
            int num13 = nact[0] + 1;
            double num18 = 0;
            for (num6 = num13;num6 <= n;num6++)
            {
                int num7 = num6 - 1;
                num18 += Math.Pow(w[ztg + num7], 2);
            }
            if ((num16 > 0) && (num18 >= num16))
            {
                if (num4 == 1)
                {
                    return ;
                }
                 
                num4 = 1;
            }
            else
            {
                num16 = num18;
                num4 = 0;
            } 
            num9 = nact[0];
            this.l_l20nf(n, m, a, ia, iact, nact, w, z, u, d, relacc, mdeg, gmnew, parnew, cgrad);
            if (nact[0] == num9)
            {
                return ;
            }
             
            par[nact[0] - 1] = 0;
        }
         
        if (num15 == 1)
        {
            for (num1 = 1;num1 <= n;num1++)
            {
                num2 = num1 - 1;
                w[gmnew + num2] = w[gm + num2];
            }
            num9 = nact[0];
        }
         
        double num19 = 0;
        int num5 = num9;
        num1 = 1;
        while (num1 <= n)
        {
            num2 = num1 - 1;
            num19 += w[(z + num5) - 1] * w[gmnew + num2];
            num5 += n;
            num1++;
        }
        num19 *= w[(u + num9) - 1];
        w[(parnew + num9) - 1] = par[num9 - 1] + num19;
        if (num9 == nact[0])
        {
            w[(parnew + num9) - 1] = Math.Min(w[(parnew + num9) - 1], 0);
        }
         
        num6 = iact[num9 - 1];
        if (num6 <= m)
        {
            for (num1 = 1;num1 <= n;num1++)
            {
                num2 = num1 - 1;
                (numArray1 = w)[(int)(ptr1 = (IntPtr)(gmnew + num2))] = numArray1[(int)ptr1] + (-num19 * a[(num2 * ia) + (num6 - 1)]);
            }
        }
        else
        {
            int num8 = num6 - m;
            if (num8 <= n)
            {
                (numArray1 = w)[(int)(ptr1 = (IntPtr)((gmnew + num8) - 1))] = numArray1[(int)ptr1] + num19;
            }
            else
            {
                (numArray1 = w)[(int)(ptr1 = (IntPtr)(((gmnew + num8) - n) - 1))] = numArray1[(int)ptr1] - num19;
            } 
        } 
        num9--;
        if (num9 > meql)
        {
            num14 = 0;
            num15 = 0;
            goto Label_0023
        }
         
        double num17 = 0;
        if (num12 < nact[0])
        {
            int num11 = nact[0] - 1;
            for (num9 = num12;num9 <= num11;num9++)
            {
                num10 = num9 - 1;
                if (w[parnew + num10] > 0)
                {
                    num17 = w[parnew + num10] / (w[parnew + num10] - par[num10]);
                    num3 = num9;
                }
                 
            }
        }
         
        double num20 = 1 - num17;
        num9 = num12;
        while (num9 <= nact[0])
        {
            num10 = num9 - 1;
            par[num10] = Math.Min((double)((num20 * w[parnew + num10]) + (num17 * par[num10])), (double)0);
            num9++;
        }
        for (num1 = 1;num1 <= n;num1++)
        {
            num2 = num1 - 1;
            w[gm + num2] = (num20 * w[gmnew + num2]) + (num17 * w[gm + num2]);
        }
        if (num17 > 0)
        {
            this.l_l15nf(n, m, a, ia, iact, nact, w, z, u, relacc, num3);
            for (num9 = num3;num9 <= nact[0];num9++)
            {
                num10 = num9 - 1;
                par[num10] = par[num10 + 1];
            }
            num14 = 0;
            num15 = 1;
            goto Label_0023
        }
         
        if (nact[0] < n)
        {
            num14 = 1;
            num15 = 1;
            goto Label_0023
        }
         
        ddotgm[0] = 0;
    }

    private void l_l18ng(int n, int m, double[] a, int ia, int[] iact, int[] nact, double[] par, double[] w, int z, int u, int d, int ztg, int gm, double relacc, double[] ddotgm, int meql, int mdeg, int gmnew, int parnew, int cgrad) throws Exception {
        int num4 = 0;
        int num3 = 0;
        int num12 = meql + 1;
        double num14 = 0;
        Label_001A:this.l_l19ng(n, nact[0], w, z, d, ztg, gm, relacc, ddotgm);
        if ((ddotgm[0] != 0) && (nact[0] != mdeg))
        {
            int num13 = nact[0] + 1;
            double num16 = 0;
            int num6 = num13;
            while (num6 <= n)
            {
                int num7 = num6 - 1;
                num16 += Math.Pow(w[ztg + num7], 2);
                num6++;
            }
            if ((num14 > 0) && (num16 >= num14))
            {
                if (num4 == 1)
                {
                    return ;
                }
                 
                num4 = 1;
            }
            else
            {
                num14 = num16;
                num4 = 0;
            } 
            int num9 = nact[0];
            this.l_l20ng(n, m, a, ia, iact, nact, w, z, u, d, relacc, mdeg, gmnew, parnew, cgrad);
            if (nact[0] != num9)
            {
                int num2;
                int num10 = new int();
                double[] numArray1 = new double[]();
                IntPtr ptr1 = new IntPtr();
                par[nact[0] - 1] = 0;
                int num1 = 1;
                while (num1 <= n)
                {
                    num2 = num1 - 1;
                    w[gmnew + num2] = w[gm + num2];
                    num1++;
                }
                num9 = nact[0];
                double num17 = 0;
                int num5 = num9;
                num1 = 1;
                while (num1 <= n)
                {
                    num2 = num1 - 1;
                    num17 += w[(z + num5) - 1] * w[gmnew + num2];
                    num5 += n;
                    num1++;
                }
                num17 *= w[(u + num9) - 1];
                w[(parnew + num9) - 1] = par[num9 - 1] + num17;
                if (num9 == nact[0])
                {
                    w[(parnew + num9) - 1] = Math.Min(w[(parnew + num9) - 1], 0);
                }
                 
                num6 = iact[num9 - 1];
                if (num6 <= m)
                {
                    for (num1 = 1;num1 <= n;num1++)
                    {
                        num2 = num1 - 1;
                        (numArray1 = w)[(int)(ptr1 = (IntPtr)(gmnew + num2))] = numArray1[(int)ptr1] + (-num17 * a[(num2 * ia) + (num6 - 1)]);
                    }
                }
                else
                {
                    int num8 = num6 - m;
                    if (num8 <= n)
                    {
                        (numArray1 = w)[(int)(ptr1 = (IntPtr)((gmnew + num8) - 1))] = numArray1[(int)ptr1] + num17;
                    }
                    else
                    {
                        (numArray1 = w)[(int)(ptr1 = (IntPtr)(((gmnew + num8) - n) - 1))] = numArray1[(int)ptr1] - num17;
                    } 
                } 
                num9--;
                if (num9 > meql)
                {
                    goto Label_001A
                }
                 
                double num15 = 0;
                if (num12 < nact[0])
                {
                    int num11 = nact[0] - 1;
                    for (num9 = num12;num9 <= num11;num9++)
                    {
                        num10 = num9 - 1;
                        if (w[parnew + num10] > 0)
                        {
                            num15 = w[parnew + num10] / (w[parnew + num10] - par[num10]);
                            num3 = num9;
                        }
                         
                    }
                }
                 
                double num18 = 1 - num15;
                num9 = num12;
                while (num9 <= nact[0])
                {
                    num10 = num9 - 1;
                    par[num10] = Math.Min((double)((num18 * w[parnew + num10]) + (num15 * par[num10])), (double)0);
                    num9++;
                }
                for (num1 = 1;num1 <= n;num1++)
                {
                    num2 = num1 - 1;
                    w[gm + num2] = (num18 * w[gmnew + num2]) + (num15 * w[gm + num2]);
                }
                if (num15 > 0)
                {
                    this.l_l15ng(n, m, a, ia, iact, nact, w, z, u, relacc, num3);
                    for (num9 = num3;num9 <= nact[0];num9++)
                    {
                        num10 = num9 - 1;
                        par[num10] = par[num10 + 1];
                    }
                    goto Label_001A
                }
                 
                if (nact[0] < n)
                {
                    goto Label_001A
                }
                 
                ddotgm[0] = 0;
            }
             
        }
         
    }

    private void l_l19nf(int n, int nact, double[] w, int z, int d, int ztg, int gm, double relacc, double[] ddotgm) throws Exception {
        ddotgm[0] = 0;
        if (nact < n)
        {
            int num1;
            int num2;
            int num3;
            int num5;
            double num7;
            double num8;
            double num9;
            int num6 = nact + 1;
            int num4 = num6;
            while (num4 <= n)
            {
                num5 = num4 - 1;
                num7 = 0;
                num8 = 0;
                num3 = num4;
                for (num1 = 1;num1 <= n;num1++)
                {
                    num2 = num1 - 1;
                    num9 = w[(z + num3) - 1] * w[gm + num2];
                    num7 += num9;
                    num8 += Math.Abs(num9);
                    num3 += n;
                }
                if (Math.Abs(num7) <= (relacc * num8))
                {
                    num7 = 0;
                }
                 
                w[ztg + num5] = num7;
                num4++;
            }
            num3 = 0;
            num1 = 1;
            while (num1 <= n)
            {
                num2 = num1 - 1;
                num7 = 0;
                num8 = 0;
                for (num4 = num6;num4 <= n;num4++)
                {
                    num5 = num4 - 1;
                    num9 = w[(z + num3) + num5] * w[ztg + num5];
                    num7 -= num9;
                    num8 += Math.Abs(num9);
                }
                if (Math.Abs(num7) <= (relacc * num8))
                {
                    num7 = 0;
                }
                 
                w[d + num2] = num7;
                num3 += n;
                num1++;
            }
            num8 = 0;
            for (num1 = 1;num1 <= n;num1++)
            {
                double[] numArray1 = new double[]();
                num2 = num1 - 1;
                num9 = w[d + num2] * w[gm + num2];
                (numArray1 = ddotgm)[0] = numArray1[0] + num9;
                num8 += Math.Abs(num9);
            }
            if ((ddotgm[0] + (relacc * num8)) >= 0)
            {
                ddotgm[0] = 0;
            }
             
        }
         
    }

    private void l_l19ng(int n, int nact, double[] w, int z, int d, int ztg, int gm, double relacc, double[] ddotgm) throws Exception {
        ddotgm[0] = 0;
        if (nact < n)
        {
            int num1;
            int num2;
            int num3;
            int num5;
            double num7;
            double num8;
            double num9;
            int num6 = nact + 1;
            int num4 = num6;
            while (num4 <= n)
            {
                num5 = num4 - 1;
                num7 = 0;
                num8 = 0;
                num3 = num4;
                for (num1 = 1;num1 <= n;num1++)
                {
                    num2 = num1 - 1;
                    num9 = w[(z + num3) - 1] * w[gm + num2];
                    num7 += num9;
                    num8 += Math.Abs(num9);
                    num3 += n;
                }
                if (Math.Abs(num7) <= (relacc * num8))
                {
                    num7 = 0;
                }
                 
                w[ztg + num5] = num7;
                num4++;
            }
            num3 = 0;
            num1 = 1;
            while (num1 <= n)
            {
                num2 = num1 - 1;
                num7 = 0;
                num8 = 0;
                for (num4 = num6;num4 <= n;num4++)
                {
                    num5 = num4 - 1;
                    num9 = w[(z + num3) + num5] * w[ztg + num5];
                    num7 -= num9;
                    num8 += Math.Abs(num9);
                }
                if (Math.Abs(num7) <= (relacc * num8))
                {
                    num7 = 0;
                }
                 
                w[d + num2] = num7;
                num3 += n;
                num1++;
            }
            num8 = 0;
            for (num1 = 1;num1 <= n;num1++)
            {
                double[] numArray1 = new double[]();
                num2 = num1 - 1;
                num9 = w[d + num2] * w[gm + num2];
                (numArray1 = ddotgm)[0] = numArray1[0] + num9;
                num8 += Math.Abs(num9);
            }
            if ((ddotgm[0] + (relacc * num8)) >= 0)
            {
                ddotgm[0] = 0;
            }
             
        }
         
    }

    private void l_l20nf(int n, int m, double[] a, int ia, int[] iact, int[] nact, double[] w, int z, int u, int d, double relacc, int mdeg, int zzdiag, int gmnew, int cgrad) throws Exception {
        int num2;
        int num5;
        int num7;
        double num15 = new double();
        double num18 = new double();
        double num19 = new double();
        double num21 = new double();
        double[] numArray1 = new double[]();
        IntPtr ptr1 = new IntPtr();
        int num3 = 0;
        double num17 = 0;
        double num16 = 0;
        int num12 = nact[0] + 1;
        int num11 = mdeg;
        int num4 = 0;
        int num1 = 1;
        while (num1 <= n)
        {
            num2 = num1 - 1;
            w[zzdiag + num2] = 0;
            for (num5 = num12;num5 <= n;num5++)
            {
                int num6 = num5 - 1;
                (numArray1 = w)[(int)(ptr1 = (IntPtr)(zzdiag + num2))] = numArray1[(int)ptr1] + Math.Pow(w[(z + num4) + num6], 2);
            }
            num4 += n;
            num1++;
        }
        Label_0091:num15 = 0;
        int num9 = num12;
        while (num9 <= num11)
        {
            double num20 = new double();
            int num10 = num9 - 1;
            num5 = iact[num10];
            if (num5 <= m)
            {
                num18 = 0;
                num19 = 0;
                num20 = 0;
                for (num1 = 1;num1 <= n;num1++)
                {
                    num2 = num1 - 1;
                    num21 = w[d + num2] * a[(num2 * ia) + (num5 - 1)];
                    num18 += num21;
                    num19 += Math.Abs(num21);
                    num20 += w[zzdiag + num2] * Math.Pow(a[(num2 * ia) + (num5 - 1)], 2);
                }
            }
            else
            {
                num7 = num5 - m;
                if (num7 <= n)
                {
                    num18 = -w[(d + num7) - 1];
                }
                else
                {
                    num7 -= n;
                    num18 = w[(d + num7) - 1];
                } 
                num19 = Math.Abs(num18);
                num20 = w[(zzdiag + num7) - 1];
            } 
            if (num18 > (relacc * num19))
            {
                double num14 = (num18 * num18) / num20;
                if (num14 > num15)
                {
                    num15 = num14;
                    num3 = num9;
                    num17 = num18;
                    num16 = num19;
                }
                 
            }
             
            num9++;
        }
        if (num15 > 0)
        {
            int num13 = 0;
            if (nact[0] == 0)
            {
                num13 = 1;
            }
             
            if (num13 == 0)
            {
                int num8 = new int();
                num13 = 1;
                num5 = iact[num3 - 1];
                if (num5 <= m)
                {
                    num8 = 0;
                    for (num1 = 1;num1 <= n;num1++)
                    {
                        num2 = num1 - 1;
                        w[gmnew + num2] = a[(num2 * ia) + (num5 - 1)];
                    }
                }
                else
                {
                    num8 = num5 - m;
                    for (num1 = 1;num1 <= n;num1++)
                    {
                        num2 = num1 - 1;
                        w[gmnew + num2] = 0;
                    }
                    if (num8 <= n)
                    {
                        w[(gmnew + num8) - 1] = -1;
                    }
                    else
                    {
                        num8 -= n;
                        w[(gmnew + num8) - 1] = 1;
                    } 
                } 
                num9 = nact[0];
                do
                {
                    num21 = 0;
                    num4 = num9;
                    num1 = 1;
                    while (num1 <= n)
                    {
                        num2 = num1 - 1;
                        num21 += w[(z + num4) - 1] * w[gmnew + num2];
                        num4 += n;
                        num1++;
                    }
                    num21 *= w[(u + num9) - 1];
                    num5 = iact[num9 - 1];
                    if (num5 <= m)
                    {
                        for (num1 = 1;num1 <= n;num1++)
                        {
                            num2 = num1 - 1;
                            (numArray1 = w)[(int)(ptr1 = (IntPtr)(gmnew + num2))] = numArray1[(int)ptr1] + (-num21 * a[(num2 * ia) + (num5 - 1)]);
                        }
                    }
                    else
                    {
                        num7 = num5 - m;
                        if (num7 <= n)
                        {
                            (numArray1 = w)[(int)(ptr1 = (IntPtr)((gmnew + num7) - 1))] = numArray1[(int)ptr1] + num21;
                        }
                        else
                        {
                            (numArray1 = w)[(int)(ptr1 = (IntPtr)(((gmnew + num7) - n) - 1))] = numArray1[(int)ptr1] - num21;
                        } 
                    } 
                    num18 = 0;
                    num19 = 0;
                    for (num1 = 1;num1 <= n;num1++)
                    {
                        num2 = num1 - 1;
                        num21 = w[d + num2] * w[gmnew + num2];
                        num18 += num21;
                        num19 += Math.Abs(num21);
                    }
                    num17 = Math.Min(num17, num18);
                    num16 = Math.Max(num16, num19);
                    num9--;
                }
                while (num9 >= 1);
                if (num8 > 0)
                {
                    w[(d + num8) - 1] = 0;
                }
                 
                if (num17 <= (relacc * num16))
                {
                    num13 = 0;
                }
                 
            }
             
            if (num13 == 1)
            {
                num9 = nact[0];
                this.l_l9onf(n, m, a, ia, iact, nact, w, z, u, relacc, num3, gmnew, cgrad);
                if (nact[0] > num9)
                {
                    return ;
                }
                 
                num3 = num12;
            }
             
            if (num12 < num11)
            {
                num9 = iact[num11 - 1];
                iact[num11 - 1] = iact[num3 - 1];
                iact[num3 - 1] = num9;
                num11--;
                goto Label_0091
            }
             
        }
         
    }

    private void l_l20ng(int n, int m, double[] a, int ia, int[] iact, int[] nact, double[] w, int z, int u, int d, double relacc, int mdeg, int zzdiag, int gmnew, int cgrad) throws Exception {
        int num2;
        int num5;
        int num7;
        double num15 = new double();
        double num18 = new double();
        double num19 = new double();
        double num21 = new double();
        double[] numArray1 = new double[]();
        IntPtr ptr1 = new IntPtr();
        int num3 = 0;
        double num17 = 0;
        double num16 = 0;
        int num12 = nact[0] + 1;
        int num11 = mdeg;
        int num4 = 0;
        int num1 = 1;
        while (num1 <= n)
        {
            num2 = num1 - 1;
            w[zzdiag + num2] = 0;
            for (num5 = num12;num5 <= n;num5++)
            {
                int num6 = num5 - 1;
                (numArray1 = w)[(int)(ptr1 = (IntPtr)(zzdiag + num2))] = numArray1[(int)ptr1] + Math.Pow(w[(z + num4) + num6], 2);
            }
            num4 += n;
            num1++;
        }
        Label_0091:num15 = 0;
        int num9 = num12;
        while (num9 <= num11)
        {
            double num20 = new double();
            int num10 = num9 - 1;
            num5 = iact[num10];
            if (num5 <= m)
            {
                num18 = 0;
                num19 = 0;
                num20 = 0;
                for (num1 = 1;num1 <= n;num1++)
                {
                    num2 = num1 - 1;
                    num21 = w[d + num2] * a[(num2 * ia) + (num5 - 1)];
                    num18 += num21;
                    num19 += Math.Abs(num21);
                    num20 += w[zzdiag + num2] * Math.Pow(a[(num2 * ia) + (num5 - 1)], 2);
                }
            }
            else
            {
                num7 = num5 - m;
                if (num7 <= n)
                {
                    num18 = -w[(d + num7) - 1];
                }
                else
                {
                    num7 -= n;
                    num18 = w[(d + num7) - 1];
                } 
                num19 = Math.Abs(num18);
                num20 = w[(zzdiag + num7) - 1];
            } 
            if (num18 > (relacc * num19))
            {
                double num14 = (num18 * num18) / num20;
                if (num14 > num15)
                {
                    num15 = num14;
                    num3 = num9;
                    num17 = num18;
                    num16 = num19;
                }
                 
            }
             
            num9++;
        }
        if (num15 > 0)
        {
            int num13 = 0;
            if (nact[0] == 0)
            {
                num13 = 1;
            }
             
            if (num13 == 0)
            {
                int num8 = new int();
                num13 = 1;
                num5 = iact[num3 - 1];
                if (num5 <= m)
                {
                    num8 = 0;
                    for (num1 = 1;num1 <= n;num1++)
                    {
                        num2 = num1 - 1;
                        w[gmnew + num2] = a[(num2 * ia) + (num5 - 1)];
                    }
                }
                else
                {
                    num8 = num5 - m;
                    for (num1 = 1;num1 <= n;num1++)
                    {
                        num2 = num1 - 1;
                        w[gmnew + num2] = 0;
                    }
                    if (num8 <= n)
                    {
                        w[(gmnew + num8) - 1] = -1;
                    }
                    else
                    {
                        num8 -= n;
                        w[(gmnew + num8) - 1] = 1;
                    } 
                } 
                num9 = nact[0];
                do
                {
                    num21 = 0;
                    num4 = num9;
                    num1 = 1;
                    while (num1 <= n)
                    {
                        num2 = num1 - 1;
                        num21 += w[(z + num4) - 1] * w[gmnew + num2];
                        num4 += n;
                        num1++;
                    }
                    num21 *= w[(u + num9) - 1];
                    num5 = iact[num9 - 1];
                    if (num5 <= m)
                    {
                        for (num1 = 1;num1 <= n;num1++)
                        {
                            num2 = num1 - 1;
                            (numArray1 = w)[(int)(ptr1 = (IntPtr)(gmnew + num2))] = numArray1[(int)ptr1] + (-num21 * a[(num2 * ia) + (num5 - 1)]);
                        }
                    }
                    else
                    {
                        num7 = num5 - m;
                        if (num7 <= n)
                        {
                            (numArray1 = w)[(int)(ptr1 = (IntPtr)((gmnew + num7) - 1))] = numArray1[(int)ptr1] + num21;
                        }
                        else
                        {
                            (numArray1 = w)[(int)(ptr1 = (IntPtr)(((gmnew + num7) - n) - 1))] = numArray1[(int)ptr1] - num21;
                        } 
                    } 
                    num18 = 0;
                    num19 = 0;
                    for (num1 = 1;num1 <= n;num1++)
                    {
                        num2 = num1 - 1;
                        num21 = w[d + num2] * w[gmnew + num2];
                        num18 += num21;
                        num19 += Math.Abs(num21);
                    }
                    num17 = Math.Min(num17, num18);
                    num16 = Math.Max(num16, num19);
                    num9--;
                }
                while (num9 >= 1);
                if (num8 > 0)
                {
                    w[(d + num8) - 1] = 0;
                }
                 
                if (num17 <= (relacc * num16))
                {
                    num13 = 0;
                }
                 
            }
             
            if (num13 == 1)
            {
                num9 = nact[0];
                this.l_l9ong(n, m, a, ia, iact, nact, w, z, u, relacc, num3, gmnew, cgrad);
                if (nact[0] > num9)
                {
                    return ;
                }
                 
                num3 = num12;
            }
             
            if (num12 < num11)
            {
                num9 = iact[num11 - 1];
                iact[num11 - 1] = iact[num3 - 1];
                iact[num3 - 1] = num9;
                num11--;
                goto Label_0091
            }
             
        }
         
    }

    private void l_l21nf(Provisdom.Optimization.NlpLinearConstraintSolver.IFunction fcn, int n, double[] xc, double fc, double[] w, int gc, int[] nfcn) throws Exception {
        double[] numArray2 = new double[1];
        double[] numArray1 = numArray2;
        double num3 = Math.Sqrt(Math.Max((double)0, (double)2.2204460492503131E-16));
        for (int num1 = 1;num1 <= n;num1++)
        {
            int[] numArray3 = new int[]();
            int num2 = num1 - 1;
            double num4 = num3 * Math.Max(Math.Abs(xc[num2]), 1);
            if (xc[num2] < 0)
            {
                num4 = -num4;
            }
             
            double num5 = xc[num2];
            xc[num2] = num5 + num4;
            numArray1[0] = fcn.F(xc);
            (numArray3 = nfcn)[0] = numArray3[0] + 1;
            xc[num2] = num5;
            w[gc + num2] = (numArray1[0] - fc) / num4;
        }
    }

    private void l_l2onf(Provisdom.Optimization.NlpLinearConstraintSolver.IFunction fcn, int n, int m, int meq, double[] a, int ia, double[] b, double[] xl, double[] xu, double[] x, double acc, int[] maxfcn, double[] sol, double[] obj, int[] nact, int[] iact, double[] par, int iprint, int[] info, double[] w) throws Exception {
        BLAS.copy(n, x, sol);
        int num3 = 1;
        int num6 = num3 + n;
        int num10 = num6 + n;
        int num7 = num10 + (n * n);
        int num8 = num7 + n;
        int num1 = num8 + n;
        int num2 = ((num1 + m) + n) + n;
        int num11 = num2 + n;
        int num4 = num11 + n;
        int num9 = num4 + n;
        int num5 = num9 + n;
        info[0] = 0;
        this.l_l3onf(fcn, n, m, meq, a, ia, b, xl, xu, sol, obj, acc, maxfcn, iact, nact, par, iprint, info, w, num3 - 1, num10 - 1, num7 - 1, num8 - 1, num6 - 1, num1 - 1, num2 - 1, num11 - 1, num4 - 1, num9 - 1, num5 - 1);
        if (((info[0] != 1) && (info[0] != 2)) && (info[0] != 3))
        {
            if (info[0] == 5)
            {
                throw new ConstraintsInconsistentException();
            }
             
            if (info[0] == 6)
            {
                throw new VarBoundsInconsistentException();
            }
             
            if (info[0] == 7)
            {
                throw new ConstraintsNotSatisfiedException();
            }
             
            if (info[0] == 9)
            {
                throw new EqualityConstraintsException();
            }
             
        }
         
    }

    private void l_l2ong(Provisdom.Optimization.NlpLinearConstraintSolver.IFunction fcn, Provisdom.Optimization.NlpLinearConstraintSolver.IGradient grad, int n, int m, int meq, double[] a, int ia, double[] b, double[] xl, double[] xu, double[] x, double acc, int[] maxfcn, double[] sol, double[] obj, int[] nact, int[] iact, double[] par, int iprint, int[] info, double[] w) throws Exception {
        BLAS.copy(n, x, sol);
        int num3 = 1;
        int num6 = num3 + n;
        int num10 = num6 + n;
        int num7 = num10 + (n * n);
        int num8 = num7 + n;
        int num1 = num8 + n;
        int num2 = ((num1 + m) + n) + n;
        int num11 = num2 + n;
        int num4 = num11 + n;
        int num9 = num4 + n;
        int num5 = num9 + n;
        info[0] = 0;
        this.l_l3ong(fcn, grad, n, m, meq, a, ia, b, xl, xu, sol, obj, acc, maxfcn, iact, nact, par, iprint, info, w, num3 - 1, num10 - 1, num7 - 1, num8 - 1, num6 - 1, num1 - 1, num2 - 1, num11 - 1, num4 - 1, num9 - 1, num5 - 1);
        if (((info[0] != 1) && (info[0] != 2)) && (info[0] != 3))
        {
            if (info[0] == 5)
            {
                throw new ConstraintsInconsistentException();
            }
             
            if (info[0] == 6)
            {
                throw new VarBoundsInconsistentException();
            }
             
            if (info[0] == 7)
            {
                throw new ConstraintsNotSatisfiedException();
            }
             
            if (info[0] == 9)
            {
                throw new EqualityConstraintsException();
            }
             
        }
         
    }

    private void l_l3onf(Provisdom.Optimization.NlpLinearConstraintSolver.IFunction fcn, int n, int m, int meq, double[] a, int ia, double[] b, double[] xl, double[] xu, double[] x, double[] obj, double acc, int[] maxfcn, int[] iact, int[] nact, double[] par, int iprint, int[] info, double[] w, int g, int z, int u, int xbig, int reskt, int bres, int d, int ztg, int gm, int xs, int gs) throws Exception {
        int[] numArray8 = new int[1];
        int[] numArray1 = numArray8;
        numArray8 = new int[1];
        int[] numArray2 = numArray8;
        numArray8 = new int[1];
        int[] numArray3 = numArray8;
        numArray8 = new int[1];
        int[] numArray4 = numArray8;
        double[] numArray9 = new double[1];
        double[] numArray5 = numArray9;
        numArray9 = new double[1];
        double[] numArray6 = numArray9;
        numArray9 = new double[1];
        double[] numArray7 = numArray9;
        numArray7[0] = -1;
        numArray3[0] = 0;
        numArray4[0] = 0;
        int num7 = maxfcn[0];
        info[0] = 4;
        this.l_l4onf(n, m, xl, xu, x, iact, numArray1, info, w, z, u, xbig, numArray5);
        numArray6[0] = Math.Max((double)0.01, (double)(10 * numArray5[0]));
        if (meq > 0)
        {
            this.l_l5onf(n, m, meq, a, ia, b, xu, iact, numArray1, info, w, z, u, numArray5[0], xs, gs);
            if (info[0] == 5)
            {
                maxfcn[0] = numArray4[0];
                return ;
            }
             
        }
         
        nact[0] = numArray1[0];
        numArray2[0] = numArray1[0];
        int num6 = nact[0];
        for (int num1 = 1;num1 <= n;num1++)
        {
            int num2 = num1 - 1;
            if (xl[num2] < xu[num2])
            {
                num6 += 2;
                iact[num6 - 2] = m + num1;
                iact[num6 - 1] = (m + n) + num1;
            }
             
        }
        this.l_l6onf(n, m, a, ia, b, xl, xu, x, iact, nact, par, info, w, g, z, u, xbig, numArray5[0], numArray6, numArray1[0], numArray2, num6, bres, d, ztg, gm, reskt, xs, gs);
        if (numArray2[0] < num6)
        {
            info[0] = 6;
            maxfcn[0] = numArray4[0];
            return ;
        }
         
        if (m > meq)
        {
            int num5 = meq + 1;
            for (int num3 = num5;num3 <= m;num3++)
            {
                int num4 = num3 - 1;
                num6++;
                iact[num6 - 1] = num3;
            }
        }
         
        Label_01CA:this.l_l6onf(n, m, a, ia, b, xl, xu, x, iact, nact, par, info, w, g, z, u, xbig, numArray5[0], numArray6, numArray1[0], numArray2, num6, bres, d, ztg, gm, reskt, xs, gs);
        if (numArray2[0] < num6)
        {
            info[0] = 7;
            maxfcn[0] = numArray4[0];
        }
        else if (numArray1[0] == n)
        {
            info[0] = 9;
            maxfcn[0] = numArray4[0];
        }
        else
        {
            this.l_l7onf(fcn, n, m, a, ia, b, xl, xu, x, obj, acc, iact, nact, par, iprint, info, w, g, z, u, xbig, numArray5[0], numArray7, numArray6[0], numArray1[0], num6, numArray3, numArray4, num7, reskt, bres, d, ztg, gm, xs, gs);
            if ((numArray6[0] > numArray5[0]) && (nact[0] > 0))
            {
                if (numArray4[0] != num7)
                {
                    this.l_l8onf(n, m, a, ia, b, xl, xu, x, iact, nact[0], w, xbig, numArray5[0], numArray6, numArray1[0]);
                    goto Label_01CA
                }
                 
                info[0] = 8;
            }
             
            maxfcn[0] = numArray4[0];
        }  
    }

    private void l_l3ong(Provisdom.Optimization.NlpLinearConstraintSolver.IFunction fcn, Provisdom.Optimization.NlpLinearConstraintSolver.IGradient grad, int n, int m, int meq, double[] a, int ia, double[] b, double[] xl, double[] xu, double[] x, double[] obj, double acc, int[] maxfcn, int[] iact, int[] nact, double[] par, int iprint, int[] info, double[] w, int g, int z, int u, int xbig, int reskt, int bres, int d, int ztg, int gm, int xs, int gs) throws Exception {
        int[] numArray8 = new int[1];
        int[] numArray1 = numArray8;
        numArray8 = new int[1];
        int[] numArray2 = numArray8;
        numArray8 = new int[1];
        int[] numArray3 = numArray8;
        numArray8 = new int[1];
        int[] numArray4 = numArray8;
        double[] numArray9 = new double[1];
        double[] numArray5 = numArray9;
        numArray9 = new double[1];
        double[] numArray6 = numArray9;
        numArray9 = new double[1];
        double[] numArray7 = numArray9;
        numArray5[0] = -1;
        numArray3[0] = 0;
        numArray4[0] = 0;
        int num7 = maxfcn[0];
        info[0] = 4;
        this.l_l4ong(n, m, xl, xu, x, iact, numArray1, info, w, z, u, xbig, numArray6);
        numArray7[0] = Math.Max((double)0.01, (double)(10 * numArray6[0]));
        if (meq > 0)
        {
            this.l_l5ong(n, m, meq, a, ia, b, xu, iact, numArray1, info, w, z, u, numArray6, xs, gs);
            if (info[0] == 5)
            {
                maxfcn[0] = numArray4[0];
                return ;
            }
             
        }
         
        nact[0] = numArray1[0];
        numArray2[0] = numArray1[0];
        int num6 = nact[0];
        for (int num1 = 1;num1 <= n;num1++)
        {
            int num2 = num1 - 1;
            if (xl[num2] < xu[num2])
            {
                num6 += 2;
                iact[num6 - 2] = m + num1;
                iact[num6 - 1] = (m + n) + num1;
            }
             
        }
        this.l_l6ong(n, m, a, ia, b, xl, xu, x, iact, nact, par, info, w, g, z, u, xbig, numArray6[0], numArray7, numArray1[0], numArray2, num6, bres, d, ztg, gm, reskt, xs, gs);
        if (numArray2[0] < num6)
        {
            info[0] = 6;
            maxfcn[0] = numArray4[0];
            return ;
        }
         
        if (m > meq)
        {
            int num5 = meq + 1;
            for (int num3 = num5;num3 <= m;num3++)
            {
                int num4 = num3 - 1;
                num6++;
                iact[num6 - 1] = num3;
            }
        }
         
        Label_01CF:this.l_l6ong(n, m, a, ia, b, xl, xu, x, iact, nact, par, info, w, g, z, u, xbig, numArray6[0], numArray7, numArray1[0], numArray2, num6, bres, d, ztg, gm, reskt, xs, gs);
        if (numArray2[0] < num6)
        {
            info[0] = 7;
            maxfcn[0] = numArray4[0];
        }
        else if (numArray1[0] == n)
        {
            info[0] = 9;
            maxfcn[0] = numArray4[0];
        }
        else
        {
            this.l_l7ong(fcn, grad, n, m, a, ia, b, xl, xu, x, obj, acc, iact, nact, par, iprint, info, w, g, z, u, xbig, numArray6[0], numArray5, numArray7[0], numArray1[0], num6, numArray3, numArray4, num7, reskt, bres, d, ztg, gm, xs, gs);
            if ((numArray7[0] > numArray6[0]) && (nact[0] > 0))
            {
                if (numArray4[0] != num7)
                {
                    this.l_l8ong(n, m, a, ia, b, xl, xu, x, iact, nact[0], w, xbig, numArray6[0], numArray7, numArray1[0]);
                    goto Label_01CF
                }
                 
                info[0] = 8;
            }
             
            maxfcn[0] = numArray4[0];
        }  
    }

    private void l_l4onf(int n, int m, double[] xl, double[] xu, double[] x, int[] iact, int[] meql, int[] info, double[] w, int z, int u, int xbig, double[] relacc) throws Exception {
        int num2;
        double num7;
        double num8;
        double num9 = 100;
        relacc[0] = 1;
        do
        {
            double[] numArray1 = new double[]();
            (numArray1 = relacc)[0] = numArray1[0] * 0.5;
            num7 = num9 + (0.5 * relacc[0]);
            num8 = num9 + relacc[0];
            if (num9 >= num7)
            {
                break;
            }
             
        }
        while (num7 < num8);
        meql[0] = 0;
        int num1 = 1;
        while (num1 <= n)
        {
            num2 = num1 - 1;
            if (xl[num2] > xu[num2])
            {
                return ;
            }
             
            if (xl[num2] == xu[num2])
            {
                int[] numArray2 = new int[]();
                (numArray2 = meql)[0] = numArray2[0] + 1;
            }
             
            num1++;
        }
        int num5 = 0;
        int num6 = n * n;
        num1 = 1;
        while (num1 <= num6)
        {
            num2 = num1 - 1;
            w[z + num2] = 0;
            num1++;
        }
        int num3 = 0;
        for (num1 = 1;num1 <= n;num1++)
        {
            int num4;
            num2 = num1 - 1;
            if (xl[num2] == xu[num2])
            {
                x[num2] = xu[num2];
                num5++;
                w[(u + num5) - 1] = 1;
                iact[num5 - 1] = (num1 + m) + n;
                num4 = num5;
            }
            else
            {
                num4 = (num1 + meql[0]) - num5;
            } 
            w[((z + num3) + num4) - 1] = 1;
            num3 += n;
            w[xbig + num2] = Math.Abs(x[num2]);
        }
        info[0] = 1;
    }

    private void l_l4ong(int n, int m, double[] xl, double[] xu, double[] x, int[] iact, int[] meql, int[] info, double[] w, int z, int u, int xbig, double[] relacc) throws Exception {
        int num2;
        double num7;
        double num8;
        double num9 = 100;
        relacc[0] = 1;
        do
        {
            double[] numArray1 = new double[]();
            (numArray1 = relacc)[0] = numArray1[0] * 0.5;
            num7 = num9 + (0.5 * relacc[0]);
            num8 = num9 + relacc[0];
            if (num9 >= num7)
            {
                break;
            }
             
        }
        while (num7 < num8);
        meql[0] = 0;
        int num1 = 1;
        while (num1 <= n)
        {
            num2 = num1 - 1;
            if (xl[num2] > xu[num2])
            {
                return ;
            }
             
            if (xl[num2] == xu[num2])
            {
                int[] numArray2 = new int[]();
                (numArray2 = meql)[0] = numArray2[0] + 1;
            }
             
            num1++;
        }
        int num5 = 0;
        int num6 = n * n;
        num1 = 1;
        while (num1 <= num6)
        {
            num2 = num1 - 1;
            w[z + num2] = 0;
            num1++;
        }
        int num3 = 0;
        for (num1 = 1;num1 <= n;num1++)
        {
            int num4;
            num2 = num1 - 1;
            if (xl[num2] == xu[num2])
            {
                x[num2] = xu[num2];
                num5++;
                w[(u + num5) - 1] = 1;
                iact[num5 - 1] = (num1 + m) + n;
                num4 = num5;
            }
            else
            {
                num4 = (num1 + meql[0]) - num5;
            } 
            w[((z + num3) + num4) - 1] = 1;
            num3 += n;
            w[xbig + num2] = Math.Abs(x[num2]);
        }
        info[0] = 1;
    }

    private void l_l5onf(int n, int m, int meq, double[] a, int ia, double[] b, double[] xu, int[] iact, int[] meql, int[] info, double[] w, int z, int u, double relacc, int am, int cgrad) throws Exception {
        for (int num7 = 1;num7 <= meq;num7++)
        {
            int num8 = num7 - 1;
            if (meql[0] < n)
            {
                int num9 = meql[0] + 1;
                iact[num9 - 1] = num7;
                this.l_l9onf(n, m, a, ia, iact, meql, w, z, u, relacc, num9, am, cgrad);
                if (meql[0] == num9)
                {
                    goto Label_0194
                }
                 
            }
             
            double num11 = b[num8];
            double num12 = Math.Abs(b[num8]);
            if (meql[0] > 0)
            {
                int num2;
                int num1 = 1;
                while (num1 <= n)
                {
                    num2 = num1 - 1;
                    w[am + num2] = a[(num2 * ia) + num8];
                    num1++;
                }
                int num6 = meql[0];
                do
                {
                    double num10;
                    double[] numArray1 = new double[]();
                    IntPtr ptr1 = new IntPtr();
                    double num13 = 0;
                    int num3 = num6;
                    num1 = 1;
                    while (num1 <= n)
                    {
                        num2 = num1 - 1;
                        num13 += w[(z + num3) - 1] * w[am + num2];
                        num3 += n;
                        num1++;
                    }
                    num13 *= w[(u + num6) - 1];
                    int num4 = iact[num6 - 1];
                    if (num4 <= m)
                    {
                        for (num1 = 1;num1 <= n;num1++)
                        {
                            num2 = num1 - 1;
                            (numArray1 = w)[(int)(ptr1 = (IntPtr)(am + num2))] = numArray1[(int)ptr1] + (-num13 * a[(num2 * ia) + (num4 - 1)]);
                        }
                        num10 = b[num4 - 1];
                    }
                    else
                    {
                        int num5 = (num4 - m) - n;
                        (numArray1 = w)[(int)(ptr1 = (IntPtr)((am + num5) - 1))] = numArray1[(int)ptr1] - num13;
                        num10 = xu[num5 - 1];
                    } 
                    num11 += -num10 * num13;
                    num12 += Math.Abs((double)(num10 * num13));
                    num6--;
                }
                while (num6 >= 1);
            }
             
            if (Math.Abs(num11) > (relacc * num12))
            {
                info[0] = 5;
                return ;
            }
             
            Label_0194:    ;
        }
    }

    private void l_l5ong(int n, int m, int meq, double[] a, int ia, double[] b, double[] xu, int[] iact, int[] meql, int[] info, double[] w, int z, int u, double[] relacc, int am, int cgrad) throws Exception {
        for (int num7 = 1;num7 <= meq;num7++)
        {
            int num8 = num7 - 1;
            if (meql[0] < n)
            {
                int num9 = meql[0] + 1;
                iact[num9 - 1] = num7;
                this.l_l9ong(n, m, a, ia, iact, meql, w, z, u, relacc[0], num9, am, cgrad);
                if (meql[0] == num9)
                {
                    goto Label_0198
                }
                 
            }
             
            double num11 = b[num8];
            double num12 = Math.Abs(b[num8]);
            if (meql[0] > 0)
            {
                int num2;
                int num1 = 1;
                while (num1 <= n)
                {
                    num2 = num1 - 1;
                    w[am + num2] = a[(num2 * ia) + num8];
                    num1++;
                }
                int num6 = meql[0];
                do
                {
                    double num10;
                    double[] numArray1 = new double[]();
                    IntPtr ptr1 = new IntPtr();
                    double num13 = 0;
                    int num3 = num6;
                    num1 = 1;
                    while (num1 <= n)
                    {
                        num2 = num1 - 1;
                        num13 += w[(z + num3) - 1] * w[am + num2];
                        num3 += n;
                        num1++;
                    }
                    num13 *= w[(u + num6) - 1];
                    int num4 = iact[num6 - 1];
                    if (num4 <= m)
                    {
                        for (num1 = 1;num1 <= n;num1++)
                        {
                            num2 = num1 - 1;
                            (numArray1 = w)[(int)(ptr1 = (IntPtr)(am + num2))] = numArray1[(int)ptr1] + (-num13 * a[(num2 * ia) + (num4 - 1)]);
                        }
                        num10 = b[num4 - 1];
                    }
                    else
                    {
                        int num5 = (num4 - m) - n;
                        (numArray1 = w)[(int)(ptr1 = (IntPtr)((am + num5) - 1))] = numArray1[(int)ptr1] - num13;
                        num10 = xu[num5 - 1];
                    } 
                    num11 += -num10 * num13;
                    num12 += Math.Abs((double)(num10 * num13));
                    num6--;
                }
                while (num6 >= 1);
            }
             
            if (Math.Abs(num11) > (relacc[0] * num12))
            {
                info[0] = 5;
                return ;
            }
             
            Label_0198:    ;
        }
    }

    private void l_l6onf(int n, int m, double[] a, int ia, double[] b, double[] xl, double[] xu, double[] x, int[] iact, int[] nact, double[] par, int[] info, double[] w, int g, int z, int u, int xbig, double relacc, double[] tol, int meql, int[] msat, int mtot, int bres, int d, int ztg, int gm, int gmnew, int parnew, int cgrad) throws Exception {
        double[] numArray4 = new double[1];
        double[] numArray1 = numArray4;
        numArray4 = new double[1];
        double[] numArray2 = numArray4;
        int[] numArray5 = new int[1];
        int[] numArray3 = numArray5;
        int num4 = 0;
        double num7 = 0;
        int num3 = 0;
        info[0] = 0;
        int num5 = 1;
        int num6 = 1;
        Label_0043:if (num5 == 1)
        {
            this.l_l10nf(n, m, a, ia, b, xl, xu, x, iact, nact, info, w, z, u, xbig, relacc, tol[0], meql);
            if (info[0] > 0)
            {
                msat[0] = nact[0];
            }
             
            if (msat[0] == mtot)
            {
                return ;
            }
             
        }
         
        if (num6 == 1)
        {
            num4 = msat[0];
            num7 = 0;
        }
         
        this.l_l11nf(n, m, a, ia, b, xl, xu, x, iact, nact, par, w, g, z, u, xbig, bres, d, ztg, relacc, tol[0], numArray1, numArray2, meql, msat, mtot, numArray3, gm, gmnew, parnew, cgrad);
        if (numArray1[0] > 0)
        {
            for (int num1 = 1;num1 <= n;num1++)
            {
                IntPtr ptr1 = new IntPtr();
                int num2 = num1 - 1;
                (numArray4 = x)[(int)(ptr1 = (IntPtr)num2)] = numArray4[(int)ptr1] + (numArray1[0] * w[d + num2]);
                w[xbig + num2] = Math.Max(w[xbig + num2], Math.Abs(x[num2]));
            }
            this.l_l9onf(n, m, a, ia, iact, nact, w, z, u, relacc, numArray3[0], gmnew, cgrad);
        }
         
        if (msat[0] < mtot)
        {
            if (numArray1[0] == 0)
            {
                if (tol[0] <= relacc)
                {
                    return ;
                }
                 
                this.l_l8onf(n, m, a, ia, b, xl, xu, x, iact, nact[0], w, xbig, relacc, tol, meql);
                num5 = 1;
                num6 = 1;
                goto Label_0043
            }
             
            if (num4 < msat[0])
            {
                num5 = 0;
                num6 = 1;
                goto Label_0043
            }
             
            if ((num7 == 0) || (numArray2[0] < num7))
            {
                num7 = numArray2[0];
                num3 = 0;
            }
             
            num3++;
            if (num3 <= 2)
            {
                num5 = 0;
                num6 = 0;
                goto Label_0043
            }
             
            if (tol[0] > relacc)
            {
                this.l_l8onf(n, m, a, ia, b, xl, xu, x, iact, nact[0], w, xbig, relacc, tol, meql);
                num5 = 1;
                num6 = 1;
                goto Label_0043
            }
             
        }
         
    }

    private void l_l6ong(int n, int m, double[] a, int ia, double[] b, double[] xl, double[] xu, double[] x, int[] iact, int[] nact, double[] par, int[] info, double[] w, int g, int z, int u, int xbig, double relacc, double[] tol, int meql, int[] msat, int mtot, int bres, int d, int ztg, int gm, int gmnew, int parnew, int cgrad) throws Exception {
        double[] numArray4 = new double[1];
        double[] numArray1 = numArray4;
        numArray4 = new double[1];
        double[] numArray2 = numArray4;
        int[] numArray5 = new int[1];
        int[] numArray3 = numArray5;
        int num4 = 0;
        double num7 = 0;
        int num3 = 0;
        info[0] = 0;
        int num5 = 1;
        int num6 = 1;
        Label_0043:if (num5 == 1)
        {
            this.l_l10ng(n, m, a, ia, b, xl, xu, x, iact, nact, info, w, z, u, xbig, relacc, tol[0], meql);
            if (info[0] > 0)
            {
                msat[0] = nact[0];
            }
             
            if (msat[0] == mtot)
            {
                return ;
            }
             
        }
         
        if (num6 == 1)
        {
            num4 = msat[0];
            num7 = 0;
        }
         
        this.l_l11ng(n, m, a, ia, b, xl, xu, x, iact, nact, par, w, g, z, u, xbig, bres, d, ztg, relacc, tol[0], numArray1, numArray2, meql, msat, mtot, numArray3, gm, gmnew, parnew, cgrad);
        if (numArray1[0] > 0)
        {
            for (int num1 = 1;num1 <= n;num1++)
            {
                IntPtr ptr1 = new IntPtr();
                int num2 = num1 - 1;
                (numArray4 = x)[(int)(ptr1 = (IntPtr)num2)] = numArray4[(int)ptr1] + (numArray1[0] * w[d + num2]);
                w[xbig + num2] = Math.Max(w[xbig + num2], Math.Abs(x[num2]));
            }
            this.l_l9ong(n, m, a, ia, iact, nact, w, z, u, relacc, numArray3[0], gmnew, cgrad);
        }
         
        if (msat[0] < mtot)
        {
            if (numArray1[0] == 0)
            {
                if (tol[0] <= relacc)
                {
                    return ;
                }
                 
                this.l_l8ong(n, m, a, ia, b, xl, xu, x, iact, nact[0], w, xbig, relacc, tol, meql);
                num5 = 1;
                num6 = 1;
                goto Label_0043
            }
             
            if (num4 < msat[0])
            {
                num5 = 0;
                num6 = 1;
                goto Label_0043
            }
             
            if ((num7 == 0) || (numArray2[0] < num7))
            {
                num7 = numArray2[0];
                num3 = 0;
            }
             
            num3++;
            if (num3 <= 2)
            {
                num5 = 0;
                num6 = 0;
                goto Label_0043
            }
             
            if (tol[0] > relacc)
            {
                this.l_l8ong(n, m, a, ia, b, xl, xu, x, iact, nact[0], w, xbig, relacc, tol, meql);
                num5 = 1;
                num6 = 1;
                goto Label_0043
            }
             
        }
         
    }

    private void l_l7onf(Provisdom.Optimization.NlpLinearConstraintSolver.IFunction fcn, int n, int m, double[] a, int ia, double[] b, double[] xl, double[] xu, double[] x, double[] obj, double acc, int[] iact, int[] nact, double[] par, int iprint, int[] info, double[] w, int g, int z, int u, int xbig, double relacc, double[] zznorm, double tol, int meql, int mtot, int[] iterc, int[] nfvals, int nfmax, int reskt, int bres, int d, int ztg, int gm, int xs, int gs) throws Exception {
        double[] numArray9 = new double[1];
        double[] numArray1 = numArray9;
        numArray9 = new double[1];
        double[] numArray2 = numArray9;
        numArray9 = new double[1];
        double[] numArray3 = numArray9;
        numArray9 = new double[1];
        double[] numArray4 = numArray9;
        numArray9 = new double[1];
        double[] numArray5 = numArray9;
        numArray9 = new double[1];
        double[] numArray6 = numArray9;
        int[] numArray10 = new int[1];
        int[] numArray7 = numArray10;
        numArray10 = new int[1];
        int[] numArray8 = numArray10;
        numArray7[0] = mtot;
        int num3 = iterc[0];
        int num6 = nfvals[0];
        numArray1[0] = 0;
        double num8 = 0;
        if ((nfvals[0] == 0) || (info[0] == 1))
        {
            numArray1[0] = fcn.F(x);
            obj[0] = numArray1[0];
            this.l_l21nf(fcn, n, x, numArray1[0], w, g, nfvals);
            (numArray10 = nfvals)[0] = numArray10[0] + 1;
        }
         
        double num9 = Math.Abs((double)((numArray1[0] + numArray1[0]) + 1));
        int num4 = -1;
        Label_00EC:this.l_l11nf(n, m, a, ia, b, xl, xu, x, iact, nact, par, w, g, z, u, xbig, bres, d, ztg, relacc, tol, numArray2, numArray3, meql, numArray7, mtot, numArray8, gm, reskt, xs, gs);
        this.l_l12nf(n, m, a, ia, iact, nact[0], par, w, g, reskt, z, u, bres, numArray4, meql, numArray5, xs, gs);
        if (numArray5[0] <= (acc * acc))
        {
            info[0] = 1;
            if (iprint != 0)
            {
                num4 = -1;
            }
             
        }
        else if (numArray3[0] >= 0)
        {
            info[0] = 2;
            if (iprint != 0)
            {
                num4 = -1;
            }
             
        }
        else
        {
            int num7 = 0;
            if (numArray1[0] >= num9)
            {
                if (((tol == relacc) || (nact[0] == 0)) && (num8 > 0))
                {
                    num7 = 1;
                }
                 
                if (num7 == 0)
                {
                    info[0] = 3;
                    if (iprint != 0)
                    {
                        num4 = -1;
                    }
                     
                    return ;
                }
                 
            }
             
            num8 = num9 - numArray1[0];
            num9 = numArray1[0];
            if (nfvals[0] == nfmax)
            {
                info[0] = 8;
                if (iprint != 0)
                {
                    num4 = -1;
                }
                 
            }
            else if (((tol > relacc) && (iterc[0] > num3)) && ((0.1 * numArray4[0]) >= Math.Max(num8, -0.5 * numArray3[0])))
            {
                if (iprint != 0)
                {
                    num4 = -1;
                }
                 
            }
            else if (num4 == iterc[0])
            {
                num4 = iterc[0] + Math.Abs(iprint);
            }
            else
            {
                int num1;
                int num2;
                (numArray10 = iterc)[0] = numArray10[0] + 1;
                this.l_l13nf(fcn, n, x, w, g, d, xs, gs, relacc, numArray2[0], numArray3[0], numArray1, numArray6, nfvals, nfmax, bres, obj);
                if (numArray6[0] == 0)
                {
                    info[0] = 3;
                    double num10 = 0;
                    for (num1 = 1;num1 <= n;num1++)
                    {
                        num2 = num1 - 1;
                        num10 += Math.Abs((double)(w[d + num2] * w[gs + num2]));
                    }
                    if ((numArray3[0] + (relacc * num10)) >= 0)
                    {
                        info[0] = 2;
                    }
                     
                    if (iprint != 0)
                    {
                        num4 = -1;
                    }
                     
                }
                else
                {
                    for (num1 = 1;num1 <= n;num1++)
                    {
                        num2 = num1 - 1;
                        w[xbig + num2] = Math.Max(w[xbig + num2], Math.Abs(x[num2]));
                    }
                    this.l_l14nf(n, x, nact[0], w, g, z, ztg, xs, gs, zznorm);
                    if (numArray6[0] == numArray2[0])
                    {
                        int num5 = iact[numArray8[0] - 1];
                        if (num5 > m)
                        {
                            num5 -= m;
                            if (num5 <= n)
                            {
                                x[num5 - 1] = xl[num5 - 1];
                            }
                            else
                            {
                                x[(num5 - n) - 1] = xu[(num5 - n) - 1];
                            } 
                        }
                         
                        this.l_l9onf(n, m, a, ia, iact, nact, w, z, u, relacc, numArray8[0], xs, gs);
                    }
                     
                    goto Label_00EC
                } 
            }   
        }  
    }

    private void l_l7ong(Provisdom.Optimization.NlpLinearConstraintSolver.IFunction fcn, Provisdom.Optimization.NlpLinearConstraintSolver.IGradient grad, int n, int m, double[] a, int ia, double[] b, double[] xl, double[] xu, double[] x, double[] obj, double acc, int[] iact, int[] nact, double[] par, int iprint, int[] info, double[] w, int g, int z, int u, int xbig, double relacc, double[] zznorm, double tol, int meql, int mtot, int[] iterc, int[] nfvals, int nfmax, int reskt, int bres, int d, int ztg, int gm, int xs, int gs) throws Exception {
        double[] numArray10 = new double[1];
        double[] numArray1 = numArray10;
        numArray10 = new double[1];
        double[] numArray2 = numArray10;
        numArray10 = new double[1];
        double[] numArray3 = numArray10;
        numArray10 = new double[1];
        double[] numArray4 = numArray10;
        numArray10 = new double[1];
        double[] numArray5 = numArray10;
        numArray10 = new double[1];
        double[] numArray6 = numArray10;
        int[] numArray11 = new int[1];
        int[] numArray7 = numArray11;
        numArray11 = new int[1];
        int[] numArray8 = numArray11;
        numArray7[0] = mtot;
        int num3 = iterc[0];
        int num6 = nfvals[0];
        numArray1[0] = 0;
        double num8 = 0;
        if ((nfvals[0] == 0) || (info[0] == 1))
        {
            numArray1[0] = fcn.F(x);
            obj[0] = numArray1[0];
            double[] numArray9 = new double[n];
            BLAS.copy(n, w, g, 1, numArray9, 0, 1);
            grad.Gradient(x, numArray9);
            BLAS.copy(n, numArray9, 0, 1, w, g, 1);
            (numArray11 = nfvals)[0] = numArray11[0] + 1;
        }
         
        double num9 = Math.Abs((double)((numArray1[0] + numArray1[0]) + 1));
        int num4 = -1;
        Label_0108:this.l_l11ng(n, m, a, ia, b, xl, xu, x, iact, nact, par, w, g, z, u, xbig, bres, d, ztg, relacc, tol, numArray2, numArray3, meql, numArray7, mtot, numArray8, gm, reskt, xs, gs);
        this.l_l12ng(n, m, a, ia, iact, nact[0], par, w, g, reskt, z, u, bres, numArray4, meql, numArray5, xs, gs);
        if (numArray5[0] <= (acc * acc))
        {
            info[0] = 1;
            if (iprint != 0)
            {
                num4 = -1;
            }
             
        }
        else if (numArray3[0] >= 0)
        {
            info[0] = 2;
            if (iprint != 0)
            {
                num4 = -1;
            }
             
        }
        else
        {
            int num7 = 0;
            if (numArray1[0] >= num9)
            {
                if (((tol == relacc) || (nact[0] == 0)) && (num8 > 0))
                {
                    num7 = 1;
                }
                 
                if (num7 == 0)
                {
                    info[0] = 3;
                    if (iprint != 0)
                    {
                        num4 = -1;
                    }
                     
                    return ;
                }
                 
            }
             
            num8 = num9 - numArray1[0];
            num9 = numArray1[0];
            if (nfvals[0] == nfmax)
            {
                info[0] = 8;
                if (iprint != 0)
                {
                    num4 = -1;
                }
                 
            }
            else if (((tol > relacc) && (iterc[0] > num3)) && ((0.1 * numArray4[0]) >= Math.Max(num8, -0.5 * numArray3[0])))
            {
                if (iprint != 0)
                {
                    num4 = -1;
                }
                 
            }
            else if (num4 == iterc[0])
            {
                num4 = iterc[0] + Math.Abs(iprint);
            }
            else
            {
                int num1;
                int num2;
                (numArray11 = iterc)[0] = numArray11[0] + 1;
                this.l_l13ng(fcn, grad, n, x, w, g, d, xs, gs, relacc, numArray2[0], numArray3[0], numArray1, numArray6, nfvals, nfmax, bres, obj);
                if (numArray6[0] == 0)
                {
                    info[0] = 3;
                    double num10 = 0;
                    for (num1 = 1;num1 <= n;num1++)
                    {
                        num2 = num1 - 1;
                        num10 += Math.Abs((double)(w[d + num2] * w[gs + num2]));
                    }
                    if ((numArray3[0] + (relacc * num10)) >= 0)
                    {
                        info[0] = 2;
                    }
                     
                    if (iprint != 0)
                    {
                        num4 = -1;
                    }
                     
                }
                else
                {
                    for (num1 = 1;num1 <= n;num1++)
                    {
                        num2 = num1 - 1;
                        w[xbig + num2] = Math.Max(w[xbig + num2], Math.Abs(x[num2]));
                    }
                    this.l_l14ng(n, x, nact[0], w, g, z, ztg, xs, gs, zznorm);
                    if (numArray6[0] == numArray2[0])
                    {
                        int num5 = iact[numArray8[0] - 1];
                        if (num5 > m)
                        {
                            num5 -= m;
                            if (num5 <= n)
                            {
                                x[num5 - 1] = xl[num5 - 1];
                            }
                            else
                            {
                                x[(num5 - n) - 1] = xu[(num5 - n) - 1];
                            } 
                        }
                         
                        this.l_l9ong(n, m, a, ia, iact, nact, w, z, u, relacc, numArray8[0], xs, gs);
                    }
                     
                    goto Label_0108
                } 
            }   
        }  
    }

    private void l_l8onf(int n, int m, double[] a, int ia, double[] b, double[] xl, double[] xu, double[] x, int[] iact, int nact, double[] w, int xbig, double relacc, double[] tol, int meql) throws Exception {
        int num1;
        int num2;
        double num10 = 0;
        if (nact > meql)
        {
            int num7 = meql + 1;
            for (int num5 = num7;num5 <= nact;num5++)
            {
                double num8;
                double num9;
                int num6 = num5 - 1;
                int num3 = iact[num6];
                if (num3 <= m)
                {
                    num8 = b[num3 - 1];
                    num9 = Math.Abs(b[num3 - 1]);
                    for (num1 = 1;num1 <= n;num1++)
                    {
                        num2 = num1 - 1;
                        num8 += -a[(num2 * ia) + (num3 - 1)] * x[num2];
                        num9 += Math.Abs((double)(a[(num2 * ia) + (num3 - 1)] * w[xbig + num2]));
                    }
                }
                else
                {
                    int num4 = num3 - m;
                    if (num4 <= n)
                    {
                        num8 = x[num4 - 1] - xl[num4 - 1];
                        num9 = w[(xbig + num4) - 1] + Math.Abs(xl[num4 - 1]);
                    }
                    else
                    {
                        num4 -= n;
                        num8 = xu[num4 - 1] - x[num4 - 1];
                        num9 = w[(xbig + num4) - 1] + Math.Abs(xu[num4 - 1]);
                    } 
                } 
                if (num8 > 0)
                {
                    num10 = Math.Max(num10, num8 / num9);
                }
                 
            }
        }
         
        tol[0] = 0.1 * Math.Min(tol[0], num10);
        if (tol[0] <= (relacc + relacc))
        {
            tol[0] = relacc;
            for (num1 = 1;num1 <= n;num1++)
            {
                num2 = num1 - 1;
                w[xbig + num2] = Math.Abs(x[num2]);
            }
        }
         
    }

    private void l_l8ong(int n, int m, double[] a, int ia, double[] b, double[] xl, double[] xu, double[] x, int[] iact, int nact, double[] w, int xbig, double relacc, double[] tol, int meql) throws Exception {
        int num1;
        int num2;
        double num10 = 0;
        if (nact > meql)
        {
            int num7 = meql + 1;
            for (int num5 = num7;num5 <= nact;num5++)
            {
                double num8;
                double num9;
                int num6 = num5 - 1;
                int num3 = iact[num6];
                if (num3 <= m)
                {
                    num8 = b[num3 - 1];
                    num9 = Math.Abs(b[num3 - 1]);
                    for (num1 = 1;num1 <= n;num1++)
                    {
                        num2 = num1 - 1;
                        num8 += -a[(num2 * ia) + (num3 - 1)] * x[num2];
                        num9 += Math.Abs((double)(a[(num2 * ia) + (num3 - 1)] * w[xbig + num2]));
                    }
                }
                else
                {
                    int num4 = num3 - m;
                    if (num4 <= n)
                    {
                        num8 = x[num4 - 1] - xl[num4 - 1];
                        num9 = w[(xbig + num4) - 1] + Math.Abs(xl[num4 - 1]);
                    }
                    else
                    {
                        num4 -= n;
                        num8 = xu[num4 - 1] - x[num4 - 1];
                        num9 = w[(xbig + num4) - 1] + Math.Abs(xu[num4 - 1]);
                    } 
                } 
                if (num8 > 0)
                {
                    num10 = Math.Max(num10, num8 / num9);
                }
                 
            }
        }
         
        tol[0] = 0.1 * Math.Min(tol[0], num10);
        if (tol[0] <= (relacc + relacc))
        {
            tol[0] = relacc;
            for (num1 = 1;num1 <= n;num1++)
            {
                num2 = num1 - 1;
                w[xbig + num2] = Math.Abs(x[num2]);
            }
        }
         
    }

    private void l_l9onf(int n, int m, double[] a, int ia, int[] iact, int[] nact, double[] w, int z, int u, double relacc, int indxbd, int ztc, int cgrad) throws Exception {
        int num1;
        int num2;
        int num6;
        int num8 = new int();
        int num9;
        int num10 = new int();
        double num12 = new double();
        double num14 = new double();
        double[] numArray1 = new double[]();
        IntPtr ptr1 = new IntPtr();
        int num5 = 0;
        int num7 = 0;
        int num11 = nact[0] + 1;
        int num3 = iact[indxbd - 1];
        iact[indxbd - 1] = iact[num11 - 1];
        iact[num11 - 1] = num3;
        if (num3 > m)
        {
            int num4 = num3 - m;
            if (num4 <= n)
            {
                num14 = -1;
            }
            else
            {
                num4 -= n;
                num14 = 1;
            } 
            num7 = (num4 * n) - n;
            for (num8 = 1;num8 <= n;num8++)
            {
                num9 = num8 - 1;
                w[ztc + num9] = num14 * w[(z + num7) + num9];
            }
        }
        else
        {
            num1 = 1;
            while (num1 <= n)
            {
                num2 = num1 - 1;
                w[cgrad + num2] = a[(num2 * ia) + (num3 - 1)];
                num1++;
            }
            for (num8 = 1;num8 <= n;num8++)
            {
                num9 = num8 - 1;
                w[ztc + num9] = 0;
                num6 = num8;
                for (num1 = 1;num1 <= n;num1++)
                {
                    num2 = num1 - 1;
                    (numArray1 = w)[(int)(ptr1 = (IntPtr)(ztc + num9))] = numArray1[(int)ptr1] + (w[(z + num6) - 1] * w[cgrad + num2]);
                    num6 += n;
                }
            }
        } 
        num8 = n;
        Label_011C:num10 = num8;
        num8--;
        if (num8 > nact[0])
        {
            if (w[(ztc + num10) - 1] != 0)
            {
                if (Math.Abs(w[(ztc + num10) - 1]) <= (relacc * Math.Abs(w[(ztc + num8) - 1])))
                {
                    num14 = Math.Abs(w[(ztc + num8) - 1]);
                }
                else if (Math.Abs(w[(ztc + num8) - 1]) <= (relacc * Math.Abs(w[(ztc + num10) - 1])))
                {
                    num14 = Math.Abs(w[(ztc + num10) - 1]);
                }
                else
                {
                    num14 = Math.Abs(w[(ztc + num10) - 1]) * Math.Sqrt(1 + Math.Pow(w[(ztc + num8) - 1] / w[(ztc + num10) - 1], 2));
                }  
                double num17 = w[(ztc + num8) - 1] / num14;
                double num19 = w[(ztc + num10) - 1] / num14;
                w[(ztc + num8) - 1] = num14;
                num6 = num8;
                if (num3 > m)
                {
                    for (num1 = 1;num1 <= n;num1++)
                    {
                        num2 = num1 - 1;
                        num14 = (num17 * w[z + num6]) - (num19 * w[(z + num6) - 1]);
                        w[(z + num6) - 1] = (num17 * w[(z + num6) - 1]) + (num19 * w[z + num6]);
                        w[z + num6] = num14;
                        num6 += n;
                    }
                    w[((z + num7) + num10) - 1] = 0;
                }
                else
                {
                    double num18 = 0;
                    num1 = 1;
                    while (num1 <= n)
                    {
                        num2 = num1 - 1;
                        double num15 = num17 * w[z + num6];
                        double num16 = num19 * w[(z + num6) - 1];
                        num14 = Math.Abs(w[cgrad + num2]) * (Math.Abs(num15) + Math.Abs(num16));
                        if (num14 > num18)
                        {
                            num18 = num14;
                            num5 = num1;
                        }
                         
                        w[(z + num6) - 1] = (num17 * w[(z + num6) - 1]) + (num19 * w[z + num6]);
                        w[z + num6] = num15 - num16;
                        num6 += n;
                        num1++;
                    }
                    num12 = 0;
                    num6 = num10;
                    for (num1 = 1;num1 <= n;num1++)
                    {
                        num2 = num1 - 1;
                        num12 += w[(z + num6) - 1] * w[cgrad + num2];
                        num6 += n;
                    }
                    if (num12 != 0)
                    {
                        num6 = ((num5 * n) - n) + num10;
                        (numArray1 = w)[(int)(ptr1 = (IntPtr)((z + num6) - 1))] = numArray1[(int)ptr1] + (-num12 / w[(cgrad + num5) - 1]);
                    }
                     
                } 
            }
             
            goto Label_011C
        }
         
        if (w[(ztc + num11) - 1] != 0)
        {
            if (num3 <= m)
            {
                num12 = 0;
                double num13 = 0;
                num6 = num11;
                for (num1 = 1;num1 <= n;num1++)
                {
                    num2 = num1 - 1;
                    num14 = w[(z + num6) - 1] * w[cgrad + num2];
                    num12 += num14;
                    num13 += Math.Abs(num14);
                    num6 += n;
                }
                if (Math.Abs(num12) <= (relacc * num13))
                {
                    return ;
                }
                 
            }
             
            w[(u + num11) - 1] = 1 / w[(ztc + num11) - 1];
            nact[0] = num11;
        }
         
    }

    private void l_l9ong(int n, int m, double[] a, int ia, int[] iact, int[] nact, double[] w, int z, int u, double relacc, int indxbd, int ztc, int cgrad) throws Exception {
        int num1;
        int num2;
        int num6;
        int num8 = new int();
        int num9;
        int num10 = new int();
        double num12 = new double();
        double num14 = new double();
        double[] numArray1 = new double[]();
        IntPtr ptr1 = new IntPtr();
        int num5 = 0;
        int num7 = 0;
        int num11 = nact[0] + 1;
        int num3 = iact[indxbd - 1];
        iact[indxbd - 1] = iact[num11 - 1];
        iact[num11 - 1] = num3;
        if (num3 > m)
        {
            int num4 = num3 - m;
            if (num4 <= n)
            {
                num14 = -1;
            }
            else
            {
                num4 -= n;
                num14 = 1;
            } 
            num7 = (num4 * n) - n;
            for (num8 = 1;num8 <= n;num8++)
            {
                num9 = num8 - 1;
                w[ztc + num9] = num14 * w[(z + num7) + num9];
            }
        }
        else
        {
            num1 = 1;
            while (num1 <= n)
            {
                num2 = num1 - 1;
                w[cgrad + num2] = a[(num2 * ia) + (num3 - 1)];
                num1++;
            }
            for (num8 = 1;num8 <= n;num8++)
            {
                num9 = num8 - 1;
                w[ztc + num9] = 0;
                num6 = num8;
                for (num1 = 1;num1 <= n;num1++)
                {
                    num2 = num1 - 1;
                    (numArray1 = w)[(int)(ptr1 = (IntPtr)(ztc + num9))] = numArray1[(int)ptr1] + (w[(z + num6) - 1] * w[cgrad + num2]);
                    num6 += n;
                }
            }
        } 
        num8 = n;
        Label_011C:num10 = num8;
        num8--;
        if (num8 > nact[0])
        {
            if (w[(ztc + num10) - 1] != 0)
            {
                if (Math.Abs(w[(ztc + num10) - 1]) <= (relacc * Math.Abs(w[(ztc + num8) - 1])))
                {
                    num14 = Math.Abs(w[(ztc + num8) - 1]);
                }
                else if (Math.Abs(w[(ztc + num8) - 1]) <= (relacc * Math.Abs(w[(ztc + num10) - 1])))
                {
                    num14 = Math.Abs(w[(ztc + num10) - 1]);
                }
                else
                {
                    num14 = Math.Abs(w[(ztc + num10) - 1]) * Math.Sqrt(1 + Math.Pow(w[(ztc + num8) - 1] / w[(ztc + num10) - 1], 2));
                }  
                double num17 = w[(ztc + num8) - 1] / num14;
                double num19 = w[(ztc + num10) - 1] / num14;
                w[(ztc + num8) - 1] = num14;
                num6 = num8;
                if (num3 > m)
                {
                    for (num1 = 1;num1 <= n;num1++)
                    {
                        num2 = num1 - 1;
                        num14 = (num17 * w[z + num6]) - (num19 * w[(z + num6) - 1]);
                        w[(z + num6) - 1] = (num17 * w[(z + num6) - 1]) + (num19 * w[z + num6]);
                        w[z + num6] = num14;
                        num6 += n;
                    }
                    w[((z + num7) + num10) - 1] = 0;
                }
                else
                {
                    double num18 = 0;
                    num1 = 1;
                    while (num1 <= n)
                    {
                        num2 = num1 - 1;
                        double num15 = num17 * w[z + num6];
                        double num16 = num19 * w[(z + num6) - 1];
                        num14 = Math.Abs(w[cgrad + num2]) * (Math.Abs(num15) + Math.Abs(num16));
                        if (num14 > num18)
                        {
                            num18 = num14;
                            num5 = num1;
                        }
                         
                        w[(z + num6) - 1] = (num17 * w[(z + num6) - 1]) + (num19 * w[z + num6]);
                        w[z + num6] = num15 - num16;
                        num6 += n;
                        num1++;
                    }
                    num12 = 0;
                    num6 = num10;
                    for (num1 = 1;num1 <= n;num1++)
                    {
                        num2 = num1 - 1;
                        num12 += w[(z + num6) - 1] * w[cgrad + num2];
                        num6 += n;
                    }
                    if (num12 != 0)
                    {
                        num6 = ((num5 * n) - n) + num10;
                        (numArray1 = w)[(int)(ptr1 = (IntPtr)((z + num6) - 1))] = numArray1[(int)ptr1] + (-num12 / w[(cgrad + num5) - 1]);
                    }
                     
                } 
            }
             
            goto Label_011C
        }
         
        if (w[(ztc + num11) - 1] != 0)
        {
            if (num3 <= m)
            {
                num12 = 0;
                double num13 = 0;
                num6 = num11;
                for (num1 = 1;num1 <= n;num1++)
                {
                    num2 = num1 - 1;
                    num14 = w[(z + num6) - 1] * w[cgrad + num2];
                    num12 += num14;
                    num13 += Math.Abs(num14);
                    num6 += n;
                }
                if (Math.Abs(num12) <= (relacc * num13))
                {
                    return ;
                }
                 
            }
             
            w[(u + num11) - 1] = 1 / w[(ztc + num11) - 1];
            nact[0] = num11;
        }
         
    }

    private void l_m1ran(int nra, int nca, double[] a, double[] atran) throws Exception {
        int num2;
        int num3;
        int num9 = 0;
        int num4 = (nra * nca) - 1;
        if (nra != nca)
        {
            num9 = 1;
        }
         
        if (num9 == 0)
        {
            if (atran != a)
            {
                BLAS.copy(nra * nca, a, atran);
            }
             
            int num5 = 2;
            int num6 = nra + 1;
            int num7 = nra - 1;
            for (num2 = nra;num2 <= num4;num2 += nra)
            {
                int num8 = num5 + num7;
                for (num3 = num5;num3 <= num2;num3++)
                {
                    double num1 = atran[num3 - 1];
                    atran[num3 - 1] = atran[num8 - 1];
                    atran[num8 - 1] = num1;
                    num8 += nra;
                }
                num5 += num6;
            }
        }
        else
        {
            double[] numArray1 = new double[]();
            if (a == atran)
            {
                numArray1 = new double[nca * nra];
            }
            else
            {
                numArray1 = atran;
            } 
            for (num2 = 0;num2 < nra;num2++)
            {
                for (num3 = 0;num3 < nca;num3++)
                {
                    numArray1[num2 + (nra * num3)] = a[num3 + (nca * num2)];
                }
            }
            if (a == atran)
            {
                BLAS.copy(nra * nca, numArray1, atran);
            }
             
        } 
    }

    public double[] getGuess() throws Exception {
        return (double[])_guess.Clone();
    }

    public void setGuess(double[] value) throws Exception {
        if (value.Length != this._var)
        {
            Object[] objArray1 = new Object[]{ "guess", value.Length, this._var };
            ExceptionThrower.ThrowArgumentException("Numerical.Math", "NotEqual", objArray1);
        }
         
        this._guess = (double[])value.Clone();
    }

    public void solve() throws Exception {
        int num1 = 0;
        double[] numArray1 = this._guess;
        int[] numArray11 = new int[1];
        int[] numArray3 = numArray11;
        double[] numArray12 = new double[1];
        double[] numArray4 = numArray12;
        double num4 = this._tol;
        numArray11 = new int[1];
        int[] numArray5 = numArray11;
        int[] numArray6 = new int[]{ this._maxObjective };
        int num5 = this._var;
        int num6 = this._con;
        int num7 = this._eq;
        double[] numArray7 = this._a;
        double[] numArray8 = this._b;
        double[] numArray9 = this._lowerBounds;
        double[] numArray10 = this._upperBounds;
        Provisdom.Optimization.NlpLinearConstraintSolver.IFunction function1 = this._objective;
        Provisdom.Optimization.NlpLinearConstraintSolver.IGradient gradient1 = this._gradient;
        int num2 = ((num5 * num5) + (11 * num5)) + num6;
        double[] numArray2 = new double[num2];
        this._value = new double[num5];
        double num3 = num4;
        if (num6 > 1)
        {
            this.l_m1ran(num6, num5, numArray7, numArray7);
        }
         
        if (this._gradient != null)
        {
            this.l_l2ong(function1, gradient1, num5, num6, num7, numArray7, num6, numArray8, numArray9, numArray10, numArray1, num3, numArray6, this._value, numArray4, numArray5, this._iactUser, this._alamdaUser, num1, numArray3, numArray2);
        }
        else
        {
            this.l_l2onf(function1, num5, num6, num7, numArray7, num6, numArray8, numArray9, numArray10, numArray1, num3, numArray6, this._value, numArray4, numArray5, this._iactUser, this._alamdaUser, num1, numArray3, numArray2);
        } 
        this._act = numArray5[0];
        this._obj = numArray4[0];
        if (num6 > 1)
        {
            this.l_m1ran(num5, num6, numArray7, numArray7);
        }
         
    }

    public int getFinalActiveConstraintsNum() throws Exception {
        return this._act;
    }

    public double getObjectiveValue() throws Exception {
        return this._obj;
    }

    public double getTolerance() throws Exception {
        return this._tol;
    }

    public void setTolerance(double value) throws Exception {
        if (value <= 0)
        {
            Object[] objArray1 = new Object[]{ "tolerance", value };
            ExceptionThrower.ThrowArgumentException("Numerical.Math", "NotPositive", objArray1);
        }
         
        this._tol = value;
    }

    private double[] _alamdaUser = new double[]();
    private int[] _iactUser = new int[]();
    private double[] _a = new double[]();
    private double[] _b = new double[]();
    private Provisdom.Optimization.NlpLinearConstraintSolver.IFunction _objective;
    private Provisdom.Optimization.NlpLinearConstraintSolver.IGradient _gradient;
    private int _maxObjective = new int();
    private int _act = new int();
    private int _con = new int();
    private int _eq = new int();
    private int _var = new int();
    private double _obj = new double();
    private double _tol = new double();
    private double[] _guess = new double[]();
    private double[] _lowerBounds = new double[]();
    private double[] _upperBounds = new double[]();
    private double[] _value = new double[]();
    public interface IFunction   
    {
        double f(double[] x) throws Exception ;
    
    }

    public interface IGradient   extends Provisdom.Optimization.NlpLinearConstraintSolver.IFunction
    {
        void gradient(double[] x, double[] g) throws Exception ;
    
    }

}


