//
// Translated by CS2J (http://www.cs2j.com): 2/18/2016 3:06:45 PM
//

package java.Provisdom.LinearAlgebra;

public class BLAS
{
    public static void add(int n, double sa, double[] sx, int offset_sx, int incx) throws Exception {
        if (n > 0)
        {
            int num1;
            if (incx != 1)
            {
                int num2 = offset_sx;
                for (num1 = 0;num1 < n;num1++)
                {
                    sx[num2] += sa;
                    num2 += incx;
                }
            }
            else
            {
                for (num1 = offset_sx;num1 < offset_sx + n;num1++)
                {
                    sx[num1] += sa;
                }
            } 
        }
         
    }

    public static void add(int n, double sa, double[][] sx, int sxRow, int offset_sx, int incx) throws Exception {
        if (n > 0)
        {
            int num1;
            if (incx != 1)
            {
                int num2 = offset_sx;
                for (num1 = 1;num1 <= n;num1++)
                {
                    sx[sxRow][num2] += sa;
                    num2 += incx;
                }
            }
            else
            {
                for (num1 = offset_sx;num1 < offset_sx + n;num1++)
                {
                    sx[sxRow][num1] += sa;
                }
            } 
        }
         
    }

    public static double asum(int n, double[] sx, int offset_sx, int incx) throws Exception {
        double num5 = 0;
        if (n > 0)
        {
            int num3;
            if (incx != 1)
            {
                int num4 = n * incx;
                num3 = 1;
                int num2 = incx;
                for (int num1 = ((((num4 - 1) + num2) / num2) > 0) ? (((num4 - 1) + num2) / num2) : 0;num1 > 0;num1--)
                {
                    num5 += Math.abs(sx[(offset_sx + num3) - 1]);
                    num3 += num2;
                }
                return num5;
            }
             
            for (num3 = 1;num3 <= n;num3++)
            {
                num5 += Math.abs(sx[(offset_sx + num3) - 1]);
            }
        }
         
        return num5;
    }

    public static void axpy(int n, double a, double[] x, double[][] y, int yRow) throws Exception {
        for (int num1 = 0;num1 < n;num1++)
        {
            y[yRow][num1] += a * x[num1];
        }
    }

    public static void axpy(int n, double a, double[] x, int xoffset, double[] y, int yoffset) throws Exception {
        for (int num1 = 0;num1 < n;num1++)
        {
            y[num1 + yoffset] += a * x[num1 + xoffset];
        }
    }

    public static void axpy(int n, double a, double[] x, int xoffset, double[][] y, int yRow, int yoffset) throws Exception {
        for (int num1 = 0;num1 < n;num1++)
        {
            y[yRow][num1 + yoffset] += a * x[num1 + xoffset];
        }
    }

    public static void axpy(int n, double a, double[][] x, int xRow, int xoffset, double[] y, int yoffset) throws Exception {
        for (int num1 = 0;num1 < n;num1++)
        {
            y[num1 + yoffset] += a * x[xRow][num1 + xoffset];
        }
    }

    public static void axpy(int n, double da, double[] dx, int dx_off, int incx, double[] dy, int dy_off, int incy) throws Exception {
        if ((n > 0) && (da != 0))
        {
            if ((incx != 1) || (incy != 1))
            {
                int num1 = dx_off;
                int num2 = dy_off;
                if (incx < 0)
                {
                    num1 = ((-n + 1) * incx) + dx_off;
                }
                 
                if (incy < 0)
                {
                    num2 = ((-n + 1) * incy) + dy_off;
                }
                 
                for (int num4 = 0;num4 < n;num4++)
                {
                    dy[num2] += da * dx[num1];
                    num1 += incx;
                    num2 += incy;
                }
            }
            else
            {
                int num3 = n % 4;
                for (int num5 = 0;num5 < num3;num5++)
                {
                    dy[num5 + dy_off] += da * dx[num5 + dx_off];
                }
                for (int num6 = num3;num6 < n;num6 += 4)
                {
                    int num7 = num6 + dx_off;
                    int num8 = num6 + dy_off;
                    dy[num8] += da * dx[num7];
                    dy[num8 + 1] += da * dx[num7 + 1];
                    dy[num8 + 2] += da * dx[num7 + 2];
                    dy[num8 + 3] += da * dx[num7 + 3];
                }
            } 
        }
         
    }

    public static void axpy(int n, double a, double[][] x, int xRow, int xoffset, double[][] y, int yRow, int yoffset) throws Exception {
        for (int num1 = 0;num1 < n;num1++)
        {
            y[yRow][num1 + yoffset] += a * x[xRow][num1 + xoffset];
        }
    }

    public static void axpy(int n, double da, double[] dx, int dx_off, int incx, double[][] dy, int dyrow, int dy_off, int incy) throws Exception {
        if ((n > 0) && (da != 0))
        {
            if ((incx != 1) || (incy != 1))
            {
                int num1 = dx_off;
                int num2 = dy_off;
                if (incx < 0)
                {
                    num1 = ((-n + 1) * incx) + dx_off;
                }
                 
                if (incy < 0)
                {
                    num2 = ((-n + 1) * incy) + dy_off;
                }
                 
                for (int num4 = 0;num4 < n;num4++)
                {
                    dy[dyrow][num2] += da * dx[num1];
                    num1 += incx;
                    num2 += incy;
                }
            }
            else
            {
                int num3 = n % 4;
                for (int num5 = 0;num5 < num3;num5++)
                {
                    dy[dyrow][num5 + dy_off] += da * dx[num5 + dx_off];
                }
                for (int num6 = num3;num6 < n;num6 += 4)
                {
                    int num7 = num6 + dx_off;
                    int num8 = num6 + dy_off;
                    dy[dyrow][num8] += da * dx[num7];
                    dy[dyrow][num8 + 1] += da * dx[num7 + 1];
                    dy[dyrow][num8 + 2] += da * dx[num7 + 2];
                    dy[dyrow][num8 + 3] += da * dx[num7 + 3];
                }
            } 
        }
         
    }

    public static void axpy(int n, double da, double[][] dx, int dxrow, int dx_off, int incx, double[] dy, int dy_off, int incy) throws Exception {
        if ((n > 0) && (da != 0))
        {
            if ((incx != 1) || (incy != 1))
            {
                int num1 = dx_off;
                int num2 = dy_off;
                if (incx < 0)
                {
                    num1 = ((-n + 1) * incx) + dx_off;
                }
                 
                if (incy < 0)
                {
                    num2 = ((-n + 1) * incy) + dy_off;
                }
                 
                for (int num4 = 0;num4 < n;num4++)
                {
                    dy[num2] += da * dx[dxrow][num1];
                    num1 += incx;
                    num2 += incy;
                }
            }
            else
            {
                int num3 = n % 4;
                for (int num5 = 0;num5 < num3;num5++)
                {
                    dy[num5 + dy_off] += da * dx[dxrow][num5 + dx_off];
                }
                for (int num6 = num3;num6 < n;num6 += 4)
                {
                    int num7 = num6 + dx_off;
                    int num8 = num6 + dy_off;
                    dy[num8] += da * dx[dxrow][num7];
                    dy[num8 + 1] += da * dx[dxrow][num7 + 1];
                    dy[num8 + 2] += da * dx[dxrow][num7 + 2];
                    dy[num8 + 3] += da * dx[dxrow][num7 + 3];
                }
            } 
        }
         
    }

    public static void axpy(int n, double da, double[][][] dx, int dxRow, int dxCol, int dx_off, int incx, double[][] dy, int dyRow, int dy_off, int incy) throws Exception {
        if ((n > 0) && (da != 0))
        {
            if ((incx != 1) || (incy != 1))
            {
                int num1 = dx_off;
                int num2 = dy_off;
                if (incx < 0)
                {
                    num1 = ((-n + 1) * incx) + dx_off;
                }
                 
                if (incy < 0)
                {
                    num2 = ((-n + 1) * incy) + dy_off;
                }
                 
                for (int num4 = 0;num4 < n;num4++)
                {
                    dy[dyRow][num2] += da * dx[dxRow][dxCol][num1];
                    num1 += incx;
                    num2 += incy;
                }
            }
            else
            {
                int num3 = n % 4;
                for (int num5 = 0;num5 < num3;num5++)
                {
                    dy[dyRow][num5 + dy_off] += da * dx[dxRow][dxCol][num5 + dx_off];
                }
                for (int num6 = num3;num6 < n;num6 += 4)
                {
                    int num7 = num6 + dx_off;
                    int num8 = num6 + dy_off;
                    dy[dyRow][num8] += da * dx[dxRow][dxCol][num7];
                    dy[dyRow][num8 + 1] += da * dx[dxRow][dxCol][num7 + 1];
                    dy[dyRow][num8 + 2] += da * dx[dxRow][dxCol][num7 + 2];
                    dy[dyRow][num8 + 3] += da * dx[dxRow][dxCol][num7 + 3];
                }
            } 
        }
         
    }

    public static void copy(int n, double[] x, double[] y) throws Exception {
        for (int num1 = 0;num1 < n;num1++)
        {
            y[num1] = x[num1];
        }
    }

    public static void copy(int n, int[] x, int[] y) throws Exception {
        for (int num1 = 0;num1 < n;num1++)
        {
            y[num1] = x[num1];
        }
    }

    public static void copy(int n, double[][] x, int xRow, double[] y) throws Exception {
        for (int num1 = 0;num1 < n;num1++)
        {
            y[num1] = x[xRow][num1];
        }
    }

    public static void copy(int n, double[] sx, int xoffset, int incx, double[] sy, int yoffset, int incy) throws Exception {
        if (n > 0)
        {
            int num2 = 0;
            int num3 = 0;
            if (incx < 0)
            {
                num2 = (-n + 1) * incx;
            }
             
            if (incy < 0)
            {
                num3 = (-n + 1) * incy;
            }
             
            for (int num1 = 0;num1 < n;num1++)
            {
                sy[yoffset + num3] = sx[xoffset + num2];
                num2 += incx;
                num3 += incy;
            }
        }
         
    }

    public static void copy(int n, int[] sx, int xoffset, int incx, int[] sy, int yoffset, int incy) throws Exception {
        if (n > 0)
        {
            int num2 = 0;
            int num3 = 0;
            if (incx < 0)
            {
                num2 = (-n + 1) * incx;
            }
             
            if (incy < 0)
            {
                num3 = (-n + 1) * incy;
            }
             
            for (int num1 = 0;num1 < n;num1++)
            {
                sy[yoffset + num3] = sx[xoffset + num2];
                num2 += incx;
                num3 += incy;
            }
        }
         
    }

    public static void copy(int n, double[][] sx, int xRow, int xoffset, int incx, double[] sy, int yoffset, int incy) throws Exception {
        if (n > 0)
        {
            int num2 = 0;
            int num3 = 0;
            if (incx < 0)
            {
                num2 = (-n + 1) * incx;
            }
             
            if (incy < 0)
            {
                num3 = (-n + 1) * incy;
            }
             
            for (int num1 = 0;num1 < n;num1++)
            {
                sy[yoffset + num3] = sx[xRow][xoffset + num2];
                num2 += incx;
                num3 += incy;
            }
        }
         
    }

    public static double dot(int n, double[] x, int xoffset, double[] y, int yoffset) throws Exception {
        double num1 = 0;
        for (int num2 = 0;num2 < n;num2++)
        {
            num1 += x[num2 + xoffset] * y[num2 + yoffset];
        }
        return num1;
    }

    public static double dot(int n, double[] x, int xoffset, double[][] y, int yRow, int yoffset) throws Exception {
        double num1 = 0;
        for (int num2 = 0;num2 < n;num2++)
        {
            num1 += x[num2 + xoffset] * y[yRow][num2 + yoffset];
        }
        return num1;
    }

    public static double dot(int n, double[][] x, int xRow, int xoffset, double[] y, int yoffset) throws Exception {
        double num1 = 0;
        for (int num2 = 0;num2 < n;num2++)
        {
            num1 += x[xRow][num2 + xoffset] * y[num2 + yoffset];
        }
        return num1;
    }

    public static double dot(int n, double[] dx, int dx_off, int incx, double[] dy, int dy_off, int incy) throws Exception {
        double num4 = 0;
        if (n > 0)
        {
            if ((incx != 1) || (incy != 1))
            {
                int num1 = dx_off;
                int num2 = dy_off;
                if (incx < 0)
                {
                    num1 = ((-n + 1) * incx) + dx_off;
                }
                 
                if (incy < 0)
                {
                    num2 = ((-n + 1) * incy) + dy_off;
                }
                 
                for (int num5 = 0;num5 < n;num5++)
                {
                    num4 += dx[num1] * dy[num2];
                    num1 += incx;
                    num2 += incy;
                }
                return num4;
            }
             
            int num3 = n % 5;
            for (int num6 = 0;num6 < num3;num6++)
            {
                num4 += dx[num6 + dx_off] * dy[num6 + dy_off];
            }
            for (int num7 = num3;num7 < n;num7 += 5)
            {
                int num8 = num7 + dx_off;
                int num9 = num7 + dy_off;
                num4 += ((((dx[num8] * dy[num9]) + (dx[num8 + 1] * dy[num9 + 1])) + (dx[num8 + 2] * dy[num9 + 2])) + (dx[num8 + 3] * dy[num9 + 3])) + (dx[num8 + 4] * dy[num9 + 4]);
            }
        }
         
        return num4;
    }

    public static double dot(int n, double[][] x, int xRow, int xoffset, double[][] y, int yRow, int yoffset) throws Exception {
        double num1 = 0;
        for (int num2 = 0;num2 < n;num2++)
        {
            num1 += x[xRow][num2 + xoffset] * y[yRow][num2 + yoffset];
        }
        return num1;
    }

    public static double dot(int n, double[][][] x, int xRow, int xCol, int xoffset, double[] y, int yoffset) throws Exception {
        double num1 = 0;
        for (int num2 = 0;num2 < n;num2++)
        {
            num1 += x[xRow][xCol][num2 + xoffset] * y[num2 + yoffset];
        }
        return num1;
    }

    public static double dot(int n, double[] dx, int dx_off, int incx, double[][] dy, int dyRow, int dy_off, int incy) throws Exception {
        double num4 = 0;
        if (n > 0)
        {
            if ((incx != 1) || (incy != 1))
            {
                int num1 = dx_off;
                int num2 = dy_off;
                if (incx < 0)
                {
                    num1 = ((-n + 1) * incx) + dx_off;
                }
                 
                if (incy < 0)
                {
                    num2 = ((-n + 1) * incy) + dy_off;
                }
                 
                for (int num5 = 0;num5 < n;num5++)
                {
                    num4 += dx[num1] * dy[dyRow][num2];
                    num1 += incx;
                    num2 += incy;
                }
                return num4;
            }
             
            int num3 = n % 5;
            for (int num6 = 0;num6 < num3;num6++)
            {
                num4 += dx[num6 + dx_off] * dy[dyRow][num6 + dy_off];
            }
            for (int num7 = num3;num7 < n;num7 += 5)
            {
                int num8 = num7 + dx_off;
                int num9 = num7 + dy_off;
                num4 += ((((dx[num8] * dy[dyRow][num9]) + (dx[num8 + 1] * dy[dyRow][num9 + 1])) + (dx[num8 + 2] * dy[dyRow][num9 + 2])) + (dx[num8 + 3] * dy[dyRow][num9 + 3])) + (dx[num8 + 4] * dy[dyRow][num9 + 4]);
            }
        }
         
        return num4;
    }

    public static double dot(int n, double[][] dx, int dxRow, int dx_off, int incx, double[] dy, int dy_off, int incy) throws Exception {
        double num4 = 0;
        if (n > 0)
        {
            if ((incx != 1) || (incy != 1))
            {
                int num1 = dx_off;
                int num2 = dy_off;
                if (incx < 0)
                {
                    num1 = ((-n + 1) * incx) + dx_off;
                }
                 
                if (incy < 0)
                {
                    num2 = ((-n + 1) * incy) + dy_off;
                }
                 
                for (int num5 = 0;num5 < n;num5++)
                {
                    num4 += dx[dxRow][num1] * dy[num2];
                    num1 += incx;
                    num2 += incy;
                }
                return num4;
            }
             
            int num3 = n % 5;
            for (int num6 = 0;num6 < num3;num6++)
            {
                num4 += dx[dxRow][num6 + dx_off] * dy[num6 + dy_off];
            }
            for (int num7 = num3;num7 < n;num7 += 5)
            {
                int num8 = num7 + dx_off;
                int num9 = num7 + dy_off;
                num4 += ((((dx[dxRow][num8] * dy[num9]) + (dx[dxRow][num8 + 1] * dy[num9 + 1])) + (dx[dxRow][num8 + 2] * dy[num9 + 2])) + (dx[dxRow][num8 + 3] * dy[num9 + 3])) + (dx[dxRow][num8 + 4] * dy[num9 + 4]);
            }
        }
         
        return num4;
    }

    public static double dot(int n, double[][] dx, int dxRow, int dx_off, int incx, double[][] dy, int dyRow, int dy_off, int incy) throws Exception {
        double num4 = 0;
        if (n > 0)
        {
            if ((incx != 1) || (incy != 1))
            {
                int num1 = dx_off;
                int num2 = dy_off;
                if (incx < 0)
                {
                    num1 = ((-n + 1) * incx) + dx_off;
                }
                 
                if (incy < 0)
                {
                    num2 = ((-n + 1) * incy) + dy_off;
                }
                 
                for (int num5 = 0;num5 < n;num5++)
                {
                    num4 += dx[dxRow][num1] * dy[dyRow][num2];
                    num1 += incx;
                    num2 += incy;
                }
                return num4;
            }
             
            int num3 = n % 5;
            for (int num6 = 0;num6 < num3;num6++)
            {
                num4 += dx[dxRow][num6 + dx_off] * dy[dyRow][num6 + dy_off];
            }
            for (int num7 = num3;num7 < n;num7 += 5)
            {
                int num8 = num7 + dx_off;
                int num9 = num7 + dy_off;
                num4 += ((((dx[dxRow][num8] * dy[dyRow][num9]) + (dx[dxRow][num8 + 1] * dy[dyRow][num9 + 1])) + (dx[dxRow][num8 + 2] * dy[dyRow][num9 + 2])) + (dx[dxRow][num8 + 3] * dy[dyRow][num9 + 3])) + (dx[dxRow][num8 + 4] * dy[dyRow][num9 + 4]);
            }
        }
         
        return num4;
    }

    public static double dot(int n, double[][][] x, int xRow, int xCol, int xoffset, double[][][] y, int yRow, int yCol, int yoffset) throws Exception {
        double num1 = 0;
        for (int num2 = 0;num2 < n;num2++)
        {
            num1 += x[xRow][xCol][num2 + xoffset] * y[yRow][yCol][num2 + yoffset];
        }
        return num1;
    }

    public static void dvcal(int n, double sa, double[] sx, int offset_sx, int incx, double[] sy, int offset_sy, int incy) throws Exception {
        if (n > 0)
        {
            int num1;
            if ((incx != 1) || (incy != 1))
            {
                int num2 = 1;
                int num3 = 1;
                for (num1 = 1;num1 <= n;num1++)
                {
                    sy[(offset_sy + num3) - 1] = sa * sx[(offset_sx + num2) - 1];
                    num2 += incx;
                    num3 += incy;
                }
            }
            else
            {
                for (num1 = 1;num1 <= n;num1++)
                {
                    sy[(offset_sy + num1) - 1] = sa * sx[(offset_sx + num1) - 1];
                }
            } 
        }
         
    }

    public static void gemm(char transa, char transb, int m, int n, int k, double alpha, double[] a, int offset_a, int lda, double[] b, int offset_b, int ldb, double beta, double[] c, int offset_c, int ldc) throws Exception {
        boolean flag1 = (transa == 'N') || (transa == 'n');
        boolean flag2 = (transb == 'N') || (transb == 'n');
        boolean flag3 = (((transa == 'T') || (transa == 't')) || (transa == 'C')) || (transa == 'c');
        boolean flag4 = (((transb == 'T') || (transb == 't')) || (transb == 'C')) || (transb == 'c');
        if (((m != 0) && (n != 0)) && (((alpha != 0) && (k != 0)) || (beta != 1)))
        {
            int num1;
            int num2;
            if (beta == 0)
            {
                for (num2 = 0;num2 < n;num2++)
                {
                    for (num1 = 0;num1 < m;num1++)
                    {
                        c[offset_c + num2 * ldc + num1] = 0;
                    }
                }
            }
            else if (beta == -1)
            {
                for (num2 = 0;num2 < n;num2++)
                {
                    for (num1 = 0;num1 < m;num1++)
                    {
                        c[offset_c + num2 * ldc + num1] = -c[offset_c + num2 * ldc + num1];
                    }
                }
            }
            else if (beta != 1)
            {
                for (num2 = 0;num2 < n;num2++)
                {
                    for (num1 = 0;num1 < m;num1++)
                    {
                        c[offset_c + num2 * ldc + num1] *= beta;
                    }
                }
            }
               
            if ((k != 0) && (alpha != 0))
            {
                int num3;
                double num4;
                if (flag2 && flag1)
                {
                    for (num3 = 0;num3 < k;num3++)
                    {
                        for (num2 = 0;num2 < n;num2++)
                        {
                            num4 = alpha * b[(offset_b + (num2 * ldb)) + num3];
                            for (num1 = 0;num1 < m;num1++)
                            {
                                c[offset_c + num2 * ldc + num1] += (num4 * a[offset_a + num3 * lda + num1]);
                            }
                        }
                    }
                }
                else if (flag2 && flag3)
                {
                    for (num3 = 0;num3 < k;num3++)
                    {
                        for (num2 = 0;num2 < n;num2++)
                        {
                            num4 = alpha * b[offset_b + num2 * ldb + num3];
                            for (num1 = 1;num1 <= m;num1++)
                            {
                                c[offset_c + num2 * ldc + num1] += num4 * a[offset_a + num1 * lda + num3];
                            }
                        }
                    }
                }
                else if (flag4 && flag3)
                {
                    for (num3 = 0;num3 < k;num3++)
                    {
                        for (num2 = 0;num2 < n;num2++)
                        {
                            num4 = alpha * b[offset_b + num3 * ldb + num2 - 1];
                            for (num1 = 0;num1 < m;num1++)
                            {
                                c[offset_c + num2 * ldc + num1] += num4 * a[offset_a + num1 * lda + num3];
                            }
                        }
                    }
                }
                else
                {
                    if (!flag4 || !flag1)
                    {
                        return ;
                    }
                     
                    for (num3 = 0;num3 < k;num3++)
                    {
                        for (num2 = 0;num2 < n;num2++)
                        {
                            num4 = alpha * b[offset_b + num3 * ldb + num2];
                            for (num1 = 0;num1 < m;num1++)
                            {
                                c[offset_c + num2 * ldc + num1] += num4 * a[offset_a + num3 * lda + num1];
                            }
                        }
                    }
                }   
            }
             
        }
         
    }

    public static void gemv(char trans, int m, int n, double alpha, double[][] a, int row_index, int col_index, double[] x, int x_off, double beta, double[] y, int y_off) throws Exception {
        boolean flag1 = false;
        if ((trans == 'N') || (trans == 'n'))
        {
            flag1 = true;
        }
         
        if (((m != 0) && (n != 0)) && ((alpha != 0) || (beta != 1))) {
            int num1;
            int num2;
            if (flag1) {
                num1 = n;
                num2 = m;
            } else {
                num1 = m;
                num2 = n;
            }
            double[] numArray1 = new double[Math.max(m, n)];
            if (beta != 1) {
                if (beta == 0) {
                    BLAS.set(num2, 0, y, y_off);
                } else {
                    BLAS.scal(num2, beta, y, y_off);
                }
            }

            if (alpha != 0) {
                if (flag1) {
                    int num3 = x_off;
                    for (int num4 = 0; num4 < n; num4++) {
                        for (int num5 = 0; num5 < m; num5++) {
                            numArray1[num5] = a[row_index + num4][col_index + num5];
                        }
                        BLAS.axpy(m, alpha * x[num3], numArray1, 0, y, y_off);
                        for (int num6 = 0; num6 < m; num6++) {
                            a[row_index + num4][col_index + num6] = numArray1[num6];
                        }
                        num3++;
                    }
                } else {
                    int num7 = y_off;
                    for (int num8 = 0; num8 < n; num8++) {
                        for (int num9 = 0; num9 < m; num9++) {
                            numArray1[num9] = a[row_index + num8][col_index + num9];
                        }
                        y[num7] += alpha * BLAS.dot(m, numArray1, 0, x, x_off);
                        for (int num10 = 0; num10 < m; num10++) {
                            a[row_index + num8][col_index + num10] = numArray1[num10];
                        }
                        num7++;
                    }
                }
            }

        }
    }

    public static void gemv(char trans, int m, int n, double alpha, double[][] a, int row_index, int col_index, double[][] x, int xRow, int x_off, double beta, double[] y, int y_off) throws Exception {
        boolean flag1 = false;
        if ((trans == 'N') || (trans == 'n'))
        {
            flag1 = true;
        }
         
        if (((m != 0) && (n != 0)) && ((alpha != 0) || (beta != 1)))
        {
            int num1;
            int num2;
            if (flag1)
            {
                num1 = n;
                num2 = m;
            }
            else
            {
                num1 = m;
                num2 = n;
            } 
            double[] numArray1 = new double[Math.max(m, n)];
            if (beta != 1)
            {
                if (beta == 0)
                {
                    BLAS.set(num2, 0, y, y_off);
                }
                else
                {
                    BLAS.scal(num2, beta, y, y_off);
                } 
            }
             
            if (alpha != 0)
            {
                if (flag1)
                {
                    int num3 = x_off;
                    for (int num4 = 0;num4 < n;num4++)
                    {
                        for (int num5 = 0;num5 < m;num5++)
                        {
                            numArray1[num5] = a[row_index + num4][col_index + num5];
                        }
                        BLAS.axpy(m, alpha * x[xRow][num3], numArray1, 0, y, y_off);
                        for (int num6 = 0;num6 < m;num6++)
                        {
                            a[row_index + num4][col_index + num6] = numArray1[num6];
                        }
                        num3++;
                    }
                }
                else
                {
                    int num7 = y_off;
                    for (int num8 = 0;num8 < n;num8++)
                    {
                        for (int num9 = 0;num9 < m;num9++)
                        {
                            numArray1[num9] = a[row_index + num8][col_index + num9];
                        }
                        y[num7] += alpha * BLAS.dot(m, numArray1, 0, x, xRow, x_off);
                        for (int num10 = 0;num10 < m;num10++)
                        {
                            a[row_index + num8][col_index + num10] = numArray1[num10];
                        }
                        num7++;
                    }
                } 
            }
             
        }
         
    }

    public static void gemv(char trans, int m, int n, double alpha, double[] a, int offset_a, int lda, double[] x, int offset_x, int incx, double beta, double[] y, int offset_y, int incy) throws Exception {
        boolean flag2 = (trans == 'N') || (trans == 'n');
        boolean flag3 = (trans == 'T') || (trans == 't');
        boolean flag1 = (trans == 'C') || (trans == 'c');
        if (((m != 0) && (n != 0)) && ((alpha != 0) || (beta != 1)))
        {
            int num6;
            int num7;
            if (flag2)
            {
                num6 = n;
                num7 = m;
            }
            else
            {
                num6 = m;
                num7 = n;
            } 
            int num2 = 1;
            int num3 = 1;
            if (incx < 0)
            {
                num2 = ((-num6 + 1) * incx) + 1;
            }
             
            if (incy < 0)
            {
                num3 = ((-num7 + 1) * incy) + 1;
            }
             
            if (beta != 1)
            {
                if (incy == 0)
                {
                    if (beta == 0)
                    {
                        y[offset_y] = 0;
                    }
                    else
                    {
                        y[offset_y] *= Math.pow(beta, (double) num7);
                    } 
                }
                else if (beta == 0)
                {
                    BLAS.set(num7, 0, y, offset_y, Math.abs(incy));
                }
                else
                {
                    BLAS.scal(num7, beta, y, offset_y, Math.abs(incy));
                }  
            }
             
            if (alpha != 0)
            {
                int num1;
                if (flag2)
                {
                    int num4 = num2;
                    for (num1 = 0;num1 < n;num1++)
                    {
                        BLAS.axpy(m, alpha * x[offset_x + num4], a, offset_a + lda * num1, 1, y, offset_y, incy);
                        num4 += incx;
                    }
                }
                else
                {
                    int num5 = num3;
                    for (num1 = 0;num1 < n;num1++)
                    {
                        y[offset_y + num5] += alpha * BLAS.dot(m, a, offset_a + lda * num1, 1, x, offset_x, incx);
                        num5 += incy;
                    }
                } 
            }
             
        }
         
    }

    public static void ger(int m, int n, double alpha, double[] x, int x_off, double[] y, int y_off, double[][] a, int row_index, int col_index) throws Exception {
        if (((m != 0) && (n != 0)) && (alpha != 0))
        {
            double[] numArray1 = new double[Math.max(m, n)];
            int num1 = y_off;
            for (int num2 = 0;num2 < n;num2++)
            {
                for (int num3 = 0;num3 < m;num3++)
                {
                    numArray1[num3] = a[row_index + num2][col_index + num3];
                }
                BLAS.axpy(m, alpha * y[num1], x, x_off, numArray1, 0);
                for (int num4 = 0;num4 < m;num4++)
                {
                    a[row_index + num2][col_index + num4] = numArray1[num4];
                }
                num1++;
            }
        }
         
    }

    public static void ger(int m, int n, double alpha, double[][] x, int xRow, int x_off, double[] y, int y_off, double[][] a, int row_index, int col_index) throws Exception {
        if (((m != 0) && (n != 0)) && (alpha != 0))
        {
            double[] numArray1 = new double[Math.max(m, n)];
            int num1 = y_off;
            for (int num2 = 0;num2 < n;num2++)
            {
                for (int num3 = 0;num3 < m;num3++)
                {
                    numArray1[num3] = a[row_index + num2][col_index + num3];
                }
                BLAS.axpy(m, alpha * y[num1], x, xRow, x_off, numArray1, 0);
                for (int num4 = 0;num4 < m;num4++)
                {
                    a[row_index + num2][col_index + num4] = numArray1[num4];
                }
                num1++;
            }
        }
         
    }

    public static void ger(int m, int n, double alpha, double[] x, int offset_x, int incx, double[] y, int offset_y, int incy, double[] a, int offset_a, int lda) throws Exception {
        if (((m != 0) && (n != 0)) && (alpha != 0))
        {
            int num2 = 1;
            if (incy < 0)
            {
                num2 = ((-n + 1) * incy) + 1;
            }
             
            int num1 = 1;
            for (int num3 = 1;num3 <= n;num3++)
            {
                BLAS.axpy(m, alpha * y[(offset_y + num2) - 1], x, offset_x, incx, a, (offset_a + num1) - 1, 1);
                num2 += incy;
                num1 += lda;
            }
        }
         
    }

    public static int iamax(int n, double[] sx, int offset_sx, int incx) throws Exception {
        int num3 = 0;
        if (n >= 1)
        {
            int num1;
            double num5;
            double num6;
            num3 = 1;
            if (n <= 1)
            {
                return num3;
            }
             
            if (incx != 1)
            {
                num5 = Math.abs(sx[offset_sx]);
                int num4 = n * incx;
                int num2 = 1;
                for (num1 = 1;num1 <= num4;num1 += incx)
                {
                    num6 = Math.abs(sx[(offset_sx + num1) - 1]);
                    if (num6 > num5)
                    {
                        num3 = num2;
                        num5 = num6;
                    }
                     
                    num2++;
                }
                return num3;
            }
             
            num5 = Math.abs(sx[offset_sx]);
            for (num1 = 2;num1 <= n;num1++)
            {
                num6 = Math.abs(sx[(offset_sx + num1) - 1]);
                if (num6 > num5)
                {
                    num3 = num1;
                    num5 = num6;
                }
                 
            }
        }
         
        return num3;
    }

    public static int iimax(int n, int[] ix, int offset_ix, int incx) throws Exception {
        int num2 = 0;
        if (n >= 1)
        {
            int num1;
            int num3;
            num2 = 1;
            if (n == 1)
            {
                return num2;
            }
             
            if (incx != 1)
            {
                int num4 = 1;
                num3 = ix[offset_ix];
                num4 += incx;
                for (num1 = 2;num1 <= n;num1++)
                {
                    if (ix[(offset_ix + num4) - 1] > num3)
                    {
                        num2 = num1;
                        num3 = ix[(offset_ix + num4) - 1];
                    }
                     
                    num4 += incx;
                }
                return num2;
            }
             
            num3 = ix[offset_ix];
            for (num1 = 2;num1 <= n;num1++)
            {
                if (ix[(offset_ix + num1) - 1] > num3)
                {
                    num2 = num1;
                    num3 = ix[(offset_ix + num1) - 1];
                }
                 
            }
        }
         
        return num2;
    }

    public static int iimin(int n, int[] ix, int offset_ix, int incx) throws Exception {
        int num2 = 0;
        if (n >= 1)
        {
            int num1;
            int num3;
            num2 = 1;
            if (n == 1)
            {
                return num2;
            }
             
            if (incx != 1)
            {
                int num4 = 1;
                num3 = ix[offset_ix];
                num4 += incx;
                for (num1 = 2;num1 <= n;num1++)
                {
                    if (ix[(offset_ix + num4) - 1] < num3)
                    {
                        num2 = num1;
                        num3 = ix[(offset_ix + num4) - 1];
                    }
                     
                    num4 += incx;
                }
                return num2;
            }
             
            num3 = ix[offset_ix];
            for (num1 = 2;num1 <= n;num1++)
            {
                if (ix[(offset_ix + num1) - 1] < num3)
                {
                    num2 = num1;
                    num3 = ix[(offset_ix + num1) - 1];
                }
                 
            }
        }
         
        return num2;
    }

    public static int ismax(int n, double[] sx, int offset_sx, int incx) throws Exception {
        int num2 = 0;
        if (n >= 1)
        {
            int num1;
            double num4;
            num2 = 1;
            if (n == 1)
            {
                return num2;
            }
             
            if (incx != 1)
            {
                int num3 = 1;
                num4 = sx[offset_sx];
                num3 += incx;
                for (num1 = 2;num1 <= n;num1++)
                {
                    if (sx[(offset_sx + num3) - 1] > num4)
                    {
                        num2 = num1;
                        num4 = sx[(offset_sx + num3) - 1];
                    }
                     
                    num3 += incx;
                }
                return num2;
            }
             
            num4 = sx[offset_sx];
            for (num1 = 2;num1 <= n;num1++)
            {
                if (sx[(offset_sx + num1) - 1] > num4)
                {
                    num2 = num1;
                    num4 = sx[(offset_sx + num1) - 1];
                }
                 
            }
        }
         
        return num2;
    }

    public static int ismax(int n, double[][] sx, int sxRow, int offset_sx, int incx) throws Exception {
        int num2 = 0;
        if (n >= 1)
        {
            int num1;
            double num4;
            num2 = 1;
            if (n == 1)
            {
                return num2;
            }
             
            if (incx != 1)
            {
                int num3 = 1;
                num4 = sx[sxRow][offset_sx];
                num3 += incx;
                for (num1 = 2;num1 <= n;num1++)
                {
                    if (sx[sxRow][(offset_sx + num3) - 1] > num4)
                    {
                        num2 = num1;
                        num4 = sx[sxRow][(offset_sx + num3) - 1];
                    }
                     
                    num3 += incx;
                }
                return num2;
            }
             
            num4 = sx[sxRow][offset_sx];
            for (num1 = 2;num1 <= n;num1++)
            {
                if (sx[sxRow][(offset_sx + num1) - 1] > num4)
                {
                    num2 = num1;
                    num4 = sx[sxRow][(offset_sx + num1) - 1];
                }
                 
            }
        }
         
        return num2;
    }

    public static int ismin(int n, double[] sx, int offset_sx, int incx) throws Exception {
        int num2 = 0;
        if (n >= 1)
        {
            int num1;
            double num4;
            num2 = 1;
            if (n == 1)
            {
                return num2;
            }
             
            if (incx != 1)
            {
                int num3 = 1;
                num4 = sx[offset_sx];
                num3 += incx;
                for (num1 = 2;num1 <= n;num1++)
                {
                    if (sx[(offset_sx + num3) - 1] < num4)
                    {
                        num2 = num1;
                        num4 = sx[(offset_sx + num3) - 1];
                    }
                     
                    num3 += incx;
                }
                return num2;
            }
             
            num4 = sx[offset_sx];
            for (num1 = 2;num1 <= n;num1++)
            {
                if (sx[(offset_sx + num1) - 1] < num4)
                {
                    num2 = num1;
                    num4 = sx[(offset_sx + num1) - 1];
                }
                 
            }
        }
         
        return num2;
    }

    public static int ismin(int n, double[][] sx, int sxRow, int offset_sx, int incx) throws Exception {
        int num2 = 0;
        if (n >= 1)
        {
            int num1;
            double num4;
            num2 = 1;
            if (n == 1)
            {
                return num2;
            }
             
            if (incx != 1)
            {
                int num3 = 1;
                num4 = sx[sxRow][offset_sx];
                num3 += incx;
                for (num1 = 2;num1 <= n;num1++)
                {
                    if (sx[sxRow][(offset_sx + num3) - 1] < num4)
                    {
                        num2 = num1;
                        num4 = sx[sxRow][(offset_sx + num3) - 1];
                    }
                     
                    num3 += incx;
                }
                return num2;
            }
             
            num4 = sx[sxRow][offset_sx];
            for (num1 = 2;num1 <= n;num1++)
            {
                if (sx[sxRow][(offset_sx + num1) - 1] < num4)
                {
                    num2 = num1;
                    num4 = sx[sxRow][(offset_sx + num1) - 1];
                }
                 
            }
        }
         
        return num2;
    }

    public static double nrm2(int n, double[] dx, int offset) throws Exception {
        double num2;
        double num1 = 0;
        double num3 = 0;
        for (int num4 = 0;num4 < n;num4++)
        {
            num1 += Math.abs(dx[num4 + offset]);
            if (num1 > BLAS._threshold[3])
            {
                num1 = BLAS._threshold[4];
                num2 = BLAS._threshold[5];
                num3 = 0;
                for (int num5 = 0;num5 < n;num5++)
                {
                    num3 += (num1 * dx[num5 + offset]) * (num1 * dx[num5 + offset]);
                }
                return (Math.sqrt(num3) * num2);
            }
             
        }
        if (num1 < (n * BLAS._threshold[0]))
        {
            num1 = BLAS._threshold[1];
            num2 = BLAS._threshold[2];
            num3 = 0;
            for (int num6 = 0;num6 < n;num6++)
            {
                num3 += (num1 * dx[num6 + offset]) * (num1 * dx[num6 + offset]);
            }
            return (Math.sqrt(num3) * num2);
        }
         
        for (int num7 = 0;num7 < n;num7++)
        {
            num3 += dx[num7 + offset] * dx[num7 + offset];
        }
        return Math.sqrt(num3);
    }

    public static double nrm2(int n, double[][] dx, int dxRow, int offset) throws Exception {
        double num2;
        double num1 = 0;
        double num3 = 0;
        for (int num4 = 0;num4 < n;num4++)
        {
            num1 += Math.abs(dx[dxRow][num4 + offset]);
            if (num1 > BLAS._threshold[3])
            {
                num1 = BLAS._threshold[4];
                num2 = BLAS._threshold[5];
                num3 = 0;
                for (int num5 = 0;num5 < n;num5++)
                {
                    num3 += (num1 * dx[dxRow][num5 + offset]) * (num1 * dx[dxRow][num5 + offset]);
                }
                return (Math.sqrt(num3) * num2);
            }
             
        }
        if (num1 < (n * BLAS._threshold[0]))
        {
            num1 = BLAS._threshold[1];
            num2 = BLAS._threshold[2];
            num3 = 0;
            for (int num6 = 0;num6 < n;num6++)
            {
                num3 += (num1 * dx[dxRow][num6 + offset]) * (num1 * dx[dxRow][num6 + offset]);
            }
            return (Math.sqrt(num3) * num2);
        }
         
        for (int num7 = 0;num7 < n;num7++)
        {
            num3 += dx[dxRow][num7 + offset] * dx[dxRow][num7 + offset];
        }
        return Math.sqrt(num3);
    }

    public static double nrm2(int n, double[] dx, int dx_off, int incx) throws Exception {
        int num2;
        double num4;
        double[] numArray1 = new double[6];
        double num6 = 0;
        numArray1[0] = 1.0010415475916E-146;
        numArray1[1] = 4.4989137945432E+161;
        numArray1[2] = 2.2227587494851E-162;
        numArray1[3] = 1.9979190722022E+146;
        numArray1[4] = 5.0104209000224E-293;
        numArray1[5] = 1.9958403095347E+292;
        double num5 = num6;
        if ((n <= 0) || (incx < 0))
        {
            return num5;
        }
         
        if (incx == 0)
        {
            num5 = n;
            return (Math.abs(dx[0]) * Math.sqrt(num5));
        }
         
        double num3 = num6;
        int num1 = 1;
        while (num1 <= (n * incx))
        {
            num3 += Math.abs(dx[(num1 - 1) + dx_off]);
            if (num3 > numArray1[3])
            {
                num3 = numArray1[4];
                num4 = numArray1[5];
                num5 = num6;
                for (num2 = 1;num2 <= (n * incx);num2++)
                {
                    num5 += (num3 * dx[(num2 - 1) + dx_off]) * (num3 * dx[(num2 - 1) + dx_off]);
                }
                return (Math.sqrt(num5) * num4);
            }
             
            num1++;
        }
        if (num3 < (n * numArray1[0]))
        {
            num3 = numArray1[1];
            num4 = numArray1[2];
            num5 = num6;
            for (num2 = 1;num2 <= (n * incx);num2++)
            {
                num5 += (num3 * dx[(num2 - 1) + dx_off]) * (num3 * dx[(num2 - 1) + dx_off]);
            }
            return (Math.sqrt(num5) * num4);
        }
         
        for (num1 = 1;num1 <= (n * incx);num1++)
        {
            num5 += dx[(num1 - 1) + dx_off] * dx[(num1 - 1) + dx_off];
        }
        return Math.sqrt(num5);
    }

    public static void rot(int n, double[] dx, int incx, double[] dy, int incy, double[] c, double[] s) throws Exception {
        if (n > 0)
        {
            for (int num1 = 0;num1 < n;num1++)
            {
                double num2 = (c[0] * dx[num1 + incx]) + (s[0] * dy[num1 + incy]);
                dy[num1 + incy] = (c[0] * dy[num1 + incy]) - (s[0] * dx[num1 + incx]);
                dx[num1 + incx] = num2;
            }
        }
         
    }

    public static void rot(int n, double[][] dx, int dxRow, int incx, double[][] dy, int dyRow, int incy, double[] c, double[] s) throws Exception {
        if (n > 0)
        {
            for (int num1 = 0;num1 < n;num1++)
            {
                double num2 = (c[0] * dx[dxRow][num1 + incx]) + (s[0] * dy[dyRow][num1 + incy]);
                dy[dyRow][num1 + incy] = (c[0] * dy[dyRow][num1 + incy]) - (s[0] * dx[dxRow][num1 + incx]);
                dx[dxRow][num1 + incx] = num2;
            }
        }
         
    }

    public static void rot(int n, double[] dx, int dx_off, int incx, double[] dy, int dy_off, int incy, double[] c, double[] s) throws Exception {
        if (n > 0)
        {
            int num1;
            double num4;
            if ((incx != 1) || (incy != 1))
            {
                int num2 = 1 + dx_off;
                int num3 = 1 + dy_off;
                if (incx < 0)
                {
                    num2 = (((-n + 1) * incx) + 1) + dx_off;
                }
                 
                if (incy < 0)
                {
                    num3 = (((-n + 1) * incy) + 1) + dy_off;
                }
                 
                for (num1 = 0;num1 < n;num1++)
                {
                    num4 = (c[0] * dx[num2]) + (s[0] * dy[num3]);
                    dy[num3] = (c[0] * dy[num3]) - (s[0] * dx[num2]);
                    dx[num2] = num4;
                    num2 += incx;
                    num3 += incy;
                }
            }
            else
            {
                for (num1 = 0;num1 < n;num1++)
                {
                    num4 = (c[0] * dx[num1 + dx_off]) + (s[0] * dy[num1 + dy_off]);
                    dy[num1 + dy_off] = (c[0] * dy[num1 + dy_off]) - (s[0] * dx[num1 + dx_off]);
                    dx[num1 + dx_off] = num4;
                }
            } 
        }
         
    }

    public static void rot(int n, double[][] dx, int dxRow, int dx_off, int incx, double[][] dy, int dyRow, int dy_off, int incy, double[] c, double[] s) throws Exception {
        if (n > 0)
        {
            int num1;
            double num4;
            if ((incx != 1) || (incy != 1))
            {
                int num2 = 1 + dx_off;
                int num3 = 1 + dy_off;
                if (incx < 0)
                {
                    num2 = (((-n + 1) * incx) + 1) + dx_off;
                }
                 
                if (incy < 0)
                {
                    num3 = (((-n + 1) * incy) + 1) + dy_off;
                }
                 
                for (num1 = 0;num1 < n;num1++)
                {
                    num4 = (c[0] * dx[dxRow][num2]) + (s[0] * dy[dyRow][num3]);
                    dy[dyRow][num3] = (c[0] * dy[dyRow][num3]) - (s[0] * dx[dxRow][num2]);
                    dx[dxRow][num2] = num4;
                    num2 += incx;
                    num3 += incy;
                }
            }
            else
            {
                for (num1 = 0;num1 < n;num1++)
                {
                    num4 = (c[0] * dx[dxRow][num1 + dx_off]) + (s[0] * dy[dyRow][num1 + dy_off]);
                    dy[dyRow][num1 + dy_off] = (c[0] * dy[dyRow][num1 + dy_off]) - (s[0] * dx[dxRow][num1 + dx_off]);
                    dx[dxRow][num1 + dx_off] = num4;
                }
            } 
        }
         
    }

    public static void rotg(double[] da, double[] db, double[] dc, double[] ds) throws Exception {
        double num2;
        double num3;
        if (Math.abs(da[0]) > Math.abs(db[0]))
        {
            num2 = da[0] + da[0];
            num3 = db[0] / num2;
            double num1 = Math.sqrt(0.25 + (num3 * num3)) * num2;
            dc[0] = da[0] / num1;
            ds[0] = num3 * (dc[0] + dc[0]);
            db[0] = ds[0];
            da[0] = num1;
        }
        else if (db[0] == 0)
        {
            dc[0] = 1;
            ds[0] = 0;
            da[0] = 0;
            db[0] = 0;
        }
        else
        {
            num2 = db[0] + db[0];
            num3 = da[0] / num2;
            da[0] = Math.sqrt(0.25 + (num3 * num3)) * num2;
            ds[0] = db[0] / da[0];
            dc[0] = num3 * (ds[0] + ds[0]);
            if (dc[0] == 0)
            {
                db[0] = 1;
            }
            else
            {
                db[0] = 1 / dc[0];
            } 
        }  
    }

    public static void rotm(int n, double[] sx, int incx, double[] sy, int incy, double[] sparam) throws Exception {
        BLAS.rotm(n, sx, 0, incx, sy, 0, incy, sparam);
    }

    public static void rotm(int n, double[] sx, int xOffset, int incx, double[] sy, int yOffset, int incy, double[] sparam) throws Exception {
        double num4 = sparam[0];
        if ((n > 0) && (num4 != -2))
        {
            int num1;
            double num5;
            double num6;
            double num7;
            double num8;
            double num9;
            double num10;
            if ((incx == 1) && (incy == 1))
            {
                if (num4 == 0)
                {
                    num6 = sparam[3];
                    num7 = sparam[2];
                    for (num1 = 1;num1 <= n;num1++)
                    {
                        num9 = sx[(xOffset + num1) - 1];
                        num10 = sy[(yOffset + num1) - 1];
                        sx[(xOffset + num1) - 1] = num9 + (num10 * num6);
                        sy[(yOffset + num1) - 1] = (num9 * num7) + num10;
                    }
                }
                else if (num4 > 0)
                {
                    num5 = sparam[1];
                    num8 = sparam[4];
                    for (num1 = 1;num1 <= n;num1++)
                    {
                        num9 = sx[(xOffset + num1) - 1];
                        num10 = sy[(yOffset + num1) - 1];
                        sx[(xOffset + num1) - 1] = (num9 * num5) + num10;
                        sy[(yOffset + num1) - 1] = -num9 + (num8 * num10);
                    }
                }
                else
                {
                    if (num4 >= 0)
                    {
                        return ;
                    }
                     
                    num5 = sparam[1];
                    num6 = sparam[3];
                    num7 = sparam[2];
                    num8 = sparam[4];
                    for (num1 = 1;num1 <= n;num1++)
                    {
                        num9 = sx[(xOffset + num1) - 1];
                        num10 = sy[(yOffset + num1) - 1];
                        sx[(xOffset + num1) - 1] = (num9 * num5) + (num10 * num6);
                        sy[(yOffset + num1) - 1] = (num9 * num7) + (num10 * num8);
                    }
                }  
            }
            else
            {
                int num2 = 1;
                int num3 = 1;
                if (incx < 0)
                {
                    num2 = 1 + ((1 - n) * incx);
                }
                 
                if (incy < 0)
                {
                    num3 = 1 + ((1 - n) * incy);
                }
                 
                if (num4 == 0)
                {
                    num6 = sparam[3];
                    num7 = sparam[2];
                    for (num1 = 1;num1 <= n;num1++)
                    {
                        num9 = sx[(xOffset + num2) - 1];
                        num10 = sy[(yOffset + num3) - 1];
                        sx[(xOffset + num2) - 1] = num9 + (num10 * num6);
                        sy[(yOffset + num3) - 1] = (num9 * num7) + num10;
                        num2 += incx;
                        num3 += incy;
                    }
                }
                else if (num4 > 0)
                {
                    num5 = sparam[1];
                    num8 = sparam[4];
                    for (num1 = 1;num1 <= n;num1++)
                    {
                        num9 = sx[(xOffset + num2) - 1];
                        num10 = sy[(yOffset + num3) - 1];
                        sx[(xOffset + num2) - 1] = (num9 * num5) + num10;
                        sy[(yOffset + num3) - 1] = -num9 + (num8 * num10);
                        num2 += incx;
                        num3 += incy;
                    }
                }
                else
                {
                    if (num4 >= 0)
                    {
                        return ;
                    }
                     
                    num5 = sparam[1];
                    num6 = sparam[3];
                    num7 = sparam[2];
                    num8 = sparam[4];
                    for (num1 = 1;num1 <= n;num1++)
                    {
                        num9 = sx[(xOffset + num2) - 1];
                        num10 = sy[(yOffset + num3) - 1];
                        sx[(xOffset + num2) - 1] = (num9 * num5) + (num10 * num6);
                        sy[(yOffset + num3) - 1] = (num9 * num7) + (num10 * num8);
                        num2 += incx;
                        num3 += incy;
                    }
                }  
            } 
        }
         
    }

    public static void rotm(int n, double[][] sx, int sxrow, int xOffset, int incx, double[] sy, int yOffset, int incy, double[] sparam) throws Exception {
        double num4 = sparam[0];
        if ((n > 0) && (num4 != -2))
        {
            int num1;
            double num5;
            double num6;
            double num7;
            double num8;
            double num9;
            double num10;
            if ((incx == 1) && (incy == 1))
            {
                if (num4 == 0)
                {
                    num6 = sparam[3];
                    num7 = sparam[2];
                    for (num1 = 1;num1 <= n;num1++)
                    {
                        num9 = sx[sxrow][(xOffset + num1) - 1];
                        num10 = sy[(yOffset + num1) - 1];
                        sx[sxrow][(xOffset + num1) - 1] = num9 + (num10 * num6);
                        sy[(yOffset + num1) - 1] = (num9 * num7) + num10;
                    }
                }
                else if (num4 > 0)
                {
                    num5 = sparam[1];
                    num8 = sparam[4];
                    for (num1 = 1;num1 <= n;num1++)
                    {
                        num9 = sx[sxrow][(xOffset + num1) - 1];
                        num10 = sy[(yOffset + num1) - 1];
                        sx[sxrow][(xOffset + num1) - 1] = (num9 * num5) + num10;
                        sy[(yOffset + num1) - 1] = -num9 + (num8 * num10);
                    }
                }
                else
                {
                    if (num4 >= 0)
                    {
                        return ;
                    }
                     
                    num5 = sparam[1];
                    num6 = sparam[3];
                    num7 = sparam[2];
                    num8 = sparam[4];
                    for (num1 = 1;num1 <= n;num1++)
                    {
                        num9 = sx[sxrow][(xOffset + num1) - 1];
                        num10 = sy[(yOffset + num1) - 1];
                        sx[sxrow][(xOffset + num1) - 1] = (num9 * num5) + (num10 * num6);
                        sy[(yOffset + num1) - 1] = (num9 * num7) + (num10 * num8);
                    }
                }  
            }
            else
            {
                int num2 = 1;
                int num3 = 1;
                if (incx < 0)
                {
                    num2 = 1 + ((1 - n) * incx);
                }
                 
                if (incy < 0)
                {
                    num3 = 1 + ((1 - n) * incy);
                }
                 
                if (num4 == 0)
                {
                    num6 = sparam[3];
                    num7 = sparam[2];
                    for (num1 = 1;num1 <= n;num1++)
                    {
                        num9 = sx[sxrow][(xOffset + num2) - 1];
                        num10 = sy[(yOffset + num3) - 1];
                        sx[sxrow][(xOffset + num2) - 1] = num9 + (num10 * num6);
                        sy[(yOffset + num3) - 1] = (num9 * num7) + num10;
                        num2 += incx;
                        num3 += incy;
                    }
                }
                else if (num4 > 0)
                {
                    num5 = sparam[1];
                    num8 = sparam[4];
                    for (num1 = 1;num1 <= n;num1++)
                    {
                        num9 = sx[sxrow][(xOffset + num2) - 1];
                        num10 = sy[(yOffset + num3) - 1];
                        sx[sxrow][(xOffset + num2) - 1] = (num9 * num5) + num10;
                        sy[(yOffset + num3) - 1] = -num9 + (num8 * num10);
                        num2 += incx;
                        num3 += incy;
                    }
                }
                else
                {
                    if (num4 >= 0)
                    {
                        return ;
                    }
                     
                    num5 = sparam[1];
                    num6 = sparam[3];
                    num7 = sparam[2];
                    num8 = sparam[4];
                    for (num1 = 1;num1 <= n;num1++)
                    {
                        num9 = sx[sxrow][(xOffset + num2) - 1];
                        num10 = sy[(yOffset + num3) - 1];
                        sx[sxrow][(xOffset + num2) - 1] = (num9 * num5) + (num10 * num6);
                        sy[(yOffset + num3) - 1] = (num9 * num7) + (num10 * num8);
                        num2 += incx;
                        num3 += incy;
                    }
                }  
            } 
        }
         
    }

    public static void rotm(int n, double[][] sx, int sxrow, int xOffset, int incx, double[][] sy, int syrow, int yOffset, int incy, double[] sparam) throws Exception {
        double num4 = sparam[0];
        if ((n > 0) && (num4 != -2))
        {
            int num1;
            double num5;
            double num6;
            double num7;
            double num8;
            double num9;
            double num10;
            if ((incx == 1) && (incy == 1))
            {
                if (num4 == 0)
                {
                    num6 = sparam[3];
                    num7 = sparam[2];
                    for (num1 = 1;num1 <= n;num1++)
                    {
                        num9 = sx[sxrow][(xOffset + num1) - 1];
                        num10 = sy[syrow][(yOffset + num1) - 1];
                        sx[sxrow][(xOffset + num1) - 1] = num9 + (num10 * num6);
                        sy[syrow][(yOffset + num1) - 1] = (num9 * num7) + num10;
                    }
                }
                else if (num4 > 0)
                {
                    num5 = sparam[1];
                    num8 = sparam[4];
                    for (num1 = 1;num1 <= n;num1++)
                    {
                        num9 = sx[sxrow][(xOffset + num1) - 1];
                        num10 = sy[syrow][(yOffset + num1) - 1];
                        sx[sxrow][(xOffset + num1) - 1] = (num9 * num5) + num10;
                        sy[syrow][(yOffset + num1) - 1] = -num9 + (num8 * num10);
                    }
                }
                else
                {
                    if (num4 >= 0)
                    {
                        return ;
                    }
                     
                    num5 = sparam[1];
                    num6 = sparam[3];
                    num7 = sparam[2];
                    num8 = sparam[4];
                    for (num1 = 1;num1 <= n;num1++)
                    {
                        num9 = sx[sxrow][(xOffset + num1) - 1];
                        num10 = sy[syrow][(yOffset + num1) - 1];
                        sx[sxrow][(xOffset + num1) - 1] = (num9 * num5) + (num10 * num6);
                        sy[syrow][(yOffset + num1) - 1] = (num9 * num7) + (num10 * num8);
                    }
                }  
            }
            else
            {
                int num2 = 1;
                int num3 = 1;
                if (incx < 0)
                {
                    num2 = 1 + ((1 - n) * incx);
                }
                 
                if (incy < 0)
                {
                    num3 = 1 + ((1 - n) * incy);
                }
                 
                if (num4 == 0)
                {
                    num6 = sparam[3];
                    num7 = sparam[2];
                    for (num1 = 1;num1 <= n;num1++)
                    {
                        num9 = sx[sxrow][(xOffset + num2) - 1];
                        num10 = sy[syrow][(yOffset + num3) - 1];
                        sx[sxrow][(xOffset + num2) - 1] = num9 + (num10 * num6);
                        sy[syrow][(yOffset + num3) - 1] = (num9 * num7) + num10;
                        num2 += incx;
                        num3 += incy;
                    }
                }
                else if (num4 > 0)
                {
                    num5 = sparam[1];
                    num8 = sparam[4];
                    for (num1 = 1;num1 <= n;num1++)
                    {
                        num9 = sx[sxrow][(xOffset + num2) - 1];
                        num10 = sy[syrow][(yOffset + num3) - 1];
                        sx[sxrow][(xOffset + num2) - 1] = (num9 * num5) + num10;
                        sy[syrow][(yOffset + num3) - 1] = -num9 + (num8 * num10);
                        num2 += incx;
                        num3 += incy;
                    }
                }
                else
                {
                    if (num4 >= 0)
                    {
                        return ;
                    }
                     
                    num5 = sparam[1];
                    num6 = sparam[3];
                    num7 = sparam[2];
                    num8 = sparam[4];
                    for (num1 = 1;num1 <= n;num1++)
                    {
                        num9 = sx[sxrow][(xOffset + num2) - 1];
                        num10 = sy[syrow][(yOffset + num3) - 1];
                        sx[sxrow][(xOffset + num2) - 1] = (num9 * num5) + (num10 * num6);
                        sy[syrow][(yOffset + num3) - 1] = (num9 * num7) + (num10 * num8);
                        num2 += incx;
                        num3 += incy;
                    }
                }  
            } 
        }
         
    }

    public static void rotmg(double[] givens, double[] sparam) throws Exception {
        double num1 = 0;
        double num2 = 0;
        double num3 = 0;
        double num4 = 0;
        double num5 = 0;
        double num6 = 0;
        double num7 = 0;
        double num8 = 0;
        double num9 = 0;
        double num10 = 0;
        double num11 = 0;
        double num12 = 0;
        double num13 = 1;
        double num14 = 2;
        double num15 = 4096;
        double num16 = 16780000;
        double num17 = 5.96E-08;
        boolean flag1 = true;
        if (givens[0] >= num12)
        {
            double[] numArray1 = new double[3];
            num7 = givens[1] * givens[3];
            if (num7 == num12)
            {
                num1 = -num14;
                sparam[0] = num1;
                return ;
            }
             
            num6 = givens[0] * givens[2];
            num9 = num7 * givens[3];
            num8 = num6 * givens[2];
            if (Math.abs(num8) > Math.abs(num9))
            {
                num4 = -givens[3] / givens[2];
                num3 = num7 / num6;
                num11 = num13 - (num3 * num4);
                if (num11 > num12)
                {
                    num1 = num12;
                    (numArray1 = givens)[0] = numArray1[0] / num11;
                    (numArray1 = givens)[1] = numArray1[1] / num11;
                    (numArray1 = givens)[2] = numArray1[2] * num11;
                    flag1 = true;
                }
                else
                {
                    flag1 = false;
                } 
            }
            else if (num9 >= num12)
            {
                num1 = num13;
                num2 = num6 / num7;
                num5 = givens[2] / givens[3];
                num11 = num13 + (num2 * num5);
                num10 = givens[1] / num11;
                givens[1] = givens[0] / num11;
                givens[0] = num10;
                givens[2] = givens[3] * num11;
                flag1 = true;
            }
            else
            {
                flag1 = false;
            }  
            if (flag1)
            {
                while (givens[0] <= num17)
                {
                    if (num1 == num12)
                    {
                        num2 = num13;
                        num5 = num13;
                        num1 = -num13;
                    }
                    else
                    {
                        num4 = -num13;
                        num3 = num13;
                        num1 = -num13;
                    } 
                    (numArray1 = givens)[0] = numArray1[0] * Math.pow(num15, 2);
                    (numArray1 = givens)[2] = numArray1[2] / num15;
                    num2 /= num15;
                    num3 /= num15;
                }
                if (givens[0] != num12)
                {
                    while (givens[0] >= num16)
                    {
                        if (num1 == num12)
                        {
                            num2 = num13;
                            num5 = num13;
                            num1 = -num13;
                        }
                        else
                        {
                            num4 = -num13;
                            num3 = num13;
                            num1 = -num13;
                        } 
                        givens[0] /= Math.pow(num15, 2);
                        (numArray1 = givens)[2] = numArray1[2] * num15;
                        num2 *= num15;
                        num3 *= num15;
                    }
                }
                 
                if (givens[1] != num12)
                {
                    while (Math.abs(givens[1]) <= num17)
                    {
                        if (num1 == num12)
                        {
                            num2 = num13;
                            num5 = num13;
                            num1 = -num13;
                        }
                        else
                        {
                            num4 = -num13;
                            num3 = num13;
                            num1 = -num13;
                        } 
                        (numArray1 = givens)[1] = numArray1[1] * Math.pow(num15, 2);
                        num4 /= num15;
                        num5 /= num15;
                    }
                    while (Math.abs(givens[1]) >= num16)
                    {
                        if (num1 == num12)
                        {
                            num2 = num13;
                            num5 = num13;
                            num1 = -num13;
                        }
                        else
                        {
                            num4 = -num13;
                            num3 = num13;
                            num1 = -num13;
                        } 
                        givens[1] /= Math.pow(num15, 2);
                        num4 *= num15;
                        num5 *= num15;
                    }
                }
                 
            }
            else
            {
                num1 = -num13;
                num2 = num12;
                num3 = num12;
                num4 = num12;
                num5 = num12;
                givens[0] = num12;
                givens[1] = num12;
                givens[2] = num12;
            } 
        }
        else
        {
            num1 = -num13;
            num2 = num12;
            num3 = num12;
            num4 = num12;
            num5 = num12;
            givens[0] = num12;
            givens[1] = num12;
            givens[2] = num12;
        } 
        if (num1 < num12)
        {
            sparam[1] = num2;
            sparam[2] = num4;
            sparam[3] = num3;
            sparam[4] = num5;
        }
        else if (num1 == num12)
        {
            sparam[2] = num4;
            sparam[3] = num3;
        }
        else
        {
            sparam[1] = num2;
            sparam[4] = num5;
        }  
        sparam[0] = num1;
    }

    public static void scal(int n, double da, double[] dx, int ix) throws Exception {
        for (int num1 = 0;num1 < n;num1++)
        {
            dx[num1 + ix] *= da;
        }
    }

    public static void scal(int n, double da, double[] dx, int dx_off, int incx) throws Exception {
        if (n > 0)
        {
            int num1;
            if (incx != 1)
            {
                int num4 = n * incx;
                for (num1 = 0;num1 < num4;num1 += incx)
                {
                    dx[num1 + dx_off] *= da;
                }
            }
            else
            {
                int num3 = n - ((n / 5) * 5);
                num1 = 0;
                while (num1 < num3)
                {
                    dx[num1 + dx_off] *= da;
                    num1++;
                }
                for (num1 = num3;num1 < n;num1 += 5)
                {
                    int num2 = num1 + dx_off;
                    dx[num2] *= da;
                    dx[num2 + 1] *= da;
                    dx[num2 + 2] *= da;
                    dx[num2 + 3] *= da;
                    dx[num2 + 4] *= da;
                }
            } 
        }
         
    }

    public static void scal(int n, double da, double[][] dx, int dxRow, int ix) throws Exception {
        for (int num1 = 0;num1 < n;num1++)
        {
            dx[dxRow][num1 + ix] *= da;
        }
    }

    public static void set(int n, double da, double[] dx, int ix) throws Exception {
        for (int num1 = 0;num1 < n;num1++)
        {
            dx[num1 + ix] = da;
        }
    }

    public static void set(int n, int da, int[] dx, int ix) throws Exception {
        for (int num1 = 0;num1 < n;num1++)
        {
            dx[num1 + ix] = da;
        }
    }

    public static void set(int n, double da, double[] dx, int ix, int incx) throws Exception {
        if (n > 0)
        {
            if (incx != 1)
            {
                if (incx < 0)
                {
                    ix += (-n + 1) * incx;
                }
                 
                for (int num1 = 0;num1 < n;num1++)
                {
                    dx[ix] = da;
                    ix += incx;
                }
            }
            else
            {
                for (int num2 = 0;num2 < n;num2++)
                {
                    dx[num2 + ix] = da;
                }
            } 
        }
         
    }

    public static void set(int n, double da, double[][] dx, int dxRow, int ix) throws Exception {
        for (int num1 = 0;num1 < n;num1++)
        {
            dx[dxRow][num1 + ix] = da;
        }
    }

    public static void set(int n, double da, double[][][] dx, int dxRow, int dxCol, int ix) throws Exception {
        for (int num1 = 0;num1 < n;num1++)
        {
            dx[dxRow][dxCol][num1 + ix] = da;
        }
    }

    public static double sum(int n, double[] sx, int offset_sx, int incx) throws Exception {
        double num3 = 0;
        if (n > 0)
        {
            int num1;
            if (incx != 1)
            {
                int num2 = n * incx;
                for (num1 = 1;num1 <= num2;num1 += incx)
                {
                    num3 += sx[(offset_sx + num1) - 1];
                }
                return num3;
            }
             
            for (num1 = 1;num1 <= n;num1++)
            {
                num3 += sx[(offset_sx + num1) - 1];
            }
        }
         
        return num3;
    }

    public static int sum(int n, int[] ix, int offset_ix, int incx) throws Exception {
        int num4 = 0;
        if (n > 0)
        {
            int num3;
            if (incx != 1)
            {
                int num5 = n * incx;
                num3 = 1;
                int num2 = incx;
                for (int num1 = ((((num5 - 1) + num2) / num2) > 0) ? (((num5 - 1) + num2) / num2) : 0;num1 > 0;num1--)
                {
                    num4 += ix[(offset_ix + num3) - 1];
                    num3 += num2;
                }
                return num4;
            }
             
            for (num3 = 1;num3 <= n;num3++)
            {
                num4 += ix[(offset_ix + num3) - 1];
            }
        }
         
        return num4;
    }

    public static double sum(int n, double[][] sx, int sxRow, int offset_sx, int incx) throws Exception {
        double num3 = 0;
        if (n > 0)
        {
            int num1;
            if (incx != 1)
            {
                int num2 = n * incx;
                for (num1 = 1;num1 <= num2;num1 += incx)
                {
                    num3 += sx[sxRow][(offset_sx + num1) - 1];
                }
                return num3;
            }
             
            for (num1 = 1;num1 <= n;num1++)
            {
                num3 += sx[sxRow][(offset_sx + num1) - 1];
            }
        }
         
        return num3;
    }

    public static void swap(int n, double[] dx, int ix, double[] dy, int iy) throws Exception {
        for (int num1 = 0;num1 < n;num1++)
        {
            double num2 = dx[num1 + ix];
            dx[num1 + ix] = dy[num1 + iy];
            dy[num1 + iy] = num2;
        }
    }

    public static void swap(int n, double[][] dx, int dxRow, int ix, double[][] dy, int dyRow, int iy) throws Exception {
        for (int num1 = 0;num1 < n;num1++)
        {
            double num2 = dx[dxRow][num1 + ix];
            dx[dxRow][num1 + ix] = dy[dyRow][num1 + iy];
            dy[dyRow][num1 + iy] = num2;
        }
    }

    public static void trsv(char uplo, char trans, char diag, int n, double[] a, int offset_a, int lda, double[] x, int offset_x, int incx) throws Exception {
        if (n != 0)
        {
            int num1;
            int num2;
            boolean flag1 = (diag == 'N') || (diag == 'n');
            boolean flag3 = (uplo == 'U') || (uplo == 'u');
            boolean flag2 = (((trans == 'T') || (trans == 't')) || (trans == 'C')) || (trans == 'c');
            if (flag3)
            {
                if (flag2)
                {
                    if (incx > 0)
                    {
                        num2 = 0;
                        for (num1 = 0;num1 < n;num1++)
                        {
                            x[offset_x + num2] -= BLAS.dot(num1 - 1, a, offset_a + lda * num1, 1, x, offset_x, incx);
                            if (flag1)
                            {
                                x[offset_x + num2] /= a[offset_a + num1 + lda * num1];
                            }
                             
                            num2 += incx;
                        }
                    }
                    else
                    {
                        num2 = ((-n + 1) * incx);
                        for (num1 = 0;num1 < n;num1++)
                        {
                            x[offset_x + num2] -= BLAS.dot(num1, a, offset_a + lda * num1, 1, x, offset_x + num2 - incx, incx);
                            if (flag1)
                            {
                                x[offset_x + num2] /= a[offset_a + num1 + lda * num1];
                            }
                             
                            num2 += incx;
                        }
                    } 
                }
                else if (incx > 0)
                {
                    num2 = ((n - 1) * incx);
                    for (num1 = n;num1 >= 1;num1--)
                    {
                        if (num1 < n)
                        {
                            x[offset_x + num2] -= BLAS.dot(n - num1, a, ((offset_a + num1) + (lda * num1)) - 1, lda, x, offset_x + num2 + incx, incx);
                        }
                         
                        if (flag1)
                        {
                            x[offset_x + num2] /= a[((offset_a + num1) + (lda * (num1 - 1))) - 1];
                        }
                         
                        num2 -= incx;
                    }
                }
                else
                {
                    num2 = 0;
                    for (num1 = n;num1 >= 1;num1--)
                    {
                        if (num1 < n)
                        {
                            x[offset_x + num2] -= BLAS.dot(n - num1, a, ((offset_a + num1) + (lda * num1)) - 1, lda, x, offset_x, incx);
                        }
                         
                        if (flag1)
                        {
                            x[offset_x + num2] /= a[((offset_a + num1) + (lda * (num1 - 1))) - 1];
                        }
                         
                        num2 -= incx;
                    }
                }  
            }
            else if (flag2)
            {
                if (incx > 0)
                {
                    num2 = ((n - 1) * incx);
                    for (num1 = n;num1 >= 1;num1--)
                    {
                        if (num1 < n)
                        {
                            x[offset_x + num2] -= BLAS.dot(n - num1, a, (offset_a + num1) + (lda * (num1 - 1)), 1, x, offset_x + num2 + incx, incx);
                        }
                         
                        if (flag1)
                        {
                            x[offset_x + num2] /= a[((offset_a + num1) + (lda * (num1 - 1))) - 1];
                        }
                         
                        num2 -= incx;
                    }
                }
                else
                {
                    num2 = 0;
                    for (num1 = n;num1 >= 1;num1--)
                    {
                        if (num1 < n)
                        {
                            x[offset_x + num2] -= BLAS.dot(n - num1, a, (offset_a + num1) + (lda * (num1 - 1)), 1, x, offset_x, incx);
                        }
                         
                        if (flag1)
                        {
                            x[offset_x + num2] /= a[((offset_a + num1) + (lda * (num1 - 1))) - 1];
                        }
                         
                        num2 -= incx;
                    }
                } 
            }
            else if (incx > 0)
            {
                num2 = 0;
                for (num1 = 0;num1 < n;num1++)
                {
                    x[offset_x + num2] -= BLAS.dot(num1, a, offset_a + num1, lda, x, offset_x, incx);
                    if (flag1)
                    {
                        x[offset_x + num2] /= a[offset_a + num1 + lda * num1];
                    }
                     
                    num2 += incx;
                }
            }
            else
            {
                num2 = ((-n + 1) * incx);
                for (num1 = 0;num1 < n;num1++)
                {
                    x[offset_x + num2] -= BLAS.dot(num1, a, offset_a + num1, lda, x, offset_x + num2 - incx, incx);
                    if (flag1)
                    {
                        x[offset_x + num2] /= a[offset_a + num1 + lda * num1];
                    }
                     
                    num2 += incx;
                }
            }   
        }
         
    }

    public static double xyz(int n, double[] sx, int incx, double[] sy, int incy, double[] sz, int incz) throws Exception {
        double num5 = 0;
        if (n > 0)
        {
            int num1;
            if (((incx != 1) || (incy != 1)) || (incz != 1))
            {
                int num2 = 1;
                int num3 = 1;
                int num4 = 1;
                if (incx < 0)
                {
                    num2 = ((-n + 1) * incx) + 1;
                }
                 
                if (incy < 0)
                {
                    num3 = ((-n + 1) * incy) + 1;
                }
                 
                if (incz < 0)
                {
                    num4 = ((-n + 1) * incz) + 1;
                }
                 
                for (num1 = 1;num1 <= n;num1++)
                {
                    num5 += (sx[num2 - 1] * sy[num3 - 1]) * sz[num4 - 1];
                    num2 += incx;
                    num3 += incy;
                    num4 += incz;
                }
                return num5;
            }
             
            for (num1 = 1;num1 <= n;num1++)
            {
                num5 += (sx[num1 - 1] * sy[num1 - 1]) * sz[num1 - 1];
            }
        }
         
        return num5;
    }

    public static double xyz(int n, double[] sx, int xoffset, int incx, double[] sy, int yoffset, int incy, double[] sz, int zoffset, int incz) throws Exception {
        double num5 = 0;
        if (n > 0)
        {
            int num1;
            if (((incx != 1) || (incy != 1)) || (incz != 1))
            {
                int num2 = 1;
                int num3 = 1;
                int num4 = 1;
                if (incx < 0)
                {
                    num2 = ((-n + 1) * incx) + 1;
                }
                 
                if (incy < 0)
                {
                    num3 = ((-n + 1) * incy) + 1;
                }
                 
                if (incz < 0)
                {
                    num4 = ((-n + 1) * incz) + 1;
                }
                 
                for (num1 = 1;num1 <= n;num1++)
                {
                    num5 += (sx[(xoffset + num2) - 1] * sy[(yoffset + num3) - 1]) * sz[(zoffset + num4) - 1];
                    num2 += incx;
                    num3 += incy;
                    num4 += incz;
                }
                return num5;
            }
             
            for (num1 = 1;num1 <= n;num1++)
            {
                num5 += (sx[(xoffset + num1) - 1] * sy[(yoffset + num1) - 1]) * sz[(zoffset + num1) - 1];
            }
        }
         
        return num5;
    }

    public static double xyz(int n, double[][] sx, int sxrow, int xoffset, int incx, double[][] sy, int syrow, int yoffset, int incy, double[][] sz, int szrow, int zoffset, int incz) throws Exception {
        double num5 = 0;
        if (n > 0)
        {
            int num1;
            if (((incx != 1) || (incy != 1)) || (incz != 1))
            {
                int num2 = 1;
                int num3 = 1;
                int num4 = 1;
                if (incx < 0)
                {
                    num2 = ((-n + 1) * incx) + 1;
                }
                 
                if (incy < 0)
                {
                    num3 = ((-n + 1) * incy) + 1;
                }
                 
                if (incz < 0)
                {
                    num4 = ((-n + 1) * incz) + 1;
                }
                 
                for (num1 = 1;num1 <= n;num1++)
                {
                    num5 += (sx[sxrow][(xoffset + num2) - 1] * sy[syrow][(yoffset + num3) - 1]) * sz[szrow][(zoffset + num4) - 1];
                    num2 += incx;
                    num3 += incy;
                    num4 += incz;
                }
                return num5;
            }
             
            for (num1 = 1;num1 <= n;num1++)
            {
                num5 += (sx[sxrow][(xoffset + num1) - 1] * sy[syrow][(yoffset + num1) - 1]) * sz[szrow][(zoffset + num1) - 1];
            }
        }
         
        return num5;
    }

    public static double xyz(int n, double[][][] sx, int sxRow, int sxCol, int xoffset, int incx, double[][][] sy, int syRow, int syCol, int yoffset, int incy, double[][] sz, int szRow, int zoffset, int incz) throws Exception {
        double num5 = 0;
        if (n > 0)
        {
            int num1;
            if (((incx != 1) || (incy != 1)) || (incz != 1))
            {
                int num2 = 1;
                int num3 = 1;
                int num4 = 1;
                if (incx < 0)
                {
                    num2 = ((-n + 1) * incx) + 1;
                }
                 
                if (incy < 0)
                {
                    num3 = ((-n + 1) * incy) + 1;
                }
                 
                if (incz < 0)
                {
                    num4 = ((-n + 1) * incz) + 1;
                }
                 
                for (num1 = 1;num1 <= n;num1++)
                {
                    num5 += (sx[sxRow][sxCol][(xoffset + num2) - 1] * sy[syRow][syCol][(yoffset + num3) - 1]) * sz[szRow][(zoffset + num4) - 1];
                    num2 += incx;
                    num3 += incy;
                    num4 += incz;
                }
                return num5;
            }
             
            for (num1 = 1;num1 <= n;num1++)
            {
                num5 += (sx[sxRow][sxCol][(xoffset + num1) - 1] * sy[syRow][syCol][(yoffset + num1) - 1]) * sz[szRow][(zoffset + num1) - 1];
            }
        }
         
        return num5;
    }

    private static final double[] _threshold = new double[]{ 1.0010415475916E-146, 4.4989137945432E+161, 2.2227587494851E-162, 1.9979190722022E+146, 5.0104209000224E-293, 1.9958403095347E+292 };
}


