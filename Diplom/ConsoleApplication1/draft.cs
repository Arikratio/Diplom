using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ConsoleApplication1
{
    class draft
    {
        //public double DipoleIntegralR(int i)
        //{
        //    double I = 0;
        //    for (int j = 0; j < K; j++)
        //    {
        //        I += GreenFuncR(pointsrho[i], pointsz[i], pointsrho[j], pointsz[j]) * DipoleFuncR(pointsrho[j], pointsz[j]) - GreenFuncC(pointsrho[i], pointsz[i], pointsrho[j], pointsz[j]) * DipoleFuncC(pointsrho[j], pointsz[j]);
        //    }
        //    return I * pointsrho[i] * (ssquare);// * pointsrho[j];
        //}

        //public double DipoleIntegralC(int i)
        //{
        //    double I = 0;
        //    for (int j = 0; j < K; j++)
        //    {
        //        I += GreenFuncC(pointsrho[i], pointsz[i], pointsrho[j], pointsz[j]) * DipoleFuncR(pointsrho[j], pointsz[j]) + GreenFuncR(pointsrho[i], pointsz[i], pointsrho[j], pointsz[j]) * DipoleFuncC(pointsrho[j], pointsz[j]);
        //    }
        //    return I * pointsrho[i] * (ssquare);// * pointsrho[j];
        //}

        //public double DipoleIntegralRS(int i)
        //{
        //    double I = 0;
        //    for (int j = 0; j < K; j++)
        //    {
        //        I += GreenFuncR(screenrho[i], screenz[i], pointsrho[j], pointsz[j]) * DipoleFuncR(pointsrho[j], pointsz[j]) - GreenFuncC(screenrho[i], screenz[i], pointsrho[j], pointsz[j]) * DipoleFuncC(pointsrho[j], pointsz[j]);
        //    }
        //    return I * screenrho[i] * (ssquare);// * pointsrho[j];
        //}

        //public double DipoleIntegralCS(int i)
        //{
        //    double I = 0;
        //    for (int j = 0; j < K; j++)
        //    {
        //        I += GreenFuncC(screenrho[i], screenz[i], pointsrho[j], pointsz[j]) * DipoleFuncR(pointsrho[j], pointsz[j]) + GreenFuncR(screenrho[i], screenz[i], pointsrho[j], pointsz[j]) * DipoleFuncC(pointsrho[j], pointsz[j]);
        //    }
        //    return I * screenrho[i] * (ssquare);// * pointsrho[j];
        //}

        //public double DipoleFuncC(double rho, double z)
        //{
        //    double R = Math.Sqrt(rho * rho + z * z);
        //    double Rr = rho / (Math.Sqrt(rho * rho + z * z));
        //    double cos = Math.Cos(k0 * R), sin = Math.Sin(k0 * R);
        //    return T * (Rr * R * k0 * cos - Rr * sin) / (R * R);
        //}

        //public double DipoleFuncR(double rho, double z)
        //{
        //    double R = Math.Sqrt(rho * rho + z * z);
        //    double Rr = rho / (Math.Sqrt(rho * rho + z * z));
        //    double cos = Math.Cos(k0 * R), sin = Math.Sin(k0 * R);
        //    return T * (-Rr * R * k0 * sin - Rr * cos) / (R * R);
        //}

        //public double DipoleIntegralR(int i)
        //{
        //    double I = 0;
        //    var lockObj = new object();
        //    Parallel.For(0, K, j =>
        //    {
        //        lock (lockObj)
        //        {
        //            I += GreenFuncR(pointsrho[i], pointsz[i], pointsrho[j], pointsz[j]) * DipoleFuncR(pointsrho[j], pointsz[j]) - GreenFuncC(pointsrho[i], pointsz[i], pointsrho[j], pointsz[j]) * DipoleFuncC(pointsrho[j], pointsz[j]);
        //        }
        //    });
        //    return I * pointsrho[i] * (ssquare);// * pointsrho[j];
        //}

        //public double DipoleIntegralC(int i)
        //{
        //    double I = 0;
        //    var lockObj = new object();
        //    Parallel.For(0, K, j =>
        //    {
        //        lock (lockObj)
        //        {
        //            I += GreenFuncC(pointsrho[i], pointsz[i], pointsrho[j], pointsz[j]) * DipoleFuncR(pointsrho[j], pointsz[j]) + GreenFuncR(pointsrho[i], pointsz[i], pointsrho[j], pointsz[j]) * DipoleFuncC(pointsrho[j], pointsz[j]);
        //        }
        //    });
        //    return I * pointsrho[i] * (ssquare);// * pointsrho[j];
        //}

        //public double DipoleIntegralRS(int i)
        //{
        //    double I = 0;
        //    var lockObj = new object();
        //    Parallel.For(0, K, j =>
        //    {
        //        lock (lockObj)
        //        {
        //            I += GreenFuncR(screenrho[i], screenz[i], pointsrho[j], pointsz[j]) * DipoleFuncR(pointsrho[j], pointsz[j]) - GreenFuncC(screenrho[i], screenz[i], pointsrho[j], pointsz[j]) * DipoleFuncC(pointsrho[j], pointsz[j]);
        //        }
        //    });
        //    return I * screenrho[i] * (ssquare);// * pointsrho[j];
        //}

        //public double DipoleIntegralCS(int i)
        //{
        //    double I = 0;
        //    var lockObj = new object();
        //    Parallel.For(0, K, j =>
        //    {
        //        lock (lockObj)
        //        {
        //            I += GreenFuncC(screenrho[i], screenz[i], pointsrho[j], pointsz[j]) * DipoleFuncR(pointsrho[j], pointsz[j]) + GreenFuncR(screenrho[i], screenz[i], pointsrho[j], pointsz[j]) * DipoleFuncC(pointsrho[j], pointsz[j]);
        //        }
        //    });
        //    return I * screenrho[i] * (ssquare);// * pointsrho[j];
        //}

        //public double GreenFuncR(double rho1, double z1, double rho2, double z2)
        //{
        //    var lockObj = new object();
        //    double step = 0.06;
        //    double IR = 0;
        //    int n_pi = (int)Math.Round(2*Math.PI / step); //количество шагов от 0 до 2*pi
        //    Parallel.For(0, n_pi, i =>
        //    {

        //        double R;
        //        double cos, ecos;
        //        cos = Math.Cos(i*step);
        //        R = Math.Sqrt(rho1 * rho1 + rho2 * rho2 - 2 * rho1 * rho2 * cos + Math.Pow(z1 - z2, 2));
        //        ecos = Math.Cos(k0 * R);
        //        if (R != 0)
        //        {
        //            lock (lockObj)
        //            {
        //                IR = IR + ecos / R * cos;
        //            }
        //        }
        //    });
        //    IR = IR * step;
        //    return IR;
        //}

        //public double GreenFuncC(double rho1, double z1, double rho2, double z2)
        //{
        //    double IC = 0;
        //    var lockObj = new object();
        //    double step = 0.06;
        //    int n_pi = (int)Math.Round(2*Math.PI / step); //количество шагов от 0 до 2*pi
        //    Parallel.For(0, n_pi, i =>
        //    {
        //        double R;
        //        double cos, esin;
        //        cos = Math.Cos(i*step);
        //        R = Math.Sqrt(rho1 * rho1 + rho2 * rho2 - 2 * rho1 * rho2 * cos + Math.Pow(z1 - z2, 2));
        //        esin = Math.Sin(k0 * R);
        //        if (R != 0)
        //        {
        //            lock (lockObj)
        //            {
        //                IC = IC + esin / R * cos;
        //            }
        //        }

        //    });
        //    IC = IC * step;
        //    return IC;
        //}

        //public Complex[] ZeidelI()
        //{
        //    Complex[] xx0 = new Complex[K], xx1 = new Complex[K];
        //    for (int i = 0; i < K; i++)
        //    {
        //        xx0[i] = coef;
        //        xx1[i] = 0;
        //    }
        //    Complex[][] MatrixC = new Complex[K][];
        //    for (int i = 0; i < K; i++)
        //    {
        //        MatrixC[i] = new Complex[K + 1];
        //    }
        //    Parallel.For(0, K, i =>
        //    {
        //        for (int j = 0; j <= K; j++)
        //        {
        //            if (i != j)
        //                MatrixC[i][j] = -IMatrix[i][j] / IMatrix[i][i];
        //            else
        //                MatrixC[i][j] = 0;
        //        }
        //    });
        //    while (true)
        //    {

        //        for (int i = 0; i < K; i++)
        //        {
        //            for (int j = 0; j < i; j++)
        //            {
        //                xx1[i] += MatrixC[i][j] * xx1[j];
        //            }
        //            for (int j = i; j < K; j++)
        //            {
        //                xx1[i] += MatrixC[i][j] * xx0[j];
        //            }
        //            xx1[i] += MatrixC[i][K];
        //        }
        //        double delta = 0;
        //        for (int i = 0; i < K; i++)
        //        {
        //            delta += Complex.Abs(xx1[i] - xx0[i]);
        //        }
        //        if (delta < 1e11)
        //        {
        //            for (int i = 0; i < K; i++)
        //            {
        //                coeff[i] = xx1[i];
        //            }
        //            break;
        //        }
        //        Console.WriteLine(delta);
        //        for (int i = 0; i < K; i++)
        //        {
        //            xx0[i] = xx1[i];
        //            xx1[i] = 0;
        //        }
        //    }
        //    System.IO.StreamWriter filec = new System.IO.StreamWriter(@"C:\Users\User\Desktop\parallel\Coeff.txt");
        //    for (int i = 0; i < K; i++)
        //    {
        //        filec.Write((coeff[i].Real).ToString().Replace(",", "."));
        //        filec.Write(" ");
        //    }
        //    for (int i = 0; i < K; i++)
        //    {
        //        filec.Write((coeff[i].Imaginary).ToString().Replace(",", "."));
        //        filec.Write(" ");
        //    }
        //    filec.Close();
        //    return coeff;
        //}

        //public Complex[] GaussI()
        //{
        //    Complex tmp;
        //    //прямой ход
        //    //Console.WriteLine(MatrixR[K - 1][K - 1]);
        //    xx = new Complex[K];
        //    xx2 = new Complex[K];
        //    Complex[][] IMatrixTemp = new Complex[K][];
        //    for (int i = 0; i < K; i++)
        //    {
        //        IMatrixTemp[i] = new Complex[K + 1];
        //    }
        //    for (int i = 0; i < K; i++)
        //    {
        //        for (int j = 0; j <= K; j++)
        //        {
        //            IMatrixTemp[i][j] = IMatrix[i][j];
        //        }
        //    }
        //    for (int i = K - 1; i >= 0; i--)
        //    {
        //        tmp = IMatrixTemp[i][i];
        //        //Console.WriteLine("Gauss tmp {0}", tmp);
        //        for (int j = 0; j < i; j++)
        //        {
        //            //Console.WriteLine("Gauss Matrix[{0}][{1}] {2}", i, j, MatrixR[i][j]);
        //            IMatrixTemp[i][j] /= tmp;
        //            if (IMatrixTemp[i][j].Real == Double.NaN)
        //            {
        //                while (true)
        //                {
        //                    Console.Write('!');
        //                }
        //            }
        //            //Console.WriteLine("Gauss Matrix[{0}][{1}] {2}", i, j, MatrixR[i][j]);
        //        }
        //        for (int j = i - 1; j >= 0; j--)
        //        {
        //            tmp = IMatrixTemp[j][i];
        //            //Console.WriteLine("Gauss tmp {0}", tmp);
        //            for (int k = 0; k <= i; k++)
        //            {
        //                IMatrixTemp[j][k] -= tmp * IMatrixTemp[i][k];
        //            }

        //        }

        //    }
        //    //for (int i = 0; i < K; i++)
        //    //{
        //    //    for (int j = 0; j <= K; j++)
        //    //    {
        //    //        Console.WriteLine("Gauss Matrix[{0}][{1}] {2}", i, j, IMatrix[i][j]);
        //    //        //Console.WriteLine(i);
        //    //        /*обратный ход*/

        //    //    }
        //    //}
        //    xx[K - 1] = IMatrixTemp[K - 1][K];
        //    for (int i = K - 2; i >= 0; i--)
        //    {
        //        xx[i] = IMatrixTemp[i][K];
        //        for (int j = i + 1; j < K; j++)
        //        { xx[i] -= IMatrixTemp[i][j] * xx[j]; }

        //    }
        //    //System.IO.StreamWriter filec = new System.IO.StreamWriter(@"C:\Users\User\Desktop\parallel\xx.txt");
        //    //for (int i = 0; i < K; i++)
        //    //{
        //    //    filec.Write((xx[i].Real).ToString().Replace(",", "."));
        //    //    filec.Write(" ");
        //    //}
        //    //for (int i = 0; i < K; i++)
        //    //{
        //    //    filec.Write((xx[i].Imaginary).ToString().Replace(",", "."));
        //    //    filec.Write(" ");
        //    //}
        //    //filec.Close();


        //    //Complex[] b1 = new Complex[K];
        //    //for (int i = 0; i < K; i++)
        //    //{
        //    //    b1[i] = 0;
        //    //    for (int j = 0; j < K; j++)
        //    //    {
        //    //        b1[i] += IMatrix[i][j] * xx[j];
        //    //    }
        //    //}
        //    //for (int i = 0; i < K; i++)
        //    //{
        //    //    for (int j = 0; j < K; j++)
        //    //    {
        //    //        IMatrixTemp[i][j] = IMatrix[i][j];
        //    //    }
        //    //    IMatrixTemp[i][K] = IMatrix[i][K]-b1[i];
        //    //}
        //    //for (int i = 0; i < K; i++)
        //    //{
        //    //    tmp = IMatrixTemp[i][i];
        //    //    //Console.WriteLine("Gauss tmp {0}", tmp);
        //    //    for (int j = K; j >= i; j--)
        //    //    {
        //    //        //Console.WriteLine("Gauss Matrix[{0}][{1}] {2}", i, j, MatrixR[i][j]);
        //    //        IMatrixTemp[i][j] /= tmp;
        //    //        if (IMatrixTemp[i][j].Real == Double.NaN)
        //    //        {
        //    //            while (true)
        //    //            {
        //    //                Console.Write('!');
        //    //            }
        //    //        }
        //    //        //Console.WriteLine("Gauss Matrix[{0}][{1}] {2}", i, j, MatrixR[i][j]);
        //    //    }
        //    //    for (int j = i + 1; j < K; j++)
        //    //    {
        //    //        tmp = IMatrixTemp[j][i];
        //    //        //Console.WriteLine("Gauss tmp {0}", tmp);
        //    //        for (int k = K; k >= i; k--)
        //    //        {
        //    //            IMatrixTemp[j][k] -= tmp * IMatrixTemp[i][k];
        //    //        }

        //    //    }



        //    //}
        //    //xx2[K - 1] = IMatrixTemp[K - 1][K];
        //    //for (int i = K - 2; i >= 0; i--)
        //    //{
        //    //    xx2[i] = IMatrixTemp[i][K];
        //    //    for (int j = i + 1; j < K; j++)
        //    //    { xx2[i] -= IMatrixTemp[i][j] * xx2[j]; }

        //    //}
        //    System.IO.StreamWriter filec = new System.IO.StreamWriter(@"C:\Users\User\Desktop\parallel\xx.txt");
        //    for (int i = 0; i < K; i++)
        //    {
        //        filec.Write(((xx[i]).Real).ToString().Replace(",", "."));
        //        filec.Write(" ");
        //    }
        //    for (int i = 0; i < K; i++)
        //    {
        //        filec.Write(((xx[i]).Imaginary).ToString().Replace(",", "."));
        //        filec.Write(" ");
        //    }
        //    filec.Close();

        //    return xx;



        //    //Complex tmp;
        //    //Complex[] xx2 = new Complex[K];
        //    //for (int i = 0; i < K; i++)
        //    //{
        //    //    xx2[i] = coeff[i];
        //    //}
        //    ////прямой ход
        //    //int Maxi = 0, p = 0;
        //    ////double Max = Complex.Abs(IMatrix[p][p]);
        //    ////for (int j = p; j < K; j++)
        //    ////{
        //    ////    for (int k = 0; k < K; k++)
        //    ////    {
        //    ////        if (Complex.Abs(IMatrix[j][k]) > Max)
        //    ////        {
        //    ////            Max = Complex.Abs(IMatrix[j][k]);
        //    ////            Maxi = j;
        //    ////        }
        //    ////    }
        //    ////}
        //    ////for (int j = 0; j <= K; j++)
        //    ////{
        //    ////    tmp = IMatrix[p][j];
        //    ////    IMatrix[p][j] = IMatrix[Maxi][j];
        //    ////    IMatrix[Maxi][j] = tmp;
        //    ////}
        //    ////p++;

        //    ////Console.WriteLine(IMatrix[i][K]);
        //    ////System.IO.StreamWriter file1 = new System.IO.StreamWriter(@"C:\Users\User\Desktop\parallel\Forward_Matrix.txt", false);
        //    ////for (int i = 0; i < K; i++)
        //    ////{
        //    ////    for (int j = 0; j < K; j++)
        //    ////    {
        //    ////        file1.Write((IMatrix[i][j]).ToString().Replace(",", "."));
        //    ////        file1.Write("    ");
        //    ////    }
        //    ////}



        //    //for (int i = 0; i < K; i++)
        //    //{


        //    //    for (int j = i + 1; j < K; j++)
        //    //    {
        //    //        //tmp = IMatrix[j][i];
        //    //        for (int k = K; k >= i; k--)
        //    //        {
        //    //            IMatrix[j][k] = IMatrix[j][k] - IMatrix[j][i]/ IMatrix[i][i]* IMatrix[i][k];
        //    //            //IMatrix[j][k] -= tmp * IMatrix[i][k];
        //    //        }

        //    //    }
        //    //    tmp = IMatrix[i][i];
        //    //    for (int j = K; j >= i; j--)
        //    //    { IMatrix[i][j] /= tmp; }
        //    //    Console.WriteLine("IMatrix last {0}", IMatrix[K - 1][K]);
        //    //}

        //    ///*обратный ход*/
        //    //coeff[K - 1] = IMatrix[K - 1][K];
        //    //for (int i = K - 2; i >= 0; i--)
        //    //{
        //    //    coeff[i] = IMatrix[i][K];
        //    //    for (int j = i + 1; j < K; j++)
        //    //    { coeff[i] -= IMatrix[i][j] * coeff[j]; }
        //    //    Console.WriteLine("GaussI coeff[{0}] {1}", i, Complex.Abs(2*Math.PI/Complex.Sqrt(coeff[i]*4*Math.PI+k0*k0)));
        //    //}

        //    //System.IO.StreamWriter filec = new System.IO.StreamWriter(@"C:\Users\User\Desktop\parallel\Coeff.txt");
        //    ////for (int i = 0; i < K; i++)
        //    ////{
        //    ////    filec.Write((2 * Math.PI / Math.Sqrt((Complex.Abs(coeff[i]) * 4 * Math.PI + k0 * k0))).ToString().Replace(",", "."));
        //    ////    filec.Write(" ");
        //    ////}
        //    //for (int i = 0; i < K; i++)
        //    //{
        //    //    filec.Write((coeff[i].Real).ToString().Replace(",", "."));
        //    //    filec.Write(" ");
        //    //}
        //    //for (int i = 0; i < K; i++)
        //    //{
        //    //    filec.Write((coeff[i].Imaginary).ToString().Replace(",", "."));
        //    //    filec.Write(" ");
        //    //}
        //    //filec.Close();

        //    //return coeff;

        //}
    }
}
