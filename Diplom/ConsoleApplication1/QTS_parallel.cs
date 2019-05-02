using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Numerics;
using MathWorks.MATLAB.NET.Utility;
using MathWorks.MATLAB.NET.Arrays;
using Screen_forward_parallel;

namespace ConsoleApplication1
{
    class QTS_parallel
    {
        public Complex[][] MatrixR, screenEQ;
        public double pi = 3.14;
        public double k1, k0;
        public double Rad = 9e-6;
        public double lambda0 = 600E-9, lambda1 = 500E-9;
        public int N, M, K;
        public double steprho = 0.05e-6, stepz = 0.05e-6;
        public double[] pointsrho, pointsz, screenrho, screenx, screeny, screenz;
        public Complex coef;
        public double omega, mu1 = 1, eps1 = 1, mu2 = 1, eps2 = 1, eps0 = 8.85E-12, mu0 = 1.25E-6;
        public double ssquare;
        public Complex[] xx,xx2, Ephi, coeff;
        public double[] jphi;
        public double T;
        public double a =  2e-6, b = 1e-6, c= 0.03e-6;
        public Complex[] screenE;
        public int numsteps = 300;
        double dist = 10e-6;
        public double[][] ScreenMatrix;
        //public Complex alpha = 1e-20;
        public double dist_dip = 5e-6;
        public Complex[][] IMatrix;
        public int dim=50;
        public Complex[][] Green, GreenScreen;
        public double sled=0.01e-6;
        public int num = 1;
        double alp;
        public Complex[][] MatrixSum(Complex[][] A, Complex[][] B, int N)
        {
            Complex[][] C;
            C = new Complex[N][];
            for (int i = 0; i < N; i++)
            {
                C[i] = new Complex[N];
            }
            for (int i = 0; i < A.Length; i++)
            {
                for (int j = 0; j < A[i].Length; j++)
                    C[i][j] = A[i][j] + B[i][j];
            }
            return C;
        }

        public Complex[][] MatrixMul(Complex[][] A, Complex[][] B, int N)
        {
            Complex[][] C;
            C = new Complex[N][];
            for (int i = 0; i < N; i++)
            {
                C[i] = new Complex[N];
            }
            for (int i = 0; i < A.Length; i++)
            {
                for (int j = 0; j < B[i].Length; j++)
                {
                    C[i][j] = 0;
                    for (int k = 0; k < A[i].Length; k++)
                    {
                        C[i][j] += A[i][k] * B[k][j];
                    }
                }
                    
            }
            return C;
        }

        public Complex MatrixDet(Complex[][] matrix)
        {
            if (matrix.Length == 2)
            {
                
                return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0];
            }
            Complex sign = 1, result = 0;
            for (int i = 0; i < matrix.GetLength(0); i++)
            {
                Complex[][] minor = GetMinor(matrix, i, 0);
                result += sign * matrix[0][i] * MatrixDet(minor);
                
                sign = -sign;
            }
            return result;
        }

        public Complex[][] GetMinor(Complex[][] matrix, int n, int m)
        {
            Complex[][] result = new Complex[matrix.GetLength(0)-1][];
            for (int i = 0; i < matrix.GetLength(0)-1; i++)
            {
                result[i] = new Complex[matrix.GetLength(0) - 1];
            }
            for (int i = 0, temp=0; i < matrix.GetLength(0); i++)
            {
                if (i == m)
                {
                    temp = 1;
                    continue;
                }
                    
                for (int j = 0, col = 0; j < matrix.GetLength(0); j++)
                {
                    if (j == n)
                        continue;
                    result[i-temp][col] = matrix[i][j];
                    col++;
                }
            }
            return result;
        }

        public Complex[][] ReverseMatrix(Complex[][] matrix)
        {
            Complex[][] result = new Complex[matrix.GetLength(0)][];
            for (int i = 0; i < matrix.GetLength(0); i++)
            {
                result[i] = new Complex[matrix.GetLength(0)];
            }
            Complex sign = 1;
            Complex oper = MatrixDet(matrix);
            for (int i = 0; i < matrix.GetLength(0); i++)
            {
                for (int j = 0; j < matrix.GetLength(0); j++)
                {
                    result[i][j] = sign * MatrixDet(GetMinor(matrix, i, j))/oper;
                    //Console.Write("{0} ", result[i][j]);
                    sign = -sign;
                }
                //Console.WriteLine();
            }

            return result;
        }

        public Complex[][] ConMatrix(Complex[][] matrix)
        {
            Complex[][] result = new Complex[matrix.GetLength(0)][];
            for (int i = 0; i < matrix.GetLength(0); i++)
            {
                result[i] = new Complex[matrix.GetLength(0)];
            }
            for (int i = 0; i < matrix.GetLength(0); i++)
            {
                for (int j = 0; j < matrix.GetLength(0); j++)
                {
                    result[i][j] = Complex.Conjugate(matrix[j][i]);
                }
                Console.WriteLine();
            }
            return result;
        }

        public double DipoleFuncC(double rho, double z)
        {
            //double R = Math.Sqrt(rho * rho + z * z);
            //double Rr = rho / R;
            //double cos = Math.Cos(k0 * R), sin = Math.Sin(k0 * R);
            //double Result = -T * (Rr * R * k0 * cos - Rr * sin) / (R * R);
            //R = Math.Sqrt((rho - sled) * (rho - sled) + z * z);
            //Rr = (rho - sled) / R;
            //cos = Math.Cos(k0 * R); sin = Math.Sin(k0 * R);
            //Result += -T * (Rr * R * k0 * cos - Rr * sin) / (R * R);
            //R = Math.Sqrt((rho + sled) * (rho + sled) + z * z);
            //Rr = (rho + sled) / R;
            //cos = Math.Cos(k0 * R); sin = Math.Sin(k0 * R);
            //Result += -T * (Rr * R * k0 * cos - Rr * sin) / (R * R);
            double Result = 0;
            for (int i = 0; i < num; i++)
            {
                double R = Math.Sqrt((rho + i * sled) * (rho + i * sled) + z * z);
                double Rr = (rho + i * sled) / R;
                double cos = Math.Cos(k0 * R), sin = Math.Sin(k0 * R);
                Result += -T * (Rr * R * k0 * cos - Rr * sin) / (R * R);
            }
            return Result;
            //return T * sin / R;
        }

        public double DipoleFuncR(double rho, double z)
        {

            //R = Math.Sqrt((rho-sled) * (rho-sled) + z * z);
            //Rr = (rho-sled) / R;
            //cos = Math.Cos(k0 * R); sin = Math.Sin(k0 * R);
            //Result += -T * (-Rr * R * k0 * sin - Rr * cos) / (R * R);
            //R = Math.Sqrt((rho + sled) * (rho + sled) + z * z);
            //Rr = (rho + sled) / R;
            //cos = Math.Cos(k0 * R); sin = Math.Sin(k0 * R);
            //Result += -T * (-Rr * R * k0 * sin - Rr * cos) / (R * R);
            double Result = 0;
            for (int i=-(num-1)/2; i<=(num-1)/2; i++)
            {
                double R = Math.Sqrt((rho + i*sled) * (rho + i*sled) + z * z);
                double Rr = (rho + i*sled) / R;
                double cos = Math.Cos(k0 * R), sin = Math.Sin(k0 * R);
                Result += -T * (-Rr * R * k0 * sin - Rr * cos) / (R * R);
            }
            return Result;
            //return T * cos / R;
        }

        

        public double GreenFuncR(double rho1, double z1, double rho2, double z2, double k)
        {
            double IR = 0;
            double R;
            double start = 1e-3;
            double step = (2 * Math.PI) / 500;
            double cos, ecos;
            double Rho1=rho1, Rho2=rho2, Z1=z1, Z2=z2;
            //if ((Math.Abs(Rho1-Rho2)<= steprho)&&(Math.Abs(Z1-Z2)<=stepz))
            //{
                Rho2+= steprho / 10;
            //    Z2 += stepz / 2;
            //}
            for (double i = 0; i <= 2 * Math.PI; i += step)
            {
                cos = Math.Cos(i);
                R = Math.Sqrt(Rho1 * Rho1 + Rho2 * Rho2 - 2 * Rho1 * Rho2 * cos + Math.Pow(Z1 - Z2, 2));
                ecos = Math.Cos(k * R);
                IR += ecos * cos / R;
            }
            IR *= step;
            return IR;
        }

        public double GreenFuncC(double rho1, double z1, double rho2, double z2, double k)
        {
            double IC = 0;
            double R;
            double start = 1e-3;
            double step = (2 * Math.PI) / 500;
            double cos, esin;
            double Rho1 = rho1, Rho2 = rho2, Z1 = z1, Z2 = z2;
            //if ((Math.Abs(Rho1 - Rho2) <= steprho) && (Math.Abs(Z1 - Z2) <= stepz))
            //{
                Rho2 += steprho / 10;
            //    Z2 += stepz / 2;
            //}
            for (double i = 0; i <= 2*Math.PI; i += step)
            {
                cos = Math.Cos(i);
                R = Math.Sqrt(Rho1 * Rho1 + Rho2 * Rho2 - 2 * Rho1 * Rho2 * cos + (Z1-Z2)* (Z1 - Z2));
                esin = Math.Sin(k * R);
                IC += esin * cos / R;
            }
            IC *= step;
            return IC;
        }

        public bool IsInCircuit(double rho, double z)
        {
            return (Math.Sqrt(rho * rho + (a - z) * (a - z)) <= a);
        }

        public bool IsInEllipse(double rho, double z)
        {
            return (Math.Sqrt(rho * rho / (a * a) + (b - z) * (b - z) / (b * b)) <= 1);
        }

        public bool IsInCilinder(double rho, double z)
        {
            return ((rho < a) && (z < 2 * b));
        }

        public void InitJphi()
        {
            jphi = new double[K];
            for (int i = 0; i < K; i++)
            {
                jphi[i] = 1e-12;
            }
        }

        public void InitMatrix()
        {
            MatrixR = new Complex[K][];
            for (int i = 0; i < K; i++)
            {
                MatrixR[i] = new Complex[K + 1];
            }
            coef = (double)(k1 * k1 - k0 * k0) / (double)(4 * Math.PI);
            Complex GR = 0; //Complex GI = 0;
            Parallel.For (0, K, i =>
            {
                for (int j = 0; j < K; j++)
                {
                    Complex GI = 0;
                    GI = coef * new Complex(GreenFuncR(pointsrho[i], pointsz[i], pointsrho[j], pointsz[j], k0), GreenFuncC(pointsrho[i], pointsz[i], pointsrho[j], pointsz[j],k0)) * ssquare;
                    if (i == j)
                    {
                        MatrixR[i][j] = 1 - pointsrho[i] * GI;
                        //Console.WriteLine("InitMatrix2[{0}][{1}] {2}", i, j, MatrixR[i][j]);
                        //Console.WriteLine("1-InitMatrix2[{0}][{1}] {2}", i, j, pointsrho[i] * GI);
                        //Console.WriteLine(MatrixR[i][j]);
                    }
                    else
                    {
                        MatrixR[i][j] = -pointsrho[i] * GI;
                    }
                }

                //MatrixR[i][K] = new Complex(DipoleIntegralR(i) / (double)(4 * Math.PI), DipoleIntegralC(i) / (double)(4 * Math.PI));
                MatrixR[i][K] = new Complex(DipoleFuncR(pointsrho[i], pointsz[i]), DipoleFuncC(pointsrho[i], pointsz[i]));
                //Console.WriteLine(MatrixR[i][K]);
            });
            System.IO.StreamWriter file1 = new System.IO.StreamWriter(@"C:\Users\User\Desktop\parallel\Matrix_Contour_0.1.txt");
            for (int i = 0; i < K; i++)
            {
                for (int j = 0; j < K; j++)
                {
                    file1.Write((MatrixR[i][j]).Real.ToString().Replace(",", "."));
                    file1.Write(" ");
                }
                for (int j = 0; j < K; j++)
                {
                    file1.Write((MatrixR[i][j]).Imaginary.ToString().Replace(",", "."));
                    file1.Write(" ");
                }
                file1.Write("\n");
            }
            file1.WriteLine();
            file1.Close();
            System.IO.StreamWriter file2 = new System.IO.StreamWriter(@"C:\Users\User\Desktop\parallel\Right_Contour_0.1.txt");
            for (int i = 0; i < K; i++)
            {
                file2.Write((MatrixR[i][K]).Real.ToString().Replace(",", "."));
                file2.Write(" ");
            }
            for (int i = 0; i < K; i++)
            {
                file2.Write((MatrixR[i][K]).Imaginary.ToString().Replace(",", "."));
                file2.Write(" ");
            }

            file2.WriteLine();
            file2.Close();
            //System.IO.StreamWriter file3 = new System.IO.StreamWriter(@"C:\Users\User\Desktop\parallel\DipoleR_Contour_0.1.txt");
            //for (int i = 0; i < K; i++)
            //{
            //    //file3.Write((Math.Sqrt(GreenFuncR(pointsrho[0], pointsz[0], pointsrho[i], pointsz[i]) * GreenFuncR(pointsrho[0], pointsz[0], pointsrho[i], pointsz[i]) + GreenFuncC(pointsrho[0], pointsz[0], pointsrho[i], pointsz[i]) * GreenFuncC(pointsrho[0], pointsz[0], pointsrho[i], pointsz[i]))).ToString().Replace(",", "."));
            //    file3.Write((Math.Sqrt(DipoleFuncR(pointsrho[i], pointsz[i]) * DipoleFuncR(pointsrho[i], pointsz[i]) + DipoleFuncC(pointsrho[i], pointsz[i]) * DipoleFuncC(pointsrho[i], pointsz[i]))).ToString().Replace(",", "."));
            //    file3.Write(" ");
            //}

            //file3.WriteLine();
            //file3.Close();

        }

        public void CountNumberMas()
        {   // что-то добавить в конструктор
            M = (int)(2*a / steprho);
            N = (int)(2 * b / stepz);
            N = Math.Max(M, N);
            //Console.WriteLine(N);
            //Console.WriteLine(M);
            k1 = 2 * pi / lambda1;
            k0 = 2 * pi / lambda0;
            ssquare = stepz * steprho;
            omega = k1 / Math.Sqrt(eps0 * mu0);
            int n = 0;
            double rhomid, zmid;
            T = 1e-9 / (4 * 3.14);
            for (int i = 0; i < N; i++)
            {
                for (int j = 0; j < N; j++)
                {
                    rhomid = i * steprho + steprho / 2;
                    zmid = j * stepz + stepz / 2;
                    //if (IsInCircuit(rhomid, zmid))
                    //    n++;
                    //if (IsInEllipse(rhomid, zmid))
                    ////    n++;
                    if (IsInCilinder(rhomid, zmid))
                        n++;
                }
            }
            K = n;
            xx = new Complex[n];
            pointsrho = new double[n];
            pointsz = new double[n];
            n = 0;
            for (int i = 0; i < M; i++)
            {
                rhomid = i * steprho + steprho / 2;
                for (int j = 0; j < N; j++)
                {

                    zmid = j * stepz + stepz / 2;
                    //if (IsInCircuit(rhomid, zmid))
                    //{
                    //    pointsrho[n] = rhomid;
                    //    zmid = zmid + dist_dip;
                    //    //pointsrho[n] = Math.Sqrt(rhomid * rhomid + zmid * zmid);
                    //    pointsz[n] = zmid;
                    //    //Console.WriteLine(pointsz[n]);
                    //    n++;
                    //}
                    //if (IsInEllipse(rhomid, zmid))
                    //{
                    //    pointsrho[n] = rhomid;
                    //    zmid = zmid + 1e-6;
                    //    //pointsrho[n] = Math.Sqrt(rhomid * rhomid + zmid * zmid);
                    //    pointsz[n] = zmid;
                    //    n++;
                    //}
                    if (IsInCilinder(rhomid, zmid))
                    {
                        pointsrho[n] = rhomid;
                        zmid = zmid + dist_dip;
                        //pointsrho[n] = Math.Sqrt(rhomid * rhomid + zmid * zmid);
                        pointsz[n] = zmid;
                        n++;
                    }
                }
            }
            System.IO.StreamWriter file2 = new System.IO.StreamWriter(@"C:\Users\User\Desktop\parallel\Contour_x_0.1.txt");
            for (int i = 0; i < K; i++)
            {
                file2.Write(pointsrho[i].ToString().Replace(",", "."));
                //file2.Write((Math.Sqrt(IQ.pointsrho[i]* IQ.pointsrho[i]- IQ.pointsz[i]* IQ.pointsz[i])).ToString().Replace(",", "."));
                file2.Write(" ");
            }
            System.IO.StreamWriter file3 = new System.IO.StreamWriter(@"C:\Users\User\Desktop\parallel\Contour_y_0.1.txt");
            file2.Close();
            for (int i = 0; i < K; i++)
            {
                file3.Write(pointsz[i].ToString().Replace(",", "."));
                file3.Write(" ");
            }
            file3.Close();

        }

        public Complex[] Gauss()
        {
            Complex tmp;
            //прямой ход
            //Console.WriteLine(MatrixR[K - 1][K - 1]);
            xx = new Complex[K];
            
            for (int i = 0; i < K; i++)
            {
                tmp = MatrixR[i][i];
                //Console.WriteLine("Gauss tmp {0}", tmp);
                for (int j = K; j >= i; j--)
                {
                    //Console.WriteLine("Gauss Matrix[{0}][{1}] {2}", i, j, MatrixR[i][j]);
                    MatrixR[i][j] /= tmp;
                    if (MatrixR[i][j].Real==Double.NaN)
                    {
                        while(true)
                        {
                            Console.Write('!');
                        }
                    }
                    //Console.WriteLine("Gauss Matrix[{0}][{1}] {2}", i, j, MatrixR[i][j]);
                }
                for (int j=i + 1; j<K; j++)
                {
                   tmp = MatrixR[j][i];
                   //Console.WriteLine("Gauss tmp {0}", tmp);
                   for (int k = K; k >= i; k--)
                   {
                       MatrixR[j][k] -= tmp * MatrixR[i][k];
                   }
                   //Console.WriteLine("Gauss Matrix[{0}][{1}] {2}", i,j, MatrixR[j][i]);
                }

                //Console.WriteLine(i);
            }
            /*обратный ход*/
            xx[K - 1] = MatrixR[K - 1][K];
            
            for (int i = K - 2; i >= 0; i--)
            {
                xx[i] = MatrixR[i][K];
                for (int j=i + 1; j<K; j++)
                    { xx[i] -= MatrixR[i][j] * xx[j]; }
                
            }
            System.IO.StreamWriter filec = new System.IO.StreamWriter(@"C:\Users\User\Desktop\parallel\xx.txt");
            for (int i = 0; i < K; i++)
            {
                filec.Write((xx[i].Real).ToString().Replace(",", "."));
                filec.Write(" ");
            }
            for (int i = 0; i < K; i++)
            {
                filec.Write((xx[i].Imaginary).ToString().Replace(",", "."));
                filec.Write(" ");
            }
            filec.Close();
            return xx;
        }


        public void InitScreen()
        {

            screenx = new double[2*numsteps];
            screeny = new double[2*numsteps];
            screenz = new double[2 * numsteps];
            screenrho = new double[2 * numsteps];

            screenE = new Complex[2 * numsteps];
            screenEQ = new Complex[2 * numsteps][];
            for (int i=0; i<2*numsteps; i++)
            {
                screenEQ[i] = new Complex[2 * numsteps];
            }

            //screenz = new double[numsteps];
            //screenrho = new double[numsteps];

            //screenE = new Complex[numsteps];


            double screendim = steprho * numsteps;
            for (int i = 0; i < numsteps; i++)
            {
                screenz[i] = dist_dip+dist;
                screenrho[i] = 8 * c * (i);
            }
            //for (int i = numsteps; i < 2*numsteps; i++)
            //{
            //    screenz[i] = -(dist_dip + dist);
            //    screenrho[i] = 8 * c * (i - 3*numsteps / 2);
            //}
            //for (int i = 0; i < numsteps; i++)
            //{
            //    screenz[i] = dist;
            //    screenrho[i] = 4 * steprho * (i);
            //}
            //for (int i = 2 * numsteps; i < 3 * numsteps; i++)
            //{
            //    screenz[i] =  dist_dip + 2*stepz * (i - 2 * numsteps);
            //    screenrho[i] = dist;
            //}
            //for (int i = 3 * numsteps; i < 4 * numsteps; i++)
            //{
            //    screenz[i] = dist_dip + stepz * (i - 3 * numsteps);
            //    screenrho[i] = -dist;
            //}
            for (int i = 0; i < 2 * numsteps; i++)
            {
                screenx[i] = screeny[i] = 8 * c * (i - numsteps);
            }
        }

        public void ScreenE()
        {
            Parallel.For(0, numsteps, i =>
            {
                Complex IR = 0, IC = 0;
                for (int j = 0; j < K; j++)
                {
                    IC = IC + (new Complex(GreenFuncR(screenrho[i], screenz[i], pointsrho[j], pointsz[j], k0), GreenFuncC(screenrho[i], screenz[i], pointsrho[j], pointsz[j], k0))/*+ new Complex(GreenFuncR(screenrho[i], screenz[i], -pointsrho[j], pointsz[j], k0), GreenFuncC(screenrho[i], screenz[i], -pointsrho[j], pointsz[j], k0))*/) * xx[j];
                }
                IC = screenrho[i] * coef * IC * ssquare;
                
                screenE[i] = new Complex(DipoleFuncR(screenrho[i], screenz[i]), DipoleFuncC(screenrho[i], screenz[i])) + IC;
                IR = IC = 0;

            });
            System.IO.StreamWriter file5 = new System.IO.StreamWriter(@"C:\Users\User\Desktop\parallel\x.txt", false);
            for (int i = 0; i < numsteps; i++)
            {
                file5.Write(screenrho[i].ToString().Replace(",", "."));
                file5.Write(" ");
            }
            file5.Close();
            System.IO.StreamWriter file4 = new System.IO.StreamWriter(@"C:\Users\User\Desktop\parallel\Screen.txt", false);
            for (int i = 0; i < numsteps; i++)
            {
                file4.Write((Math.Sqrt(screenE[i].Real * screenE[i].Real + screenE[i].Imaginary * screenE[i].Imaginary)).ToString().Replace(",", "."));
                file4.Write(" ");
                //Console.WriteLine((Math.Sqrt(screenE[i].Real * screenE[i].Real + screenE[i].Imaginary * screenE[i].Imaginary)).ToString().Replace(",", "."));
            }
            file4.Close();



            //Parallel.For(0, numsteps, i =>
            //{
            //    Complex IR = 0, IC = 0;
            //    for (int j = 0; j < numsteps; j++)
            //    {
            //        double screenxy = Math.Sqrt(screenx[i + numsteps] * screenx[i + numsteps] + screeny[j + numsteps] * screeny[j + numsteps]);
            //        for (int k = 0; k < K; k++)
            //        {
            //            IC = IC + (new Complex(GreenFuncR(screenxy, screenz[i + numsteps], pointsrho[k], pointsz[k], k0), GreenFuncC(screenxy, screenz[i + numsteps], pointsrho[k], pointsz[k], k0))/* + new Complex(GreenFuncR(screenxy, screenz[i], -pointsrho[k], pointsz[k], k0), GreenFuncC(screenxy, screenz[i], -pointsrho[k], pointsz[k], k0))*/) * xx[k];
            //        }
            //        IC = screenxy * coef * IC * ssquare;

            //        screenEQ[i + numsteps][j + numsteps] = screenEQ[i + numsteps][numsteps - j - 1] = screenEQ[numsteps - i - 1][j + numsteps] = screenEQ[numsteps - i - 1][numsteps - j - 1] = new Complex(DipoleFuncR(screenxy, screenz[i]), DipoleFuncC(screenxy, screenz[i])) + IC;
            //        IR = IC = 0;
            //    }


            //});
            //file5 = new System.IO.StreamWriter(@"C:\Users\User\Desktop\parallel\xy.txt", false);
            //for (int i = 0; i < 2 * numsteps; i++)
            //{
            //    file5.Write(screenx[i].ToString().Replace(",", "."));
            //    file5.Write(" ");
            //}
            //file5.Close();
            //file4 = new System.IO.StreamWriter(@"C:\Users\User\Desktop\parallel\ScreenQ.txt", false);
            //for (int i = 0; i < 2 * numsteps; i++)
            //{
            //    for (int j = 0; j < 2 * numsteps; j++)
            //    {
            //        file4.Write((Math.Sqrt(screenEQ[i][j].Real * screenEQ[i][j].Real + screenEQ[i][j].Imaginary * screenEQ[i][j].Imaginary)).ToString().Replace(",", "."));
            //        file4.Write(" ");
            //    }
            //    file4.Write("\n");
            //    //Console.WriteLine((Math.Sqrt(screenE[i].Real * screenE[i].Real + screenE[i].Imaginary * screenE[i].Imaginary)).ToString().Replace(",", "."));
            //}
            //file4.Close();

        }


        //обратная задача
        public void IterationProcessKOLD(double alpha)
        {
            double steprho_f = (Rad) / (dim);
            double stepz_f = (Rad) / dim;
            double rhomid, zmid;
            int n = 0;
            steprho = steprho_f;
            stepz = stepz_f;
            pointsrho = new double[2 * dim * dim];
            pointsz = new double[2 * dim * dim];
            //double[] E_phi_0_R = new double[dim * dim], E_phi_0_C = new double[dim * dim];
            //double[] E_phi_R = new double[dim * dim], E_phi_C = new double[dim * dim];
            double[] K_1 = new double[2 * dim * dim], K_0 = new double[2 * dim * dim], ksi = new double[2 * dim * dim];
            for (int i = 0; i < dim; i++)
            {
                rhomid = i * steprho_f + steprho_f / 2;
                for (int j = 0; j < dim; j++)
                {
                    zmid = j * stepz_f + stepz_f / 2;
                    pointsrho[n] = rhomid;
                    //zmid = zmid; //+ dist_dip;
                    pointsz[n] = dist_dip /*5*stepz_f*/ + zmid;
                    n++;
                }
            }
            K = n;
            System.IO.StreamWriter file3 = new System.IO.StreamWriter(@"C:\Users\User\Desktop\parallel\Forward_x.txt",false);
            for (int i = 0; i < n; i++)
            {
                file3.Write(pointsrho[i].ToString().Replace(",", "."));
                //file2.Write((Math.Sqrt(IQ.pointsrho[i]* IQ.pointsrho[i]- IQ.pointsz[i]* IQ.pointsz[i])).ToString().Replace(",", "."));
                file3.Write(" ");
            }
            System.IO.StreamWriter file4 = new System.IO.StreamWriter(@"C:\Users\User\Desktop\parallel\Forward_y.txt",false);
            file3.Close();
            for (int i = 0; i < n; i++)
            {
                file4.Write(pointsz[i].ToString().Replace(",", "."));
                file4.Write(" ");
            }
            file4.Close();
            //numsteps = numsteps * 2;
            MatrixR = new Complex[numsteps][];
            for (int i = 0; i < numsteps; i++)
            {
                MatrixR[i] = new Complex[K + 1];
            }
            Complex[][] MatrixRInter = new Complex[K][];
            for (int i = 0; i < K; i++)
            {
                MatrixRInter[i] = new Complex[numsteps];
            }
            Green = new Complex[K][];
            for (int i = 0; i < K; i++)
            {
                Green[i] = new Complex[K];
            }
            GreenScreen = new Complex[numsteps][];
            for (int i = 0; i < numsteps; i++)
            {
                GreenScreen[i] = new Complex[K + 1];
            }
            Parallel.For(0, numsteps, i =>
            {
                Complex GR = 0, GI = 0;
                double DR = 0, DC = 0;
                double k = k0;
                for (int j = 0; j < K; j++)
                {
                    GreenScreen[i][j] = new Complex(GreenFuncR(screenrho[i], screenz[i], pointsrho[j], pointsz[j], k), GreenFuncC(screenrho[i], screenz[i], pointsrho[j], pointsz[j], k));
                    
                }
            });
            Parallel.For(0, K, i =>
            {
                Complex GR = 0, GI = 0;
                double DR = 0, DC = 0;
                double k = k0;
                for (int j = 0; j < K; j++)
                {
                    Green[i][j] = new Complex(GreenFuncR(pointsrho[i], pointsz[i], pointsrho[j], pointsz[j], k), GreenFuncC(pointsrho[i], pointsz[i], pointsrho[j], pointsz[j], k));

                }
            });
            


            
            Parallel.For(0, numsteps, i =>
            {
                double square = steprho_f * stepz_f;
                Complex GR = 0, GI = 0;
                double DR = 0, DC = 0;
                for (int j = 0; j < K; j++)
                {
                    GI = screenrho[i] * GreenScreen[i][j]*square;
                    DR = DipoleFuncR(pointsrho[j], pointsz[j]);
                    DC = DipoleFuncC(pointsrho[j], pointsz[j]);
                    MatrixR[i][j] = GI * new Complex(DR, DC);
                    MatrixRInter[j][i] = Complex.Conjugate(MatrixR[i][j]);
                    //Console.WriteLine(MatrixR[i][j]);
                }
                //MatrixR[i][K] = -new Complex(DipoleIntegralRS(i), DipoleIntegralCS(i)) / (double)(4 * Math.PI) - screenE[i];
                MatrixR[i][K] = (-new Complex(DipoleFuncR(screenrho[i], screenz[i]), DipoleFuncC(screenrho[i], screenz[i])) + screenE[i]);
                //MatrixRInter[K][i] = Complex.Conjugate(MatrixR[i][K]);
            });
            IMatrix = new Complex[K][];
            for (int i = 0; i < K; i++)
            {
                IMatrix[i] = new Complex[K + 1];
            }
            for (int i = 0; i < K; i++)
            {
                for (int j = 0; j < K; j++)
                {
                    IMatrix[i][j] = 0;
                    for (int k = 0; k < numsteps; k++)
                    {
                        IMatrix[i][j] += MatrixRInter[i][k] * MatrixR[k][j];
                    }

                }
            }

            for (int i = 0; i < K; i++)
            {
                IMatrix[i][K] = 0;
                for (int k = 0; k < numsteps; k++)
                {
                    IMatrix[i][K] += MatrixRInter[i][k]*MatrixR[k][K];
                }
            }
            alp= alpha * Complex.Abs(IMatrix[0][0]);
            for (int i = 0; i < K; i++)
            {
                IMatrix[i][i] = IMatrix[i][i] + alpha * Complex.Abs(IMatrix[0][0]);
            }
            
            //Complex min = IMatrix[0][0];
            //for (int i = 0; i < K; i++)
            //{
            //    for (int j = 0; j < K; j++)
            //    {
            //        if (Complex.Abs(IMatrix[i][j]) > Complex.Abs(min))
            //            min = IMatrix[i][j];
            //    }
            //}
            //for (int i = 0; i < K; i++)
            //{
            //    for (int j = 0; j <= K; j++)
            //    {
            //        IMatrix[i][j] /= min;
            //    }
            //}
            System.IO.StreamWriter file1 = new System.IO.StreamWriter(@"C:\Users\User\Desktop\parallel\Forward_Matrix.txt",false);
            for (int i = 0; i < K; i++)
            {
                for (int j = 0; j < K; j++)
                {
                    file1.Write((IMatrix[i][j]).Real.ToString().Replace(",", "."));
                    file1.Write(" ");
                }
                for (int j = 0; j < K; j++)
                {
                    file1.Write((IMatrix[i][j]).Imaginary.ToString().Replace(",", "."));
                    file1.Write(" ");
                }
                file1.Write("\n");
            }
            file1.WriteLine();
            file1.Close();
            System.IO.StreamWriter file2 = new System.IO.StreamWriter(@"C:\Users\User\Desktop\parallel\Forward_Right.txt",false);
            for (int i = 0; i < K; i++)
            {
                file2.Write((IMatrix[i][K].Real).ToString().Replace(",", "."));
                file2.Write(" ");
            }
            for (int i = 0; i < K; i++)
            {
                file2.Write((IMatrix[i][K].Imaginary).ToString().Replace(",", "."));
                file2.Write(" ");
            }

            file2.WriteLine();
            file2.Close();

            coeff = new Complex[K];
            for (int i = 0; i < K; i++)
            {
                coeff[i] = 0;
            }
        }

        public Complex[] MatlabGauss()
        {
            xx = new Complex[K];
            MWArray[] res = null; //выходной массив метода plane
            MWNumericArray descriptor = null; //массив возвращаемого параметра 
            Class1 obj_plane = new Class1();
            res = obj_plane.plot_dielectr_parallel(1);
            descriptor = (MWNumericArray)res[0]; //выбор первого элемента из массива MWArray и преобразование в числовой тип MWNumericArray
            double[,] d_r = (double[,])descriptor.ToArray(MWArrayComponent.Real);//преобразование массива MWNUmericArray  к масииву типа double
            double[,] d_c = (double[,])descriptor.ToArray(MWArrayComponent.Imaginary);//преобразование массива MWNUmericArray  к масииву типа double
            for (int i = 0; i < K; i++)
            {
                xx[i] = new Complex(d_r[i, 0], d_c[i, 0]);
                //Console.WriteLine(xx[i]);
            }
            System.IO.StreamWriter filec = new System.IO.StreamWriter(@"C:\Users\User\Desktop\parallel\xx.txt");
            for (int i = 0; i < K; i++)
            {
                filec.Write((xx[i].Real).ToString().Replace(",", "."));
                filec.Write(" ");
            }
            for (int i = 0; i < K; i++)
            {
                filec.Write((xx[i].Imaginary).ToString().Replace(",", "."));
                filec.Write(" ");
            }
            filec.Close();
            return xx;
        }

        public Complex[] MatlabGaussI()
        {
            MWArray[] res = null; //выходной массив метода plane
            MWNumericArray descriptor = null; //массив возвращаемого параметра 
            Class1 obj_plane = new Class1();
            res = obj_plane.Screen_forward_parallel(1,dim);
            descriptor = (MWNumericArray)res[0]; //выбор первого элемента из массива MWArray и преобразование в числовой тип MWNumericArray
            double[,] d_r = (double[,])descriptor.ToArray(MWArrayComponent.Real);//преобразование массива MWNUmericArray  к масииву типа double
            double[,] d_c = (double[,])descriptor.ToArray(MWArrayComponent.Imaginary);//преобразование массива MWNUmericArray  к масииву типа double
            for (int i= 0;i< K; i++)
            {
                coeff[i] = new Complex(d_r[i, 0], d_c[i, 0]);
                //Console.WriteLine("coeff {0}",coeff[i]);
            }
            return coeff;
        }

        


        public void InitMatrix2()
        {
            MatrixR = new Complex[K][];
            for (int i = 0; i < K; i++)
            {
                MatrixR[i] = new Complex[K + 1];
            }
            //coef = (double)(k1 * k1 - k0 * k0) / (double)(4 * Math.PI);
            
            Parallel.For(0, K, i =>
            {
                Complex GR = 0, GI = 0;
                for (int j = 0; j < K; j++)
                {
                    //Console.WriteLine("GreenFunc [{0}][{1}] {2}", i, j, GR);
                    GI = coeff[j] * Green[i][j] * ssquare;
                    if (i == j)
                    {
                        MatrixR[i][j] = 1 - pointsrho[i] * GI;
                        //Console.WriteLine("InitMatrix2[{0}][{1}] {2}", i, j, MatrixR[i][j]);
                        //Console.WriteLine("1-InitMatrix2[{0}][{1}] {2}", i, j, pointsrho[i] * GI);
                    }
                    else
                    {
                        MatrixR[i][j] = -pointsrho[i] * GI;
                        //Console.WriteLine("InitMatrix2[{0}][{1}] {2}", i, j, MatrixR[i][j]);
                    }
                    //Console.ReadLine();
                    //Console.WriteLine("InitMatrix2[{0}][{1}] {2}", i,j,MatrixR[i][j]);
                }

                MatrixR[i][K] = new Complex(DipoleFuncR(pointsrho[i], pointsz[i]), DipoleFuncC(pointsrho[i], pointsz[i]));
                //Console.WriteLine(MatrixR[i][K]);
            });
            System.IO.StreamWriter file1 = new System.IO.StreamWriter(@"C:\Users\User\Desktop\parallel\Matrix_Contour_0.1.txt");
            for (int i = 0; i < K; i++)
            {
                for (int j = 0; j < K; j++)
                {
                    file1.Write((MatrixR[i][j]).Real.ToString().Replace(",", "."));
                    file1.Write(" ");
                }
                for (int j = 0; j < K; j++)
                {
                    file1.Write((MatrixR[i][j]).Imaginary.ToString().Replace(",", "."));
                    file1.Write(" ");
                }
                file1.Write("\n");
            }
            file1.WriteLine();
            file1.Close();
            System.IO.StreamWriter file2 = new System.IO.StreamWriter(@"C:\Users\User\Desktop\parallel\Right_Contour_0.1.txt");
            for (int i = 0; i < K; i++)
            {
                file2.Write((MatrixR[i][K]).Real.ToString().Replace(",", "."));
                file2.Write(" ");
            }
            for (int i = 0; i < K; i++)
            {
                file2.Write((MatrixR[i][K]).Imaginary.ToString().Replace(",", "."));
                file2.Write(" ");
            }

            file2.WriteLine();
            file2.Close();
            System.IO.StreamWriter file4 = new System.IO.StreamWriter(@"C:\Users\User\Desktop\parallel\Contour_x_0.1.txt");
            for (int i = 0; i < K; i++)
            {
                file4.Write(pointsrho[i].ToString().Replace(",", "."));
                //file2.Write((Math.Sqrt(IQ.pointsrho[i]* IQ.pointsrho[i]- IQ.pointsz[i]* IQ.pointsz[i])).ToString().Replace(",", "."));
                file4.Write(" ");
            }
            System.IO.StreamWriter file3 = new System.IO.StreamWriter(@"C:\Users\User\Desktop\parallel\Contour_y_0.1.txt");
            file4.Close();
            for (int i = 0; i < K; i++)
            {
                file3.Write(pointsz[i].ToString().Replace(",", "."));
                file3.Write(" ");
            }
            file3.Close();
            //filec.WriteLine();
            //filec.Close();


        }

        public void IterationProcess2()
        {
            double steprho_f = (Rad) / (dim);
            double stepz_f = (Rad) / dim;
            int n = 0;
            MatrixR = new Complex[numsteps][];
            for (int i = 0; i < numsteps; i++)
            {
                MatrixR[i] = new Complex[K + 1];
            }
            
            
            Parallel.For(0, numsteps, i =>
            {
                double square = steprho_f * stepz_f;
                Complex GR = 0, GI = 0;
                Complex DR = 0, DC = 0;
                for (int j = 0; j < K; j++)
                {
                    GI = screenrho[i] * GreenScreen[i][j] * square;
                    DR = xx[j];
                    //DC = DipoleIntegralC(j) / (double)(4 * Math.PI);
                    MatrixR[i][j] = GI * DR;
                    //Console.WriteLine(DR);
                }
                MatrixR[i][K] =  (-new Complex(DipoleFuncR(screenrho[i], screenz[i]), DipoleFuncC(screenrho[i], screenz[i])) + screenE[i]);
            });
            IMatrix = new Complex[K][];
            for (int i = 0; i < K; i++)
            {
                IMatrix[i] = new Complex[K + 1];
            }
            for (int i = 0; i < K; i++)
            {
                for (int j = 0; j < K; j++)
                {
                    IMatrix[i][j] = 0;
                    for (int k = 0; k < numsteps; k++)
                    {
                        IMatrix[i][j] += Complex.Conjugate(MatrixR[k][i]) * MatrixR[k][j];

                    }
                }
            }

            for (int i = 0; i < K; i++)
            {
                IMatrix[i][K] = 0;
                for (int k = 0; k < numsteps; k++)
                {
                    IMatrix[i][K] += Complex.Conjugate(MatrixR[k][i]) * MatrixR[k][K];
                }
            }
            for (int i = 0; i < K; i++)
            {
                IMatrix[i][i] += alp;
            }
            System.IO.StreamWriter file1 = new System.IO.StreamWriter(@"C:\Users\User\Desktop\parallel\Forward_Matrix.txt", false);
            for (int i = 0; i < K; i++)
            {
                for (int j = 0; j < K; j++)
                {
                    file1.Write((IMatrix[i][j]).Real.ToString().Replace(",", "."));
                    file1.Write(" ");
                }
                for (int j = 0; j < K; j++)
                {
                    file1.Write((IMatrix[i][j]).Imaginary.ToString().Replace(",", "."));
                    file1.Write(" ");
                }
                file1.Write("\n");
            }
            file1.WriteLine();
            file1.Close();
            System.IO.StreamWriter file2 = new System.IO.StreamWriter(@"C:\Users\User\Desktop\parallel\Forward_Right.txt", false);
            for (int i = 0; i < K; i++)
            {
                file2.Write((IMatrix[i][K].Real).ToString().Replace(",", "."));
                file2.Write(" ");
            }
            for (int i = 0; i < K; i++)
            {
                file2.Write((IMatrix[i][K].Imaginary).ToString().Replace(",", "."));
                file2.Write(" ");
            }

            file2.WriteLine();
            file2.Close();
            xx = new Complex[K];

        }



    }
}
