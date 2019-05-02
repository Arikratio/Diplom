using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Numerics;

namespace ConsoleApplication1
{
    class QTS
    {
        public Complex[][] MatrixR;
        public double pi = 3.14;
        public double k1, k0;
        public double Rad = 3E-6;
        public double lambda0 = 600E-9, lambda1 = 450E-9;
        public int N, M, K;
        public double steprho = 0.1E-6, stepz = 0.1E-6;
        public double[] pointsrho, pointsz, screenrho, screenx, screeny, screenz;
        public Complex coef;
        public double omega, mu1 = 1, eps1 = 1, mu2 = 1, eps2 = 1, eps0 = 8.85E-12, mu0 = 1.25E-6;
        public double ssquare;
        public Complex[] xx, coeff;
        public double[] jphi;
        public double T;
        public double a = 3e-6, b = 3e-6;
        public Complex[] screenE;
        public int numsteps = 15;
        double dist = 7.5E-6;
        public double[][] ScreenMatrix;
        public double alpha = 1e-50;
        public double dist_dip = 1e-6;
        public Complex[][] IMatrix;


        public double DipoleFuncC(double rho, double z)
        {
            double R = Math.Sqrt(rho * rho + z * z);
            double Rr = rho / (Math.Sqrt(rho * rho + z * z));
            double cos = Math.Cos(k0 * R), sin = Math.Sin(k0 * R);
            return T * (Rr * R * k0 * cos - Rr * sin) / (R * R);
        }

        public double DipoleFuncR(double rho, double z)
        {
            double R = Math.Sqrt(rho * rho + z * z);
            double Rr = rho / (Math.Sqrt(rho * rho + z * z));
            double cos = Math.Cos(k0 * R), sin = Math.Sin(k0 * R);
            return T * (-Rr * R * k0 * sin - Rr * cos) / (R * R);
        }

        public double DipoleIntegralR(int i)
        {
            double I = 0;
            for (int j = 0; j < K; j++)
            {
                I += GreenFuncR(pointsrho[i], pointsz[i], pointsrho[j], pointsz[j]) * DipoleFuncR(pointsrho[j], pointsz[j]) - GreenFuncC(pointsrho[i], pointsz[i], pointsrho[j], pointsz[j]) * DipoleFuncC(pointsrho[j], pointsz[j]);
            }
            return I * pointsrho[i] * (ssquare);// * pointsrho[j];
        }

        public double DipoleIntegralC(int i)
        {
            double I = 0;
            for (int j = 0; j < K; j++)
            {
                I += GreenFuncC(pointsrho[i], pointsz[i], pointsrho[j], pointsz[j]) * DipoleFuncR(pointsrho[j], pointsz[j]) + GreenFuncR(pointsrho[i], pointsz[i], pointsrho[j], pointsz[j]) * DipoleFuncC(pointsrho[j], pointsz[j]);
            }
            return I * pointsrho[i] * (ssquare);// * pointsrho[j];
        }

        public double DipoleIntegralRS(int i)
        {
            double I = 0;
            for (int j = 0; j < K; j++)
            {
                I += GreenFuncR(screenrho[i], screenz[i], pointsrho[j], pointsz[j]) * DipoleFuncR(pointsrho[j], pointsz[j]) - GreenFuncC(screenrho[i], screenz[i], pointsrho[j], pointsz[j]) * DipoleFuncC(pointsrho[j], pointsz[j]);
            }
            return I * screenrho[i] * (ssquare);// * pointsrho[j];
        }

        public double DipoleIntegralCS(int i)
        {
            double I = 0;
            for (int j = 0; j < K; j++)
            {
                I += GreenFuncC(screenrho[i], screenz[i], pointsrho[j], pointsz[j]) * DipoleFuncR(pointsrho[j], pointsz[j]) + GreenFuncR(screenrho[i], screenz[i], pointsrho[j], pointsz[j]) * DipoleFuncC(pointsrho[j], pointsz[j]);
            }
            return I * screenrho[i] * (ssquare);// * pointsrho[j];
        }

        public double GreenFuncR(double rho1, double z1, double rho2, double z2)
        {
            double IR = 0;
            double R;
            double step = 0.06;
            double cos, ecos;
            for (double i = 0; i < 6.28; i = i + step)
            {
                cos = Math.Cos(i);
                R = Math.Sqrt(rho1 * rho1 + rho2 * rho2 - 2 * rho1 * rho2 * cos + Math.Pow(z1 - z2, 2));
                ecos = Math.Cos(k0 * R);
                if (R != 0)
                    IR = IR + ecos / R * cos;
            }
            IR = IR * step;
            return IR;
        }

        public double GreenFuncC(double rho1, double z1, double rho2, double z2)
        {
            double IC = 0;
            double R;
            double step = 0.06;
            double cos, esin;
            for (double i = 0; i < 6.28; i += step)
            {
                cos = Math.Cos(i);
                R = Math.Sqrt(rho1 * rho1 + rho2 * rho2 - 2 * rho1 * rho2 * cos + Math.Pow(z1 - z2, 2));
                esin = Math.Sin(k0 * R);
                if (R != 0)
                    IC = IC + esin / R * cos;
            }
            IC = IC * step;
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
            Complex GR = 0, GI = 0;
            for (int i = 0; i < K; i++)
            {
                for (int j = 0; j < K; j++)
                {

                    GI = coef * new Complex(GreenFuncR(pointsrho[i], pointsz[i], pointsrho[j], pointsz[j]), GreenFuncC(pointsrho[i], pointsz[i], pointsrho[j], pointsz[j])) * ssquare;
                    if (i == j)
                    {
                        MatrixR[i][j] = 1 - pointsrho[i]*GI;
                        //Console.WriteLine(MatrixR[i][j]);
                    }
                    else
                    {
                        MatrixR[i][j] = -pointsrho[i]*GI;
                    }
                }

                MatrixR[i][K] = -new Complex(DipoleIntegralR(i) / (double)(4 * Math.PI), DipoleIntegralC(i) / (double)(4 * Math.PI));
                //Console.WriteLine(MatrixR[i][K]);
            }
            System.IO.StreamWriter file1 = new System.IO.StreamWriter(@"C:\Users\User\Desktop\данные\Matrix_Contour_0.1.txt");
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
            System.IO.StreamWriter file2 = new System.IO.StreamWriter(@"C:\Users\User\Desktop\данные\Right_Contour_0.1.txt");
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
            System.IO.StreamWriter file3 = new System.IO.StreamWriter(@"C:\Users\User\Desktop\данные\DipoleR_Contour_0.1.txt");
            for (int i = 0; i < K; i++)
            {
                //file3.Write((Math.Sqrt(GreenFuncR(pointsrho[0], pointsz[0], pointsrho[i], pointsz[i]) * GreenFuncR(pointsrho[0], pointsz[0], pointsrho[i], pointsz[i]) + GreenFuncC(pointsrho[0], pointsz[0], pointsrho[i], pointsz[i]) * GreenFuncC(pointsrho[0], pointsz[0], pointsrho[i], pointsz[i]))).ToString().Replace(",", "."));
                file3.Write((Math.Sqrt(DipoleFuncR(pointsrho[i], pointsz[i]) * DipoleFuncR(pointsrho[i], pointsz[i]) + DipoleFuncC(pointsrho[i], pointsz[i]) * DipoleFuncC(pointsrho[i], pointsz[i]))).ToString().Replace(",", "."));
                file3.Write(" ");
            }

            file3.WriteLine();
            file3.Close();

        }

        public void CountNumberMas()
        {   // что-то добавить в конструктор
            M = (int)(a / steprho);
            N = (int)(2 * b / stepz);
            //Console.WriteLine(N);
            //Console.WriteLine(M);
            k1 = 2 * pi / lambda1;
            k0 = 2 * pi / lambda0;
            ssquare = stepz * steprho;
            omega = k1 / Math.Sqrt(eps0 * mu0);
            int n = 0;
            double rhomid, zmid;
            T = 1e-12 / (4 * 3.14);
            for (int i = 0; i < M; i++)
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
                    //    Console.WriteLine(pointsz[n]);
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
            System.IO.StreamWriter file2 = new System.IO.StreamWriter(@"C:\Users\User\Desktop\данные\Contour_x_0.1.txt");
            for (int i = 0; i < K; i++)
            {
                file2.Write(pointsrho[i].ToString().Replace(",", "."));
                //file2.Write((Math.Sqrt(IQ.pointsrho[i]* IQ.pointsrho[i]- IQ.pointsz[i]* IQ.pointsz[i])).ToString().Replace(",", "."));
                file2.Write(" ");
            }
            System.IO.StreamWriter file3 = new System.IO.StreamWriter(@"C:\Users\User\Desktop\данные\Contour_y_0.1.txt");
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
            for (int i = 0; i < K; i++)
            {
                tmp = MatrixR[i][i];
                for (int j = K; j >= i; j--)
                    MatrixR[i][j] /= tmp;
                for (int j = i + 1; j < K; j++)
                {
                    tmp = MatrixR[j][i];
                    for (int k = K; k >= i; k--)
                    {
                        MatrixR[j][k] -= tmp * MatrixR[i][k];
                    }
                }
            }
            /*обратный ход*/
            xx[K - 1] = MatrixR[K - 1][K];

            for (int i = K - 2; i >= 0; i--)
            {
                xx[i] = MatrixR[i][K];
                for (int j = i + 1; j < K; j++)
                    xx[i] -= MatrixR[i][j] * xx[j];
                //Console.WriteLine(xx[i]);
            }
            return xx;
        }


        public void InitScreen()
        {

            screenx = new double[numsteps];
            screeny = new double[numsteps];
            screenz = new double[2 * numsteps];
            screenrho = new double[2 * numsteps];

            screenE = new Complex[2 * numsteps];


            //screenz = new double[numsteps];
            //screenrho = new double[numsteps];

            //screenE = new Complex[numsteps];


            double screendim = steprho * numsteps;
            for (int i = 0; i < 2 * numsteps; i++)
            {
                screenz[i] = dist;
                screenrho[i] = 8 * steprho * (i - numsteps);
            }
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
            for (int i = 0; i < numsteps; i++)
            {
                screenx[i] = screeny[i] = (i - numsteps / 2) * steprho;
            }
        }

        public void ScreenE()
        {
            Complex IR = 0, IC = 0;
            for (int i = 0; i < numsteps; i++)
            {
                for (int j = 0; j < K; j++)
                {
                    IC += (new Complex(GreenFuncR(Math.Abs(screenrho[i + numsteps]), screenz[i + numsteps], pointsrho[j], pointsz[j]), GreenFuncC(Math.Abs(screenrho[i + numsteps]), screenz[i + numsteps], pointsrho[j], pointsz[j]))+ new Complex(GreenFuncR(Math.Abs(screenrho[i + numsteps]), screenz[i + numsteps], -pointsrho[j], pointsz[j]), GreenFuncC(Math.Abs(screenrho[i + numsteps]), screenz[i + numsteps], -pointsrho[j], pointsz[j]))) * xx[j];
                }
                IC = Math.Abs(screenrho[i + numsteps]) * coef * IC * ssquare;
                //Console.WriteLine(coef);
                //for (int j = 0; j < K; j++)
                //{
                //    IC += new Complex(GreenFuncR(Math.Abs(screenrho[i]), screenz[i], pointsrho[j], pointsz[j]), GreenFuncC(Math.Abs(screenrho[i]), screenz[i], pointsrho[j], pointsz[j])) * xx[j];
                //}
                //IC = Math.Abs(screenrho[i]) * coef * IC * ssquare;
                //Console.WriteLine(IC);
                //screenE[i] = new Complex(DipoleIntegralRS(i) / (double)(4 * Math.PI), DipoleIntegralCS(i) / (double)(4 * Math.PI)) - IC;

                screenE[numsteps + i] = screenE[numsteps - 1 - i] = -new Complex(DipoleIntegralRS(i + numsteps) / (double)(4 * Math.PI), DipoleIntegralCS(i + numsteps) / (double)(4 * Math.PI)) + IC;
                //Console.WriteLine(Math.Sqrt(screenE[i + numsteps].Real * screenE[i + numsteps].Real + screenE[i + numsteps].Imaginary * screenE[i + numsteps].Imaginary));
                IR = IC = 0;

            }
            System.IO.StreamWriter file5 = new System.IO.StreamWriter(@"C:\Users\User\Desktop\данные\x.txt");
            for (int i = 0; i < 2*numsteps; i++)
            {
                file5.Write(screenrho[i].ToString().Replace(",", "."));
                file5.Write(" ");
            }
            file5.Close();
            System.IO.StreamWriter file4 = new System.IO.StreamWriter(@"C:\Users\User\Desktop\данные\Screen.txt");
            for (int i = 0; i < 2*numsteps; i++)
            {
                file4.Write((Math.Sqrt(screenE[i].Real * screenE[i].Real + screenE[i].Imaginary * screenE[i].Imaginary)).ToString().Replace(",", "."));
                file4.Write(" ");
                //Console.WriteLine((Math.Sqrt(screenE[i].Real * screenE[i].Real + screenE[i].Imaginary * screenE[i].Imaginary)).ToString().Replace(",", "."));
            }
            file4.Close();

        }

        //обратная задача
        public void IterationProcess()
        {
            int dim = 15;
            double steprho_f =  (2*Rad) / dim;
            double stepz_f =  (2*Rad) / dim;
            double rhomid, zmid;
            int n = 0;
            pointsrho = new double[2 * dim * dim];
            pointsz = new double[2 * dim * dim];
            //double[] E_phi_0_R = new double[dim * dim], E_phi_0_C = new double[dim * dim];
            //double[] E_phi_R = new double[dim * dim], E_phi_C = new double[dim * dim];
            double[] K_1 = new double[2 * dim * dim], K_0 = new double[2 * dim * dim], ksi = new double[2 * dim * dim];
            for (int i = -dim; i < dim; i++)
            {
                rhomid = i * steprho_f + steprho_f / 2;
                for (int j = 0; j < dim; j++)
                {
                    zmid = j * stepz_f + stepz_f / 2;
                    pointsrho[n] = rhomid;
                    //zmid = zmid; //+ dist_dip;
                    pointsz[n] = zmid;
                    K_1[n] = k0;
                    K_0[n] = k0;
                    ksi[n] = (K_1[n] - K_0[n]) / (4 * Math.PI);
                    n++;
                }
            }
            K = n;
            System.IO.StreamWriter file3 = new System.IO.StreamWriter(@"C:\Users\User\Desktop\данные\Forward_x.txt");
            for (int i = 0; i < n; i++)
            {
                file3.Write(pointsrho[i].ToString().Replace(",", "."));
                //file2.Write((Math.Sqrt(IQ.pointsrho[i]* IQ.pointsrho[i]- IQ.pointsz[i]* IQ.pointsz[i])).ToString().Replace(",", "."));
                file3.Write(" ");
            }
            System.IO.StreamWriter file4 = new System.IO.StreamWriter(@"C:\Users\User\Desktop\данные\Forward_y.txt");
            file3.Close();
            for (int i = 0; i < n; i++)
            {
                file4.Write(pointsz[i].ToString().Replace(",", "."));
                file4.Write(" ");
            }
            file4.Close();
            numsteps = numsteps * 2;
            MatrixR = new Complex[numsteps][];
            for (int i = 0; i < numsteps; i++)
            {
                MatrixR[i] = new Complex[K + 1];
            }
            Complex[][] MatrixRInter = new Complex[numsteps][];
            for (int i = 0; i < numsteps; i++)
            {
                MatrixRInter[i] = new Complex[K + 1];
            }
            Complex GR = 0, GI = 0;
            double DR = 0, DC = 0;
            ssquare = steprho_f * stepz_f;
            for (int i = 0; i < numsteps; i++)
            {
                for (int j = 0; j < K; j++)
                {
                    GI = screenrho[i]*new Complex(GreenFuncR(screenrho[i], screenz[i], pointsrho[j], pointsz[j]), GreenFuncC(screenrho[i], screenz[i],pointsrho[j], pointsz[j])) * ssquare;
                    DR = -DipoleIntegralR(j) / (double)(4 * Math.PI);
                    DC = -DipoleIntegralC(j) / (double)(4 * Math.PI);
                    MatrixR[i][j] = GI * new Complex(DR, DC);
                    //Console.WriteLine(MatrixR[i][j]);
                }
                MatrixR[i][K] = new Complex(DipoleIntegralRS(i), DipoleIntegralCS(i)) / (double)(4 * Math.PI) + screenE[i];
            }
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
                IMatrix[i][i] += alpha;
            }
            System.IO.StreamWriter file1 = new System.IO.StreamWriter(@"C:\Users\User\Desktop\данные\Forward_Matrix.txt");
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
            System.IO.StreamWriter file2 = new System.IO.StreamWriter(@"C:\Users\User\Desktop\данные\Forward_Right.txt");
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
            System.IO.StreamWriter filec = new System.IO.StreamWriter(@"C:\Users\User\Desktop\данные\Coeff.txt");
            for (int i = 0; i < K; i++)
            {
                filec.Write((coeff[i].Real).ToString().Replace(",", "."));
                filec.Write(" ");
            }
            for (int i = 0; i < K; i++)
            {
                filec.Write((coeff[i].Imaginary).ToString().Replace(",", "."));
                filec.Write(" ");
            }

            filec.WriteLine();
            filec.Close();
        }

        public void IterationProcess2()
        {
            int dim = 15;
            double steprho_f = (2 * Rad) / dim;
            double stepz_f = (2 * Rad) / dim;
            double rhomid, zmid;
            int n = 0;
            pointsrho = new double[2 * dim * dim];
            pointsz = new double[2 * dim * dim];
            //double[] E_phi_0_R = new double[dim * dim], E_phi_0_C = new double[dim * dim];
            //double[] E_phi_R = new double[dim * dim], E_phi_C = new double[dim * dim];
            double[] K_1 = new double[2 * dim * dim], K_0 = new double[2 * dim * dim], ksi = new double[2 * dim * dim];
            for (int i = -dim; i < dim; i++)
            {
                rhomid = i * steprho_f + steprho_f / 2;
                for (int j = 0; j < dim; j++)
                {
                    zmid = j * stepz_f + stepz_f / 2;
                    pointsrho[n] = rhomid;
                    //zmid = zmid; //+ dist_dip;
                    pointsz[n] = zmid;
                    K_1[n] = k0;
                    K_0[n] = k0;
                    ksi[n] = (K_1[n] - K_0[n]) / (4 * Math.PI);
                    n++;
                }
            }
            K = n;
            System.IO.StreamWriter file3 = new System.IO.StreamWriter(@"C:\Users\User\Desktop\данные\Forward_x.txt");
            for (int i = 0; i < n; i++)
            {
                file3.Write(pointsrho[i].ToString().Replace(",", "."));
                //file2.Write((Math.Sqrt(IQ.pointsrho[i]* IQ.pointsrho[i]- IQ.pointsz[i]* IQ.pointsz[i])).ToString().Replace(",", "."));
                file3.Write(" ");
            }
            System.IO.StreamWriter file4 = new System.IO.StreamWriter(@"C:\Users\User\Desktop\данные\Forward_y.txt");
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
            Complex[][] MatrixRInter = new Complex[numsteps][];
            for (int i = 0; i < numsteps; i++)
            {
                MatrixRInter[i] = new Complex[K + 1];
            }
            Complex GR = 0, GI = 0;
            Complex DR = 0, DC = 0;
            ssquare = steprho_f * stepz_f;
            for (int i = 0; i < numsteps; i++)
            {
                for (int j = 0; j < K; j++)
                {
                    GI = screenrho[i] * new Complex(GreenFuncR(screenrho[i], screenz[i], pointsrho[j], pointsz[j]), GreenFuncC(screenrho[i], screenz[i], pointsrho[j], pointsz[j])) * ssquare;
                    DR = xx[j];
                    //DC = DipoleIntegralC(j) / (double)(4 * Math.PI);
                    MatrixR[i][j] = GI * DR;
                    //Console.WriteLine(MatrixR[i][j]);
                }
                MatrixR[i][K] = new Complex(DipoleIntegralRS(i) / (double)(4 * Math.PI), DipoleIntegralCS(i) / (double)(4 * Math.PI)) + screenE[i];
            }
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
                IMatrix[i][i] += alpha;
            }
            System.IO.StreamWriter file1 = new System.IO.StreamWriter(@"C:\Users\User\Desktop\данные\Forward_Matrix.txt");
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
            System.IO.StreamWriter file2 = new System.IO.StreamWriter(@"C:\Users\User\Desktop\данные\Forward_Right.txt");
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

        public void InitMatrix2()
        {
            MatrixR = new Complex[K][];
            for (int i = 0; i < K; i++)
            {
                MatrixR[i] = new Complex[K + 1];
            }
            //coef = (double)(k1 * k1 - k0 * k0) / (double)(4 * Math.PI);
            Complex GR = 0, GI = 0;
            for (int i = 0; i < K; i++)
            {
                for (int j = 0; j < K; j++)
                {

                    GI = xx[j] * new Complex(GreenFuncR(pointsrho[i], pointsz[i], pointsrho[j], pointsz[j]), GreenFuncC(pointsrho[i], pointsz[i], pointsrho[j], pointsz[j])) * ssquare;
                    if (i == j)
                    {
                        MatrixR[i][j] = 1 - pointsrho[i] * GI;
                        //Console.WriteLine(MatrixR[i][j]);
                    }
                    else
                    {
                        MatrixR[i][j] = -pointsrho[i] * GI;
                    }
                }

                MatrixR[i][K] = -new Complex(DipoleIntegralR(i) / (double)(4 * Math.PI), DipoleIntegralC(i) / (double)(4 * Math.PI));
                //Console.WriteLine(MatrixR[i][K]);
            }
            coeff = new Complex[K];
            for (int i = 0; i < K; i++)
            {
                coeff[i] = xx[i];
            }
            System.IO.StreamWriter filec = new System.IO.StreamWriter(@"C:\Users\User\Desktop\данные\Coeff.txt");
            for (int i = 0; i < K; i++)
            {
                filec.Write((coeff[i].Real).ToString().Replace(",", "."));
                filec.Write(" ");
            }
            for (int i = 0; i < K; i++)
            {
                filec.Write((coeff[i].Imaginary).ToString().Replace(",", "."));
                filec.Write(" ");
            }

            filec.WriteLine();
            filec.Close();
            

        }

        public Complex[] GaussI()
        {
            Complex tmp;
            //прямой ход
            Complex[] xx2 = new Complex[K];
            for (int i=0; i<K; i++)
            {
                xx2[i] = coeff[i];
            }
            Console.WriteLine(IMatrix[K - 1][K - 1]);
            for (int i = 0; i < K; i++)
            {
                tmp = IMatrix[i][i];
                for (int j = K; j >= i; j--)
                    IMatrix[i][j] /= tmp;
                for (int j = i + 1; j < K; j++)
                {
                    tmp = IMatrix[j][i];
                    for (int k = K; k >= i; k--)
                    {
                        IMatrix[j][k] -= tmp * IMatrix[i][k];
                    }
                }
            }
            /*обратный ход*/
            xx[K - 1] = IMatrix[K - 1][K];

            for (int i = K - 2; i >= 0; i--)
            {
                xx[i] = IMatrix[i][K];
                for (int j = i + 1; j < K; j++)
                    xx[i] -= IMatrix[i][j] * xx[j];
                Console.WriteLine(xx[i]);
            }
            //for (int i = 0; i < K; i++)
            //{
            //    xx[i] = xx2[i]-xx[i];
            //}
            return xx;

        }

    }
}
