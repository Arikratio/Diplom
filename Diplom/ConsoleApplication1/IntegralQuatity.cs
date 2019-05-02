using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading;
using System.Threading.Tasks;

namespace ConsoleApplication1
{
    class IntegralQuantity
    {
        public double[][] MatrixR;
        public double pi = 3.14;
        public double k1, k0;
        public double Rad=3E-6;
        public double lambda0=600E-9, lambda1=400E-9;
        public int N, M, K;
        public double steprho=0.15E-6, stepz=0.15E-6;
        public double[] pointsrho, pointsz, screenrho, screenx, screeny, screenz;
        public double coef;
        public double omega, mu1=1, eps1=1, mu2 = 1, eps2 = 1, eps0 = 8.85E-12, mu0 = 1.25E-6;
        public double ssquare;
        public double[] xx;
        public double[] jphi;
        public double T;
        public double a=3e-6, b=3e-6;
        public double[] screenE;
        public int numsteps = 20;
        double dist = 7.5E-6;
        public double[][] ScreenMatrix;
        public double alpha = 1e-50;
        public double dist_dip = 1e-6;
        public double[][] IMatrix;


        public double DipoleFuncC(double rho, double z)
        {
            double R = Math.Sqrt(rho * rho + z * z);
            double Rr=rho/(Math.Sqrt(rho * rho + z * z));
            double cos = Math.Cos(k0*R), sin = Math.Sin(k0 * R);
            return T*(Rr * R * k0 * cos - Rr * sin) / (R * R);
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
            return I * Math.Abs(pointsrho[i])*(ssquare);// * pointsrho[j];
        }

        public double DipoleIntegralC(int i)
        {
            double I = 0;
            for (int j = 0; j < K; j++)
            {
                I += GreenFuncC(pointsrho[i], pointsz[i], pointsrho[j], pointsz[j]) * DipoleFuncR(pointsrho[j], pointsz[j]) + GreenFuncR(pointsrho[i], pointsz[i], pointsrho[j], pointsz[j]) * DipoleFuncC(pointsrho[j], pointsz[j]);
            }
            return I * Math.Abs(pointsrho[i])*(ssquare);// * pointsrho[j];
        }

        public double DipoleIntegralRS(int i)
        {
            double I = 0;
            for (int j = 0; j < K; j++)
            {
                I += GreenFuncR(screenrho[i], screenz[i], pointsrho[j], pointsz[j]) * DipoleFuncR(pointsrho[j], pointsz[j]) - GreenFuncC(screenrho[i], screenz[i], pointsrho[j], pointsz[j]) * DipoleFuncC(pointsrho[j], pointsz[j]);
            }
            return I * Math.Abs(screenrho[i]) * (ssquare);// * pointsrho[j];
        }

        public double DipoleIntegralCS(int i)
        {
            double I = 0;
            for (int j = 0; j < K; j++)
            {
                I += GreenFuncC(screenrho[i], screenz[i], pointsrho[j], pointsz[j]) * DipoleFuncR(pointsrho[j], pointsz[j]) + GreenFuncR(screenrho[i], screenz[i], pointsrho[j], pointsz[j]) * DipoleFuncC(pointsrho[j], pointsz[j]);
            }
            return I * Math.Abs(screenrho[i]) * (ssquare);// * pointsrho[j];
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
                ecos = Math.Cos(k1 * R);
                if (R!=0)
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
                esin = Math.Sin(k1 * R);
                if (R != 0)
                    IC = IC + esin / R * cos;
            }
            IC = IC * step;
            return IC;
        }

        public bool IsInCircuit(double rho, double z)
        {
            return (Math.Sqrt(rho * rho + (a-z) * (a-z)) <= a);
        }

        public bool IsInEllipse(double rho, double z)
        {
            return (Math.Sqrt(rho * rho/(a*a) + (b - z) * (b - z)/(b*b)) <= 1);
        }

        public bool IsInCilinder(double rho, double z)
        {
            return ((rho<a)&&(z<2*b));
        }

        public void InitJphi()
        {
            jphi = new double[K];
            for (int i = 0; i < K; i++)
            {
                jphi[i] = 1e-12;
            }
        }

        //public void InitDipole()
        //{
        //    Dipolefunc = new double[K];
        //    for (int i = 0; i < K; i++)
        //    {
        //        Dipolefunc[i] = DipoleFunc(jphi[i]);
        //    }
        //}

        //public void InitMatrixRScreen()
        //{
        //    MatrixR = new double[2 * K][];
        //    for (int i = 0; i < 2 * K; i++)
        //    {
        //        MatrixR[i] = new double[2 * K + 1];
        //    }
        //    coef = (double)(k1 * k1 - k0 * k0) / (double)(4 * Math.PI);
        //    double GR = 0, GI = 0;
        //    for (int i = 0; i < K; i++)
        //    {
        //        for (int j = i; j < K; j++)
        //        {
        //            GR = coef * GreenFuncR(pointsrho[i], pointsz[i], pointsrho[j], pointsz[j]) * ssquare;
        //            GI = coef * GreenFuncC(pointsrho[i], pointsz[i], pointsrho[j], pointsz[j]) * ssquare;
        //            if (i == j)
        //            {
        //                MatrixR[i][j] = 1 + pointsrho[i] * GR;
        //                MatrixR[i + K][j + K] = 1 + pointsrho[i] * GR;
        //            }
        //            else
        //            {
        //                MatrixR[i][j] = pointsrho[i] * GR;
        //                MatrixR[i + K][j + K] = pointsrho[i] * GR;
        //                MatrixR[j][i] = pointsrho[j] * GR;
        //                MatrixR[j + K][i + K] = pointsrho[j] * GR;
        //            }
        //            MatrixR[i][j + K] = -pointsrho[i] * GI;
        //            MatrixR[i + K][j] = pointsrho[i] * GI;
        //            MatrixR[j][i + K] = -pointsrho[j] * GI;
        //            MatrixR[j + K][i] = pointsrho[j] * GI;

        //        }

        //        MatrixR[i][2 * K] = DipoleIntegralR(i) / (double)(4 * Math.PI);
        //        MatrixR[i + K][2 * K] = DipoleIntegralC(i) / (double)(4 * Math.PI);
        //        //Console.WriteLine(i);
        //    }
        //    //Console.WriteLine(MatrixR[0][2 * K]);
        //    System.IO.StreamWriter file1 = new System.IO.StreamWriter(@"C:\Users\tsybr\Desktop\данные\Matrix_Contour_0.1.txt");
        //    for (int i = 0; i < 2 * K; i++)
        //    {
        //        for (int j = 0; j < 2 * K; j++)
        //        {
        //            file1.Write((MatrixR[i][j]).ToString().Replace(",", "."));
        //            file1.Write(" ");
        //        }
        //        file1.Write("\n");
        //    }
        //    file1.WriteLine();
        //    file1.Close();
        //    System.IO.StreamWriter file2 = new System.IO.StreamWriter(@"C:\Users\tsybr\Desktop\данные\Right_Contour_0.1.txt");
        //    for (int i = 0; i < 2 * K; i++)
        //    {
        //        file2.Write((MatrixR[i][2 * K]).ToString().Replace(",", "."));
        //        file2.Write(" ");
        //    }

        //    file2.WriteLine();
        //    file2.Close();
        //    System.IO.StreamWriter file3 = new System.IO.StreamWriter(@"C:\Users\tsybr\Desktop\данные\DipoleR_Contour_0.1.txt");
        //    for (int i = 0; i < K; i++)
        //    {
        //        //file3.Write((Math.Sqrt(GreenFuncR(pointsrho[0], pointsz[0], pointsrho[i], pointsz[i]) * GreenFuncR(pointsrho[0], pointsz[0], pointsrho[i], pointsz[i]) + GreenFuncC(pointsrho[0], pointsz[0], pointsrho[i], pointsz[i]) * GreenFuncC(pointsrho[0], pointsz[0], pointsrho[i], pointsz[i]))).ToString().Replace(",", "."));
        //        file3.Write((Math.Sqrt(DipoleFuncR(pointsrho[i], pointsz[i]) * DipoleFuncR(pointsrho[i], pointsz[i]) + DipoleFuncC(pointsrho[i], pointsz[i]) * DipoleFuncC(pointsrho[i], pointsz[i]))).ToString().Replace(",", "."));
        //        file3.Write(" ");
        //    }

        //    file3.WriteLine();
        //    file3.Close();

        //}

        public void InitMatrixR()
        {
            MatrixR = new double[2*K][];
            for (int i = 0; i < 2*K; i++)
            {
                MatrixR[i] = new double[2*K+1];
            }
            coef =(double)(k1 * k1 - k0 * k0) / (double)(4 * Math.PI);
            double GR=0, GI=0;
            for (int i = 0; i < K; i++)
            {
                for (int j = i; j < K; j++)
                {
                    GR= coef * GreenFuncR(pointsrho[i], pointsz[i], pointsrho[j], pointsz[j]) * ssquare;
                    GI= coef * GreenFuncC(pointsrho[i], pointsz[i], pointsrho[j], pointsz[j]) * ssquare;
                    if (i == j) 
                    {
                        MatrixR[i][j] = 1 + pointsrho[i]*GR;
                        MatrixR[i+K][j+K] = 1 + pointsrho[i]*GR;
                    }
                    else
                    {
                        MatrixR[i][j] = pointsrho[i]*GR;
                        MatrixR[i + K][j + K] = pointsrho[i]*GR;
                        MatrixR[j][i] = pointsrho[j] * GR;
                        MatrixR[j + K][i + K] = pointsrho[j] * GR;
                    }
                    MatrixR[i][j+K] = -pointsrho[i]*GI;
                    MatrixR[i+K][j] = pointsrho[i]*GI;
                    MatrixR[j][i + K] = -pointsrho[j] * GI;
                    MatrixR[j + K][i] = pointsrho[j] * GI;
                    
                }
               
                MatrixR[i][2 * K] = DipoleIntegralR(i) / (double)(4 * Math.PI);
                MatrixR[i+K][2 * K] = DipoleIntegralC(i) / (double)(4 * Math.PI);
                //Console.WriteLine(i);
            }
            //Console.WriteLine(MatrixR[0][2 * K]);
            System.IO.StreamWriter file1 = new System.IO.StreamWriter(@"C:\Users\tsybr\Desktop\данные\Matrix_Contour_0.1.txt");
            for (int i = 0; i < 2*K; i++)
            {
                for (int j = 0; j < 2 * K; j++)
                {
                    file1.Write((MatrixR[i][j]).ToString().Replace(",", "."));
                    file1.Write(" ");
                }
                file1.Write("\n");
            }
            file1.WriteLine();
            file1.Close();
            System.IO.StreamWriter file2 = new System.IO.StreamWriter(@"C:\Users\tsybr\Desktop\данные\Right_Contour_0.1.txt");
                for (int i = 0; i < 2 * K; i++)
                {
                    file2.Write((MatrixR[i][2*K]).ToString().Replace(",", "."));
                    file2.Write(" ");
                }

            file2.WriteLine();
            file2.Close();
            System.IO.StreamWriter file3 = new System.IO.StreamWriter(@"C:\Users\tsybr\Desktop\данные\DipoleR_Contour_0.1.txt");
            for (int i = 0; i < K; i++)
            {
                //file3.Write((Math.Sqrt(GreenFuncR(pointsrho[0], pointsz[0], pointsrho[i], pointsz[i]) * GreenFuncR(pointsrho[0], pointsz[0], pointsrho[i], pointsz[i]) + GreenFuncC(pointsrho[0], pointsz[0], pointsrho[i], pointsz[i]) * GreenFuncC(pointsrho[0], pointsz[0], pointsrho[i], pointsz[i]))).ToString().Replace(",", "."));
                file3.Write((Math.Sqrt(DipoleFuncR(pointsrho[i], pointsz[i])* DipoleFuncR(pointsrho[i], pointsz[i])+ DipoleFuncC(pointsrho[i], pointsz[i])* DipoleFuncC(pointsrho[i], pointsz[i]))).ToString().Replace(",", "."));
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
                    if (IsInCircuit(rhomid, zmid))
                        n++;
                    //if (IsInEllipse(rhomid, zmid))
                    ////    n++;
                    //if (IsInCilinder(rhomid, zmid))
                    //    n++;
                }
            K = n;
            xx = new double[2*n];
            pointsrho = new double[n];
            pointsz = new double[n];
            n = 0;
            for (int i = 0; i < M; i++)
            {
                rhomid = i * steprho + steprho / 2;
                for (int j = 0; j < N; j++)
                {

                    zmid = j * stepz + stepz / 2;
                    if (IsInCircuit(rhomid, zmid))
                    {
                        pointsrho[n] = rhomid;
                        zmid = zmid + dist_dip;
                        //pointsrho[n] = Math.Sqrt(rhomid * rhomid + zmid * zmid);
                        pointsz[n] = zmid;
                        Console.WriteLine(pointsz[n]);
                        n++;
                    }
                    //if (IsInEllipse(rhomid, zmid))
                    //{
                    //    pointsrho[n] = rhomid;
                    //    zmid = zmid + 1e-6;
                    //    //pointsrho[n] = Math.Sqrt(rhomid * rhomid + zmid * zmid);
                    //    pointsz[n] = zmid;
                    //    n++;
                    //}
                    //if (IsInCilinder(rhomid, zmid))
                    //{
                    //    pointsrho[n] = rhomid;
                    //    zmid = zmid + dist_dip;
                    //    //pointsrho[n] = Math.Sqrt(rhomid * rhomid + zmid * zmid);
                    //    pointsz[n] = zmid;
                    //    n++;
                    //}
                }
            }
            System.IO.StreamWriter file2 = new System.IO.StreamWriter(@"C:\Users\tsybr\Desktop\данные\Contour_x_0.1.txt");
            for (int i = 0; i < K; i++)
            {
                file2.Write(pointsrho[i].ToString().Replace(",", "."));
                //file2.Write((Math.Sqrt(IQ.pointsrho[i]* IQ.pointsrho[i]- IQ.pointsz[i]* IQ.pointsz[i])).ToString().Replace(",", "."));
                file2.Write(" ");
            }
            System.IO.StreamWriter file3 = new System.IO.StreamWriter(@"C:\Users\tsybr\Desktop\данные\Contour_y_0.1.txt");
            file2.Close();
            for (int i = 0; i < K; i++)
            {
                file3.Write(pointsz[i].ToString().Replace(",", "."));
                file3.Write(" ");
            }
            file3.Close();

        }

        public double[] Gauss()
        {
            double tmp;
            //прямой ход
            Console.WriteLine(MatrixR[2 * K - 1][2 * K-1]);
            for (int i = 0; i < 2*K; i++)
            {
                tmp = MatrixR[i][i];
                for (int j = 2*K; j >= i; j--)
                    MatrixR[i][j] /= tmp;
                for (int j = i + 1; j < 2*K; j++)
                {
                    tmp = MatrixR[j][i];
                    for (int k = 2*K; k >= i; k--)
                    {
                        MatrixR[j][k] -= tmp * MatrixR[i][k];
                    }
                }
            }
            /*обратный ход*/
            xx[2*K - 1] = MatrixR[2*K-1][2*K];
            
            for (int i = 2*K - 2; i >= 0; i--)
            {
                xx[i] = MatrixR[i][2*K];
                for (int j = i + 1; j < 2*K; j++)
                    xx[i] -= MatrixR[i][j] * xx[j];
            }
            return xx;
        }

        public void InitScreen()
        {
           
            screenx = new double[numsteps];
            screeny = new double[numsteps];
            screenz = new double[2 * numsteps];
            screenrho = new double[2 * numsteps];

            screenE = new double[4*numsteps];
            ScreenMatrix= new double[2*numsteps][];
            for (int i = 0; i < 2*numsteps; i++)
            {
                ScreenMatrix[i] = new double[numsteps];
            }
            double screendim = steprho * numsteps;
            for (int i=0; i<2*numsteps; i++)
            {
                screenz[i] = dist;
                screenrho[i] = 2*steprho * (i-numsteps);
            }
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
                screenx[i] = screeny[i] = (i - numsteps/2) * steprho ;
            }
        }

        public void ScreenE()
        {
            double IR = 0, IC = 0;
            for (int i=0; i<numsteps; i++)
            {
                for (int j=0; j<K; j++)
                {
                    IR += GreenFuncR(Math.Abs(screenrho[i+ numsteps]), screenz[i + numsteps], pointsrho[j], pointsz[j]) * xx[j] - GreenFuncC(screenrho[i + numsteps], screenz[i + numsteps], pointsrho[j], pointsz[j]) * xx[j + K];
                    IC += GreenFuncC(Math.Abs(screenrho[i + numsteps]), screenz[i + numsteps], pointsrho[j], pointsz[j]) * xx[j] + GreenFuncR(screenrho[i + numsteps], screenz[i + numsteps], pointsrho[j], pointsz[j]) * xx[j + K];
                }
                IR = Math.Abs(screenrho[i + numsteps])*coef * IR * ssquare;
                IC = Math.Abs(screenrho[i + numsteps])*coef * IC * ssquare;
                screenE[numsteps+i] = screenE[numsteps-1 - i] =DipoleIntegralRS(i) / (double)(4 * Math.PI) - IR;;
                screenE[i+3*numsteps] = screenE[-i -1 + 3 * numsteps]= DipoleIntegralCS(i) / (double)(4 * Math.PI) - IC;
                Console.WriteLine(Math.Sqrt(screenE[i+numsteps] * screenE[i + numsteps] + screenE[i + numsteps + 2 * numsteps] * screenE[i + numsteps + 2 * numsteps]));
                IR = IC = 0;

            }
            //for (int i = 2*numsteps; i < 3*numsteps; i++)
            //{
            //    for (int j = 0; j < K; j++)
            //    {
            //        IR += GreenFuncR(screenrho[i], screenz[i], pointsrho[j], pointsz[j]) * xx[j] - GreenFuncC(screenrho[i], screenz[i], pointsrho[j], pointsz[j]) * xx[j + K];
            //        IC += GreenFuncC(screenrho[i], screenz[i], pointsrho[j], pointsz[j]) * xx[j] + GreenFuncR(screenrho[i], screenz[i], pointsrho[j], pointsz[j]) * xx[j + K];
            //    }
            //    IR = Math.Abs(screenrho[i])*coef * IR * ssquare;
            //    IC = Math.Abs(screenrho[i])*coef * IC * ssquare;
            //    screenE[i] = screenE[i+numsteps] = DipoleIntegralRS(i) / (double)(4 * Math.PI) - IR;
            //    screenE[i+4*numsteps] = screenE[i + numsteps + 4 * numsteps] = DipoleIntegralCS(i) / (double)(4 * Math.PI) - IC;
            //    Console.WriteLine(Math.Sqrt(screenE[i] * screenE[i] + screenE[i + 4*numsteps] * screenE[i + 4*numsteps]));
            //    IR = IC = 0;

            //}
            System.IO.StreamWriter file5 = new System.IO.StreamWriter(@"C:\Users\tsybr\Desktop\данные\x.txt");
            for (int i = 0; i < 2*numsteps; i++)
            {
                file5.Write(screenrho[i].ToString().Replace(",", "."));
                file5.Write(" ");
            }
            file5.Close();
            System.IO.StreamWriter file4 = new System.IO.StreamWriter(@"C:\Users\tsybr\Desktop\данные\Screen.txt");
            for (int i = 0; i < 2* numsteps; i++)
            {
                file4.Write((Math.Sqrt(screenE[i] * screenE[i] + screenE[i + 2*numsteps] * screenE[i + 2*numsteps])).ToString().Replace(",", "."));
                file4.Write(" ");
            }
            file4.Close();
            
        }

        //public void ScreenEMatrix()
        //{
        //    double IR = 0, IC = 0;
        //    for (int i = 0; i < numsteps; i++)
        //    {
        //        for (int j = 0; j < numsteps; j++)
        //        {
        //            screenrho[i] = Math.Sqrt(screenx[i] * screenx[i] + screeny[j] * screeny[j]);
        //            for (int k = 0; k < K; k++)
        //            {
        //                IR += GreenFuncR(screenrho[i], screenz[i], pointsrho[k], pointsz[k]) * xx[k] - GreenFuncC(screenrho[i], screenz[i], pointsrho[k], pointsz[k]) * xx[k + K];
        //                IC += GreenFuncC(screenrho[i], screenz[i], pointsrho[k], pointsz[k]) * xx[k] + GreenFuncR(screenrho[i], screenz[i], pointsrho[k], pointsz[k]) * xx[k + K];
        //            }
        //            IR = coef * IR * ssquare;
        //            IC = coef * IC * ssquare;
        //            ScreenMatrix[i][j] = DipoleIntegralRS(i) / (double)(4 * Math.PI) - IR;
        //            ScreenMatrix[i + numsteps][j] = DipoleIntegralCS(i) / (double)(4 * Math.PI) - IC;
        //            Console.WriteLine(Math.Sqrt(ScreenMatrix[i][j] * ScreenMatrix[i][j] + ScreenMatrix[i + numsteps][j] * ScreenMatrix[i + numsteps][j]));
        //            IR = IC = 0;
        //        }
                
        //    }

        //}
        //обратная задача
        public void IterationProcess()
        {
            int dim = 30;
            double steprho_f = (2 * Rad + 3e-6) / dim;
            double stepz_f = (2 * Rad+3e-6) / dim;
            double rhomid, zmid;
            int n = 0;
            pointsrho = new double[2*dim * dim];
            pointsz = new double[2*dim * dim];
            //double[] E_phi_0_R = new double[dim * dim], E_phi_0_C = new double[dim * dim];
            //double[] E_phi_R = new double[dim * dim], E_phi_C = new double[dim * dim];
            double[] K_1 = new double[2*dim * dim], K_0 = new double[2*dim * dim], ksi= new double[2*dim * dim];
            for (int i = -dim; i < dim; i++)
            {
                rhomid = i * steprho_f + steprho_f / 2;
                for (int j = 0; j < dim; j++)
                {
                    zmid = j * stepz_f + stepz_f / 2;
                    pointsrho[n] = rhomid;
                    zmid = zmid; //+ dist_dip;
                    pointsz[n] = zmid;
                    K_1[n] = k0;
                    K_0[n] = k0;
                    ksi[n] = (K_1[n] - K_0[n]) / (4*Math.PI);
                    n++;
                }
            }
            K = n;
            System.IO.StreamWriter file3 = new System.IO.StreamWriter(@"C:\Users\tsybr\Desktop\данные\Forward_x.txt");
            for (int i = 0; i < n; i++)
            {
                file3.Write(pointsrho[i].ToString().Replace(",", "."));
                //file2.Write((Math.Sqrt(IQ.pointsrho[i]* IQ.pointsrho[i]- IQ.pointsz[i]* IQ.pointsz[i])).ToString().Replace(",", "."));
                file3.Write(" ");
            }
            System.IO.StreamWriter file4 = new System.IO.StreamWriter(@"C:\Users\tsybr\Desktop\данные\Forward_y.txt");
            file3.Close();
            for (int i = 0; i < n; i++)
            {
                file4.Write(pointsz[i].ToString().Replace(",", "."));
                file4.Write(" ");
            }
            file4.Close();
            numsteps = numsteps * 2;
            MatrixR = new double[2 * numsteps][];
            for (int i = 0; i < 2 * numsteps; i++)
            {
                MatrixR[i] = new double[2 * K + 1];
            }
            double[][] MatrixRInter = new double[2 * numsteps][];
            for (int i = 0; i < 2 * numsteps; i++)
            {
                MatrixRInter[i] = new double[2 * K + 1];
            }
            double GR = 0, GI = 0;
            double DR = 0, DC = 0;
            ssquare = steprho_f * stepz_f;
            for (int i = 0; i < numsteps; i++)
            {
                for (int j = 0; j < K; j++)
                {
                    GR = GreenFuncR(Math.Abs(screenrho[i]), screenz[i], Math.Abs(pointsrho[j]), pointsz[j]) * ssquare;
                    GI = GreenFuncC(Math.Abs(screenrho[i]), screenz[i], Math.Abs(pointsrho[j]), pointsz[j]) * ssquare;
                    DR = DipoleIntegralR(j) / (double)(4 * Math.PI);
                    DC = DipoleIntegralC(j) / (double)(4 * Math.PI);
                    MatrixR[i][j] = Math.Abs(screenrho[i]) * (GR * DR-GI*DC) ;
                    MatrixR[i + numsteps][j + K] = Math.Abs(screenrho[i]) * (GR * DR - GI * DC);
                    MatrixR[i][j + K] = -Math.Abs(screenrho[i]) * (GR * DC+GI*DC);
                    MatrixR[i + numsteps][j] = Math.Abs(screenrho[i]) * (GR * DC + GI * DC);
                    MatrixRInter[i][j] = Math.Abs(screenrho[i]) * (GR * DR - GI * DC);
                    MatrixRInter[i + numsteps][j + K] = Math.Abs(screenrho[i]) * (GR * DR - GI * DC);
                    MatrixRInter[i][j + K] = Math.Abs(screenrho[i]) * (GR * DC + GI * DR);
                    MatrixRInter[i + numsteps][j] = -Math.Abs(screenrho[i]) * (GR * DC + GI * DR);
                }

                MatrixR[i][2 * K] = DipoleIntegralRS(i)-screenE[i];
                MatrixR[i + numsteps][2 * K] = DipoleIntegralCS(i) - screenE[i+numsteps];
            }
            IMatrix = new double[2 * K][];
            for (int i = 0; i < 2 * K; i++)
            {
                IMatrix[i] = new double[2 * K + 1];
            }
            for (int i = 0; i < 2 * K; i++)
            {
                for (int j = 0; j < 2 * K; j++)
                {
                    IMatrix[i][j] = 0;
                    for (int k = 0; k < 2 * numsteps; k++)
                    {
                        IMatrix[i][j] += MatrixRInter[k][i] * MatrixR[k][j];
                        
                    }
                }
            }
            
            for (int i = 0; i < 2 * K; i++)
            {
                IMatrix[i][2 * K] = 0;
                for (int k = 0; k < 2 * numsteps; k++)
                {
                    IMatrix[i][2 * K] += MatrixRInter[k][i] * MatrixR[k][2 * K];
                }
            }
            for (int i = 0; i < 2 * K; i++)
            {
                IMatrix[i][i] += alpha;
            }
            System.IO.StreamWriter file1 = new System.IO.StreamWriter(@"C:\Users\tsybr\Desktop\данные\Forward_Matrix.txt");
            for (int i = 0; i < 2 * K; i++)
            {
                for (int j = 0; j < 2 * K; j++)
                {
                    file1.Write((IMatrix[i][j]).ToString().Replace(",", "."));
                    file1.Write(" ");
                }
                file1.Write("\n");
            }
            file1.WriteLine();
            file1.Close();
            System.IO.StreamWriter file2 = new System.IO.StreamWriter(@"C:\Users\tsybr\Desktop\данные\Forward_Right.txt");
            for (int i = 0; i < 2 * K; i++)
            {
                file2.Write((IMatrix[i][2 * K]).ToString().Replace(",", "."));
                file2.Write(" ");
            }

            file2.WriteLine();
            file2.Close();
            xx = new double[2 * K];
        }

        public double[] GaussI()
        {
            double tmp;
            //прямой ход
            Console.WriteLine(IMatrix[2 * K - 1][2 * K - 1]);
            for (int i = 0; i < 2 * K; i++)
            {
                tmp = IMatrix[i][i];
                for (int j = 2 * K; j >= i; j--)
                    IMatrix[i][j] /= tmp;
                for (int j = i + 1; j < 2 * K; j++)
                {
                    tmp = IMatrix[j][i];
                    for (int k = 2 * K; k >= i; k--)
                    {
                        IMatrix[j][k] -= tmp * IMatrix[i][k];
                    }
                }
            }
            /*обратный ход*/
            xx[2 * K - 1] = IMatrix[2 * K - 1][2 * K];

            for (int i = 2 * K - 2; i >= 0; i--)
            {
                xx[i] = IMatrix[i][2 * K];
                for (int j = i + 1; j < 2 * K; j++)
                    xx[i] -= IMatrix[i][j] * xx[j];
            }
            
            for (int i=0; i<K; i++)
            {

                Console.WriteLine(-Math.Sqrt(xx[i] * xx[i] + xx[i + K] * xx[i + K])*4*Math.PI+k0*k0);
            }
            Console.WriteLine("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
            Console.WriteLine(k0);
            return xx;
        }

        //обратная задача
        public void IterationProcess2()
        {
            int dim = 30;
            double steprho_f = 2 * Rad / dim;
            double stepz_f = (2 * Rad + 1e-6) / dim;
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
                    zmid = zmid; //+ dist_dip;
                    pointsz[n] = zmid;
                    K_1[n] = k0;
                    K_0[n] = k0;
                    ksi[n] = (K_1[n] - K_0[n]) / (4 * Math.PI);
                    n++;
                }
            }
            K = n;
            System.IO.StreamWriter file3 = new System.IO.StreamWriter(@"C:\Users\tsybr\Desktop\данные\Forward_x.txt");
            for (int i = 0; i < n; i++)
            {
                file3.Write(pointsrho[i].ToString().Replace(",", "."));
                //file2.Write((Math.Sqrt(IQ.pointsrho[i]* IQ.pointsrho[i]- IQ.pointsz[i]* IQ.pointsz[i])).ToString().Replace(",", "."));
                file3.Write(" ");
            }
            System.IO.StreamWriter file4 = new System.IO.StreamWriter(@"C:\Users\tsybr\Desktop\данные\Forward_y.txt");
            file3.Close();
            for (int i = 0; i < n; i++)
            {
                file4.Write(pointsz[i].ToString().Replace(",", "."));
                file4.Write(" ");
            }
            file4.Close();
            MatrixR = new double[2 * numsteps][];
            for (int i = 0; i < 2 * numsteps; i++)
            {
                MatrixR[i] = new double[2 * K + 1];
            }
            double[][] MatrixRInter = new double[2 * numsteps][];
            for (int i = 0; i < 2 * numsteps; i++)
            {
                MatrixRInter[i] = new double[2 * K + 1];
            }
            double GR = 0, GI = 0;
            double DR = 0, DC = 0;
            ssquare = steprho_f * stepz_f;
            for (int i = 0; i < numsteps; i++)
            {
                for (int j = 0; j < K; j++)
                {
                    GR = GreenFuncR(Math.Abs(screenrho[i]), screenz[i], Math.Abs(pointsrho[j]), pointsz[j]) * ssquare;
                    GI = GreenFuncC(Math.Abs(screenrho[i]), screenz[i], Math.Abs(pointsrho[j]), pointsz[j]) * ssquare;
                    DR = xx[j];
                    DC = xx[j+K];
                    MatrixR[i][j] = Math.Abs(screenrho[i]) * (GR * DR - GI * DC);
                    MatrixR[i + numsteps][j + K] = Math.Abs(screenrho[i]) * (GR * DR - GI * DC);
                    MatrixR[i][j + K] = -Math.Abs(screenrho[i]) * (GR * DC + GI * DR);
                    MatrixR[i + numsteps][j] = Math.Abs(screenrho[i]) * (GR * DC + GI * DR);
                    MatrixRInter[i][j] = Math.Abs(screenrho[i]) * (GR * DR - GI * DC);
                    MatrixRInter[i + numsteps][j + K] = Math.Abs(screenrho[i]) * (GR * DR - GI * DC);
                    MatrixRInter[i][j + K] = Math.Abs(screenrho[i]) * (GR * DC + GI * DR);
                    MatrixRInter[i + numsteps][j] = -Math.Abs(screenrho[i]) * (GR * DC + GI * DR);
                    //Console.WriteLine(Math.Sqrt(MatrixR[i][j] * MatrixR[i][j] + MatrixR[i][j + K] * MatrixR[i][j + K]));
                }
                //Console.WriteLine("&&&&&&&&&");
                MatrixR[i][2 * K] = DipoleIntegralRS(i) - screenE[i];
                MatrixR[i + numsteps][2 * K] = DipoleIntegralCS(i) - screenE[i + numsteps];
                //Console.WriteLine(Math.Sqrt(MatrixR[i][2 * K] * MatrixR[i][2 * K] + MatrixR[i + numsteps][2 * K] * MatrixR[i + numsteps][2 * K]));
            }
            IMatrix = new double[2 * K][];
            for (int i = 0; i < 2 * K; i++)
            {
                IMatrix[i] = new double[2 * K + 1];
            }
            for (int i = 0; i < 2 * K; i++)
            {
                for (int j = 0; j < 2 * K; j++)
                {
                    IMatrix[i][j] = 0;
                    for (int k = 0; k < 2 * numsteps; k++)
                    {
                        IMatrix[i][j] += MatrixRInter[k][i] * MatrixR[k][j];
                    }
                }
            }

            for (int i = 0; i < 2 * K; i++)
            {
                IMatrix[i][2 * K] = 0;
                for (int k = 0; k < 2 * numsteps; k++)
                {
                    IMatrix[i][2 * K] += MatrixRInter[k][i] * MatrixR[k][2 * K];
                }
            }
            for (int i = 0; i < 2 * K; i++)
            {
                IMatrix[i][i] += alpha;
            }
            System.IO.StreamWriter file1 = new System.IO.StreamWriter(@"C:\Users\tsybr\Desktop\данные\Forward_Matrix.txt");
            for (int i = 0; i < 2 * K; i++)
            {
                for (int j = 0; j < 2 * K; j++)
                {
                    file1.Write((IMatrix[i][j]).ToString().Replace(",", "."));
                    file1.Write(" ");
                }
                file1.Write("\n");
            }
            file1.WriteLine();
            file1.Close();
            System.IO.StreamWriter file2 = new System.IO.StreamWriter(@"C:\Users\tsybr\Desktop\данные\Forward_Right.txt");
            for (int i = 0; i < 2 * K; i++)
            {
                file2.Write((IMatrix[i][2 * K]).ToString().Replace(",", "."));
                file2.Write(" ");
            }

            file2.WriteLine();
            file2.Close();
            xx = new double[2 * K];
        }

        public void InitMatrixR2()
        {
            MatrixR = new double[2 * K][];
            for (int i = 0; i < 2 * K; i++)
            {
                MatrixR[i] = new double[2 * K + 1];
            }
            //coef = (double)(k1 * k1 - k0 * k0) / (double)(4 * Math.PI);
            double GR = 0, GI = 0;
            for (int i = 0; i < K; i++)
            {
                for (int j = 0; j < K; j++)
                {
                    GR = (double)(Math.Sqrt(xx[j] * xx[j] + xx[j + K] * xx[j + K]))* GreenFuncR(pointsrho[i], pointsz[i], pointsrho[j], pointsz[j]) * ssquare;
                    GI = (double)(Math.Sqrt(xx[j] * xx[j] + xx[j + K] * xx[j + K]))* GreenFuncC(pointsrho[i], pointsz[i], pointsrho[j], pointsz[j]) * ssquare;
                    if (i == j)
                    {
                        MatrixR[i][j] = 1 + Math.Abs(pointsrho[i]) * GR;
                        MatrixR[i + K][j + K] = 1 + Math.Abs(pointsrho[i]) * GR;
                    }
                    else
                    {
                        MatrixR[i][j] = Math.Abs(pointsrho[i]) * GR;
                        MatrixR[i + K][j + K] = Math.Abs(pointsrho[i]) * GR;
                        //MatrixR[j][i] = Math.Abs(pointsrho[j]) * GR;
                        //MatrixR[j + K][i + K] = Math.Abs(pointsrho[j]) * GR;
                    }
                    MatrixR[i][j + K] = -Math.Abs(pointsrho[i]) * GI;
                    MatrixR[i + K][j] = Math.Abs(pointsrho[i]) * GI;
                    //MatrixR[j][i + K] = -Math.Abs(pointsrho[j]) * GI;
                    //MatrixR[j + K][i] = Math.Abs(pointsrho[j]) * GI;

                }

                MatrixR[i][2 * K] = DipoleIntegralR(i) / (double)(4 * Math.PI);
                MatrixR[i + K][2 * K] = DipoleIntegralC(i) / (double)(4 * Math.PI);
                //Console.WriteLine(i);
            }
            //Console.WriteLine(MatrixR[0][2 * K]);
            System.IO.StreamWriter file1 = new System.IO.StreamWriter(@"C:\Users\tsybr\Desktop\данные\Matrix_Contour_0.1.txt");
            for (int i = 0; i < 2 * K; i++)
            {
                for (int j = 0; j < 2 * K; j++)
                {
                    file1.Write((MatrixR[i][j]).ToString().Replace(",", "."));
                    file1.Write(" ");
                }
                file1.Write("\n");
            }
            file1.WriteLine();
            file1.Close();
            System.IO.StreamWriter file2 = new System.IO.StreamWriter(@"C:\Users\tsybr\Desktop\данные\Right_Contour_0.1.txt");
            for (int i = 0; i < 2 * K; i++)
            {
                file2.Write((MatrixR[i][2 * K]).ToString().Replace(",", "."));
                file2.Write(" ");
            }

            file2.WriteLine();
            file2.Close();
            System.IO.StreamWriter file3 = new System.IO.StreamWriter(@"C:\Users\tsybr\Desktop\данные\DipoleR_Contour_0.1.txt");
            for (int i = 0; i < K; i++)
            {
                //file3.Write((Math.Sqrt(GreenFuncR(pointsrho[0], pointsz[0], pointsrho[i], pointsz[i]) * GreenFuncR(pointsrho[0], pointsz[0], pointsrho[i], pointsz[i]) + GreenFuncC(pointsrho[0], pointsz[0], pointsrho[i], pointsz[i]) * GreenFuncC(pointsrho[0], pointsz[0], pointsrho[i], pointsz[i]))).ToString().Replace(",", "."));
                file3.Write((Math.Sqrt(DipoleFuncR(pointsrho[i], pointsz[i]) * DipoleFuncR(pointsrho[i], pointsz[i]) + DipoleFuncC(pointsrho[i], pointsz[i]) * DipoleFuncC(pointsrho[i], pointsz[i]))).ToString().Replace(",", "."));
                file3.Write(" ");
            }

            file3.WriteLine();
            file3.Close();



        }
    }
}
