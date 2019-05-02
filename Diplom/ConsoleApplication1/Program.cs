using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Numerics;
using System.Text;
using System.Threading.Tasks;

namespace ConsoleApplication1
{
    class Program
    {
        static IntegralQuantity IQ = new IntegralQuantity();
        static QTS IQS = new QTS();
        static QTS_parallel IQS_parallel = new QTS_parallel();
        static void Main(string[] args)
        {
            Process.GetCurrentProcess().PriorityClass = ProcessPriorityClass.Idle;
            var stopwatch = Stopwatch.StartNew();
            //IQS.CountNumberMas();
            //IQS.InitJphi();
            //IQS.InitMatrix();
            //IQS.Gauss();
            //IQS.InitScreen();
            //IQS.ScreenE();
            //IQS.IterationProcess();
            //stopwatch.Stop();
            //Console.WriteLine("не параллельно  " + stopwatch.Elapsed);
            //stopwatch.Restart();

            //Random rand = new Random();
            //Complex[][] C;
            //C = new Complex[3][];
            //for (int i = 0; i <3; i++)
            //{
            //    C[i] = new Complex[3];
            //}
            //for (int i = 0; i < 3; i++)
            //{
            //    for (int j = 0; j < 3; j++)
            //    {
            //        C[i][j] = rand.Next(0, 10);
            //        Console.Write("{0} ",C[i][j]);
            //    }
            //    Console.WriteLine();
            //}
            //Console.WriteLine(C.GetLength(0));
            //Console.WriteLine(IQS_parallel.MatrixDet(C));
            //Console.WriteLine(IQS_parallel.ReverseMatrix(C));

            Console.BufferHeight = 30000;
            IQS_parallel.CountNumberMas();
            IQS_parallel.InitJphi();
            IQS_parallel.InitMatrix();
            IQS_parallel.MatlabGauss();
            IQS_parallel.InitScreen();
            IQS_parallel.ScreenE();
            //IQS_parallel.InitIterProc(50);

            IQS_parallel.IterationProcessKOLD(0.1);
            stopwatch.Stop();
            Console.WriteLine("параллельно  " + stopwatch.Elapsed);
            

            for (int i = 0; i < 5; i++)
            {
                //IQS_parallel.GaussI();
                IQS_parallel.MatlabGaussI();
                IQS_parallel.InitMatrix2();
                IQS_parallel.MatlabGauss();
                IQS_parallel.IterationProcess2();
            }
            IQS_parallel.MatlabGaussI();
            //IQS_parallel.ZeidelI();
            //IQS_parallel.GaussI();
            //IQ.CountNumberMas();
            //Console.WriteLine(IQ.K);
            //IQ.InitJphi();
            //IQ.InitMatrixR();
            //IQ.Gauss();
            //IQ.InitScreen();
            //IQ.ScreenE();
            //IQ.IterationProcess();
            //for (int i = 0; i < 10; i++)
            //{
            //    IQ.GaussI();
            //    IQ.InitMatrixR2();
            //    IQ.Gauss();
            //    IQ.IterationProcess2();
            //}

            //System.IO.StreamWriter file1 = new System.IO.StreamWriter(@"J:\Diplom\Diplom\ConsoleApplication1\WriteLines1.txt");
            //for (int i=0; i<IQ.K; i++)
            //{
            //    file1.Write(Math.Sqrt(IQ.xx[i]* IQ.xx[i] + IQ.xx[i+IQ.K] * IQ.xx[i + IQ.K]).ToString().Replace(",", "."));
            //    file1.Write(" ");
            //}
            ////for (int i = IQ.K; i < 2*IQ.K; i++)
            ////{
            ////    file1.Write(IQ.xx[i].ToString().Replace(",", "."));
            ////    file1.Write("; ");
            ////}
            //file1.Close();

            //Console.WriteLine(IQ.GreenFuncR(1,1,2,1));
            //System.IO.StreamWriter file4 = new System.IO.StreamWriter(@"C:\Users\tsybr\Desktop\данные\Screen.txt");
            //for (int i = 0; i < IQ.numsteps; i++)
            //{
            //    for (int j = 0; j < IQ.numsteps; j++)
            //    {
            //        file4.Write((Math.Sqrt(IQ.ScreenMatrix[i][j] * IQ.ScreenMatrix[i][j] + IQ.ScreenMatrix[i + IQ.numsteps][j] * IQ.ScreenMatrix[i + IQ.numsteps][j])).ToString().Replace(",", "."));
            //        file4.Write(" ");
            //    }
            //}
            //file4.Close();

            //System.IO.StreamWriter file5 = new System.IO.StreamWriter(@"C:\Users\tsybr\Desktop\данные\x.txt");
            //for (int i = 0; i < IQ.numsteps; i++)
            //{
            //    file5.Write(IQ.screenrho[i].ToString().Replace(",", "."));
            //    file5.Write(" ");
            //}
            //file5.Close();
            //System.IO.StreamWriter file6 = new System.IO.StreamWriter(@"C:\Users\tsybr\Desktop\данные\y.txt");
            //for (int i = 0; i < IQ.numsteps; i++)
            //{
            //    file6.Write(IQ.screeny[i].ToString().Replace(",", "."));
            //    file6.Write(" ");
            //}
            //file6.Close();
            Console.WriteLine("Complete!");
            Console.ReadLine();
        }
    }
}