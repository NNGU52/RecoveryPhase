using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;
using System.Windows.Forms.DataVisualization.Charting;

using System.Numerics;

namespace RecoveryPhase
{
    public partial class Form1 : Form
    {
        int length;
        int calc;
        double TAU, Er, Eps;
        PointD[] signal;            // массив для точек сигнала
        PointD[] signal_recover;    // массив для точек сигнала
        double[] spectr_phase;      //массив для точек спектра(типа double)
        double[] spectr_magnitude;  //массив для точек спектра(типа double)
        double[] recover;           //массив для точек спектра(типа double)
        Cmplx[] F_recover;          //комплексный массив
        double modul;
        bool Timer1, init, stop_Timer1; // переменные к анимации
        bool Timer2, stop_Timer2;

        // структура, для хранения точек сигнала типв double
        public struct PointD 
        {
            public double X;
            public double Y;
        }

        // структура для комплексного числа
        public struct Cmplx
        {
            public double real;
            public double image;
        }

        public Form1()
        {
            InitializeComponent();

            comboBoxLength.Text = "512";
            init = true;
            Timer1 = true;
            Timer2 = true;
            calc = 0;
        }

        void Init()
        {
            Recover.Enabled = true;
            Shift.Enabled = false;
            Stop.Enabled = false;

        }

        // функция рисования (принимающая массив типа PoitD)
        private void SendChart(Chart ChartI, PointD[] data, string p)
        {
            ChartI.Series.Clear(); // очистка коллекции Series

            double Maxval = 0.0;
            Series NewSeries = new Series(p);
            for (int i = 0; i < data.Length; i++)
            {
                NewSeries.Points.AddXY(i, data[i].Y);

                if (data[i].Y > Maxval)
                {
                    Maxval = data[i].Y;
                }
            }
            NewSeries.ChartType = SeriesChartType.Line;
            NewSeries.BorderWidth = 2;

            ChartI.Series.Add(NewSeries);

            ChartI.ChartAreas[0].AxisX.Minimum = 0;
            ChartI.ChartAreas[0].AxisX.Maximum = data.Length;
            ChartI.ChartAreas[0].AxisY.Maximum = Maxval + 1;

            ChartI.Invalidate();  // метод, вызывающий перерисовку chart
        }

        // функция рисования (принимающая массив типа double)
        private void SendChart_(Chart ChartI, double[] data, string p)
        {
            ChartI.Series.Clear(); // очистка коллекции Series

            double Maxval = 0.0;
            Series NewSeries = new Series(p);
            for (int i = 0; i < data.Length; i++)
            {
                NewSeries.Points.AddXY(i, data[i]);

                if (data[i] > Maxval)
                {
                    Maxval = data[i];
                }
            }
            NewSeries.ChartType = SeriesChartType.Line;
            NewSeries.BorderWidth = 1;

            ChartI.Series.Add(NewSeries);

            ChartI.ChartAreas[0].AxisX.Minimum = 0;
            ChartI.ChartAreas[0].AxisX.Maximum = data.Length;
            ChartI.ChartAreas[0].AxisY.Maximum = Maxval + 1;

            ChartI.Invalidate();  // метод, вызывающий перерисовку chart
        }

        private void SecondaryChart(Chart ChartI, PointD[] Data, string p)
        {
            Double Maxval = ChartI.ChartAreas[0].AxisY.Maximum;

            Series NewSeries = new Series(p);
            for (UInt16 i = 0; i < Data.Length; i++)
            {
                NewSeries.Points.AddXY(i, Data[i].Y);

                if (Data[i].Y > Maxval)
                {
                    Maxval = Data[i].Y;
                }
            }
            NewSeries.ChartType = SeriesChartType.Line;
            NewSeries.BorderWidth = 3;

            ChartI.Series.Add(NewSeries);

            ChartI.ChartAreas[0].BorderWidth = 3;
            ChartI.ChartAreas[0].AxisY.Maximum = Maxval * 1.05;
            ChartI.Invalidate();
        }

        // функуция генерации значений сигнала по двум осям
        PointD[] GeneratingPointsSignal(double A, double M, double D, int length, bool b)  
        {
            PointD[] points = new PointD[length];

            for (int i = 0; i < points.Length; i++)
            {
                points[i].X = Convert.ToDouble(i);

                if (b)
                    points[i].Y = A * Math.Exp(-1 * Math.Pow((points[i].X - M), 2) / Math.Pow(D, 2));
                else
                    points[i].Y = 0;
            }

            return points;
        }

        // функция суммирования трех сигналов
        PointD[] SignalSummation(PointD[] signal1, PointD[] signal2, PointD[] signal3, PointD[] signal4, PointD[] signal5)  
        {
            PointD[] sum = new PointD[signal1.Length];

            for (int i = 0; i < sum.Length; i++)
            {
                sum[i].X = signal1[i].X + signal2[i].X + signal3[i].X + signal4[i].X + signal5[i].X;
                sum[i].Y = signal1[i].Y + signal2[i].Y + signal3[i].Y + signal4[i].Y + signal5[i].Y;
            }

            return sum;
        }

        // функция для Фурье - преобразований, s=1 - обратное, s=-1 -прямое
        public static void Fourea(Cmplx[] sig, int n, int s)
        {
            int i, j, istep;
            int m, mmax;
            double r, r1, theta, w_r, w_i, temp_r, temp_i;
            double pi = 3.1415926f;

            r = pi * s;
            j = 0;
            for (i = 0; i < n; i++)
            {
                if (i < j)
                {

                    temp_r = sig[j].real;
                    temp_i = sig[j].image;
                    sig[j].real = sig[i].real;
                    sig[j].image = sig[i].image;
                    sig[i].real = temp_r;
                    sig[i].image = temp_i;
                }
                m = n >> 1;
                while (j >= m) { j -= m; m = (m + 1) / 2; }
                j += m;
            }
            mmax = 1;
            while (mmax < n)
            {
                istep = mmax << 1;
                r1 = r / (double)mmax;
                for (m = 0; m < mmax; m++)
                {
                    theta = r1 * m;
                    w_r = (double)Math.Cos((double)theta);
                    w_i = (double)Math.Sin((double)theta);
                    for (i = m; i < n; i += istep)
                    {

                        j = i + mmax;
                        temp_r = w_r * sig[j].real - w_i * sig[j].image;
                        temp_i = w_r * sig[j].image + w_i * sig[j].real;
                        sig[j].real = sig[i].real - temp_r;
                        sig[j].image = sig[i].image - temp_i;
                        sig[i].real += temp_r;
                        sig[i].image += temp_i;

                    }
                }
                mmax = istep;
            }
            if (s > 0)
                for (i = 0; i < n; i++)
                {
                    sig[i].real /= (Convert.ToDouble(n));
                    sig[i].image /= (Convert.ToDouble(n));
                }

        }

        // амплитудный спектр
        public static double[] FuncSpectrMagnitude(int length, PointD[] Signal)
        {
            Cmplx[] Cmplx = new Cmplx[length + 1];
            double[] Spectr = new double[length + 1];

            for (int i = 0; i < length; i++)
            {
                Cmplx[i].real = (float)Signal[i].Y;
                Cmplx[i].image = 0;
            }
            Fourea(Cmplx, length, -1);
            for (int i = 0; i < length; i++)
            {
                Spectr[i] = Math.Sqrt(Cmplx[i].real * Cmplx[i].real + Cmplx[i].image * Cmplx[i].image);
            }
            return Spectr;
        }

        // фазовый спектр
        public static double[] FuncSpectrPhase(int length, PointD[] Signal)
        {
            Cmplx[] Cmplx = new Cmplx[length + 1];
            double[] Spectr = new double[length + 1];

            for (int i = 0; i < length; i++)
            {
                Cmplx[i].real = (float)Signal[i].Y;
                Cmplx[i].image = 0;
            }
            Fourea(Cmplx, length, -1);
            for (int i = 0; i < length; i++)
            {
                Spectr[i] = Math.Atan(Cmplx[i].real / Cmplx[i].image);
            }
            return Spectr;
        }

        public static void InitRec(int len, int Prec, double[] spectr, Cmplx[] F_recover, double[] recover, PointD[] signal_recover)
        {
            var min = 0.0;                  // диапазон для генерации чисел
            var max = 2 * Math.PI;          // в интервале от 0 до 2Пи
            double phase = 0;
            var rnd = new Random();         // элемент класса Random
            F_recover = new Cmplx[len];


            for (int i = 0; i < len; i++)
            {
                phase = rnd.NextDouble() * (max - min) + min; // генерируем фазу от 0 до 2пи
                F_recover[i].real = spectr[i] * Math.Cos(phase);
                F_recover[i].image = spectr[i] * Math.Sin(phase);
            }

            Fourea(F_recover, len, 1);  // выполняем обратное преобразование Фурье


            for (int i = 0; i < len; i++)
            {
                recover[i] = F_recover[i].real;
                if (recover[i] < 0)
                    recover[i] = 0;
                signal_recover[i].Y = recover[i];
            }

        }

        private void Generating_Click(object sender, EventArgs e)
        {
            Init();
            // первая вкладка
            signal = SignalSummation(GeneratingPointsSignal(Convert.ToDouble(textBox1A.Text), Convert.ToDouble(textBox1M.Text), Convert.ToDouble(textBox1D.Text), Convert.ToInt32(comboBoxLength.Text), checkBox1.Checked),
                GeneratingPointsSignal(Convert.ToDouble(textBox2A.Text), Convert.ToDouble(textBox2M.Text), Convert.ToDouble(textBox2D.Text), Convert.ToInt32(comboBoxLength.Text), checkBox2.Checked),
                GeneratingPointsSignal(Convert.ToDouble(textBox3A.Text), Convert.ToDouble(textBox3M.Text), Convert.ToDouble(textBox3D.Text), Convert.ToInt32(comboBoxLength.Text), checkBox3.Checked),
                GeneratingPointsSignal(Convert.ToDouble(textBox4A.Text), Convert.ToDouble(textBox4M.Text), Convert.ToDouble(textBox4D.Text), Convert.ToInt32(comboBoxLength.Text), checkBox4.Checked),
                GeneratingPointsSignal(Convert.ToDouble(textBox5A.Text), Convert.ToDouble(textBox5M.Text), Convert.ToDouble(textBox5D.Text), Convert.ToInt32(comboBoxLength.Text), checkBox5.Checked));

                
            for (int i = 0; i < signal.Length; i++)
            {
                signal[i].X = Convert.ToDouble(i);
                signal[i].Y = signal[i].Y;
            }

            spectr_magnitude = FuncSpectrMagnitude(Convert.ToInt32(comboBoxLength.Text), signal);   //вычисляем модуль спектра
            spectr_phase = FuncSpectrPhase(Convert.ToInt32(comboBoxLength.Text), signal);

            SendChart(chart1, signal, "Исходный сигнал");

            // вторая вкладка
            SendChart_(chart2, spectr_magnitude, "Амплитудный спектр исходного сигнала");

            // третья вкладка
            SendChart_(chart3, spectr_phase, "Фазовый спектр исходного сигнала");

        }

        private void Shift_Click(object sender, EventArgs e)
        {
            if (Timer2 == true)
            {
                stop_Timer2 = false;
                timer_shift.Enabled = true;
                Timer2 = false;
                Shift.Text = "Стоп";
            }
            else
            {
                timer_shift.Enabled = false;
                Timer2 = true;
                Shift.Text = "Сдвиг";
            }
        }

        private void Stop_Click(object sender, EventArgs e)
        {
            timer_data.Stop();
            Timer1 = true;
            MessageBox.Show(
            "The process of recovery have been cancelled after " + calc + " iterations.",
            "Abort", MessageBoxButtons.OK, MessageBoxIcon.Exclamation);
        }

        private void Fault_Click(object sender, EventArgs e)
        {
            calc = 0;
            Er = 0;
            Eps = 0;

            chart1.Series[1].Points.Clear();

            init = true;
            timer_data.Enabled = false;
            timer_shift.Enabled = false;
            Timer1 = true;
            stop_Timer1 = false;
            Timer2 = true;
            stop_Timer2 = false;
            Recover.Enabled = true;
            Generating.Enabled = true;

            Shift.Text = "Сдвиг";
            Shift.Enabled = false;
            Stop.Enabled = false;

        }

        private void Recover_Click(object sender, EventArgs e)
        {
            timer_shift.Enabled = false;
            Recover.Enabled = true;
            timer_data.Enabled = true;
            int Prec = UInt16.Parse(textBox_correct.Text);  //степень точности вычислений
            length = Convert.ToInt32(comboBoxLength.Text);  //длина сигнала, вынимаем из текстбокса
            signal_recover = new PointD[length];            //наш восстановленный сигнал
            F_recover = new Cmplx[length];                  //комплексный массив для восстановления
            recover = new Double[length];
            if (init == true)
            {
                InitRec(length, Prec, spectr_magnitude, F_recover, recover, signal_recover);    // Инцилизация для выполнения алгоритма
                init = false;
            }
            Stop.Enabled = true;
        }

        private void timer_data_Tick(object sender, EventArgs e)
        {
            signal = SignalSummation(GeneratingPointsSignal(Convert.ToDouble(textBox1A.Text), Convert.ToDouble(textBox1M.Text), Convert.ToDouble(textBox1D.Text), Convert.ToInt32(comboBoxLength.Text), checkBox1.Checked),
                GeneratingPointsSignal(Convert.ToDouble(textBox2A.Text), Convert.ToDouble(textBox2M.Text), Convert.ToDouble(textBox2D.Text), Convert.ToInt32(comboBoxLength.Text), checkBox2.Checked),
                GeneratingPointsSignal(Convert.ToDouble(textBox3A.Text), Convert.ToDouble(textBox3M.Text), Convert.ToDouble(textBox3D.Text), Convert.ToInt32(comboBoxLength.Text), checkBox3.Checked),
                GeneratingPointsSignal(Convert.ToDouble(textBox4A.Text), Convert.ToDouble(textBox4M.Text), Convert.ToDouble(textBox4D.Text), Convert.ToInt32(comboBoxLength.Text), checkBox4.Checked),
                GeneratingPointsSignal(Convert.ToDouble(textBox5A.Text), Convert.ToDouble(textBox5M.Text), Convert.ToDouble(textBox5D.Text), Convert.ToInt32(comboBoxLength.Text), checkBox5.Checked));


            for (int i = 0; i < signal.Length; i++)
            {
                signal[i].X = Convert.ToDouble(i);
                signal[i].Y = signal[i].Y;
            }

            for (int i = 0; i < length; i++)
            {
                F_recover[i].real = recover[i];
                F_recover[i].image = 0;
            }

            Fourea(F_recover, length, -1);  // прямое преобразование Фурье

            for (int i = 0; i < length; i++)
            {
                modul = Math.Sqrt(F_recover[i].real * F_recover[i].real + F_recover[i].image * F_recover[i].image); //новый модуль
                F_recover[i].real = F_recover[i].real * spectr_magnitude[i] / modul;                                // меняем модуль комплексного числа
                F_recover[i].image = F_recover[i].image * spectr_magnitude[i] / modul;                              // на модуль исходного сигнала
            }

            Fourea(F_recover, length, 1);       //обратное преобразование Фурье


            for (int i = 0; i < length; i++)    //зануляем отрицательную часть воостановленого сигнала
            {
                recover[i] = F_recover[i].real;
                if (recover[i] < 0)
                    recover[i] = 0;

            }

            Er = 0;
            //считаем СКО
            for (int i = 0; i < length; i++)
            {
                Er += (signal_recover[i].Y - recover[i]) * (signal_recover[i].Y - recover[i]) / length;
            }

            textBox_error.Text = String.Format("{0:E}", Er); //Выводит отклонение на экран

            for (int i = 0; i < length; i++)                //запоминаем предыдущее значение, чтобы использовать на следующей итерации
            {
                signal_recover[i].Y = recover[i];
            }

            Int32 Prec = Int32.Parse(textBox_corect.Text);  //степень точности вычислений, вынимаем из текстбокса
            TAU = Math.Pow(10.0, (-(double)Prec));          // Точность вычислений    
            if (Er < TAU)                                   //пока наше СКО отклонение не станет меньше задданой точности, то продолжаем процесс
            {
                stop_Timer1 = true;
            }


            calc++; //счетчик итераций данного цикла
 
            SendChart(chart1, signal, " Исходный сигнал ");
            SecondaryChart(chart1, signal_recover, " Восстановленный сигнал ");
            //SendChart_(chart2, spectr, " Спектр исходного сигнала ");
            //SecondaryChart_(chart2, modul, " Восстановленный сигнал ");
            if (stop_Timer1 == true)
            {
                timer_data.Enabled = false;
                Recover.Enabled = false;
                Shift.Enabled = true;
                Stop.Enabled = false;

                //Stop.Text = "Остановить";
                MessageBox.Show(
                   "The process of recovery have finished working after " + calc + " iterations.",
                   "Success", MessageBoxButtons.OK, MessageBoxIcon.Information);
            }
        }

        private void timer_shift_Tick(object sender, EventArgs e)
        {
            length = Convert.ToInt32(comboBoxLength.Text);
            chart1.Series[1].Points.Clear();

            // циклический сдвиг
            for (int i = 0; i < length; i++)
            {
                if (i + 1 < length)
                    signal_recover[i].Y = signal_recover[i + 1].Y;
                else
                    signal_recover[i].Y = signal_recover[length - i + 1].Y;
            }

            for (int i = 1; i < length; i++)
            {
                chart1.Series[1].Points.AddXY(i, signal_recover[i].Y);  //отрисовка
            }

            Eps = 0;
            /*for (int i = 0; i < length; i++)
            {
                Eps += (signal[i].Y - signal_recover[i].Y) * (signal[i].Y - signal_recover[i].Y) / length;//СКО
            }*/

            if (stop_Timer2 == true)
            {
                timer_shift.Enabled = false;
                Timer2 = false;
                //Shift.Text = "Сдвиг";
            }

            //textBox_error.Text = String.Format("{0:E}", Eps);//Выводит невязки на экран
        } 

        private void Form1_Load(object sender, EventArgs e)
        {
            Fault_Click(sender, e);
        }

        private void Close_Click(object sender, EventArgs e)
        {
            Close();
        }

        private void checkBox_reflex_CheckedChanged(object sender, EventArgs e)
        {
            chart1.Series[1].Points.Clear();

            double buf;
            for (int i = 0; i < length / 2; i++)
            {
                buf = signal_recover[i].Y;
                signal_recover[i].Y = signal_recover[length - i - 1].Y;
                signal_recover[length - i - 1].Y = buf;
            }

            for (int i = 1; i < length; i++)
            {
                chart1.Series[1].Points.AddXY(i, signal_recover[i].Y);
            }
        }
    }
}
