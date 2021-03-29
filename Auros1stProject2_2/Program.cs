using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Numerics;
using System.Text;
using System.Threading.Tasks;
using static System.Console;
using static System.Math;

namespace Auros1stProject2_1
{
    class Program
    {
        static void Main(string[] args)
        {
            //
            // "SiO2 1000nm_on_Si.dat" 파일 로딩 후
            // 측정 스펙트럼 데이터를 alpha, beta 로 변환한다.
            //
            // 2021.03.24 이지원.
            //
            #region psi, delta -> alpha, beta

            List<string> MeasurementSpectrumData = new List<string>();  // 측정 스펙트럼 데이터 저장할 배열. (한 줄씩 저장)
            string[] SingleLineData;                                    // 한 줄의 스펙트럼 데이터를 임시로 저장할 배열.

            // "SiO2 2nm_on_Si.dat" 파일 읽기. (한 줄씩)
            MeasurementSpectrumData.AddRange(File.ReadAllLines("SiO2 1000nm_on_Si.dat"));

            // 무의미한 공백 행을 제거한다.
            int lenSpectrumData = MeasurementSpectrumData.Count;
            string Blank = "";
            for (int i = 0; i < lenSpectrumData; i++)
            {
                if (MeasurementSpectrumData[0] == Blank)
                    MeasurementSpectrumData.RemoveAt(0);
                else
                    break;
            }

            // wavelength : 350 ~ 980(nm)인 측정 스펙트럼 데이터를 담을 리스트 선언.
            List<double> wavelength_exp = new List<double>();   // 파장 데이터 리스트.
            List<double> AOI_exp = new List<double>();          // 입사각 데이터 리스트.
            List<double> alpha_exp = new List<double>();        // Psi 데이터 리스트.
            List<double> beta_exp = new List<double>();         // Delta 데이터 리스트.

            // 데이터의 첫번째 줄은 column 명이다.
            // 이를 제외하기 위해 반복문을 1부터 시작한다.
            int StartIndex = 1;
            int LenData = MeasurementSpectrumData.Count;
            for (int i = StartIndex; i < LenData; i++)
            {
                // tsv 형식의 데이터를 SingleLineData에 저장한다.
                SingleLineData = MeasurementSpectrumData[i].Split((char)0x09);  // 0x09 : 수평 탭.

                // 파장이 350 ~ 980(nm) 이내인 데이터만 저장한다.
                if (Convert.ToDouble(SingleLineData[0]) >= 350.0 &&
                    Convert.ToDouble(SingleLineData[0]) <= 980.0)
                {
                    // 각 컬럼에 해당하는 데이터를 저장한다.
                    wavelength_exp.Add(Double.Parse(SingleLineData[0]));
                    AOI_exp.Add(Double.Parse(SingleLineData[1]));
                    alpha_exp.Add(Double.Parse(SingleLineData[2]));
                    beta_exp.Add(Double.Parse(SingleLineData[3]));
                }
            }

            // psi, delta -> alpha, beta 변환.

            // degree, radian 변환 인라인 함수 정의.
            double degree2radian(double angle) => ((angle * (PI)) / 180.0);
            // double radian2degree(double angle) => (angle * (180.0 / PI));

            // Polarizer offset 각도. (45도)
            double PolarizerRadian = degree2radian(45.0);

            // psi, delta 데이터를 alpha, beta 로 변환한다.
            LenData = wavelength_exp.Count;
            for (int i = 0; i < LenData; i++)
            {
                // psi, delta 값을 radian 으로 변환한다.
                double PsiRadian = degree2radian(alpha_exp[i]);
                double DeltaRadian = degree2radian(beta_exp[i]);

                // psi, delta 데이터를 alpha, beta 로 갱신한다.
                alpha_exp[i] = (
                    (Pow(Tan(PsiRadian), 2.0) - Pow(Tan(PolarizerRadian), 2.0))
                    / (Pow(Tan(PsiRadian), 2.0) + Pow(Tan(PolarizerRadian), 2.0)));
                beta_exp[i] = (
                    (2.0 * Tan(PsiRadian) * Tan(PolarizerRadian) * Cos(DeltaRadian))
                    / (Pow(Tan(PsiRadian), 2.0) + Pow(Tan(PolarizerRadian), 2.0)));
            }

            #endregion

            //
            // "Si_new.txt", "SiO2_new.txt" 파일 물성값 로딩.
            //
            // 2021.03.24 이지원.
            //
            #region MyRegion

            // "Si_new.txt" 파일 읽기.
            string[] Si_new = File.ReadAllLines("Si_new.txt");  // Si 기판 물성값 저장.(한 줄씩)

            // 데이터의 첫번째 줄은 column 명이다.
            // 이를 제외하고 데이터를 받기 위해 LenData 변수를 선언한다.
            LenData = Si_new.Length - 1;
            double[] wavelength_Si = new double[LenData];
            double[] n_Si = new double[LenData];
            double[] k_Si = new double[LenData];

            // Si_new 에 받은 데이터를 각 컬럼별로 저장한다.
            LenData = Si_new.Length;
            for (int i = StartIndex; i < LenData; i++)
            {
                // tsv 형식의 데이터를 SingleLineData에 저장한다.
                SingleLineData = Si_new[i].Split((char)0x09);  // 0x09 : 수평 탭.

                // 각 컬럼에 해당하는 데이터를 저장한다.
                wavelength_Si[i - 1] = Double.Parse(SingleLineData[0]);
                n_Si[i - 1] = Double.Parse(SingleLineData[1]);
                k_Si[i - 1] = Double.Parse(SingleLineData[2]);
            }


            // "SiO2_new.txt" 파일 읽기.
            string[] SiO2_new = File.ReadAllLines("SiO2_new.txt");  // Si 기판 물성값 저장.(한 줄씩)

            // 데이터의 첫번째 줄은 column 명이다.
            // 이를 제외하고 데이터를 받기 위해 LenData 변수를 선언한다.
            LenData = SiO2_new.Length - 1;
            double[] wavelength_SiO2 = new double[LenData];
            double[] n_SiO2 = new double[LenData];
            double[] k_SiO2 = new double[LenData];

            // SiO2_new 에 받은 데이터를 각 컬럼별로 저장한다.
            LenData = SiO2_new.Length;
            for (int i = StartIndex; i < LenData; i++)
            {
                // tsv 형식의 데이터를 SingleLineData에 저장한다.
                SingleLineData = SiO2_new[i].Split((char)0x09);  // 0x09 : 수평 탭.

                // 각 컬럼에 해당하는 데이터를 저장한다.
                wavelength_SiO2[i - 1] = Double.Parse(SingleLineData[0]);
                n_SiO2[i - 1] = Double.Parse(SingleLineData[1]);
                k_SiO2[i - 1] = Double.Parse(SingleLineData[2]);
            }

            #region Si_new, SiO2_new 데이터 출력 (Test)
            /*LenData = wavelength_Si.Length;
            for (int i = 0; i < LenData; i++)
                WriteLine(wavelength_Si[i] + "\t" + n_Si[i] + "\t" + k_Si[i]);
            WriteLine("============================================");
            for (int i = 0; i < LenData; i++)
                WriteLine(wavelength_SiO2[i] + "\t" + n_SiO2[i] + "\t" + k_SiO2[i]);*/
            #endregion
            #endregion

            //
            // "Si_new.txt", "SiO2_new.txt" 의 n, k 를 사용하여
            // 각 계면에서의 반사, 투과계수를 계산한다.
            //
            // 2021.03.24 이지원.
            //
            #region 각 계면에서의 반사, 투과계수 계산

            LenData = wavelength_Si.Length;

            // 반사계수를 담을 배열.
            Complex[] r12p = new Complex[LenData],
                      r12s = new Complex[LenData],
                      r01p = new Complex[LenData],
                      r01s = new Complex[LenData];
            // 투과계수를 담을 배열.
            Complex[] t12p = new Complex[LenData],
                      t12s = new Complex[LenData],
                      t01p = new Complex[LenData],
                      t01s = new Complex[LenData];

            double AOI_air = degree2radian(65.0);  // 입사각. (라디안) 
            Complex N_air = new Complex(1.0, 0);    // 공기의 굴절률.

            // 반사, 투과계수를 계산한다.
            for (int i = 0; i < LenData; i++)
            {
                // 파장에 대한 물질의 복소굴절률을 구한다.
                Complex N_SiO2 = new Complex(n_SiO2[i], -k_SiO2[i]);
                Complex N_Si = new Complex(n_Si[i], -k_Si[i]);

                // air, SiO2 경계면에서의 굴절각을 구한다. (스넬의 법칙)
                Complex Sintheta_j = new Complex(Sin(AOI_air), 0);
                Complex Costheta_j = new Complex(Cos(AOI_air), 0);
                Complex Sintheta_k = (N_air / N_SiO2) * Sintheta_j;
                Complex theta_k = Complex.Asin(Sintheta_k);
                // air, SiO2 경계면에서의 굴절각.
                Complex Costheta_k = Complex.Cos(theta_k);

                // air, SiO2 경계면에서의 반사계수를 구한다.
                r01p[i] = ((N_SiO2 * Costheta_j) - (N_air * Costheta_k)) /
                               ((N_SiO2 * Costheta_j) + (N_air * Costheta_k));

                r01s[i] = ((N_air * Costheta_j) - (N_SiO2 * Costheta_k)) /
                               ((N_air * Costheta_j) + (N_SiO2 * Costheta_k));

                // air, SiO2 경계면에서의 투과계수를 구한다.
                t01p[i] = (N_air * Costheta_j * 2.0) /
                               ((N_SiO2 * Costheta_j) + (N_air * Costheta_k));

                t01s[i] = (N_air * Costheta_j * 2.0) /
                               ((N_air * Costheta_j) + (N_SiO2 * Costheta_k));

                // SiO2, Si 경계면에서의 굴절각을 구한다. (스넬의 법칙)
                Sintheta_j = Complex.Sin(theta_k);
                Costheta_j = Complex.Cos(theta_k);
                Sintheta_k = (N_SiO2 / N_Si) * Sintheta_j;
                theta_k = Complex.Asin(Sintheta_k);
                Costheta_k = Complex.Cos(theta_k);

                // SiO2, Si 경계면에서의 반사계수를 구한다.
                r12p[i] = ((N_Si * Costheta_j) - (N_SiO2 * Costheta_k)) /
                             ((N_Si * Costheta_j) + (N_SiO2 * Costheta_k));

                r12s[i] = ((N_SiO2 * Costheta_j) - (N_Si * Costheta_k)) /
                             ((N_SiO2 * Costheta_j) + (N_Si * Costheta_k));

                // SiO2, Si 경계면에서의 투과계수를 구한다.
                t12p[i] = (N_SiO2 * Costheta_j * 2.0) /
                             ((N_Si * Costheta_j) + (N_SiO2 * Costheta_k));

                t12s[i] = (N_SiO2 * Costheta_j * 2.0) /
                             ((N_SiO2 * Costheta_j) + (N_Si * Costheta_k));
            }

            #region 위에서 구한 반사, 투과계수 출력 (Test)
            /*WriteLine("====== air, SiO2 경계 ======");
            for (int i = 0; i < LenData; i++)
            {
                WriteLine(
                    r01p[i] + " " +
                    r01s[i] + " " +
                    t01p[i] + " " +
                    t01s[i]);
            }
            WriteLine("====== SiO2, Si 경계 ======");
            for (int i = 0; i < LenData; i++)
            {
                WriteLine(
                    r12p[i] + " " +
                    r12s[i] + " " +
                    t12p[i] + " " +
                    t12s[i]);
            }*/
            #endregion

            #endregion

            //
            // 위상 두께를 구하고 위에서 구한 반사계수를 통해
            // 무한등비급수 수렴식을 계산한다. => 총 반사계수를 구한다.
            //
            // 2021.03.24 이지원. 
            //
            #region 총 반사계수를 구한다.

            // 두께 범위와 두께 간격을 설정한다.
            double StartThickness = 700.0;
            double EndThickness = 1300.0;
            double gap = 6.0;

            // MSE 와 두께를 담을 배열을 선언, 초기화한다.
            double numMSE = (EndThickness - StartThickness) / gap + 1;
            double[] MSEs        = new double[(int)numMSE];
            double[] thicknesses = new double[(int)numMSE];

            // 두께별 MSE 를 계산해서 MSEs 배열에 저장한다.
            int idx = 0;
            for (double thickness = StartThickness; thickness <= EndThickness; thickness += gap)
            {
                // 총 반사계수를 저장할 배열 선언.
                Complex[] Rp = new Complex[LenData],
                          Rs = new Complex[LenData];

                for (int i = 0; i < LenData; i++)
                {
                    // SiO2의 복소 굴절률.
                    Complex N_SiO2 = new Complex(n_SiO2[i], -k_SiO2[i]);

                    // air, SiO2 경계면에서의 굴절각을 구한다. (스넬의 법칙)
                    Complex Sintheta_j = new Complex(Sin((double)AOI_air), 0);
                    Complex Sintheta_k = (N_air / N_SiO2) * Sintheta_j;
                    Complex theta_k = Complex.Asin(Sintheta_k);         // air, SiO2 경계면에서의 굴절각.
                    Complex Costheta_k = Complex.Cos(theta_k);

                    // 위상 두께를 구한다.
                    Complex PhaseThickness = ((double)thickness * Math.PI * 2.0 / wavelength_SiO2[i]) * N_SiO2 * Costheta_k;

                    //WriteLine(PhaseThickness);

                    // 총 반사계수를 구한다.
                    Complex E = Complex.Exp(PhaseThickness * new Complex(0, -2.0));

                    Rp[i] = (r01p[i] + r12p[i] * E) /
                            (1 + r01p[i] * r12p[i] * E);

                    Rs[i] = (r01s[i] + r12s[i] * E) /
                            (1 + r01s[i] * r12s[i] * E);

                }

                #region 총 반사계수 출력 (Test)

                /*for (int i = 0; i < LenData; i++)
                    WriteLine(Rp[i] + "\t" + Rs[i]);*/

                #endregion

                #endregion

                //
                // 위에서 구한 총 반사계수로부터 alpha, beta 를 구한다.
                //
                // 2021.03.24 이지원.
                //
                #region 총 반사계수로부터 alpha, beta 도출.

                // alpha, beta 이론값을 담을 배열 선언.
                double[] alpha_cal = new double[LenData],
                         beta_cal = new double[LenData];

                // Polarizer 오프셋 각.
                double polarizerAngle = degree2radian(45.0);

                for (int i = 0; i < LenData; i++)
                {
                    // 총 반사계수비. (복소반사계수비)
                    Complex rho = Rp[i] / Rs[i];

                    // Psi, Delta.
                    double Psi = Atan(rho.Magnitude);
                    double Delta = rho.Phase;


                    alpha_cal[i] = (Pow(Tan(Psi), 2.0) - Pow(Tan(polarizerAngle), 2.0)) /
                                           (Pow(Tan(Psi), 2.0) + Pow(Tan(polarizerAngle), 2.0));

                    beta_cal[i] = (2.0 * Tan(Psi) * Cos(Delta) * Tan(polarizerAngle)) /
                                           (Pow(Tan(Psi), 2.0) + Pow(Tan(polarizerAngle), 2.0));
                }

                #region alpha, beta 이론값 출력 (Test)

                /*for (int i = 0; i < LenData; i++)
                    WriteLine(alpha_cal[i] + " " 
                    + beta_cal[i]);*/

                #endregion

                #endregion

                //
                // 측정값과 이론값의 MSE 를 계산한다.
                //
                // 2021.03.24 이지원.
                //
                #region 측정값과 이론값의 MSE 계산.

                double sum = 0;
                for (int i = 0; i < LenData; i++)
                {
                    double difference_MSE =
                         Pow((alpha_exp[i] - alpha_cal[i]), 2.0) +
                         Pow((beta_exp[i] - beta_cal[i]), 2.0);
                    sum += difference_MSE;

                }

                MSEs[idx] = sum / LenData;
                thicknesses[idx] = thickness;
                ++idx;

                #endregion
            }
            #region MSE 배열 출력

            /*foreach (var item in MSEs)
                WriteLine(item);
            WriteLine("==========================");
            foreach (var item in thicknesses)
                WriteLine(item);*/

            #endregion

            //
            // 두께별로 계산된 MSE 를 통해
            // global minimum 일 때의 두께, MSE 값을 찾는다.
            //
            // 2021.03.26 이지원.
            //
            #region global minimum 탐색.

            // MSEs 에서 global minimum 에서의 MSE, 두께값을 구한다.
            int idxGlobalMinimum = 0;           // global minimum 에서의 index.
            double GlobalMinimum = MSEs.Min();  // global minimum 값.

            // global minimum 에서의 index 를 찾는다.
            LenData = MSEs.Length;
            for (int i = 0; i < LenData; i++)
            {
                if (MSEs[i] == GlobalMinimum)
                {
                    idxGlobalMinimum = i;
                    break;
                }
            }
            // WriteLine($"{MSEs[idxGlobalMinimum]}  {idxGlobalMinimum}");

            // global minimum 에서의 두께 값.
            double d0 = thicknesses[idxGlobalMinimum];

            #endregion

            //
            // d0 근처에 있는 "실제" global minimum 이 되는 두께 d_sol 을 찾는다.
            // MSE 값의 변화가 10^-5 이하가 되면, 최적화 수행을 멈춘다.
            //
            // 2021.03.29 이지원.
            //
            #region 두께 범위를 재설정하고 MSE를 구한다.

            double d_sol = 0.0; // "실제" global minimum 이 되는 두께.
            int cnt = 0;        // global minimum 을 찾기 위한 연산 반복 횟수.
            //double beforeGM = GlobalMinimum;
            while (true)
            {
                ++cnt;
                WriteLine($"==== {cnt} 회차 ====");

                // 두께 범위와 두께 간격을 정한다. 
                StartThickness = d0 - gap;
                EndThickness = d0 + gap;
                gap = gap * 0.7;

                // MSEs 배열을 초기화한다.
                // double minus = EndThickness - StartThickness;
                numMSE = (EndThickness - StartThickness) / gap + 1;
                MSEs = new double[(int)numMSE];
                thicknesses = new double[(int)numMSE];

                // 새로 지정된 두께 범위에서 MSE 값들을 계산한다.
                idx = 0;
                for (double thickness = StartThickness; thickness <= EndThickness; thickness += gap)
                {
                    #region 총 반사계수를 구한다.
                    // 총 반사계수를 저장할 배열 선언.
                    Complex[] Rp = new Complex[LenData],
                              Rs = new Complex[LenData];

                    for (int i = 0; i < LenData; i++)
                    {
                        // SiO2의 복소 굴절률.
                        Complex N_SiO2 = new Complex(n_SiO2[i], -k_SiO2[i]);

                        // air, SiO2 경계면에서의 굴절각을 구한다. (스넬의 법칙)
                        Complex Sintheta_j = new Complex(Sin((double)AOI_air), 0);
                        Complex Sintheta_k = (N_air / N_SiO2) * Sintheta_j;
                        Complex theta_k = Complex.Asin(Sintheta_k);
                        // air, SiO2 경계면에서의 굴절각.
                        Complex Costheta_k = Complex.Cos(theta_k);

                        // 위상 두께를 구한다.
                        Complex PhaseThickness = ((double)thickness * Math.PI * 2.0 / wavelength_SiO2[i]) * N_SiO2 * Costheta_k;

                        //WriteLine(PhaseThickness);

                        // 총 반사계수를 구한다.
                        Complex E = Complex.Exp(PhaseThickness * new Complex(0, -2.0));

                        Rp[i] = (r01p[i] + r12p[i] * E) /
                                (1 + r01p[i] * r12p[i] * E);

                        Rs[i] = (r01s[i] + r12s[i] * E) /
                                (1 + r01s[i] * r12s[i] * E);

                    }
                    #endregion
                    #region 총 반사계수로부터 alpha, beta 를 계산한다.
                    // alpha, beta 이론값을 담을 배열 선언.
                    double[] alpha_cal = new double[LenData],
                             beta_cal = new double[LenData];

                    // Polarizer 오프셋 각.
                    double polarizerAngle = degree2radian(45.0);

                    for (int i = 0; i < LenData; i++)
                    {
                        // 총 반사계수비. (복소반사계수비)
                        Complex rho = Rp[i] / Rs[i];

                        // Psi, Delta.
                        double Psi = Atan(rho.Magnitude);
                        double Delta = rho.Phase;


                        alpha_cal[i] = (Pow(Tan(Psi), 2.0) - Pow(Tan(polarizerAngle), 2.0)) /
                                               (Pow(Tan(Psi), 2.0) + Pow(Tan(polarizerAngle), 2.0));

                        beta_cal[i] = (2.0 * Tan(Psi) * Cos(Delta) * Tan(polarizerAngle)) /
                                               (Pow(Tan(Psi), 2.0) + Pow(Tan(polarizerAngle), 2.0));
                    }
                    #endregion
                    #region MSE 를 계산한다.

                    double sum = 0;
                    for (int i = 0; i < LenData; i++)
                    {
                        double difference_MSE =
                             Pow((alpha_exp[i] - alpha_cal[i]), 2.0) +
                             Pow((beta_exp[i] - beta_cal[i]), 2.0);
                        sum += difference_MSE;

                    }

                    MSEs[idx] = sum / LenData;
                    thicknesses[idx] = thickness;
                    ++idx;

                    #endregion
                }
                // MSEs 출력.
                foreach (var item in MSEs)
                    WriteLine(item);

                // global minimum 에서의 MSE, 두께 값을 구한다.
                idxGlobalMinimum = 0;               // global minimum 에서의 index.
                GlobalMinimum = MSEs.Min();      // global minimum 값.

                // global minimum 에서의 index 를 찾는다.
                LenData = MSEs.Length;
                for (int i = 0; i < LenData; i++)
                {
                    if (MSEs[i] == GlobalMinimum)
                    {
                        idxGlobalMinimum = i;
                        break;
                    }
                }

                d0 = thicknesses[idxGlobalMinimum]; // global minimum 에서의 두께 값.


                // MSE 값의 간격을 확인한다.
                // 1. 양 옆.
                for (int i = 1; i < LenData; i++)
                {
                    // MSE 의 변화율이 10^-5 이하이면 while 문을 탈출한다.
                    if (Math.Abs((MSEs[i] - MSEs[i - 1])) <= 0.00001)
                    {
                        d_sol = d0;
                        goto Find_d_sol;
                    }
                }

                // 2. global minimum 끼리 비교.
                /*for (int i = 1; i < LenData; i++)
                {
                    // MSE 의 변화율이 10^-5 이하이면 while 문을 탈출한다.
                    if (Math.Abs((beforeGM - GlobalMinimum)) <= 0.00001)
                    {
                        d_sol = d0;
                        goto Find_d_sol;
                    }
                }
                beforeGM = GlobalMinimum;*/

            }
            Find_d_sol: WriteLine($"d_sol: {d_sol}    MSE: {MSEs[idxGlobalMinimum]}    idxGlobalMinimum: {idxGlobalMinimum}");

            #endregion

        }

    }
}
