using System;
using System.Text;
using System.Collections.Generic;
using System.Linq;

namespace MyApp
{
    class Program
    {
        static void Main(string[] args)
        {
            SimulationParameters p = new SimulationParameters();
            MontecarloSimimulator euro = new MontecarloSimimulator();
            Console.Write("What is the Underlying Price: ");
            p.S0 = Convert.ToDouble(Console.ReadLine());
            Console.Write("What is the Strike Price: ");
            p.Strike = Convert.ToDouble(Console.ReadLine());
            Console.Write("How much years the Option expired: ");
            p.Tenor = Convert.ToDouble(Console.ReadLine());
            Console.Write("What is the risk free rate: ");
            p.r = Convert.ToDouble(Console.ReadLine());
            Console.Write("What is the Volatility: ");
            p.Volatility = Convert.ToDouble(Console.ReadLine());
            Console.Write("Is Call? (true or false) ");
            euro.IsCall = Convert.ToBoolean(Console.ReadLine());
            Console.Write("Do you want Antithetic method for montecarlo? (true or false) ");
            euro.IsAntithetic = Convert.ToBoolean(Console.ReadLine());
            Console.Write("How many simulations you want ");
            p.Simulations = Convert.ToInt32(Console.ReadLine());
            Console.Write("How many steps you want ");
            p.Steps = Convert.ToInt32(Console.ReadLine());
            if (p.Steps == 1)
            {   
                VDCparameters vp = new VDCparameters();
                Console.Write("Do you want Van der Corput functionality? (true or false)");
                vp.IsVDP = Convert.ToBoolean(Console.ReadLine());
                if (vp.IsVDP == true)
                {
                    Console.Write("What base you want for your Van der Corput functionality? ");
                    vp.Base = Convert.ToInt32(Console.ReadLine());
                    OptionResult VDCresult = new OptionResult();
                    VDCfunc VDCs = new VDCfunc();
                    VDCresult.Price = VDCs.VDCOptionPrice(vp.Base, p, euro.IsCall);
                    (VDCresult.Delta, VDCresult.Gamma) = VDCs.DeltaandGamma(vp.Base, p, euro.IsCall);
                    VDCresult.Theta = VDCs.ThetaValue(vp.Base, p, euro.IsCall);
                    VDCresult.Vega = VDCs.VegaValue(vp.Base, p, euro.IsCall);
                    VDCresult.Rho = VDCs.RhoValue(vp.Base, p, euro.IsCall);
                    Console.WriteLine("Option Price: ");
                    Console.WriteLine(VDCresult.Price);
                    Console.WriteLine("Option Delta: ");
                    Console.WriteLine(VDCresult.Delta);
                    Console.WriteLine("Option Gamma: ");
                    Console.WriteLine(VDCresult.Gamma);
                    Console.WriteLine("Option Vega: ");
                    Console.WriteLine(VDCresult.Vega);
                    Console.WriteLine("Option Theta: ");
                    Console.WriteLine(VDCresult.Theta);
                    Console.WriteLine("Option Rho: ");
                    Console.WriteLine(VDCresult.Rho);
                }
            }
            if (euro.IsAntithetic == true)
            {
                OptionResult AntitheticResult = new OptionResult();
                AntitheticResult.Price = euro.OptionPandStd(p, euro.IsCall, euro.IsAntithetic).VP;
                (AntitheticResult.Delta, AntitheticResult.Gamma) = euro.DeltaandGamma(p, euro.IsCall, euro.IsAntithetic);
                AntitheticResult.Theta = euro.ThetaValue(p, euro.IsCall, euro.IsAntithetic);
                AntitheticResult.Vega = euro.VegaValue(p, euro.IsCall, euro.IsAntithetic);
                AntitheticResult.Rho = euro.RhoValue(p, euro.IsCall, euro.IsAntithetic);
                AntitheticResult.StdErrorNorm = euro.StandardError(p, euro.IsCall, false);
                AntitheticResult.StdErrorReduce = euro.StandardError(p, euro.IsCall, true); 
                Console.WriteLine("Option Price: ");
                Console.WriteLine(AntitheticResult.Price);
                Console.WriteLine("Option Delta: ");
                Console.WriteLine(AntitheticResult.Delta);
                Console.WriteLine("Option Gamma: ");
                Console.WriteLine(AntitheticResult.Gamma);
                Console.WriteLine("Option Vega: ");
                Console.WriteLine(AntitheticResult.Vega);
                Console.WriteLine("Option Theta: ");
                Console.WriteLine(AntitheticResult.Theta);
                Console.WriteLine("Option Rho: ");
                Console.WriteLine(AntitheticResult.Rho);
                Console.WriteLine("Regular Standard Error");
                Console.WriteLine(AntitheticResult.StdErrorNorm);
                Console.WriteLine("Antithetic Standard Error");
                Console.WriteLine(AntitheticResult.StdErrorReduce);




            }
            else
            {
                OptionResult Result = new OptionResult();
                Result.Price = euro.OptionPandStd(p, euro.IsCall, euro.IsAntithetic).VP;
                (Result.Delta, Result.Gamma) = euro.DeltaandGamma(p, euro.IsCall, euro.IsAntithetic);
                Result.Theta = euro.ThetaValue(p, euro.IsCall, euro.IsAntithetic);
                Result.Vega = euro.VegaValue(p, euro.IsCall, euro.IsAntithetic);
                Result.Rho = euro.RhoValue(p, euro.IsCall, euro.IsAntithetic);
                Result.StdErrorNorm = euro.StandardError(p, euro.IsCall, false);
                Console.WriteLine("Option Price: ");
                Console.WriteLine(Result.Price);
                Console.WriteLine("Option Delta: ");
                Console.WriteLine(Result.Delta);
                Console.WriteLine("Option Gamma: ");
                Console.WriteLine(Result.Gamma);
                Console.WriteLine("Option Vega: ");
                Console.WriteLine(Result.Vega);
                Console.WriteLine("Option Theta: ");
                Console.WriteLine(Result.Theta);
                Console.WriteLine("Option Rho: ");
                Console.WriteLine(Result.Rho);
                Console.WriteLine("Regular Standard Error");
                Console.WriteLine(Result.StdErrorNorm);
            }

            
        }
    
        class SimulationParameters
        {
            public double S0{get; set;}
            public double r{get; set;}
            public int Steps{get; set;}
            public int Simulations{get; set;}
            public double Tenor {get; set;}
            public double Volatility{get; set;}

            public double Strike{get; set;}
        }
        class VDCparameters
        {
            public bool IsVDP{get; set;}
            public int Base{get; set;}
        }
        class OptionResult 
        {
            public double Price {get; set;}
            public double Delta{get; set;}
            public double Theta{get; set;}
            public double Gamma{get; set;}
            public double Vega{get; set;}
            public double Rho{get; set;}
            public double StdErrorNorm{get; set;}
            public double StdErrorReduce{get; set;}
        }
        class GaussianRandoms
        {
            
            public double [,] PopulatedNrands(int rows, int cols)
            {   
                double [,] NRands = new double [rows, cols];
                Random r_1 = new Random(8);
                Random r_2 = new Random(11);
                for(int i = 0; i <rows; i++)
                {
                    for (int j = 0; j< cols; j++)
                    {
                        double x1 = r_1.NextDouble();
                        double x2 = r_2.NextDouble();
                        double z1 = Math.Sqrt(-2 * Math.Log(x1))*Math.Cos(-2*Math.PI*x2);
                        double z2 = Math.Sqrt(-2 * Math.Log(x1))*Math.Sin(-2*Math.PI*x2);
                        NRands[i,j] = z1;
                    }

                }
                return NRands;
            }

            public double[,] ReverseNrands(int rows, int cols)
            {
                double [, ]NRands = PopulatedNrands(rows, cols);
                double [, ] NegNrands = new double [rows, cols];
                for(int i = 0; i <rows; i++)
                {
                    for (int j = 0; j< cols; j++)
                    {
                        NegNrands[i,j] = NRands[i,j] *-1;
                    }

                }
                return NegNrands;
                
            }
        }

        class SimulationResult
        {
            public double[,] SimulatedPaths {get; set;}
        }
        class MontecarloSimimulator
        {   
            public bool IsCall {get; set;}
            public bool IsAntithetic {get; set;}
            public static SimulationResult GeneratePaths(SimulationParameters p)
            {
                GaussianRandoms Nrands = new GaussianRandoms();
                double [, ]R1 = Nrands.PopulatedNrands(p.Simulations,p.Steps);

                SimulationResult results = new SimulationResult();
                results.SimulatedPaths = new double [p.Simulations,p.Steps];
                double dt = p.Tenor/p.Steps;
                
                for (int i = 0; i < p.Simulations; i++)
                {   
                    results.SimulatedPaths[i,0] = p.S0;
                    for (int j = 1; j < p.Steps; j++)
                    {
                
                        results.SimulatedPaths[i,j] = results.SimulatedPaths[i, j-1] * Math.Exp(((p.r - 0.5 * Math.Pow(p.Volatility,2))* dt) + (p.Volatility * Math.Sqrt(dt)*R1[i,j]));
        
                    }
                }

                return results;
            }

            public static SimulationResult AntetheticPaths(SimulationParameters p )
            {
                GaussianRandoms Nrands2 = new GaussianRandoms();
                double [,] R2 = Nrands2.ReverseNrands(p.Simulations,p.Steps);
                SimulationResult reverseresult = new SimulationResult();
                reverseresult.SimulatedPaths = new double [p.Simulations, p.Steps];
                double dt = p.Tenor/p.Steps;
                
                for (int i = 0; i < p.Simulations; i++)
                {   
                    reverseresult.SimulatedPaths[i,0] = p.S0;
                    for (int j = 1; j < p.Steps; j++)
                    {
                
                        reverseresult.SimulatedPaths[i,j] = reverseresult.SimulatedPaths[i, j-1] * Math.Exp(((p.r - 0.5 * Math.Pow(p.Volatility,2))* dt) + (p.Volatility * Math.Sqrt(dt)*R2[i,j]));
        
                    }
                }
                return reverseresult;
            }

            public (double VP, double std) OptionPandStd(SimulationParameters p, bool IsCall, bool IsAntithetic)
            {
                double VP= new double();
                double VNDC = new double();
                double std = new double();
                double var = new double();
                double total1 = new double();
                double total2 = new double();
                double discontfator = new double();
                discontfator = Math.Exp(-p.r*p.Tenor);
                SimulationResult SpathNorm = GeneratePaths(p);
                SimulationResult SpathAnti = AntetheticPaths(p);
                List<double> eachOption = new List<double>();
                List<double> reverseOption = new List<double>();
                 
                for (int i = 0; i < p.Simulations; i++)
                {
                    if (IsCall == true)
                    {
                        eachOption.Add(Math.Max(SpathNorm.SimulatedPaths[i, p.Steps - 1]-p.Strike, 0));
                    }
                
                    else if (IsCall == false)
                    {
                        eachOption.Add(Math.Max(p.Strike - SpathNorm.SimulatedPaths[i, p.Steps - 1], 0));
                    }

                } 
                total1 = eachOption.Sum();
                if (IsAntithetic == true)
                {
                    for (int i = 0; i < p.Simulations; i++)
                    {
                       if (IsCall == true)
                        {
                            reverseOption.Add(Math.Max(SpathAnti.SimulatedPaths[i, p.Steps - 1]-p.Strike, 0));
                        }
                
                        else if (IsCall == false)
                        {   
                            reverseOption.Add(Math.Max(p.Strike - SpathAnti.SimulatedPaths[i, p.Steps - 1], 0));
                        } 
                    }
                    total2 = reverseOption.Sum();
                    VNDC  = (total1 + total2)/(2*p.Simulations);
                    VP = discontfator*(total1 + total2)/(2*p.Simulations);
                    List<double> varlistOne = new List<double>();
                    for (int i = 0; i< p.Simulations; i++)
                    {
                        var = Math.Pow(eachOption[i] - VNDC,2) + Math.Pow(reverseOption[i] - VNDC, 2);
                        varlistOne.Add(var);
                    }
                    std = Math.Sqrt(varlistOne.Sum() / (p.Simulations * 2));
                    return (VP,std);
                }
                else
                {
                    VP = discontfator*total1/p.Simulations;
                    VNDC = total1/p.Simulations;
                    List<double> varlistTwo = new List<double>();
                    for (int i = 0; i < p.Simulations; i++)
                    {
                        var = Math.Pow(eachOption[i]- VNDC, 2);
                        varlistTwo.Add(var);
                    }
                    std = Math.Sqrt(varlistTwo.Sum() / p.Simulations);
                    return (VP,std);
                }    
            } 
            public double StandardError(SimulationParameters p, bool IsCall, bool IsAntithetic)
            {
                double SE = new double();
                double mu = OptionPandStd(p, IsCall, IsAntithetic).VP;
                double std = OptionPandStd(p, IsCall, IsAntithetic).std;
                if (IsAntithetic == true)
                {
                    SE = std/Math.Sqrt(p.Simulations * 2);
                }
                else 
                {
                    SE = std/Math.Sqrt(p.Simulations);
                }
                return SE;
            }

            public (double Delta, double Gamma) DeltaandGamma(SimulationParameters p, bool IsCall, bool IsAntithetic)
            {
                double deltaS = p.S0*0.1;
                double Delta = new double();
                double Gamma = new double();
                double VP = OptionPandStd(p, IsCall, IsAntithetic).VP;
                SimulationParameters p1 = new SimulationParameters();
                SimulationParameters p2 = new SimulationParameters();
                

                p1.S0 = p.S0 + deltaS;
                p1.r = p.r;
                p1.Steps = p.Steps;
                p1.Simulations = p.Simulations;
                p1.Tenor = p.Tenor;
                p1.Volatility= p.Volatility ;
                p2.S0 = p.S0 -deltaS;
                p2.r = p.r;
                p2.Steps = p.Steps;
                p2.Simulations = p.Simulations;
                p2.Tenor = p.Tenor;
                p2.Volatility = p.Volatility;

                double VP1 = OptionPandStd(p1, IsCall, IsAntithetic).VP;
                double VP2 = OptionPandStd(p2, IsCall, IsAntithetic).VP;
                Delta = (VP1 - VP2)/ (2 * deltaS);
                Gamma = (VP1 + VP2 - 2*VP)/ (Math.Pow(deltaS,2));
                if (IsCall == true)
                {
                    return(Delta, Gamma);
                }
                else
                {
                    return(Delta - 1, Gamma);
                }                 
            }
            public double VegaValue(SimulationParameters p, bool IsCall, bool IsAntithetic)
            {
                double deltaSig = p.Volatility*0.01;
                double Vega = new double();
                SimulationParameters p1 = new SimulationParameters();
                SimulationParameters p2 = new SimulationParameters();
                
                p1.S0 = p.S0;
                p1.r = p.r;
                p1.Steps = p.Steps;
                p1.Simulations = p.Simulations;
                p1.Tenor = p.Tenor;
                p1.Volatility= p.Volatility + deltaSig;
                p2.S0 = p.S0;
                p2.r = p.r;
                p2.Steps = p.Steps;
                p2.Simulations = p.Simulations;
                p2.Tenor = p.Tenor;
                p2.Volatility = p.Volatility - deltaSig;
                double VP1 = OptionPandStd(p1, IsCall, IsAntithetic).VP;
                double VP2 = OptionPandStd(p2, IsCall, IsAntithetic).VP;
                Vega = (VP1 - VP2) / (2 * deltaSig);
                return Vega;
            }
            public double ThetaValue(SimulationParameters p, bool IsCall, bool IsAntithetic)
            {
                double deltaT = p.Tenor * 0.1;
                double Theta = new double();
                double VP = OptionPandStd(p, IsCall, IsAntithetic).VP;
                SimulationParameters p1 = new SimulationParameters();
                p1.S0 = p.S0;
                p1.r = p.r;
                p1.Steps = p.Steps;
                p1.Simulations = p.Simulations;
                p1.Volatility= p.Volatility;
                p1.Tenor = p.Tenor + deltaT;
                double VP1 = OptionPandStd(p1, IsCall, IsAntithetic).VP;
                Theta = (VP1 - VP)/ (deltaT);
                return Theta;
            }
            public double RhoValue(SimulationParameters p, bool IsCall, bool IsAntithetic)
            {
                double deltar = p.r*0.1;
                double Rho = new double();
                SimulationParameters p1 = new SimulationParameters();
                SimulationParameters p2 = new SimulationParameters();
                
                
                p1.S0 = p.S0;
                p1.r= p.r + deltar;
                p1.Steps = p.Steps;
                p1.Simulations = p.Simulations;
                p1.Tenor = p.Tenor;
                p1.Volatility= p.Volatility;
                p2.S0 = p.S0;
                p2.r = p.r - deltar;
                p2.Steps = p.Steps;
                p2.Simulations = p.Simulations;
                p2.Tenor = p.Tenor;
                double VP1 = OptionPandStd(p1, IsCall, IsAntithetic).VP;
                double VP2 = OptionPandStd(p2, IsCall, IsAntithetic).VP;
                Rho = (VP1 - VP2)/(2 * deltar);
                return Rho;

            }


        }
        class VDCfunc
        {
            public static List<double> VDC(int Base, int n)
            {
                List<double> list = new List<double>();
                long p = 0, q = 0;
                double r = 0;
                int i = 0;
                while( q <= Base)
                {
                    i++;
                    p = n%Base;
                    q = (n-p) / Base;
                    r += Math.Pow(Base, -i);
                }
                list.Add(r);
                return list;
            }
            public static List<double> VDCRands (int Base, int Simulations)
            {
                List<double> NRands = new List<double>();
                List<double> list1 = VDC(Base, Simulations);
                List<double> list2 = VDC(Base + 2, Simulations);
                for (int i = 0; i <  list1.Count() - 1; i++)
                {
                    NRands[i] = Math.Sqrt(-2 * Math.Log(list2[i])) * Math.Cos(-2*Math.PI*list1[i]);
                }
                return NRands;
            }
            public static List<double> SSimulation(int Base, SimulationParameters p)
            {
                List<double> SPs = new List<double>();
                List<double> rands  = VDCRands(Base,p.Simulations);
                for (int i = 0; i < rands.Count() - 1; i++)
                {
                    SPs[i] = p.S0 * Math.Exp(((p.r - 0.5 * Math.Pow(p.Volatility,2))* p.Tenor) + (p.Volatility * Math.Sqrt(p.Tenor)*rands[i]));

                }
                return SPs;
            }
            public double VDCOptionPrice (int Base, SimulationParameters p, bool IsCall)
            {
                List<double> SPS = SSimulation(Base,p);
                List<double> eachOption = new List<double>();
                double discontfator = new double();
                discontfator = Math.Exp(-p.r*p.Tenor);
                double VP= new double();
                for (int i = 0; i < SPS.Count() - 1; i++)
                {
                    if (IsCall == true)
                    {
                        eachOption.Add(Math.Max(SPS[i]-p.Strike, 0));
                    }
                
                    else if (IsCall == false)
                    {
                        eachOption.Add(Math.Max(p.Strike - SPS[i], 0));
                    }

                }
                VP = discontfator*(eachOption.Sum() / p.Simulations);
                return VP;
            }
            public (double Delta, double Gamma) DeltaandGamma(int Base, SimulationParameters p, bool IsCall)
            {
                double deltaS = p.S0*0.01;
                double Delta = new double();
                double Gamma = new double();
                double VP = VDCOptionPrice(Base, p, IsCall);
                SimulationParameters p1 = new SimulationParameters();
                SimulationParameters p2 = new SimulationParameters();
                p1.S0 = p.S0 + deltaS;
                p1.r = p.r;
                p1.Steps = p.Steps;
                p1.Simulations = p.Simulations;
                p1.Tenor = p.Tenor;
                p1.Volatility= p.Volatility ;
                p2.S0 = p.S0 -deltaS;
                p2.r = p.r;
                p2.Steps = p.Steps;
                p2.Simulations = p.Simulations;
                p2.Tenor = p.Tenor;
                p2.Volatility = p.Volatility;
                double VP1 = VDCOptionPrice(Base, p1, IsCall);
                double VP2 = VDCOptionPrice(Base, p2, IsCall);
                Delta = (VP1 - VP2)/ (2 * deltaS);
                Gamma = (VP1 + VP2 - 2*VP)/ (Math.Pow(deltaS,2));
                if (IsCall == true)
                {
                    return(Delta, Gamma);
                }
                else
                {
                    return(Delta - 1, Gamma);
                }                 
            }
            public double VegaValue(int Base, SimulationParameters p, bool IsCall)
            {
                double deltaSig = p.Volatility *0.1;
                double Vega = new double();
                SimulationParameters p1 = new SimulationParameters();
                SimulationParameters p2 = new SimulationParameters();
                p1.S0 = p.S0;
                p1.r = p.r;
                p1.Steps = p.Steps;
                p1.Simulations = p.Simulations;
                p1.Tenor = p.Tenor;
                p1.Volatility= p.Volatility + deltaSig;
                p2.S0 = p.S0;
                p2.r = p.r;
                p2.Steps = p.Steps;
                p2.Simulations = p.Simulations;
                p2.Tenor = p.Tenor;
                p2.Volatility = p.Volatility - deltaSig;
                double VP1 = VDCOptionPrice(Base, p1, IsCall);
                double VP2 = VDCOptionPrice(Base, p2, IsCall);
                Vega = (VP1 - VP2) / (2 * deltaSig);
                return Vega;
            }
            public double ThetaValue(int Base, SimulationParameters p, bool IsCall)
            {
                double deltaT = p.Tenor/p.Strike;
                double Theta = new double();
                double VP = VDCOptionPrice(Base, p, IsCall);
                SimulationParameters p1 = new SimulationParameters();
                p1.S0 = p.S0;
                p1.r = p.r;
                p1.Steps = p.Steps;
                p1.Simulations = p.Simulations;
                p1.Volatility= p.Volatility;
                p1.Tenor = p.Tenor + deltaT;
                double VP1 = VDCOptionPrice(Base, p1, IsCall);
                Theta = (VP1 - VP)/ (deltaT);
                return Theta;
            }
            public double RhoValue(int Base, SimulationParameters p, bool IsCall)
            {
                double deltar = p.r * 0.1;
                double Rho = new double();
                SimulationParameters p1 = new SimulationParameters();
                SimulationParameters p2 = new SimulationParameters();
                p1.S0 = p.S0;
                p1.r= p.r + deltar;
                p1.Steps = p.Steps;
                p1.Simulations = p.Simulations;
                p1.Tenor = p.Tenor;
                p1.Volatility= p.Volatility;
                p2.S0 = p.S0;
                p2.r = p.r - deltar;
                p2.Steps = p.Steps;
                p2.Simulations = p.Simulations;
                p2.Tenor = p.Tenor;
                double VP1 = VDCOptionPrice(Base, p1, IsCall);
                double VP2 = VDCOptionPrice(Base, p2, IsCall);
                Rho = (VP1 - VP2)/(2 * deltar);
                return Rho; 
            }

        }

    }
}
/// Colab with Hu and Alex