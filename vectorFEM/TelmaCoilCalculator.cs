using ReactiveUI;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using TelmaConsole;
using TelmaQuasarCommon.Core.ContextInterfaces;
using TelmaQuasarCommon.ProblemBase;
using TelmaQuasarCommon.Problems;
using TelmaQuasarCommon.ProjectSystem;

namespace Core
{
    public class TelmaCoilCalculator : IDisposable
    {
        static TelmaCoilCalculator()
        {
            LogManager.Init();
            ContextManager.Context.Config.Load();

            //App.splash = new GUI.SplashScreen();
            //App.splash.Show();
            //            AppDomain.CurrentDomain.AssemblyResolve += AssemblyResolver;
            AppDomain.CurrentDomain.UnhandledException += ExceptionHandler;
            System.Globalization.CultureInfo? culture =
                Thread.CurrentThread.CurrentCulture.Clone() as System.Globalization.CultureInfo;
            //            culture.NumberFormat.NumberDecimalSeparator = ".";
            culture!.NumberFormat = System.Globalization.CultureInfo.InvariantCulture.NumberFormat;
            Thread.CurrentThread.CurrentCulture = culture;
            System.Globalization.CultureInfo.DefaultThreadCurrentCulture = culture;
            System.Globalization.CultureInfo.DefaultThreadCurrentUICulture = culture;

        }

        public async Task Load(string path, string problemName)
        {
            var project = await Project.Load(path) ?? throw new Exception($"Can't read file {path}");
            problem = await project.BuildProblem(problemName) as ElectroMagneticVectorScalarProblem ?? throw new ArgumentNullException();
        }

        ElectroMagneticVectorScalarProblem problem;

        public Core.Vector3D Hext(Core.Vector3D point, double t)
        {
            TelmaQuasarCommon.Vector3D p = new(point.AsArray());

            var H = problem.Coil.CoilCalculatorH.Value(p);

            H *= problem.SourceTimeMultiplier.Value(default, t);

            return new Core.Vector3D(H.AsArray()) ;
        }
        
        static void ExceptionHandler(object sender, UnhandledExceptionEventArgs e)
        {
            if (e.ExceptionObject is Exception)
            {

                var ex = e.ExceptionObject as Exception;
                ContextManager.Context.Logger.Log(LogLevel.Critical, ex);
                Console.WriteLine(string.Format("Exception: {0}", ex!.Message), "Unhandled Exception in Telma! :(");
            }
            else
            {
                //if(ProgressScreenManager.IsActive)
                //            ProgressScreenManager.ProcessException(new Exception("Unknown exception"));

                ContextManager.Context.Logger.Log(LogLevel.Critical, "Unknown exception (maybe in unmanaged code).");
                Console.WriteLine("Unknown exception (maybe in unmanaged code)." + "Unhandled Exception in Telma! :(");
            }
        }

        private bool disposedValue;

        protected virtual void Dispose(bool disposing)
        {
            if (!disposedValue)
            {
                if (disposing)
                {
                    problem.Dispose();
                }

                // TODO: освободить неуправляемые ресурсы (неуправляемые объекты) и переопределить метод завершения
                // TODO: установить значение NULL для больших полей
                disposedValue = true;
            }
        }

        // // TODO: переопределить метод завершения, только если "Dispose(bool disposing)" содержит код для освобождения неуправляемых ресурсов
        // ~TelmaCoilCalculator()
        // {
        //     // Не изменяйте этот код. Разместите код очистки в методе "Dispose(bool disposing)".
        //     Dispose(disposing: false);
        // }

        void IDisposable.Dispose()
        {
            // Не изменяйте этот код. Разместите код очистки в методе "Dispose(bool disposing)".
            Dispose(disposing: true);
            GC.SuppressFinalize(this);
        }
    }
}
