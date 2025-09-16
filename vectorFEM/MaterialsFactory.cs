using Core;
using FEM;
using Microsoft.CodeAnalysis.CSharp.Scripting;
using Microsoft.CodeAnalysis.Scripting;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Reflection;
using System.Text;
using System.Threading.Tasks;

namespace Core2
{
    public static class MaterialsFactory
    {
        static ScriptOptions options = ScriptOptions.Default.AddReferences(typeof(Program).Assembly);

        public static IDictionary<string, IMaterial> CreateMaterials(string path)
        {
            if (!Path.Exists(path)) throw new ArgumentException("The path is incorrect");

            using var reader = new StreamReader(Path.Combine(path, "materialsFunctions.txt"));

            int n = int.Parse(reader.ReadLine());
            IDictionary<string, IMaterial> materials = new Dictionary<string, IMaterial>(n);

            for (int i = 0; i < n; ++i)
            {
                string[] values = reader.ReadLine().Split(' ', StringSplitOptions.RemoveEmptyEntries);
                int numOfFunc = int.Parse(values[2]);
                MaterialType type = (MaterialType)int.Parse(values[1]);
                Material material = new(type);

                for (int f = 0; f < numOfFunc; ++f)
                {
                    string funcName = reader.ReadLine().Trim();
                    string funcBody = reader.ReadLine().Trim();

                    AddFunctionInMaterial(material, funcName, funcBody);
                }

                materials.Add(values[0], material);
            }

            return materials;
        }

        static void AddFunctionInMaterial(Material material, string funcName, string funcBody)
        {
            switch (funcName)
            {
                case "Lambda":
                    material.Lambda = Compile<Func<Vector3D, double>>(funcBody, "p => ");
                    break;
                case "Sigma":
                    material.Sigma = Compile<Func<Vector3D, double>>(funcBody, "p => ");
                    break;
                case "Epsilon":
                    material.Epsilon = Compile<Func<Vector3D, double>>(funcBody, "p => ");
                    break;
                case "Mu":
                    material.Mu = Compile<Func<Vector3D, double>>(funcBody, "p => ");
                    break;
                case "Betta":
                    material.Betta = Compile<Func<Vector3D, double>>(funcBody, "p => ");
                    break;
                case "UBetta":
                    material.UBetta = Compile<Func<Vector3D, double, double>>(funcBody, "(p, t) => ");
                    break;
                case "Theta":
                    material.Theta = Compile<Func<Vector3D, double, double>>(funcBody, "(p, t) => ");
                    break;
                case "Ug":
                    material.Ug = Compile<Func<Vector3D, double, double>>(funcBody, "(p, t) => ");
                    break;
                case "F":
                    material.F = Compile<Func<Vector3D, double, double>>(funcBody, "(p, t) => ");
                    break;
                case "Htheta":
                    material.Htheta = Compile<Func<Vector3D, double, Vector3D>>(funcBody, "(p, t) => "); 
                    break;
                case "Hext":
                    material.Hext = Compile<Func<Vector3D, double, Vector3D>>(funcBody, "(p, t) => ");
                    break;
                case "Ag":
                    material.Ag = Compile<Func<Vector3D, double, Vector3D>>(funcBody, "(p, t) => ");
                    break;
                case "Fv":
                    material.Fv = Compile<Func<Vector3D, double, Vector3D>>(funcBody, "(p, t) => ");
                    break;

                default:
                    throw new ArgumentException("Incorrect function name");
            }
        }

        static T Compile<T>(string body, string prefix) => 
            CSharpScript.EvaluateAsync<T>(prefix + body, options).GetAwaiter().GetResult();
    }
}
