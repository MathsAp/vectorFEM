using FEM;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Core
{
    public class VectorScalarParabolicProblem : IProblem
    {
        public IDictionary<string, IMaterial> Materials => throw new NotImplementedException();

        public void Prepare()
        {
            throw new NotImplementedException();
        }

        public void Solve(ISolution result)
        {
            throw new NotImplementedException();
        }
    }
}
