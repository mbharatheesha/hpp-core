//
// Copyright (c) 2014 CNRS
// Authors: Mukunda Bharatheesha 
//
// This file is part of hpp-core
// hpp-core is free software: you can redistribute it
// and/or modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation, either version
// 3 of the License, or (at your option) any later version.
//
// hpp-core is distributed in the hope that it will be
// useful, but WITHOUT ANY WARRANTY; without even the implied warranty
// of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Lesser Public License for more details.  You should have
// received a copy of the GNU Lesser General Public License along with
// hpp-core  If not, see
// <http://www.gnu.org/licenses/>.

#ifndef HPP_CORE_SIMPLEPENDULUM_HH
# define HPP_CORE_SIMPLEPENDULUM_HH

# include <iostream>
# include <hpp/core/fwd.hh>
# include <hpp/core/system-dynamics.hh>
# include <hpp/core/integrate-dynamics.hh>
# include <eigen3/Eigen/Core>

#define GRAVITY 9.81
namespace hpp {
    namespace core {
        class HPP_CORE_DLLAPI SimplePendulum : public SystemDynamics, protected IntegrateDynamics {
            protected:
                double length_;
                double mass_;
                double friction_;
            public:
                /// Constructor
                SimplePendulum () : SystemDynamics (), IntegrateDynamics (), length_ (1), mass_ (1), friction_ (0.01)
            {
            }
                // Abstract class for system dynamics
                void setProblemDimension (int n);
                void setParameters (void);
                VectorXd computeStateDerivative (double time, VectorXd stateVector, VectorXd control);
                VectorXd getControl (VectorXd stateVector);
                MatrixXd simulateDynamics (VectorXd timeVector, VectorXd initState);
                MatrixXd simulateDynamics (VectorXd timeVector, VectorXd initState, MatrixXd control);
                VectorXd integrateRK4 (double time, VectorXd state, VectorXd control, double timeStep);
                VectorXd integrateEuler (double time, VectorXd state, VectorXd control, double timeStep);
        };
    }
}
#endif
