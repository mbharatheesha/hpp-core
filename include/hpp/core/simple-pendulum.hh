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
        class HPP_CORE_DLLAPI SimplePendulum : public SystemDynamics, public IntegrateDynamics {
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
                vectorOut_t computeStateDerivative (double time, vectorIn_t stateVector, vectorIn_t control);
                vectorOut_t getControl (vectorIn_t stateVector);
                matrixOut_t simulateDynamics (vectorIn_t timeVector, vectorIn_t initState);
                matrixOut_t simulateDynamics (vectorIn_t timeVector, vectorIn_t initState, matrixIn_t control);
                vectorOut_t integrateRK4 (double time, vectorIn_t state, vectorIn_t control, double timeStep);
                vectorOut_t integrateEuler (double time, vectorIn_t state, vectorIn_t control, double timeStep);
        };
    }
}
#endif
