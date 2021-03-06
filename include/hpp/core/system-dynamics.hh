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

#ifndef HPP_CORE_SYSTEMDYNAMICS_HH
# define HPP_CORE_SYSTEMDYNAMICS_HH

# include <iostream>
# include <hpp/core/fwd.hh>
# include <eigen3/Eigen/Core>
#include <hpp/core/config.hh>

using namespace Eigen;
namespace hpp {
    namespace core {
        class HPP_CORE_DLLAPI SystemDynamics {
            protected:
                int nDOF_;
                VectorXd dynAutonomous_;
                VectorXd dynForced_;
                VectorXd state_;

            public:
                // Abstract class for system dynamics
                virtual void setProblemDimension (int) = 0;
                virtual vectorOut_t computeStateDerivative (
                        double time, vectorIn_t state, vectorIn_t control) = 0;
                virtual matrixOut_t simulateDynamics (vectorIn_t timeVec,
                        vectorIn_t initState) = 0;
                virtual matrixOut_t simulateDynamics (vectorIn_t timeVec,
                        vectorIn_t initState, matrixIn_t controlFeedFwd) = 0;
        };
    }
}
#endif
