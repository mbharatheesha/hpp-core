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

#ifndef HPP_CORE_STEERING_ILQR
# define HPP_CORE_STEERING_ILQR

# include <iostream>
# include <hpp/core/fwd.hh>
# include <eigen3/Eigen/Core>
# include <hpp/core/steering-statespace.hh>

namespace hpp {
    namespace core {
        class HPP_CORE_DLLAPI SteeringILQR : public SteeringStateSpace {
            protected:
                int numIter_;
                double DT_;
                MatrixXd Qstate_;
                MatrixXd Qfinal_;
                MatrixXd Rcontrol_;

            public:
                /// Constructor
                SteeringILQR () : SteeringStateSpace (), numIter_ (20), DT_ (0.005),
                Qstate_ (), Qfinal_ (), Rcontrol_ ()
            {
            }
                /// Set the time discretization parameter
                void setTimeDiscretization (double);
                
                /// Set the weighting matrices for computing cost
                void setWeightingMatrices (VectorXd, VectorXd, VectorXd);
                
                /// Cost Function
                double computeCost (MatrixXd, MatrixXd);
                
                /// Forward Pass
                MatrixXd forwardPass (MatrixXd (*fdyn) (VectorXd, VectorXd));

                /// Linearize nonlinear dynamics around a given trajectory and input
                void linearizeSystem (MatrixXd, MatrixXd);

                /// Forward Finite Difference Linearization
                MatrixXd linearizeFDJacobian (VectorXd, VectorXd, MatrixXd&);

        };
    }

}
#endif
