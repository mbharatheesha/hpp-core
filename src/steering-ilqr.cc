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

#include <hpp/util/debug.hh>
#include <hpp/core/steering-ilqr.hh>

namespace hpp {
    namespace core {
        /// Set the time discretization parameter
        void SteeringILQR::setTimeDiscretization (double dt)
        {
            DT_ = dt;
        }
        /// Set the weighting matrices for computing cost
        void SteeringILQR::setWeightingMatrices (VectorXd stateWeights, VectorXd finalStateWeights,
                VectorXd controlWeights)
        {
            Qstate_ = stateWeights.asDiagonal ();
            Qfinal_ = finalStateWeights.asDiagonal ();
            Rcontrol_ = controlWeights.asDiagonal ();
        }

        /// Cost Function
        double SteeringILQR::computeCost (MatrixXd stateTraj, MatrixXd control)
        {
            double costFinal = finalState_.transpose () * Qfinal_ * finalState_;
            double costTraj = 0;

            for (int i=0; i < stateTraj.cols (); i++)
            {
                costTraj = costTraj + 
                    (stateTraj.col (i).transpose () * Qstate_ * stateTraj.col (i)) +
                    (control.col (i).transpose () * Rcontrol_ * control.col (i));
            }

            return (0.5 * (costTraj + costFinal));
        }

        /// Forward Pass
        MatrixXd SteeringILQR::forwardPass (MatrixXd (*fdyn) (VectorXd, VectorXd))
        {
            return (*fdyn) (timeVec_, initState_);
        }

        void SteeringILQR::linearizeSystem (MatrixXd traj, MatrixXd control)
        {
            //Create an array of matrices to store the linearizations around
            //the trajectory
            MatrixXd *linMatrices = new MatrixXd [traj.cols ()];

            for (int i = 0; i < traj.cols (); i++)
            {
                MatrixXd &linMat = linMatrices [i];
                linearizeFDJacobian (traj.col (i), control.col (i), linMat);
            }
        }

        /// Forward Finite Difference Linearization
        MatrixXd SteeringILQR::linearizeFDJacobian (VectorXd state, VectorXd control,
                MatrixXd& lMat)
        {
        }
    }
}
        
