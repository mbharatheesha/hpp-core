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
#include <hpp/core/simple-pendulum.hh>

namespace hpp {
    namespace core {
        //Dynamics of Simple Pendulum
        void SimplePendulum::setProblemDimension (int n)
        {
            nDOF_ = n;
        }
        void SimplePendulum::setParameters (void)
        {
            length_ = 1.0; // m
            mass_ = 1.0;   // kg
            friction_ = 0.1; // kg/s
        }
        vectorOut_t SimplePendulum::computeStateDerivative
            (double time, vectorIn_t stateVector, vectorIn_t control)
            {
                VectorXd stateDot (stateVector.size ());

                dynAutonomous_.resize (stateVector.size ());

                dynForced_.resize (stateVector.size ());

                dynAutonomous_ (0) = stateVector (1);
                dynAutonomous_ (1) = -(GRAVITY / length_) * sin (stateVector (0)) -
                    ((friction_ * stateVector (1))/(mass_ * length_*length_));

                dynForced_ << 0,1;

                stateDot = dynAutonomous_ + (dynForced_ * control);

                return stateDot;
            }
        vectorOut_t SimplePendulum::getControl (vectorIn_t stateVector)
        {
            VectorXd control;

            control.resize (nDOF_);

            control(0) = 0.5 * 0;

            return control;
        }

        /// Used for feedback control
        matrixOut_t SimplePendulum::simulateDynamics (vectorIn_t tVec,

                vectorIn_t initState)
        {
            int trajLen = tVec.size();
            MatrixXd stateTraj;
            double stepSize;

            stateTraj.resize (initState.size (), trajLen);

            stateTraj.col (0) = initState;


            for (int i = 1; i < trajLen ; i++)
            {
                stepSize = tVec (i) - tVec (i-1);
                VectorXd newState = integrateRK4 (tVec (i-1), stateTraj.col (i-1),
                        getControl (stateTraj.col(i-1)), stepSize);
                stateTraj.col (i) = newState;
            }

            return stateTraj;
        }

        /// Used for feedforward control
        matrixOut_t SimplePendulum::simulateDynamics (vectorIn_t tVec,
                vectorIn_t initState, matrixIn_t control)
        {
            int trajLen = tVec.size();
            MatrixXd stateTraj;
            double stepSize;

            stateTraj.resize (initState.size (), trajLen);

            stateTraj.col (0) = initState;


            for (int i = 1; i < trajLen ; i++)
            {
                stepSize = tVec (i) - tVec (i-1);
                VectorXd newState = integrateEuler (tVec (i-1), stateTraj.col (i-1),
                        control.col (i), stepSize);
                stateTraj.col (i) = newState;
            }

            return stateTraj;
        }
        
        vectorOut_t SimplePendulum::integrateRK4 (double t, vectorIn_t state, vectorIn_t u, double h)
        {
            VectorXd st1 = computeStateDerivative (t, state, u);
            VectorXd st2 = computeStateDerivative (t + (0.5 * h), state + (0.5 * h * st1), u);
            VectorXd st3 = computeStateDerivative (t + (0.5 * h), state + (0.5 * h * st2),  u);
            VectorXd st4 = computeStateDerivative (t + h, state + (h * st3), u);

            VectorXd stNew = state + ((1/6.0) * h * (st1 + 2.0*st2 + 2.0*st3 + st4));

            return (stNew);
        }
        
        vectorOut_t SimplePendulum::integrateEuler (double t, vectorIn_t state, vectorIn_t u, double h)
        {
            VectorXd st = computeStateDerivative (t, state, u);
            VectorXd stNew = state + h*st;
            return (stNew);
        }
    }
}
