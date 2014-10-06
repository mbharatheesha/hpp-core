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
        void SteeringILQR::setILQRParams (int nIter, int numCtl, int uLen, double dt)
        {
            numIter_ = nIter;
            numControl_ = numCtl;
            controlLen_ = uLen;
            linMatrices_.resize (uLen); 
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
        double SteeringILQR::computeCost ()
        {
            double costFinal = (stateTraj_->col (controlLen_-1) - finalState_).transpose () * Qfinal_
                * (stateTraj_->col (controlLen_-1) - finalState_);
            double costTraj = 0;

            for (int i=0; i < stateTraj_->cols () - 1; i++)
            {
                costTraj += ((stateTraj_->col (i) - finalState_).transpose () * Qstate_ * (stateTraj_->col (i) - finalState_)) +
                    (controlSeq_->col (i).transpose () * Rcontrol_ * controlSeq_->col (i));
            }

            return (0.5 * (costTraj + costFinal));
        }

        /// Forward Pass
        void SteeringILQR::forwardPass (const MatrixPtr_t& ctl)
        {
            VectorXd timeSeq (VectorXd::LinSpaced (controlLen_, 0.0, (DT_*(controlLen_-1))));
            *stateTraj_ = sysDyn_->simulateDynamics (timeSeq, initState_, *ctl);
        }

        /// Linearize the system around the trajectory
        void SteeringILQR::linearizeSystem ()
        {
            for (int i = 0; i < stateTraj_->cols (); i++)
            {
                linearizeFDJacobian (stateTraj_->col (i), controlSeq_->col (i), linMatrices_.at (i));
            }
        }

        /// Forward Finite Difference Linearization
        void SteeringILQR::linearizeFDJacobian (VectorXd state, VectorXd control,
                MatrixXd& lMat)
        {
            lMat.resize (state.rows (), state.rows () + control.rows ());
            int idxState;
            MatrixXd discreteIdentity (MatrixXd::Identity (state.rows (), state.rows ()));

            VectorXd nomFnVal = sysDyn_-> computeStateDerivative (0.0, state, control);
            for (idxState = 0; idxState < state.size (); idxState++)
            {
                VectorXd pertState = VectorXd::Zero (state.size ());

                pertState (idxState) = HPP_CORE_FDSTEP_SIZE;


                VectorXd pertFnVal = sysDyn_-> computeStateDerivative (0.0, state + pertState,
                        control);

                lMat.col (idxState) = discreteIdentity.col (idxState) +
                    ((pertFnVal - nomFnVal) * HPP_CORE_FDSTEP_SIZE_INV * DT_ );
            }
            
            for (int idxCtl = 0; idxCtl < control.size (); idxCtl++)
            {

                VectorXd pertControl = VectorXd::Zero (control.size ());

                pertControl (idxCtl) = HPP_CORE_FDSTEP_SIZE;

                VectorXd pertFnVal = sysDyn_-> computeStateDerivative (0.0, state,
                        control + pertControl);

                lMat.col (idxState + idxCtl) = (pertFnVal - nomFnVal) * HPP_CORE_FDSTEP_SIZE_INV * DT_;
            }
        }

        MatrixXd* SteeringILQR::deltaOptControl ()
        {
            MatrixXd Alin; 
            MatrixXd Blin;
            MatrixXd commonTerm2;
            Alin = (linMatrices_.at (controlLen_-1)).leftCols (initState_.size ());
            Blin = (linMatrices_.at (controlLen_-1)).rightCols (numControl);
            MatrixXd Snext = Qfinal_;
            MatrixXd valFun = MatrixXd::Zero (initState_.size (), controlLen_);

            valFun.col (controlLen_-1) = Snext *
                (stateTraj_->col (controlLen_-1) - finalState_);
            commonTerm2 = Blin * Rcontrol_.inverse () * Blin.transpose ();
            MatrixXd dStateNext = (MatrixXd::Identity (initState_.size (), initState_.size ()) + 
                    ( commonTerm2 * Snext).inverse ()) *
                ((Alin * (*dStateTraj_) - commonTerm2 * valFun - Blin * (*controlSeq_));

            VectorXd valNext = valFun.col (controlLen_-1);

            MatrixXd K = MatrixXd::Zero (numControl_, initState_.size ());
            MatrixXd Kv = MatrixXd::Zero (numControl_, initState_.size ());
            MatrixXd Ku = MatrixXd::Zero (numControl_, numControl_);
            MatrixXd commonTerm;
            MatrixXd Snew;
            MatrixXd closedLoopDyn;
            for (int i = controlLen_-2; i >= 0; i--)
            {
                Alin = (linMatrices_.at (i)).leftCols (initState_.size ());
                Blin = (linMatrices_.at (i)).rightCols (numControl_);

                commonTerm = ((Blin.transpose () * Snext * Blin) + Rcontrol_ ).inverse (); 

                K = commonTerm * Blin.transpose () * Snext * Alin;


                Kv = commonTerm * Blin.transpose ();

                Ku = commonTerm * Rcontrol_;

                closedLoopDyn = Alin - (Blin * K);

                Snew = Alin.transpose () * Snext * closedLoopDyn + Qstate_;

                valFun.col (i) = closedLoopDyn.transpose () * valNext -
                        (K.transpose () * Rcontrol_ * controlSeq_->col (i)) +
                        (Qstate_ * stateTraj_->col (i));

                dControlSeq_->col (i) =- ((K * dStateTraj_->col (i)) + (Kv * valNext) +
                        (Ku * controlSeq_->col (i)));
                

                valNext = valFun.col(i);
                Snext = Snew;
            }

            return dControlSeq_;
        }

        MatrixPtr_t SteeringILQR::steerState (VectorXd init, VectorXd final)
        {
            initState_ = init;
            finalState_ = final;
            std::vector <double> costVec (numIter_);

            initializeMemory ();

            /// Do the forward pass with initial control input and compute the cost
            forwardPass (controlSeq_);
            double curCost = computeCost ();

            /// Iteratively find the control improvement along the trajectory
            for (int i=0; i < numIter_; i++)
            {
                /// Linearize the system around the forward pass trajectory
                linearizeSystem ();

                /// Compute the control improvement via backward pass
                dControlSeq_ = deltaOptControl ();

                *controlSeq_ = *controlSeq_ + (0.25 * (*dControlSeq_));
                forwardPass (controlSeq_);
                costVec.at (i) = computeCost ();

                curCost = costVec.at (i);

                std::cout << curCost << std::endl;
            }

            std::ofstream file1;
            file1.open ("control_profile.dat");

            file1 << *controlSeq_;
            file1.close ();
            
            std::cout << "ILQR control computed successfully!" << std::endl;
            delete dStateTraj_;
            delete dControlSeq_;
            return stateTraj_;
        } 
    }
}
