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
            std::cout <<"It looks okay till here!"<<std::endl;
        }

        /// Cost Function
        double SteeringILQR::computeCost ()
        {
            double costFinal = (stateTraj_->col (controlLen_-1) - finalState_).transpose () * Qfinal_
                * (stateTraj_->col (controlLen_-1) - finalState_);
            double costTraj = 0;
            std::cout << "state trajectory: " << stateTraj_->col (controlLen_-1) << std::endl;

            for (int i=0; i < stateTraj_->cols () - 1; i++)
            {
                costTraj = costTraj + 
                    (stateTraj_->col (i).transpose () * Qstate_ * stateTraj_->col (i)) +
                    (controlSeq_->col (i).transpose () * Rcontrol_ * controlSeq_->col (i));
            }

            std::cout <<"Cost computed successfully: "<< costTraj << std::endl;
            return (0.5 * (costTraj + costFinal));
        }

        /// Forward Pass
        void SteeringILQR::forwardPass (const MatrixPtr_t& ctl)
        {
            VectorXd timeSeq (VectorXd::LinSpaced (controlLen_,0.0,controlLen_-1));
            *stateTraj_ = sysDyn_->simulateDynamics (timeSeq, initState_, *ctl);
            *dStateTraj_ = MatrixXd::Zero (2,100);
        }

        void SteeringILQR::linearizeSystem ()
        {
            for (int i = 0; i < stateTraj_->cols () - 1; i++)
            {
                linearizeFDJacobian (stateTraj_->col (i), controlSeq_->col (i), linMatrices_.at (i));
            }
            std::cout << "Linearizations finished successfully" << std::endl;
        }

        /// Forward Finite Difference Linearization
        void SteeringILQR::linearizeFDJacobian (VectorXd state, VectorXd control,
                MatrixXd& lMat)
        {
            lMat.resize (state.rows (), state.rows () + control.rows ());
            int idxState;
            MatrixXd discreteIdentity (MatrixXd::Identity (state.rows (), state.rows ()));

            for (idxState = 0; idxState < state.size (); idxState++)
            {
                VectorXd nomFnVal = sysDyn_-> computeStateDerivative (0.0, state, control);

                VectorXd pertState = VectorXd::Zero (state.size ());

                pertState (idxState) = HPP_CORE_FDSTEP_SIZE;


                VectorXd pertFnVal = sysDyn_-> computeStateDerivative (0.0, state + pertState,
                        control);

                lMat.col (idxState) = discreteIdentity.col (idxState) +
                    ((pertFnVal - nomFnVal) * HPP_CORE_FDSTEP_SIZE_INV * DT_);
            }
            for (int idxCtl = 0; idxCtl < control.size (); idxCtl++)
            {
                VectorXd nomFnVal = sysDyn_-> computeStateDerivative (0.0, state, control);

                VectorXd pertControl = VectorXd::Zero (control.size ());

                pertControl (idxCtl) = HPP_CORE_FDSTEP_SIZE;

                VectorXd pertFnVal = sysDyn_-> computeStateDerivative (0.0, state,
                        control + pertControl);

                lMat.col (idxState + idxCtl) = (pertFnVal - nomFnVal) * HPP_CORE_FDSTEP_SIZE_INV
                    * DT_;
            }
        }

        MatrixPtr_t SteeringILQR::deltaOptControl ()
        {
            MatrixXd Snext = Qfinal_;
            MatrixXd valFun = MatrixXd::Zero (initState_.size (), controlLen_);

            valFun.col (controlLen_-1) = Snext *
                (stateTraj_->col (controlLen_-1) - finalState_);

            VectorXd valNext = valFun.col (controlLen_-1);

            MatrixXd K = MatrixXd::Zero (numControl_, initState_.size ());
            MatrixXd Kv = MatrixXd::Zero (numControl_, initState_.size ());
            MatrixXd Ku = MatrixXd::Zero (numControl_, numControl_);
            MatrixXd Alin; 
            MatrixXd Blin; 
            for (int i = controlLen_-2; i > 0; i--)
            {
                Alin = (linMatrices_.at (i)).leftCols (initState_.size ());
                Blin = (linMatrices_.at (i)).rightCols (numControl_);



                K = ((Blin.transpose () * Snext * Blin) + Rcontrol_ ).inverse () *
                    Blin.transpose () * Snext * Alin;


                Kv = ((Blin.transpose () * Snext * Blin) + Rcontrol_ ).inverse () *
                    Blin.transpose ();

                Ku = ((Blin.transpose () * Snext * Blin) + Rcontrol_ ).inverse () *
                    Rcontrol_;

                MatrixXd Snew = Alin.transpose () * Snext * (Alin - (Blin * K)) + Qstate_;

                valFun.col (i-1) = ((Alin - (Blin * K)).transpose () * valNext -
                        (K.transpose () * Rcontrol_ * controlSeq_->col (i)) +
                        (Qstate_ * stateTraj_->col (i)));

                std::cout << "K: " << K << std::endl;
                std::cout << "stateTraj: " << dStateTraj_->col (i) << std::endl;
                std::cout << "Kv: " << Kv << std::endl;
                std::cout << "valNext: " << valNext << std::endl;
                std::cout << "Ku: " << Ku << std::endl;
                std::cout << "controlseq: " << controlSeq_->col (i) << std::endl;
                dControlSeq_->col (i) =- ((K * dStateTraj_->col (i)) + (Kv * valNext) +
                        (Ku * controlSeq_->col (i)));
                std::cout<< "Here is the issue with Matrix Multiplication" << std::endl;

                valNext = valFun.col(i-1);
                Snext = Snew;
            }

            return dControlSeq_;
        }

        MatrixPtr_t SteeringILQR::steerState (VectorXd init, VectorXd final)
        {
           initState_ = init;
           finalState_ = final;
           std::vector <double> costVec;

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
           }

           return controlSeq_;
        } 
    }
}
