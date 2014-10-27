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
# include <hpp/core/system-dynamics.hh>
# include <hpp/core/steering-statespace.hh>

# define HPP_CORE_FDSTEP_SIZE 1.0e-17
# define HPP_CORE_FDSTEP_SIZE_INV 1/HPP_CORE_FDSTEP_SIZE

namespace hpp {
    namespace core {
        class HPP_CORE_DLLAPI SteeringILQR : public SteeringStateSpace {
            protected:
                int numIter_;
                int numControl_;
                int controlLen_;
                double DT_;
                MatrixXd Qstate_;
                MatrixXd Qfinal_;
                MatrixXd Rcontrol_;
                MatrixXd valFun_;
                MatrixPtr_t controlSeq_;
                MatrixXd* dControlSeq_;
                MatrixPtr_t stateTraj_;
                MatrixXd* dStateTraj_;
                std::vector<MatrixXd> linMatrices_;
                SystemDynamicsPtr_t sysDyn_;

            public:
                /// Constructor
                SteeringILQR () : SteeringStateSpace (), numIter_ (20), numControl_ (1),
                controlLen_ (100), DT_ (0.005), Qstate_ (), Qfinal_ (), Rcontrol_ (), valFun_ (),
                linMatrices_ (), sysDyn_ ()
            {
            }
                /// Get the dynamical System Information
                void setSysDynPtr (const SystemDynamicsPtr_t& sysDynPtr)
                {
                    sysDyn_ = sysDynPtr;
                }
                
                /// Create Shared pointers for the state and control sequences
                /// Allocate memory sizes for all data
                void initializeMemory ( void )
                {
                    MatrixXd* ptrCtl = new MatrixXd (numControl_, controlLen_);
                    MatrixPtr_t shPtrCtl (ptrCtl);
                    controlSeq_ = shPtrCtl;

                    MatrixXd* ptrState = new MatrixXd (initState_.size (), controlLen_);
                    MatrixPtr_t shPtrState (ptrState);
                    stateTraj_ =  shPtrState;

                    dStateTraj_ = new MatrixXd ();
                    dStateTraj_->resize (initState_.size (), controlLen_);
                    dStateTraj_->setZero ();

                    dControlSeq_ = new MatrixXd ();
                    dControlSeq_->resize (numControl_, controlLen_);
                    dControlSeq_->setConstant (1.0e-6);

                    valFun_ = MatrixXd::Zero (initState_.size (), controlLen_);
                } 


                /// Set the time discretization parameter
                void setILQRParams (int, int, int, double);
                
                /// Set the weighting matrices for computing cost
                void setWeightingMatrices (vectorIn_t, vectorIn_t, vectorIn_t);
                
                /// Cost Function
                double computeCost ();
                
                /// Forward Pass
                void forwardPass (const MatrixPtr_t&); 

                /// Linearize nonlinear dynamics around a given trajectory and input
                void linearizeSystem ();

                /// Forward Finite Difference Linearization
                void linearizeFDJacobian (vectorIn_t, vectorIn_t, matrixOut_t);

                /// Find optimal control improvement by iteratively solving LQR
                MatrixXd* deltaOptControl ();

                /// Steer from initial state to final state
                MatrixPtr_t steerState (vectorIn_t, vectorIn_t);
        };
    }

}
#endif
