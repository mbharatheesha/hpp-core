// Copyright (C) 2014 LAAS-CNRS
// Author: Mathieu Geisert
//
// This file is part of the hpp-core.
//
// hpp-core is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// test-hpp is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with hpp-core.  If not, see <http://www.gnu.org/licenses/>.

#include <boost/assign.hpp>

#include <hpp/util/debug.hh>
#include <hpp/model/device.hh>
#include <hpp/model/joint.hh>
#include <hpp/core/fwd.hh>
#include <hpp/core/roadmap.hh>
#include <hpp/core/k-d-tree.hh>
#include <hpp/core/weighed-distance.hh>
#include "hpp/core/basic-configuration-shooter.hh"
#include <hpp/core/connected-component.hh>
#include <hpp/core/node.hh>
#include <hpp/model/joint-configuration.hh>

#include <hpp/core/steering-method-straight.hh>
#include <hpp/core/weighed-distance.hh>

#include <hpp/core/system-dynamics.hh>
#include <hpp/core/simple-pendulum.hh>
#include <hpp/core/steering-statespace.hh>
#include <hpp/core/steering-ilqr.hh>

#define BOOST_TEST_MODULE roadmap-1
#include <boost/test/included/unit_test.hpp>

using hpp::model::Configuration_t;
using hpp::core::ConfigurationPtr_t;
using hpp::model::JointPtr_t;
using hpp::model::Device;
using hpp::model::DevicePtr_t;
using hpp::model::JointTranslation;
using hpp::core::SteeringMethodStraight;
using hpp::core::RoadmapPtr_t;
using hpp::core::Roadmap;
using hpp::core::NodePtr_t;
using hpp::core::WeighedDistance;
using hpp::core::SystemDynamics;
using hpp::core::SystemDynamicsPtr_t;
using hpp::core::MatrixPtr_t;
using hpp::core::SimplePendulum;
using hpp::core::SteeringStateSpace;
using hpp::core::SteeringILQR;

BOOST_AUTO_TEST_SUITE( test_hpp_core )

BOOST_AUTO_TEST_CASE (Roadmap1) {
  // Build robot
  DevicePtr_t robot = Device::create("robot");
  JointPtr_t xJoint = new JointTranslation(fcl::Transform3f());
  xJoint->isBounded(0,1);
  xJoint->lowerBound(0,-3.);
  xJoint->upperBound(0,3.);
  JointPtr_t yJoint = new JointTranslation
    (fcl::Transform3f(fcl::Quaternion3f (sqrt (2)/2, 0, 0, sqrt(2)/2)));
  yJoint->isBounded(0,1);
  yJoint->lowerBound(0,-3.);
  yJoint->upperBound(0,3.);

  robot->rootJoint (xJoint);
  xJoint->addChildJoint (yJoint);

  // Create steering method
  SteeringMethodStraight sm (robot);
  // create roadmap
  RoadmapPtr_t r = Roadmap::create (WeighedDistance::create
				    (robot, boost::assign::list_of (1)(1)),
				    robot);
  ConfigurationPtr_t q_init (new Configuration_t (robot->configSize ()));
  (*q_init) [0] = 0; (*q_init) [1] = 0;
  // Set init node
  r->initNode (q_init);
  NodePtr_t n_init = r->initNode ();

  ConfigurationPtr_t q_goal (new Configuration_t (robot->configSize ()));
  (*q_goal) [0] = 1; (*q_goal) [1] = 0;
  NodePtr_t n_goal = r->addNode (q_goal);
  r->addGoalNode (q_goal);

  ConfigurationPtr_t q_randnode (new Configuration_t (robot->configSize ()));
  (*q_randnode)[0] = 0.5; (*q_randnode) [1] = 0.9;
  NodePtr_t n_randnode = r->addNode (q_randnode);

  ConfigurationPtr_t q_randnode1 (new Configuration_t (robot->configSize ()));
  (*q_randnode1)[0] = -0.1; (*q_randnode1) [1] = -0.9;
  NodePtr_t n_randnode1 = r->addNode (q_randnode1);

  ConfigurationPtr_t q_randnode2 (new Configuration_t (robot->configSize ()));
  (*q_randnode2)[0] = 1.5; (*q_randnode2) [1] = 2.9;
  NodePtr_t n_randnode2 = r->addNode (q_randnode2);



  r->addEdge (n_init, n_goal, sm (*q_init, *q_goal));
  r->addEdge (n_goal, n_init, sm (*q_goal, *q_init));
  r->addEdge (n_goal, n_randnode1, sm (*q_goal, *q_randnode1));
  r->addEdge (n_init, n_randnode, sm (*q_init, *q_randnode));
  r->addEdge (n_randnode, n_randnode1, sm (*q_randnode, *q_randnode1));
  r->addEdge (n_randnode1, n_randnode2, sm (*q_randnode1, *q_randnode2));
  r->addEdge (n_randnode2, n_randnode, sm (*q_randnode2, *q_randnode));
  r->addEdge (n_randnode, n_goal, sm (*q_randnode, *q_goal));

  std::cout << *r << std::endl;

  SystemDynamics* ptr = new SimplePendulum ();
  
//  SimplePendulum pendulum;
//
//  ptr = &pendulum;

  SystemDynamicsPtr_t shPtr (ptr);

  
  int numDOF = 1;
  VectorXd initState;
  VectorXd dState;
  VectorXd tVec = VectorXd::LinSpaced (1000,0,5);
  MatrixXd stateTraj;

  shPtr->setProblemDimension (numDOF);

  initState = 0*(VectorXd::Ones (2*numDOF));

  initState(1) = 0.9;

  static_cast<SimplePendulum*>(ptr)->setParameters ();

  stateTraj = MatrixXd::Zero (initState.size(),tVec.size());

  stateTraj = shPtr->simulateDynamics (tVec, initState);

  std::cout << "Dynamics of SimplePendulum were integrated successfully" << std::endl;

  std::ofstream file_out;
  file_out.open("test.dat");  

  //std::cout << stateTraj;
  file_out << stateTraj;
  file_out.close();


  /// Testing for the ILQR functionality
  SteeringILQR iLQRsteer;

  int ilqrIter = 40;
  int numIn = 1;
  int lenIn = 450;
  double dt = 0.05;
  MatrixPtr_t statePath;

  VectorXd stateWeight (numDOF*2);
  VectorXd finalStateWeight (numDOF*2);
  VectorXd controlWeight (numIn);

  stateWeight.fill (0.0001);
  finalStateWeight.fill (1);
  controlWeight.fill (1);

  initState (0) = 0;
  initState (1) = 0;

  VectorXd finalState (numDOF*2);
  finalState (0) = 3.1416;
  finalState (1) = 0;

  iLQRsteer.setILQRParams (ilqrIter, numIn, lenIn, dt);
  iLQRsteer.setWeightingMatrices (stateWeight, finalStateWeight, controlWeight);
  iLQRsteer.setSysDynPtr (shPtr);

  statePath = iLQRsteer.steerState (initState, finalState);
  
  std::ofstream file_path;
  file_path.open("state_path.dat");  

  file_path << *statePath;
  file_path.close();

  BOOST_CHECK (r->pathExists ());
}
BOOST_AUTO_TEST_SUITE_END()



