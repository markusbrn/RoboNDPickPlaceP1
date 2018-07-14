#!/usr/bin/env python

# Copyright (C) 2017 Udacity Inc.
#
# This file is part of Robotic Arm: Pick and Place project for Udacity
# Robotics nano-degree program
#
# All Rights Reserved.

# Author: Harsh Pandya

# import modules
import rospy
import tf
from kuka_arm.srv import *
from trajectory_msgs.msg import JointTrajectory, JointTrajectoryPoint
from geometry_msgs.msg import Pose
from mpmath import *
from sympy import *

########################################################################################
class IK_Service:
    def __init__(self):
        # Create symbols
        self.q1, self.q2, self.q3, self.q4, self.q5, self.q6, self.qG = symbols('q1:8')
        self.d1, self.d2, self.d3, self.d4, self.d5, self.d6, self.dG = symbols('d1:8')
        self.a0, self.a1, self.a2, self.a3, self.a4, self.a5, self.a6 = symbols('a0:7')
        self.alpha0, self.alpha1, self.alpha2, self.alpha3, self.alpha4, self.alpha5, self.alpha6 = symbols('alpha0:7')
        # Create Modified DH parameters
        self.s = {self.alpha0:     0,   self.a0:      0,   self.d1:  0.75,
                  self.alpha1: -pi/2,   self.a1:   0.35,   self.d2:     0,   self.q2: self.q2-pi/2,
                  self.alpha2:     0,   self.a2:   1.25,   self.d3:     0,
                  self.alpha3: -pi/2,   self.a3: -0.054,   self.d4:  1.50,
                  self.alpha4:  pi/2,   self.a4:      0,   self.d5:     0,
                  self.alpha5: -pi/2,   self.a5:      0,   self.d6:     0,
                  self.alpha6:     0,   self.a6:      0,   self.dG: 0.303,   self.qG: 0}
        # Define Modified DH Transformation matrix
        def Matrix_DH(alpha, a, d, q):
            return Matrix([[           cos(q),           -sin(q),           0,             a],
                           [sin(q)*cos(alpha), cos(q)*cos(alpha), -sin(alpha), -sin(alpha)*d],
                           [sin(q)*sin(alpha), cos(q)*sin(alpha),  cos(alpha),  cos(alpha)*d],
                           [                0,                 0,           0,             1]])
        # Create individual transformation matrices
        T_0_1 = Matrix_DH(self.alpha0, self.a0, self.d1, self.q1).subs(self.s)
        T_1_2 = Matrix_DH(self.alpha1, self.a1, self.d2, self.q2).subs(self.s)
        T_2_3 = Matrix_DH(self.alpha2, self.a2, self.d3, self.q3).subs(self.s)
        T_3_4 = Matrix_DH(self.alpha3, self.a3, self.d4, self.q4).subs(self.s)
        T_4_5 = Matrix_DH(self.alpha4, self.a4, self.d5, self.q5).subs(self.s)
        T_5_6 = Matrix_DH(self.alpha5, self.a5, self.d6, self.q6).subs(self.s)
        T_6_G = Matrix_DH(self.alpha6, self.a6, self.dG, self.qG).subs(self.s)
        self.T_0_3 = (T_0_1 * T_1_2 * T_2_3)
        # Compensate for rotation discrepancy between DH parameters and Gazebo
        T_z = Matrix([[cos(pi),   -sin(pi),   0,   0],
                      [sin(pi),    cos(pi),   0,   0],
                      [         0,       0,   1,   0],
                      [         0,       0,   0,   1]])

	T_y = Matrix([[ cos(-pi/2),   0,   sin(-pi/2),   0],
                      [          0,   1,            0,   0],                      
                      [-sin(-pi/2),   0,   cos(-pi/2),   0],
                      [          0,   0,            0,   1]])

        self.T_U_G = T_z * T_y #transformation from gripper to urdf COS
        self.T_3_U = T_3_4 * T_4_5 * T_5_6 * T_6_G * self.T_U_G**-1

    def calculate_IK(self, req):
        rospy.loginfo("Received %s eef-poses from the plan" % len(req.poses))
        if len(req.poses) < 1:
            print "No valid poses received"
            return -1
        else:
            # Initialize service response
            joint_trajectory_list = []
            for x in xrange(0, len(req.poses)):
                req_pose = req.poses[x]
                # IK code starts here
                joint_trajectory_point = JointTrajectoryPoint()        

                # Initialize service response
                # IK code starts here
                # Extract end-effector position and orientation from request
                # px,py,pz = end-effector position
                # roll, pitch, yaw = end-effector orientation
                px = req_pose.position.x
                py = req_pose.position.y
                pz = req_pose.position.z

                (roll, pitch, yaw) = tf.transformations.euler_from_quaternion(
                    [req_pose.orientation.x, req_pose.orientation.y,
                        req_pose.orientation.z, req_pose.orientation.w])

                ### Your IK code here
                # Calculate joint angles using Geometric IK method
                #1.) calculate rotation matrix from 0 to G with the recieved robot pose message orientations
                R_roll = Matrix([[ 1,              0,          0],
                                 [ 0,      cos(roll), -sin(roll)],
                                 [ 0,      sin(roll),  cos(roll)]])

                R_pitch = Matrix([[ cos(pitch),        0,  sin(pitch)],
                                  [          0,        1,           0],
                                  [-sin(pitch),        0,  cos(pitch)]])

                R_yaw = Matrix([[ cos(yaw), -sin(yaw),        0],
                                [ sin(yaw),  cos(yaw),        0],
                                [        0,         0,        1]])

                R_0_U = R_yaw * R_pitch * R_roll
                R_0_G = R_0_U * self.T_U_G[0:3,0:3]

                #2.) calculate the wrist center position
                self.WC = Matrix([[px-self.dG*R_0_G[0,2]],
                                  [py-self.dG*R_0_G[1,2]],
                                  [pz-self.dG*R_0_G[2,2]]]).evalf(subs = self.s)

                #3.) compute q1-3
                theta1 = atan2(self.WC[1],self.WC[0])

                A = self.d4.evalf(subs = self.s)
                B = (sqrt((sqrt(self.WC[0]**2+self.WC[1]**2)-self.a1)**2 + (self.WC[2]-self.d1)**2)).evalf(subs = self.s)
                C = self.a2.evalf(subs = self.s)

                alpha = acos((B**2+C**2-A**2)/(2*B*C))
                beta = acos((A**2+C**2-B**2)/(2*A*C))
                gamma = acos((A**2+B**2-C**2)/(2*A*B))

                theta2 = pi/2 - alpha - (atan2(self.WC[2]-self.d1, sqrt(self.WC[0]**2+self.WC[1]**2)-self.a1)).evalf(subs = self.s)
                theta3 = pi/2 - (beta + 0.036)

                #4.) Calculate R_3_0
                R_3_0 = self.T_0_3[:3,:3].evalf(subs = {self.q1:theta1, self.q2:theta2, self.q3:theta3}).transpose()

                R_3_G = R_3_0 * R_0_G
            
                #5.) Compute q4-6
                theta5 = atan2(sqrt(R_3_G[2,2]**2+R_3_G[0,2]**2),R_3_G[1,2])
                
                #handle multiple solutions (I found this solution in a slack posting :-)
                if sin(theta5) < 0:
                    theta4 = atan2(-R_3_G[2,2],R_3_G[0,2])
                    theta6 = atan2(R_3_G[1,1],-R_3_G[1,0])
                else:
                    theta4 = atan2(R_3_G[2,2],-R_3_G[0,2])
                    theta6 = atan2(-R_3_G[1,1],R_3_G[1,0])

                # Populate response for the IK request
                # In the next line replace theta1,theta2...,theta6 by your joint angle variables
	        joint_trajectory_point.positions = [theta1, theta2, theta3, theta4, theta5, theta6]
	        joint_trajectory_list.append(joint_trajectory_point)

            rospy.loginfo("length of Joint Trajectory List: %s" % len(joint_trajectory_list))
            return CalculateIKResponse(joint_trajectory_list)
########################################################################################

def IK_server():
    # initialize node and declare calculate_ik service
    rospy.init_node('IK_server')
    IK = IK_Service()
    s = rospy.Service('calculate_ik', CalculateIK, IK.calculate_IK)
    print "Ready to receive an IK request"
    rospy.spin()

if __name__ == "__main__":
    IK_server()
