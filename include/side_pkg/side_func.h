/* Author: 	Emanuele Sansebastiano
 * Date:	August 2017
 *
 * Desc:	Library to encapsulate some useful side functions of everyday coding
 *
 */

#ifndef SIDE_PKG_SIDE_FUNC_H
#define SIDE_PKG_SIDE_FUNC_H

// ROS
#include <ros/ros.h>

// C++
#include <geometry_msgs/Vector3.h>
#include <geometry_msgs/PoseStamped.h>
#include <boost/math/constants/constants.hpp>
////#include <geometry_msgs/PointStamped.h>
#include <tf_conversions/tf_kdl.h>

//////////////////////////////////////////////////////////////////////////////////////
// VALUES MODIFIABLE BY THE USER \\

// Common define values
#define	std_time					0.01

//////////////////////////////////////////////////////////////////////////////////////

/*HOW TO USE THIS LIBRARY
   At the beginning of you program you must insert the following lines:
 * 	ros::init(argc, argv, "your_node_name");
 *	ros::NodeHandle node_handle("~");
 *	ros::AsyncSpinner spinner(1);
 *  spinner.start();
 */

namespace basic_side_classes
{
  //brief: Class to sum up all the most used char and string
  class UsefulCharString
  {
  	  public:
	  	  UsefulCharString()
  	  	  {};
	  	  ~UsefulCharString()
	  	  {};
	  	  //brief: set of function to get common chars to manage (write and read) .txt files
	  	  char get_newline_char(void);
	  	  char get_tab_char(void);
	  	  char get_comment_char(void);

	  	  //brief: set of function to get common strings to communicate with the user
	  	  std::string get_invalid_input_str(void);

  	  protected:
	  	  //brief: set of common chars used to generate and read .txt files
	  	  char newline_ = '\n';
	  	  char tab_ = '\t';
	  	  char comment_ = '%';

	  	  //brief: set of common strings to communicate with the user
	      std::string invalid_input_str_ = "Error: An invalid input has been inserted!";

  };

// End namespace "basic_side_classes"
}


namespace basic_side_functions
{
  //brief: Function to generate a waiting countdown
  void countdown_sec(int n_sec = 3);

  //brief: Standard sleep time - default value [seconds]
  void standardSleep(double sleep_time = std_time);

  //brief: Function to get the current PC's time in milliseconds
  long int getCurrentTime();

  //brief: Convert radiant to degree
  double rad2deg(double rad);

  //brief: Convert degree to radiant
  double deg2rad(double deg);

  //brief: Absolute function to convert double values
  double abs_f(double val);

  //brief: Round a value to the closest integer
  //       examples: 2.4 --> 2; 2.7 --> 3; 2.0 --> 2; 2.5 --> 3;
  int round_f(double val);

  //brief: Function to shift decimal values and cut off the non-integer part
  int decimal_shift(double val2shift, int decimal_considered);

  //brief: Function to cluster to vector in one vector
  std::vector<double> vector_cluster_double(std::vector<double> left, std::vector<double> right);

  //brief: Function to print on video the matrix mA
  void matrix2DPRINT(std::vector <std::vector<double> > mA);

  //brief: Set of function to Sum or Subtract two matrixes having 2 or 3 dimensions
  bool matrix2DSUM(std::vector <std::vector<double> > mA, std::vector <std::vector<double> > mB, std::vector <std::vector<double> > &mRes);
  bool matrix2DSUBTRACTION(std::vector <std::vector<double> > mA, std::vector <std::vector<double> > mB, std::vector <std::vector<double> > &mRes);
  bool matrix3DSUM(std::vector<std::vector <std::vector<double> > > mA, std::vector<std::vector <std::vector<double> > > mB, std::vector<std::vector <std::vector<double> > > &mRes);
  bool matrix3DSUBTRACTION(std::vector<std::vector <std::vector<double> > > mA, std::vector<std::vector <std::vector<double> > > mB, std::vector<std::vector <std::vector<double> > > &mRes);

  //brief: Function to make the product of two 2D matrixes defined as vector of vectors of double
  bool matrix2DPROD(std::vector <std::vector <double> > mA, std::vector <std::vector<double> > mB, std::vector <std::vector<double> > &mProd);

  //brief: Function to transpose mA (not just squared ones)
  bool matrix2DTRAN(std::vector <std::vector <double> > &mA);

  //brief: Function to transform mA to its upper triangular matrix (square matrix)
  bool matrix2DTRIANup(std::vector <std::vector <double> > &mA);

  //brief: Function to calculate the determinant of a square 2D matrix
  bool matrix2DDET(std::vector <std::vector <double> > mA, double &det);

  //brief: Function to transform mA to its inverse matrix (square matrix)
  bool matrix2DINVERT(std::vector <std::vector <double> > &mA);

  //brief: This function compares two vectors of double by a threshold
  bool VectorEquivalence_theshold(std::vector<double> A, std::vector<double> B, double threshold = 0.0);

  //brief: This function checks if a topic exist or not
  bool CheckTopicExistence(std::string topic2check);

// End namespace "basic_side_functions"
}

namespace geometry_side_functions
{
  //brief: Make a Vector3 data in just one line
  geometry_msgs::Vector3 makeVector3(double x, double y, double z);

  //brief: Make a Quaternion data in just one line
 geometry_msgs::Quaternion makeQuat(double w, double x, double y, double z);

  //brief: Function to convert RPY to quaternion angles
  //       If RPY is NOT in radiant unit change RPY_rad to false
  //       If you do not want to arrange the angle in to interval (-pi ; pi] change turn_corr to false
  geometry_msgs::Quaternion RPY2Quat(geometry_msgs::Vector3 RPY, bool RPY_rad = true, bool turn_corr = true);

  //brief: Function to convert quaternion into RPY angles
  //       If you want RPY NOT in radiant unit change RPY_rad to false
  geometry_msgs::Vector3 Quat2RPY(geometry_msgs::Quaternion Quat, bool RPY_rad = true);

  //brief: Function to make a Pose message from quaternions
  geometry_msgs::Pose makePose(geometry_msgs::Quaternion orientation, geometry_msgs::Vector3 XYZ_location);
  //brief: Function to make a Pose message from eulerian angles in degree
  geometry_msgs::Pose makePose(geometry_msgs::Vector3 RPY_orientation, geometry_msgs::Vector3 XYZ_location);

  //brief: Convert a Pose data to a PoseStamped data
  geometry_msgs::PoseStamped Pose2PoseStamped( geometry_msgs::Pose old_pose);
  //brief: Convert a Pose vector to a PoseStamped vector
  std::vector<geometry_msgs::PoseStamped> Pose2PoseStamped( std::vector<geometry_msgs::Pose> old_pose);

  //brief: Convert a PoseStamped data to a Pose data
  geometry_msgs::Pose PoseStamped2Pose( geometry_msgs::PoseStamped old_pose_s);
  //brief: Convert a PoseStamped vector to a Pose vector
  std::vector<geometry_msgs::Pose> PoseStamped2Pose( std::vector<geometry_msgs::PoseStamped> old_pose_s);

  //brief: Function to convert easily a Vector3 data into a Pose.position data
  void vector32posePosition(geometry_msgs::Vector3 vector, geometry_msgs::Pose &pose2update);
  void vector32posePosition(geometry_msgs::Vector3 vector, geometry_msgs::PoseStamped &poseStamped2update);

  //brief: Function to convert easily a Pose.position data into a Vector3 data
  void posePosition2vector3(geometry_msgs::Pose pose, geometry_msgs::Vector3 &vector2update);
  void posePosition2vector3(geometry_msgs::PoseStamped poseStamped, geometry_msgs::Vector3 &vector2update);

  //brief: This function returns the sum of the absolute difference of every term of two quaternions
  double Quat_diff(geometry_msgs::Quaternion a, geometry_msgs::Quaternion b);

  //brief: This function compares two poses until a specific decimal value
  //       'decimal_considered' is the number of decimals you want to consider - negative values mean 'all the available values'
  bool PoseEquivalence_decimal_value(geometry_msgs::Pose A, geometry_msgs::Pose B, int decimal_considered = -1);

  //brief: This function compares two poses by a threshold (true output if the values are equivalent)
  //       Pose.position check
  bool PoseEquivalence_XYZ(geometry_msgs::Pose A, geometry_msgs::Pose B, double threshold_XYZ = 0.0);
  //       Pose.orientation check
  bool PoseEquivalence_Quat(geometry_msgs::Pose A, geometry_msgs::Pose B, double threshold_Quat = 0.0);
  //       Pose full check
  bool PoseEquivalence_theshold(geometry_msgs::Pose A, geometry_msgs::Pose B, double threshold_XYZ = 0.0, double threshold_Quat = 0.0);

  /*//brief:: Function to easy generate point matrix representing a point in the 3D
  std::vector <std::vector <double> > point2matrix_gen(geometry_msgs::Vector3 point);

  //brief: Function to easy convert a matrix representing a point in 3D into a vector3 data
  geometry_msgs::Vector3 matrix2point_gen(std::vector <std::vector <double> > matrix);*/

  //brief: Function to generate the translation matrix from the translation vector
  std::vector <std::vector <double> > translation_matrix(geometry_msgs::Vector3 translation);

  //brief:: Function to extract the translation values from the full transfer matrix
  geometry_msgs::Vector3 transMatrix2XYZ(std::vector <std::vector <double> > matrix);

  //brief: Set of function to generate the rotational matrix around every single axis
  //       this rotation corresponds to the RPY
  //       if the angles are NOT in radiant put 'rad = false' to convert them
  std::vector <std::vector <double> > rotational_matrix_X(double angle_X, bool rad = true);
  std::vector <std::vector <double> > rotational_matrix_Y(double angle_Y, bool rad = true);
  std::vector <std::vector <double> > rotational_matrix_Z(double angle_Z, bool rad = true);

  //brief:: Function to generate the full rotational matrix from the RPY axis rotation
  //	    if the angle unit is NOT radiant put 'rad = false' to convert
  std::vector <std::vector <double> > RPY2rotMatrix(geometry_msgs::Vector3 angles, bool rad = true);
  //brief:: Function to generate the full rotational matrix from the RPY axis rotation starting from quaternions
  std::vector <std::vector <double> > Quaternion2rotMatrix(geometry_msgs::Quaternion quat);

  //brief:: Function to extract the RPY axis rotation from the full transfer matrix
  //	    if the angle unit you want is NOT radiant put 'rad = false' to convert
  geometry_msgs::Vector3 rotMatrix2RPY(std::vector<std::vector <double> > matrix, bool rad = true);

  //brief:: Function to extract the quaternion of the axis rotation from the full transfer matrix
  geometry_msgs::Quaternion rotMatrix2Quaternion( std::vector <std::vector <double> > matrix);

// End namespace "geometry_side_functions"
}

#endif /* SIDE_PKG_SIDE_FUNC_H */

