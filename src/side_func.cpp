/* Author: 	Emanuele Sansebastiano
 * Date:	August 2017
 *
 * Desc:	Library to encapsulate some useful side functions of everyday coding
 *
 */

// this pkg
#include <side_pkg/side_func.h>

namespace bsc = basic_side_classes;
char bsc::UsefulCharString::get_newline_char()
{
	  return newline_;
}
char bsc::UsefulCharString::get_tab_char()
{
	  return tab_;
}
char bsc::UsefulCharString::get_comment_char()
{
	  return comment_;
}

std::string bsc::UsefulCharString::get_invalid_input_str()
{
	  return invalid_input_str_;
}
// End namespace "basic_side_classes"


namespace basic_side_functions
{
  void countdown_sec(int n_sec)
  {
	std::cout << "countdown: " << n_sec << std::endl;
	for(int i = 0; i < n_sec; i++)
    {
		ros::Duration(1.0).sleep();
		std::cout << "countdown: " << n_sec-i-1 << std::endl;
	}
  }

  void standardSleep(double sleep_time)
  {
	  ros::Duration(sleep_time).sleep();
  }

  long int getCurrentTime(void)
  {
	  struct timeval tp;
	  gettimeofday(&tp, NULL);
	  return tp.tv_sec * 1000 + tp.tv_usec / 1000;
  }

  double rad2deg(double rad)
  {
	  double pi = boost::math::constants::pi<double>();
	  int deg_pi = 180;
	  double deg = deg_pi*rad/pi;
	  return deg;
  }

  double deg2rad(double deg)
  {
	  double pi = boost::math::constants::pi<double>();
	  int deg_pi = 180;
	  double rad = pi*deg/deg_pi;
	  return rad;
  }

  double abs_f(double val)
  {
	  double abs_v;
	  if (val >= 0)
		  abs_v = val;
	  else
		  abs_v = -val;
	  return abs_v;
  }

  int round_f(double val)
  {
	  int temp = (int) val;
	  double temp1 = val - temp;
	  // if val is positive
	  if (abs_f(temp1) >= 0.5 && temp1 >= 0)
		  temp = val +1;
	  // if val is negative
	  else if (abs_f(temp1) >= 0.5 && temp1 < 0)
		  temp =  val -1;
	  //else {temp is already the integer part of val}
	  return temp;
  }

  int decimal_shift(double val2shift, int decimal_considered)
  {
	  int val2return = val2shift*std::pow(10,decimal_considered);
	  return val2return;
  }

  std::vector<double> vector_cluster_double(std::vector<double> left, std::vector<double> right)
  {
	  int final_size = right.size() + left.size();
	  std::vector<double> final_vector; final_vector.resize(final_size);

	  for (int i = 0; i < left.size(); i++)
	  		  final_vector[i] = left[i];
	  for (int i = 0; i < right.size(); i++)
	  		  final_vector[i + left.size()] = right[i];

	  return final_vector;
  }

  bool matrix2DSUM(std::vector <std::vector<double> > mA, std::vector <std::vector<double> > mB, std::vector <std::vector<double> > &mRes)
  {
	  int temp_int;

	  //function check
	  //matrix size check
	  temp_int = mA[0].size();
	  for(int i = 0; i < mA.size(); i++)
	  {
		  if(mA[i].size() != temp_int){
			  std::cout << "Error: the raws of the first matrix have NOT constant size!" << std::endl;
			  goto return_0;
		  }
	  }
	  temp_int = mB[0].size();
	  for(int i = 0; i < mB.size(); i++)
	  {
		  if(mB[i].size() != temp_int){
			  std::cout << "Error: the raws of the second matrix have NOT constant size!" << std::endl;
			  goto return_0;
		  }
	  }
	  //proper dimensions check
	  if(mA.size() != mB.size() || mA[0].size() != mB[0].size())
	  {
		  std::cout << "Error: the function is trying to make the sum of matrixes having the following dimensions:" << std::endl;
		  std::cout << mA.size() << " x " << mA[0].size() << "  and  " << mB.size() << " x " << mB[0].size() << std::endl;
		  goto return_0;
	  }

	  //main program
	  mRes = mA;
	  for(int i = 0; i < mA.size(); i++)
	  {
		  for(int j = 0; j < mA[0].size(); j++)
		  {
			  mRes[i][j] = mA[i][j] + mB[i][j];
		  }
	  }

	  return true;

	  //exit return in case something is not correct
	  return_0:
	  return false;
  }

  bool matrix2DSUBTRACTION(std::vector <std::vector<double> > mA, std::vector <std::vector<double> > mB, std::vector <std::vector<double> > &mRes)
  {
	  int temp_int;

	  //function check
	  //matrix size check
	  temp_int = mA[0].size();
	  for(int i = 0; i < mA.size(); i++)
	  {
		  if(mA[i].size() != temp_int){
			  std::cout << "Error: the raws of the first matrix have NOT constant size!" << std::endl;
			  goto return_0;
		  }
	  }
	  temp_int = mB[0].size();
	  for(int i = 0; i < mB.size(); i++)
	  {
		  if(mB[i].size() != temp_int){
			  std::cout << "Error: the raws of the second matrix have NOT constant size!" << std::endl;
			  goto return_0;
		  }
	  }
	  //proper dimensions check
	  if(mA.size() != mB.size() || mA[0].size() != mB[0].size())
	  {
		  std::cout << "Error: the function is trying to make the subtraction of matrixes having the following dimensions:" << std::endl;
		  std::cout << mA.size() << " x " << mA[0].size() << "  and  " << mB.size() << " x " << mB[0].size() << std::endl;
		  goto return_0;
	  }

	  //main program
	  mRes = mA;
	  for(int i = 0; i < mA.size(); i++)
	  {
		  for(int j = 0; j < mA[0].size(); j++)
		  {
			  mRes[i][j] = mA[i][j] - mB[i][j];
		  }
	  }

	  return true;

	  //exit return in case something is not correct
	  return_0:
	  return false;
  }

  bool matrix3DSUM(std::vector<std::vector <std::vector<double> > > mA, std::vector<std::vector <std::vector<double> > > mB, std::vector<std::vector <std::vector<double> > > &mRes)
  {
	  int temp_int;

	  //function check
	  //matrix size check
	  temp_int = mA[0].size();
	  for(int i = 0; i < mA.size(); i++)
	  {
		  if(mA[i].size() != temp_int){
			  std::cout << "Error: the raws of the first matrix have NOT constant size!" << std::endl;
			  goto return_0;
		  }
	  }
	  temp_int = mA[0][0].size();
	  for(int i = 0; i < mA.size(); i++)
	  {
		  for(int j = 0; j < mA[0].size(); j++)
		  {
			  if(mA[i][j].size() != temp_int){
				  std::cout << "Error: the depth lines of the first matrix have NOT constant size!" << std::endl;
				  goto return_0;
			  }
		  }
	  }
	  temp_int = mB[0].size();
	  for(int i = 0; i < mB.size(); i++)
	  {
		  if(mB[i].size() != temp_int){
			  std::cout << "Error: the raws of the second matrix have NOT constant size!" << std::endl;
			  goto return_0;
		  }
	  }
	  temp_int = mB[0][0].size();
	  for(int i = 0; i < mB.size(); i++)
	  {
		  for(int j = 0; j < mB[0].size(); j++)
		  {
			  if(mB[i][j].size() != temp_int){
				  std::cout << "Error: the depth lines of the second matrix have NOT constant size!" << std::endl;
				  goto return_0;
			  }
		  }
	  }
	  //proper dimensions check
	  if(mA.size() != mB.size() || mA[0].size() != mB[0].size() || mA[0][0].size() != mB[0][0].size())
	  {
		  std::cout << "Error: the function is trying to make the sum of matrixes having the following dimensions:" << std::endl;
		  std::cout << mA.size() << " x " << mA[0].size() << " x " << mA[0][0].size() << "  and  " << mB.size() << " x " << mB[0].size() << " x " << mB[0][0].size() << std::endl;
		  goto return_0;
	  }

	  //main program
	  mRes = mA;
	  for(int i = 0; i < mA.size(); i++)
	  {
		  for(int j = 0; j < mA[0].size(); j++)
		  {
			  for(int z = 0; z < mA[0][0].size(); z++)
			  {
				  mRes[i][j][z] = mA[i][j][z] + mB[i][j][z];
			  }
		  }
	  }

	  return true;

	  //exit return in case something is not correct
	  return_0:
	  return false;
  }

  bool matrix3DSUBTRACTION(std::vector<std::vector <std::vector<double> > > mA, std::vector<std::vector <std::vector<double> > > mB, std::vector<std::vector <std::vector<double> > > &mRes)
  {
	  int temp_int;

	  //function check
	  //matrix size check
	  temp_int = mA[0].size();
	  for(int i = 0; i < mA.size(); i++)
	  {
		  if(mA[i].size() != temp_int){
			  std::cout << "Error: the raws of the first matrix have NOT constant size!" << std::endl;
			  goto return_0;
		  }
	  }
	  temp_int = mA[0][0].size();
	  for(int i = 0; i < mA.size(); i++)
	  {
		  for(int j = 0; j < mA[0].size(); j++)
		  {
			  if(mA[i][j].size() != temp_int){
				  std::cout << "Error: the depth lines of the first matrix have NOT constant size!" << std::endl;
				  goto return_0;
			  }
		  }
	  }
	  temp_int = mB[0].size();
	  for(int i = 0; i < mB.size(); i++)
	  {
		  if(mB[i].size() != temp_int){
			  std::cout << "Error: the raws of the second matrix have NOT constant size!" << std::endl;
			  goto return_0;
		  }
	  }
	  temp_int = mB[0][0].size();
	  for(int i = 0; i < mB.size(); i++)
	  {
		  for(int j = 0; j < mB[0].size(); j++)
		  {
			  if(mB[i][j].size() != temp_int){
				  std::cout << "Error: the depth lines of the second matrix have NOT constant size!" << std::endl;
				  goto return_0;
			  }
		  }
	  }
	  //proper dimensions check
	  if(mA.size() != mB.size() || mA[0].size() != mB[0].size() || mA[0][0].size() != mB[0][0].size())
	  {
		  std::cout << "Error: the function is trying to make the subtraction of matrixes having the following dimensions:" << std::endl;
		  std::cout << mA.size() << " x " << mA[0].size() << " x " << mA[0][0].size() << "  and  " << mB.size() << " x " << mB[0].size() << " x " << mB[0][0].size() << std::endl;
		  goto return_0;
	  }

	  //main program
	  mRes = mA;
	  for(int i = 0; i < mA.size(); i++)
	  {
		  for(int j = 0; j < mA[0].size(); j++)
		  {
			  for(int z = 0; z < mA[0][0].size(); z++)
			  {
				  mRes[i][j][z] = mA[i][j][z] - mB[i][j][z];
			  }
		  }
	  }

	  return true;

	  //exit return in case something is not correct
	  return_0:
	  return false;
  }

  bool matrix2Dprod(std::vector <std::vector<double> > mA, std::vector <std::vector<double> > mB, std::vector <std::vector<double> > &mProd)
  {
	  int temp_int;
	  double temp_double;
	  mProd.clear();

	  //function check
	  //matrix size check
	  temp_int = mA[0].size();
	  for(int i = 0; i < mA.size(); i++)
	  {
		  if(mA[i].size() != temp_int){
			  std::cout << "Error: the raws of the first matrix have NOT constant size!" << std::endl;
			  goto return_0;
		  }
	  }
	  temp_int = mB[0].size();
	  for(int i = 0; i < mB.size(); i++)
	  {
		  if(mB[i].size() != temp_int){
			  std::cout << "Error: the raws of the second matrix have NOT constant size!" << std::endl;
			  goto return_0;
		  }
	  }
	  //proper dimensions check
	  if(mA[0].size() != mB.size())
	  {
		  std::cout << "Error: the function is trying to make the product of matrixes having the following dimensions:" << std::endl;
		  std::cout << mA.size() << " x " << mA[0].size() << "  and  " << mB.size() << " x " << mB[0].size() << std::endl;
		  goto return_0;
	  }

	  //output definition
	  mProd.resize(mA.size());
	  for(int i = 0; i < mA.size(); i++)
	  {
		  mProd[i].resize(mB[0].size());
	  }

	  //main program
	  for(int j = 0; j < mB[0].size(); j++)
	  {
		  for(int i = 0; i < mA.size(); i++)
		  {
			  temp_double = 0.0;
			  for(int z = 0; z < mB.size(); z++)
			  {
				  temp_double += mA[i][z]*mB[z][j];
			  }
			  mProd[i][j] = temp_double;
		  }
  	  }

	  return true;

	  //exit return in case something is not correct
	  return_0:
	  return false;
  }

  bool VectorEquivalence_theshold(std::vector<double> A, std::vector<double> B, double threshold)
  {
	  bool equal = false;
	  //size check
	  if(A.size() != B.size()){
		  perror("The two vector has different size, a false value has been returned\n");
		  return false;
	  }
  	  if(threshold < 0.0)
  	  {
  		 std::cout << "Warning: The XYZ threshold must be positive (it has been converted to the positive value)" << std::endl;
  		 threshold = abs_f(threshold);
  	  }

  	  equal = true;
  	  for(int i = 0; i < A.size(); i++)
  	  {
  		  if(A[i] > B[i] + threshold || A[i] < B[i] - threshold)
  	  		  equal = false;
  	  }
  	  return equal;
  }

  bool CheckTopicExistence(std::string topic2check)
  {
	  bool success = false;
	  //get all topics in the system running
	  ros::master::V_TopicInfo master_topics;
	  ros::master::getTopics(master_topics);
	  for (ros::master::V_TopicInfo::iterator it = master_topics.begin() ; it != master_topics.end(); it++) {
	      ros::master::TopicInfo& info = *it;
	      if(info.name == topic2check)
	      {
	    	  success = true;
	    	  break;
	      }
	  }

	  return success;
  }

// End namespace "basic_side_functions"
}

namespace geometry_side_functions
{
  //Other namespaces:
  namespace bsf = basic_side_functions;

  geometry_msgs::Vector3 makeVector3(double x, double y, double z)
  {
	  geometry_msgs::Vector3 vector;
	  vector.x = x;
	  vector.y = y;
	  vector.z = z;

	  return vector;
  }

  geometry_msgs::Quaternion makeQuat(double w, double x, double y, double z)
  {
	  geometry_msgs::Quaternion quat;
	  quat.w = w;
	  quat.x = x;
	  quat.y = y;
	  quat.z = z;

	  return quat;
  }

  /* Angles names:
   * Tilt - Bank - Roll - Psi - around x
   * Elevation - Attitude - Pitch - Phi - around Y
   * Azimuth - Heading - Yaw - Theta - around Z
  */
  geometry_msgs::Quaternion RPY2Quat(geometry_msgs::Vector3 RPY, bool RPY_rad, bool turn_corr)
  {
	  geometry_msgs::Quaternion Quat;
	  const int arr_len = 3;

	  double v[arr_len] = {0,0,0};
	  v[0] = RPY.x; v[1] = RPY.y; v[2] = RPY.z;

	  // convert degree to radiant
	  if (!RPY_rad){
		  for (int i = 0; i < arr_len; i++)
		  {
			  v[i] = bsf::deg2rad(v[i]);
		  }
	  }

	  // keep the angles values in (-pi ; pi]
	  if (turn_corr)
	  {
		  double pi = boost::math::constants::pi<double>();
		  for (int i = 0; i < arr_len; i++)
		  {
			  if (v[i] == -pi){
				  v[i] = pi;
			  }
			  // absolute converter
			  if (bsf::abs_f(v[i]) > pi){
				  v[i] = v[i] - bsf::round_f(v[i]/(2*pi))*(2*pi);
			  }
		  }
	  }

	  /*
	  // pre-definition
	  double c1 = cos(v[2]/2);
	  double c2 = cos(v[1]/2);
	  double c3 = cos(v[0]/2);
	  double s1 = sin(v[2]/2);
	  double s2 = sin(v[1]/2);
	  double s3 = sin(v[0]/2);

	  Quat.w = c1*c2*c3 - s1*s2*s3;
	  Quat.x = s1*s2*c3 + c1*c2*s3;
	  Quat.y = c1*s2*c3 - s1*c2*s3;
	  Quat.z = s1*c2*c3 + c1*s2*s3;
    */

	  tf::Quaternion t;
	  t.setRPY(v[0], v[1], v[2]);

	  Quat.w = t.getW();
	  Quat.x = t.getX();
	  Quat.y = t.getY();
	  Quat.z = t.getZ();

	  return Quat;
  }

  geometry_msgs::Vector3 Quat2RPY(geometry_msgs::Quaternion Quat, bool RPY_rad)
  {
	  geometry_msgs::Vector3 RPY2return;

	  KDL::Rotation k;
	  tf::Quaternion t;
	  t.setW(Quat.w); t.setX(Quat.x); t.setY(Quat.y); t.setZ(Quat.z);

	  tf::QuaternionTFToKDL(t,k);

	  k.GetRPY(RPY2return.x, RPY2return.y, RPY2return.z);

	  // convert radiant to degree
	  if (!RPY_rad){
		  RPY2return.x = bsf::rad2deg(RPY2return.x);
		  RPY2return.y = bsf::rad2deg(RPY2return.y);
		  RPY2return.z = bsf::rad2deg(RPY2return.z);
	  }

	  return RPY2return;
  }

  geometry_msgs::Pose makePose(geometry_msgs::Quaternion orientation, geometry_msgs::Vector3 XYZ_location)
  {
	  geometry_msgs::Pose pose;
	  pose.orientation = orientation;
	  pose.position.x = XYZ_location.x; pose.position.y = XYZ_location.y; pose.position.z = XYZ_location.z;

	  return pose;
  }
  //RPY_orientation in degree
  geometry_msgs::Pose makePose(geometry_msgs::Vector3 RPY_orientation, geometry_msgs::Vector3 XYZ_location)
  {
	  geometry_msgs::Pose pose;
	  pose.orientation = RPY2Quat(RPY_orientation, false);
	  pose.position.x = XYZ_location.x; pose.position.y = XYZ_location.y; pose.position.z = XYZ_location.z;

	  return pose;
  }

  geometry_msgs::PoseStamped Pose2PoseStamped( geometry_msgs::Pose old_pose)
  {
	  geometry_msgs::PoseStamped new_pose_s;
	  new_pose_s.pose.position.x = old_pose.position.x;
	  new_pose_s.pose.position.y = old_pose.position.y;
	  new_pose_s.pose.position.z = old_pose.position.z;
	  new_pose_s.pose.orientation.w = old_pose.orientation.w;
	  new_pose_s.pose.orientation.x = old_pose.orientation.x;
	  new_pose_s.pose.orientation.y = old_pose.orientation.y;
	  new_pose_s.pose.orientation.z = old_pose.orientation.z;

	  return new_pose_s;
  }

  std::vector<geometry_msgs::PoseStamped> Pose2PoseStamped( std::vector<geometry_msgs::Pose> old_pose)
  {
	  std::vector<geometry_msgs::PoseStamped> new_pose_s;
	  new_pose_s.resize(old_pose.size());

	  for(int i = 0; i < old_pose.size(); i++)
	  {
		  new_pose_s[i].pose.position.x = old_pose[i].position.x;
		  new_pose_s[i].pose.position.y = old_pose[i].position.y;
		  new_pose_s[i].pose.position.z = old_pose[i].position.z;
		  new_pose_s[i].pose.orientation.w = old_pose[i].orientation.w;
		  new_pose_s[i].pose.orientation.x = old_pose[i].orientation.x;
		  new_pose_s[i].pose.orientation.y = old_pose[i].orientation.y;
		  new_pose_s[i].pose.orientation.z = old_pose[i].orientation.z;
	  }
	  return new_pose_s;
  }

  geometry_msgs::Pose PoseStamped2Pose( geometry_msgs::PoseStamped old_pose_s)
  {
	  geometry_msgs::Pose new_pose;
	  new_pose.position.x = old_pose_s.pose.position.x;
	  new_pose.position.y = old_pose_s.pose.position.y;
	  new_pose.position.z = old_pose_s.pose.position.z;
	  new_pose.orientation.w = old_pose_s.pose.orientation.w;
	  new_pose.orientation.x = old_pose_s.pose.orientation.x;
	  new_pose.orientation.y = old_pose_s.pose.orientation.y;
	  new_pose.orientation.z = old_pose_s.pose.orientation.z;

	  return new_pose;
  }

  std::vector<geometry_msgs::Pose> PoseStamped2Pose( std::vector<geometry_msgs::PoseStamped> old_pose_s)
  {
	  std::vector<geometry_msgs::Pose> new_pose;
	  new_pose.resize(old_pose_s.size());

	  for(int i = 0; i < old_pose_s.size(); i++)
	  {
		  new_pose[i].position.x = old_pose_s[i].pose.position.x;
		  new_pose[i].position.y = old_pose_s[i].pose.position.y;
		  new_pose[i].position.z = old_pose_s[i].pose.position.z;
		  new_pose[i].orientation.w = old_pose_s[i].pose.orientation.w;
		  new_pose[i].orientation.x = old_pose_s[i].pose.orientation.x;
		  new_pose[i].orientation.y = old_pose_s[i].pose.orientation.y;
		  new_pose[i].orientation.z = old_pose_s[i].pose.orientation.z;
	  }
	  return new_pose;
  }

  void vector32posePosition(geometry_msgs::Vector3 vector, geometry_msgs::Pose &pose2update)
  {
	  pose2update.position.x = vector.x;
	  pose2update.position.y = vector.y;
	  pose2update.position.z = vector.z;
  }
  void vector32posePosition(geometry_msgs::Vector3 vector, geometry_msgs::PoseStamped &poseStamped2update)
  {
	  poseStamped2update.pose.position.x = vector.x;
	  poseStamped2update.pose.position.y = vector.y;
	  poseStamped2update.pose.position.z = vector.z;
  }

  void posePosition2vector3(geometry_msgs::Pose pose, geometry_msgs::Vector3 &vector2update)
  {
	  vector2update.x = pose.position.x;
	  vector2update.y = pose.position.y;
	  vector2update.z = pose.position.z;
  }
  void posePosition2vector3(geometry_msgs::PoseStamped poseStamped, geometry_msgs::Vector3 &vector2update)
  {
	  vector2update.x = poseStamped.pose.position.x;
	  vector2update.y = poseStamped.pose.position.y;
	  vector2update.z = poseStamped.pose.position.z;
  }

  double Quat_diff(geometry_msgs::Quaternion a, geometry_msgs::Quaternion b)
  {
	  const int val_quat = 4;
	  double temp[val_quat];
	  double sum = 0.0;

	  temp[0] = bsf::abs_f(a.x - b.x);
	  temp[1] = bsf::abs_f(a.y - b.y);
	  temp[2] = bsf::abs_f(a.z - b.z);
	  temp[3] = bsf::abs_f(a.w - b.w);

	  for (int i; i < val_quat; i++)
	  {
		  sum += temp[i];
	  }

	  return sum;
  }

  bool PoseEquivalence_decimal_value(geometry_msgs::Pose A, geometry_msgs::Pose B, int decimal_considered)
  {
	  bool equal = false;
	  if(decimal_considered < 0)
	  {
		  if(A.position.x == B.position.x &&
			 A.position.y == B.position.y &&
			 A.position.z == B.position.z &&
			 A.orientation.x == B.orientation.x &&
			 A.orientation.y == B.orientation.y &&
			 A.orientation.z == B.orientation.z &&
			 A.orientation.w == B.orientation.w)

			  equal = true;
	  }else{
		  int dc = decimal_considered;
		  if(bsf::decimal_shift(A.position.x, dc) == bsf::decimal_shift(B.position.x, dc) &&
			 bsf::decimal_shift(A.position.y, dc) == bsf::decimal_shift(B.position.y, dc) &&
			 bsf::decimal_shift(A.position.z, dc) == bsf::decimal_shift(B.position.z, dc) &&
			 bsf::decimal_shift(A.orientation.x, dc) == bsf::decimal_shift(B.orientation.x, dc) &&
			 bsf::decimal_shift(A.orientation.y, dc) == bsf::decimal_shift(B.orientation.y, dc) &&
			 bsf::decimal_shift(A.orientation.z, dc) == bsf::decimal_shift(B.orientation.z, dc) &&
			 bsf::decimal_shift(A.orientation.w, dc) == bsf::decimal_shift(B.orientation.w, dc))

			  equal = true;
	  }

	  return equal;
  }

  bool PoseEquivalence_XYZ(geometry_msgs::Pose A, geometry_msgs::Pose B, double threshold_XYZ)
  {
	  bool equal = false;
	  if(threshold_XYZ < 0.0)
	  {
		 std::cout << "Warning: The XYZ threshold must be positive (it has been converted to the positive value)" << std::endl;
		 threshold_XYZ = bsf::abs_f(threshold_XYZ);
	  }

	  double point_dist = std::sqrt(std::pow((A.position.x - B.position.x),2) + std::pow((A.position.y - B.position.y),2) + std::pow((A.position.z - B.position.z),2));
	  if(point_dist <= threshold_XYZ)
		  equal = true;

	  return equal;
  }
  bool PoseEquivalence_Quat(geometry_msgs::Pose A, geometry_msgs::Pose B, double threshold_Quat)
  {
	  bool equal = false;
	  if(threshold_Quat < 0.0)
	  {
		 std::cout << "Warning: The quaternion threshold must be positive (it has been converted to the positive value)" << std::endl;
		 threshold_Quat = bsf::abs_f(threshold_Quat);
	  }

	  double orientation_dist = std::sqrt(std::pow((A.orientation.x - B.orientation.x),2) + std::pow((A.orientation.y - B.orientation.y),2) + std::pow((A.orientation.z - B.orientation.z),2) + std::pow((A.orientation.w - B.orientation.w),2));
	  if(orientation_dist <= threshold_Quat)
		  equal = true;

	  return equal;
  }
  bool PoseEquivalence_theshold(geometry_msgs::Pose A, geometry_msgs::Pose B, double threshold_XYZ, double threshold_Quat)
  {
	  bool equal = false;
	  if(threshold_XYZ < 0.0)
	  {
		 std::cout << "Warning: The XYZ threshold must be positive (it has been converted to the positive value)" << std::endl;
		 threshold_XYZ = bsf::abs_f(threshold_XYZ);
	  }
	  if(threshold_Quat < 0.0)
	  {
		 std::cout << "Warning: The XYZ threshold must be positive (it has been converted to the positive value)" << std::endl;
		 threshold_Quat = bsf::abs_f(threshold_Quat);
	  }

	  if(PoseEquivalence_XYZ(A,B,threshold_XYZ) && PoseEquivalence_Quat(A,B,threshold_Quat))
		  equal = true;

	  return equal;
  }

  std::vector <std::vector <double> > point2matrix_gen(geometry_msgs::Vector3 point)
    {
  	  //initialization
  	  std::vector <std::vector <double> > matrixP; matrixP.resize(4);
  	  std::vector <double> temp_line; temp_line.resize(1);
  	  for(int i = 0; i < matrixP.size(); i++)
  		  matrixP[i] = temp_line;

  	  //main prog
  	  matrixP[0][0] = point.x;
  	  matrixP[1][0] = point.y;
  	  matrixP[2][0] = point.z;
  	  matrixP[3][0] = 1.0;

  	  return matrixP;
    }

  geometry_msgs::Vector3 matrix2point_gen(std::vector <std::vector <double> > matrix)
    {
	  geometry_msgs::Vector3 point;

  	  //initialization check
	  if(matrix.size() != 4)
	  {
		  std::cout << "Error: the number of raws are supposed to be 4!" << std::endl;
		  std::cout << "A default null vector has been returned..." << std::endl;
		  return point;
	  }
	  for(int i = 0; i < 4; i++)
	  {
		  if(matrix[i].size() != 1)
		  {
			  std::cout << "Error: the raw number " << i << " is supposed to be composed by just 1 element!" << std::endl;
			  std::cout << "A default null vector has been returned..." << std::endl;
			  return point;
		  }
	  }
	  if(matrix[3][0] != 1.0)
	  {
		  std::cout << "Error: the last element is supposed to be equal to 1.0!" << std::endl;
		  std::cout << "A default null vector has been returned..." << std::endl;
		  return point;
	  }

  	  //main prog
  	  point.x = matrix[0][0];
  	  point.y = matrix[1][0];
  	  point.z = matrix[2][0];

  	  return point;
    }

  std::vector <std::vector <double> > translation_matrix(geometry_msgs::Vector3 translation)
  {
	  //initialization
	  std::vector <std::vector <double> > matrixT; matrixT.resize(4);
	  std::vector <double> temp_line; temp_line.resize(4);
	  for(int i = 0; i < matrixT.size(); i++)
		  matrixT[i] = temp_line;

	  //main prog
	  matrixT[0][0] = 1.0;
	  matrixT[1][1] = 1.0;
	  matrixT[2][2] = 1.0;
	  matrixT[3][3] = 1.0;
	  matrixT[0][3] = translation.x;
	  matrixT[1][3] = translation.y;
	  matrixT[2][3] = translation.z;

	  return matrixT;
  }

  std::vector <std::vector <double> > rotational_matrix_X(double angle_X, bool rad)
  {
	  //initialization
	  std::vector <std::vector <double> > matrixX; matrixX.resize(4);
	  std::vector <double> temp_line; temp_line.resize(4);
	  for(int i = 0; i < matrixX.size(); i++)
		  matrixX[i] = temp_line;

	  //conversion to radiant (if needed)
	  if(!rad)
		  angle_X = bsf::deg2rad(angle_X);

	  //main prog
	  matrixX[0][0] = 1.0;
	  matrixX[3][3] = 1.0;
	  matrixX[1][1] = cos(angle_X);
	  matrixX[2][2] = cos(angle_X);
	  matrixX[1][2] = -sin(angle_X);
	  matrixX[2][1] = sin(angle_X);

	  return matrixX;
  }

  std::vector <std::vector <double> > rotational_matrix_Y(double angle_Y, bool rad)
  {
	  //initialization
	  std::vector <std::vector <double> > matrixY; matrixY.resize(4);
	  std::vector <double> temp_line; temp_line.resize(4);
	  for(int i = 0; i < matrixY.size(); i++)
		  matrixY[i] = temp_line;

	  //conversion to radiant (if needed)
	  if(!rad)
		  angle_Y = bsf::deg2rad(angle_Y);

	  //main prog
	  matrixY[1][1] = 1.0;
	  matrixY[3][3] = 1.0;
	  matrixY[0][0] = cos(angle_Y);
	  matrixY[2][2] = cos(angle_Y);
	  matrixY[0][2] = sin(angle_Y);
	  matrixY[2][0] = -sin(angle_Y);

	  return matrixY;
  }

  std::vector <std::vector <double> > rotational_matrix_Z(double angle_Z, bool rad)
  {
	  //initialization
	  std::vector <std::vector <double> > matrixZ; matrixZ.resize(4);
	  std::vector <double> temp_line; temp_line.resize(4);
	  for(int i = 0; i < matrixZ.size(); i++)
		  matrixZ[i] = temp_line;

	  //conversion to radiant (if needed)
	  if(!rad)
		  angle_Z = bsf::deg2rad(angle_Z);

	  //main prog
	  matrixZ[2][2] = 1.0;
	  matrixZ[3][3] = 1.0;
	  matrixZ[0][0] = cos(angle_Z);
	  matrixZ[1][1] = cos(angle_Z);
	  matrixZ[0][1] = -sin(angle_Z);
	  matrixZ[1][0] = sin(angle_Z);

	  return matrixZ;
  }

  std::vector <std::vector <double> > RPY2rotMatrix(geometry_msgs::Vector3 angles, bool rad)
  {
	  //initialization
	  std::vector <std::vector <double> > matrixR; matrixR.resize(4);
	  std::vector <double> temp_line; temp_line.resize(4);
	  for(int i = 0; i <  matrixR.size(); i++)
		  matrixR[i] = temp_line;

	  //conversion to radiant (if needed)
	  if(!rad)
	  {
		  angles.x = bsf::deg2rad(angles.x);
		  angles.y = bsf::deg2rad(angles.y);
		  angles.z = bsf::deg2rad(angles.z);
	  }

	  //main prog
	  std::vector <std::vector <double> > matrixX = rotational_matrix_X(angles.x);
	  std::vector <std::vector <double> > matrixY = rotational_matrix_Y(angles.y);
	  std::vector <std::vector <double> > matrixZ = rotational_matrix_Z(angles.z);

	  bsf::matrix2Dprod(matrixZ,matrixY,matrixR);
	  bsf::matrix2Dprod(matrixR,matrixX,matrixR);

	  return matrixR;
  }
  std::vector <std::vector <double> > Quaternion2rotMatrix(geometry_msgs::Quaternion quat)
  {
	  //initialization
	  std::vector <std::vector <double> > matrixR; matrixR.resize(4);
	  std::vector <double> temp_line; temp_line.resize(4);
	  for(int i = 0; i < matrixR.size(); i++)
		  matrixR[i] = temp_line;

	  geometry_msgs::Vector3 angles = Quat2RPY(quat);

	  //main prog
	  std::vector <std::vector <double> > matrixX = rotational_matrix_X(angles.x);
	  std::vector <std::vector <double> > matrixY = rotational_matrix_Y(angles.y);
	  std::vector <std::vector <double> > matrixZ = rotational_matrix_Z(angles.z);

	  bsf::matrix2Dprod(matrixZ,matrixY,matrixR);
	  bsf::matrix2Dprod(matrixR,matrixX,matrixR);

	  return matrixR;
  }

  geometry_msgs::Vector3 rotMatrix2RPY( std::vector <std::vector <double> > matrix, bool rad)
  {
	  geometry_msgs::Vector3 rot_vector;

  	  //initialization check
	  if(matrix.size() != 4)
	  {
		  std::cout << "Error: the number of raws are supposed to be 4!" << std::endl;
		  std::cout << "A default null vector has been returned..." << std::endl;
		  return rot_vector;
	  }
	  for(int i = 0; i < 4; i++)
	  {
		  if(matrix[i].size() != 4)
		  {
			  std::cout << "Error: the raw number " << i << " is supposed to be composed by just 1 element!" << std::endl;
			  std::cout << "A default null vector has been returned..." << std::endl;
			  return rot_vector;
		  }
	  }
	  if(matrix[3][3] != 1.0)
	  {
		  std::cout << "Error: the last element is supposed to be equal to 1.0!" << std::endl;
		  std::cout << "A default null vector has been returned..." << std::endl;
		  return rot_vector;
	  }
	  for(int i = 0; i < 3; i++)
	  {
		  if(matrix[3][i] != 0.0 || matrix[i][3] != 0.0)
		  {
			  std::cout << "Error: the values contained in the cells [" << i << ",3] and [3," << i << "] is supposed to be ZERO!" << std::endl;
			  std::cout << "A default null vector has been returned..." << std::endl;
			  return rot_vector;
		  }
	  }

  	  //main prog
	  rot_vector.x = atan2(matrix[2][1], matrix[2][2]);
	  rot_vector.y = atan2(-matrix[2][0], sqrt((matrix[2][1]*matrix[2][1]) + (matrix[2][2]*matrix[2][2])));
	  rot_vector.z = atan2(matrix[1][0], matrix[0][0]);

	  //conversion to radiant (if needed)
	  if(!rad)
	  {
		  rot_vector.x = bsf::rad2deg(rot_vector.x);
		  rot_vector.y = bsf::rad2deg(rot_vector.y);
		  rot_vector.z = bsf::rad2deg(rot_vector.z);
	  }

  	  return rot_vector;
  }

  geometry_msgs::Quaternion rotMatrix2Quaternion(std::vector<std::vector <double> > matrix)
  {
	  geometry_msgs::Vector3 rot_vector;
	  geometry_msgs::Quaternion quat_vec;

  	  //initialization check
	  if(matrix.size() != 4)
	  {
		  std::cout << "Error: the number of raws are supposed to be 4!" << std::endl;
		  std::cout << "A default null vector has been returned..." << std::endl;
		  return quat_vec;
	  }
	  for(int i = 0; i < 4; i++)
	  {
		  if(matrix[i].size() != 4)
		  {
			  std::cout << "Error: the raw number " << i << " is supposed to be composed by just 1 element!" << std::endl;
			  std::cout << "A default null vector has been returned..." << std::endl;
			  return quat_vec;
		  }
	  }
	  if(matrix[3][3] != 1.0)
	  {
		  std::cout << "Error: the last element is supposed to be equal to 1.0!" << std::endl;
		  std::cout << "A default null vector has been returned..." << std::endl;
		  return quat_vec;
	  }
	  for(int i = 0; i < 3; i++)
	  {
		  if(matrix[3][i] != 0.0 || matrix[i][3] != 0.0)
		  {
			  std::cout << "Error: the values contained in the cells [" << i <<",3] and [3," << i << "] is supposed to be ZERO!" << std::endl;
			  std::cout << "A default null vector has been returned..." << std::endl;
			  return quat_vec;
		  }
	  }

  	  //main prog
	  rot_vector.x = atan2(matrix[2][1], matrix[2][2]);
	  rot_vector.y = atan2(-matrix[2][0], sqrt((matrix[2][1]*matrix[2][1]) + (matrix[2][2]*matrix[2][2])));
	  rot_vector.z = atan2(matrix[1][0], matrix[0][0]);

	  //conversion to quaternions
	  quat_vec = RPY2Quat(rot_vector);

	  return quat_vec;
  }

// End namespace "geometry_side_functions"
}

