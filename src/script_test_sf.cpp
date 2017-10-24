 /*Author: Emanuele Sansebastiano */

// libs
#include <ros/ros.h>
#include <side_pkg/side_func.h>

int main(int argc, char **argv)
{
	ros::init(argc, argv, "script_test_sf");
	ros::NodeHandle node_handle("~");

	ros::AsyncSpinner spinner(1);

	spinner.start();

	namespace bsc = basic_side_classes;
	namespace bsf = basic_side_functions;
	namespace gsf = geometry_side_functions;

	bsf::countdown_sec(3);
	//real program init


	ros::shutdown();
	return 0;
}
