PROJECT(SRF-Odometry)

CMAKE_MINIMUM_REQUIRED(VERSION 2.4)
if(COMMAND cmake_policy)
      cmake_policy(SET CMP0003 NEW)  # Required by CMake 2.7+
endif(COMMAND cmake_policy)


FIND_PACKAGE(MRPT REQUIRED base gui opengl nav obs maps)



ADD_LIBRARY(srf_lib
	laser_odometry_v1.cpp
	laser_odometry_v1.h
	laser_odometry_standard.cpp
	laser_odometry_standard.h
	laser_odometry_nosym.cpp
	laser_odometry_nosym.h
	laser_odometry_3scans.cpp
	laser_odometry_3scans.h
	laser_odometry_refscans.cpp
	laser_odometry_refscans.h
)



ADD_EXECUTABLE(Laser-odometry-randomnav 
	main_laserodo_randomnav.cpp
	laserodo_randomnav.h
	map.xpm
	map_lines_rf2o.xpm
	map_lab_rf2o.xpm
	map_lab.xpm
	)
	
TARGET_LINK_LIBRARIES(Laser-odometry-randomnav
		${MRPT_LIBS}
		srf_lib)




ADD_EXECUTABLE(Rawlog-groundtruth  
	main_rawlog_gt.cpp
	rawlog_gt.h
	polar_match.h
	polar_match.cpp
	)
	
TARGET_LINK_LIBRARIES(Rawlog-groundtruth 
		${MRPT_LIBS})		

 
# Set optimized building:
IF(CMAKE_COMPILER_IS_GNUCXX AND NOT CMAKE_BUILD_TYPE MATCHES "Debug")
	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -mtune=native")
ENDIF(CMAKE_COMPILER_IS_GNUCXX AND NOT CMAKE_BUILD_TYPE MATCHES "Debug")

