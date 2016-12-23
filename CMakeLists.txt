PROJECT(SRF-Odometry)

CMAKE_MINIMUM_REQUIRED(VERSION 2.4)
if(COMMAND cmake_policy)
      cmake_policy(SET CMP0003 NEW)  # Required by CMake 2.7+
endif(COMMAND cmake_policy)


FIND_PACKAGE(MRPT REQUIRED base gui opengl nav obs maps)

FIND_PACKAGE(PCL REQUIRED)
INCLUDE_DIRECTORIES(${PCL_INCLUDE_DIRS})
LINK_DIRECTORIES(${PCL_LIBRARY_DIRS})
ADD_DEFINITIONS(${PCL_DEFINITIONS})

# Require we have pkgconfig installed to find CSM
find_package ( PkgConfig REQUIRED )
pkg_check_modules( CSM REQUIRED csm )

IF (${CSM_FOUND})
	MESSAGE("CSM_LIBRARY_DIRS: ${CSM_LIBRARY_DIRS}")
	MESSAGE("CSM_LIBRARIES: ${CSM_LIBRARIES}")
	MESSAGE("CSM_INCLUDE_DIRS: ${CSM_INCLUDE_DIRS}")
	
	INCLUDE_DIRECTORIES( ${CSM_INCLUDE_DIRS} ) # important !
	LINK_DIRECTORIES( ${CSM_LIBRARY_DIRS} )	 # important !
ELSE (${CSM_FOUND})
	MESSAGE( FATAL_ERROR " CSM not found . Check that the environment variable PKG_CONFIG_PATH includes the path #containing the file 'csm.pc'.")
ENDIF (${CSM_FOUND})


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


ADD_EXECUTABLE(Laser-odometry-comparison  
	main_laserodo_comparison.cpp
	laserodo_comparison.h
	map.xpm
	map_lines_rf2o.xpm
	map_lab_rf2o.xpm
	map_lab.xpm
	polar_match.h
	polar_match.cpp
)
	
TARGET_LINK_LIBRARIES(Laser-odometry-comparison 
		${MRPT_LIBS}
		${CSM_LIBRARIES}
		${PCL_LIBRARIES} 
		srf_lib)


ADD_EXECUTABLE(Laser-odometry-experiments 
	main_laserodo_experiments.cpp
	laserodo_experiments.h
	map.xpm
	map_lines_rf2o.xpm
	map_lab_rf2o.xpm
	map_lab.xpm
	polar_match.h
	polar_match.cpp
)
	
TARGET_LINK_LIBRARIES(Laser-odometry-experiments 
		${MRPT_LIBS}
		${CSM_LIBRARIES}
		srf_lib)


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


ADD_EXECUTABLE(Laser-odometry-rawlog  
	main_laserodo_rawlog.cpp
	laserodo_rawlog.h
	map_lab_big.xpm
	map_lab.xpm
	polar_match.h
	polar_match.cpp
	)
	
TARGET_LINK_LIBRARIES(Laser-odometry-rawlog
		${MRPT_LIBS}
		${CSM_LIBRARIES}
		srf_lib)


ADD_EXECUTABLE(Laser-odometry-targets  
	main_laserodo_targets.cpp
	laserodo_targets.h
	map_lab_big.xpm
	map_lab.xpm
	polar_match.h
	polar_match.cpp
	)
	
TARGET_LINK_LIBRARIES(Laser-odometry-targets
		${MRPT_LIBS}
		${CSM_LIBRARIES}
		srf_lib)


ADD_EXECUTABLE(Laser-odometry-video-raw  
	main_laserodo_raw_video.cpp
	laserodo_raw_video.h
	polar_match.h
	polar_match.cpp
	)
	
TARGET_LINK_LIBRARIES(Laser-odometry-video-raw
		${MRPT_LIBS}
		${CSM_LIBRARIES}
		srf_lib)

ADD_EXECUTABLE(Laser-odometry-video-sim  
	main_laserodo_sim_video.cpp
	laserodo_sim_video.h
	polar_match.h
	polar_match.cpp
	)
	
TARGET_LINK_LIBRARIES(Laser-odometry-video-sim
		${MRPT_LIBS}
		${CSM_LIBRARIES}
		srf_lib)


ADD_EXECUTABLE(Laser-odometry-video-multirobot 
	main_laserodo_video_multirobot.cpp
	laserodo_video_multirobot.h
	polar_match.h
	polar_match.cpp
	)
	
TARGET_LINK_LIBRARIES(Laser-odometry-video-multirobot
		${MRPT_LIBS}
		${CSM_LIBRARIES}
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

