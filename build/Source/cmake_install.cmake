# Install script for directory: /Users/admin/Documents/GitHub/QUIT/Source

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Release")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/Users/admin/Documents/GitHub/QUIT/build/Source/qi")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/qi" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/qi")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/Library/Developer/CommandLineTools/usr/bin/strip" -u -r "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/qi")
    endif()
  endif()
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("/Users/admin/Documents/GitHub/QUIT/build/Source/Core/cmake_install.cmake")
  include("/Users/admin/Documents/GitHub/QUIT/build/Source/ImageIO/cmake_install.cmake")
  include("/Users/admin/Documents/GitHub/QUIT/build/Source/Sequences/cmake_install.cmake")
  include("/Users/admin/Documents/GitHub/QUIT/build/Source/CoreProgs/cmake_install.cmake")
  include("/Users/admin/Documents/GitHub/QUIT/build/Source/B1/cmake_install.cmake")
  include("/Users/admin/Documents/GitHub/QUIT/build/Source/MT/cmake_install.cmake")
  include("/Users/admin/Documents/GitHub/QUIT/build/Source/Perfusion/cmake_install.cmake")
  include("/Users/admin/Documents/GitHub/QUIT/build/Source/Relaxometry/cmake_install.cmake")
  include("/Users/admin/Documents/GitHub/QUIT/build/Source/RUFIS/cmake_install.cmake")
  include("/Users/admin/Documents/GitHub/QUIT/build/Source/Stats/cmake_install.cmake")
  include("/Users/admin/Documents/GitHub/QUIT/build/Source/Susceptibility/cmake_install.cmake")
  include("/Users/admin/Documents/GitHub/QUIT/build/Source/Utils/cmake_install.cmake")

endif()

