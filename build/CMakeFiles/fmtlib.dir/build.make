# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.17

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Disable VCS-based implicit rules.
% : %,v


# Disable VCS-based implicit rules.
% : RCS/%


# Disable VCS-based implicit rules.
% : RCS/%,v


# Disable VCS-based implicit rules.
% : SCCS/s.%


# Disable VCS-based implicit rules.
% : s.%


.SUFFIXES: .hpux_make_needs_suffix_list


# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/local/Cellar/cmake/3.17.3/bin/cmake

# The command to remove a file.
RM = /usr/local/Cellar/cmake/3.17.3/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/admin/Documents/GitHub/QUIT

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/admin/Documents/GitHub/QUIT/build

# Include any dependencies generated for this target.
include CMakeFiles/fmtlib.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/fmtlib.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/fmtlib.dir/flags.make

CMakeFiles/fmtlib.dir/External/include/fmt/format.cc.o: CMakeFiles/fmtlib.dir/flags.make
CMakeFiles/fmtlib.dir/External/include/fmt/format.cc.o: ../External/include/fmt/format.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/admin/Documents/GitHub/QUIT/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/fmtlib.dir/External/include/fmt/format.cc.o"
	/Library/Developer/CommandLineTools/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/fmtlib.dir/External/include/fmt/format.cc.o -c /Users/admin/Documents/GitHub/QUIT/External/include/fmt/format.cc

CMakeFiles/fmtlib.dir/External/include/fmt/format.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/fmtlib.dir/External/include/fmt/format.cc.i"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/admin/Documents/GitHub/QUIT/External/include/fmt/format.cc > CMakeFiles/fmtlib.dir/External/include/fmt/format.cc.i

CMakeFiles/fmtlib.dir/External/include/fmt/format.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/fmtlib.dir/External/include/fmt/format.cc.s"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/admin/Documents/GitHub/QUIT/External/include/fmt/format.cc -o CMakeFiles/fmtlib.dir/External/include/fmt/format.cc.s

CMakeFiles/fmtlib.dir/External/include/fmt/posix.cc.o: CMakeFiles/fmtlib.dir/flags.make
CMakeFiles/fmtlib.dir/External/include/fmt/posix.cc.o: ../External/include/fmt/posix.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/admin/Documents/GitHub/QUIT/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/fmtlib.dir/External/include/fmt/posix.cc.o"
	/Library/Developer/CommandLineTools/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/fmtlib.dir/External/include/fmt/posix.cc.o -c /Users/admin/Documents/GitHub/QUIT/External/include/fmt/posix.cc

CMakeFiles/fmtlib.dir/External/include/fmt/posix.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/fmtlib.dir/External/include/fmt/posix.cc.i"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/admin/Documents/GitHub/QUIT/External/include/fmt/posix.cc > CMakeFiles/fmtlib.dir/External/include/fmt/posix.cc.i

CMakeFiles/fmtlib.dir/External/include/fmt/posix.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/fmtlib.dir/External/include/fmt/posix.cc.s"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/admin/Documents/GitHub/QUIT/External/include/fmt/posix.cc -o CMakeFiles/fmtlib.dir/External/include/fmt/posix.cc.s

# Object files for target fmtlib
fmtlib_OBJECTS = \
"CMakeFiles/fmtlib.dir/External/include/fmt/format.cc.o" \
"CMakeFiles/fmtlib.dir/External/include/fmt/posix.cc.o"

# External object files for target fmtlib
fmtlib_EXTERNAL_OBJECTS =

libfmtlib.a: CMakeFiles/fmtlib.dir/External/include/fmt/format.cc.o
libfmtlib.a: CMakeFiles/fmtlib.dir/External/include/fmt/posix.cc.o
libfmtlib.a: CMakeFiles/fmtlib.dir/build.make
libfmtlib.a: CMakeFiles/fmtlib.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/admin/Documents/GitHub/QUIT/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Linking CXX static library libfmtlib.a"
	$(CMAKE_COMMAND) -P CMakeFiles/fmtlib.dir/cmake_clean_target.cmake
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/fmtlib.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/fmtlib.dir/build: libfmtlib.a

.PHONY : CMakeFiles/fmtlib.dir/build

CMakeFiles/fmtlib.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/fmtlib.dir/cmake_clean.cmake
.PHONY : CMakeFiles/fmtlib.dir/clean

CMakeFiles/fmtlib.dir/depend:
	cd /Users/admin/Documents/GitHub/QUIT/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/admin/Documents/GitHub/QUIT /Users/admin/Documents/GitHub/QUIT /Users/admin/Documents/GitHub/QUIT/build /Users/admin/Documents/GitHub/QUIT/build /Users/admin/Documents/GitHub/QUIT/build/CMakeFiles/fmtlib.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/fmtlib.dir/depend
