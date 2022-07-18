# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.23

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

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/local/Cellar/cmake/3.23.0/bin/cmake

# The command to remove a file.
RM = /usr/local/Cellar/cmake/3.23.0/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/starman/Desktop/black_hole

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/starman/Desktop/black_hole/build

# Include any dependencies generated for this target.
include src/ray_tracing/CMakeFiles/ray_tracing.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include src/ray_tracing/CMakeFiles/ray_tracing.dir/compiler_depend.make

# Include the progress variables for this target.
include src/ray_tracing/CMakeFiles/ray_tracing.dir/progress.make

# Include the compile flags for this target's objects.
include src/ray_tracing/CMakeFiles/ray_tracing.dir/flags.make

src/ray_tracing/CMakeFiles/ray_tracing.dir/ray_tracing.cpp.o: src/ray_tracing/CMakeFiles/ray_tracing.dir/flags.make
src/ray_tracing/CMakeFiles/ray_tracing.dir/ray_tracing.cpp.o: ../src/ray_tracing/ray_tracing.cpp
src/ray_tracing/CMakeFiles/ray_tracing.dir/ray_tracing.cpp.o: src/ray_tracing/CMakeFiles/ray_tracing.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/starman/Desktop/black_hole/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object src/ray_tracing/CMakeFiles/ray_tracing.dir/ray_tracing.cpp.o"
	cd /Users/starman/Desktop/black_hole/build/src/ray_tracing && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT src/ray_tracing/CMakeFiles/ray_tracing.dir/ray_tracing.cpp.o -MF CMakeFiles/ray_tracing.dir/ray_tracing.cpp.o.d -o CMakeFiles/ray_tracing.dir/ray_tracing.cpp.o -c /Users/starman/Desktop/black_hole/src/ray_tracing/ray_tracing.cpp

src/ray_tracing/CMakeFiles/ray_tracing.dir/ray_tracing.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/ray_tracing.dir/ray_tracing.cpp.i"
	cd /Users/starman/Desktop/black_hole/build/src/ray_tracing && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/starman/Desktop/black_hole/src/ray_tracing/ray_tracing.cpp > CMakeFiles/ray_tracing.dir/ray_tracing.cpp.i

src/ray_tracing/CMakeFiles/ray_tracing.dir/ray_tracing.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/ray_tracing.dir/ray_tracing.cpp.s"
	cd /Users/starman/Desktop/black_hole/build/src/ray_tracing && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/starman/Desktop/black_hole/src/ray_tracing/ray_tracing.cpp -o CMakeFiles/ray_tracing.dir/ray_tracing.cpp.s

# Object files for target ray_tracing
ray_tracing_OBJECTS = \
"CMakeFiles/ray_tracing.dir/ray_tracing.cpp.o"

# External object files for target ray_tracing
ray_tracing_EXTERNAL_OBJECTS =

src/ray_tracing/libray_tracing.a: src/ray_tracing/CMakeFiles/ray_tracing.dir/ray_tracing.cpp.o
src/ray_tracing/libray_tracing.a: src/ray_tracing/CMakeFiles/ray_tracing.dir/build.make
src/ray_tracing/libray_tracing.a: src/ray_tracing/CMakeFiles/ray_tracing.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/starman/Desktop/black_hole/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX static library libray_tracing.a"
	cd /Users/starman/Desktop/black_hole/build/src/ray_tracing && $(CMAKE_COMMAND) -P CMakeFiles/ray_tracing.dir/cmake_clean_target.cmake
	cd /Users/starman/Desktop/black_hole/build/src/ray_tracing && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/ray_tracing.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
src/ray_tracing/CMakeFiles/ray_tracing.dir/build: src/ray_tracing/libray_tracing.a
.PHONY : src/ray_tracing/CMakeFiles/ray_tracing.dir/build

src/ray_tracing/CMakeFiles/ray_tracing.dir/clean:
	cd /Users/starman/Desktop/black_hole/build/src/ray_tracing && $(CMAKE_COMMAND) -P CMakeFiles/ray_tracing.dir/cmake_clean.cmake
.PHONY : src/ray_tracing/CMakeFiles/ray_tracing.dir/clean

src/ray_tracing/CMakeFiles/ray_tracing.dir/depend:
	cd /Users/starman/Desktop/black_hole/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/starman/Desktop/black_hole /Users/starman/Desktop/black_hole/src/ray_tracing /Users/starman/Desktop/black_hole/build /Users/starman/Desktop/black_hole/build/src/ray_tracing /Users/starman/Desktop/black_hole/build/src/ray_tracing/CMakeFiles/ray_tracing.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/ray_tracing/CMakeFiles/ray_tracing.dir/depend

