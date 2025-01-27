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
include src/ray_tracing/math_formula/CMakeFiles/math_formula.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include src/ray_tracing/math_formula/CMakeFiles/math_formula.dir/compiler_depend.make

# Include the progress variables for this target.
include src/ray_tracing/math_formula/CMakeFiles/math_formula.dir/progress.make

# Include the compile flags for this target's objects.
include src/ray_tracing/math_formula/CMakeFiles/math_formula.dir/flags.make

src/ray_tracing/math_formula/CMakeFiles/math_formula.dir/math_formula.cpp.o: src/ray_tracing/math_formula/CMakeFiles/math_formula.dir/flags.make
src/ray_tracing/math_formula/CMakeFiles/math_formula.dir/math_formula.cpp.o: ../src/ray_tracing/math_formula/math_formula.cpp
src/ray_tracing/math_formula/CMakeFiles/math_formula.dir/math_formula.cpp.o: src/ray_tracing/math_formula/CMakeFiles/math_formula.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/starman/Desktop/black_hole/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object src/ray_tracing/math_formula/CMakeFiles/math_formula.dir/math_formula.cpp.o"
	cd /Users/starman/Desktop/black_hole/build/src/ray_tracing/math_formula && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT src/ray_tracing/math_formula/CMakeFiles/math_formula.dir/math_formula.cpp.o -MF CMakeFiles/math_formula.dir/math_formula.cpp.o.d -o CMakeFiles/math_formula.dir/math_formula.cpp.o -c /Users/starman/Desktop/black_hole/src/ray_tracing/math_formula/math_formula.cpp

src/ray_tracing/math_formula/CMakeFiles/math_formula.dir/math_formula.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/math_formula.dir/math_formula.cpp.i"
	cd /Users/starman/Desktop/black_hole/build/src/ray_tracing/math_formula && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/starman/Desktop/black_hole/src/ray_tracing/math_formula/math_formula.cpp > CMakeFiles/math_formula.dir/math_formula.cpp.i

src/ray_tracing/math_formula/CMakeFiles/math_formula.dir/math_formula.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/math_formula.dir/math_formula.cpp.s"
	cd /Users/starman/Desktop/black_hole/build/src/ray_tracing/math_formula && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/starman/Desktop/black_hole/src/ray_tracing/math_formula/math_formula.cpp -o CMakeFiles/math_formula.dir/math_formula.cpp.s

# Object files for target math_formula
math_formula_OBJECTS = \
"CMakeFiles/math_formula.dir/math_formula.cpp.o"

# External object files for target math_formula
math_formula_EXTERNAL_OBJECTS =

src/ray_tracing/math_formula/libmath_formula.a: src/ray_tracing/math_formula/CMakeFiles/math_formula.dir/math_formula.cpp.o
src/ray_tracing/math_formula/libmath_formula.a: src/ray_tracing/math_formula/CMakeFiles/math_formula.dir/build.make
src/ray_tracing/math_formula/libmath_formula.a: src/ray_tracing/math_formula/CMakeFiles/math_formula.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/starman/Desktop/black_hole/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX static library libmath_formula.a"
	cd /Users/starman/Desktop/black_hole/build/src/ray_tracing/math_formula && $(CMAKE_COMMAND) -P CMakeFiles/math_formula.dir/cmake_clean_target.cmake
	cd /Users/starman/Desktop/black_hole/build/src/ray_tracing/math_formula && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/math_formula.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
src/ray_tracing/math_formula/CMakeFiles/math_formula.dir/build: src/ray_tracing/math_formula/libmath_formula.a
.PHONY : src/ray_tracing/math_formula/CMakeFiles/math_formula.dir/build

src/ray_tracing/math_formula/CMakeFiles/math_formula.dir/clean:
	cd /Users/starman/Desktop/black_hole/build/src/ray_tracing/math_formula && $(CMAKE_COMMAND) -P CMakeFiles/math_formula.dir/cmake_clean.cmake
.PHONY : src/ray_tracing/math_formula/CMakeFiles/math_formula.dir/clean

src/ray_tracing/math_formula/CMakeFiles/math_formula.dir/depend:
	cd /Users/starman/Desktop/black_hole/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/starman/Desktop/black_hole /Users/starman/Desktop/black_hole/src/ray_tracing/math_formula /Users/starman/Desktop/black_hole/build /Users/starman/Desktop/black_hole/build/src/ray_tracing/math_formula /Users/starman/Desktop/black_hole/build/src/ray_tracing/math_formula/CMakeFiles/math_formula.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/ray_tracing/math_formula/CMakeFiles/math_formula.dir/depend

