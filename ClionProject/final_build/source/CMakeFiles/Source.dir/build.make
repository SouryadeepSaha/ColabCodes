# CMAKE generated file: DO NOT EDIT!
# Generated by "MinGW Makefiles" Generator, CMake Version 3.25

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

SHELL = cmd.exe

# The CMake executable.
CMAKE_COMMAND = "C:\Program Files\CMake\bin\cmake.exe"

# The command to remove a file.
RM = "C:\Program Files\CMake\bin\cmake.exe" -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = C:\Users\Saha\Documents\GitHub\ColabCodes\ClionProject

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = C:\Users\Saha\Documents\GitHub\ColabCodes\ClionProject\final_build

# Include any dependencies generated for this target.
include source/CMakeFiles/Source.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include source/CMakeFiles/Source.dir/compiler_depend.make

# Include the progress variables for this target.
include source/CMakeFiles/Source.dir/progress.make

# Include the compile flags for this target's objects.
include source/CMakeFiles/Source.dir/flags.make

source/CMakeFiles/Source.dir/Grid.cpp.obj: source/CMakeFiles/Source.dir/flags.make
source/CMakeFiles/Source.dir/Grid.cpp.obj: C:/Users/Saha/Documents/GitHub/ColabCodes/ClionProject/source/Grid.cpp
source/CMakeFiles/Source.dir/Grid.cpp.obj: source/CMakeFiles/Source.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=C:\Users\Saha\Documents\GitHub\ColabCodes\ClionProject\final_build\CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object source/CMakeFiles/Source.dir/Grid.cpp.obj"
	cd /d C:\Users\Saha\Documents\GitHub\ColabCodes\ClionProject\final_build\source && C:\MinGW\bin\g++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT source/CMakeFiles/Source.dir/Grid.cpp.obj -MF CMakeFiles\Source.dir\Grid.cpp.obj.d -o CMakeFiles\Source.dir\Grid.cpp.obj -c C:\Users\Saha\Documents\GitHub\ColabCodes\ClionProject\source\Grid.cpp

source/CMakeFiles/Source.dir/Grid.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/Source.dir/Grid.cpp.i"
	cd /d C:\Users\Saha\Documents\GitHub\ColabCodes\ClionProject\final_build\source && C:\MinGW\bin\g++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E C:\Users\Saha\Documents\GitHub\ColabCodes\ClionProject\source\Grid.cpp > CMakeFiles\Source.dir\Grid.cpp.i

source/CMakeFiles/Source.dir/Grid.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/Source.dir/Grid.cpp.s"
	cd /d C:\Users\Saha\Documents\GitHub\ColabCodes\ClionProject\final_build\source && C:\MinGW\bin\g++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S C:\Users\Saha\Documents\GitHub\ColabCodes\ClionProject\source\Grid.cpp -o CMakeFiles\Source.dir\Grid.cpp.s

source/CMakeFiles/Source.dir/Poisson_FD.cpp.obj: source/CMakeFiles/Source.dir/flags.make
source/CMakeFiles/Source.dir/Poisson_FD.cpp.obj: C:/Users/Saha/Documents/GitHub/ColabCodes/ClionProject/source/Poisson_FD.cpp
source/CMakeFiles/Source.dir/Poisson_FD.cpp.obj: source/CMakeFiles/Source.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=C:\Users\Saha\Documents\GitHub\ColabCodes\ClionProject\final_build\CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object source/CMakeFiles/Source.dir/Poisson_FD.cpp.obj"
	cd /d C:\Users\Saha\Documents\GitHub\ColabCodes\ClionProject\final_build\source && C:\MinGW\bin\g++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT source/CMakeFiles/Source.dir/Poisson_FD.cpp.obj -MF CMakeFiles\Source.dir\Poisson_FD.cpp.obj.d -o CMakeFiles\Source.dir\Poisson_FD.cpp.obj -c C:\Users\Saha\Documents\GitHub\ColabCodes\ClionProject\source\Poisson_FD.cpp

source/CMakeFiles/Source.dir/Poisson_FD.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/Source.dir/Poisson_FD.cpp.i"
	cd /d C:\Users\Saha\Documents\GitHub\ColabCodes\ClionProject\final_build\source && C:\MinGW\bin\g++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E C:\Users\Saha\Documents\GitHub\ColabCodes\ClionProject\source\Poisson_FD.cpp > CMakeFiles\Source.dir\Poisson_FD.cpp.i

source/CMakeFiles/Source.dir/Poisson_FD.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/Source.dir/Poisson_FD.cpp.s"
	cd /d C:\Users\Saha\Documents\GitHub\ColabCodes\ClionProject\final_build\source && C:\MinGW\bin\g++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S C:\Users\Saha\Documents\GitHub\ColabCodes\ClionProject\source\Poisson_FD.cpp -o CMakeFiles\Source.dir\Poisson_FD.cpp.s

source/CMakeFiles/Source.dir/Poisson_GMG.cpp.obj: source/CMakeFiles/Source.dir/flags.make
source/CMakeFiles/Source.dir/Poisson_GMG.cpp.obj: C:/Users/Saha/Documents/GitHub/ColabCodes/ClionProject/source/Poisson_GMG.cpp
source/CMakeFiles/Source.dir/Poisson_GMG.cpp.obj: source/CMakeFiles/Source.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=C:\Users\Saha\Documents\GitHub\ColabCodes\ClionProject\final_build\CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object source/CMakeFiles/Source.dir/Poisson_GMG.cpp.obj"
	cd /d C:\Users\Saha\Documents\GitHub\ColabCodes\ClionProject\final_build\source && C:\MinGW\bin\g++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT source/CMakeFiles/Source.dir/Poisson_GMG.cpp.obj -MF CMakeFiles\Source.dir\Poisson_GMG.cpp.obj.d -o CMakeFiles\Source.dir\Poisson_GMG.cpp.obj -c C:\Users\Saha\Documents\GitHub\ColabCodes\ClionProject\source\Poisson_GMG.cpp

source/CMakeFiles/Source.dir/Poisson_GMG.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/Source.dir/Poisson_GMG.cpp.i"
	cd /d C:\Users\Saha\Documents\GitHub\ColabCodes\ClionProject\final_build\source && C:\MinGW\bin\g++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E C:\Users\Saha\Documents\GitHub\ColabCodes\ClionProject\source\Poisson_GMG.cpp > CMakeFiles\Source.dir\Poisson_GMG.cpp.i

source/CMakeFiles/Source.dir/Poisson_GMG.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/Source.dir/Poisson_GMG.cpp.s"
	cd /d C:\Users\Saha\Documents\GitHub\ColabCodes\ClionProject\final_build\source && C:\MinGW\bin\g++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S C:\Users\Saha\Documents\GitHub\ColabCodes\ClionProject\source\Poisson_GMG.cpp -o CMakeFiles\Source.dir\Poisson_GMG.cpp.s

source/CMakeFiles/Source.dir/ProblemDesc.cpp.obj: source/CMakeFiles/Source.dir/flags.make
source/CMakeFiles/Source.dir/ProblemDesc.cpp.obj: C:/Users/Saha/Documents/GitHub/ColabCodes/ClionProject/source/ProblemDesc.cpp
source/CMakeFiles/Source.dir/ProblemDesc.cpp.obj: source/CMakeFiles/Source.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=C:\Users\Saha\Documents\GitHub\ColabCodes\ClionProject\final_build\CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object source/CMakeFiles/Source.dir/ProblemDesc.cpp.obj"
	cd /d C:\Users\Saha\Documents\GitHub\ColabCodes\ClionProject\final_build\source && C:\MinGW\bin\g++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT source/CMakeFiles/Source.dir/ProblemDesc.cpp.obj -MF CMakeFiles\Source.dir\ProblemDesc.cpp.obj.d -o CMakeFiles\Source.dir\ProblemDesc.cpp.obj -c C:\Users\Saha\Documents\GitHub\ColabCodes\ClionProject\source\ProblemDesc.cpp

source/CMakeFiles/Source.dir/ProblemDesc.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/Source.dir/ProblemDesc.cpp.i"
	cd /d C:\Users\Saha\Documents\GitHub\ColabCodes\ClionProject\final_build\source && C:\MinGW\bin\g++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E C:\Users\Saha\Documents\GitHub\ColabCodes\ClionProject\source\ProblemDesc.cpp > CMakeFiles\Source.dir\ProblemDesc.cpp.i

source/CMakeFiles/Source.dir/ProblemDesc.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/Source.dir/ProblemDesc.cpp.s"
	cd /d C:\Users\Saha\Documents\GitHub\ColabCodes\ClionProject\final_build\source && C:\MinGW\bin\g++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S C:\Users\Saha\Documents\GitHub\ColabCodes\ClionProject\source\ProblemDesc.cpp -o CMakeFiles\Source.dir\ProblemDesc.cpp.s

source/CMakeFiles/Source.dir/GMG_solver.cpp.obj: source/CMakeFiles/Source.dir/flags.make
source/CMakeFiles/Source.dir/GMG_solver.cpp.obj: C:/Users/Saha/Documents/GitHub/ColabCodes/ClionProject/source/GMG_solver.cpp
source/CMakeFiles/Source.dir/GMG_solver.cpp.obj: source/CMakeFiles/Source.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=C:\Users\Saha\Documents\GitHub\ColabCodes\ClionProject\final_build\CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building CXX object source/CMakeFiles/Source.dir/GMG_solver.cpp.obj"
	cd /d C:\Users\Saha\Documents\GitHub\ColabCodes\ClionProject\final_build\source && C:\MinGW\bin\g++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT source/CMakeFiles/Source.dir/GMG_solver.cpp.obj -MF CMakeFiles\Source.dir\GMG_solver.cpp.obj.d -o CMakeFiles\Source.dir\GMG_solver.cpp.obj -c C:\Users\Saha\Documents\GitHub\ColabCodes\ClionProject\source\GMG_solver.cpp

source/CMakeFiles/Source.dir/GMG_solver.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/Source.dir/GMG_solver.cpp.i"
	cd /d C:\Users\Saha\Documents\GitHub\ColabCodes\ClionProject\final_build\source && C:\MinGW\bin\g++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E C:\Users\Saha\Documents\GitHub\ColabCodes\ClionProject\source\GMG_solver.cpp > CMakeFiles\Source.dir\GMG_solver.cpp.i

source/CMakeFiles/Source.dir/GMG_solver.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/Source.dir/GMG_solver.cpp.s"
	cd /d C:\Users\Saha\Documents\GitHub\ColabCodes\ClionProject\final_build\source && C:\MinGW\bin\g++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S C:\Users\Saha\Documents\GitHub\ColabCodes\ClionProject\source\GMG_solver.cpp -o CMakeFiles\Source.dir\GMG_solver.cpp.s

# Object files for target Source
Source_OBJECTS = \
"CMakeFiles/Source.dir/Grid.cpp.obj" \
"CMakeFiles/Source.dir/Poisson_FD.cpp.obj" \
"CMakeFiles/Source.dir/Poisson_GMG.cpp.obj" \
"CMakeFiles/Source.dir/ProblemDesc.cpp.obj" \
"CMakeFiles/Source.dir/GMG_solver.cpp.obj"

# External object files for target Source
Source_EXTERNAL_OBJECTS =

source/libSource.a: source/CMakeFiles/Source.dir/Grid.cpp.obj
source/libSource.a: source/CMakeFiles/Source.dir/Poisson_FD.cpp.obj
source/libSource.a: source/CMakeFiles/Source.dir/Poisson_GMG.cpp.obj
source/libSource.a: source/CMakeFiles/Source.dir/ProblemDesc.cpp.obj
source/libSource.a: source/CMakeFiles/Source.dir/GMG_solver.cpp.obj
source/libSource.a: source/CMakeFiles/Source.dir/build.make
source/libSource.a: source/CMakeFiles/Source.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=C:\Users\Saha\Documents\GitHub\ColabCodes\ClionProject\final_build\CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Linking CXX static library libSource.a"
	cd /d C:\Users\Saha\Documents\GitHub\ColabCodes\ClionProject\final_build\source && $(CMAKE_COMMAND) -P CMakeFiles\Source.dir\cmake_clean_target.cmake
	cd /d C:\Users\Saha\Documents\GitHub\ColabCodes\ClionProject\final_build\source && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles\Source.dir\link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
source/CMakeFiles/Source.dir/build: source/libSource.a
.PHONY : source/CMakeFiles/Source.dir/build

source/CMakeFiles/Source.dir/clean:
	cd /d C:\Users\Saha\Documents\GitHub\ColabCodes\ClionProject\final_build\source && $(CMAKE_COMMAND) -P CMakeFiles\Source.dir\cmake_clean.cmake
.PHONY : source/CMakeFiles/Source.dir/clean

source/CMakeFiles/Source.dir/depend:
	$(CMAKE_COMMAND) -E cmake_depends "MinGW Makefiles" C:\Users\Saha\Documents\GitHub\ColabCodes\ClionProject C:\Users\Saha\Documents\GitHub\ColabCodes\ClionProject\source C:\Users\Saha\Documents\GitHub\ColabCodes\ClionProject\final_build C:\Users\Saha\Documents\GitHub\ColabCodes\ClionProject\final_build\source C:\Users\Saha\Documents\GitHub\ColabCodes\ClionProject\final_build\source\CMakeFiles\Source.dir\DependInfo.cmake --color=$(COLOR)
.PHONY : source/CMakeFiles/Source.dir/depend

