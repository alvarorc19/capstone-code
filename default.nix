{ stdenv, cmake, boost, tomlplusplus, hdf5-cpp, catch2, ... }:

stdenv.mkDerivation {
  # 1. Package Metadata
  pname = "capstone";
  version = "1.0.0"; # Use an actual version number

  # 2. Source
  src = ./.; # Tells Nix to use the current directory as the source

  # 3. Build Tool Dependencies
  # These are the tools needed to perform the build (e.g., build system, compiler, linker)
  # Nix's stdenv recognizes 'cmake' here and automatically runs the configure/build steps.
  nativeBuildInputs = [
    cmake
  ];

  # 4. Runtime/Linker Dependencies (Libraries)
  # These are the C++ libraries your project links against (e.g., Boost, HDF5).
  # Nix will automatically set up the paths so CMake can find them.
  buildInputs = [
    boost
    tomlplusplus
    hdf5-cpp
    catch2 # Assuming you need this at link time for tests
  ];

  # 5. CMake Specific Configuration
  # Since your CMakeLists.txt doesn't use the default install commands,
  # we need to ensure the final executable ends up in the correct path ($out/bin).
  # We use postInstall to move the executable from the build directory to the installation directory.
  # The executable name 'ising' comes from your add_executable(ising ...) in CMake.
  postInstall = ''
    mkdir -p $out/bin
    # The executable is in the build directory because of RUNTIME_OUTPUT_DIRECTORY
    mv ising $out/bin/capstone-ising # Rename it for clarity in the install path
    # If your CMakeLists.txt had an 'install' step, this would not be necessary.
  '';
}
