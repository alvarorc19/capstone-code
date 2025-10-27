{
  stdenv,
  cmake,
  boost,
  tomlplusplus,
  hdf5-cpp,
  catch2,
  ...
}:
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
}
