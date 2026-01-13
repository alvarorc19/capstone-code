{
  description = "Developing in C++ and python plotting";

  inputs = {
    flake-parts.url = "github:hercules-ci/flake-parts";
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";
  };

  outputs = inputs @ {flake-parts, ...}:
    flake-parts.lib.mkFlake {inherit inputs;} {
      systems = ["x86_64-linux" "aarch64-linux" "aarch64-darwin" "x86_64-darwin"];
      perSystem = {
        config,
        self',
        inputs',
        pkgs,
        system,
        ...
      }: {
        packages.default = pkgs.callPackage ./default.nix {};
        devShells.default = pkgs.mkShell {
          name = "c++ with pixi";
          nativeBuildInputs = with pkgs; [
            clang
          ];
          buildInputs = with pkgs; [
            ffmpeg # for video generation
            pyright
            clang-tools
            # hdfview
          ];
          packages = with pkgs; [
            boost
            doxygen
            tomlplusplus
            # hdf5-cpp
            highfive
            catch2
            cmake
            pixi
          ];
          shellHook = ''
            echo "pixi shell and c++ environment initialised"
            clear
            exec pixi shell
            exec zsh
          '';
        };
      };
    };
}
