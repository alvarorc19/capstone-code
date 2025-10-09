{
  description = "Developing in C++";

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
        devShells.default = pkgs.mkShell {
          name = "c++ with pixi";
          nativeBuildInputs = with pkgs; [
            clang
            clang-tools
          ];
          buildInputs = with pkgs; [
            ffmpeg # for video generation
          ];
          packages = with pkgs; [
            boost
            catch2
            cmake
            pixi
            pyright
          ];
        };
      };
    };
}
