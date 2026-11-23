{
  description = "Julia environment for L314 DSP Project";

  inputs.nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";

  outputs = { nixpkgs, ... }:
    let
      systems = [ "x86_64-linux" "aarch64-linux" ];
      forAllSystems = nixpkgs.lib.genAttrs systems;
    in {
      devShells = forAllSystems (system:
        let
          pkgs = nixpkgs.legacyPackages.${system};
          julia = pkgs.julia-bin.withPackages [
            "DSP"
            "FFTW"
            "Plots"
            "Pluto"
            "Images"
            "ImageShow"
          ];
        in {
          default = pkgs.mkShell {
            packages = [
              julia
              pkgs.poppler-utils
            ];
          };
        });
    };
}
