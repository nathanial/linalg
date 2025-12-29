import Lake
open Lake DSL

package linalg where
  version := v!"0.1.0"
  leanOptions := #[
    ⟨`autoImplicit, false⟩,
    ⟨`relaxedAutoImplicit, false⟩
  ]

require crucible from git "https://github.com/nathanial/crucible" @ "v0.0.1"

@[default_target]
lean_lib Linalg where
  roots := #[`Linalg]

lean_lib LinalgTests where
  globs := #[.submodules `LinalgTests]

@[test_driver]
lean_exe linalg_tests where
  root := `LinalgTests.Main
