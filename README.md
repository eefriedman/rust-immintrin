# immintrin

This crate implements the names defined by Intel for SIMD and related instructions on x86
(functions with an `_mm_` prefix). These are useful for porting code, and because people
are familiar with them.  See https://github.com/rust-lang/rfcs/issues/1639 for discussion.

Currently contains signatures for all SSE/SSE2/SSE2/SSSE3 intrinsics. Intrinsics involving
MMX and non-temporal stores are not yet implemented; everything else should work. Requires a
nightly compiler for SIMD and various unstable intrinsics.

If anyone is looking to improve this crate, probably the most important thing at the moment
is some sort of testing infrastructure; currently, this crate has no tests.  Not precisely
sure what the tests would look like; probably some combination of checking the generated LLVM IR,
checking the generated assembly, and runtime tests.

It's not clear whether the signatures for intrinsics which take a constant integer should stay
the way they are, or be changed somehow.  Currently, we just accept that the parameter might
not be constant, and just make sure we generate code which will be optimized well if
the parameter is constant (for example, the match statement in `_mm_slli_si128` will fold down
to a single instruction).

Licensed under the Apache License, Version 2.0
<LICENSE-APACHE or
http://www.apache.org/licenses/LICENSE-2.0> or the MIT
license <LICENSE-MIT or http://opensource.org/licenses/MIT>,
at your option. All files in the project carrying such
notice may not be copied, modified, or distributed except
according to those terms.
